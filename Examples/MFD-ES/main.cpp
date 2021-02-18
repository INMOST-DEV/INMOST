#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace INMOST;

#ifndef M_PI
#define M_PI 3.141592653589
#endif

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

bool report_residual = false;
bool write_steps = false;
bool two_point = false;

//shortcuts
typedef Storage::bulk bulk;
typedef Storage::real real;
typedef Storage::integer integer;
typedef Storage::enumerator enumerator;
typedef Storage::real_array real_array;
typedef Storage::var_array var_array;


variable refU(real xyz[3], real nu, real t)
{
	unknown x(xyz[0],0);
	unknown y(xyz[1],1);
	unknown z(xyz[2],2);
	const real a = M_PI/4.0;
	const real d = M_PI/2.0;
	return -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-nu*d*d*t);
}

variable refV(real xyz[3], real nu, real t)
{
	unknown x(xyz[0],0);
	unknown y(xyz[1],1);
	unknown z(xyz[2],2);
	const real a = M_PI/4.0;
	const real d = M_PI/2.0;
	return -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-nu*d*d*t);
}

variable refW(real xyz[3], real nu, real t)
{
	unknown x(xyz[0],0);
	unknown y(xyz[1],1);
	unknown z(xyz[2],2);
	const real a = M_PI/4.0;
	const real d = M_PI/2.0;
	return -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x))*exp(-nu*d*d*t);
}

real refP(real xyz[3], real nu, real t)
{
	real x = xyz[0];
	real y = xyz[1];
	real z = xyz[2];
	const real a = M_PI/4.0;
	const real d = M_PI/2.0;
	return -a*a/2.0*(exp(2*a*x)+exp(2*a*y)+exp(2*a*z)+
					 2*sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))+
					 2*sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(z+x))+
					 2*sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y)))*exp(-2*nu*d*d*t);
}

// (done)todo: adapt for parallel use
void CheckResidual(Mesh & m, const TagReal & tag_P, const TagReal & tag_nU, real visc, real T, MarkerType mrk)
{
	{
		TagReal      tag_refP = m.CreateTag("refP" ,DATA_REAL,CELL|FACE,NONE,1);
		real p_ref_min = 1.0e20, p_ref_max = -1.0e20;
		real p_min = 1.0e20, p_max = -1.0e20;
		for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype) )
			for(integer it = 0; it < m.LastLocalID(etype); ++it) if( m.isValidElement(etype,it) )
			{
				Element e = m.ElementByLocalID(etype,it);
				if( e.GetStatus() == Element::Ghost ) continue;
				real x[3],p;
				e.Centroid(x);
				p = refP(x,visc,T);
				if( p < p_ref_min ) p_ref_min = p;
				if( p > p_ref_max ) p_ref_max = p;
				p = tag_P[e];
				if( p < p_min ) p_min = p;
				if( p > p_max ) p_max = p;
			}
		p_max = m.AggregateMax(p_max);
		p_min = m.AggregateMin(p_min);
		p_ref_max = m.AggregateMax(p_ref_max);
		p_ref_min = m.AggregateMin(p_ref_min);
		for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype) )
			for(integer it = 0; it < m.LastLocalID(etype); ++it) if( m.isValidElement(etype,it) )
			{
				Element e = m.ElementByLocalID(etype,it);
				if( e.GetStatus() == Element::Ghost ) continue;
				real x[3];
				e.Centroid(x);
				tag_refP[e] = refP(x,visc,T) - p_ref_min;
			}
		if( m.GetProcessorRank() == 0 )
		{
			std::cout << "p_ref " << p_ref_min << ":" << p_ref_max << std::endl;
			std::cout << "p     " << p_min << ":" << p_max << std::endl;
		}
		TagReal      tag_E   = m.CreateTag("ERROR_P" ,DATA_REAL,CELL|FACE,NONE,1);
		TagRealArray tag_EU  = m.CreateTag("ERROR_U" ,DATA_REAL,CELL,NONE,3);
		TagRealArray tag_U  =  m.CreateTag("U" ,DATA_REAL,CELL,NONE,3);
		TagRealArray tag_refU  =  m.CreateTag("refU" ,DATA_REAL,CELL,NONE,3);
		TagReal      tag_EnU = m.CreateTag("ERROR_nU",DATA_REAL,FACE,NONE,1);
		real C, L2, volume, Cu[3], L2u[3], Cnu, L2nu;
		C = L2 = volume = 0.0;
		Cu[0] = Cu[1] = Cu[2] = L2u[0] = L2u[1] = L2u[2] = 0;
		for(integer q = 0; q < m.CellLastLocalID(); ++q ) if( m.isValidCell(q) )
		{
			Cell c = m.CellByLocalID(q);
			if( c.GetStatus() == Element::Ghost ) continue;
			real x[3], z[3], n[3];
			c.Centroid(x);
			real err = tag_P[c] - p_min - tag_refP[c];
			real vol = c.Volume();
			if( C < fabs(err) ) C = fabs(err);
			L2 += err*err*vol;
			volume += vol;
			tag_E[c] = err;
			
			real V[3],Q[3]; //restored velocity
			V[0] = V[1] = V[2] = 0;
			Q[0] = Q[1] = Q[2] = 0;
			real a,s,v;
			ElementArray<Face> cfaces = c.getFaces();
			for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt)
			{
				//cfaces[kt].FixNormalOrientation();
				//if( !cfaces[kt].FixNormalOrientation() ) std::cout << "Orientation is bad" << std::endl;
				cfaces[kt].Centroid(z);
				cfaces[kt].UnitNormal(n);
				a = cfaces[kt].Area();
				s = cfaces[kt].FaceOrientedOutside(c) ? 1.0 : -1.0;
				s*= cfaces[kt].GetMarker(mrk) ? -1.0 : 1.0;
				v = tag_nU[cfaces[kt]]*s*a;
				//v = n[0]*refU(z,visc,T).GetValue()
				//  + n[1]*refV(z,visc,T).GetValue()
				//  + n[2]*refW(z,visc,T).GetValue();
				//v*= a*s;
				V[0] += v*(z[0]-x[0]);
				V[1] += v*(z[1]-x[1]);
				V[2] += v*(z[2]-x[2]);
			}
			V[0] /= c.Volume();
			V[1] /= c.Volume();
			V[2] /= c.Volume();
			
			tag_U[c][0] = V[0];
			tag_U[c][1] = V[1];
			tag_U[c][2] = V[2];
			
			tag_refU[c][0] = refU(x,visc,T).GetValue();
			tag_refU[c][1] = refV(x,visc,T).GetValue();
			tag_refU[c][2] = refW(x,visc,T).GetValue();
			
			err = V[0] - tag_refU[c][0];
			if( Cu[0] < fabs(err) ) Cu[0] = fabs(err);
			L2u[0] += err*err*vol;
			tag_EU[c][0] = err;
			
			err = V[1] - tag_refU[c][1];
			if( Cu[1] < fabs(err) ) Cu[1] = fabs(err);
			L2u[1] += err*err*vol;
			tag_EU[c][1] = err;
			
			err = V[2] - tag_refU[c][2];
			if( Cu[2] < fabs(err) ) Cu[2] = fabs(err);
			L2u[2] += err*err*vol;
			tag_EU[c][2] = err;
		}
		m.ExchangeData(tag_U,CELL,0);
		m.ExchangeData(tag_refU,CELL,0);
		volume = m.Integrate(volume);
		L2 = sqrt(m.Integrate(L2)/volume);
		m.Integrate(L2u,3);
		L2u[0] = sqrt(L2u[0]/volume);
		L2u[1] = sqrt(L2u[1]/volume);
		L2u[2] = sqrt(L2u[2]/volume);
		C = m.AggregateMax(C);
		m.AggregateMax(Cu,3);
		if( m.GetProcessorRank() == 0 )
		{
			std::cout << "Pressure error on cells, C-norm " << C << " L2-norm " << L2 << std::endl;
			std::cout << "Reconstructed U error on cells, C-norm " << Cu[0] << " L2-norm " << L2u[0] << std::endl;
			std::cout << "Reconstructed V error on cells, C-norm " << Cu[1] << " L2-norm " << L2u[1] << std::endl;
			std::cout << "Reconstructed W error on cells, C-norm " << Cu[2] << " L2-norm " << L2u[2] << std::endl;
		}
		
		C = L2 = volume = 0.0;
		Cnu = L2nu = 0.0;
		for( int q = 0; q < m.FaceLastLocalID(); ++q ) if( m.isValidFace(q) )
		{
			Face f = m.FaceByLocalID(q);
			
			if( f.GetStatus() == Element::Ghost ) continue;
			
			real x[3],y[3],n[3];
			f.Centroid(x);
			f.UnitNormal(n);
			f.BackCell().Centroid(y);
			real H = n[0]*(x[0]-y[0])+n[1]*(x[1]-y[1])+n[2]*(x[2]-y[2]);
			if( f.FrontCell().isValid() )
			{
				f.FrontCell().Centroid(y);
				H += n[0]*(y[0]-x[0])+n[1]*(y[1]-x[1])+n[2]*(y[2]-x[2]);
			}
			real err = tag_P[f] - p_min - tag_refP[f];
			real vol = f.Area()*H;
			if( C < fabs(err) ) C = fabs(err);
			L2 += err*err*vol;
			volume += vol;
			tag_E[f] = err;
			
			real V[3], nV;
			V[0] = refU(x,visc,T).GetValue();
			V[1] = refV(x,visc,T).GetValue();
			V[2] = refW(x,visc,T).GetValue();
			nV = (n[0]*V[0]+n[1]*V[1]+n[2]*V[2])*(f.GetMarker(mrk)?-1:1);
			err = tag_nU[f] - nV;
			
			if( Cnu < fabs(err) ) Cnu = fabs(err);
			L2nu += err*err*vol;
			volume += vol;
			tag_EnU[f] = err;
		}
		volume = m.Integrate(volume);
		L2 = sqrt(m.Integrate(L2)/volume);
		L2nu = sqrt(m.Integrate(L2nu)/volume);
		C = m.AggregateMax(C);
		Cnu = m.AggregateMax(Cnu);
		if( m.GetProcessorRank() == 0 )
		{
			if( !two_point ) std::cout << "Pressure error on faces, C-norm " << C << " L2-norm " << L2 << std::endl;
			std::cout << "Normal velocity error on faces, C-norm " << Cnu << " L2-norm " << L2nu << std::endl;
		}
		
		real div_max = -1.0e20, div_min = 1.0e20, div_int = 0;
		for(integer it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) ) //loop over cells
		{
			Cell c = m.CellByLocalID(it);
			if( c.GetStatus() == Element::Ghost ) continue;
			ElementArray<Face> cfaces = c.getFaces(); //obtain faces of the cell
			real div = 0, a, s;
			for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt)
			{
				a = cfaces[kt].Area();
				s = cfaces[kt].FaceOrientedOutside(c)?1.0:-1.0;
				s*= cfaces[kt].GetMarker(mrk) ? -1.0 : 1.0;
				div += tag_nU[cfaces[kt]]*a*s;
			}
			if( div > div_max ) div_max = div;
			if( div < div_min ) div_min = div;
			div_int += div;
		}
		div_max = m.AggregateMax(div_max);
		div_min = m.AggregateMin(div_min);
		div_int = m.Integrate(div_int);
		if( m.GetProcessorRank() == 0 )
		{
			std::cout << "Divergence of velocity: " << div_min << ":" << div_max << " integral " << div_int << std::endl;
		}
	}
}
		
int main(int argc,char ** argv)
{
    Solver::Initialize(&argc,&argv,""); // Initialize the solver and MPI activity
#if defined(USE_PARTITIONER)
    Partitioner::Initialize(&argc,&argv); // Initialize the partitioner activity
#endif
    if( argc > 1 )
    {
		
		//parameters
		
		real dT = 1.0/50.0;
		real T = 1.0/10.0;
		real visc = 0.0;
		
		if( argc > 2 ) dT = atof(argv[2]);
		if( argc > 3 ) T = atof(argv[3]);
		if( argc > 4 ) visc = atof(argv[4]);
		
		if( visc > 0 )
		{
			std::cout << "Viscosity is not calculated right now." << std::endl;
			visc = 0;
		}
		
        bool repartition = false; // Is it required to redistribute the mesh?
        Mesh m; // Create an empty mesh
        { // Load the mesh
            m.SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
            if( m.GetProcessorRank() == 0 ) std::cout << "Processors: " << m.GetProcessorsNumber() << std::endl;
            if( m.isParallelFileFormat(argv[1]) ) //The format is
            {
                m.Load(argv[1]); // Load mesh from the parallel file format
                repartition = true; // Ask to repartition the mesh
            }
            else if( m.GetProcessorRank() == 0 ) m.Load(argv[1]); // Load mesh from the serial file format
        }

		
#if defined(USE_PARTITIONER)
        if (m.GetProcessorsNumber() > 1 )//&& !repartition) // Currently only non-distributed meshes are supported by Inner_RCM partitioner
        {
            { // Compute mesh partitioning
                Partitioner p(&m); //Create Partitioning object
                p.SetMethod(Partitioner::INNER_KMEANS,repartition ? Partitioner::Repartition : Partitioner::Partition); // Specify the partitioner
                p.Evaluate(); // Compute the partitioner and store new processor ID in the mesh
            }

            { //Distribute the mesh
                m.Redistribute(); // Redistribute the mesh data
                m.ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
            }
        }
#endif
		//map grid into [-1,1] box
		real xmax[3] = {-1.0e20,-1.0e20,-1.0e20};
		real xmin[3] = {+1.0e20,+1.0e20,+1.0e20};
		for(integer it = 0; it < m.NodeLastLocalID(); ++it) if( m.isValidNode(it) )
		{
			real_array c = m.NodeByLocalID(it).Coords();
			for(int k = 0; k < 3; ++k)
			{
				if( xmin[k] > c[k] ) xmin[k] = c[k];
				if( xmax[k] < c[k] ) xmax[k] = c[k];
			}
		}
		m.AggregateMax(xmax,3);
		m.AggregateMin(xmin,3);
		
		for(integer it = 0; it < m.NodeLastLocalID(); ++it) if( m.isValidNode(it) )
		{
			real_array c = m.NodeByLocalID(it).Coords();
			for(int k = 0; k < 3; ++k)
				c[k] = 2.0*(c[k]-xmin[k])/(xmax[k]-xmin[k])-1.0;
		}
		
		
		
        { // prepare geometrical data on the mesh
            Mesh::GeomParam table;
            table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
            table[NORMAL]      = FACE;        //Compute normals
            table[ORIENTATION] = FACE;        //Check and fix normal orientation
            table[MEASURE]     = CELL | FACE; //Compute volumes and areas
            //table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
			m.RemoveGeometricData(table);
            m.PrepareGeometricData(table); //Ask to precompute the data
        }
		
		
        // data tags for
		TagReal tag_nU; // Normal component of velocity
		TagReal tag_nUold; // on previous step
        TagReal tag_P;  // Pressure
		TagReal tag_Q;  // Pressure update
        TagReal tag_F;  // Forcing term
        TagRealArray tag_W; // Gradient matrix
		TagInteger tag_i; // row number for current face in back cell matrix
		
		
        if( m.GetProcessorsNumber() > 1 ) //skip for one processor job
        { // Exchange ghost cells
            m.ExchangeGhost(1,FACE); // Produce layer of ghost cells
        }
		
		
		MarkerType nrm_ornt = m.CreateMarker();
		MarkerType boundary = m.CreateMarker();
		
		m.MarkNormalOrientation(nrm_ornt);
		m.MarkBoundaryFaces(boundary);
		
        { //initialize data
            tag_P = m.CreateTag("PRESSURE",DATA_REAL,CELL|FACE,NONE,1); // Create a new tag for the pressure
			tag_W = m.CreateTag("W",DATA_REAL,CELL,NONE);
			tag_i = m.CreateTag("row",DATA_INTEGER,FACE,NONE,1);
			tag_F = m.CreateTag("FORCE",DATA_REAL,CELL,NONE,1);
			tag_Q = m.CreateTag("PUPDATE",DATA_REAL,CELL|FACE,NONE,1);
			tag_nU = m.CreateTag("NORMAL_VELOCITY",DATA_REAL,FACE,NONE,1);
			tag_nUold = m.CreateTag("NORMAL_VELOCITY_OLD",DATA_REAL,FACE,NONE,1);
			//Assemble gradient matrix W on cells
#if defined(USE_OMP)
#pragma omp parallel for
#endif
            for(integer it = 0; it < m.CellLastLocalID(); ++it ) if( m.isValidCell(it) )
            {
                Cell c = m.CellByLocalID(it);
				ElementArray<Face> faces = c.getFaces(); //obtain faces of the cell
				Matrix<real,real_array> W(tag_W[c]);
				W.Resize(faces.size(),faces.size());
				real xP[3]; //center of the cell
				real yF[3]; //center of the face
				real nF[3]; //normal to the face
				real aF; //area of the face
				c.Centroid(xP);
				rMatrix N(faces.size(),3), R(faces.size(),3); //big gradient matrix, co-normals, directions
				for(INMOST_DATA_ENUM_TYPE k = 0; k < faces.size(); ++k) //loop over faces
				{
					aF = faces[k].Area();
					faces[k].Centroid(yF);
					faces[k].OrientedUnitNormal(c,nF);
					// assemble matrix of directions
					R(k,0) = (yF[0]-xP[0])*aF;
					R(k,1) = (yF[1]-xP[1])*aF;
					R(k,2) = (yF[2]-xP[2])*aF;
					// assemble matrix of co-normals
					N(k,0) = nF[0];
					N(k,1) = nF[1];
					N(k,2) = nF[2];
					// remember position
					if( faces[k].BackCell() == c )
						tag_i[faces[k]] = k;
				} //end of loop over faces
				W = N*(N.Transpose()*R).Invert()*N.Transpose(); //stability part
				W += (rMatrix::Unit(faces.size()) - R*(R.Transpose()*R).Invert()*R.Transpose())*(4.0/(faces.size())*W.Trace());
            } //end of loop over cells
			
			//initialize normal velocity
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(integer it = 0; it < m.FaceLastLocalID(); ++it) if( m.isValidFace(it) )
			{
				Face f = m.FaceByLocalID(it);
				real x[3], n[3], u[3];
				f.UnitNormal(n);
				f.Centroid(x);
				u[0] = refU(x,visc,0).GetValue();
				u[1] = refV(x,visc,0).GetValue();
				u[2] = refW(x,visc,0).GetValue();
				tag_nUold[f] = tag_nU[f] = (u[0]*n[0]+u[1]*n[1]+u[2]*n[2])*(f.GetMarker(nrm_ornt)?-1:1);
				//pressure
				tag_P[f] = refP(x,visc,0);
			}
			
			
			//initialize pressure
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(integer it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
			{
				Cell c = m.CellByLocalID(it);
				real x[3];
				c.Centroid(x);
				tag_P[c] = refP(x,visc,0);
			}
			
			//CheckResidual(m,tag_P,tag_nU,visc,0,nrm_ornt);
			
			m.ExchangeData(tag_nU,FACE,0);
			m.ExchangeData(tag_nUold,FACE,0);
			
			//CheckResidual(m,tag_P,tag_nU,visc,0,nrm_ornt);
			
        } //end of initialize data
		
		

        if( m.GetProcessorRank() == 0 ) std::cout << "Initialization done" << std::endl;
        { //Main loop for problem solution
            Automatizator aut; // declare class to help manage unknowns
			
			dynamic_variable Q(aut,aut.RegisterTag(tag_Q,two_point ? CELL : (CELL|FACE))); //register pressure as primary unknown
            aut.EnumerateEntries(); //enumerate all primary variables
			
			std::cout << "Enumeration done, size " << aut.GetLastIndex() - aut.GetFirstIndex() << std::endl;

            Residual R("",aut.GetFirstIndex(),aut.GetLastIndex());
            Sparse::LockService Locks(aut.GetFirstIndex(),aut.GetLastIndex());
            Sparse::AnnotationService Text(aut.GetFirstIndex(),aut.GetLastIndex());
            Sparse::Vector Update  ("",aut.GetFirstIndex(),aut.GetLastIndex()); //vector for update
            {//Annotate matrix
                for( int q = 0; q < m.CellLastLocalID(); ++q ) if( m.isValidCell(q) )
                {
                    Cell cell = m.CellByLocalID(q);
                    if( cell.GetStatus() != Element::Ghost )
                        Text.SetAnnotation(Q.Index(cell),"Cell-centered pressure value");
                }
                if( !two_point ) for( int q = 0; q < m.FaceLastLocalID(); ++q ) if( m.isValidFace(q) )
                {
                    Face face = m.FaceByLocalID(q);
                    if( face.GetStatus() != Element::Ghost )
                    {
                        //if( tag_BC.isValid() && face.HaveData(tag_BC) )
                        //    Text.SetAnnotation(P.Index(face),"Pressure guided by boundary condition");
                        //else
                            Text.SetAnnotation(Q.Index(face),"Interface pressure");
                    }
                }
            }

            if( m.GetProcessorRank() == 0 )  std::cout << "Matrix was annotated" << std::endl;
			
			//checks
			/*
			for(integer it = 0; it < m.FaceLastLocalID(); ++it) if( m.isValidFace(it) )
			{
				Face f = m.FaceByLocalID(it);
				Cell c = f.BackCell();
				ElementArray<Face> cfaces = c.getFaces();
				Matrix<real,real_array> W(tag_W[c],cfaces.size(),cfaces.size());
				real ngrad[3] = {0,0,0};
				real n[3], xf[3], xc[3];
				f.UnitNormal(n);
				c.Centroid(xc);
				for(integer kt = 0; kt < cfaces.size(); ++kt)
				{
					f.Centroid(xf);
					ngrad[0] += W(tag_i[f],kt)*(xf[0] - xc[0])*cfaces[kt].Area();
					ngrad[1] += W(tag_i[f],kt)*(xf[1] - xc[1])*cfaces[kt].Area();
					ngrad[2] += W(tag_i[f],kt)*(xf[2] - xc[2])*cfaces[kt].Area();
				}
				std::cout << "1/area: " << 1.0/f.Area() << std::endl;
				std::cout << "normal: " << n[0] << " " << n[1] << " " << n[2] << std::endl;
				std::cout << "ngrad:  " << ngrad[0] << " " << ngrad[1] << " " << ngrad[2] << std::endl;
			}
			 */
			
			
			int step = 0;
			real alpha = 1, beta = -1, gamma = 0;
			for(real t = 0; t < T-1.0e-9; t+=dT)
			{
				//if( step == 1 ) break;
				if( t > T ) t = T;
				
				if( report_residual )
				{
					if( m.GetProcessorRank() == 0 )  std::cout << "Before velocity update " << std::endl;
					CheckResidual(m,tag_P,tag_nU,visc,t,nrm_ornt);
				}
				
				
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer it = 0; it < m.FaceLastLocalID(); ++it) if( m.isValidFace(it) )
				{
					Face f = m.FaceByLocalID(it);
					if( f.GetStatus() == Element::Ghost ) continue;
					real x[3], n[3], V[3], nV, U, Uold, s;
					s = f.GetMarker(nrm_ornt) ? -1 : 1;
					f.UnitNormal(n);
					f.Centroid(x);
					U = tag_nU[f];
					Uold = tag_nUold[f];
					if( f.GetMarker(boundary) ) //Boundary condition
					{
						f.Centroid(x);
						V[0] = refU(x,visc,t+dT).GetValue();
						V[1] = refV(x,visc,t+dT).GetValue();
						V[2] = refW(x,visc,t+dT).GetValue();
						nV = (n[0]*V[0]+n[1]*V[1]+n[2]*V[2])*s;
						tag_nU[f] = nV;
						tag_nUold[f] = U;
					}
					else //Update with advection and pressure gradient
					{
						real ngradP = 0;
						if( two_point )
						{
							if( f.FrontCell().isValid() )
							{
								real xf[3],xb[3],n[3],d;
								f.UnitNormal(n);
								f.FrontCell().Centroid(xf);
								f.BackCell().Centroid(xb);
								d = n[0]*(xf[0]-xb[0]) + n[1]*(xf[1]-xb[1]) + n[2]*(xf[2]-xb[2]);
								ngradP = (tag_P[f.FrontCell()] - tag_P[f.BackCell()])/d;
							}
						}
						else
						{
							Cell c = f.BackCell();
							ElementArray<Face> cfaces = c.getFaces();
							Matrix<real,real_array> W(tag_W[c],cfaces.size(),cfaces.size());
							for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt)
								ngradP += W(tag_i[f],kt)*(tag_P[cfaces[kt]] - tag_P[c])*cfaces[kt].Area();
						}
						variable u = refU(x,visc,t+dT);
						variable v = refV(x,visc,t+dT);
						variable w = refW(x,visc,t+dT);
						V[0] = u.GetValue()*u.GetRow()[0] + v.GetValue()*u.GetRow()[1] + w.GetValue()*u.GetRow()[2];
						V[1] = u.GetValue()*v.GetRow()[0] + v.GetValue()*v.GetRow()[1] + w.GetValue()*v.GetRow()[2];
						V[2] = u.GetValue()*w.GetRow()[0] + v.GetValue()*w.GetRow()[1] + w.GetValue()*w.GetRow()[2];
						nV = (n[0]*V[0]+n[1]*V[1]+n[2]*V[2])*s;
						tag_nU[f] = -(beta*U + gamma*Uold + (nV + ngradP)*dT)/alpha;
						tag_nUold[f] = U;
					}
				}
				
				m.ExchangeData(tag_nU,FACE,0);
				m.ExchangeData(tag_nUold,FACE,0);
				
				if( report_residual )
				{
					if( m.GetProcessorRank() == 0 )  std::cout << "After velocity update " << std::endl;
					CheckResidual(m,tag_P,tag_nU,visc,t+dT,nrm_ornt);
				}
				
				R.Clear(); //clean up the residual
				
				//clean up Q just in case, check does this affect the result
				//for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype) )
				//	for(integer it = 0; it < m.LastLocalID(etype); ++it) if( m.isValidElement(etype,it) )
				//		tag_Q[m.ElementByLocalID(etype,it)] = 0;
				
				
				Automatizator::MakeCurrent(&aut);
				
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) ) //loop over cells
				{
					Cell c = m.CellByLocalID(it);
					ElementArray<Face> cfaces = c.getFaces(); //obtain faces of the cell
					vMatrix FLUX(cfaces.size(),1); //computed flux on faces
					if( two_point )
					{
						for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt)
						{
							Face f = cfaces[kt];
							if( f.FrontCell().isValid() )
							{
								real xf[3],xb[3],n[3],d,s;
								f.UnitNormal(n);
								f.FrontCell().Centroid(xf);
								f.BackCell().Centroid(xb);
								d = n[0]*(xf[0]-xb[0]) + n[1]*(xf[1]-xb[1]) + n[2]*(xf[2]-xb[2]);
								s = f.FaceOrientedOutside(c) ? 1.0 : -1.0;
								FLUX(kt,0) = s*(Q[f.FrontCell()] - Q[f.BackCell()])/d;
							}
						}
					}
					else
					{
						Matrix<real,real_array> W(tag_W[c],cfaces.size(),cfaces.size()); //Matrix for gradient
						vMatrix dQ(cfaces.size(),1); //vector of pressure differences on faces
						for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt)
							dQ(kt,0) = (Q[cfaces[kt]] - Q[c])*cfaces[kt].Area();
						FLUX = W*dQ; //fluxes on faces
					}
					
					if( c.GetStatus() != Element::Ghost )
					{
						real div = 0, a, s;
						for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt)
						{
							a = cfaces[kt].Area();
							s = cfaces[kt].FaceOrientedOutside(c)?1.0:-1.0;
							s*= cfaces[kt].GetMarker(nrm_ornt) ? -1 : 1;
							div += tag_nU[cfaces[kt]]*a*s;
						}
						tag_F[c] = div/dT;
						R[Q.Index(c)] = alpha*div/dT;
						for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt) //loop over faces of
							R[Q.Index(c)] += FLUX(kt,0)*cfaces[kt].Area();
					}
					if( !two_point ) for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt) //loop over faces of current cell
					{
						if( cfaces[kt].GetStatus() == Element::Ghost ) continue;
						int index = Q.Index(cfaces[kt]);
						Locks.Lock(index);
						R[index] += FLUX(kt,0);
						Locks.UnLock(index);
					}
				} //end of loop over cells
				
				{
					//set one of the boundary pressure corrections to be zero
					//(done) todo: in parallel processors should calculate number of local boundary interfaces
					//if all are neumann then the processor with minimum number should set one of the boundary
					//faces to zero
					//
					// todo: in future detect dirichlet bc, if they are present this is not needed
					if( m.GetProcessorRank() == 0 )
					{
						for(integer it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
						{
							Cell c = m.CellByLocalID(it);
							if( c.GetStatus() != Element::Ghost  )
							{
								R[Q.Index(c)] = Q[c] - 0.0;
								break;
							}
						}
					}
				}
				
				Automatizator::RemoveCurrent();
				
				//std::cout << "residual: " << R.Norm() << std::endl;
				
				//Solver S(Solver::SUPERLU);
				//Solver S(Solver::INNER_ILU2);
				Solver S(Solver::INNER_MPTILUC);
				//Solver S(Solver::INNER_DDPQILUC);
				S.SetParameter("relative_tolerance", "1.0e-14");
				S.SetParameter("absolute_tolerance", "1.0e-12");
				S.SetParameter("drop_tolerance", "5.0e-2");
				S.SetParameter("reuse_tolerance", "1.0e-3");
				S.SetParameter("schwartz_overlap", "2");
				S.SetMatrix(R.GetJacobian());;
				if( S.Solve(R.GetResidual(),Update) )
				{
					
					//real Q_min = 1.0e54, Q_max = -1.0e54;
					for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype) )
					{
						if( two_point && etype == FACE ) continue;
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer it = 0; it < m.LastLocalID(etype); ++it) if( m.isValidElement(etype,it) )
						{
							Element e = m.ElementByLocalID(etype,it);
							if( e.GetStatus() == Element::Ghost ) continue;
							tag_Q[e] -= Update[Q.Index(e)];
							//if( tag_Q[e] < Q_min ) Q_min = tag_Q[e];
							//if( tag_Q[e] > Q_max ) Q_max = tag_Q[e];
						}
					}
					
					//std::cout << "Q: " << Q_min << ":" << Q_max << std::endl;
					/*
					std::cout << "Q_min:" << Q_min << std::endl;
					for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype) )
						for(integer it = 0; it < m.LastLocalID(etype); ++it) if( m.isValidElement(etype,it) )
							tag_Q[m.ElementByLocalID(etype,it)] -= Q_min;
					 */
					
					m.ExchangeData(tag_Q, CELL|FACE, 0);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for(integer it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
					{
						Cell c = m.CellByLocalID(it);
						tag_P[c] -= tag_Q[c];
					}
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for(integer it = 0; it < m.FaceLastLocalID(); ++it) if( m.isValidFace(it) )
					{
						Face f = m.FaceByLocalID(it);
						//update pressure
						tag_P[f] -= tag_Q[f];
						//update for velocity field
						real ngradP = 0;
						if( two_point )
						{
							if( f.FrontCell().isValid() )
							{
								real xf[3],xb[3],n[3],d;
								f.UnitNormal(n);
								f.FrontCell().Centroid(xf);
								f.BackCell().Centroid(xb);
								d = n[0]*(xf[0]-xb[0]) + n[1]*(xf[1]-xb[1]) + n[2]*(xf[2]-xb[2]);
								ngradP = (tag_Q[f.FrontCell()] - tag_Q[f.BackCell()])/d;
							}
						}
						else
						{
							Cell c = f.BackCell();
							ElementArray<Face> cfaces = c.getFaces();
							Matrix<real,real_array> W(tag_W[c],cfaces.size(),cfaces.size());
							for(INMOST_DATA_ENUM_TYPE kt = 0; kt < cfaces.size(); ++kt)
								ngradP += W(tag_i[f],kt)*(tag_Q[cfaces[kt]] - tag_Q[c])*cfaces[kt].Area();
						}
						tag_nU[f] += ngradP*dT/alpha;
						
					}
					m.ExchangeData(tag_P, CELL|FACE, 0);
					m.ExchangeData(tag_nU, FACE, 0);
					
					if( report_residual )
					{
						if( m.GetProcessorRank() == 0 )  std::cout << "After projection" << std::endl;
						CheckResidual(m,tag_P,tag_nU,visc,t+dT,nrm_ornt);
					}
					
					if( write_steps )
					{
						std::stringstream str;
						str << "step" << step;
						if( m.GetProcessorsNumber() == 1 )
							str << ".vtk";
						else
							str << ".pvtk";
						m.Save(str.str());
					}
					
				}
				else
				{
					if( m.GetProcessorRank() == 0 )  std::cout << "Unable to solve: " << S.ReturnReason() << std::endl;
					break;
				}
				
				if( m.GetProcessorRank() == 0 )  std::cout << "step " << step << " time " << t+dT << "/" << T << std::endl;
				step++;
				alpha = 3.0/2.0;
				beta = -4.0/2.0;
				gamma = 1.0/2.0;
			} //end time step
			
		}
		
		
		CheckResidual(m,tag_P,tag_nU,visc,T,nrm_ornt);
		

        if( m.GetProcessorsNumber() == 1 )
            m.Save("out.vtk");
        else
            m.Save("out.pvtk");
		
		
		m.ReleaseMarker(boundary,FACE);
		m.ReleaseMarker(nrm_ornt,FACE);
    }
    else
    {
        std::cout << argv[0] << " mesh_file [dt=1/50] [T=1/10] [visc=0]" << std::endl;
    }

#if defined(USE_PARTITIONER)
    Partitioner::Finalize(); // Finalize the partitioner activity
#endif
    Solver::Finalize(); // Finalize solver and close MPI activity
    return 0;
}
