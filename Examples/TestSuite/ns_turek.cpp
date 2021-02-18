#include "inmost.h"

using namespace INMOST;

// use ADMFD or ADVDIFF to solve for phi

//setup BC
//setup on mesh:
//
// FORCE - 3 entries
//
// BOUNDARY_CONDITION_VELOCITY - 7 entries
// r = (a5,a6,a7)
// n.(a1*u + a2*t) = n.r
// (I - nn.)(a3*u + a4*t) = (I-nn.)r
//
// BOUNDARY_CONDITION_PRESSURE - 1 entry
//
// REFERENCE_VELOCITY - 3 entries
//
// REFERENCE_PRESSURE - 1 entry
//
// BOUNDARY_CONDTION - 3 entries (for calculation of phi and measures of drag and lift)


typedef Storage::real real;
const double eps = 1.0e-5;


int main(int argc, char ** argv)
{
	
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid_out.pmf] [Umax=2.25] [fix_cylinder=0 (0:no fix; 1:project nodes; 2: project faces)]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try
	{
		m->Load(argv[1]);
	}
	catch(...)
	{
		std::cout << "Cannot load the mesh " << argv[1] << std::endl;
		return -1;
	}
	
	std::string fout = "grid_out.pmf";
	double Umax = 2.25;
	int fix_cylinder = 0;
	if( argc > 2 ) fout = std::string(argv[2]);
	if( argc > 3 ) Umax = atof(argv[3]);
	if( argc > 4 ) fix_cylinder = atoi(argv[4]);
	
	
	
	
	double cmax[3] = {-1.0e20,-1.0e20,-1.0e20}, cmin[3] = {1.0e20,1.0e20,1.0e20};
	for(Mesh::iteratorNode n = m->BeginNode(); n != m->EndNode(); ++n)
	{
		Storage::real_array c = n->Coords();
		for(int k = 0; k < 3; ++k)
		{
			if( cmax[k] < c[k] ) cmax[k] = c[k];
			if( cmin[k] > c[k] ) cmin[k] = c[k];
		}
	}
	for(int d = 0; d < 3; ++d)
		std::cout << d << " " << cmin[d] << ":" << cmax[d] << std::endl;
			
	
	
	//TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray  bc      = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	TagRealArray  bcphi   = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,3);
	TagReal       bcp     = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagRealArray  uvw     = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3);
	TagReal       p       = m->CreateTag("P",DATA_REAL,CELL,NONE,1);
	m->self().Real(m->CreateTag("Umax",DATA_REAL,MESH,NONE,1)) = Umax;
	
	
	//this should not be needed?
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
		it->FixNormalOrientation();
	
	{ // prepare geometrical data on the mesh
		Mesh::GeomParam table;
		table[CENTROID]    = CELL | FACE; //Compute averaged center of mass
		table[NORMAL]      = FACE;        //Compute normals
		table[MEASURE]     = CELL | FACE; //Compute volumes and areas
		table[BARYCENTER]  = CELL | FACE; //Compute volumetric center of mass
		m->RemoveGeometricData(table); //Ask to precompute the data
	}
	
	if( fix_cylinder )
	{
		MarkerType cylinder = m->CreateMarker();
		
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
		{
			Storage::real  n[3], c[3];
			it->UnitNormal(n);
			it->Centroid(c);
			if( fabs(n[0]-1) < 1.0e-3 && c[0] > 2.5-eps) // outflow
			{
			}
			else if(  fabs(n[0]+1) < 1.0e-3 && c[0] < 0.0+eps) //inflow
			{
			}
			else //no-slip walls
			{
				if( fabs(n[2]) < 0.7 && c[0] > 0.45-eps && c[0] < 0.55+eps && c[1] > 0.15-eps && c[1] < 0.25+eps )
					it->SetMarker(cylinder);
			}
		}
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
			if( it->nbAdjElements(FACE,cylinder) ) it->SetMarker(cylinder);
			
		
		if( fix_cylinder != 3 )
		{
			std::cout << "project boundary nodes onto cylinder " << std::endl;
			for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) if( it->GetMarker(cylinder) )
			{
				Storage::real x = it->Coords()[0], y = it->Coords()[1], dx, dy;
				Storage::real r = sqrt((x-0.5)*(x-0.5) + (y-0.2)*(y-0.2));
				if( r )
				{
					dx = (0.05/r-1.0)*(x-0.5);
					dy = (0.05/r-1.0)*(y-0.2);
					//std::cout << "at " << x << "," << y << " r " << r << " dx " << dx << " dy " << dy << std::endl;
					it->Coords()[0] += dx;
					it->Coords()[1] += dy;
				}
				else std::cout << "node at center of cylinder: " << x << "," << y << std::endl;
			}
		}
		else fix_cylinder = 2;
		if( fix_cylinder == 2 )
		{
			TagRealArray vec_t = m->CreateTag("vec_t",DATA_REAL,FACE,FACE,2);
			std::cout << "project centers of boundary faces onto cylinder " << std::endl;
			int iter = 0;
			while(iter < 5000)
			{
				Storage::real A = 0, err = 0, fA;
				for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->GetMarker(cylinder) )
				{
					Storage::real cnt[3], dx, dy;
					it->Centroid(cnt);
					Storage::real r = sqrt((cnt[0]-0.5)*(cnt[0]-0.5) + (cnt[1]-0.2)*(cnt[1]-0.2));
					if( r )
					{
						dx = (0.05/r-1.0)*(cnt[0]-0.5);
						dy = (0.05/r-1.0)*(cnt[1]-0.2);
						//std::cout << "at " << cnt[0] << "," << cnt[1] << " r " << r << " dx " << dx << " dy " << dy << std::endl;
						vec_t[*it][0] = dx;
						vec_t[*it][1] = dy;
						fA = it->Area();
						err += sqrt(dx*dx+dy*dy)*fA;
						A += fA;
					}
					else std::cout << " face center is at center of cylinder: " << cnt[0] << "," << cnt[1] << std::endl;
				}
				err = err/A;
				if( iter % 50 == 0 )
					std::cout << "iter " << iter << " error: " << err << " area " << A << std::endl;
				if( err < 1.0e-8 ) break;
				for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) if( it->GetMarker(cylinder) )
				{
					ElementArray<Face> faces = it->getFaces(cylinder);
					Storage::real dxy[2] = {0,0}, dxyA = 0;
					for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
					{
						fA = jt->Area();
						dxy[0] += vec_t[*jt][0]*fA;
						dxy[1] += vec_t[*jt][1]*fA;
						dxyA += fA;
					}
					dxy[0] /= dxyA;
					dxy[1] /= dxyA;
					it->Coords()[0] += dxy[0];
					it->Coords()[1] += dxy[1];
				}
				iter++;
			}
			m->DeleteTag(vec_t);
		}
		m->ReleaseMarker(cylinder,FACE|NODE);
	}
	
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		Storage::real  n[3], c[3];
		it->UnitNormal(n);
		it->Centroid(c);
		bcphi[*it][0] = 1;
		bcphi[*it][1] = 0;
		bcphi[*it][2] = 0;
		if( fabs(n[0]-1) < 1.0e-3 && c[0] > 2.5-eps) // outflow
		{
			bcp[*it] = 0;
		}
		else if(  fabs(n[0]+1) < 1.0e-3 && c[0] < 0.0+eps) //inflow
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 16*Umax*c[1]*c[2]*(0.41-c[1])*(0.41-c[2])/pow(0.41,4);
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
		else //no-slip walls
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
			if( fabs(n[2]) < 0.7 && c[0] > 0.45-eps && c[0] < 0.55+eps && c[1] > 0.15-eps && c[1] < 0.25+eps )
				bcphi[*it][2] = 1;
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
