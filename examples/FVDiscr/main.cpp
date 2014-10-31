//g++ main.cpp -lpetsc -L/usr/X11R6/lib -lX11 -ldmumps -lmumps_common -lmpi_f77 -lscalapack -lpord -lblacs -lparmetis -lmetis -lmpi -lHYPRE -lgfortran -lblas -llapack ../../INMOST.a
//mpicxx main.cpp -lpetsc -L/usr/X11R6/lib -lX11 -ldmumps -lmumps_common -lmpi_f77 -lscalapack -lpord -lblacs -lparmetis -lmetis -lmpi -lHYPRE -lgfortran -lblas -llapack ../../INMOST.a -lzoltan -lparmetis -lmetis
// run ./a.out grids/rezultMesh.vtk MATERIALS -ksp_monitor -ksp_view -mat_type aij -vec_type standard

#include "../../inmost.h"
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.141592653589
#endif
using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

void make_vec(Storage::real v1[3], Storage::real v2[3], Storage::real out[3])
{
	out[0] = v1[0] - v2[0];
	out[1] = v1[1] - v2[1];
	out[2] = v1[2] - v2[2];
}

Storage::real dot_prod(Storage::real v1[3], Storage::real v2[3])
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

Storage::real transmissibility(Storage::real vec[3], Storage::real K, Storage::real normal_face[3])
{
	Storage::real Kn[3];
	Kn[0] = K * normal_face[0], Kn[1] = K * normal_face[1], Kn[2] = K * normal_face[2];
	return dot_prod(vec,Kn);
}

Storage::real func(Storage::real x[3], Storage::real tmp)
{
//  	return x[0] + 2 * x[1] + 3 * x[2];
	return sin (M_PI * x[0]) * sin (M_PI * x[1]) * sin (M_PI * x[2]);
	(void) tmp;
}

Storage::real func_rhs(Storage::real x[3], Storage::real tmp)
{
//  	return 0;
	return -3 * tmp * M_PI * M_PI * sin (M_PI * x[0]) * sin (M_PI * x[1]) * sin (M_PI * x[2]);
}

int main(int argc,char ** argv)
{
	Solver::Initialize(&argc,&argv,"");
#if defined(USE_PARTITIONER)
	Partitioner::Initialize(&argc,&argv);
#endif
	if( argc > 1 )
	{
		Tag phi, tensor_K, id;
		Mesh * m = new Mesh();
		double ttt = Timer();
		bool repartition = false;
		m->SetCommunicator(INMOST_MPI_COMM_WORLD);
		if( m->GetProcessorRank() == 0 ) std::cout << argv[0] << std::endl;
		
		if( m->isParallelFileFormat(argv[1]) ) 
		{
			m->Load(argv[1]);
			repartition = true;
		}
		else
		{
			if( m->GetProcessorRank() == 0 ) 
				m->Load(argv[1]);
		}
		BARRIER
		
		if( m->GetProcessorRank() == 0 ) std::cout << "Processors: " << m->GetProcessorsNumber() << std::endl;
		if( m->GetProcessorRank() == 0 ) std::cout << "Load(MPI_File): " << Timer()-ttt << std::endl;
		
		//~ double ttt2 = Timer();
		//~ Mesh t;
		//~ t.SetCommunicator(INMOST_MPI_COMM_WORLD);
		//~ t.SetParallelFileStrategy(0);
		//~ t.Load(argv[1]);
		//~ BARRIER
		//~ if( m->GetProcessorRank() == 0 ) std::cout << "Load(MPI_Scatter): " << Timer()-ttt2 << std::endl;
		
#if defined(USE_PARTITIONER)
		ttt = Timer();
		Partitioner * p = new Partitioner(m);
		p->SetMethod(Partitioner::Inner_RCM,Partitioner::Partition);
		p->Evaluate();
		delete p;
		BARRIER
	
		if( m->GetProcessorRank() == 0 ) std::cout << "Evaluate: " << Timer()-ttt << std::endl;
		
		ttt = Timer();
		m->Redistribute();
		m->ReorderEmpty(CELL|FACE|EDGE|NODE);
		BARRIER
	
		if( m->GetProcessorRank() == 0 ) std::cout << "Redistribute: " << Timer()-ttt << std::endl;
#endif
		
		ttt = Timer();
		m->AssignGlobalID(CELL | EDGE | FACE | NODE);
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Assign id: " << Timer()-ttt << std::endl;		
		id = m->GlobalIDTag();
		
		phi = m->CreateTag("Solution",DATA_REAL,CELL,NONE,1);
		tensor_K = m->CreateTag("K",DATA_REAL,CELL,NONE,1);
		
		for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
		if( cell->GetStatus() != Element::Ghost )
			cell->Real(tensor_K) = 1.0;
		
		ttt = Timer();
		m->ExchangeGhost(1,FACE);
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;
		
		ttt = Timer();
		Solver S(Solver::INNER_ILU2);
		Solver::Matrix A;
		Solver::Vector x,b;
		
		std::map<GeometricData,ElementType> table;
		
		table[MEASURE] = CELL | FACE;
		table[CENTROID] = CELL | FACE;
		table[BARYCENTER] = CELL | FACE;
		table[NORMAL] = FACE;
		table[ORIENTATION] = FACE;
		m->PrepareGeometricData(table);
		//~ BARRIER
		//~ if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
		
		unsigned idmax = 0, idmin = UINT_MAX;
		for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
		if( cell->GetStatus() != Element::Ghost )
		{
			unsigned pid = cell->Integer(id);
			unsigned pid2 = cell->GlobalID();
			if( pid < idmin ) idmin = pid;
			if( pid+1 > idmax ) idmax = pid+1;
		}
		
		A.SetInterval(idmin,idmax);
		x.SetInterval(idmin,idmax);
		b.SetInterval(idmin,idmax);
		//~ std::cout << m->GetProcessorRank() << " A,x,b interval " << idmin << ":" << idmax << " size " << idmax-idmin << std::endl;
		
		//solve \nabla \cdot \nabla phi = f equation
		for(Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face)
		{
			//~ std::cout << face->LocalID() << " / " << m->NumberOfFaces() << std::endl;
			Element::Status s1,s2;
			Cell * r1 = face->BackCell();
			Cell * r2 = face->FrontCell();
			if( ((r1 == NULL || (s1 = r1->GetStatus()) == Element::Ghost)?0:1)+ 
			    ((r2 == NULL || (s2 = r2->GetStatus()) == Element::Ghost)?0:1) == 0) continue;
			Storage::real f_nrm[3], r1_cnt[3], r2_cnt[3], f_cnt[3], d1[3], Coef;
			Storage::real f_area = face->Area();
			Storage::real vol1 = r1->Volume(), vol2;
			Storage::integer id1 = r1->Integer(id), id2;
			Storage::real K1 = r1->Real(tensor_K), K2, Kav;
			face->Normal(f_nrm);
			f_nrm[0] /= f_area;
			f_nrm[1] /= f_area;
			f_nrm[2] /= f_area;
			r1->Barycenter(r1_cnt);
			face->Barycenter(f_cnt);
			if( r2 == NULL ) //boundary condition
			{
				Storage::real bnd_pnt[3], dist;
				make_vec(f_cnt,r1_cnt,d1);
				dist = dot_prod(f_nrm,d1) / dot_prod(f_nrm,f_nrm);
				// bnd_pnt is a projection of the cell center to the face
				bnd_pnt[0] = r1_cnt[0] + dist * f_nrm[0], 
				bnd_pnt[1] = r1_cnt[1] + dist * f_nrm[1], 
				bnd_pnt[2] = r1_cnt[2] + dist * f_nrm[2];
				Coef = K1 * f_area / dist;
				A[id1][id1] += -Coef;
				b[id1] += -Coef * func(bnd_pnt, 0);
			}
			else
			{
				vol2 = r2->Volume();
				K2 = r2->Real(tensor_K);
				id2 = r2->Integer(id);
				r2->Barycenter(r2_cnt);
// 				Kav = 0.5 * ( K1 + K2 ); // Arithmetic mean
				Kav = 2.0 / ( 1.0 / K1 + 1.0 / K2 ); // Harmonic mean
				make_vec(r2_cnt,r1_cnt,d1);
				Coef = transmissibility(d1,Kav,f_nrm)/dot_prod(d1,d1) * f_area;
				
				if( s1 != Element::Ghost )
				{
					A[id1][id1] += -Coef;
					A[id1][id2] += Coef;
				}
				if( s2 != Element::Ghost )
				{
					A[id2][id1] += Coef;
					A[id2][id2] += -Coef;
				}
			}
		}
		for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
		if( cell->GetStatus() != Element::Ghost )
			b[cell->Integer(id)] += cell->Mean(func_rhs, cell->Real(tensor_K)) * cell->Volume();
			
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Matrix assemble: " << Timer()-ttt << std::endl;
		
		m->RemoveGeometricData(table);
		
		//~ ttt = Timer();
		//~ A.Save("A.mtx");
		//~ b.Save("b.rhs");
		//~ BARRIER
		//~ if( m->GetProcessorRank() == 0 ) std::cout << "Save matrix and RHS: " << Timer()-ttt << std::endl;
		
		ttt = Timer();
		
		S.SetMatrix(A);
		S.Solve(b,x);
		
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Solve system: " << Timer()-ttt << std::endl;

		ttt = Timer();
		
		Storage::real err_C = 0.0, err_L2 = 0.0;
		for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); cell++)
		if( cell->GetStatus() != Element::Ghost )
		{
			Storage::real exact = cell->Mean(func, 0);
			Storage::real err = fabs (x[cell->Integer(id)] - exact);
			if (err > err_C)
				err_C = err;
			err_L2 += err * err * cell->Volume();
// 			x[cell->Integer(id)] = err;
		}
		err_C = m->AggregateMax(err_C);
		err_L2 = sqrt(m->Integrate(err_L2));
		if( m->GetProcessorRank() == 0 ) std::cout << "err_C = " << err_C << std::endl;
		if( m->GetProcessorRank() == 0 ) std::cout << "err_L2 = " << err_L2 << std::endl;
		
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Compute residual: " << Timer()-ttt << std::endl;

		ttt = Timer();
		
		for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); cell++)
		if( cell->GetStatus() != Element::Ghost )
			cell->Real(phi) = x[cell->Integer(id)];
		
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Retrive data: " << Timer()-ttt << std::endl;
		
		ttt = Timer();
		m->ExchangeData(phi,CELL);
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Exchange phi: " << Timer()-ttt << std::endl;
		
		ttt = Timer();
		if (m->GetProcessorsNumber() == 1 )
			m->Save("grid.vtk");
		else
			m->Save("grid.pvtk");
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Save pvtk: " << Timer()-ttt << std::endl;
		
		delete m;
	}
	else
	{
		std::cout << argv[0] << " [mesh_file]" << std::endl;
	}
#if defined(USE_PARTITIONER)
	Partitioner::Finalize();
#endif
	Solver::Finalize();
	return 0;
}
