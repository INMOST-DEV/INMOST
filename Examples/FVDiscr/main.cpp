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

void make_vec(double v1[3], double v2[3], double out[3])
{
	out[0] = v1[0] - v2[0];
	out[1] = v1[1] - v2[1];
	out[2] = v1[2] - v2[2];
}

double dot_prod(double v1[3], double v2[3])
{
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double transmissibility(double vec[3], double K, double normal_face[3])
{
	double Kn[3];
	Kn[0] = K * normal_face[0], Kn[1] = K * normal_face[1], Kn[2] = K * normal_face[2];
	return dot_prod(vec,Kn);
}

double func(double x[3], double tmp)
{
//  	return x[0] + 2 * x[1] + 3 * x[2];
	return sin (M_PI * x[0]) * sin (M_PI * x[1]) * sin (M_PI * x[2]);
	(void) tmp;
}

double func_rhs(double x[3], double tmp)
{
//  	return 0;
	return -3 * tmp * M_PI * M_PI * sin (M_PI * x[0]) * sin (M_PI * x[1]) * sin (M_PI * x[2]);
}

int main(int argc,char ** argv)
{
	Solver::Initialize(&argc,&argv,""); // Initialize the solver and MPI activity
#if defined(USE_PARTITIONER)
	Partitioner::Initialize(&argc,&argv); // Initialize the partitioner activity
#endif
	if( argc > 1 )
	{
		Tag phi, tensor_K, id;
		Mesh * m = new Mesh(); // Create an empty mesh
		double ttt = Timer();
		bool repartition = false;
		m->SetCommunicator(INMOST_MPI_COMM_WORLD); // Set the MPI communicator for the mesh
		if( m->GetProcessorRank() == 0 ) // If the current process is the master one
			std::cout << argv[0] << std::endl;

		if( m->isParallelFileFormat(argv[1]) )
		{
			m->Load(argv[1]); // Load mesh from the parallel file format
			repartition = true;
		}
		else
		{
			if( m->GetProcessorRank() == 0 )
				m->Load(argv[1]); // Load mesh from the serial file format
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
		if (m->GetProcessorsNumber() > 1 && !repartition)
    { // currently only non-distributed meshes are supported by Inner_RCM partitioner
			ttt = Timer();
			Partitioner * p = new Partitioner(m);
			p->SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition); // Specify the partitioner
			p->Evaluate(); // Compute the partitioner and store new processor ID in the mesh
			delete p;
			BARRIER

			if( m->GetProcessorRank() == 0 ) std::cout << "Evaluate: " << Timer()-ttt << std::endl;

			ttt = Timer();
			m->Redistribute(); // Redistribute the mesh data
			m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
			BARRIER

			if( m->GetProcessorRank() == 0 ) std::cout << "Redistribute: " << Timer()-ttt << std::endl;
		}
#endif

		ttt = Timer();
		m->AssignGlobalID(CELL | EDGE | FACE | NODE);
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Assign id: " << Timer()-ttt << std::endl;
		id = m->GlobalIDTag(); // Get the tag of the global ID

		phi = m->CreateTag("Solution",DATA_REAL,CELL,NONE,1); // Create a new tag for the solution phi
		tensor_K = m->CreateTag("K",DATA_REAL,CELL,NONE,1); // Create a new tag for K tensor

		for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell ) // Loop over mesh cells
			if( cell->GetStatus() != Element::Ghost ) // If the cell is an own one
				cell->Real(tensor_K) = 1.0; // Store the tensor K value into the tag

		ttt = Timer();
		m->ExchangeGhost(1,FACE);
		m->ExchangeData(tensor_K,CELL,0); // Exchange the tensor_K data over processors
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;

		ttt = Timer();
		Solver S("inner_ilu2"); // Specify the linear solver to ASM+ILU2+BiCGStab one
		S.SetParameter("absolute_tolerance", "1e-8");
    Sparse::LockService L;
		Sparse::Matrix A; // Declare the matrix of the linear system to be solved
		Sparse::Vector x,b; // Declare the solution and the right-hand side vectors

		Mesh::GeomParam table;

		table[CENTROID] = CELL | FACE;
		table[NORMAL] = FACE;
		table[ORIENTATION] = FACE;
		table[MEASURE] = CELL | FACE;
		table[BARYCENTER] = CELL | FACE;
		m->PrepareGeometricData(table);
		//~ BARRIER
		//~ if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;

		unsigned idmax = 0, idmin = UINT_MAX;
		for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell )
			if( cell->GetStatus() != Element::Ghost )
			{
				unsigned pid = cell->Integer(id);
				//unsigned pid2 = cell->GlobalID();
				if( pid < idmin ) idmin = pid;
				if( pid+1 > idmax ) idmax = pid+1;
			}

		// Set the indeces intervals for the matrix and vectors
    L.SetInterval(idmin,idmax);
		A.SetInterval(idmin,idmax);
		x.SetInterval(idmin,idmax);
		b.SetInterval(idmin,idmax);
		//~ std::cout << m->GetProcessorRank() << " A,x,b interval " << idmin << ":" << idmax << " size " << idmax-idmin << std::endl;

		// Solve \nabla \cdot \nabla phi = f equation
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		//for( Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face )
    for(Storage::integer fi = 0; fi < m->FaceLastLocalID(); ++fi) if( m->isValidElement(FACE,fi) )
		{
      Face face(m,ComposeHandle(FACE,fi));
			//~ std::cout << face->LocalID() << " / " << m->NumberOfFaces() << std::endl;
			Element::Status s1,s2;
			Cell r1 = face->BackCell();
			Cell r2 = face->FrontCell();
			if( ((!r1->isValid() || (s1 = r1->GetStatus()) == Element::Ghost)?0:1) +
			    ((!r2->isValid() || (s2 = r2->GetStatus()) == Element::Ghost)?0:1) == 0) continue;
			double f_nrm[3], r1_cnt[3], r2_cnt[3], f_cnt[3], d1[3], Coef;
			double f_area = face->Area(); // Get the face area
			Storage::integer id1 = r1->Integer(id), id2;
			double K1 = r1->Real(tensor_K), K2, Kav;
			face->Normal(f_nrm); // Get the face normal
			f_nrm[0] /= f_area;
			f_nrm[1] /= f_area;
			f_nrm[2] /= f_area;
			r1->Barycenter(r1_cnt);  // Get the barycenter of the cell
			face->Barycenter(f_cnt); // Get the barycenter of the face
			if( !r2->isValid() ) // boundary condition
			{
				double bnd_pnt[3], dist;
				make_vec(f_cnt,r1_cnt,d1);
				dist = dot_prod(f_nrm,d1) / dot_prod(f_nrm,f_nrm);
				// bnd_pnt is a projection of the cell center to the face
				bnd_pnt[0] = r1_cnt[0] + dist * f_nrm[0];
				bnd_pnt[1] = r1_cnt[1] + dist * f_nrm[1];
				bnd_pnt[2] = r1_cnt[2] + dist * f_nrm[2];
				Coef = K1 * f_area / dist;
        L.Lock(id1);
				A[id1][id1] += -Coef;
				b[id1] += -Coef * func(bnd_pnt, 0);
        L.UnLock(id1);
			}
			else
			{
				K2 = r2->Real(tensor_K);
				id2 = r2->Integer(id);
				r2->Barycenter(r2_cnt);
// 				Kav = 0.5 * ( K1 + K2 ); // Arithmetic mean
				Kav = 2.0 / ( 1.0 / K1 + 1.0 / K2 ); // Harmonic mean
				make_vec(r2_cnt,r1_cnt,d1);
				Coef = transmissibility(d1,Kav,f_nrm)/dot_prod(d1,d1) * f_area;

				if( s1 != Element::Ghost )
				{
          L.Lock(id1);
					A[id1][id1] += -Coef;
					A[id1][id2] += Coef;
          L.UnLock(id2);
				}
				if( s2 != Element::Ghost )
				{
          L.Lock(id2);
					A[id2][id1] += Coef;
					A[id2][id2] += -Coef;
          L.UnLock(id2);
				}
			}
		}


		//for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell )
#if defined(USE_OMP)
#pragma omp parallel for
#endif
    for(Storage::integer ci = 0; ci < m->CellLastLocalID(); ++ci) if( m->isValidElement(CELL,ci) )
    {
      Cell cell(m,ComposeHandle(CELL,ci));
			if( cell->GetStatus() != Element::Ghost )
				b[cell->Integer(id)] += cell->Mean(func_rhs, cell->Real(tensor_K)) * cell->Volume();
    }

		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Matrix assemble: " << Timer()-ttt << std::endl;

		m->RemoveGeometricData(table); // Clean the computed geometric data

		if( argc > 3 ) // Save the matrix and RHS if required
		{
			ttt = Timer();
			A.Save(std::string(argv[2])); // "A.mtx"
			b.Save(std::string(argv[3])); // "b.rhs"
			BARRIER
			if( m->GetProcessorRank() == 0 ) std::cout << "Save matrix \"" << argv[2] << "\" and RHS \"" << argv[3] << "\": " << Timer()-ttt << std::endl;
		}

		ttt = Timer();

		S.SetMatrix(A); // Compute the preconditioner for the original matrix
		S.Solve(b,x);   // Solve the linear system with the previously computted preconditioner

		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "\nSolve system: " << Timer()-ttt << std::endl;

		ttt = Timer();

    Tag error = m->CreateTag("error",DATA_REAL,CELL,NONE,1);

		double err_C = 0.0, err_L2 = 0.0;
		//for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell )
#if defined(USE_OMP)
#pragma omp parallel for
#endif
    for(Storage::integer ci = 0; ci < m->CellLastLocalID(); ++ci) if( m->isValidElement(CELL,ci) )
    {
      Cell cell(m,ComposeHandle(CELL,ci));
			if( cell->GetStatus() != Element::Ghost )
			{
				double exact = cell->Mean(func, 0); // Compute the mean value of the function over the cell
				double err = fabs (x[cell->Integer(id)] - exact);
				if (err > err_C)
					err_C = err;
				err_L2 += err * err * cell->Volume();
				cell->Real(error) = err;
// 				x[cell->Integer(id)] = err;
			}
    }
		err_C = m->AggregateMax(err_C); // Compute the maximal C norm for the error
		err_L2 = sqrt(m->Integrate(err_L2)); // Compute the global L2 norm for the error
		if( m->GetProcessorRank() == 0 ) std::cout << "err_C  = " << err_C << std::endl;
		if( m->GetProcessorRank() == 0 ) std::cout << "err_L2 = " << err_L2 << std::endl;

		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Compute true residual: " << Timer()-ttt << std::endl;

		ttt = Timer();
		for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell )
			if( cell->GetStatus() != Element::Ghost )
				cell->Real(phi) = x[cell->Integer(id)];
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Retrieve data: " << Timer()-ttt << std::endl;

		ttt = Timer();
		m->ExchangeData(phi,CELL,0); // Data exchange over processors
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Exchange phi: " << Timer()-ttt << std::endl;

		std::string filename = "result";
		if( m->GetProcessorsNumber() == 1 )
			filename += ".vtk";
		else
			filename += ".pvtk";
		ttt = Timer();
		m->Save(filename);
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Save \"" << filename << "\": " << Timer()-ttt << std::endl;

		delete m;
	}
	else
	{
		std::cout << argv[0] << " mesh_file [A.mtx b.rhs]" << std::endl;
	}

#if defined(USE_PARTITIONER)
	Partitioner::Finalize(); // Finalize the partitioner activity
#endif
	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
