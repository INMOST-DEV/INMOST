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
		TagReal phi;
		TagReal tag_F;
		TagRealArray tag_K;
		TagRealArray tag_BC;
		TagReal phi_ref;
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
		BARRIER;
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
		if (m->GetProcessorsNumber() > 1)
		{ // currently only non-distributed meshes are supported by Inner_RCM partitioner
			ttt = Timer();
			Partitioner * p = new Partitioner(m);
			p->SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition); // Specify the partitioner
			p->Evaluate(); // Compute the partitioner and store new processor ID in the mesh
			delete p;
			BARRIER;

			if( m->GetProcessorRank() == 0 ) std::cout << "Evaluate: " << Timer()-ttt << std::endl;

			ttt = Timer();
			m->Redistribute(); // Redistribute the mesh data
			m->ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
			BARRIER;

			if( m->GetProcessorRank() == 0 ) std::cout << "Redistribute: " << Timer()-ttt << std::endl;
		}
#endif

		ttt = Timer();
		phi = m->CreateTag("Solution",DATA_REAL,CELL,NONE,1); // Create a new tag for the solution phi
		
		bool makerefsol = true;
		
		if( m->HaveTag("PERM" ) )
		{
			tag_K = m->GetTag("PERM");
			makerefsol = false;
		}
		else
		{
			std::cout << "Set perm" << std::endl;
			tag_K = m->CreateTag("PERM",DATA_REAL,CELL,NONE,1); // Create a new tag for K tensor
			for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell ) // Loop over mesh cells
				tag_K[*cell][0] = 1.0; // Store the tensor K value into the tag
		}
		
		if( m->HaveTag("FORCE") )
		{
			tag_F = m->GetTag("FORCE");
			makerefsol = false;
		}
		else
		{
			std::cout << "Set rhs" << std::endl;
			tag_F = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,1); // Create a new tag for external force
			double x[3];
			for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell ) // Loop over mesh cells
			{
				cell->Centroid(x);
				tag_F[*cell] = func_rhs(x,1);  //cell->Mean(func_rhs, 1);
			}
		}
		
		if( m->HaveTag("BOUNDARY_CONDITION") )
		{
			tag_BC = m->GetTag("BOUNDARY_CONDITION");
			makerefsol = false;
		}
		else
		{
			std::cout << "Set boundary conditions" << std::endl;
			double x[3];
			tag_BC = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,NONE,3);
			for( Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face )
				if( face->Boundary() && !(face->GetStatus() == Element::Ghost) )
				{
					face->Centroid(x);
					tag_BC[*face][0] = 1; //dirichlet
					tag_BC[*face][1] = 0; //neumann
					tag_BC[*face][2] = func(x,0);//face->Mean(func, 0);
				}
		}
		
		if(m->HaveTag("REFERENCE_SOLUTION") )
			phi_ref = m->GetTag("REFERENCE_SOLUTION");
		else if( makerefsol )
		{
			phi_ref = m->CreateTag("REFRENCE_SOLUTION",DATA_REAL,CELL,NONE,1);
			double x[3];
			for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell )
			{
				cell->Centroid(x);
				phi_ref[*cell] = func(x,0);//cell->Mean(func, 0);
			}
		}

		ttt = Timer();
		m->ExchangeGhost(1,FACE);
		m->ExchangeData(tag_K,CELL,0); // Exchange the tag_K data over processors
		BARRIER;
		if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;



		ttt = Timer();
		Solver S("inner_ilu2"); // Specify the linear solver to ASM+ILU2+BiCGStab one
		S.SetParameter("absolute_tolerance", "1e-8");
		S.SetParameter("schwartz_overlap", "2");
    	Residual R; // Residual vector
    	Sparse::LockService Locks;
		Sparse::Vector Update; // Declare the solution and the right-hand side vectors

		{
			Mesh::GeomParam table;
			table[CENTROID] = CELL | FACE;
			table[NORMAL] = FACE;
			table[ORIENTATION] = FACE;
			table[MEASURE] = CELL | FACE;
			table[BARYCENTER] = CELL | FACE;
			m->PrepareGeometricData(table);
		}
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;
		
		{
			Automatizator aut;
			Automatizator::MakeCurrent(&aut);
			INMOST_DATA_ENUM_TYPE iphi = aut.RegisterTag(phi,CELL);
			aut.EnumerateEntries();

			// Set the indeces intervals for the matrix and vectors
			R.SetInterval(aut.GetFirstIndex(),aut.GetLastIndex());
			R.InitLocks();
			Update.SetInterval(aut.GetFirstIndex(),aut.GetLastIndex());
			
			dynamic_variable Phi(aut,iphi);
			// Solve \nabla \cdot \nabla phi = f equation
			//for( Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face )
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				variable flux; //should be more efficient to define here to avoid multiple memory allocations if storage for variations should be expanded
				rMatrix x1(3,1), x2(3,1), xf(3,1), n(3,1);
				double d1, d2, k1, k2, area, T, a, b, c;
#if defined(USE_OMP)
#pragma omp for
#endif
				for(Storage::integer iface = 0; iface < m->FaceLastLocalID(); ++iface ) if( m->isValidFace(iface) )
				{
					Face face = Face(m,ComposeFaceHandle(iface));
					Element::Status s1,s2;
					Cell r1 = face->BackCell();
					Cell r2 = face->FrontCell();
					if( ((!r1->isValid() || (s1 = r1->GetStatus()) == Element::Ghost)?0:1) +
						((!r2->isValid() || (s2 = r2->GetStatus()) == Element::Ghost)?0:1) == 0) continue;
					
					
					area = face->Area(); // Get the face area
					face->UnitNormal(n.data()); // Get the face normal
					face->Centroid(xf.data()); // Get the barycenter of the face
					r1->Centroid(x1.data());  // Get the barycenter of the cell
					k1 = n.DotProduct(rMatrix::FromTensor(tag_K[r1].data(),
														  tag_K[r1].size(),3)*n);
					d1 = fabs(n.DotProduct(xf-x1));
					if( !r2->isValid() ) // boundary condition
					{
						// bnd_pnt is a projection of the cell center to the face
						// a*pb + bT(pb-p1) = c
						// F = T(pb-p1)
						// pb = (c + bTp1)/(a+bT)
						// F = T/(a+bT)(c - ap1)
						T = k1/d1;
						a = 0;
						b = 1;
						c = 0;
						if( tag_BC.isValid() )
						{
							a = tag_BC[face][0];
							b = tag_BC[face][1];
							c = tag_BC[face][2];
						}
						R.Lock(Phi.Index(r1));
						R[Phi.Index(r1)] -=  T/(a + b*T) * area * (c - a*Phi(r1));
						R.UnLock(Phi.Index(r1));
					}
					else
					{
						r2->Centroid(x2.data());
						k2 = n.DotProduct(rMatrix::FromTensor(tag_K[r2].data(),
															  tag_K[r2].size(),3)*n);
						d2 = fabs(n.DotProduct(x2-xf));
						T = 1.0/(d1/k1 + d2/k2);
						flux = T* area * (Phi(r2) - Phi(r1));
						if( s1 != Element::Ghost )
						{
							R.Lock(Phi.Index(r1));
							R[Phi.Index(r1)] -= flux;
							R.UnLock(Phi.Index(r1));
						}
						if( s2 != Element::Ghost )
						{
							R.Lock(Phi.Index(r2));
							R[Phi.Index(r2)] += flux;
							R.UnLock(Phi.Index(r2));
						}
					}
				}
			}
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for( Storage::integer icell = 0; icell < m->CellLastLocalID(); ++icell ) if( m->isValidCell(icell) )
			{
				Cell cell = Cell(m,ComposeCellHandle(icell));
				if( cell->GetStatus() != Element::Ghost )
					R[Phi.Index(cell)] += tag_F[cell] * cell->Volume();
			}
			BARRIER;
			if( m->GetProcessorRank() == 0 ) std::cout << "Matrix assemble: " << Timer()-ttt << std::endl;

			//m->RemoveGeometricData(table); // Clean the computed geometric data

			if( argc > 3 ) // Save the matrix and RHS if required
			{
				ttt = Timer();
				R.GetJacobian().Save(std::string(argv[2])); // "A.mtx"
				R.GetResidual().Save(std::string(argv[3])); // "b.rhs"
				BARRIER;
				if( m->GetProcessorRank() == 0 ) std::cout << "Save matrix \"" << argv[2] << "\" and RHS \"" << argv[3] << "\": " << Timer()-ttt << std::endl;
			}

			ttt = Timer();

			S.SetMatrix(R.GetJacobian()); // Compute the preconditioner for the original matrix
			S.Solve(R.GetResidual(),Update);   // Solve the linear system with the previously computted preconditioner

			BARRIER;
			if( m->GetProcessorRank() == 0 ) 
			{
				std::cout << S.Residual() << " " << S.Iterations() << " " << S.ReturnReason() << std::endl;
				std::cout << "Solve system: " << Timer()-ttt << std::endl;
			}

			ttt = Timer();

			
			if( phi_ref.isValid() )
			{
				Tag error = m->CreateTag("error",DATA_REAL,CELL,NONE,1);
				double err_C = 0.0, err_L2 = 0.0;
#if defined(USE_OMP)
#pragma omp parallel
#endif
				{
					double local_err_C = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:err_L2)
#endif
					for( Storage::integer icell = 0; icell < m->CellLastLocalID(); ++icell )
					{
						Cell cell = Cell(m,ComposeCellHandle(icell));
						if( cell->GetStatus() != Element::Ghost )
						{
							double old = phi[cell];
							double exact = phi_ref[cell];
							double res = Update[Phi.Index(cell)];
							double sol = old-res;
							double err = fabs (sol - exact);
							if (err > local_err_C) local_err_C = err;
							err_L2 += err * err * cell->Volume();
							cell->Real(error) = err;
							phi[cell] = sol;
						}
					}
#if defined(USE_OMP)
#pragma omp critical
#endif
					{
						if( local_err_C > err_C ) err_C = local_err_C;
					}
				}
				err_C = m->AggregateMax(err_C); // Compute the maximal C norm for the error
				err_L2 = sqrt(m->Integrate(err_L2)); // Compute the global L2 norm for the error
				if( m->GetProcessorRank() == 0 ) std::cout << "err_C  = " << err_C << std::endl;
				if( m->GetProcessorRank() == 0 ) std::cout << "err_L2 = " << err_L2 << std::endl;
			}
		}
		BARRIER;
		if( m->GetProcessorRank() == 0 ) std::cout << "Compute true residual: " << Timer()-ttt << std::endl;

		ttt = Timer();
		m->ExchangeData(phi,CELL,0); // Data exchange over processors
		BARRIER;
		if( m->GetProcessorRank() == 0 ) std::cout << "Exchange phi: " << Timer()-ttt << std::endl;

		std::string filename = "result";
		if( m->GetProcessorsNumber() == 1 )
			filename += ".vtk";
		else
			filename += ".pvtk";
		ttt = Timer();
		m->Save(filename);
		m->Save("result.pmf");
		BARRIER;
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
