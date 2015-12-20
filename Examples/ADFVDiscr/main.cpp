#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "inmost.h"
using namespace INMOST;

#ifndef M_PI
#define M_PI 3.141592653589
#endif

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


Storage::real func(Storage::real x[3], Storage::real tmp)
{
//  	return x[0] + 2 * x[1] + 3 * x[2];
  double s0 = sin (M_PI * x[0]); 
  double s1 = sin (M_PI * x[1]);
  double s2 = sin (M_PI * x[2]);
	return s0 * s1 * s2;
	(void) tmp;
}

Storage::real func_rhs(Storage::real x[3], Storage::real tmp)
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
		if (!repartition) 
    { // currently only non-distributed meshes are supported by Inner_RCM partitioner
			ttt = Timer();
			Partitioner * p = new Partitioner(m);
			p->SetMethod(Partitioner::Inner_RCM,Partitioner::Partition); // Specify the partitioner
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
		//m->Save("solution_check_0.vtk");
		phi = m->CreateTag("Solution",DATA_REAL,CELL,NONE,1); // Create a new tag for the solution phi
		tensor_K = m->CreateTag("K",DATA_REAL,CELL,NONE,1); // Create a new tag for K tensor
    //m->Save("solution_check_1.vtk");

		
		for( Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell ) // Loop over mesh cells
			if( cell->GetStatus() != Element::Ghost ) // If the cell is an own one
				cell->Real(tensor_K) = 1.0; // Store the tensor K value into the tag
		
		ttt = Timer();
		m->ExchangeGhost(1,FACE);
		m->ExchangeData(tensor_K,CELL,0); // Exchange the tensor_K data over processors
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;

    
		
		ttt = Timer();
		Solver S(Solver::INNER_ILU2); // Specify the linear solver to ASM+ILU2+BiCGStab one
		S.SetParameterReal("absolute_tolerance",1e-8);
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
		
    {
		  Automatizator aut(m);
      Automatizator::MakeCurrent(&aut);
      INMOST_DATA_ENUM_TYPE iphi = aut.RegisterDynamicTag(phi,CELL);
      aut.EnumerateDynamicTags();

		  // Set the indeces intervals for the matrix and vectors
      A.SetInterval(aut.GetFirstIndex(),aut.GetLastIndex());
      x.SetInterval(aut.GetFirstIndex(),aut.GetLastIndex());
      b.SetInterval(aut.GetFirstIndex(),aut.GetLastIndex());
		  //~ std::cout << m->GetProcessorRank() << " A,x,b interval " << idmin << ":" << idmax << " size " << idmax-idmin << std::endl;
      dynamic_variable Phi(aut,iphi);
		  // Solve \nabla \cdot \nabla phi = f equation
		  //for( Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face )
#if defined(USE_OMP)
#pragma omp parallel
#endif
      {
        variable flux; //should be more efficient to define here to avoid multiple memory allocations if storage for variations should be expanded
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
          Storage::integer i1 = aut.GetDynamicIndex(r1,iphi), i2;
			    Storage::real f_nrm[3], r1_cnt[3], r2_cnt[3], f_cnt[3], d1, d2, D, v[3], T;
			    Storage::real f_area = face->Area(); // Get the face area
			    face->UnitNormal(f_nrm); // Get the face normal
			    r1->Centroid(r1_cnt);  // Get the barycenter of the cell
			    face->Centroid(f_cnt); // Get the barycenter of the face
          Sparse::RowMerger & r = aut.GetMerger();
			    if( !r2->isValid() ) // boundary condition
			    {
				    Storage::real bnd_pnt[3], dist;
				    make_vec(f_cnt,r1_cnt,v);
				    dist = dot_prod(f_nrm,v);
				    // bnd_pnt is a projection of the cell center to the face
				    bnd_pnt[0] = r1_cnt[0] + dist * f_nrm[0];
				    bnd_pnt[1] = r1_cnt[1] + dist * f_nrm[1];
				    bnd_pnt[2] = r1_cnt[2] + dist * f_nrm[2];
            T = r1->Real(tensor_K) * f_area / dist;
            //flux =  T * (func(bnd_pnt,0) - variable(aut,r1,iphi));
            flux = T * (func(bnd_pnt,0) - Phi(r1));
            A[i1].Lock();
            r.PushRow(1.0,A[i1]);
            r.AddRow(-1.0,flux.GetRow());
            r.RetriveRow(A[i1]);
            r.Clear();
            b[i1] -= flux.GetValue();
            A[i1].Unlock();
			    }
			    else
			    {
            i2 = aut.GetDynamicIndex(r2,iphi);
				    r2->Centroid(r2_cnt);
				    D = dot_prod(f_nrm,f_cnt);
            d1 = fabs(dot_prod(r1_cnt,f_nrm) - D);
            d2 = fabs(dot_prod(r2_cnt,f_nrm) - D);
            T = 1.0 / (d1/r1->Real(tensor_K) + d2/r2->Real(tensor_K)) * f_area;
				    //flux = T * (variable(aut,r2,iphi) - variable(aut,r1,iphi));//(unknown(aut,r2,iphi) - unknown(aut,r1,iphi));
            flux = T * (Phi(r2) - Phi(r1));
				    if( s1 != Element::Ghost )
				    {
              A[i1].Lock();
              r.PushRow(1.0,A[i1]);
              r.AddRow(-1.0,flux.GetRow());
              r.RetriveRow(A[i1]);
              r.Clear();
              b[i1] -= flux.GetValue();
              A[i1].Unlock();
				    }
				    if( s2 != Element::Ghost )
				    {
              A[i2].Lock();
              r.PushRow(1.0,A[i2]);
              r.AddRow(1.0,flux.GetRow());
              r.RetriveRow(A[i2]);
              r.Clear();
              b[i2] += flux.GetValue();
              A[i2].Unlock();
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
		  if( m->GetProcessorRank() == 0 ) 
      {
        std::cout << S.Residual() << " " << S.Iterations() << " " << S.GetReason() << std::endl;
        std::cout << "Solve system: " << Timer()-ttt << std::endl;
      }

		  ttt = Timer();

      Tag error = m->CreateTag("error",DATA_REAL,CELL,NONE,1);

		  Storage::real err_C = 0.0, err_L2 = 0.0;
#if defined(USE_OMP)
#pragma omp parallel
#endif
      {
        Storage::real local_err_C = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:err_L2)
#endif
		    for( Storage::integer icell = 0; icell < m->CellLastLocalID(); ++icell )
        {
          Cell cell = Cell(m,ComposeCellHandle(icell));
			    if( cell->GetStatus() != Element::Ghost )
			    {
            Storage::real old = cell->Real(phi);
				    Storage::real exact = cell->Mean(func, 0); // Compute the mean value of the function over the cell
            Storage::real res = x[aut.GetDynamicIndex(cell->self(),iphi)];
            Storage::real sol = old-res;
            Storage::real err = fabs (sol - exact);
  		      if (err > local_err_C) local_err_C = err;
            err_L2 += err * err * cell->Volume();
            cell->Real(error) = err;
            cell->Real(phi) = sol;
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
		BARRIER
		if( m->GetProcessorRank() == 0 ) std::cout << "Compute true residual: " << Timer()-ttt << std::endl;

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
    m->Save("result.pmf");
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
