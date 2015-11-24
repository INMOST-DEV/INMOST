#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>

#include "inmost.h"
using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

int main(int argc, char ** argv)
{
	if( argc < 3 || argc > 1 && ( atoi(argv[1]) < 0 || atoi(argv[1]) > 11 ) )
	{
		std::cout << "Usage: " << argv[0] << " method_number<0:INNER_ILU2,1:INNER_DDPQILUC,2:INNER_MPTILUC,3:INNER_MPTILU2,4:Trilinos_Aztec,5:Trilinos_Belos,6:Trilinos_ML,7:Trilinos_Ifpack,8:PETSc,9:ANI,10:FCBIILU2,11:K3BIILU2> matrix.mtx [right_hand_side.rhs] [exact_solution] [solver_options.txt]" << std::endl;
		std::cout << "Example: " << argv[0] << "  0 a.mtx b.rhs" << std::endl;
		std::cout << "or just: " << argv[0] << " 11 a.mtx - 1 database.txt" << std::endl;
		return -1;
	}
	Solver::Type type;
	switch(atoi(argv[1])) // Actually: type = atoi(argv[1])
	{
		case  0: type = Solver::INNER_ILU2;      break;
		case  1: type = Solver::INNER_DDPQILUC;  break;
		case  2: type = Solver::INNER_MPTILUC;   break;
		case  3: type = Solver::INNER_MPTILU2;   break;
		case  4: type = Solver::Trilinos_Aztec;  break;
		case  5: type = Solver::Trilinos_Belos;  break;
		case  6: type = Solver::Trilinos_ML;     break;
		case  7: type = Solver::Trilinos_Ifpack; break;
		case  8: type = Solver::PETSc;           break;
		case  9: type = Solver::ANI;             break;
		case 10: type = Solver::FCBIILU2;        break;
		case 11: type = Solver::K3BIILU2;        break;
	}
	Solver::Initialize(&argc,&argv,argc > 5 ? argv[5] : NULL); // Initialize the linear solver in accordance with args
	{
		int rank, procs;
#if defined(USE_MPI)
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);  // Get the rank of the current process
		MPI_Comm_size(MPI_COMM_WORLD,&procs); // Get the total number of processors used
#else
		rank = 0;
		procs = 1;
#endif
    if( rank == 0 )
      std::cout << "Solving with " << Solver::TypeName(type) << std::endl;
		//std::cout << rank << "/" << procs << " " << argv[0] << std::endl;
		Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
		Sparse::Vector b("rhs"); // Declare the right-hand side vector
		Sparse::Vector x("sol"); // Declare the solution vector
		//std::cout << rank << " load matrix from " << std::string(argv[2]) << " ..." << std::endl;
		double t = Timer(), tt = Timer(), tsol = 0;
		mat.Load(std::string(argv[2])); //if interval parameters not set, matrix will be divided automatically
		BARRIER
		if( !rank ) std::cout << "load matrix:    " << Timer() - t << std::endl;
		//mat.Save("test.mtx");
		t = Timer();
		if( argc > 4 && std::string(argv[3]) == std::string("-")
		             && std::string(argv[4]) == std::string("1") )
		{
			if( !rank ) std::cout << "[" << rank << "]: set RHS = A * [1]" << std::endl;
			Solver::Vector ex("exact");  // Declare the temporary exact solution vector
			INMOST_DATA_ENUM_TYPE mbeg,mend,k;
			mat.GetInterval(mbeg,mend);
			ex.SetInterval(mbeg,mend);
			for( k = mbeg; k < mend; ++k ) ex[k] = 1.0;
			Solver::Matrix a = mat; // Declare a temporary copy of the original matrix 
			Solver::OrderInfo info;
			info.PrepareMatrix(a,0);
			info.PrepareVector(ex);
			info.Update(ex);
			a.MatVec(1.0,ex,0.0,b); // Multiply the original matrix by a vector
			//b.Save("b.rhs");      // Save the computed RHS if required
		}
		else if( argc > 3 && std::string(argv[3]) != std::string("1") )
		{
			if( !rank ) std::cout << "[" << rank << "]: load vector from " << std::string(argv[3]) << std::endl;
			b.Load(std::string(argv[3])); // Load RHS vector
		}
		else // Set local RHS to 1 if it was not specified
		{
			if( !rank ) std::cout << "[" << rank << "]: set RHS=1" << std::endl;
			INMOST_DATA_ENUM_TYPE mbeg,mend,k;
			mat.GetInterval(mbeg,mend);
			b.SetInterval(mbeg,mend);
			for( k = mbeg; k < mend; ++k ) b[k] = 1.0;
		}
		BARRIER
		if( !rank ) std::cout << "load vector:    " << Timer() - t << std::endl;
		tt = Timer();
		bool success = false;
		int iters;
		double resid, realresid = 0;
		std::string reason;
		{
			Solver s(type); // Declare the linear solver by specified type

			s.SetParameterEnum("gmres_substeps",4);
			s.SetParameterReal("relative_tolerance",1.0e-6);
			s.SetParameterReal("absolute_tolerance",1.0e-16);

			s.SetParameterEnum("reorder_nonzeros",0);
			s.SetParameterEnum("rescale_iterations",8);
			s.SetParameterEnum("adapt_ddpq_tolerance",0);
			
			s.SetParameterReal("drop_tolerance",3.0e-3);
			s.SetParameterReal("reuse_tolerance",1.0e-5);
			s.SetParameterReal("ddpq_tolerance",0.0);

			s.SetParameterEnum("condition_estimation",1);
      s.SetParameterEnum("schwartz_overlap",3);
			

			t = Timer();
			s.SetMatrix(mat); // Compute the preconditioner for the original matrix
			BARRIER
      tsol += Timer()-t;
			if( !rank ) std::cout << "preconditioner: " << Timer() - t << std::endl;
			t = Timer();
			success = s.Solve(b,x); // Solve the linear system with the previously computted preconditioner
			BARRIER
			if( !rank ) std::cout << "iterations:     " << Timer() - t << std::endl;
			iters  = s.Iterations(); // Get the number of iterations performed
			resid  = s.Residual();   // Get the final residual achieved
			reason  = s.GetReason(); // Get the convergence reason
			//x.Save("output.sol");  // Save the solution if required
		}
		tt = Timer() - tt;

		{ // Compute the true residual
			double aresid = 0, bresid = 0;
			Sparse::Vector test;
			t = Timer();
			Solver::OrderInfo info;
			info.PrepareMatrix(mat,0);
			info.PrepareVector(x);
			info.Update(x);

			mat.MatVec(1.0,x,0.0,test); // Multiply the original matrix by a vector

			{
				INMOST_DATA_ENUM_TYPE mbeg,mend,k;
				info.GetLocalRegion(info.GetRank(),mbeg,mend);
				for( k = mbeg; k < mend; ++k )
				{
					aresid += (test[k]-b[k])*(test[k]-b[k]);
					bresid += b[k]*b[k];
				}
			}
			double temp[2] = {aresid,bresid}, recv[2] = {aresid,bresid};
#if defined(USE_MPI)
			MPI_Reduce(temp,recv,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			if( info.GetRank() == 0 ) std::cout << "||Ax-b||=" << sqrt(recv[0]) << " ||b||=" << sqrt(recv[1]) << " ||Ax-b||/||b||=" << sqrt(recv[0]/(recv[1]+1.0e-100)) << std::endl;
#endif
			realresid = sqrt(recv[0]/(recv[1]+1.0e-100));
			//realresid = sqrt(realresid);

			info.RestoreVector(x);
			if( !rank ) std::cout << "norms: " << Timer() - t << std::endl;
		}

		{
			if( !rank )
			{
				std::cout << procs << " processors for Solver::type=" << type;
				if( success )
				{
					std::cout << " solved in " << tt << " secs (without reading " << tsol << " secs)";
					std::cout << " with " << iters << " iterations to " << resid << " norm";
				}
				else std::cout << " failed to solve";
				std::cout  << " matrix \"" << argv[2] << "\"";
				if( argc > 3 ) std::cout  << " vector \"" << argv[3] << "\"";
				std::cout << " true residual ||Ax-b||/||b||=" << realresid;
				std::cout << std::endl;
				std::cout << "reason: " << reason << std::endl;
			}
		}

		{ // Compare solution with exact one if available
			if( argc > 4 )
			{
				Solver::Vector ex("exact");  // Declare the exact solution vector
				Solver::Vector err("error"); // Declare the solution error vector
				INMOST_DATA_ENUM_TYPE mbeg,mend,k;
				mat.GetInterval(mbeg,mend);
				err.SetInterval(mbeg,mend);
				if( std::string(argv[4]) == std::string("1"))
				{
					ex.SetInterval(mbeg,mend);
					for( k = mbeg; k < mend; ++k ) ex[k] = 1.0;
				}
				else
				{
					//std::cout << "[" << rank << "]: load exact solution vector from " << std::string(argv[4]) << std::endl;
					ex.Load(std::string(argv[4])); // Load exact solution vector
				}
				BARRIER
				double dif1 = 0, dif2 = 0, dif8 = 0, norm = 0;
				//std::cout << "[" << rank << "]: mbeg=" << mbeg << " mend=" << mend << std::endl;
				for( k = mbeg; k < mend; ++k )
				{
					double dloc = err[k] = fabs(x[k] - ex[k]);
					dif1 += dloc;
					dif2 += dloc*dloc;
					dif8 = (dif8 > dloc) ? dif8 : dloc;
					norm += fabs(ex[k]);
				}
#if defined(USE_MPI)
				if( procs > 1 )
				{
					double temp = dif1;
					MPI_Reduce(&temp,&dif1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					temp = dif2;
					MPI_Reduce(&temp,&dif2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
					temp = dif8;
					MPI_Reduce(&temp,&dif8,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
					temp = norm;
					MPI_Reduce(&temp,&norm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
				}
#endif
				dif2 = sqrt(dif2);
				norm += 1.0e-100;
				if( !rank )
				{
					std::cout << "Difference with exact solution \"" << std::string(argv[4]) << "\": " << std::scientific << std::setprecision(6) << std::endl;
					std::cout << "dif1 = " << dif1 << "  dif2 = " << dif2 << "  dif8 = " << dif8 << "  ||ex||_1 = " << norm << std::endl;
					std::cout << "rel1 = " << dif1/norm << "  rel2 = " << dif2/norm << "  rel8 = " << dif8/norm << std::endl;
				}
				//err.Save("error.sol");
			}
		}
	}

	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
