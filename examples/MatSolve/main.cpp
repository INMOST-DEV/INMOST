#include <string>
#include <iostream>

#include "../../inmost.h"
using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

int main(int argc, char ** argv)
{
	int rank,procs;
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " method_number<0:INNER_ILU2,1:INNER_MLILUC,2:PETSc,3:Trilinos_Aztec,4:Trilinos_Belos,5:Trilinos_Ifpack,6:Trilinos_ML,7:ANI> matrix.mtx [right_hand_side.rhs] [solver_options.txt]" << std::endl;
		return -1;
	}
	Solver::Type type;
	switch(atoi(argv[1]))
	{
		case 0: type = Solver::INNER_ILU2; break;
		case 1: type = Solver::INNER_DDPQILUC; break;
		case 2: type = Solver::PETSc; break;
		case 3: type = Solver::Trilinos_Aztec; break;
		case 4: type = Solver::Trilinos_Belos; break;
		case 5: type = Solver::Trilinos_Ifpack; break;
		case 6: type = Solver::Trilinos_ML; break;
    case 7: type = Solver::ANI; break;
		case 8: type = Solver::INNER_MPTILUC; break;
	}
	Solver::Initialize(&argc,&argv,argc > 4 ? argv[4] : NULL); // Initialize the linear solver in accordance with args
	{
#if defined(USE_MPI)
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);  // Get the rank of the current process
		MPI_Comm_size(MPI_COMM_WORLD,&procs); // Get the total number of processors used
#else
		rank = 0;
		procs = 1;
#endif
		//std::cout << rank << "/" << procs << " " << argv[0] << std::endl;
		Solver::Matrix mat("A"); // Declare the matrix of the linear system to be solved
		Solver::Vector b("rhs"); // Declare the right-hand side vector
		Solver::Vector x("sol"); // Declare the solution vector
		//std::cout << rank << " load matrix from " << std::string(argv[2]) << " ..." << std::endl;
		double t = Timer(), tt = Timer();
		mat.Load(std::string(argv[2])); //if interval parameters not set, matrix will be divided automatically
		BARRIER
		if( !rank ) std::cout << "load matrix: " << Timer() - t << std::endl;
		//mat.Save("test.mtx");
		t = Timer();
		if( argc > 3 )
		{
			//std::cout << rank << " load vector from " << std::string(argv[3]) << std::endl;
			b.Load(std::string(argv[3])); // Load RHS vector
		}
		else // Set local RHS to 1 if it was not specified
		{
			INMOST_DATA_ENUM_TYPE mbeg,mend,k;
			mat.GetInterval(mbeg,mend);
			b.SetInterval(mbeg,mend);
			for( k = mbeg; k < mend; ++k ) b[k] = 1.0;
		}
		BARRIER
		if( !rank ) std::cout << "load vector: " << Timer() - t << std::endl;
		bool success = false;
		int iters;
		double resid, realresid = 0;
		std::string reason;
		{
			Solver s(type); // Declare the linear solver by specified type

			s.SetParameterEnum("gmres_substeps",4);
			s.SetParameterReal("relative_tolerance",1.0e-9);
			s.SetParameterReal("absolute_tolerance",1.0e-16);

			s.SetParameterEnum("reorder_nonzeros",0);
			s.SetParameterEnum("rescale_iterations",8);
			s.SetParameterEnum("adapt_ddpq_tolerance",0);
			
			s.SetParameterReal("drop_tolerance",1.0e-4);
			s.SetParameterReal("reuse_tolerance",1.0e-8);
			s.SetParameterReal("ddpq_tolerance",0.8);

			s.SetParameterEnum("condition_estimation",1);
			

			t = Timer();
			s.SetMatrix(mat); // Compute the preconditioner for the original matrix
			BARRIER
			if( !rank ) std::cout << "preconditioner: " << Timer() - t << std::endl;
			t = Timer();
			success = s.Solve(b,x); // Solve the linear system with the previously computted preconditioner
			BARRIER
			if( !rank ) std::cout << "solver: " << Timer() - t << "\t\t\t" << std::endl;
			iters = s.Iterations(); // Get the number of iterations performed
			resid = s.Residual();   // Get the final residual achieved
			reason = s.GetReason();
			//x.Save("output.rhs");
		}
		tt = Timer() - tt;

		{ // Compute the true residual
			double aresid = 0, bresid = 0;
			Solver::Vector test;
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
			if( info.GetRank() == 0 ) std::cout << "||Ax-b|| " << sqrt(recv[0]) << " ||b|| " << sqrt(recv[1]) << " ||Ax-b||/||b|| " << sqrt(recv[0]/(recv[1]+1.0e-100)) << std::endl;
#endif
			realresid = sqrt(recv[0]/(recv[1]+1.0e-100));
			//realresid = sqrt(realresid);

			info.RestoreVector(x);
			if( !rank ) std::cout << "norms: " << Timer() - t << std::endl;
		}

		{
			if( rank == 0 )
			{
				std::cout << procs << " processors";
				if( success )
				{
					std::cout << " solved in " << tt << " secs";
					std::cout << " with " << iters << " iterations to " << resid << " norm";
				}
				else std::cout << " failed to solve";
				std::cout  << " matrix \"" << argv[2] << "\"";
				if( argc > 3 ) std::cout  << " vector \"" << argv[3] << "\"";
				std::cout << " true residual ||Ax-b||/||b|| " << realresid;
				std::cout << std::endl;
				std::cout << "reason: " << reason << std::endl;
			}
		}
	}

	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}

