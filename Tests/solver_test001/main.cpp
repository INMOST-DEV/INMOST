#include <string>
#include <iostream>

#include "inmost.h"
using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

int main(int argc, char ** argv)
{
	int rank,procs;
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " matrix.mtx [solver_type] [right_hand_side.rhs]" << std::endl;
		std::cout << "solver_type: inner_ilu2 (default), inner_mptilu2, inner_mptiluc, inner_ddpqiluc." << std::endl;
		return -1;
	}
	Solver::Type type = Solver::INNER_ILU2;
	Solver::Initialize(&argc,&argv,NULL); // Initialize the linear solver in accordance with args
	{
#if defined(USE_MPI)
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);  // Get the rank of the current process
		MPI_Comm_size(MPI_COMM_WORLD,&procs); // Get the total number of processors used
#else
		rank = 0;
		procs = 1;
#endif
		if( argc > 2 )
		{
			std::string stype(argv[2]);
			for(size_t k = 0; k < stype.size(); ++k) stype[k] = tolower(stype[k]);
			/*
			if( stype == "inner_mptilu2" )
				type = Solver::INNER_MPTILU2;
			else if( stype == "inner_mptiluc" )
				type = Solver::INNER_MPTILUC;
			else if( stype == "inner_ddpqiluc" )
				type = Solver::INNER_DDPQILUC;
			else
			 */
			if( Solver::isSolverAvailable(stype) )
				type = stype;
			else if( stype != "inner_ilu2" && stype != "*" )
			{
				std::cout << "Unknown type of the solver " << stype << std::endl;
				return -1;
			}
		}
		//std::cout << rank << "/" << procs << " " << argv[0] << std::endl;
		Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
		Sparse::Vector b("rhs"); // Declare the right-hand side vector
		Sparse::Vector x("sol"); // Declare the solution vector
		//std::cout << rank << " load matrix from " << std::string(argv[2]) << " ..." << std::endl;
		double t = Timer(), tt = Timer();
		mat.Load(std::string(argv[1])); //if interval parameters not set, matrix will be divided automatically
		BARRIER
		if( !rank ) std::cout << "load matrix: " << Timer() - t << std::endl;
		
		/*
		b.SetInterval(0,mat.Size());
		x.SetInterval(0,mat.Size());
		for(int k = 0; k < (int)mat.Size(); ++k)
			x[k] = (1+k)*10;
		mat.MatVec(1.0,x,0.0,b);
		for(int k = 0; k < (int)mat.Size(); ++k)
			x[k] = 0;
		*/
		
		//mat.Save("test.mtx");
		t = Timer();
		if( argc > 3 )
		{
			//std::cout << rank << " load vector from " << std::string(argv[3]) << std::endl;
			b.Load(std::string(argv[3])); // Load RHS vector
		}
		else if( b.Empty() )// Set local RHS to 1 if it was not specified
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
		std::string reason;
		double resid, realresid = 0;

		{
			Solver s(type); // Declare the linear solver by specified type

			s.SetParameterEnum("gmres_substeps",4);
			s.SetParameterEnum("condition_estimation",1);
			s.SetParameterReal("relative_tolerance",1.0e-9);
			s.SetParameterReal("absolute_tolerance",1.0e-16);
			s.SetParameterEnum("rescale_iterations",16);
			s.SetParameterReal("drop_tolerance",1.0e-2);
			s.SetParameterReal("reuse_tolerance",1.0e-4);
			//s.SetParameterReal("drop_tolerance",0);
			//s.SetParameterReal("reuse_tolerance",0);
			/*
			
			s.SetParameterEnum("reorder_nonzeros",0);
			s.SetParameterEnum("adapt_ddpq_tolerance",0);
			s.SetParameterReal("ddpq_tolerance",0.0);
			*/

			t = Timer();
			s.SetMatrix(mat); // Compute the preconditioner for the original matrix
			BARRIER
			if( !rank ) std::cout << "preconditioner: " << Timer() - t << std::endl;
			t = Timer();
			success = s.Solve(b,x); // Solve the linear system with the previously computted preconditioner
			BARRIER
			if( !rank ) std::cout << "solver: " << Timer() - t << std::endl;
			iters = s.Iterations(); // Get the number of iterations performed
			resid = s.Residual();   // Get the final residual achieved
			reason = s.GetReason();
			//x.Save("sol.mtx");
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
				//Sparse::Vector tmp("",mbeg,mend);
				for( k = mbeg; k < mend; ++k )
				{
					//tmp[k] = test[k]-b[k];
					aresid += (test[k]-b[k])*(test[k]-b[k]);
					bresid += b[k]*b[k];
				}
				//tmp.Save("diff.mtx");
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
					std::cout << " with " << iters << " iterations to " << resid << " norm" << std::endl;
				}
				else std::cout << " failed to solve" << std::endl;
				std::cout  << " matrix \"" << argv[1] << "\"" << std::endl;
				if( argc > 3 ) std::cout  << " vector \"" << argv[3] << "\"" << std::endl;
				else std::cout  << " unit right hand side" << std::endl;
				std::cout  << " solver \"" << type << "\"" << std::endl;
				std::cout << " true residual ||Ax-b||/||b|| " << realresid << std::endl;;
				std::cout << "converged reason: " << reason << std::endl;
			}
		}

		if( !success || (success && realresid > 1.0e-3) )
		{
#if defined(USE_MPI)
			MPI_Abort(MPI_COMM_WORLD,-1);
#else
			exit(-1);
#endif
		}
	}

	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
