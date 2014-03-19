
#include "../../inmost.h"
#include <string>
#include <iostream>

using namespace INMOST;
#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

int main(int argc, char ** argv)
{
	/*
	int i;
	for (int j = 0; j < 10000; j++) 
	{
		Solver::Matrix A;
		double t = Timer();
		for(int i = 0; i < 10000; i++)
			A[i][i] = 0;
		printf("%d %g\n",j,Timer()-t);
	}
	return 0;
	*/
	int rank,procs;
	if( argc < 3 )
	{
		std::cout << "usage: " << argv[0] << " method_number<0:Inner,1:ANI3D,2:PETSc,3:Trilinos> matrix.mtx [right_hand_side.rhs] [petsc_options.txt]" << std::endl;
		return -1;
	}
	Solver::Type type;
	switch(atoi(argv[1]))
	{
		case 0: type = Solver::INNER_ILU2; break;
		case 1: type = Solver::INNER_MLILUC; break;
		case 2: type = Solver::ANI; break;
		case 3: type = Solver::PETSC; break;
		
	}
	Solver::Initialize(&argc,&argv,argc > 4 ? argv[4] : NULL);
	{
#if defined(USE_MPI)
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&procs);
#else
		rank = 0;
		procs = 1;
#endif
		//std::cout << rank << "/" << procs << " " << argv[0] << std::endl;
		Solver::Matrix mat("A");
		Solver::Vector b("rhs");
		Solver::Vector x("sol");
		//std::cout << rank << " load matrix from " << std::string(argv[2]) << std::endl;
		long double t = Timer(), tt = Timer();
		mat.Load(std::string(argv[2])); //if interval parameters not set, matrix will be divided automatically
		BARRIER
		if( !rank ) std::cout << "load matrix: " << Timer() - t << std::endl;
		//mat.Save("test.mtx");
		t = Timer();
		if( argc > 3 ) 
		{
			//std::cout << rank << " load vector from " << std::string(argv[3]) << std::endl;
			b.Load(std::string(argv[3]));
		}
		else
		{
			INMOST_DATA_ENUM_TYPE mbeg,mend,k;
			mat.GetInterval(mbeg,mend);
			b.SetInterval(mbeg,mend);
			for(k = mbeg; k < mend; k++) b[k] = 1.0;
		}
		BARRIER
		if( !rank ) std::cout << "load vector: " << Timer() - t << std::endl;
		bool success = false;
		int iters;
		double resid, realresid = 0;
		
		{
			Solver s(type);
			t = Timer();
			s.SetMatrix(mat);
			BARRIER
			if( !rank ) std::cout << "preconditioner: " << Timer() - t << std::endl;
			t = Timer();
			success = s.Solve(b,x);
			BARRIER
			if( !rank ) std::cout << "solver: " << Timer() - t << std::endl;
			iters = s.Iterations();
			resid = s.Residual();
			//x.Save("output.rhs");
		}
		tt = Timer() - tt;
		
		{
			double aresid = 0, bresid = 0;
			Solver::Vector test;
			t = Timer();
			Solver::OrderInfo info;
			info.PrepareMatrix(mat,0);
			info.PrepareVector(x);
			info.Update(x);
			
			
			mat.MatVec(1.0,x,0.0,test);
		
			{
				INMOST_DATA_ENUM_TYPE mbeg,mend,k;
				info.GetLocalRegion(info.GetRank(),mbeg,mend);
				for(k = mbeg; k < mend; k++)
				{
					aresid += (test[k]-b[k])*(test[k]-b[k]);
					bresid += b[k]*b[k];
				}
			}
			double temp[2] = {aresid,bresid}, recv[2] = {aresid,bresid};
#if defined(USE_MPI)
			MPI_Reduce(temp,recv,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			if( info.GetRank() == 0 ) std::cout << "||Ax-b|| " << sqrt(recv[0]) << " ||b|| " << sqrt(recv[1]) << " ||Ax-b||/||b|| " << sqrt(recv[0]/recv[1]) << std::endl;
#endif
			realresid = sqrt(recv[0]/recv[1]);
			//realresid = sqrt(realresid);
			
			info.RestoreVector(x);
			if( !rank ) std::cout << "norms: " << Timer() - t << std::endl;
		}
		
		
		{
			if( rank == 0 ) 
			{
				std::cout  << procs << " processors";
				if( success )
				{
					std::cout  << " solved in " << tt << " secs";
					std::cout  << " with " << iters << " iterations to " << resid << " norm";
				}
				else std::cout << " failed to solve";
				std::cout  << " matrix " << argv[2];
				if( argc > 3 ) std::cout  << " vector " << argv[3];
				std::cout << " real residual ||Ax-b||/||b|| " << realresid;
				std::cout << std::endl;
			}
		}
	}
	Solver::Finalize();
	return 0;
}
