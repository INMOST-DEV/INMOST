#include <cstdio>
#include <cmath>

#include "../../inmost.h"
using namespace INMOST;

int main(int argc,char ** argv)
{
	int permut = 0;
	int solver = 0; // 0 - INNER_ILU2, 1 - INNER_MLILUC, 2 - PETSc
	// 3 - Trilinos_Aztec, 4 - Trilinos_Ifpack,
	// 5 - Trilinos_ML, 6 - Trilinos_Belos, 7 - ANI
	int rank,procs,newrank;

	Solver::Initialize(&argc,&argv,"database.txt"); // Initialize the solver and MPI activity
#if defined(USE_MPI)
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);  // Get the rank of the current process
	MPI_Comm_size(MPI_COMM_WORLD,&procs); // Get the total number of processors used
#else
	rank = 0;
	procs = 1;
#endif

	if (argc > 1)  permut = atoi(argv[1]);
	if (argc > 2)  solver = atoi(argv[2]);

	if (permut < procs)  newrank = (rank + permut) % procs;
	else newrank = (permut - rank) % procs;

	std::cout << rank <<  " -> " << newrank << std::endl;

	Solver::Type type;
	switch(solver)
	{
	case 0: type = Solver::INNER_ILU2; break;
	case 1: type = Solver::INNER_DDPQILUC; break;
	case 2: type = Solver::PETSc; break;
	case 3: type = Solver::Trilinos_Aztec; break;
	case 4: type = Solver::Trilinos_Ifpack; break;
	case 5: type = Solver::Trilinos_ML; break;
	case 6: type = Solver::Trilinos_Belos; break;
	case 7: type = Solver::ANI; break;
	case 8: type = Solver::INNER_MPTILUC; break;
	case 9: type = Solver::INNER_MPTILU2; break;
	}

	{
		Solver S(type); // Specify the linear solver


		Solver::Matrix A; // Declare the matrix of the linear system to be solved
		Solver::Vector x,b; // Declare the solution and the right-hand side vectors


		INMOST_DATA_ENUM_TYPE mbeg, mend;

		mbeg = newrank * 10;
		mend = (newrank+1) * 10;
		A.SetInterval(mbeg,mend);
		x.SetInterval(mbeg,mend);
		b.SetInterval(mbeg,mend);

		for( INMOST_DATA_ENUM_TYPE i = mbeg; i !=  mend; i++ )
		{
			A[i][i] = 10.0;
			b[i] = i;
			x[i] = 0.1;
		}
		if (rank==0)  std::cout << "next call S.SetMatrix(A);" << std::endl;
		S.SetMatrix(A); // Compute the preconditioner for the original matrix
		if (rank==0)  std::cout << "next call S.Solve(b,x);" << std::endl;
		if( !S.Solve(b,x) )   // Solve the linear system with the previously computted preconditioner
		{
			if( rank == 0 )
				std::cout << S.GetReason() << std::endl;
#if defined(USE_MPI)
			MPI_Abort(MPI_COMM_WORLD,-1);
#else
			return -1;
#endif

		}
		

		double err = 0;
		for( INMOST_DATA_ENUM_TYPE i = mbeg; i !=  mend; i++ )
		{
			err += fabs(x[i] - static_cast<double>(i)/10.0);
		}
#if defined(USE_MPI)
		double tmp;
		MPI_Allreduce(&err,&tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		err = tmp;
#endif

		if (rank==0)  
		{
			std::cout << "done, error " << err << "\t\t\t" << std::endl;
		}

		if( err > 1.0e-5 || err!=err)
		{
#if defined(USE_MPI)
			MPI_Abort(MPI_COMM_WORLD,-1);
#else
			return -1;
#endif
		}
		//no pollution with files
		//x.Save("output.rhs");
		//b.Save("b.rhs");
		//A.Save("A.mtx");
	}

	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
