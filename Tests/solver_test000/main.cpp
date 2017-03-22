#include <cstdio>
#include <cmath>

#include "inmost.h"
using namespace INMOST;

int main(int argc,char ** argv)
{
	int permut = 0;
	int copy_test = 0;
	std::string solver = "inner_ilu2";
	int rank,procs,newrank;

	Solver::Initialize(&argc,&argv,"database.xml"); // Initialize the solver and MPI activity
#if defined(USE_MPI)
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);  // Get the rank of the current process
	MPI_Comm_size(MPI_COMM_WORLD,&procs); // Get the total number of processors used
#else
	rank = 0;
	procs = 1;
#endif

	if (argc > 1)  permut = atoi(argv[1]);
	if (argc > 2)  solver = std::string(argv[2]);
	if (argc > 3)  copy_test = atoi(argv[3]);

	if (permut < procs)  newrank = (rank + permut) % procs;
	else newrank = (permut - rank) % procs;

	std::cout << rank <<  " -> " << newrank << std::endl;

	{
		Solver S(solver); // Specify the linear solver
        Solver &actualSolver = S;

        //Copy test
        if (copy_test == 1) {
            Solver copySolver(S);
            actualSolver = copySolver;
        }

		Sparse::Matrix A; // Declare the matrix of the linear system to be solved
		Sparse::Vector x,b; // Declare the solution and the right-hand side vectors


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
		if (rank==0)  std::cout << "next call SetMatrix(A);" << std::endl;
		actualSolver.SetMatrix(A); // Compute the preconditioner for the original matrix
		if (rank==0)  std::cout << "next call Solve(b,x);" << std::endl;

        //Assign test
        if (copy_test == 2) {
            Solver assignSolver(solver);
            assignSolver = Solver(actualSolver);
            actualSolver = assignSolver;
        }

		if( !actualSolver.Solve(b,x) )   // Solve the linear system with the previously computted preconditioner
		{
			if( rank == 0 )
				std::cout << actualSolver.ReturnReason() << std::endl;
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
#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
	Solver::Finalize(); // Finalize solver and close MPI activity
	return 0;
}
