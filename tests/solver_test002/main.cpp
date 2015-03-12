#include <string>
#include <iostream>

#include "../../inmost.h"
using namespace INMOST;

/***********************************************************************

(*) 5402 3D Poisson equation solver

(*) Test solvers on 3D Poisson equation.

This test is located in Tests/solver_test002

(*) Brief

Test solvers and tune parameters for on-the-fly generated matrices for 3D
Poisson equation.

(*) Description

This test will run solvers in both serial and parallel modes with NP 
processes for a model problem of 3D Poisson equation.
The coefficient matrix and vectors for a parallel run will be generated
directly at the target processor, so no external reordering is required
to be installed.
The artificial right-hand side rhs=(1,1,...,1) is used.
The specific solver is defined by a user.
User may also provide options file to alter default solver options.

Main purpose of this test is to assess robustness of internal or external
solvers during development.
Another purpose is to check the behaviour of the liner solver for large
and extremely large test problems without taking into account the disk
memory requirements.

(*) Arguments

Usage: ./solver_test002 method_number<0:INNER_ILU2,1:INNER_MLILUC,2:PETSc,3:Trilinos_Aztec,4:Trilinos_Belos,5:Trilinos_Ifpack,6:Trilinos_ML,7:ANI> N<for NxNxN problem> [solver_options.txt]


    * First parameter is the Solver type:
        0 – INNER_ILU2, inner Solver based on BiCGStab(L) solver with second
order ILU factorization as preconditioner;
        1 – INNER_MLILUC, inner Solver based on BiCGStab(L) solver with second
order Crout-ILU with inversed-based condition estimation and unsymmetric
reordering for diagonal dominance as preconditioner;
        2 – PETSc, external Solver AztecOO from Trilinos package;
        3 – Trilinos_Aztec, external Solver Belos from Trilinos package,
currently without preconditioner;
        4 – Trilinos_Belos, external Solver AztecOO with ML preconditioner;
        5 – Trilinos_Ifpack, external Solver AztecOO with Ifpack
preconditioner;
        6 – Trilinos_ML, external Solver PETSc;
        7 – ANI, external Solver from ANI3D based on ILU2 (sequential Fortran
version).
    * Second parameter is the dimension N of the 3D Poisson problem for NxNxN
mesh.
    * Third optional parameter is the file with solver parameters, see
examples/MatSolve/database.txt as example.

(*) Running test

You can run the test directly from the command line.
For example, you can specify the 100x100x100 test case and solve it by the 
internal ILU2 based solver with the default parameters on 4 processors:

$ cd tests/solver_test002
$ mpirun -np 4 ./solver_test002 0 100 

(*) Source

    * Source code is adopted from examples/MatSolve

***********************************************************************/

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

void Poisson3D(int n1, int n2, int n3, Solver::Matrix & mat); // Generate 3D Poisson matrix

int main(int argc, char ** argv)
{
	int rank,procs;
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " method_number<0:INNER_ILU2,1:INNER_MLILUC,2:PETSc,3:Trilinos_Aztec,4:Trilinos_Belos,5:Trilinos_Ifpack,6:Trilinos_ML,7:ANI> N<for NxNxN problem> [solver_options.txt]" << std::endl;
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
	int n = atoi(argv[2]);
	Solver::Initialize(&argc,&argv,argc > 3 ? argv[3] : NULL); // Initialize the linear solver in accordance with args
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
		int n1=n, n2=n, n3=n;
		if( !rank ) std::cout << "Poisson equation: " << n1 << " x " << n2 << " x " << n3 << std::endl;
		BARRIER
		double tt = Timer(), t = Timer();
		Poisson3D(n1,n2,n3,mat);
		BARRIER
		if( !rank ) std::cout << "Generate matrix: " << Timer() - t << std::endl;
		//mat.Save("test.mtx");
		// Set local RHS to 1 if it was not specified
		INMOST_DATA_ENUM_TYPE mbeg,mend,k;
		mat.GetInterval(mbeg,mend);
		b.SetInterval(mbeg,mend);
		for( k = mbeg; k < mend; ++k ) b[k] = 1.0;
		bool success = false;
		int iters;
		double resid, realresid = 0;
		std::string reason;
		{
			Solver s(type); // Declare the linear solver by specified type

			s.SetParameterEnum("gmres_substeps",3);
			
			s.SetParameterEnum("reorder_nonzeros",0);
			s.SetParameterEnum("rescale_iterations",8);
			s.SetParameterEnum("adapt_ddpq_tolerance",0);
			
			s.SetParameterReal("drop_tolerance",0.001);
			s.SetParameterReal("reuse_tolerance",0.00001);
			s.SetParameterReal("ddpq_tolerance",0.7);
			

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

// Generate 3D Poisson matrix
void Poisson3D(int n1, int n2, int n3, Solver::Matrix & A)
{
    int myid=0, nproc=1;
#if defined(USE_MPI)
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
#endif

    int n = n1 * n2 * n3;
    unsigned idmax, idmin, block = n / nproc;
    idmin = myid * block;
    idmax = idmin + block;
    if (myid == nproc-1) idmax = n;
    A.SetInterval(idmin,idmax);

    static const unsigned ndiag = 7;
    int id[ndiag] = {  0, 0,-1, 0, 1, 0, 0 }; // !    !    -1    !    !
    int jd[ndiag] = {  0,-1, 0, 0, 0, 1, 0 }; // ! -1 ! -1  6 -1 ! -1 !
    int kd[ndiag] = { -1, 0, 0, 0, 0, 0, 1 }; // !    !    -1    !    !
    int ad[ndiag] = { -1,-1,-1, 6,-1,-1,-1 };

    for (int k=0; k<n3; k++) {
        for (int j=0; j<n2; j++) {
            for (int i=0; i<n1; i++) {
                unsigned iii = i + j*n1 + k*n1*n2;
                if (iii >= idmin && iii < idmax) {
                    for (unsigned d=0; d<ndiag; d++) {
                        int ii = i + id[d];
                        int jj = j + jd[d];
                        int kk = k + kd[d];
                        if (ii >= 0 && ii < n1 &&
                            jj >= 0 && jj < n2 &&
                            kk >= 0 && kk < n3) {
                            int jjj = ii + jj*n1 + kk*n1*n2;
                            A[iii][jjj] = ad[d];
                        }
                    }
                }
            }
        }
    }

}
