#include "inmost.h"
#include <cstdio>
#include <cmath>
using namespace INMOST;

int main(int argc,char ** argv)
{
	int test = 0;
	if (argc > 1)  test = atoi(argv[1]);

	Storage::real err = 1;

	if( test == 0 ) //A*B
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		Storage::real B[] =
		{
			4, 2,
			3, 5,
			9, 8,
			1, 7
		};
		Storage::real C[] =
		{
			36,  68,
			55,  69,
			103,  90
		};
		err = (raMatrixMake(A,3,4)*raMatrixMake(B,4,2)-raMatrixMake(C,3,2)).FrobeniusNorm();
	}
	else if( test == 1 ) // A/B == B^{-1}*A
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		Storage::real B[] =
		{
			4, 2,
			9, 8,
			2, 1
		};
		Storage::real C[] =
		{
			0.942857142857144,  0.685714285714285,  2.685714285714290,  2.142857142857140,
			-0.685714285714285,  0.228571428571429, -2.771428571428570, -2.285714285714290
		};
		err = (raMatrixMake(A,3,4)/raMatrixMake(B,3,2)-raMatrixMake(C,2,4)).FrobeniusNorm();
	}
	else if( test == 2 ) // A+B
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		Storage::real B[] =
		{
			2, 6, 9, 1,
			9, 1, 5, 3,
			8, 21, 91, 9
		};
		Storage::real C[] =
		{
			3, 9, 11, 6,
			12, 9, 7, 4,
			12, 23, 100, 9
		};
		err = (raMatrixMake(A,3,4)+raMatrixMake(B,3,4)-raMatrixMake(C,3,4)).FrobeniusNorm();
	}
	else if( test == 3 ) // A-B
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		Storage::real B[] =
		{
			2, 6, 9, 1,
			9, 1, 5, 3,
			8, 21, 91, 9
		};
		Storage::real C[] =
		{
			-1, -3, -7, 4,
			-6, 7, -3, -2,
			-4, -19, -82, -9
		};
		err = (raMatrixMake(A,3,4)-raMatrixMake(B,3,4)-raMatrixMake(C,3,4)).FrobeniusNorm();
	}
	else if( test == 4 ) // svd square
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0,
			2, 1, 10, 9
		};
		rMatrix U, S, V;
		raMatrix mA = raMatrixMake(A,4,4);
		mA.SVD(U,S,V);
		err = (U*S*V.Transpose() - mA).FrobeniusNorm();
	}
	else if( test == 5 ) // svd n > m
	{
		Storage::real A[] =
		{
			1, 3, 2,
			3, 8, 2,
			4, 2, 9,
			2, 1, 10
		};
		rMatrix U, S, V;
		raMatrix mA = raMatrixMake(A,4,3);
		mA.SVD(U,S,V);
		err = (U*S*V.Transpose() - mA).FrobeniusNorm();
	}
	else if( test == 6 ) // svd n < m
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		rMatrix U, S, V;
		raMatrix mA = raMatrixMake(A,3,4);
		mA.SVD(U,S,V);
		err = (U*S*V.Transpose() - mA).FrobeniusNorm();
	}
	else if( test == 7 ) // inv square
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0,
			2, 1, 10, 9
		};
		raMatrix mA = raMatrixMake(A,4,4);
		err = (mA.Invert()*mA-rMatrix::Unit(4)).FrobeniusNorm();
		err+= (mA*mA.Invert()-rMatrix::Unit(4)).FrobeniusNorm();
	}
	else if( test == 8 ) // inv n > m
	{
		Storage::real A[] =
		{
			1, 3, 2,
			3, 8, 2,
			4, 2, 9,
			2, 1, 10
		};
		raMatrix mA = raMatrixMake(A,4,3);
		err = (mA.Invert()*mA-rMatrix::Unit(3)).FrobeniusNorm();
	}
	else if( test == 9 ) // inv n < m
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		raMatrix mA = raMatrixMake(A,3,4);
		err = (mA*mA.Invert()-rMatrix::Unit(3)).FrobeniusNorm();
	}
	else if( test == 10 ) // cholinv square
	{
		Storage::real A[] =
		{
			10, 3, 2, 5,
			3, 8, 2, 1,
			2, 2, 9, 0,
			5, 1, 0, 9
		};
		raMatrix mA = raMatrixMake(A,4,4);
		err = (mA.CholeskyInvert()*mA-rMatrix::Unit(4)).FrobeniusNorm();
		err+= (mA*mA.CholeskyInvert()-rMatrix::Unit(4)).FrobeniusNorm();
	}
	else if( test == 11 ) // pinv square
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0,
			2, 1, 10, 9
		};
		raMatrix mA = raMatrixMake(A,4,4);
		err = (mA.PseudoInvert()*mA-rMatrix::Unit(4)).FrobeniusNorm();
		err+= (mA*mA.PseudoInvert()-rMatrix::Unit(4)).FrobeniusNorm();
	}
	else if( test == 12 ) // pinv n > m
	{
		Storage::real A[] =
		{
			1, 3, 2,
			3, 8, 2,
			4, 2, 9,
			2, 1, 10
		};
		raMatrix mA = raMatrixMake(A,4,3);
		err = (mA.PseudoInvert()*mA-rMatrix::Unit(3)).FrobeniusNorm();
	}
	else if( test == 13 ) // pinv n < m
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		raMatrix mA = raMatrixMake(A,3,4);
		err = (mA*mA.PseudoInvert()-rMatrix::Unit(3)).FrobeniusNorm();
	}
	else if( test == 14 ) // sol square
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0,
			2, 1, 10, 9
		};
		Storage::real B[] =
		{
			1, 2,
			3, 9,
			4, 1,
			2, 5
		};
		Storage::real C[] =
		{
			1.000000000000000, -8.333333333333330,
			0.000000000000000,  3.666666666666670,
			-0.000000000000001,  3.000000000000000,
			0.000000000000001, -1.333333333333330
		};
		raMatrix mA = raMatrixMake(A,4,4);
		raMatrix mB = raMatrixMake(B,4,2);
		raMatrix mC = raMatrixMake(C,4,2);
		err = (mA.Solve(mB)-mC).FrobeniusNorm();
	}
	else if( test == 15 ) // sol n > m
	{
		Storage::real A[] =
		{
			1, 3, 2,
			3, 8, 2,
			4, 2, 9,
			2, 1, 10
		};
		Storage::real B[] =
		{
			1, 2,
			3, 9,
			4, 1,
			2, 5
		};
		raMatrix mA = raMatrixMake(A,4,3);
		raMatrix mB = raMatrixMake(B,4,2);
		//raMatrix mC = raMatrixMake(C,4,2);
		//err = (mA.PseudoSolve(mB)-mC).FrobeniusNorm();
		err = (mA.Transpose()*(mA*mA.Solve(mB)-mB)).FrobeniusNorm();
	}
	else if( test == 16 ) // sol n < m
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		Storage::real B[] =
		{
			1, 2,
			3, 9,
			4, 1,
		};
		raMatrix mA = raMatrixMake(A,3,4);
		raMatrix mB = raMatrixMake(B,3,2);
		//raMatrix mC = raMatrixMake(C,4,2);
		//err = (mA.PseudoSolve(mB)-mC).FrobeniusNorm();
		err = (mA*mA.Solve(mB)-mB).FrobeniusNorm();
	}
	else if( test == 17 ) // cholsol square
	{
		Storage::real A[] =
		{
			10, 3, 2, 5,
			3, 8, 2, 1,
			2, 2, 9, 0,
			5, 1, 0, 9
		};
		Storage::real B[] =
		{
			1, 2,
			3, 9,
			4, 1,
			2, 5
		};
		Storage::real C[] =
		{
			-0.241562583045445, -0.514217379750200,
			0.318628753654002,  1.242625564709010,
			0.427318628753653, -0.050757374435302,
			0.321020462397020,  0.703162370449109
		};
		raMatrix mA = raMatrixMake(A,4,4);
		raMatrix mB = raMatrixMake(B,4,2);
		raMatrix mC = raMatrixMake(C,4,2);
		err = (mA.CholeskySolve(mB)-mC).FrobeniusNorm();
	}
	else if( test == 18 ) // psol square
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0,
			2, 1, 10, 9
		};
		Storage::real B[] =
		{
			1, 2,
			3, 9,
			4, 1,
			2, 5
		};
		Storage::real C[] =
		{
			1.000000000000000, -8.333333333333330,
			0.000000000000000,  3.666666666666670,
			-0.000000000000001,  3.000000000000000,
			0.000000000000001, -1.333333333333330
		};
		raMatrix mA = raMatrixMake(A,4,4);
		raMatrix mB = raMatrixMake(B,4,2);
		raMatrix mC = raMatrixMake(C,4,2);
		err = (mA.PseudoSolve(mB)-mC).FrobeniusNorm();
		//err = (mA*mA.PseudoSolve(mB)-mB).FrobeniusNorm();
	}
	else if( test == 19 ) // psol n > m
	{
		Storage::real A[] =
		{
			1, 3, 2,
			3, 8, 2,
			4, 2, 9,
			2, 1, 10
		};
		Storage::real B[] =
		{
			1, 2,
			3, 9,
			4, 1,
			2, 5
		};
		raMatrix mA = raMatrixMake(A,4,3);
		raMatrix mB = raMatrixMake(B,4,2);
		//raMatrix mC = raMatrixMake(C,4,2);
		//err = (mA.PseudoSolve(mB)-mC).FrobeniusNorm();
		err = (mA.Transpose()*(mA*mA.PseudoSolve(mB)-mB)).FrobeniusNorm();
	}
	else if( test == 20 ) // psol n < m
	{
		Storage::real A[] =
		{
			1, 3, 2, 5,
			3, 8, 2, 1,
			4, 2, 9, 0
		};
		Storage::real B[] =
		{
			1, 2,
			3, 9,
			4, 1,
		};
		raMatrix mA = raMatrixMake(A,3,4);
		raMatrix mB = raMatrixMake(B,3,2);
		//raMatrix mC = raMatrixMake(C,4,2);
		//err = (mA.PseudoSolve(mB)-mC).FrobeniusNorm();
		err = (mA*mA.PseudoSolve(mB)-mB).FrobeniusNorm();
	}
	else if( test == 21 ) // transform
	{
		Storage::real A[] = {2,4,9};
		Storage::real B[] = {-1,5,-21};
		raMatrix mA = raMatrixMake(A,3,1);
		raMatrix mB = raMatrixMake(B,3,1);
		rMatrix Q = mA.Transform(mB);
		std::cout << "Q:" << std::endl;
		Q.Print();
		std::cout << "Q*mA:" << std::endl;
		(Q*mA).Print();
		std::cout << "Q*mB:" << std::endl;
		(Q*mB).Print();
		err = (mA.Transform(mB)*mA-mB).FrobeniusNorm();
	}
	else if (test == 22)
	{
		Storage::real A[] = { 1,2,3 };
		Storage::real B[] = { 2,3,4 };
		Storage::real C[] = { -1,2,-1 };
		raMatrix mA = raMatrixMake(A, 3, 1);
		raMatrix mB = raMatrixMake(B, 3, 1);
		raMatrix mC = raMatrixMake(C, 3, 1);
		err = (mA.CrossProduct(mB) - mC).FrobeniusNorm();
	}
	else if (test == 23)
	{
		Storage::real A[] = { 1,2,3 };
		Storage::real B[] = { 2,3,4 };
		Storage::real C[] = { -1,2,-1 };
		rMatrix Q = rMatrix::CrossProductMatrix(A);
		raMatrix mB = raMatrixMake(B, 3, 1);
		raMatrix mC = raMatrixMake(C, 3, 1);
		err = (Q*mB - mC).FrobeniusNorm();
	}
	else if (test == 24)  //check memory_pool with large matrix (was shown to fail in hydr_frac problem)
	{
#if defined(USE_OMP)
#pragma omp parallel
#endif
		{
			rMatrix A(224, 12), B(224, 1);
			for (int k = 0; k < 224; ++k)
			{
				for (int l = 0; l < 12; ++l)
					A(k, l) = 2.0 * (rand() / (1. * RAND_MAX)) - 1.0;
				B(k, 0) = 2.0 * (rand() / (1. * RAND_MAX)) - 1.0;
			}
#if defined(USE_OMP)
#pragma omp critical
#endif
			(A.PseudoSolve(B)).Transpose().Print();
		}
	}
#if defined(USE_FP64)
	if( fabs(err) > 1.0e-10 )
#else
	if( fabs(err) > 1.0e-5 )
#endif
	{
		std::cout << "error is " << err << std::endl;
		return -1;
	}
	else std::cout << "no error" << std::endl;
	return 0;
}
