
#ifndef __SOLVER_BCGS__
#define __SOLVER_BCGS__

//\todo
// 1. comply solvers with Method prototype, after TODO in solver_prototypes.hpp is done
// 2. Implement tricks from Read/solver/bcgsl/download.pdf with convex update and true residual correction
// 3. Detect numerical accuracy breakdown - when preconditioned residual is too far from true residual (probably 2 will fix).

#include "inmost_solver.h"

#define PSEUDOINVERSE  // same trick as in petsc with pseudoinverse
//#define USE_LAPACK_SVD // use lapack's dgesvd routine instead of built-in svdnxn

//#if !defined(NDEBUG)
#define REPORT_RESIDUAL
//#endif

namespace INMOST
{
	//lapack svd
#if defined(PSEUDOINVERSE)
#if defined(USE_LAPACK_SVD)
	extern "C"
	{
		void dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
	}
#else // SVD adopted from http://stackoverflow.com/questions/3856072/svd-implementation-c answer by Dhairya Malhotra
	void GivensL(INMOST_DATA_REAL_TYPE * S, const int N, int m, INMOST_DATA_REAL_TYPE a, INMOST_DATA_REAL_TYPE b)
	{
		INMOST_DATA_REAL_TYPE r = sqrt(a*a+b*b);
		INMOST_DATA_REAL_TYPE c = a/r;
		INMOST_DATA_REAL_TYPE s = -b/r;
		for(int i=0;i<N;i++)
		{
			INMOST_DATA_REAL_TYPE S0 = S[(m+0)*N+i];
			INMOST_DATA_REAL_TYPE S1 = S[(m+1)*N+i];
			S[(m+0)*N + i] += S0*(c-1);
			S[(m+0)*N + i] += S1*( -s);
			S[(m+1)*N + i] += S0*(s  );
			S[(m+1)*N + i] += S1*(c-1);
		}
	}

	void GivensR(INMOST_DATA_REAL_TYPE * S, const int N, int m, INMOST_DATA_REAL_TYPE a, INMOST_DATA_REAL_TYPE b)
	{
		INMOST_DATA_REAL_TYPE r = sqrt(a*a+b*b);
		INMOST_DATA_REAL_TYPE c = a/r;
		INMOST_DATA_REAL_TYPE s = -b/r;
		for(int i=0;i<N;i++)
		{
			INMOST_DATA_REAL_TYPE S0 = S[i*N+(m+0)];
			INMOST_DATA_REAL_TYPE S1 = S[i*N+(m+1)];
			S[i*N+(m+0)] += S0*(c-1);
			S[i*N+(m+0)] += S1*( -s);
			S[i*N+(m+1)] += S0*(s  );
			S[i*N+(m+1)] += S1*(c-1);
		}
	}

	void svdnxn(INMOST_DATA_REAL_TYPE * A, INMOST_DATA_REAL_TYPE * U, INMOST_DATA_REAL_TYPE * S,  INMOST_DATA_REAL_TYPE * V, const int N)
	{
		memset(S,0,sizeof(INMOST_DATA_REAL_TYPE)*N*N);
		memset(U,0,sizeof(INMOST_DATA_REAL_TYPE)*N*N);
		memset(V,0,sizeof(INMOST_DATA_REAL_TYPE)*N*N);
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
				S[i*N+j]=A[i*N+j];
			U[i*N+i] = 1;
			V[i*N+i] = 1;
		}
		INMOST_DATA_REAL_TYPE eps = -1;
		{ // Bi-diagonalization
			std::vector<INMOST_DATA_REAL_TYPE> house_vec(N);
			for(int i=0;i<N;i++)
			{
				// Column Householder
				{
					INMOST_DATA_REAL_TYPE x1= S[i*N+i];
					if(x1<0) x1=-x1;

					INMOST_DATA_REAL_TYPE x_inv_norm=0;
					for(int j=i;j<N;j++) x_inv_norm+= S[j*N+i]*S[j*N+i];
					x_inv_norm=1/sqrt(x_inv_norm);

					INMOST_DATA_REAL_TYPE alpha=sqrt(1+x1*x_inv_norm);
					INMOST_DATA_REAL_TYPE beta=x_inv_norm/alpha;

					house_vec[i]=-alpha;
					for(int j=i+1;j<N;j++) house_vec[j]=-beta*S[j*N+i];
					if(S[i*N+i]<0) for(int j=i+1;j<N;j++) house_vec[j]=-house_vec[j];
				}
				
				for(int k=i;k<N;k++)
				{
					INMOST_DATA_REAL_TYPE dot_prod=0;
					for(int j=i;j<N;j++) dot_prod+=S[j*N+k]*house_vec[j];
					
					for(int j=i;j<N;j++) S[j*N+k]-=dot_prod*house_vec[j];
				}
				
				for(int k=0;k<N;k++)
				{
					INMOST_DATA_REAL_TYPE dot_prod=0;
					for(int j=i;j<N;j++) dot_prod+=U[k*N+j]*house_vec[j];
					for(int j=i;j<N;j++) U[k*N+j]-=dot_prod*house_vec[j];
				}

				// Row Householder
				if(i>=N-1) continue;

				{
					INMOST_DATA_REAL_TYPE x1=S[i*N+(i+1)];
					if(x1<0) x1=-x1;

					INMOST_DATA_REAL_TYPE x_inv_norm=0;
					for(int j=i+1;j<N;j++) x_inv_norm+=S[i*N+j]*S[i*N+j];
					x_inv_norm=1/sqrt(x_inv_norm);

					INMOST_DATA_REAL_TYPE alpha=sqrt(1+x1*x_inv_norm);
					INMOST_DATA_REAL_TYPE beta=x_inv_norm/alpha;

					house_vec[i+1]=-alpha;
					for(int j=i+2;j<N;j++) house_vec[j]=-beta*S[i*N+j];
					if(S[i*N+(i+1)]<0) for(int j=i+2;j<N;j++) house_vec[j]=-house_vec[j];
				}
				
				for(int k=i;k<N;k++)
				{
					INMOST_DATA_REAL_TYPE dot_prod=0;
					for(int j=i+1;j<N;j++) dot_prod+=S[k*N+j]*house_vec[j];
					for(int j=i+1;j<N;j++) S[k*N+j] -= dot_prod*house_vec[j];
				}
				
				for(int k=0;k<N;k++)
				{
					INMOST_DATA_REAL_TYPE dot_prod=0;
					for(int j=i+1;j<N;j++) dot_prod+=V[j*N+k]*house_vec[j];
					for(int j=i+1;j<N;j++) V[j*N+k]-=dot_prod*house_vec[j];
				}
			}
		}

		int k0=0;
		if(eps<0)
		{
			eps=1.0;
			while(eps+(INMOST_DATA_REAL_TYPE)1.0>1.0) eps*=0.5;
			eps*=64.0;
		}
		while(k0<N-1)
		{ // Diagonalization
			INMOST_DATA_REAL_TYPE S_max=0.0;
			for(int i=0;i<N;i++) S_max=(S_max > S[i*N+i] ? S_max : S[i*N+i]);
			while(k0<N-1 && fabs(S[k0*N+(k0+1)])<=eps*S_max) k0++;
			int k=k0;
			int n=k0+1;
			while(n<N && fabs(S[(n-1)*N+n])>eps*S_max) n++;

			INMOST_DATA_REAL_TYPE mu=0;
			{ // Compute mu
				INMOST_DATA_REAL_TYPE C[2][2];
				C[0][0]=S[(n-2)*N+(n-2)]*S[(n-2)*N+(n-2)]+S[(n-3)*N+(n-2)]*S[(n-3)*N+(n-2)]; 
				C[0][1]=S[(n-2)*N+(n-2)]*S[(n-2)*N+(n-1)];
				C[1][0]=S[(n-2)*N+(n-2)]*S[(n-2)*N+(n-1)]; 
				C[1][1]=S[(n-1)*N+(n-1)]*S[(n-1)*N+(n-1)]+S[(n-2)*N+(n-1)]*S[(n-2)*N+(n-1)];
				INMOST_DATA_REAL_TYPE b =-(C[0][0]+C[1][1])/2;
				INMOST_DATA_REAL_TYPE c =  C[0][0]*C[1][1] - C[0][1]*C[1][0];
				INMOST_DATA_REAL_TYPE d = sqrt(b*b-c);
				INMOST_DATA_REAL_TYPE lambda1 = -b+d;
				INMOST_DATA_REAL_TYPE lambda2 = -b-d;
				INMOST_DATA_REAL_TYPE d1 = lambda1-C[1][1]; d1 = (d1<0?-d1:d1);
				INMOST_DATA_REAL_TYPE d2 = lambda2-C[1][1]; d2 = (d2<0?-d2:d2);
				mu = (d1<d2?lambda1:lambda2);
			}

			INMOST_DATA_REAL_TYPE alpha = S[k*N+k] * S[k*N+k] - mu;
			INMOST_DATA_REAL_TYPE beta  = S[k*N+k] * S[k*N+(k+1)];

			for(;k<N-1;k++)
			{
				GivensR(S,N,k,alpha,beta);
				GivensL(V,N,k,alpha,beta);

				alpha = S[k*N+k];
				beta  = S[(k+1)*N+k];
				GivensL(S,N,k,alpha,beta);
				GivensR(U,N,k,alpha,beta);

				alpha = S[k*N+(k+1)];
				beta  = S[k*N+(k+2)];
			}
		}
		
		for(int i=0;i<N;i++)
		{
			INMOST_DATA_REAL_TYPE temp;
			for(int j=i+1;j<N;j++)
			{
				temp = U[i+N*j];
				U[i+N*j] = U[j+N*i];
				U[j+N*i] = temp;

				temp = V[i+N*j];
				V[i+N*j] = V[j+N*i];
				V[j+N*i] = temp;
			}
		}
		

		for(int i=0;i<N;i++) if( S[i*N+i] < 0.0 )
		{
			for(int j=0;j<N;j++)
			{
				U[j+N*i] *= -1;
			}
			S[i*N+i] *= -1;
		}
		
	}
#endif //USE_LAPACK_SVD
#endif //PSEUDOINVERSE
	int solvenxn(INMOST_DATA_REAL_TYPE * A, INMOST_DATA_REAL_TYPE * x, INMOST_DATA_REAL_TYPE * b, int n, int * order)
	{
		INMOST_DATA_REAL_TYPE temp, max;
		int temp2;
		for(int i = 0; i < n; i++) order[i] = i;
		for(int i = 0; i < n; i++)
		{
			int maxk = i, maxq = i;
			max = fabs(A[maxk*n+maxq]);
			//Find best pivot
			for(int q = i; q < n; q++) // over columns
			{
				for(int k = i; k < n; k++) // over rows
				{
					if( fabs(A[k*n+q]) > max )
					{
						max = fabs(A[k*n+q]);
						maxk = k;
						maxq = q;
					}
				}
			}
			//Exchange rows
			if( maxk != i ) 
			{
				for(int q = 0; q < n; q++)
				{
					temp = A[maxk*n+q];
					A[maxk*n+q] = A[i*n+q];
					A[i*n+q] = temp;
				}
				//exchange rhs
				{
					temp = b[maxk];
					b[maxk] = b[i];
					b[i] = temp;
				}
			}
			//Exchange columns
			if( maxq != i ) 
			{
				for(int k = 0; k < n; k++)
				{
					temp = A[k*n+maxq];
					A[k*n+maxq] = A[k*n+i];
					A[k*n+i] = temp;
				}
				//remember order in sol
				{
					temp2 = order[maxq];
					order[maxq] = order[i];
					order[i] = temp2;
				}
			}
			if( fabs(b[i]/A[i*n+i]) > 1.0e+100 )
				return i+1;
		
			for(int k = i+1; k < n; k++)
			{
				A[i*n+k] /= A[i*n+i];
				A[k*n+i] /= A[i*n+i];
			}
			for(int k = i+1; k < n; k++)
			for(int q = i+1; q < n; q++)
			{
				A[k*n+q] -= A[k*n+i] * A[i*n+i] * A[i*n+q];
			}
			for(int j = i+1; j < n; j++) //iterate over columns of L
			{
				b[j] -= b[i] * A[j*n+i];
			}
			b[i] /= A[i*n+i];
		}

		for(int i = n-1; i >= 0; i--) //iterate over rows of U
			for(int j = i+1; j < n; j++) 
			{
				b[i] -= b[j] * A[i*n+j];
			}
		for(int i = 0; i < n; i++)
			x[order[i]] = b[i];
	
		return 0;
	}

	class BCGSL_solver : public IterativeMethod
	{
		INMOST_DATA_REAL_TYPE rtol, atol, divtol, last_resid;
		INMOST_DATA_ENUM_TYPE iters, maxits, l, last_it;
		INMOST_DATA_REAL_TYPE resid;
		INMOST_DATA_REAL_TYPE * tau, * sigma, * gamma, *theta1, * theta2, * theta3;
		Solver::Vector r_tilde, x0, t, * u, * r;
		Solver::Matrix * Alink;
		Method * prec;
		std::string reason;
		Solver::OrderInfo * info;
		bool init;
	public:
		INMOST_DATA_ENUM_TYPE GetIterations() {return last_it;}
		INMOST_DATA_REAL_TYPE GetResidual() {return last_resid;}
		INMOST_DATA_REAL_TYPE & RealParameter(std::string name)
		{
			if (name[0] == ':')
			{
				if (prec != NULL) return prec->RealParameter(name.substr(1, name.size() - 1));
			}
			if (name == "rtol") return rtol;
			else if (name == "atol") return atol;
			else if (name == "divtol") return divtol;
			else if (prec != NULL) return prec->RealParameter(name);
			throw - 1;
		}
		INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name)
		{
			if (name[0] == ':')
			{
				if (prec != NULL) return prec->EnumParameter(name.substr(1, name.size() - 1));
			}
			if (name == "maxits") return maxits;
			else if (name == "levels") 
			{
				if( init ) throw - 1; //solver was already initialized, value should not be changed
				return l;
			}
			else if (prec != NULL) return prec->EnumParameter(name);
			throw - 1;
		}
		BCGSL_solver(Method * prec, Solver::OrderInfo & info)
			:rtol(1e-8), atol(1e-9), divtol(1e+40), maxits(1500),l(2),prec(prec),info(&info)
		{
			Alink = NULL;
			init = false;
		}
		bool Initialize()
		{
			if (isInitialized()) Finalize();
			if (prec != NULL && !prec->isInitialized()) prec->Initialize();
			info->PrepareVector(r_tilde);
			info->PrepareVector(x0);
			info->PrepareVector(t);
			tau = new INMOST_DATA_REAL_TYPE[l * 5 + l*l];
			sigma = tau + l*l;
			gamma = sigma + l;
			theta1 = gamma + l;
			theta2 = theta1 + l;
			theta3 = theta2 + l;
			u = new Solver::Vector[l * 2 + 2];
			r = u + l + 1;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < l + 1; k++)
			{
				info->PrepareVector(r[k]);
				info->PrepareVector(u[k]);
			}
			init = true;
			return true;
		}
		bool isInitialized() { return init && (prec == NULL || prec->isInitialized()); }
		bool Finalize()
		{
			if (isInitialized())
			{
				if (!prec->isFinalized()) prec->Finalize();
				delete[] u;
				delete[] tau;
				init = false;
			}
			return true;
		}
		bool isFinalized() { return !init && (prec == NULL || prec->isFinalized()); }
		void Copy(const Method * other)
		{
			const BCGSL_solver * b = dynamic_cast<const BCGSL_solver *>(other);
			assert(b != NULL);
			rtol = b->rtol;
			atol = b->atol;
			divtol = b->divtol;
			last_resid = b->last_resid;
			iters = b->iters;
			maxits = b->maxits;
			l = b->l;
			last_it = b->last_it;
			resid = b->resid;
			Alink = b->Alink;
			info = b->info;
			if (init) Finalize();
			if (b->prec != NULL)
			{
				if (prec == NULL) prec = b->prec->Duplicate();
				else prec->Copy(b->prec);
			}
			if (b->init) Initialize();
		}
		BCGSL_solver(const BCGSL_solver & other) :IterativeMethod(other)
		{
			Copy(&other);
		}
		BCGSL_solver & operator =(BCGSL_solver const & other)
		{
			Copy(&other);
			return *this;
		}
		~BCGSL_solver()
		{
			if (!isFinalized()) Finalize();
			if (prec != NULL) delete prec;
		}
		void ApplyOperator(Solver::Vector & Input, Solver::Vector & Output)
		{
			if (prec != NULL) //right preconditioning here! for left preconditioner have to reverse order
			{
				prec->Solve(Input, t); 
				info->Update(t);
				Alink->MatVec(1.0,t,0,Output);
				info->Update(Output);
			}
			else
			{
				Alink->MatVec(1.0,t,0,Output);
				info->Update(Output);
			}
		}
		bool Solve(Solver::Vector & RHS, Solver::Vector & SOL)
		{
			assert(isInitialized());
			INMOST_DATA_ENUM_TYPE vbeg,vend, vlocbeg, vlocend;
			INMOST_DATA_REAL_TYPE rho0 = 1, rho1, alpha = 0, beta, omega = 1, eta;
			INMOST_DATA_REAL_TYPE resid0, resid, rhs_norm;//, temp[2];
			iters = 0;
			info->PrepareVector(SOL);
			info->PrepareVector(RHS);
			info->Update(SOL);
			info->Update(RHS);
			if( prec != NULL ) prec->ReplaceSOL(SOL);
			if( prec != NULL ) prec->ReplaceRHS(RHS);
			info->GetLocalRegion(info->GetRank(),vlocbeg,vlocend);
			info->GetVectorRegion(vbeg,vend);

			//rhs_norm = info->ScalarProd(RHS,RHS,vlocbeg,vlocend);
			rhs_norm = 1;
			//r[0] = b
			std::copy(RHS.Begin(),RHS.End(),r[0].Begin());
			{
				// r[0] = r[0] - A x
				Alink->MatVec(-1,SOL,1,r[0]); //global multiplication, r probably needs an update
				info->Update(r[0]); // r is good
				std::copy(SOL.Begin(),SOL.End(),x0.Begin()); //x0 = x
				std::fill(SOL.Begin(),SOL.End(),0.0); //x = 0
			}
			std::copy(r[0].Begin(),r[0].End(),r_tilde.Begin()); // r_tilde = r[0]
			std::fill(u[0].Begin(),u[0].End(),0); // u[0] = 0
			resid = info->ScalarProd(r[0],r[0],vlocbeg,vlocend); //resid = dot(r[0],r[0])
			for(INMOST_DATA_ENUM_TYPE k = vbeg; k != vend; k++) // r_tilde = r[0] / dot(r[0],r[0])
				r_tilde[k] /= resid;
			last_resid = resid = resid0 = sqrt(resid/rhs_norm); //resid = sqrt(dot(r[0],r[0])
			last_it = 0;
#if defined(REPORT_RESIDUAL)
			if( info->GetRank() == 0 ) 
			{
				//std::cout << "iter " << last_it << " residual " << resid << std::endl;
				//std::cout << "iter " << last_it << " resid " << resid << "\r";
				//printf("iter %3d resid %12g | %12g relative %12g | %12g\r", last_it, resid, atol, resid / resid0, rtol);
				printf("iter %3d resid %12g | %g\r", last_it, resid, atol);
				fflush(stdout);
			}
#endif
			INMOST_DATA_ENUM_TYPE i = 0;

			if( last_resid < atol || last_resid < rtol*resid0 ) 
			{
				reason = "initial solution satisfy tolerances";
				goto exit;
			}

			long double tt, ts, tp;
			while( true )
			{
				ts = tp = 0;
				rho0 = -omega*rho0;
				
				tt = Timer();
				for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++)
				{
					rho1 = info->ScalarProd(r[j],r_tilde,vlocbeg,vlocend); // rho1 = dot(r[j],r_tilde)
					beta = alpha * (rho1/rho0);

					if( fabs(beta) > 1.0e+100 ) 
					{
						//std::cout << "alpha " << alpha << " rho1 " << rho1 << " rho0 " << rho0 << " beta " << beta << std::endl;
						reason = "multiplier(1) is too large";
						goto exit;
					}

					if( beta != beta )
					{
						reason = "multiplier(1) is NaN";
						goto exit;
					}

					rho0 = rho1;
					for(INMOST_DATA_ENUM_TYPE i = 0; i < j+1; i++)
						for(INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k)
							u[i][k] = r[i][k] - beta*u[i][k];

					ApplyOperator(u[j],u[j+1]); // u[j+1] = A*R*u[j]
					eta = info->ScalarProd(u[j+1],r_tilde,vlocbeg,vlocend); //eta = dot(u[j+1],r_tilde)
					
					alpha = rho0 / eta;

					if( fabs(alpha) > 1.0e+100 ) 
					{
						reason = "multiplier(2) is too large";
						goto exit;
					}
					if( alpha != alpha )
					{
						reason = "multiplier(2) is NaN";
						goto exit;
					}

					for(INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k)
						SOL[k] += alpha*u[0][k];

					for(INMOST_DATA_ENUM_TYPE i = 0; i < j+1; i++)
						for(INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k) //r[i] = r[i] - alpha * u[i+1]
							r[i][k] -= alpha*u[i+1][k];

					
					resid = info->ScalarProd(r[0],r[0],vlocbeg,vlocend); // resid = dot(r[j],r[j])
					resid = sqrt(resid/rhs_norm); // resid = sqrt(dot(r[j],r[j]))

					
					if( resid < atol || resid < rtol*resid0 ) 
					{
						reason = "early exit in bi-cg block";
						last_resid = resid;
						goto exit;
					}
					

					ApplyOperator(r[j],r[j+1]); // r[j+1] = A*R*r[j]

					
				}
				
				for(INMOST_DATA_ENUM_TYPE j = 1; j < l+1; j++)
				{
					for(INMOST_DATA_ENUM_TYPE m = 1; m < j+1; m++)
					{
						tau[(m-1) + (j-1)*l] = 0;
						for(INMOST_DATA_ENUM_TYPE k = vlocbeg; k < vlocend; ++k)
							tau[(m-1) + (j-1)*l] += r[j][k]*r[m][k];
						tau[(j-1) + (m-1)*l] = tau[(m-1) + (j-1)*l];
					}
					sigma[j-1] = 0;
					for(INMOST_DATA_ENUM_TYPE k = vlocbeg; k < vlocend; ++k)
						sigma[j-1] += r[0][k]*r[j][k];
				}
				info->Integrate(tau,l*l+l); //sigma is updated with tau

#if defined(PSEUDOINVERSE)
				{
					int dgesvd_info = 0;


#if defined(USE_LAPACK_SVD)
					char c = 'A';
					INMOST_DATA_REAL_TYPE U[128*128], V[128*128], w[128];
					INMOST_DATA_REAL_TYPE work[5*128];
					int lwork = 5*128;
					int n = l;
					dgesvd_(&c,&c,&n,&n,tau,&n,w,U,&n,V,&n,work,&lwork,&dgesvd_info);
#else
					
					INMOST_DATA_REAL_TYPE U[128*128], V[128*128], S[128*128], w[128];
					svdnxn(tau,U,S,V,l);
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++) w[j] = S[j*l+j];
#endif		
					/*
					printf("w ");
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++) printf("%20g ",w[j]);
					printf("\n");

					printf("U\n");
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l*l; j++) 
					{
						printf("%20g ",U[j]);
						if( (j+1) % l == 0 ) printf("\n");
					}
					printf("\n");

					printf("VT\n");
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l*l; j++) 
					{
						printf("%20g ",V[j]);
						if( (j+1) % l == 0 ) printf("\n");
					}
					printf("\n");
					*/
					if( dgesvd_info != 0 )
					{
						printf("(%s:%d) dgesvd %d\n",__FILE__,__LINE__,dgesvd_info);
						exit(-1);
					}
					
					INMOST_DATA_REAL_TYPE maxw = w[0], tol;
					for(INMOST_DATA_ENUM_TYPE j = 1; j < l; j++) if(w[j]>maxw) maxw = w[j];
					tol = l*maxw*1.0e-14;
					memset(gamma,0,sizeof(INMOST_DATA_REAL_TYPE)*l);
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++)
					{
						if( w[j] > tol )
						{
							INMOST_DATA_REAL_TYPE sum = 0;
							for(INMOST_DATA_ENUM_TYPE k = 0; k < l; ++k)
								sum += sigma[k]*U[j*l+k];
							for(INMOST_DATA_ENUM_TYPE k = 0; k < l; ++k)
								gamma[k] += sum/w[j]*V[k*l+j];
						}
					}
				}

				//svdnxn(tau,U,S,V,l);
				//INMOST_DATA_REAL_TYPE inv_tau[64];
				//pseudoinverse(tau,inv_tau,l);
				//matmul(inv_tau,sigma,gamma,l,l,1);
#else
				int order[128];
				int row = solvenxn(tau,gamma,sigma,l,order);
				if( row != 0 )
				{
					std::cout << "breakdown on row " << row << std::endl;
					reason = "breakdown in matrix inversion in polynomial part";
					break;
				}
#endif
				omega = gamma[l-1];
				if( fabs(omega) > 1.0e+100 )
				{
					reason = "multiplier(3) is too large";
					goto exit;
				}
				if( omega != omega )
				{
					reason = "multiplier(3) is NaN";
					goto exit;
				}
				for(INMOST_DATA_ENUM_TYPE j = 1; j < l+1; ++j)
				{
					for(INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k)
					{
						u[0][k] -= gamma[j-1]*u[j][k];
						SOL[k]  += gamma[j-1]*r[j-1][k];
						r[0][k] -= gamma[j-1]*r[j][k];
					}
				}
				
				
				/*
				for(INMOST_DATA_ENUM_TYPE j = 1; j < l+1; j++)
				{
					for(INMOST_DATA_ENUM_TYPE i = 1; i < j; i++)
					{
						tau[i-1 + (j-1)*l] = 0;
						for(INMOST_DATA_ENUM_TYPE k = vlocbeg; k < vlocend; ++k)
							tau[i-1 + (j-1)*l] += r[j][k]*r[i][k];
						info->Integrate(&tau[i-1 + (j-1)*l],1);
						tau[i-1 + (j-1)*l] /= sigma[i-1];
						for(INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k)
							r[j][k] -= tau[i-1 + (j-1)*l]*r[i][k];
					}
					INMOST_DATA_REAL_TYPE temp[2] = {0,0};
					for(INMOST_DATA_ENUM_TYPE k = vlocbeg; k < vlocend; ++k)
					{
						temp[0] += r[j][k]*r[j][k];
						temp[1] += r[0][k]*r[j][k];
					}
					info->Integrate(temp,2);
					sigma[j-1] = temp[0];//+1.0e-35; //REVIEW
					theta2[j-1] = temp[1]/sigma[j-1];
				}
				omega = theta1[l-1] = theta2[l-1];
				for(INMOST_DATA_ENUM_TYPE j = l-1; j > 0; j--)
				{
					eta = 0;
					for(INMOST_DATA_ENUM_TYPE i = j+1; i < l+1; i++)
						eta += tau[j-1 + (i-1)*l] * theta1[i-1];
					theta1[j-1] = theta2[j-1] - eta;
				}
				for(INMOST_DATA_ENUM_TYPE j = 1; j < l; j++)
				{
					eta = 0;
					for(INMOST_DATA_ENUM_TYPE i = j+1; i < l; i++)
						eta += tau[j-1 + (i-1)*l] * theta1[i];
					theta3[j-1] = theta1[j] + eta;
				}
				for(INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k)
				{
					SOL[k] += theta1[0]*r[0][k];
					r[0][k] -= theta2[l-1]*r[l][k];
					u[0][k] -= theta1[l-1]*u[l][k];
				}
				for(INMOST_DATA_ENUM_TYPE j = 1; j < l; j++)
				{
					for(INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k)
					{
						u[0][k] -= theta1[j-1]*u[j][k];
						SOL[k] += theta3[j-1]*r[j][k];
						r[0][k] -= theta2[j-1]*r[j][k];
					}
				}
				*/
				last_it = i+1;
				{
					resid = info->ScalarProd(r[0],r[0],vlocbeg,vlocend);
					resid = sqrt(resid/rhs_norm);
				}
				tt = Timer() - tt;
#if defined(REPORT_RESIDUAL)
				if( info->GetRank() == 0 ) 
				{
					//std::cout << "iter " << last_it << " residual " << resid << " time " << tt << " matvec " << ts*0.5/l << " precond " << tp*0.5/l << std::endl;
					//std::cout << "iter " << last_it << " resid " << resid << "\r";
					//printf("iter %3d resid %12g | %12g relative %12g | %12g\r", last_it, resid, atol, resid / resid0, rtol);
					printf("iter %3d resid %12g | %g\r", last_it, resid, atol);
					fflush(stdout);
				}
#endif
				last_resid = resid;
				if( resid != resid )
				{
					reason = "residual is NAN";
					break;
				}
				if( resid < atol )
				{
					reason = "converged due to absolute tolerance";
					break;
				}
				if( resid < rtol*resid0 )
				{
					reason = "converged due to relative tolerance";
					break;
				}
				if( resid > divtol )
				{
					reason = "diverged due to divergence tolerance";
					break;
				}
				if( i == maxits )
				{
					reason = "reached maximum iteration number";
					break;
				}
				i++;
			}
exit:
			if (prec != NULL)
			{
				prec->Solve(SOL, r_tilde); //undo right preconditioner
				std::copy(r_tilde.Begin(), r_tilde.End(), SOL.Begin());
			}
			for(INMOST_DATA_ENUM_TYPE k = vlocbeg; k < vlocend; ++k) //undo shift
				SOL[k] += x0[k];
			//info->RestoreMatrix(A);
			info->RestoreVector(SOL);
			info->RestoreVector(RHS);
			if( last_resid < atol || last_resid < rtol*resid0 ) return true;
			return false;
		}
		bool ReplaceMAT(Solver::Matrix & A) { if (isInitialized()) Finalize(); if (prec != NULL) prec->ReplaceMAT(A);  Alink = &A; return true; }
		bool ReplaceRHS(Solver::Vector & RHS) { (void) RHS; return true; }
		bool ReplaceSOL(Solver::Vector & SOL) { (void) SOL; return true; }
		Method * Duplicate() { return new BCGSL_solver(*this);}
		std::string GetReason() {return reason;}
	};


	class BCGS_solver : public IterativeMethod
	{
		INMOST_DATA_REAL_TYPE rtol, atol, divtol, last_resid;
		INMOST_DATA_ENUM_TYPE iters, maxits, last_it;
		INMOST_DATA_REAL_TYPE resid;
		Solver::Vector r0, p, y, s, t, z, r, v;
		Solver::Matrix * Alink;
		Method * prec;
		Solver::OrderInfo * info;
		bool init;
		std::string reason;
	public:
		INMOST_DATA_ENUM_TYPE GetIterations() {return last_it;}
		INMOST_DATA_REAL_TYPE GetResidual() {return last_resid;}
		INMOST_DATA_REAL_TYPE & RealParameter(std::string name)
		{
			if (name[0] == ':')
			{
				if (prec != NULL) return prec->RealParameter(name.substr(1, name.size() - 1));
			}
			if (name == "rtol") return rtol;
			else if (name == "atol") return atol;
			else if (name == "divtol") return divtol;
			else if( prec != NULL ) return prec->RealParameter(name);
			throw - 1;
		}
		INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name)
		{
			if (name[0] == ':')
			{
				if (prec != NULL) return prec->EnumParameter(name.substr(1, name.size() - 1));
			}
			if (name == "maxits") return maxits;
			else if (prec != NULL) return prec->EnumParameter(name);
			throw - 1;
		}
		BCGS_solver(Method * prec, Solver::OrderInfo & info)
			:rtol(1e-8), atol(1e-11), divtol(1e+40), iters(0), maxits(1500),prec(prec),info(&info)
		{
			init = false;
		}
		bool Initialize()
		{
			assert(Alink != NULL);
			if (isInitialized()) Finalize();
			if (prec != NULL && !prec->isInitialized()) prec->Initialize();
			info->PrepareVector(r);
			info->PrepareVector(v);
			info->PrepareVector(p);
			info->PrepareVector(y);
			info->PrepareVector(s);
			info->PrepareVector(t);
			info->PrepareVector(z);
			info->PrepareVector(r0);
			init = true;
			return true;
		}
		bool isInitialized() { return init && (prec == NULL || prec->isInitialized()); }
		bool Finalize()
		{
			if (prec != NULL && !prec->isFinalized()) prec->Finalize();
			init = false;
			return true;
		}
		bool isFinalized() { return !init && (prec == NULL || prec->isFinalized()); }
		void Copy(const Method * other)
		{
			const BCGS_solver * b = dynamic_cast<const BCGS_solver *>(other);
			assert(b != NULL);
			info = b->info;
			rtol = b->rtol;
			atol = b->atol;
			divtol = b->divtol;
			maxits = b->maxits;
			last_resid = b->last_resid;
			iters = b->iters;
			last_it = b->last_it;
			resid = b->resid;
			Alink = b->Alink;
			if (b->prec != NULL)
			{
				if (prec == NULL) prec = b->prec->Duplicate();
				else prec->Copy(b->prec);
			}
			if (b->init) Initialize();
		}
		BCGS_solver(const BCGS_solver & other) : IterativeMethod(other)
		{
			Copy(&other);
		}
		BCGS_solver & operator =(BCGS_solver const & other)
		{
			Copy(&other);
			return *this;
		}
		~BCGS_solver()
		{
			if (!isFinalized()) Finalize();
			if (prec != NULL) delete prec;
		}
		bool Solve(Solver::Vector & RHS, Solver::Vector & SOL)
		{
			assert(isInitialized());
			INMOST_DATA_REAL_TYPE tempa = 0.0, tempb=0.0;
			INMOST_DATA_ENUM_TYPE vbeg,vend, vlocbeg, vlocend;
			INMOST_DATA_INTEGER_TYPE ivbeg,ivend, ivlocbeg, ivlocend;
			INMOST_DATA_REAL_TYPE rho = 1, alpha = 1, beta, omega = 1;
			INMOST_DATA_REAL_TYPE resid0, resid, temp[2];
			bool is_parallel = info->GetSize() > 1;
			info->PrepareVector(SOL);
			info->PrepareVector(RHS);
			if( is_parallel ) info->Update(SOL);
			if( is_parallel ) info->Update(RHS);
			if (prec != NULL)prec->ReplaceSOL(SOL);
			if (prec != NULL)prec->ReplaceRHS(RHS);
			info->GetLocalRegion(info->GetRank(),vlocbeg,vlocend);
			info->GetVectorRegion(vbeg,vend);
			std::copy(RHS.Begin(),RHS.End(),r.Begin());
			{
				Alink->MatVec(-1,SOL,1,r); //global multiplication, r probably needs an update
				info->Update(r); // r is good
			}
			std::copy(r.Begin(),r.End(),r0.Begin());
			std::fill(v.Begin(),v.End(),0.0);
			std::fill(p.Begin(),p.End(),0.0);
			{
				resid = 0;
				for(INMOST_DATA_ENUM_TYPE k = vlocbeg; k != vlocend; k++) 
					resid += r[k]*r[k];
				if( is_parallel ) info->Integrate(&resid,1);
			}
			last_resid = resid = resid0 = sqrt(resid);
			last_it = 0;
			ivbeg = vbeg;
			ivend = vend;
			ivlocbeg = vlocbeg;
			ivlocend = vlocend;
#if defined(REPORT_RESIDUAL)
			if( info->GetRank() == 0 ) 
			{
				//std::cout << "iter " << last_it << " residual " << resid << std::endl;
				//std::cout << "iter " << last_it << " resid " << resid << "\r";
				//printf("iter %3d resid %12g | %12g relative %12g | %12g\r",last_it,resid,atol,resid/resid0,rtol);
				printf("iter %3d resid %12g | %g\r", last_it, resid, atol);
				fflush(stdout);
			}
#endif
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				long double tt, ts, tp, ttt;
				INMOST_DATA_ENUM_TYPE i = 0;
				while(true)
				{
					ts = tp = 0;
					tt = Timer();
					{
						/*
						if( fabs(rho) < 1.0e-31 )
						{
							std::cout << "rho " << rho << " alpha " << alpha << " omega " << omega << " beta " << 1.0 /rho * alpha / omega << std::endl;
							reason = "denominator(1) is zero";
							break;
						}
						if( fabs(omega) < 1.0e-31 )
						{
							std::cout << "rho " << rho << " alpha " << alpha << " omega " << omega << " beta " << 1.0 /rho * alpha / omega << std::endl;
							reason = "denominator(2) is zero";
							break;
						}
						*/
						//std::cout << "rho " << rho << " alpha " << alpha << " omega " << omega << " beta " << 1.0 /rho * alpha / omega << std::endl;
						beta = 1.0 /rho * alpha / omega;
						rho = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:rho)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; k++) 
							rho += r0[k]*r[k];
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Integrate(&rho,1);
						}
						beta *= rho;

						if( fabs(beta) > 1.0e+100 )
						{
							//std::cout << "rho " << rho << " alpha " << alpha << " omega " << omega << " beta " << 1.0 /rho * alpha / omega << std::endl;
							reason = "multiplier(1) is too large";
							break;
						}
						if( beta != beta )
						{
							reason = "multiplier(1) is NaN";
							break;
						}
					}
					{
#if defined(USE_OMP)
#pragma omp for
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k) 
							p[k] = r[k] + beta*(p[k] - omega*v[k]); //global indexes r, p, v
					}
					{
						ttt = Timer();
						if (prec != NULL)
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							prec->Solve(p, y);
						}
						tp += Timer() - ttt;
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Update(y);
						}
						ttt = Timer();
						Alink->MatVec(1,y,0,v); // global multiplication, y should be updated, v probably needs an update
						ts += Timer() - ttt;
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Update(v);
						}
					}
					{
						alpha = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:alpha)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; k++)  
							alpha += r0[k]*v[k];
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Integrate(&alpha,1);
						}

						if( alpha == 0 && rho == 0 ) 
							alpha = 0;
						else
							alpha = rho / alpha; //local indexes, r0, v

						if( fabs(alpha) > 1.0e+100 )
						{
							reason = "multiplier(2) is too large";
							//std::cout << "alpha " << alpha << " rho " << rho << std::endl;
							break;
						}
						if( alpha != alpha )
						{
							reason = "multiplier(2) is NaN";
							//std::cout << "alpha " << alpha << " rho " << rho << std::endl;
							break;
						}
						
					}
					{
#if defined(USE_OMP)
#pragma omp for
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k) 
							s[k] = r[k] - alpha * v[k]; //global indexes r, v
					}
				
					{
						ttt = Timer();
						if (prec != NULL)
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							prec->Solve(s, z);
						}
						tp += Timer() - ttt;
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Update(z);
						}
						ttt = Timer();
						Alink->MatVec(1.0,z,0,t); // global multiplication, z should be updated, t probably needs an update
						ts += Timer() - ttt;
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Update(t);
						}
					}
					{
						
						temp[0] = temp[1] = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:tempa) reduction(+:tempb)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; k++)
						{
							tempa += t[k]*s[k];
							tempb += t[k]*t[k];
						}
						temp[0] = tempa;
						temp[1] = tempb;
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Integrate(temp,2);
						}
						/*
						if (fabs(temp[1]) < 1.0e-35)
						{
							std::cout << "a " << temp[0] << " b " << temp[1] << " omega " << temp[0]/temp[1] << std::endl;
						}
						*/
						//omega = temp[0] / (temp[1] + (temp[1] < 0.0 ? -1.0e-10 : 1.0e-10)); //local indexes t, s
						if( temp[0] == 0 && temp[1] == 0 )
							omega = 0;
						else
							omega = temp[0] / temp[1];

						if( fabs(omega) > 1.0e+100 )
						{
							reason = "multiplier(3) is too large";
							//std::cout << "alpha " << alpha << " rho " << rho << std::endl;
							break;
						}
						if( omega != omega )
						{
							reason = "multiplier(3) is NaN";
							//std::cout << "alpha " << alpha << " rho " << rho << std::endl;
							break;
						}
					}
					{
#if defined(USE_OMP)
#pragma omp for
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k) 
							SOL[k] += alpha * y[k] + omega * z[k]; // global indexes SOL, y, z
					}
					{
#if defined(USE_OMP)
#pragma omp for
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k) 
							r[k] = s[k] - omega * t[k]; // global indexes r, s, t
					}
					last_it = i+1;
					{
						resid = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:resid)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; k++) 
							resid += r[k]*r[k];
						if( is_parallel ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							info->Integrate(&resid,1);
						}
						resid = sqrt(resid);
					}
					tt = Timer() - tt;
#if defined(REPORT_RESIDUAL)
					if( info->GetRank() == 0 ) 
					{
						//std::cout << "iter " << last_it << " residual " << resid << " time " << tt << " matvec " << ts*0.5 << " precond " << tp*0.5 << std::endl;
						//std::cout << "iter " << last_it << " resid " << resid << "\r";
						//printf("iter %3d resid %12g | %12g relative %12g | %12g\r", last_it, resid, atol, resid / resid0, rtol);
#if defined(USE_OMP)
#pragma omp single
#endif
						{
							printf("iter %3d resid %12g | %g\r", last_it, resid, atol);
							fflush(stdout);
						}
					}
#endif
					last_resid = resid;
					if( resid != resid )
					{
						reason = "residual is NAN";
						break;
					}
					if( resid > divtol )
					{
						reason = "diverged due to divergence tolerance";
						break;
					}
					if( resid < atol )
					{
						reason = "converged due to absolute tolerance";
						break;
					}
					if( resid < rtol*resid0 )
					{
						reason = "converged due to relative tolerance";
						break;
					}
					if( i == maxits )
					{
						reason = "reached maximum iteration number";
						break;
					}
					i++;
				}
			}
			//info->RestoreMatrix(A);
			info->RestoreVector(SOL);
			info->RestoreVector(RHS);
			if( last_resid < atol || last_resid < rtol*resid0 ) return true;
			return false;
		}
		bool ReplaceMAT(Solver::Matrix & A) { if (isInitialized()) Finalize();  if (prec != NULL) prec->ReplaceMAT(A);  Alink = &A; return true; }
		bool ReplaceRHS(Solver::Vector & RHS) {(void)RHS; return true; }
		bool ReplaceSOL(Solver::Vector & SOL) {(void)SOL; return true; }
		Method * Duplicate() { return new BCGS_solver(*this);}
		std::string GetReason() {return reason;}
	};
}


#endif
