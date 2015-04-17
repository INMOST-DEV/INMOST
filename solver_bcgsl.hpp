
#ifndef __SOLVER_BCGS__
#define __SOLVER_BCGS__

//\todo
// 1. comply solvers with Method prototype, after TODO in solver_prototypes.hpp is done
// 2. Implement tricks from Read/solver/bcgsl/download.pdf with convex update and true residual correction
// 3. Detect numerical accuracy breakdown - when preconditioned residual is too far from true residual (probably 2 will fix).

#include "inmost_solver.h"

//#define CONVEX_COMBINATION
#define PSEUDOINVERSE  // same trick as in petsc with pseudoinverse
//#define USE_LAPACK_SVD // use lapack's dgesvd routine instead of built-in svdnxn

//#if !defined(NDEBUG)
#define REPORT_RESIDUAL
//#endif
//#define USE_OMP


namespace INMOST
{
	//lapack svd
#if defined(PSEUDOINVERSE)
#if defined(USE_LAPACK_SVD)
	extern "C"
	{
		void dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
	}
#else
	//SVD adopted from http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
	static INMOST_DATA_REAL_TYPE SIGNFUNC(INMOST_DATA_REAL_TYPE a, INMOST_DATA_REAL_TYPE b) { return b >= 0.0 ? fabs(a) : -fabs(a); }
	static INMOST_DATA_REAL_TYPE MAXFUNC(INMOST_DATA_REAL_TYPE x,INMOST_DATA_REAL_TYPE y) { return x > y? x : y; }
	static INMOST_DATA_REAL_TYPE PYTHAG(INMOST_DATA_REAL_TYPE a, INMOST_DATA_REAL_TYPE b)
	{
		INMOST_DATA_REAL_TYPE at = fabs(a), bt = fabs(b), ct, result;
		if (at > bt)       { ct = bt / at; result = sqrt(at) * sqrt(at + ct * bt); }
		else if (bt > 0.0) { ct = at / bt; result = sqrt(bt) * sqrt(bt + ct * at); }
		else result = 0.0;
		return(result);
	}
	int svdnxn(INMOST_DATA_REAL_TYPE * pa,INMOST_DATA_REAL_TYPE * pu, INMOST_DATA_REAL_TYPE *pw, INMOST_DATA_REAL_TYPE * pv, const int n)
	{
		shell<INMOST_DATA_REAL_TYPE> a(pa,n*n);
		shell<INMOST_DATA_REAL_TYPE> u(pu,n*n);
		shell<INMOST_DATA_REAL_TYPE> v(pv,n*n);
		shell<INMOST_DATA_REAL_TYPE> w(pw,n);
		std::copy(a.begin(),a.end(),u.begin());
		//memcpy(u,a,sizeof(INMOST_DATA_REAL_TYPE)*n*n);
		int flag, i, its, j, jj, k, l, nm;
		INMOST_DATA_REAL_TYPE c, f, h, s, x, y, z;
		INMOST_DATA_REAL_TYPE anorm = 0.0, g = 0.0, scale = 0.0;
		static std::vector<INMOST_DATA_REAL_TYPE> rv1(n);
		// Householder reduction to bidiagonal form
		for (i = 0; i < n; i++) 
		{
			// left-hand reduction
			l = i + 1;
			rv1[i] = scale * g;
			g = s = scale = 0.0;
			if (i < n) 
			{
				for (k = i; k < n; k++) 
					scale += fabs(u[k*n+i]);
				if (scale) 
				{
					for (k = i; k < n; k++) 
					{
						u[k*n+i] = u[k*n+i]/scale;
						s += (u[k*n+i] * u[k*n+i]);
					}
					f = u[i*n+i];
					g = -SIGNFUNC(sqrt(s), f);
					h = f * g - s;
					u[i*n+i] = (f - g);
					if (i != n - 1) 
					{
						for (j = l; j < n; j++) 
						{
							for (s = 0.0, k = i; k < n; k++) 
								s += (u[k*n+i] * u[k*n+j]);
							f = s / h;
							for (k = i; k < n; k++) 
								u[k*n+j] += (f * u[k*n+i]);
						}
					}
					for (k = i; k < n; k++) 
						u[k*n+i] = (u[k*n+i]*scale);
				}
			}
			w[i] = (scale * g);

			// right-hand reduction 
			g = s = scale = 0.0;
			if (i < n && i != n - 1) 
			{
				for (k = l; k < n; k++) 
					scale += fabs(u[i*n+k]);
				if (scale) 
				{
					for (k = l; k < n; k++) 
					{
						u[i*n+k] = (u[i*n+k]/scale);
						s += (u[i*n+k] * u[i*n+k]);
					}
					f = u[i*n+l];
					g = -SIGNFUNC(sqrt(s), f);
					h = f * g - s;
					u[i*n+l] = (f - g);
					for (k = l; k < n; k++) 
						rv1[k] = u[i*n+k] / h;
					if (i != n - 1) 
					{
						for (j = l; j < n; j++) 
						{
							for (s = 0.0, k = l; k < n; k++) 
								s += (u[j*n+k] * u[i*n+k]);
							for (k = l; k < n; k++) 
								u[j*n+k] += (s * rv1[k]);
						}
					}
					for (k = l; k < n; k++) 
						u[i*n+k] = (u[i*n+k]*scale);
				}
			}
			anorm = MAXFUNC(anorm, (fabs(w[i]) + fabs(rv1[i])));
		}

		// accumulate the right-hand transformation 
		for (i = n - 1; i >= 0; i--) 
		{
			if (i < n - 1) 
			{
				if (g) 
				{
					for (j = l; j < n; j++)
						v[j*n+i] = ((u[i*n+j] / u[i*n+l]) / g);
					// double division to avoid underflow 
					for (j = l; j < n; j++) 
					{
						for (s = 0.0, k = l; k < n; k++) 
							s += (u[i*n+k] * v[k*n+j]);
						for (k = l; k < n; k++) 
							v[k*n+j] += (s * v[k*n+i]);
					}
				}
				for (j = l; j < n; j++) 
					v[i*n+j] = v[j*n+i] = 0.0;
			}
			v[i*n+i] = 1.0;
			g = rv1[i];
			l = i;
		}

		// accumulate the left-hand transformation
		for (i = n - 1; i >= 0; i--) 
		{
			l = i + 1;
			g = w[i];
			if (i < n - 1) 
				for (j = l; j < n; j++) 
					u[i*n+j] = 0.0;
			if (g) 
			{
				g = 1.0 / g;
				if (i != n - 1) 
				{
					for (j = l; j < n; j++) 
					{
						for (s = 0.0, k = l; k < n; k++) 
							s += (u[k*n+i] * u[k*n+j]);
						f = (s / u[i*n+i]) * g;
						for (k = i; k < n; k++) 
							u[k*n+j] += (f * u[k*n+i]);
					}
				}
				for (j = i; j < n; j++) 
					u[j*n+i] = (u[j*n+i]*g);
			}
			else 
			{
				for (j = i; j < n; j++) 
					u[j*n+i] = 0.0;
			}
			++u[i*n+i];
		}

		// diagonalize the bidiagonal form
		for (k = n - 1; k >= 0; k--) 
		{                             
			// loop over singular values 
			for (its = 0; its < 30; its++) 
			{                         
				// loop over allowed iterations
				flag = 1;
				for (l = k; l >= 0; l--) 
				{                     
					// test for splitting 
					nm = l - 1;
					if (fabs(rv1[l]) + anorm == anorm) 
					{
						flag = 0;
						break;
					}
					if (fabs(w[nm]) + anorm == anorm) 
						break;
				}
				if (flag) 
				{
					c = 0.0;
					s = 1.0;
					for (i = l; i <= k; i++) 
					{
						f = s * rv1[i];
						if (fabs(f) + anorm != anorm) 
						{
							g = w[i];
							h = PYTHAG(f, g);
							w[i] = h; 
							h = 1.0 / h;
							c = g * h;
							s = (- f * h);
							for (j = 0; j < n; j++) 
							{
								y = u[j*n+nm];
								z = u[j*n+i];
								u[j*n+nm] = (y * c + z * s);
								u[j*n+i] = (z * c - y * s);
							}
						}
					}
				}
				z = w[k];
				if (l == k) 
				{                  
					// convergence
					if (z < 0.0) 
					{              
						// make singular value nonnegative
						w[k] = (-z);
						for (j = 0; j < n; j++) 
							v[j*n+k] = (-v[j*n+k]);
					}
					break;
				}
				if (its >= 30) 
				{
					fprintf(stderr, "No convergence after 30,000! iterations \n");
					return 1;
				}

				// shift from bottom 2 x 2 minor
				x = w[l];
				nm = k - 1;
				y = w[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = PYTHAG(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + SIGNFUNC(g, f))) - h)) / x;

				// next QR transformation
				c = s = 1.0;
				for (j = l; j <= nm; j++) 
				{
					i = j + 1;
					g = rv1[i];
					y = w[i];
					h = s * g;
					g = c * g;
					z = PYTHAG(f, h);
					rv1[j] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y = y * c;
					for (jj = 0; jj < n; jj++) 
					{
						x = v[jj*n+j];
						z = v[jj*n+i];
						v[jj*n+j] = (x * c + z * s);
						v[jj*n+i] = (z * c - x * s);
					}
					z = PYTHAG(f, h);
					w[j] = z;
					if (z) 
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = (c * g) + (s * y);
					x = (c * y) - (s * g);
					for (jj = 0; jj < n; jj++) 
					{
						y = u[jj*n+j];
						z = u[jj*n+i];
						u[jj*n+j] = (y * c + z * s);
						u[jj*n+i] = (z * c - y * s);
					}
				}
				rv1[l] = 0.0;
				rv1[k] = f;
				w[k] = x;
			}
		}

		for(int i = 0; i < n; i++)
		{
			INMOST_DATA_REAL_TYPE temp;
			for(int j = i + 1; j < n; j++)
			{
				temp = u[i+n*j];
				u[i+n*j] = u[j+n*i];
				u[j+n*i] = temp;
			}
		}
		for(i = 0; i < n; i++)
		{
			k = i;
			for(j = i+1; j < n; ++j)
				if( w[k] < w[j] ) k = j;
			INMOST_DATA_REAL_TYPE temp;
			if( w[k] > w[i] )
			{
				temp = w[k];
				w[k] = w[i];
				w[i] = temp;
				for(int j = 0; j < n; ++j)
				{
					temp = u[k*n+j];
					u[k*n+j] = u[i*n+j];
					u[i*n+j] = temp;
					temp = v[j*n+k];
					v[j*n+k] = v[j*n+i];
					v[j*n+i] = temp;
				}
			}
		}
		return 0;
	}
#endif //USE_LAPACK_SVD
#else //PSEUDOINVERSE
	static int solvenxn(INMOST_DATA_REAL_TYPE * A, INMOST_DATA_REAL_TYPE * x, INMOST_DATA_REAL_TYPE * b, int n, int * order)
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
#endif //PSEUDOINVERSE
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
			tau = new INMOST_DATA_REAL_TYPE[l * 3 + (l+1)*(l+1) + (l+1)*2];
			sigma = tau + (l+1)*(l+1);
			gamma = sigma + l+1;
			theta1 = gamma + l+1;
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
				Alink->MatVec(1.0,Input,0,Output);
				info->Update(Output);
			}
		}
		bool Solve(Solver::Vector & RHS, Solver::Vector & SOL)
		{
			assert(isInitialized());
			INMOST_DATA_ENUM_TYPE vbeg,vend, vlocbeg, vlocend;
			INMOST_DATA_INTEGER_TYPE ivbeg, ivend, ivlocbeg, ivlocend;
			INMOST_DATA_REAL_TYPE rho0 = 1, rho1, alpha = 0, beta, omega = 1, eta;
			INMOST_DATA_REAL_TYPE resid0, resid, rhs_norm, tau_sum, sigma_sum;//, temp[2];
			iters = 0;
			info->PrepareVector(SOL);
			info->PrepareVector(RHS);
			info->Update(SOL);
			info->Update(RHS);
			if( prec != NULL ) prec->ReplaceSOL(SOL);
			if( prec != NULL ) prec->ReplaceRHS(RHS);
			info->GetLocalRegion(info->GetRank(),vlocbeg,vlocend);
			info->GetVectorRegion(vbeg,vend);
			ivbeg = vbeg;
			ivend = vend;
			ivlocbeg = vlocbeg;
			ivlocend = vlocend;
			//info->ScalarProd(RHS,RHS,vlocbeg,vlocend,rhs_norm);
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
			resid = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:resid)
#endif
			for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k)
				resid += r[0][k]*r[0][k];
			info->Integrate(&resid,1);
			//info->ScalarProd(r[0],r[0],vlocbeg,vlocend,resid); //resid = dot(r[0],r[0])
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k) // r_tilde = r[0] / dot(r[0],r[0])
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
			bool halt = false;
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				long double tt, ts, tp;
				while( !halt )
				{
					ts = tp = 0;

#if defined(USE_OMP)
#pragma omp single
#endif
					{
						rho0 = -omega*rho0;
					}
				
					tt = Timer();
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++)
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						rho1 = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:rho1)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k) 
							rho1+= r[j][k]*r_tilde[k];
						info->Integrate(&rho1,1);
						//info->ScalarProd(r[j],r_tilde,vlocbeg,vlocend,rho1); // rho1 = dot(r[j],r_tilde)
#if defined(USE_OMP)
#pragma omp single
#endif
						{
							beta = alpha * (rho1/rho0);
						}

						if( fabs(beta) > 1.0e+100 ) 
						{
							//std::cout << "alpha " << alpha << " rho1 " << rho1 << " rho0 " << rho0 << " beta " << beta << std::endl;
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "multiplier(1) is too large";
							halt = true;
							break;
						}

						if( beta != beta )
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "multiplier(1) is NaN";
							halt = true;
							break;
						}
#if defined(USE_OMP)
#pragma omp single
#endif
						{
							rho0 = rho1;
						}
						for(INMOST_DATA_ENUM_TYPE i = 0; i < j+1; i++)
						{
#if defined(USE_OMP)
#pragma omp for
#endif
							for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k)
								u[i][k] = r[i][k] - beta*u[i][k];
						}


						ApplyOperator(u[j],u[j+1]); // u[j+1] = A*R*u[j]
#if defined(USE_OMP)
#pragma omp single
#endif
						eta = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:eta)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k) 
							eta += u[j+1][k]*r_tilde[k];
						info->Integrate(&eta,1);
						//info->ScalarProd(u[j+1],r_tilde,vlocbeg,vlocend,eta); //eta = dot(u[j+1],r_tilde)

#if defined(USE_OMP)
#pragma omp single
#endif
						alpha = rho0 / eta;

						if( fabs(alpha) > 1.0e+100 ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "multiplier(2) is too large";
							halt = true;
							break;
						}
						if( alpha != alpha )
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "multiplier(2) is NaN";
							halt = true;
							break;
						}

#if defined(USE_OMP)
#pragma omp for
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k)
							SOL[k] += alpha*u[0][k];

						for(INMOST_DATA_ENUM_TYPE i = 0; i < j+1; i++)
						{
#if defined(USE_OMP)
#pragma omp for
#endif
							for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k) //r[i] = r[i] - alpha * u[i+1]
								r[i][k] -= alpha*u[i+1][k];
						}

					
#if defined(USE_OMP)
#pragma omp single
#endif
						resid = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:resid)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k) 
							resid += r[0][k]*r[0][k];
						info->Integrate(&resid,1);
						//info->ScalarProd(r[0],r[0],vlocbeg,vlocend,resid); // resid = dot(r[j],r[j])
#if defined(USE_OMP)
#pragma omp single
#endif
						resid = sqrt(resid/rhs_norm); // resid = sqrt(dot(r[j],r[j]))

					
						if( resid < atol || resid < rtol*resid0 ) 
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "early exit in bi-cg block";
							last_resid = resid;
							halt = true;
							break;
						}
						ApplyOperator(r[j],r[j+1]); // r[j+1] = A*R*r[j]		
					}

					if( halt ) break;
					INMOST_DATA_ENUM_TYPE size = l;
#if defined(CONVEX_COMBINATION)
					size = l+1;
#endif
					// Penalization for convex combination for update below
					for(INMOST_DATA_ENUM_TYPE j = 1; j < l+1; j++)
					{
						for(INMOST_DATA_ENUM_TYPE m = 1; m < j+1; m++)
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							tau_sum = 0.0;
#if defined(USE_OMP)
#pragma omp for reduction(+:tau_sum)
#endif
							for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k)
								tau_sum += r[j][k]*r[m][k];
#if defined(USE_OMP)
#pragma omp single
#endif
							tau[(j-1) + (m-1)*size] = tau[(m-1) + (j-1)*size] = tau_sum;
						}
#if defined(USE_OMP)
#pragma omp single
#endif
						sigma_sum = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:sigma_sum)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k)
							sigma_sum += r[0][k]*r[j][k];
#if defined(USE_OMP)
#pragma omp single
#endif
						sigma[j-1] = sigma_sum;
					}
#if defined(CONVEX_COMBINATION)
					INMOST_DATA_REAL_TYPE lagrangian = 0.0;
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++) lagrangian += tau[j+size*j];
					sigma[l] = lagrangian;
					tau[(l+1)*(l+1)-1] = 0.0;
					for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++)
					{
						tau[l + j*(l+1)] = -lagrangian;
						tau[l*(l+1)+j] = lagrangian;
					}
#endif
					info->Integrate(tau,size*size+size); //sigma is updated with tau

#if defined(USE_OMP)
#pragma omp single
#endif
#if defined(PSEUDOINVERSE)
					{
						int dgesvd_info = 0;
#if defined(USE_LAPACK_SVD)
						char c = 'A';
						INMOST_DATA_REAL_TYPE U[128*128], V[128*128], w[128];
						INMOST_DATA_REAL_TYPE work[5*128];
						int lwork = 5*128;
						int n = static_cast<int>(size);
						dgesvd_(&c,&c,&n,&n,tau,&n,w,U,&n,V,&n,work,&lwork,&dgesvd_info);
#else
						/*
						char c = 'A';
						INMOST_DATA_REAL_TYPE U2[128*128], V2[128*128], w2[128], tau2[128*128];
						INMOST_DATA_REAL_TYPE work[5*128];
						int lwork = 5*128;
						int n = l;
						memcpy(tau2,tau,sizeof(INMOST_DATA_REAL_TYPE)*l*l);
						dgesvd_(&c,&c,&n,&n,tau2,&n,w2,U2,&n,V2,&n,work,&lwork,&dgesvd_info);
						printf("dgesvd\n");
						printf("w\n");
						for(int q = 0; q < l; ++q) printf("%g ",w2[q]);
						printf("\nU\n");
						for(int q = 0; q < l*l; ++q) 
						{
							printf("%g ",U2[q]);
							if( (q+1)%l == 0 ) printf("\n");
						}
						printf("V\n");
						for(int q = 0; q < l*l; ++q) 
						{
							printf("%g ",V2[q]);
							if( (q+1)%l == 0 ) printf("\n");
						}
						*/
						INMOST_DATA_REAL_TYPE U[128*128], V[128*128], w[128];
						dgesvd_info = svdnxn(tau,U,w,V,size);
						//for(INMOST_DATA_ENUM_TYPE j = 0; j < l; j++) w[j] = S[j*l+j];
						/*
						printf("svdnxn\n");
						printf("w\n");
						for(int q = 0; q < l; ++q) printf("%g ",w[q]);
						printf("\nU\n");
						for(int q = 0; q < l*l; ++q) 
						{
							printf("%g ",U[q]);
							if( (q+1)%l == 0 ) printf("\n");
						}
						printf("V\n");
						for(int q = 0; q < l*l; ++q) 
						{
							printf("%g ",V[q]);
							if( (q+1)%l == 0 ) printf("\n");
						}
						*/
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
						for(INMOST_DATA_ENUM_TYPE j = 1; j < size; j++) if(w[j]>maxw) maxw = w[j];
						tol = size*maxw*1.0e-14;
						memset(gamma,0,sizeof(INMOST_DATA_REAL_TYPE)*size);
						for(INMOST_DATA_ENUM_TYPE j = 0; j < size; j++)
						{
							if( w[j] > tol )
							{
								INMOST_DATA_REAL_TYPE sum = 0;
								for(INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k)
									sum += sigma[k]*U[j*size+k];
								for(INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k)
									gamma[k] += sum/w[j]*V[k*size+j];
							}
						}
					}

					//svdnxn(tau,U,S,V,l);
					//INMOST_DATA_REAL_TYPE inv_tau[64];
					//pseudoinverse(tau,inv_tau,l);
					//matmul(inv_tau,sigma,gamma,l,l,1);
#else
					{
						int order[128];
						int row = solvenxn(tau,gamma,sigma,size,order);
						/*
						double sum = 0.0;
						for(int j = 0; j < l; ++j) 
						{
							sum += gamma[j];
							std::cout << gamma[j] << " ";
						}
						std::cout << "sum: " << sum;
						//std::cout << " lagrangian: " << gamma[l];
						std::cout << std::endl;
						*/
						if( row != 0 )
						{
							std::cout << "breakdown on row " << row << std::endl;
							reason = "breakdown in matrix inversion in polynomial part";
							break;
						}
					}
#endif
#if defined(USE_OMP)
#pragma omp single
#endif
					omega = gamma[l-1];
					if( fabs(omega) > 1.0e+100 )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "multiplier(3) is too large";
						halt = true;
						break;
					}
					if( omega != omega )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "multiplier(3) is NaN";
						halt = true;
						break;
					}
					for(INMOST_DATA_ENUM_TYPE j = 1; j < l+1; ++j)
					{
#if defined(USE_OMP)
#pragma omp for
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivbeg; k < ivend; ++k)
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
#if defined(USE_OMP)
#pragma omp single
#endif
						resid = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:resid)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k)
							resid += r[0][k]*r[0][k];
						info->Integrate(&resid,1);
						//info->ScalarProd(r[0],r[0],vlocbeg,vlocend,resid);
#if defined(USE_OMP)
#pragma omp single
#endif
						resid = sqrt(resid/rhs_norm);
					}
					tt = Timer() - tt;
#if defined(REPORT_RESIDUAL)
					if( info->GetRank() == 0 ) 
					{
						//std::cout << "iter " << last_it << " residual " << resid << " time " << tt << " matvec " << ts*0.5/l << " precond " << tp*0.5/l << std::endl;
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
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "residual is NAN";
						break;
					}
					if( resid < atol )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "converged due to absolute tolerance";
						break;
					}
					if( resid < rtol*resid0 )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "converged due to relative tolerance";
						break;
					}
					if( resid > divtol )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "diverged due to divergence tolerance";
						break;
					}
					if( i == maxits )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "reached maximum iteration number";
						break;
					}
					i++;
				}
			}
exit:
			if (prec != NULL)
			{
				prec->Solve(SOL, r_tilde); //undo right preconditioner
				std::copy(r_tilde.Begin(), r_tilde.End(), SOL.Begin());
			}
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k) //undo shift
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
			info->PrepareVector(SOL);
			info->PrepareVector(RHS);
			info->Update(SOL);
			info->Update(RHS);
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; k++) 
					resid += r[k]*r[k];
				info->Integrate(&resid,1);
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
				INMOST_DATA_ENUM_TYPE i = 0;
				while(true)
				{
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
#if defined(USE_OMP)
#pragma omp single
#endif
						{
							beta = 1.0 /rho * alpha / omega;
							rho = 0;
						}
#if defined(USE_OMP)
#pragma omp for reduction(+:rho)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k) 
							rho += r0[k]*r[k];
						info->Integrate(&rho,1);
						//info->ScalarProd(r0,r,ivlocbeg,ivlocend,rho);
#if defined(USE_OMP)
#pragma omp single
#endif
						{
							beta *= rho;
						}

						if( fabs(beta) > 1.0e+100 )
						{
							//std::cout << "rho " << rho << " alpha " << alpha << " omega " << omega << " beta " << 1.0 /rho * alpha / omega << std::endl;
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "multiplier(1) is too large";
							break;
						}
						if( beta != beta )
						{
#if defined(USE_OMP)
#pragma omp single
#endif
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
						if (prec != NULL) 
						{
							prec->Solve(p, y);
							info->Update(y);
							Alink->MatVec(1,y,0,v); // global multiplication, y should be updated, v probably needs an update
							info->Update(v);
						}
						else
						{
							Alink->MatVec(1,p,0,v); // global multiplication, y should be updated, v probably needs an update
							info->Update(v);
						}
					}
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						alpha = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:alpha)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k) 
							alpha += r0[k]*v[k];
						info->Integrate(&alpha,1);
						//info->ScalarProd(r0,v,ivlocbeg,ivlocend,alpha);
#if defined(USE_OMP)
#pragma omp single
#endif
						{
							if( alpha == 0 && rho == 0 ) 
								alpha = 0;
							else
								alpha = rho / alpha; //local indexes, r0, v
						}

						if( fabs(alpha) > 1.0e+100 )
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "multiplier(2) is too large";
							//std::cout << "alpha " << alpha << " rho " << rho << std::endl;
							break;
						}
						if( alpha != alpha )
						{
#if defined(USE_OMP)
#pragma omp single
#endif
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
						if (prec != NULL) 
						{
							prec->Solve(s, z);
							info->Update(z);
							Alink->MatVec(1.0,z,0,t); // global multiplication, z should be updated, t probably needs an update
							info->Update(t);
						}
						else
						{
							Alink->MatVec(1.0,s,0,t); // global multiplication, z should be updated, t probably needs an update
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
						info->Integrate(temp,2);
						/*
						if (fabs(temp[1]) < 1.0e-35)
						{
							std::cout << "a " << temp[0] << " b " << temp[1] << " omega " << temp[0]/temp[1] << std::endl;
						}
						*/
						//omega = temp[0] / (temp[1] + (temp[1] < 0.0 ? -1.0e-10 : 1.0e-10)); //local indexes t, s
#if defined(USE_OMP)
#pragma omp single
#endif
						{
							if( temp[0] == 0 && temp[1] == 0 )
								omega = 0;
							else
								omega = temp[0] / temp[1];
						}

						if( fabs(omega) > 1.0e+100 )
						{
#if defined(USE_OMP)
#pragma omp single
#endif
							reason = "multiplier(3) is too large";
							//std::cout << "alpha " << alpha << " rho " << rho << std::endl;
							break;
						}
						if( omega != omega )
						{
#if defined(USE_OMP)
#pragma omp single
#endif
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
					//info->ScalarProd(r,r,ivlocbeg,ivlocend,resid);
#if defined(USE_OMP)
#pragma omp single
#endif
					resid = 0;
#if defined(USE_OMP)
#pragma omp for reduction(+:resid)
#endif
					for(INMOST_DATA_INTEGER_TYPE k = ivlocbeg; k < ivlocend; ++k)
						resid += r[k]*r[k];
#if defined(USE_OMP)
#pragma omp single
#endif
					{
						info->Integrate(&resid,1);
						resid = sqrt(resid);
					}
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
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "residual is NAN";
						break;
					}
					if( resid > divtol )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "diverged due to divergence tolerance";
						break;
					}
					if( resid < atol )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "converged due to absolute tolerance";
						break;
					}
					if( resid < rtol*resid0 )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
						reason = "converged due to relative tolerance";
						break;
					}
					if( i == maxits )
					{
#if defined(USE_OMP)
#pragma omp single
#endif
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
