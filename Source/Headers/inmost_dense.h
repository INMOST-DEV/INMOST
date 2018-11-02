#ifndef INMOST_DENSE_INCLUDED
#define INMOST_DENSE_INCLUDED
#include "inmost_common.h"
#include "inmost_expression.h"
#include <iomanip>

#define pool_array_t pool_array
//#define pool_array_t array

namespace INMOST
{
	
	/// Structure that selects desired class, depending on the operation.
	template<class A, class B> struct Promote;
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_INTEGER_TYPE> {typedef INMOST_DATA_INTEGER_TYPE type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, INMOST_DATA_INTEGER_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
#if defined(USE_AUTODIFF)
	//For INMOST_DATA_INTEGER_TYPE
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, unknown>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, variable>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, multivar_expression_reference>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, hessian_variable>  {typedef hessian_variable type;};
	//For INMOST_DATA_REAL_TYPE
	template<> struct Promote<INMOST_DATA_REAL_TYPE, unknown>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, variable>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, multivar_expression_reference>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, hessian_variable>  {typedef hessian_variable type;};
	//for unknown
	//For INMOST_DATA_INTEGER_TYPE
	template<> struct Promote<unknown, INMOST_DATA_INTEGER_TYPE>  {typedef variable type;};
	template<> struct Promote<unknown, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
	template<> struct Promote<unknown, unknown>  {typedef variable type;};
	template<> struct Promote<unknown, variable>  {typedef variable type;};
	template<> struct Promote<unknown, multivar_expression_reference>  {typedef variable type;};
	template<> struct Promote<unknown, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<unknown, hessian_variable>  {typedef hessian_variable type;};
	//For multivar_expression_reference
	template<> struct Promote<multivar_expression_reference, INMOST_DATA_INTEGER_TYPE>  {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, unknown> {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, variable> {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, multivar_expression_reference> {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, hessian_multivar_expression_reference> {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, hessian_variable> {typedef variable type;};
	//For variable
	template<> struct Promote<variable, INMOST_DATA_INTEGER_TYPE>  {typedef variable type;};
	template<> struct Promote<variable, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
	template<> struct Promote<variable, unknown> {typedef variable type;};
	template<> struct Promote<variable, variable> {typedef variable type;};
	template<> struct Promote<variable, multivar_expression_reference> {typedef variable type;};
	template<> struct Promote<variable, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<variable, hessian_variable>  {typedef hessian_variable type;};
	//For hessian_multivar_expression_reference
	template<> struct Promote<hessian_multivar_expression_reference, INMOST_DATA_INTEGER_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, INMOST_DATA_REAL_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, unknown>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, variable>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, hessian_multivar_expression_reference> {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, hessian_variable> {typedef hessian_variable type;};
	//For hessian_variable
	template<> struct Promote<hessian_variable, INMOST_DATA_INTEGER_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, INMOST_DATA_REAL_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, unknown>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, variable>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, hessian_multivar_expression_reference> {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, hessian_variable> {typedef hessian_variable type;};
#endif
	
	
	
	/// Abstract class for a matrix,
	/// used to abstract away all the data storage and access
	/// and provide common implementation of the algorithms.
	template<typename Var>
	class AbstractMatrix
	{
	private:
		/// Sign function for SVD algorithm.
		static Var sign_func(const Var & a, const Var & b) {return (b >= 0.0 ? fabs(a) : -fabs(a));}
		/// Max function for SVD algorithm.
		static INMOST_DATA_REAL_TYPE max_func(INMOST_DATA_REAL_TYPE x, INMOST_DATA_REAL_TYPE y) { return x > y ? x : y; }
		/// Function for QR rotation in SVD algorithm.
		static Var pythag(const Var & a, const Var & b)
		{
			Var at = fabs(a), bt = fabs(b), ct, result;
			if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
			else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
			else result = 0.0;
			return result;
		}
	public:
		typedef unsigned enumerator;
		/// Construct empty matrix.
		AbstractMatrix() {}
		/// Construct from matrix of the same type.
		/// @param other Another matrix of the same type.
		AbstractMatrix(const AbstractMatrix & b)
		{
			Resize(b.Rows(),b.Cols());
			for(enumerator i = 0; i < b.Rows(); ++i)
				for(enumerator j = 0; j < b.Cols(); ++j)
					(*this)(i,j) = b(i,j);
		}
		/// Construct from matrix of another type.
		/// @param other Another matrix of different type.
		template<typename typeB>
		AbstractMatrix(const AbstractMatrix<typeB> & b)
		{
			Resize(b.Rows(),b.Cols());
			for(enumerator i = 0; i < b.Rows(); ++i)
				for(enumerator j = 0; j < b.Cols(); ++j)
					assign((*this)(i,j),b(i,j));
		}
		/// Assign matrix of the same type.
		/// @param other Another matrix of the same type.
		/// @return Reference to matrix.
		AbstractMatrix & operator =(AbstractMatrix const & other)
		{
			Resize(other.Rows(),other.Cols());
			for(enumerator i = 0; i < other.Rows(); ++i)
				for(enumerator j = 0; j < other.Cols(); ++j)
					assign((*this)(i,j),other(i,j));
			return *this;
		}
		/// Assign matrix of another type.
		/// @param other Another matrix of different type.
		/// @return Reference to matrix.
		template<typename typeB>
		AbstractMatrix & operator =(AbstractMatrix<typeB> const & other)
		{
			Resize(other.Rows(),other.Cols());
			for(enumerator i = 0; i < other.Rows(); ++i)
				for(enumerator j = 0; j < other.Cols(); ++j)
					assign((*this)(i,j),other(i,j));
			return *this;
		}
		/// Assign value to all entries of the matrix.
		/// @param b Assigned value.
		/// @return Reference to matrix.
		AbstractMatrix & operator =(Var const & b)
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					assign((*this)(i,j),b);
			return *this;
		}
		/// Obtain number of rows.
		/// @return Number of rows.
		virtual enumerator Rows() const = 0;
		/// Obtain number of columns.
		/// @return Number of columns.
		virtual enumerator Cols() const = 0;
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		virtual Var & operator () (enumerator i, enumerator j) = 0;
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		virtual const Var & operator () (enumerator i, enumerator j) const = 0;
		/// Resize the matrix into different size.
		/// @param nrows New number of rows.
		/// @param ncols New number of columns.
		virtual void Resize(enumerator rows, enumerator cols) = 0;
		/// Check all matrix entries for nans.
		/// Also checks derivatives for matrices of variables.
		bool CheckNans() const
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					if( check_nans((*this)(i,j)) ) return true;
			return false;
		}
		bool CheckInfs() const
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					if( check_infs((*this)(i,j)) ) return true;
			return false;
		}
		bool CheckNansInfs() const
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					if( check_nans_infs((*this)(i,j)) ) return true;
			return false;
		}
		/// Maximum product transversal.
		/// Computes unsymmetric reordering that maximizes product on diagonal.
		/// Returns reordering matrix P and scaling matrix S that transforms matrix into I-dominant matrix.
		/// @param Perm Array for reordering, size of columns of the matrix.
		/// @param SL Diagonal for rescaling matrix from left, size of columns of the matrix.
		/// @param SR Diagonal for rescaling matrix from right, size of rows of the matrix.
		/// \todo 
		/// 1. Test rescaling.
		/// 2. Test on non-square matrices.
		void MPT(INMOST_DATA_ENUM_TYPE * Perm, INMOST_DATA_REAL_TYPE * SL = NULL, INMOST_DATA_REAL_TYPE * SR = NULL) const;
		/// Singular value decomposition.
		/// Reconstruct matrix: A = U*Sigma*V.Transpose().
		/// Source http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
		/// With few corrections for robustness.
		/// @param U Left unitary matrix, U^T U = I.
		/// @param Sigma Diagonal matrix with singular values.
		/// @param V Right unitary matrix, not transposed.
		/// @param order_singular_values Return singular values in descending order.
		/// @param nonnegative Change sign of singular values.
		/// \warning Somehow result differ in auto-differential mode.
		/// \todo Test implementation for auto-differentiation.
		bool SVD(AbstractMatrix<Var> & U, AbstractMatrix<Var> & Sigma, AbstractMatrix<Var> & V, bool order_singular_values = true, bool nonnegative = true) const
		{
			int flag, i, its, j, jj, k, l, nm;
			int n = Rows();
			int m = Cols();
			Var c, f, h, s, x, y, z;
			Var g = 0.0, scale = 0.0;
			INMOST_DATA_REAL_TYPE anorm = 0.0;
			if (n >= m)
			{
				U = (*this);
				Sigma.Resize(m,m);
				Sigma.Zero();
				V.Resize(m,m);
			}
			else // n < m
			{
				U.Resize(n,n);
				Sigma.Resize(n,n);
				Sigma.Zero();
				V.Resize(m,n);
			}
			if (n < m)
			{
				bool success = Transpose().SVD(V,Sigma,U);
				if( success )
				{
					//U.Swap(V);
					//U = U.Transpose();
					//V = V.Transpose();
					return true;
				}
				else return false;
			} //m <= n
			dynarray<Var,128> rv1(m);
			//array<Var> _rv1(m);
			//shell<Var> rv1(_rv1);
			std::swap(n,m); //this how original algorithm takes it
			// Householder reduction to bidiagonal form
			for (i = 0; i < n; i++)
			{
				// left-hand reduction
				l = i + 1;
				rv1[i] = scale * g;
				g = s = scale = 0.0;
				if (i < m)
				{
					for (k = i; k < m; k++) scale += fabs(U(k,i));
					if (get_value(scale))
					{
						for (k = i; k < m; k++)
						{
							U(k,i) /= scale;
							s += U(k,i) * U(k,i);
						}
						f = U(i,i);
						g = -sign_func(sqrt(s), f);
						h = f * g - s;
						U(i,i) = f - g;
						if (i != n - 1)
						{
							for (j = l; j < n; j++)
							{
								for (s = 0.0, k = i; k < m; k++) s += U(k,i) * U(k,j);
								f = s / h;
								for (k = i; k < m; k++) U(k,j) += f * U(k,i);
							}
						}
						for (k = i; k < m; k++) U(k,i) *= scale;
					}
				}
				Sigma(i,i) = scale * g;
				// right-hand reduction
				g = s = scale = 0.0;
				if (i < m && i != n - 1)
				{
					for (k = l; k < n; k++) scale += fabs(U(i,k));
					if (get_value(scale))
					{
						for (k = l; k < n; k++)
						{
							U(i,k) = U(i,k)/scale;
							s += U(i,k) * U(i,k);
						}
						f = U(i,l);
						g = -sign_func(sqrt(s), f);
						h = f * g - s;
						U(i,l) = f - g;
						for (k = l; k < n; k++) rv1[k] = U(i,k) / h;
						if (i != m - 1)
						{
							for (j = l; j < m; j++)
							{
								for (s = 0.0, k = l; k < n; k++) s += U(j,k) * U(i,k);
								for (k = l; k < n; k++) U(j,k) += s * rv1[k];
							}
						}
						for (k = l; k < n; k++) U(i,k) *= scale;
					}
				}
				anorm = max_func(anorm,fabs(get_value(Sigma(i,i))) + fabs(get_value(rv1[i])));
			}
			
			// accumulate the right-hand transformation
			for (i = n - 1; i >= 0; i--)
			{
				if (i < (n - 1))
				{
					if (get_value(g))
					{
						for (j = l; j < n; j++) V(j,i) = ((U(i,j) / U(i,l)) / g);
						// double division to avoid underflow
						for (j = l; j < n; j++)
						{
							for (s = 0.0, k = l; k < n; k++) s += U(i,k) * V(k,j);
							for (k = l; k < n; k++) V(k,j) += s * V(k,i);
						}
					}
					for (j = l; j < n; j++) V(i,j) = V(j,i) = 0.0;
				}
				V(i,i) = 1.0;
				g = rv1[i];
				l = i;
			}
			
			// accumulate the left-hand transformation
			for (i = n - 1; i >= 0; i--)
			{
				l = i + 1;
				g = Sigma(i,i);
				if (i < (n - 1))
					for (j = l; j < n; j++)
						U(i,j) = 0.0;
				if (get_value(g))
				{
					g = 1.0 / g;
					if (i != n - 1)
					{
						for (j = l; j < n; j++)
						{
							for (s = 0.0, k = l; k < m; k++) s += (U(k,i) * U(k,j));
							f = (s / U(i,i)) * g;
							for (k = i; k < m; k++) U(k,j) += f * U(k,i);
						}
					}
					for (j = i; j < m; j++) U(j,i) = U(j,i)*g;
				}
				else for (j = i; j < m; j++) U(j,i) = 0.0;
				U(i,i) += 1;
			}
			
			// diagonalize the bidiagonal form
			for (k = n - 1; k >= 0; k--)
			{// loop over singular values
				for (its = 0; its < 30; its++)
				{// loop over allowed iterations
					flag = 1;
					for (l = k; l >= 0; l--)
					{// test for splitting
						nm = l - 1;
						if (fabs(get_value(rv1[l])) + anorm == anorm)
						{
							flag = 0;
							break;
						}
						if( l == 0 ) break;
						if (fabs(get_value(Sigma(nm,nm))) + anorm == anorm)
							break;
					}
					if (flag && l != 0)
					{
						c = 0.0;
						s = 1.0;
						for (i = l; i <= k; i++)
						{
							f = s * rv1[i];
							if (fabs(get_value(f)) + anorm != anorm)
							{
								g = Sigma(i,i);
								h = pythag(f, g);
								Sigma(i,i) = h;
								h = 1.0 / h;
								c = g * h;
								s = (- f * h);
								for (j = 0; j < m; j++)
								{
									y = U(j,nm);
									z = U(j,i);
									U(j,nm) = (y * c + z * s);
									U(j,i) = (z * c - y * s);
								}
							}
						}
					}
					z = Sigma(k,k);
					if (l == k)
					{// convergence
						if (z < 0.0 && nonnegative)
						{// make singular value nonnegative
							Sigma(k,k) = -z;
							for (j = 0; j < n; j++) V(j,k) = -V(j,k);
						}
						break;
					}
					if (its >= 30)
					{
						std::cout << "No convergence after " << its << " iterations" << std::endl;
						std::swap(n,m);
						return false;
					}
					// shift from bottom 2 x 2 minor
					x = Sigma(l,l);
					nm = k - 1;
					y = Sigma(nm,nm);
					g = rv1[nm];
					h = rv1[k];
					f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
					g = pythag(f, 1.0);
					f = ((x - z) * (x + z) + h * ((y / (f + sign_func(g, f))) - h)) / x;
					// next QR transformation
					c = s = 1.0;
					for (j = l; j <= nm; j++)
					{
						i = j + 1;
						g = rv1[i];
						y = Sigma(i,i);
						h = s * g;
						g = c * g;
						z = pythag(f, h);
						rv1[j] = z;
						c = f / z;
						s = h / z;
						f = x * c + g * s;
						g = g * c - x * s;
						h = y * s;
						y = y * c;
						for (jj = 0; jj < n; jj++)
						{
							x = V(jj,j);
							z = V(jj,i);
							V(jj,j) = (x * c + z * s);
							V(jj,i) = (z * c - x * s);
						}
						z = pythag(f, h);
						Sigma(j,j) = z;
						if (get_value(z))
						{
							z = 1.0 / z;
							c = f * z;
							s = h * z;
						}
						f = (c * g) + (s * y);
						x = (c * y) - (s * g);
						for (jj = 0; jj < m; jj++)
						{
							y = U(jj,j);
							z = U(jj,i);
							U(jj,j) = (y * c + z * s);
							U(jj,i) = (z * c - y * s);
						}
					}
					rv1[l] = 0.0;
					rv1[k] = f;
					Sigma(k,k) = x;
				}
			}
			//CHECK THIS!
			if( order_singular_values )
			{
				Var temp;
				for(i = 0; i < n; i++)
				{
					k = i;
					for(j = i+1; j < n; ++j)
						if( Sigma(k,k) < Sigma(j,j) ) k = j;
					if( Sigma(k,k) > Sigma(i,i) )
					{
						temp       = Sigma(k,k);
						Sigma(k,k) = Sigma(i,i);
						Sigma(i,i) = temp;
						// U is m by n
						for(int j = 0; j < m; ++j)
						{
							temp   = U(j,k);
							U(j,k) = U(j,i);
							U(j,i) = temp;
						}
						// V is n by n
						for(int j = 0; j < n; ++j)
						{
							temp   = V(j,k);
							V(j,k) = V(j,i);
							V(j,i) = temp;
						}
					}
				}
			}
			
			std::swap(n,m);
			return true;
		}
		/// Set all the elements of the matrix to zero.
		void Zero()
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					(*this)(i,j) = 0.0;
		}
		/// Transpose the current matrix.
		/// @return Transposed matrix.
		Matrix<Var, pool_array_t<Var> > Transpose() const;
		///Exchange contents of two matrices.
		virtual void Swap(AbstractMatrix<Var> & b);
		/// Unary minus. Change sign of each element of the matrix.
		Matrix<Var, pool_array_t<Var> > operator-() const;
		/// Cross-product operation for a vector.
		/// Both right hand side and left hand side should be a vector
		/// @param other The right hand side of cross product.
		/// @return The cross product of current and right hand side vector.
		/// \warning Works for 3x1 vector and 3xm m-vectors as right hand side.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		CrossProduct(const AbstractMatrix<typeB> & other) const;
		/// Transformation matrix from current vector to provided vector.
		/// @param other Vector to transform to.
		/// @return A sqaure (rotation) matrix that transforms current vector into right hand side vector.
		/// \warning Works only for 3x1 vectors.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		Transform(const AbstractMatrix<typeB> & other) const;
		/// Subtract a matrix.
		/// @param other The matrix to be subtracted.
		/// @return Difference of current matrix with another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		operator-(const AbstractMatrix<typeB> & other) const;
		/// Subtract a matrix and store result in the current.
		/// @param other The matrix to be subtracted.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator-=(const AbstractMatrix<typeB> & other);
		/// Add a matrix.
		/// @param other The matrix to be added.
		/// @return Sum of current matrix with another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		operator+(const AbstractMatrix<typeB> & other) const;
		/// Add a matrix and store result in the current.
		/// @param other The matrix to be added.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator+=(const AbstractMatrix<typeB> & other);
		/// Multiply the matrix by another matrix.
		/// @param other Matrix to be multiplied from right.
		/// @return Matrix multipled by another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		operator*(const AbstractMatrix<typeB> & other) const;
		/// Multiply matrix with another matrix in-place.
		/// @param B Another matrix to the right in multiplication.
		/// @return Reference to current matrix.
		template<typename typeB>
		AbstractMatrix & operator*=(const AbstractMatrix<typeB> & B);
		/// Multiply the matrix by a coefficient.
		/// @param coef Coefficient.
		/// @return Matrix multiplied by the coefficient.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		operator*(typeB coef) const;
		/// Multiply the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator*=(typeB coef);
		/// Divide the matrix by a coefficient of a different type.
		/// @param coef Coefficient.
		/// @return Matrix divided by the coefficient.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		operator/(typeB coef) const;
		/// Divide the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix &
		operator/=(typeB coef);
		/// Performs B^{-1}*A, multiplication by inverse matrix from left.
		/// Throws exception if matrix is not invertable. See Mesh::PseudoSolve for
		/// singular matrices.
		/// @param other Matrix to be inverted and multiplied from left.
		/// @return Multiplication of current matrix by inverse of another
		/// @see Matrix::PseudoInvert.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		operator/(const AbstractMatrix<typeB> & other) const
		{
			Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(other.Cols(),Cols());
			ret = other.Solve(*this);
			return ret;
		}
		/// Kronecker product, latex symbol \otimes.
		/// @param other Matrix on the right of the Kronecker product.
		/// @return Result of the Kronecker product of current and another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		Kronecker(const AbstractMatrix<typeB> & other) const;
		/// Inverts matrix using Crout-LU decomposition with full pivoting for
		/// maximum element. Works well for non-singular matrices, for singular
		/// matrices look see Matrix::PseudoInvert.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr equal to a positive integer that represents the row
		///             on which the failure occured (starting from 1), in case of no failure *ierr = 0.
		/// @return Inverse matrix.
		/// @see Matrix::PseudoInvert.
		/// \todo (test) Activate and test implementation with Solve.
		Matrix<Var, pool_array_t<Var> > Invert(int * ierr = NULL) const;
		/// Inverts symmetric positive-definite matrix using Cholesky decomposition.
		Matrix<Var, pool_array_t<Var> > CholeskyInvert(int * ierr = NULL) const;
		/// Finds X in A*X=B, where A and B are general matrices.
		/// Converts system into A^T*A*X=A^T*B.
		/// Inverts matrix A^T*A using Crout-LU decomposition with full pivoting for
		/// maximum element. Works well for non-singular matrices, for singular
		/// matrices see Matrix::PseudoInvert.
		/// @param B Right hand side matrix.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr equal to a positive integer that represents the row
		///             on which the failure occured (starting from 1), in case of no failure *ierr = 0.
		/// @return Inverse matrix,
		/// @see Matrix::PseudoInvert.
		/// \warning Number of rows in matrices A and B should match.
		/// \todo
		/// 1. Test implementation.
		/// 2. Maximum product transversal + block pivoting instead of pivoting
		///    by maximum element.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		Solve(const AbstractMatrix<typeB> & B, int * ierr = NULL) const;
		/// Finds X in A*X=B, where A is a square symmetric positive definite matrix.
		/// Uses Cholesky decomposition algorithm.
		/// @param B Right hand side matrix.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr equal to a positive integer that represents the row
		///             on which the failure occured (starting from 1), in case of no failure *ierr = 0.
		/// @see Matrix::PseudoInvert.
		/// @return Inverse matrix,
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		CholeskySolve(const AbstractMatrix<typeB> & B, int * ierr = NULL) const;
		/// Calculate sum of the diagonal elements of the matrix.
		/// @return Trace of the matrix.
		Var Trace() const
		{
			assert(Cols() == Rows());
			Var ret = 0.0;
			for(enumerator i = 0; i < Rows(); ++i) ret += (*this)(i,i);
			return ret;
		}
		/// Output matrix to screen.
		/// Does not output derivatices.
		/// @param threshold Elements smaller then the number are considered to be zero.
		void Print(INMOST_DATA_REAL_TYPE threshold = 1.0e-10) const
		{
			for(enumerator k = 0; k < Rows(); ++k)
			{
				for(enumerator l = 0; l < Cols(); ++l)
				{
					if( __isinf__(get_value((*this)(k,l))) )
						std::cout << std::setw(12) << "inf";
					else if( std::isnan(get_value((*this)(k,l))) )
						std::cout << std::setw(12) << "nan";
					else if( fabs(get_value((*this)(k,l))) > threshold )
						std::cout << std::setw(12) << get_value((*this)(k,l));
					else
						std::cout << std::setw(12) << 0;
					std::cout << " ";
				}
				std::cout << std::endl;
			}
		}
		/// Check if the matrix is symmetric.
		/// @return Returns true if the matrix is symmetric, otherwise false.
		bool isSymmetric() const
		{
			if( Rows() != Cols() ) return false;
			for(enumerator k = 0; k < Rows(); ++k)
			{
				for(enumerator l = k+1; l < Rows(); ++l)
					if( fabs((*this)(k,l)-(*this)(l,k)) > 1.0e-7 )
						return false;
			}
			return true;
		}
		/// Computes dot product by summing up multiplication of entries with the
		/// same indices in the current and the provided matrix.
		/// @param other Provided matrix.
		/// @return Dot product of two matrices.
		template<typename typeB>
		typename Promote<Var,typeB>::type
		DotProduct(const AbstractMatrix<typeB> & other) const
		{
			assert(Cols() == other.Cols());
			assert(Rows() == other.Rows());
			typename Promote<Var,typeB>::type ret = 0.0;
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret += ((*this)(i,j))*other(i,j);
			return ret;
		}
		/// Computes dot product by summing up multiplication of entries with the
		/// same indices in the current and the provided matrix.
		/// @param other Provided matrix.
		/// @return Dot product of two matrices.
		template<typename typeB>
		typename Promote<Var,typeB>::type operator ^(const AbstractMatrix<typeB> & other) const
		{
			return DotProduct(other);
		}
		/// Computes frobenious norm of the matrix.
		/// @return Frobenius norm of the matrix.
		Var FrobeniusNorm() const
		{
			Var ret = 0;
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret += (*this)(i,j)*(*this)(i,j);
			return sqrt(ret);
		}
		/// Calculates Moore-Penrose pseudo-inverse of the matrix.
		/// @param tol Thershold for singular values. Singular values smaller
		///            then threshold are considered to be zero.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr = 1, in case of no failure *ierr = 0.
		/// @return A pair of pseudo-inverse matrix and boolean. If boolean is true,
		///         then the matrix was inverted successfully.
		Matrix<Var, pool_array_t<Var> > PseudoInvert(INMOST_DATA_REAL_TYPE tol = 0, int * ierr = NULL) const;
		/// Solves the system of equations of the form A*X=B, with A and B matrices.
		/// Uses Moore-Penrose pseudo-inverse of the matrix A and calculates X = A^+*B.
		/// @param B Matrix at the right hand side.
		/// @param tol Thershold for singular values. Singular values smaller
		///            then threshold are considered to be zero.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr = 1, in case of no failure *ierr = 0.
		/// @return A pair of the solution matrix X and boolean. If boolean is true,
		///         then the matrix was inverted successfully.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
		PseudoSolve(const AbstractMatrix<typeB> & B, INMOST_DATA_REAL_TYPE tol = 0, int * ierr = NULL) const;
		/// Extract submatrix of a matrix.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param ibeg Starting row in the original matrix.
		/// @param iend Last row (excluded) in the original matrix.
		/// @param jbeg Starting column in the original matrix.
		/// @param jend Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		Matrix<Var, pool_array_t<Var> > ExtractSubMatrix(enumerator ibeg, enumerator iend, enumerator jbeg, enumerator jend) const;
		/// Change representation of the matrix into matrix of another size.
		/// Useful to change representation from matrix into vector and back.
		/// Replaces original number of columns and rows with a new one.
		/// @return Matrix with same entries and provided number of rows and columns.
		Matrix<Var, pool_array_t<Var> > Repack(enumerator rows, enumerator cols) const;
	};
	
	
	template<typename Var, typename storage_type>
	class SymmetricMatrix : public AbstractMatrix<Var>
	{
	public:
		typedef unsigned enumerator; //< Integer type for indexes.
	protected:
		storage_type space; //< Array of row-wise stored elements.
		enumerator n; //< Number of rows.
	public:
		///Exchange contents of two matrices.
		void Swap(AbstractMatrix<Var> & b)
		{
#if defined(_CPPRTTI) || defined(__GXX_RTTI)
			SymmetricMatrix<Var,storage_type> * bb = dynamic_cast<SymmetricMatrix<Var,storage_type> *>(&b);
			if( bb != NULL )
			{
				space.swap((*bb).space);
				std::swap(n,(*bb).n);
			}
			else AbstractMatrix<Var>::Swap(b);
#else //_CPPRTTI
			AbstractMatrix<Var>::Swap(b);
#endif //_CPPRTTI
		}
		/// Construct empty matrix.
		SymmetricMatrix() : space(), n(0) {}
		/// Construct the matrix from provided array and sizes.
		/// @param pspace Array of elements of the matrix, stored in row-wise format.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		SymmetricMatrix(const Var * pspace, enumerator pn) : space(pspace,pspace+pn*(pn+1)/2), n(pn) {}
		/// Construct the matrix with the provided storage with known size.
		/// Could be used to wrap existing array.
		/// \warning The size of the provided container is assumed to be pn*pm.
		/// @param pspace Storage of elements of the matrix, stored in row-wise format.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// \todo Do we need reference for pspace or just pspace?
		SymmetricMatrix(const storage_type & pspace, enumerator pn) : space(pspace), n(pn) {}
		/// Construct the matrix with the provided storage and unknown size.
		/// Could be used to wrap existing array.
		/// \warning Have to call Resize afterwards.
		/// @param pspace Storage of elements of the matrix, stored in row-wise format.
		/// \todo Do we need reference for pspace or just pspace?
		SymmetricMatrix(const storage_type & pspace) : space(pspace), n(0) {}
		/// Construct a matrix with provided sizes.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// \warning The matrix does not necessery have zero entries.
		SymmetricMatrix(enumerator pn) : space(pn*(pn+1)/2), n(pn) {}
		/// Construct a matrix with provided sizes and fills with value.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// @param c Value to fill the matrix.
		SymmetricMatrix(enumerator pn, const Var & c) : space(pn*(pn+1)/2,c), n(pn) {}
		/// Copy matrix.
		/// @param other Another matrix of the same type.
		SymmetricMatrix(const SymmetricMatrix & other) : space(other.space), n(other.n)
		{
			//for(enumerator i = 0; i < n*m; ++i)
			//	space[i] = other.space[i];
		}
		/// Construct matrix from matrix of different type.
		/// Function assumes that the other matrix is square and symmetric.
		/// Copies only top-right triangular part.
		/// Uses assign function declared in inmost_expression.h.
		/// Copies derivative information if possible.
		/// @param other Another matrix of different type.
		template<typename typeB>
		SymmetricMatrix(const AbstractMatrix<typeB> & other) : space(other.Rows()*(other.Rows()+1)/2), n(other.Rows())
		{
			assert(other.Rows() == other.Cols());
			for(enumerator i = 0; i < n; ++i)
				for(enumerator j = i; j < n; ++j)
					assign((*this)(i,j),other(i,j));
		}
		/// Delete matrix.
		~SymmetricMatrix() {}
		/// Resize the matrix into different size.
		/// Number of rows must match number of columns for symmetric matrix.
		/// @param nrows New number of rows.
		/// @param ncols New number of columns.
		void Resize(enumerator nrows, enumerator ncols)
		{
			assert(nrows == ncols);
			if( space.size() != (nrows+1)*nrows/2 )
				space.resize((nrows+1)*nrows/2);
			n = nrows;
		}
		/// Assign matrix of the same type.
		/// @param other Another matrix of the same type.
		/// @return Reference to matrix.
		SymmetricMatrix & operator =(SymmetricMatrix const & other)
		{
			if( this != &other )
			{
				if( n != other.n )
					space.resize(other.n*(other.n+1)/2);
				for(enumerator i = 0; i < other.n*(other.n+1)/2; ++i)
					space[i] = other.space[i];
				n = other.n;
			}
			return *this;
		}
		/// Assign matrix of another type.
		/// Function assumes that the other matrix is square and symmetric.
		/// Copies only top-right triangular part.
		/// @param other Another matrix of different type.
		/// @return Reference to matrix.
		template<typename typeB>
		SymmetricMatrix & operator =(AbstractMatrix<typeB> const & other)
		{
			assert(other.Rows() == other.Cols());
			if( Rows() != other.Rows() )
				space.resize((other.Rows()+1)*other.Rows()+1);
			for(enumerator i = 0; i < other.Rows(); ++i)
				for(enumerator j = i; j < other.Cols(); ++j)
					assign((*this)(i,j),other(i,j));
			n = other.Rows();
			return *this;
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		Var & operator()(enumerator i, enumerator j)
		{
			assert(i >= 0 && i < n);
			assert(j >= 0 && j < n);
			if( i > j ) std::swap(i,j);
			return space[j+n*i-i*(i+1)/2];
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		const Var & operator()(enumerator i, enumerator j) const
		{
			assert(i >= 0 && i < n);
			assert(j >= 0 && j < n);
			if( i > j ) std::swap(i,j);
			return space[j+n*i-i*(i+1)/2];
		}
		
		/// Return raw pointer to matrix data, stored in row-wise format.
		/// @return Pointer to data.
		Var * data() {return space.data();}
		/// Return raw pointer to matrix data without right of change,
		/// stored in row-wise format.
		/// @return Pointer to constant data.
		const Var * data() const {return space.data();}
		/// Obtain number of rows.
		/// @return Number of rows.
		enumerator Rows() const {return n;}
		/// Obtain number of columns.
		/// @return Number of columns.
		enumerator Cols() const {return n;}
		/// Obtain number of rows.
		/// @return Reference to number of rows.
		enumerator & Rows() {return n;}
		/// Obtain number of rows.
		/// @return Reference to number of columns.
		enumerator & Cols() {return n;}
		/// Convert values in array into square matrix.
		/// Supports the following representation, depending on the size
		/// of input array and size of side of final tensors' matrix:
		///
		/// representation | (array size, tensor size)
		///
		/// scalar         | (1,1), (1,2), (1,3), (1,6)
		///
		/// diagonal       | (2,2), (3,3), (6,6)
		///
		/// symmetric      | (3,2), (6,3), (21,6)
		///
		/// For symmetric matrix elements in array are enumerated row by
		/// row starting from diagonal.
		/// @param K Array of elements to be converted into tensor.
		/// @param size Size of the input array.
		/// @param matsize Size of the final tensor.
		/// @return Matrix of the tensor of size matsize by matsize.
		static SymmetricMatrix<Var, pool_array_t<Var> > FromTensor(const Var * K, enumerator size, enumerator matsize = 3)
		{
			SymmetricMatrix<Var, pool_array_t<Var> > Kc(matsize,matsize);
			if( matsize == 1 )
			{
				assert(size == 1);
				Kc(0,0) = K[0];
			}
			if( matsize == 2 )
			{
				assert(size == 1 || size == 2 || size == 3 || size == 4);
				switch(size)
				{
					case 1: //scalar
						Kc(0,0) = Kc(1,1) = K[0];
						break;
					case 2: //diagonal
						Kc(0,0) = K[0]; // KXX
						Kc(1,1) = K[1]; // KYY
						break;
					case 3: //symmetric
						Kc(0,0) = K[0]; // KXX
						Kc(0,1) = K[1]; //KXY
						Kc(1,1) = K[2]; //KYY
						break;
				}
			}
			else if( matsize == 3 )
			{
				assert(size == 1 || size == 3 || size == 6 || size == 9);
				switch(size)
				{
					case 1: //scalar permeability tensor
						Kc(0,0) = Kc(1,1) = Kc(2,2) = K[0];
						break;
					case 3: //diagonal permeability tensor
						Kc(0,0) = K[0]; //KXX
						Kc(1,1) = K[1]; //KYY
						Kc(2,2) = K[2]; //KZZ
						break;
					case 6: //symmetric permeability tensor
						Kc(0,0) = K[0]; //KXX
						Kc(0,1) = K[1]; //KXY
						Kc(0,2) = K[2]; //KXZ
						Kc(1,1) = K[3]; //KYY
						Kc(1,2) = K[4]; //KYZ
						Kc(2,2) = K[5]; //KZZ
						break;
				}
			}
			else if( matsize == 6 )
			{
				assert(size == 1 || size == 6 || size == 21 || size == 36);
				switch(size)
				{
					case 1: //scalar elasticity tensor
						Kc(0,0) = Kc(1,1) = Kc(2,2) = Kc(3,3) = Kc(4,4) = Kc(5,5) = K[0];
						break;
					case 6: //diagonal elasticity tensor
						Kc(0,0) = K[0]; //KXX
						Kc(1,1) = K[1]; //KYY
						Kc(2,2) = K[2]; //KZZ
						break;
					case 21: //symmetric elasticity tensor (note - diagonal first, then off-diagonal rows)
					{
						Kc(0,0) = K[0]; //c11
						Kc(0,1) = K[1]; //c12
						Kc(0,2) = K[2]; //c13
						Kc(0,3) = K[3]; //c14
						Kc(0,4) = K[4]; //c15
						Kc(0,5) = K[5]; //c16
						Kc(1,1) = K[6]; //c22
						Kc(1,2) = K[7]; //c23
						Kc(1,3) = K[8]; //c24
						Kc(1,4) = K[9]; //c25
						Kc(1,5) = K[10]; //c26
						Kc(2,2) = K[11]; //c33
						Kc(2,3) = K[12]; //c34
						Kc(2,4) = K[13]; //c35
						Kc(2,5) = K[14]; //c36
						Kc(3,3) = K[15]; //c44
						Kc(3,4) = K[16]; //c45
						Kc(3,5) = K[17]; //c46
						Kc(4,4) = K[18]; //c55
						Kc(4,5) = K[19]; //c56
						Kc(5,5) = K[20]; //c66
						break;
					}
				}
			}
			return Kc;
		}
		/// Create diagonal matrix from array
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by array, other elements are zero.
		static SymmetricMatrix<Var, pool_array_t<Var> > FromDiagonal(const Var * r, enumerator size)
		{
			SymmetricMatrix<Var, pool_array_t<Var> > ret(size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
			return ret;
		}
		/// Create diagonal matrix from array of values that have to be inversed.
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by inverse of array elements.
		static SymmetricMatrix<Var, pool_array_t<Var> > FromDiagonalInverse(const Var * r, enumerator size)
		{
			SymmetricMatrix<Var, pool_array_t<Var> > ret(size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = 1.0/r[k];
			return ret;
		}
		/// Unit matrix. Creates a square matrix of size pn by pn
		/// and fills the diagonal with c.
		/// @param pn Number of rows and columns in the matrix.
		/// @param c Value to put onto diagonal.
		/// @return Returns a unit matrix.
		static SymmetricMatrix<Var, pool_array_t<Var> > Unit(enumerator pn, const Var & c = 1.0)
		{
			SymmetricMatrix<Var, pool_array_t<Var> > ret(pn,0.0);
			for(enumerator i = 0; i < pn; ++i) ret(i,i) = c;
			return ret;
		}
		
		/// Extract submatrix of a matrix for in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param first_row Starting row in the original matrix.
		/// @param last_row Last row (excluded) in the original matrix.
		/// @param first_col Starting column in the original matrix.
		/// @param last_col Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		//::INMOST::SubMatrix<Var,storage_type> SubMatrix(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);
		
		/// Extract submatrix of a matrix for in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param first_row Starting row in the original matrix.
		/// @param last_row Last row (excluded) in the original matrix.
		/// @param first_col Starting column in the original matrix.
		/// @param last_col Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		::INMOST::SubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);
		
		::INMOST::ConstSubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const;
	};
	
	/// Class for linear algebra operations on dense matrices.
	/// Matrix with n rows and m columns.
	///
	///   __m__
	///  |     |
	/// n|     |
	///  |_____|
	///
	/// \todo:
	/// 1. expression templates for operations
	///    (???) how to for multiplication?
	///    efficient multiplication would require all the
	///    matrix elements to be precomputed.
	///    consider number 5 instead.
	/// 2. (ok) template matrix type for AD variables
	/// 3. (ok,test) template container type for data storage.
	/// 4. (ok,test) option for wrapper container around provided data storage.
	///    (to perform matrix operations with existing data)
	/// 5. consider multi-threaded stack to get space for
	///    matrices for local operations and returns.
	/// 6. class SubMatrix for fortran-like access to matrix.
	/// 7. Uniform implementation of algorithms for Matrix and Submatrix.
	///    to achieve: make abdstract class with abstract element access
	///    operator, make matrix and submatrix ancestors of that class
	template<typename Var, typename storage_type>
	class Matrix : public AbstractMatrix<Var>
	{
	public:
		typedef unsigned enumerator; //< Integer type for indexes.
	protected:
		//array<Var> space; //< Array of row-wise stored elements.
		storage_type space; //< Array of row-wise stored elements.
		enumerator n; //< Number of rows.
		enumerator m; //< Number of columns.
	public:
		/// Erase single row.
		/// @param row Position of row, should be from 0 to Matrix::Rows()-1.
		void RemoveRow(enumerator row)
		{
			assert(row < n);
			for(enumerator k = row+1; k < n; ++k)
			{
				for(enumerator l = 0; l < m; ++l)
					(*this)(k-1,l) = (*this)(k,l);
			}
			space.resize((n-1)*m);
			--n;
		}
		/// Erase multiple rows. Assumes first <= last.
		/// If first == last, then no row is deleted.
		/// @param first Position of the first row, should be from 0 to Matrix::Rows()-1.
		/// @param last Position behind the last row, should be from 0 to Matrix::Rows().
		/// \todo check
		void RemoveRows(enumerator first, enumerator last)
		{
			assert(first < n);
			assert(last <= n);
			assert(first <= last);
			enumerator shift = last-first;
			for(enumerator k = last; k < n; ++k)
			{
				for(enumerator l = 0; l < m; ++l)
					(*this)(k-shift,l) = (*this)(k,l);
			}
			space.resize((n-shift)*m);
			n-=shift;
		}
		/// Erase single column.
		/// @param row Position of column, should be from 0 to Matrix::Cols()-1.
		void RemoveColumn(enumerator col)
		{
			assert(col < m);
			Matrix<Var> tmp(n,m-1);
			for(enumerator k = 0; k < n; ++k)
			{
				for(enumerator l = 0; l < col; ++l)
					tmp(k,l) = (*this)(k,l);
				for(enumerator l = col+1; l < m; ++l)
					tmp(k,l-1) = (*this)(k,l);
			}
			this->Swap(tmp);
		}
		/// Erase multiple columns. Assumes first <= last.
		/// If first == last, then no column is deleted.
		/// @param first Position of the first column, should be from 0 to Matrix::Cols()-1.
		/// @param last Position behind the last column, should be from 0 to Matrix::Cols().
		/// \todo check
		void RemoveColumns(enumerator first, enumerator last)
		{
			assert(first < m);
			assert(last <= m);
			assert(first <= last);
			enumerator shift = last-first;
			Matrix<Var> tmp(n,m-shift);
			for(enumerator k = 0; k < n; ++k)
			{
				for(enumerator l = 0; l < first; ++l)
					tmp(k,l) = (*this)(k,l);
				for(enumerator l = last; l < m; ++l)
					tmp(k,l-shift) = (*this)(k,l);
			}
			this->Swap(tmp);
		}
		/// Erase part of the matrix.
		/// @param firstrow Position of the first row, should be from 0 to Matrix::Rows()-1.
		/// @param lastrow Position behind the last row, should be from 0 to Matrix::Rows().
		/// @param firstcol Position of the first column, should be from 0 to Matrix::Cols()-1.
		/// @param lastcol Position behind the last column, should be from 0 to Matrix::Cols().
		void RemoveSubset(enumerator firstrow, enumerator lastrow, enumerator firstcol, enumerator lastcol)
		{
			enumerator shiftrow = lastrow-firstrow;
			enumerator shiftcol = lastcol-firstcol;
			Matrix<Var> tmp(n-shiftrow, m-shiftcol);
			for(enumerator k = 0; k < firstrow; ++k)
			{
				for(enumerator l = 0; l < firstcol; ++l)
					tmp(k,l) = (*this)(k,l);
				for(enumerator l = lastcol; l < m; ++l)
					tmp(k,l-shiftcol) = (*this)(k,l);
			}
			for(enumerator k = lastrow; k < n; ++k)
			{
				for(enumerator l = 0; l < firstcol; ++l)
					tmp(k-shiftrow,l) = (*this)(k,l);
				for(enumerator l = lastcol; l < m; ++l)
					tmp(k-shiftrow,l-shiftcol) = (*this)(k,l);
			}
			this->Swap(tmp);
		}
		///Exchange contents of two matrices.
		void Swap(AbstractMatrix<Var> & b)
		{
#if defined(_CPPRTTI) || defined(__GXX_RTTI)
			Matrix<Var,storage_type> * bb = dynamic_cast<Matrix<Var,storage_type> *>(&b);
			if( bb != NULL )
			{
				space.swap((*bb).space);
				std::swap(n,(*bb).n);
				std::swap(m,(*bb).m);
			}
			else AbstractMatrix<Var>::Swap(b);
#else //_CPPRTTI
			AbstractMatrix<Var>::Swap(b);
#endif //_CPPRTTI	
		}
		/// Construct empty matrix.
		Matrix() : space(), n(0), m(0) {}
		/// Construct the matrix from provided array and sizes.
		/// @param pspace Array of elements of the matrix, stored in row-wise format.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		Matrix(const Var * pspace, enumerator pn, enumerator pm) : space(pspace,pspace+pn*pm), n(pn), m(pm) {}
		/// Construct the matrix with the provided storage with known size.
		/// Could be used to wrap existing array.
		/// \warning The size of the provided container is assumed to be pn*pm.
		/// @param pspace Storage of elements of the matrix, stored in row-wise format.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// \todo Do we need reference for pspace or just pspace?
		Matrix(const storage_type & pspace, enumerator pn, enumerator pm) : space(pspace), n(pn), m(pm) {}
		/// Construct the matrix with the provided storage and unknown size.
		/// Could be used to wrap existing array.
		/// \warning Have to call Resize afterwards.
		/// @param pspace Storage of elements of the matrix, stored in row-wise format.
		/// \todo Do we need reference for pspace or just pspace?
		Matrix(const storage_type & pspace) : space(pspace), n(0), m(0) {}
		/// Construct a matrix with provided sizes.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// \warning The matrix does not necessery have zero entries.
		Matrix(enumerator pn, enumerator pm) : space(pn*pm), n(pn), m(pm) {}
		/// Construct a matrix with provided sizes and fills with value.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// @param c Value to fill the matrix.
		Matrix(enumerator pn, enumerator pm, const Var & c) : space(pn*pm,c), n(pn), m(pm) {}
		/// Copy matrix.
		/// @param other Another matrix of the same type.
		Matrix(const Matrix & other) : space(other.space), n(other.n), m(other.m)
		{
			//for(enumerator i = 0; i < n*m; ++i)
			//	space[i] = other.space[i];
		}
		/// Construct matrix from matrix of different type.
		/// Uses assign function declared in inmost_expression.h.
		/// Copies derivative information if possible.
		/// @param other Another matrix of different type.
		template<typename typeB>
		Matrix(const AbstractMatrix<typeB> & other) : space(other.Cols()*other.Rows()), n(other.Rows()), m(other.Cols())
		{
			for(enumerator i = 0; i < n; ++i)
				for(enumerator j = 0; j < m; ++j)
					assign((*this)(i,j),other(i,j));
		}
		/// Delete matrix.
		~Matrix() {}
		/// Resize the matrix into different size.
		/// @param nrows New number of rows.
		/// @param ncols New number of columns.
		void Resize(enumerator nrows, enumerator mcols)
		{
			if( space.size() != mcols*nrows )
				space.resize(mcols*nrows);
			n = nrows;
			m = mcols;
		}
		/// Assign matrix of the same type.
		/// @param other Another matrix of the same type.
		/// @return Reference to matrix.
		Matrix & operator =(Matrix const & other)
		{
			if( this != &other )
			{
				if( n*m != other.n*other.m )
					space.resize(other.n*other.m);
				for(enumerator i = 0; i < other.n*other.m; ++i)
					space[i] = other.space[i];
				n = other.n;
				m = other.m;
			}
			return *this;
		}
		/// Assign matrix of another type.
		/// @param other Another matrix of different type.
		/// @return Reference to matrix.
		template<typename typeB>
		Matrix & operator =(AbstractMatrix<typeB> const & other)
		{
			if( Cols()*Rows() != other.Cols()*other.Rows() )
				space.resize(other.Cols()*other.Rows());
			for(enumerator i = 0; i < other.Rows(); ++i)
				for(enumerator j = 0; j < other.Cols(); ++j)
					assign((*this)(i,j),other(i,j));
			n = other.Rows();
			m = other.Cols();
			return *this;
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		Var & operator()(enumerator i, enumerator j)
		{
			assert(i >= 0 && i < n);
			assert(j >= 0 && j < m);
			assert(i*m+j < n*m); //overflow check?
			return space[i*m+j];
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		const Var & operator()(enumerator i, enumerator j) const
		{
			assert(i >= 0 && i < n);
			assert(j >= 0 && j < m);
			assert(i*m+j < n*m); //overflow check?
			return space[i*m+j];
		}
		
		/// Return raw pointer to matrix data, stored in row-wise format.
		/// @return Pointer to data.
		Var * data() {return space.data();}
		/// Return raw pointer to matrix data without right of change,
		/// stored in row-wise format.
		/// @return Pointer to constant data.
		const Var * data() const {return space.data();}
		/// Obtain number of rows.
		/// @return Number of rows.
		enumerator Rows() const {return n;}
		/// Obtain number of columns.
		/// @return Number of columns.
		enumerator Cols() const {return m;}
		/// Obtain number of rows.
		/// @return Reference to number of rows.
		enumerator & Rows() {return n;}
		/// Obtain number of rows.
		/// @return Reference to number of columns.
		enumerator & Cols() {return m;}
		/// Construct row permutation matrix from array of new positions for rows.
		/// Row permutation matrix multiplies matrix from left.
		/// Column permutation matrix is obtained by transposition and is multiplied
		/// from the right.
		/// Argument Perm is filled with the new position for rows, i.e.
		/// i-th row takes new position Perm[i]
		/// @param Perm Array with new positions for rows.
		/// @param size Size of the array and the resulting matrix.
		/// @return Permutation matrix.
		static Matrix<Var, pool_array_t<Var> > Permutation(const INMOST_DATA_ENUM_TYPE * Perm, enumerator size)
		{
			Matrix<Var, pool_array_t<Var> > Ret(size,size);
			Ret.Zero();
			for(enumerator k = 0; k < size; ++k)
				Ret(k,Perm[k]) = 1;
			return Ret;
		}
		/// Convert values in array into square matrix.
		/// Supports the following representation, depending on the size
		/// of input array and size of side of final tensors' matrix:
		///
		/// representation | (array size, tensor size)
		///
		/// scalar         | (1,1), (1,2), (1,3), (1,6)
		///
		/// diagonal       | (2,2), (3,3), (6,6)
		///
		/// symmetric      | (3,2), (6,3), (21,6)
		///
		/// full           | (4,2), (9,3), (36,6)
		///
		/// For full matrix elements in array are enumerated row by row.
		/// For symmetric matrix elements in array are enumerated row by
		/// row starting from diagonal.
		/// @param K Array of elements to be converted into tensor.
		/// @param size Size of the input array.
		/// @param matsize Size of the final tensor.
		/// @return Matrix of the tensor of size matsize by matsize.
		static Matrix<Var, pool_array_t<Var> > FromTensor(const Var * K, enumerator size, enumerator matsize = 3)
		{
			Matrix<Var, pool_array_t<Var> > Kc(matsize,matsize);
			if( matsize == 1 )
			{
				assert(size == 1);
				Kc(0,0) = K[0];
			}
			if( matsize == 2 )
			{
				assert(size == 1 || size == 2 || size == 3 || size == 4);
				switch(size)
				{
					case 1: //scalar
						Kc(0,0) = Kc(1,1) = K[0];
						break;
					case 2: //diagonal
						Kc(0,0) = K[0]; // KXX
						Kc(1,1) = K[1]; // KYY
						break;
					case 3: //symmetric
						Kc(0,0) = K[0]; // KXX
						Kc(0,1) = Kc(1,0) = K[1]; //KXY
						Kc(1,1) = K[2]; //KYY
						break;
					case 4: //full
						Kc(0,0) = K[0]; //KXX
						Kc(0,1) = K[1]; //KXY
						Kc(1,0) = K[2]; //KYX
						Kc(1,1) = K[3]; //KYY
						break;
				}
			}
			else if( matsize == 3 )
			{
				assert(size == 1 || size == 3 || size == 6 || size == 9);
				switch(size)
				{
					case 1: //scalar permeability tensor
						Kc(0,0) = Kc(1,1) = Kc(2,2) = K[0];
						break;
					case 3: //diagonal permeability tensor
						Kc(0,0) = K[0]; //KXX
						Kc(1,1) = K[1]; //KYY
						Kc(2,2) = K[2]; //KZZ
						break;
					case 6: //symmetric permeability tensor
						Kc(0,0) = K[0]; //KXX
						Kc(0,1) = Kc(1,0) = K[1]; //KXY
						Kc(0,2) = Kc(2,0) = K[2]; //KXZ
						Kc(1,1) = K[3]; //KYY
						Kc(1,2) = Kc(2,1) = K[4]; //KYZ
						Kc(2,2) = K[5]; //KZZ
						break;
					case 9: //full permeability tensor
						Kc(0,0) = K[0]; //KXX
						Kc(0,1) = K[1]; //KXY
						Kc(0,2) = K[2]; //KXZ
						Kc(1,0) = K[3]; //KYX
						Kc(1,1) = K[4]; //KYY
						Kc(1,2) = K[5]; //KYZ
						Kc(2,0) = K[6]; //KZX
						Kc(2,1) = K[7]; //KZY
						Kc(2,2) = K[8]; //KZZ
						break;
				}
			}
			else if( matsize == 6 )
			{
				assert(size == 1 || size == 6 || size == 21 || size == 36);
				switch(size)
				{
					case 1: //scalar elasticity tensor
						Kc(0,0) = Kc(1,1) = Kc(2,2) = Kc(3,3) = Kc(4,4) = Kc(5,5) = K[0];
						break;
					case 6: //diagonal elasticity tensor
						Kc(0,0) = K[0]; //KXX
						Kc(1,1) = K[1]; //KYY
						Kc(2,2) = K[2]; //KZZ
						break;
					case 21: //symmetric elasticity tensor (note - diagonal first, then off-diagonal rows)
					{
						Kc(0,0) = K[0]; //c11
						Kc(0,1) = Kc(1,0) = K[1]; //c12
						Kc(0,2) = Kc(2,0) = K[2]; //c13
						Kc(0,3) = Kc(3,0) = K[3]; //c14
						Kc(0,4) = Kc(4,0) = K[4]; //c15
						Kc(0,5) = Kc(5,0) = K[5]; //c16
						Kc(1,1) = K[6]; //c22
						Kc(1,2) = Kc(2,1) = K[7]; //c23
						Kc(1,3) = Kc(3,1) = K[8]; //c24
						Kc(1,4) = Kc(4,1) = K[9]; //c25
						Kc(1,5) = Kc(5,1) = K[10]; //c26
						Kc(2,2) = K[11]; //c33
						Kc(2,3) = Kc(3,2) = K[12]; //c34
						Kc(2,4) = Kc(4,2) = K[13]; //c35
						Kc(2,5) = Kc(5,2) = K[14]; //c36
						Kc(3,3) = K[15]; //c44
						Kc(3,4) = Kc(4,3) = K[16]; //c45
						Kc(3,5) = Kc(5,3) = K[17]; //c46
						Kc(4,4) = K[18]; //c55
						Kc(4,5) = Kc(5,4) = K[19]; //c56
						Kc(5,5) = K[20]; //c66
						break;
					}
					case 36: //full elasticity tensor
						for(int i = 0; i < 6; ++i)
							for(int j = 0; j < 6; ++j)
								Kc(i,j) = K[6*i+j];
						break;
				}
			}
			return Kc;
		}
		/// Create column-vector in matrix form from array.
		/// @param r Array of elements of the vector.
		/// @param size Size of the vector.
		/// @return Vector with contents of the array.
		static Matrix<Var, pool_array_t<Var> > FromVector(const Var * r, enumerator size)
		{
			return Matrix<Var, pool_array_t<Var> >(r,size,1);
		}
		/// Create diagonal matrix from array
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by array, other elements are zero.
		static Matrix<Var, pool_array_t<Var> > FromDiagonal(const Var * r, enumerator size)
		{
			Matrix<Var, pool_array_t<Var> > ret(size,size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
			return ret;
		}
		/// Create diagonal matrix from array of values that have to be inversed.
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by inverse of array elements.
		static Matrix<Var, pool_array_t<Var> > FromDiagonalInverse(const Var * r, enumerator size)
		{
			Matrix<Var, pool_array_t<Var> > ret(size,size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = 1.0/r[k];
			return ret;
		}
		/// Cross-product matrix. Converts an array of 3 elements
		/// representing a vector into matrix, helps replace
		/// a cross product of two vectors by multiplication of matrix
		/// and vector. For a x b equivalent is CrossProduct(a)*b.
		/// @param vec Array of elements representing a vector.
		/// @return A matrix representing cross product.
		static Matrix<Var, pool_array_t<Var> > CrossProduct(const Var vec[3])
		{
			// |  0  -z   y |
			// |  z   0  -x |
			// | -y   x   0 |
			Matrix<Var, pool_array_t<Var> > ret(3,3);
			ret(0,0) = 0.0;
			ret(0,1) = -vec[2]; //-z
			ret(0,2) = vec[1]; //y
			ret(1,0) = vec[2]; //z
			ret(1,1) = 0;
			ret(1,2) = -vec[0]; //-x
			ret(2,0) = -vec[1]; //-y
			ret(2,1) = vec[0]; //x
			ret(2,2) = 0;
			return ret;
		}
		/// Unit matrix. Creates a square matrix of size pn by pn
		/// and fills the diagonal with c.
		/// @param pn Number of rows and columns in the matrix.
		/// @param c Value to put onto diagonal.
		/// @return Returns a unit matrix.
		static Matrix<Var, pool_array_t<Var> > Unit(enumerator pn, const Var & c = 1.0)
		{
			Matrix<Var, pool_array_t<Var> > ret(pn,pn);
			ret.Zero();
			for(enumerator i = 0; i < pn; ++i) ret(i,i) = c;
			return ret;
		}
		/// Matix with 1 row, Create a matrix of size 1 by pn and
		/// fills it with c.
		/// @param pn Number of columns.
		/// @param c Value to fill the matrix.
		/// @return Returns a matrix with 1 row.
		static Matrix<Var, pool_array_t<Var> > Row(enumerator pn, const Var & c = 1.0)
		{
			return Matrix<Var, pool_array_t<Var> >(1,pn,c);
		}
		/// Matix with 1 column, Create a matrix of size pn by 1 and
		/// fills it with c.
		/// @param pn Number of rows.
		/// @param c Value to fill the matrix.
		/// @return Returns a matrix with 1 column.
		static Matrix<Var, pool_array_t<Var> > Col(enumerator pn, const Var & c = 1.0)
		{
			return Matrix<Var, pool_array_t<Var> >(pn, 1, c);
		}
		/// Concatenate B matrix as columns of current matrix.
		/// Assumes that number of rows of current matrix is
		/// equal to number of rows of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatRows
		Matrix<Var, pool_array_t<Var> > ConcatCols(const Matrix & B) const
		{
			assert(Rows() == B.Rows());
			Matrix<Var, pool_array_t<Var> > ret(Rows(),Cols()+B.Cols());
			const Matrix & A = *this;
			for(enumerator i = 0; i < Rows(); ++i)
			{
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i,j) = A(i,j);
				for(enumerator j = 0; j < B.Cols(); ++j)
					ret(i,j+Cols()) = B(i,j);
			}
			return ret;
		}
		/// Concatenate B matrix as rows of current matrix.
		/// Assumes that number of colums of current matrix is
		/// equal to number of columns of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatCols
		Matrix<Var, pool_array_t<Var> > ConcatRows(const Matrix & B) const
		{
			assert(Cols() == B.Cols());
			Matrix<Var, pool_array_t<Var> > ret(Rows()+B.Rows(),Cols());
			const Matrix & A = *this;
			for(enumerator i = 0; i < Rows(); ++i)
			{
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i,j) = A(i,j);
			}
			for(enumerator i = 0; i < B.Rows(); ++i)
			{
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i+Rows(),j) = B(i,j);
			}
			return ret;
		}
		
		/// Joint diagonalization algorithm by Cardoso.
		/// Source http://perso.telecom-paristech.fr/~cardoso/Algo/Joint_Diag/joint_diag_r.m
		/// Current matrix should have size n by n*m
		/// And represent concatination of m n by n matrices.
		/// Current matrix is replaced by diagonalized matrices.
		/// For correct result requires that input matrices are
		/// exectly diagonalizable, otherwise the result may be approximate.
		/// @param threshold Optional small number.
		/// @return A unitary n by n matrix V used to diagonalize array of
		/// initial matrices. Current matrix is replaced by concatination of
		/// V^T*A_i*V, a collection of diagonalized matrices.
		Matrix<Var, pool_array_t<Var> > JointDiagonalization(INMOST_DATA_REAL_TYPE threshold = 1.0e-7)
		{
			enumerator N = Rows();
			enumerator M = Cols() / Rows();
			Matrix<Var, pool_array_t<Var> > V(Matrix<Var, pool_array_t<Var> >::Unit(m));
			Matrix<Var, pool_array_t<Var> > R(2,M);
			Matrix<Var, pool_array_t<Var> > G(2,2);
			Matrix & A = *this;
			Var ton, toff, theta, c, s, Ap, Aq, Vp, Vq;
			bool repeat;
			do
			{
				repeat = false;
				for(enumerator p = 0; p < N-1; ++p)
				{
					for(enumerator q = p+1; q < N; ++q)
					{
						for(enumerator k = 0; k < M; ++k)
						{
							R(0,k) = A(p,p + k*N) - A(q,q + k*N);
							R(1,k) = A(p,q + k*N) + A(q,p + k*N);
						}
						G = R*R.Transpose();
						Var ton  = G(0,0) - G(1,1);
						Var toff = G(0,1) + G(1,0);
						Var theta = 0.5 * atan2( toff, ton + sqrt(ton*ton + toff*toff) );
						Var c = cos(theta);
						Var s = sin(theta);
						if( fabs(s) > threshold )
						{
							//std::cout << "p,q: " << p << "," << q << " c,s: " << c << "," << s << std::endl;
							repeat = true;
							for(enumerator k = 0; k < M; ++k)
							{
								for(enumerator i = 0; i < N; ++i)
								{
									Ap = A(i,p + k*N);
									Aq = A(i,q + k*N);
									A(i,p + k*N) = Ap*c + Aq*s;
									A(i,q + k*N) = Aq*c - Ap*s;
								}
							}
							for(enumerator k = 0; k < M; ++k)
							{
								for(enumerator j = 0; j < N; ++j)
								{
									Ap = A(p,j + k*N);
									Aq = A(q,j + k*N);
									A(p,j + k*N) = Ap*c + Aq*s;
									A(q,j + k*N) = Aq*c - Ap*s;
								}
							}
							for(enumerator i = 0; i < N; ++i)
							{
								Vp = V(i,p);
								Vq = V(i,q);
								V(i,p) = Vp*c + Vq*s;
								V(i,q) = Vq*c - Vp*s;
							}
						}
					}
				}
				//Print();
			} while( repeat );
			return V;
		}
		/// Extract submatrix of a matrix for in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param first_row Starting row in the original matrix.
		/// @param last_row Last row (excluded) in the original matrix.
		/// @param first_col Starting column in the original matrix.
		/// @param last_col Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		//::INMOST::SubMatrix<Var,storage_type> SubMatrix(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);
		
		/// Extract submatrix of a matrix for in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param first_row Starting row in the original matrix.
		/// @param last_row Last row (excluded) in the original matrix.
		/// @param first_col Starting column in the original matrix.
		/// @param last_col Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		::INMOST::SubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);

		::INMOST::ConstSubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const;
	};
	/// This class allows for in-place operations on submatrix of the matrix elements.
	template<typename Var>
	class SubMatrix : public AbstractMatrix<Var>
	{
	public:
		typedef unsigned enumerator; //< Integer type for indexes.
	private:
		AbstractMatrix<Var> * M;
		enumerator brow; //< First row in matrix M.
		enumerator erow; //< Last row in matrix M.
		enumerator bcol; //< First column in matrix M.
		enumerator ecol; //< Last column in matrix M.
	public:
		/// Number of rows in submatrix.
		/// @return Number of rows.
		enumerator Rows() const {return erow-brow;}
		/// Number of columns in submatrix.
		/// @return Number of columns.
		enumerator Cols() const {return ecol-bcol;}
		/// Create submatrix for a matrix.
		/// @param rM Reference to the matrix that stores elements.
		/// @param first_row First row in the matrix.
		/// @param last_row Last row in the matrix.
		/// @param first_column First column in the matrix.
		/// @param last_column Last column in the matrix.
		SubMatrix(AbstractMatrix<Var> & rM, enumerator first_row, enumerator last_row, enumerator first_column, enumerator last_column) : M(&rM), brow(first_row), erow(last_row), bcol(first_column), ecol(last_column)
		{}
		SubMatrix(const SubMatrix & b) : M(b.M), brow(b.brow), erow(b.erow), bcol(b.bcol), ecol(b.ecol) {}
		/// Assign matrix of another type to submatrix.
		/// @param other Another matrix of different type.
		/// @return Reference to current submatrix.
		template<typename typeB>
		SubMatrix & operator =(AbstractMatrix<typeB> const & other)
		{
			assert( Cols() == other.Cols() );
			assert( Rows() == other.Rows() );
			for(enumerator i = 0; i < other.Rows(); ++i)
				for(enumerator j = 0; j < other.Cols(); ++j)
					assign((*this)(i,j),other(i,j));
			return *this;
		}
		/// Assign submatrix of another type to submatrix.
		/// @param other Another submatrix of different type.
		/// @return Reference to current submatrix.
		template<typename typeB>
		SubMatrix & operator =(SubMatrix<typeB> const & other)
		{
			assert( Cols() == other.Cols() );
			assert( Rows() == other.Rows() );
			for(enumerator i = 0; i < other.Rows(); ++i)
				for(enumerator j = 0; j < other.Cols(); ++j)
					assign((*this)(i,j),other(i,j));
			return *this;
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		Var & operator()(enumerator i, enumerator j)
		{
			assert(i >= 0 && i < Rows());
			assert(j >= 0 && j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return (*M)(i+brow,j+bcol);
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		const Var & operator()(enumerator i, enumerator j) const
		{
			assert(i >= 0 && i < Rows());
			assert(j >= 0 && j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return (*M)(i+brow,j+bcol);
		}
		/// Convert submatrix into matrix.
		/// Note, that modifying returned matrix does
		/// not affect elements of the submatrix or original matrix
		/// used to create submatrix.
		/// @return Matrix with same entries as submatrix.
		::INMOST::Matrix<Var, pool_array_t<Var> > MakeMatrix()
		{
			::INMOST::Matrix<Var, pool_array_t<Var> > ret(Rows(),Cols());
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i,j) = (*this)(i,j);
			return ret;
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. SubMatrix cannot change it's size,
		/// since it just points to a part of the original matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
            (void)cols; (void)rows;
		}
	};
	
	/// This class allows for in-place operations on submatrix of the matrix elements.
	template<typename Var>
	class ConstSubMatrix : public AbstractMatrix<Var>
	{
	public:
		typedef unsigned enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrix<Var> * M;
		enumerator brow; //< First row in matrix M.
		enumerator erow; //< Last row in matrix M.
		enumerator bcol; //< First column in matrix M.
		enumerator ecol; //< Last column in matrix M.
	public:
		/// Number of rows in submatrix.
		/// @return Number of rows.
		enumerator Rows() const {return erow-brow;}
		/// Number of columns in submatrix.
		/// @return Number of columns.
		enumerator Cols() const {return ecol-bcol;}
		/// Create submatrix for a matrix.
		/// @param rM Reference to the matrix that stores elements.
		/// @param first_row First row in the matrix.
		/// @param last_row Last row in the matrix.
		/// @param first_column First column in the matrix.
		/// @param last_column Last column in the matrix.
		ConstSubMatrix(const AbstractMatrix<Var> & rM, enumerator first_row, enumerator last_row, enumerator first_column, enumerator last_column) : M(&rM), brow(first_row), erow(last_row), bcol(first_column), ecol(last_column)
		{}
		ConstSubMatrix(const ConstSubMatrix & b) : M(b.M), brow(b.brow), erow(b.erow), bcol(b.bcol), ecol(b.ecol) {}
		
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		Var & operator()(enumerator i, enumerator j)
		{
			assert(i >= 0 && i < Rows());
			assert(j >= 0 && j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return const_cast<Var &>((*M)(i+brow,j+bcol));
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		const Var & operator()(enumerator i, enumerator j) const
		{
			assert(i >= 0 && i < Rows());
			assert(j >= 0 && j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return (*M)(i+brow,j+bcol);
		}
		/// Convert submatrix into matrix.
		/// Note, that modifying returned matrix does
		/// not affect elements of the submatrix or original matrix
		/// used to create submatrix.
		/// @return Matrix with same entries as submatrix.
		::INMOST::Matrix<Var, pool_array_t<Var> > MakeMatrix()
		{
			::INMOST::Matrix<Var, pool_array_t<Var> > ret(Rows(),Cols());
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i,j) = (*this)(i,j);
			return ret;
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. SubMatrix cannot change it's size,
		/// since it just points to a part of the original matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
		}
	};
	
	template<typename Var>
	Matrix<Var, pool_array_t<Var> >
	AbstractMatrix<Var>::Transpose() const
	{
		Matrix<Var, pool_array_t<Var> > ret(Cols(),Rows());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(j,i) = (*this)(i,j);
		return ret;
	}
	
	template<typename Var>
	void
	AbstractMatrix<Var>::Swap(AbstractMatrix<Var> & b)
	{
		Matrix<Var, pool_array_t<Var> > tmp = b;
		b = (*this);
		(*this) = tmp;
	}
	
	template<typename Var>
	Matrix<Var, pool_array_t<Var> >
	AbstractMatrix<Var>::operator-() const
	{
		Matrix<Var, pool_array_t<Var> > ret(Rows(),Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = -(*this)(i,j);
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::CrossProduct(const AbstractMatrix<typeB> & B) const
	{
		const AbstractMatrix<Var> & A = *this;
		assert(A.Rows() == 3);
		assert(A.Cols() == 1);
		assert(B.Rows() == 3);
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(3,B.Cols()); //check RVO
		for(unsigned k = 0; k < B.Cols(); ++k)
		{
			ret(0,k) = (A(1,0)*B(2,k) - A(2,0)*B(1,k));
			ret(1,k) = (A(2,0)*B(0,k) - A(0,0)*B(2,k));
			ret(2,k) = (A(0,0)*B(1,k) - A(1,0)*B(0,k));
		}
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::Transform(const AbstractMatrix<typeB> & B) const
	{
		const AbstractMatrix<Var> & A = *this;
		assert(A.Rows() == 3);
		assert(A.Cols() == 1);
		assert(B.Rows() == 3);
		assert(B.Cols() == 1);
		typedef typename Promote<Var,typeB>::type var_t;
		Matrix<var_t, pool_array_t<var_t> > Q(3,3);
		var_t x,y,z,w,n,l;
		
		x = -(B(1,0)*A(2,0) - B(2,0)*A(1,0));
		y = -(B(2,0)*A(0,0) - B(0,0)*A(2,0));
		z = -(B(0,0)*A(1,0) - B(1,0)*A(0,0));
		w = A.FrobeniusNorm()*B.FrobeniusNorm() + A.DotProduct(B);
		
		l = B.FrobeniusNorm()/A.FrobeniusNorm();
		n = 2*l/(x*x+y*y+z*z+w*w);
		
		Q(0,0) = l - (y*y + z*z)*n;
		Q(0,1) = (x*y - z*w)*n;
		Q(0,2) = (x*z + y*w)*n;
		
		Q(1,0) = (x*y + z*w)*n;
		Q(1,1) = l - (x*x + z*z)*n;
		Q(1,2) = (y*z - x*w)*n;
		
		Q(2,0) = (x*z - y*w)*n;
		Q(2,1) = (y*z + x*w)*n;
		Q(2,2) = l - (x*x + y*y)*n;
		
		return Q;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::operator-(const AbstractMatrix<typeB> & other) const
	{
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(Rows(),Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = (*this)(i,j)-other(i,j);
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	AbstractMatrix<Var> &
	AbstractMatrix<Var>::operator-=(const AbstractMatrix<typeB> & other)
	{
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				assign((*this)(i,j),(*this)(i,j)-other(i,j));
		return *this;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::operator+(const AbstractMatrix<typeB> & other) const
	{
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(Rows(),Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = (*this)(i,j)+other(i,j);
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	AbstractMatrix<Var> &
	AbstractMatrix<Var>::operator+=(const AbstractMatrix<typeB> & other)
	{
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				assign((*this)(i,j),(*this)(i,j)+other(i,j));
		return *this;
	}
	
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::operator*(const AbstractMatrix<typeB> & other) const
	{
		assert(Cols() == other.Rows());
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				typename Promote<Var,typeB>::type tmp = 0.0;
				for(enumerator k = 0; k < Cols(); ++k)
					tmp += (*this)(i,k)*other(k,j);
				ret(i,j) = tmp;
			}
		}
		return ret;
	}
	
	
	/*
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var, typeB>::type>
	AbstractMatrix<Var>::operator*(const AbstractMatrix<typeB> & other) const
	{
		const enumerator b = 32;
		assert(Cols() == other.Rows());
		Matrix<typename Promote<Var, typeB>::type> ret(Rows(), other.Cols()); //check RVO
		for (enumerator i0 = 0; i0 < Rows(); i0+=b) //loop rows
		{
			for (enumerator j0 = 0; j0 < other.Cols(); j0 += b) //loop columns
			{
				enumerator i0e = std::min(i0+b,Rows());
				enumerator j0e = std::min(j0+b,other.Cols());
				for (enumerator i = i0; i < i0e; ++i)
					for(enumerator j = j0; j < j0e; ++j)
						ret(i,j) = 0.0;
				for (enumerator k0 = 0; k0 < Cols(); k0+=b)
				{
					enumerator k0e = std::min(k0+b,Cols());
					for (enumerator i = i0; i < i0e; ++i)
					{
						for (enumerator k = k0; k < k0e; ++k)
						{
							for(enumerator j = j0; j < j0e; ++j)
								assign(ret(i, j),ret(i, j)+(*this)(i, k)*other(k, j));
						}
					}
				}
			}
		}
		return ret;
	}
	*/
	
	template<typename Var>
	template<typename typeB>
	AbstractMatrix<Var> &
	AbstractMatrix<Var>::operator*=(const AbstractMatrix<typeB> & B)
	{
		(*this) = (*this)*B;
		return *this;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::operator*(typeB coef) const
	{
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(Rows(),Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = (*this)(i,j)*coef;
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	AbstractMatrix<Var> &
	AbstractMatrix<Var>::operator*=(typeB coef)
	{
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				assign((*this)(i,j),(*this)(i,j)*coef);
		return *this;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::operator/(typeB coef) const
	{
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(Rows(),Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = (*this)(i,j)/coef;
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	AbstractMatrix<Var> &
	AbstractMatrix<Var>::operator/=(typeB coef)
	{
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				assign((*this)(i,j),(*this)(i,j)/coef);
		return *this;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::Kronecker(const AbstractMatrix<typeB> & other) const
	{
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(Rows()*other.Rows(),Cols()*other.Cols());
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Rows(); ++j) //loop columns
			{
				for(enumerator k = 0; k < Cols(); ++k)
				{
					for(enumerator l = 0; l < other.Cols(); ++l)
					{
						ret(i*other.Rows()+j,k*other.Cols()+l) = other(j,l)*(*this)(i,k);
					}
				}
			}
		}
		return ret;
	}
	
	template<typename Var>
	Matrix<Var, pool_array_t<Var> >
	AbstractMatrix<Var>::Invert(int * ierr) const
	{
		Matrix<Var, pool_array_t<Var> > ret(Cols(),Rows());
		ret = Solve(Matrix<Var, pool_array_t<Var> >::Unit(Rows()),ierr);
		return ret;
	}
	
	template<typename Var>
	Matrix<Var, pool_array_t<Var> >
	AbstractMatrix<Var>::CholeskyInvert(int * ierr) const
	{
		Matrix<Var, pool_array_t<Var> > ret(Rows(),Rows());
		ret = CholeskySolve(Matrix<Var, pool_array_t<Var> >::Unit(Rows()),ierr);
		return ret;
	}
	
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::CholeskySolve(const AbstractMatrix<typeB> & B, int * ierr) const
	{
		const AbstractMatrix<Var> & A = *this;
		assert(A.Rows() == A.Cols());
		assert(A.Rows() == B.Rows());
		enumerator n = A.Rows();
		enumerator l = B.Cols();
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(B);
		SymmetricMatrix<Var, pool_array_t<Var> > L(A);
		
		//SAXPY
		/*
		for(enumerator i = 0; i < n; ++i)
		{
			for(enumerator k = 0; k < i; ++k)
				for(enumerator j = i; j < n; ++j)
					L(i,j) -= L(i,k)*L(j,k);
			if( L(i,i) < 0.0 )
			{
				ret.second = false;
				if( print_fail ) std::cout << "Negative diagonal pivot " << get_value(L(i,i)) << std::endl;
				return ret;
			}
			
			L(i,i) = sqrt(L(i,i));
			
			if( fabs(L(i,i)) < 1.0e-24 )
			{
				ret.second = false;
				if( print_fail ) std::cout << "Diagonal pivot is too small " << get_value(L(i,i)) << std::endl;
				return ret;
			}
			
			for(enumerator j = i+1; j < n; ++j)
				L(i,j) = L(i,j)/L(i,i);
		}
		*/
		//Outer product
		for(enumerator k = 0; k < n; ++k)
		{
			if( L(k,k) < 0.0 )
			{
				if( ierr )
				{
					if( *ierr == -1 ) std::cout << "Negative diagonal pivot " << get_value(L(k,k)) << " row " << k << std::endl;
					*ierr = k+1;
				}
				else throw MatrixCholeskySolveFail;
				return ret;
			}
			
			L(k,k) = sqrt(L(k,k));
			
			if( fabs(L(k,k)) < 1.0e-24 )
			{
				if( ierr )
				{
					if( *ierr == -1 ) std::cout << "Diagonal pivot is too small " << get_value(L(k,k)) << " row " << k << std::endl;
					*ierr = k+1;
				}
				else throw MatrixCholeskySolveFail;
				return ret;
			}
			
			for(enumerator i = k+1; i < n; ++i)
				L(i,k) = L(i,k)/L(k,k);
			
			for(enumerator j = k+1; j < n; ++j)
			{
				for(enumerator i = j; i < n; ++i)
					L(i,j) -= L(i,k)*L(j,k);
			}
		}
		// LY=B
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > & Y = ret;
		for(enumerator i = 0; i < n; ++i)
		{
			for(enumerator k = 0; k < l; ++k)
			{
				for(enumerator j = 0; j < i; ++j)
					Y(i,k) -= Y(j,k)*L(j,i);
				Y(i,k) /= L(i,i);
			}
		}
		// L^TX = Y
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > & X = ret;
		for(enumerator it = n; it > 0; --it)
		{
			enumerator i = it-1;
			for(enumerator k = 0; k < l; ++k)
			{
				for(enumerator jt = n; jt > it; --jt)
				{
					enumerator j = jt-1;
					X(i,k) -= X(j,k)*L(i,j);
				}
				X(i,k) /= L(i,i);
			}
		}
		if( ierr ) *ierr = 0;
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::Solve(const AbstractMatrix<typeB> & B, int * ierr) const
	{
		// for A^T B
		assert(Rows() == B.Rows());
		/*
		if( Rows() != Cols() )
		{
			Matrix<Var> At = this->Transpose(); //m by n matrix
			return (At*(*this)).Solve(At*B,print_fail);
		}
		Matrix<typeB> AtB = B; //m by l matrix
		Matrix<Var> AtA = (*this); //m by m matrix
		 */
		enumerator l = B.Cols();
		enumerator m = Cols();
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(m,l);
		Matrix<Var, pool_array_t<Var> > At = this->Transpose(); //m by n matrix
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > AtB = At*B; //m by l matrix
		Matrix<Var, pool_array_t<Var> > AtA = At*(*this); //m by m matrix
		assert(l == AtB.Cols());
		//enumerator l = AtB.Cols();
        //enumerator n = Rows();
		
		dynarray<enumerator,128> order(m);
		
		//Var temp;
		INMOST_DATA_REAL_TYPE max,v;
		//typeB tempb;
		for(enumerator i = 0; i < m; ++i) order[i] = i;
		for(enumerator i = 0; i < m; i++)
		{
			enumerator maxk = i, maxq = i;//, temp2;
			max = fabs(get_value(AtA(maxk,maxq)));
			//Find best pivot
			//if( max < 1.0e-8 )
			{
				for(enumerator k = i; k < m; k++) // over rows
				{
					for(enumerator q = i; q < m; q++) // over columns
					{
						v = fabs(get_value(AtA(k,q)));
						if( v > max )
						{
							max = v;
							maxk = k;
							maxq = q;
						}
					}
				}
				//Exchange rows
				if( maxk != i )
				{
					for(enumerator q = 0; q < m; q++) // over columns of A
					{
						std::swap(AtA(maxk,q),AtA(i,q));
						//temp = AtA(maxk,q);
						//AtA(maxk,q) = AtA(i,q);
						//AtA(i,q) = temp;
					}
					//exchange rhs
					for(enumerator q = 0; q < l; q++) // over columns of B
					{
						std::swap(AtB(maxk,q),AtB(i,q));
						//tempb = AtB(maxk,q);
						//AtB(maxk,q) = AtB(i,q);
						//AtB(i,q) = tempb;
					}
				}
				//Exchange columns
				if( maxq != i )
				{
					for(enumerator k = 0; k < m; k++) //over rows
					{
						std::swap(AtA(k,maxq),AtA(k,i));
						//temp = AtA(k,maxq);
						//AtA(k,maxq) = AtA(k,i);
						//AtA(k,i) = temp;
					}
					//remember order in sol
					{
						std::swap(order[maxq],order[i]);
						//temp2 = order[maxq];
						//order[maxq] = order[i];
						//order[i] = temp2;
					}
				}
			}
			//Check small entry
			if( fabs(get_value(AtA(i,i))) < 1.0e-54 )
			{
				bool ok = true;
				for(enumerator k = 0; k < l; k++) // over columns of B
				{
					if( fabs(get_value(AtB(i,k))/1.0e-54) > 1 )
					{
						ok = false;
						break;
					}
				}
				if( ok ) AtA(i,i) = AtA(i,i) < 0.0 ? - 1.0e-12 : 1.0e-12;
				else
				{
					if( ierr )
					{
						if( *ierr == -1 ) std::cout << "Failed to invert matrix" << std::endl;
						*ierr = i+1;
					}
					else throw MatrixSolveFail;
					return ret;
				}
			}
			//Divide row and column by diagonal values
			for(enumerator k = i+1; k < m; k++)
			{
				AtA(i,k) /= AtA(i,i);
				AtA(k,i) /= AtA(i,i);
			}
			//Elimination step for matrix
			for(enumerator k = i+1; k < m; k++)
				for(enumerator q = i+1; q < m; q++)
				{
					AtA(k,q) -= AtA(k,i) * AtA(i,i) * AtA(i,q);
				}
			//Elimination step for right hand side
			for(enumerator k = 0; k < l; k++)
			{
				for(enumerator j = i+1; j < m; j++) //iterate over columns of L
				{
					AtB(j,k) -= AtB(i,k) * AtA(j,i);
				}
				AtB(i,k) /= AtA(i,i);
			}
		}
		//Back substitution step
		for(enumerator k = 0; k < l; k++)
		{
			for(enumerator i = m; i-- > 0; ) //iterate over rows of U
				for(enumerator j = i+1; j < m; j++)
				{
					AtB(i,k) -= AtB(j,k) * AtA(i,j);
				}
			for(enumerator i = 0; i < m; i++)
				ret(order[i],k) = AtB(i,k);
		}
		//delete [] order;
		if( ierr ) *ierr = 0;
		return ret;
	}
	
	template<typename Var>
	Matrix<Var, pool_array_t<Var> >
	AbstractMatrix<Var>::ExtractSubMatrix(enumerator ibeg, enumerator iend, enumerator jbeg, enumerator jend) const
	{
		assert(ibeg < Rows());
		assert(iend < Rows());
		assert(jbeg < Cols());
		assert(jend < Cols());
		Matrix<Var, pool_array_t<Var> > ret(iend-ibeg,jend-jbeg);
		for(enumerator i = ibeg; i < iend; ++i)
		{
			for(enumerator j = jbeg; j < jend; ++j)
				ret(i-ibeg,j-jbeg) = (*this)(i,j);
		}
		return ret;
	}
	
	template<typename Var>
	Matrix<Var, pool_array_t<Var> >
	AbstractMatrix<Var>::Repack(enumerator rows, enumerator cols) const
	{
		assert(Cols()*Rows()==rows*cols);
		Matrix<Var, pool_array_t<Var> > ret(*this);
		ret.Rows() = rows;
		ret.Cols() = cols;
		return ret;
	}
	
	template<typename Var>
	Matrix<Var, pool_array_t<Var> >
	AbstractMatrix<Var>::PseudoInvert(INMOST_DATA_REAL_TYPE tol, int * ierr) const
	{
		Matrix<Var, pool_array_t<Var> > ret(Cols(),Rows());
		Matrix<Var, pool_array_t<Var> > U,S,V;
		bool success = SVD(U,S,V);
		if( !success )
		{
			if( ierr )
			{
				if( *ierr == -1 ) std::cout << "Failed to compute Moore-Penrose inverse of the matrix" << std::endl;
				*ierr = 1;
				return ret;
			}
			else throw MatrixPseudoSolveFail;
		}
        for(INMOST_DATA_ENUM_TYPE k = 0; k < S.Cols(); ++k)
		{
			if( S(k,k) > tol )
				S(k,k) = 1.0/S(k,k);
			else
				S(k,k) = 0.0;
		}
		ret = V*S*U.Transpose();
		if( ierr ) *ierr = 0;
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> >
	AbstractMatrix<Var>::PseudoSolve(const AbstractMatrix<typeB> & B, INMOST_DATA_REAL_TYPE tol, int * ierr) const
	{
		Matrix<typename Promote<Var,typeB>::type, pool_array_t<typename Promote<Var,typeB>::type> > ret(Cols(),B.Cols());
		Matrix<Var, pool_array_t<Var> > U,S,V;
		bool success = SVD(U,S,V);
		if( !success )
		{
			if( ierr )
			{
				if( *ierr == -1 ) std::cout << "Failed to compute Moore-Penrose inverse of the matrix" << std::endl;
				*ierr = 1;
			}
			else throw MatrixPseudoSolveFail;
		}
		for(int k = 0; k < S.Cols(); ++k)
		{
			if( S(k,k) > tol )
				S(k,k) = 1.0/S(k,k);
			else
				S(k,k) = 0.0;
		}
		ret = V*S*U.Transpose()*B;
		if( ierr ) *ierr = 0;
		return ret;
	}
	
	//template<typename Var,typename storage_type>
	//SubMatrix<Var,storage_type> Matrix<Var,storage_type>::SubMatrix(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col)
	//{
	//	return ::INMOST::SubMatrix<Var,storage_type>(*this,first_row,last_row,first_col,last_col);
	//}
	
	template<typename Var,typename storage_type>
	SubMatrix<Var> SymmetricMatrix<Var,storage_type>::operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col)
	{
		return ::INMOST::SubMatrix<Var>(*this,first_row,last_row,first_col,last_col);
	}
	template<typename Var, typename storage_type>
	ConstSubMatrix<Var> SymmetricMatrix<Var, storage_type>::operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const
	{
		return ::INMOST::ConstSubMatrix<Var>(*this, first_row, last_row, first_col, last_col);
	}
	
	template<typename Var,typename storage_type>
	SubMatrix<Var> Matrix<Var,storage_type>::operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col)
	{
		return ::INMOST::SubMatrix<Var>(*this,first_row,last_row,first_col,last_col);
	}
	template<typename Var, typename storage_type>
	ConstSubMatrix<Var> Matrix<Var, storage_type>::operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const
	{
		return ::INMOST::ConstSubMatrix<Var>(*this, first_row, last_row, first_col, last_col);
	}
	
	
	template<class T> struct make_integer;
	template<> struct make_integer<float> {typedef int type;};
	template<> struct make_integer<double> {typedef long long type;};

	__INLINE static bool compare(INMOST_DATA_REAL_TYPE * a, INMOST_DATA_REAL_TYPE * b)
	{
		return (*reinterpret_cast< make_integer<INMOST_DATA_REAL_TYPE>::type * >(a)) <=
			(*reinterpret_cast< make_integer<INMOST_DATA_REAL_TYPE>::type * >(b));
	}
	
	class BinaryHeapDense
	{
		INMOST_DATA_ENUM_TYPE size_max, size;
		std::vector<INMOST_DATA_ENUM_TYPE> heap;
		std::vector<INMOST_DATA_ENUM_TYPE> index;
		INMOST_DATA_REAL_TYPE * keys;
		
		void swap(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
		{
			INMOST_DATA_ENUM_TYPE t = heap[i];
			heap[i] = heap[j];
			heap[j] = t;
			index[heap[i]] = i;
			index[heap[j]] = j;
		}
		void bubbleUp(INMOST_DATA_ENUM_TYPE k)
		{
			while(k > 1 && keys[heap[k/2]] > keys[heap[k]])
			{
				swap(k, k/2);
				k = k/2;
			}
		}
		void bubbleDown(INMOST_DATA_ENUM_TYPE k)
		{
			INMOST_DATA_ENUM_TYPE j;
			while(2*k <= size)
			{
				j = 2*k;
				if(j < size && keys[heap[j]] > keys[heap[j+1]])
					j++;
				if(keys[heap[k]] <= keys[heap[j]])
					break;
				swap(k, j);
				k = j;
			}
		}
	public:
		BinaryHeapDense(INMOST_DATA_REAL_TYPE * pkeys, INMOST_DATA_ENUM_TYPE len)
		{
			size_max = len;
			keys = pkeys;
			size = 0;
			heap.resize(size_max+1);
			index.resize(size_max+1);
			for(INMOST_DATA_ENUM_TYPE i = 0; i <= size_max; i++)
				index[i] = ENUMUNDEF;
		}
		void PushHeap(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			size++;
			index[i] = size;
			heap[size] = i;
			keys[i] = key;
			bubbleUp(size);
		}
		INMOST_DATA_ENUM_TYPE PopHeap()
		{
			if( size == 0 ) return ENUMUNDEF;
			INMOST_DATA_ENUM_TYPE min = heap[1];
			swap(1, size--);
			bubbleDown(1);
			index[min] = ENUMUNDEF;
			heap[size+1] = ENUMUNDEF;
			return min;
		}
		void DecreaseKey(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			keys[i] = key;
			bubbleUp(index[i]);
		}
		void Clear()
		{
			while( size ) keys[PopHeap()] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
		}
		bool Contains(INMOST_DATA_ENUM_TYPE i)
		{
			return index[i] != ENUMUNDEF;
		}
	};
	
	template<typename Var>
	void AbstractMatrix<Var>::MPT(INMOST_DATA_ENUM_TYPE * Perm, INMOST_DATA_REAL_TYPE * SL, INMOST_DATA_REAL_TYPE * SR) const
	{
		const INMOST_DATA_ENUM_TYPE EOL = ENUMUNDEF-1;
		int n = Rows();
		int m = Cols();
		INMOST_DATA_REAL_TYPE u, l;
		array<INMOST_DATA_REAL_TYPE> Cmax(m,0.0);
		array<INMOST_DATA_REAL_TYPE> U(m,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		array<INMOST_DATA_REAL_TYPE> V(n,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		std::fill(Perm,Perm+m,ENUMUNDEF);
		//array<INMOST_DATA_ENUM_TYPE> Perm(m,ENUMUNDEF);
		array<INMOST_DATA_ENUM_TYPE> IPerm(std::max(n,m),ENUMUNDEF);
		array<INMOST_DATA_ENUM_TYPE> ColumnList(m,ENUMUNDEF);
		array<INMOST_DATA_ENUM_TYPE> ColumnPosition(n,ENUMUNDEF);
		array<INMOST_DATA_ENUM_TYPE> Parent(n,ENUMUNDEF);
		array<INMOST_DATA_REAL_TYPE> Dist(m,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		const AbstractMatrix<Var> & A = *this;
		Matrix<INMOST_DATA_REAL_TYPE> C(n,m);
		INMOST_DATA_ENUM_TYPE Li, Ui;
		INMOST_DATA_ENUM_TYPE ColumnBegin, PathEnd, Trace, IPermPrev;
		INMOST_DATA_REAL_TYPE ShortestPath, AugmentPath;
		BinaryHeapDense Heap(&Dist[0],m);
		//Initial LOG transformation to dual problem and initial extreme match
		for(int k = 0; k < n; ++k)
		{
			for(int i = 0; i < m; ++i)
			{
				C(k,i) = fabs(get_value(A(k,i)));
				if( Cmax[i] < C(k,i) ) Cmax[i] = C(k,i);
			}
		}
		for(int k = 0; k < n; ++k)
		{
			for (int i = 0; i < m; ++i)
			{
				if( Cmax[i] == 0 || C(k,i) == 0 )
					C(k,i) = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
				else
				{
					C(k,i) = log(Cmax[i])-log(C(k,i));
					if( C(k,i) < U[i] ) U[i] = C(k,i);
				}
			}
		}
		for(int k = 0; k < n; ++k)
		{
			for (int i = 0; i < m; ++i)
			{
				u = C(k,i) - U[i];
				if( u < V[k] ) V[k] = u;
			}
		}
		// Update cost and match 
		for(int k = 0; k < n; ++k)
		{
			for (int i = 0; i < m; ++i)
			{
				u = fabs(C(k,i) - V[k] - U[i]);
				if( u < 1.0e-30 && Perm[i] == ENUMUNDEF && IPerm[k] == ENUMUNDEF )
				{
					 Perm[i] = k;
					 IPerm[k] = i;
					 ColumnPosition[k] = i;
				}
			}
		}
		//1-step augmentation 
		for(int k = 0; k < n; ++k)
		{
			if( IPerm[k] == ENUMUNDEF ) //unmatched row
			{
				for (int i = 0; i < m && IPerm[k] == ENUMUNDEF; ++i)
				{
					u = fabs(C(k,i) - V[k] - U[i]);
					if( u <= 1.0e-30 )
					{
						Li = Perm[i];
						assert(Li != ENUMUNDEF);
						// Search other row in C for 0
						for (int Lit = 0; Lit < m; ++Lit)
						{
							u = fabs(C(Li,Lit) - V[Li] - U[Lit]);
							if( u <= 1.0e-30 && Perm[Lit] == ENUMUNDEF )
							{
								Perm[i] = k;
								IPerm[k] = i;
								ColumnPosition[k] = i;
								Perm[Lit] = Li;
								IPerm[Li] = Lit;
								ColumnPosition[Li] = Lit;
								break;
							}
						}
					}
				}
			}
		}
		// Weighted bipartite matching 
		for(int k = 0; k < n; ++k)
		{
			if( IPerm[k] != ENUMUNDEF )
				continue;
			Li = k;
			ColumnBegin = EOL;
			Parent[Li] = ENUMUNDEF;
			PathEnd = ENUMUNDEF;
			Trace = k;
			ShortestPath = 0;
			AugmentPath = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
			while(true)
			{
				for (int Lit = 0; Lit < m; ++Lit)
				{
					if( ColumnList[Lit] != ENUMUNDEF ) continue;
					l = ShortestPath + C(Li,Lit) - V[Li] - U[Lit];
					if( l < 0.0 && fabs(l) < 1.0e-12 ) l = 0;
					if( l < AugmentPath )
					{
						if( Perm[Lit] == ENUMUNDEF )
						{
							PathEnd = Lit;
							Trace = Li;
							AugmentPath = l;
						}
						else if( l < Dist[Lit] )
						{
							Parent[Perm[Lit]] = Li;
							if( Heap.Contains(Lit) )
								Heap.DecreaseKey(Lit,l);
							else 
								Heap.PushHeap(Lit,l);
						}
					}
				}

				INMOST_DATA_ENUM_TYPE pop_heap_pos = Heap.PopHeap();
				if( pop_heap_pos == ENUMUNDEF ) break;
			
				Ui = pop_heap_pos;
				ShortestPath = Dist[Ui];

				if( AugmentPath <= ShortestPath ) 
				{
					Dist[Ui] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					break;
				}


				ColumnList[Ui] = ColumnBegin;
				ColumnBegin = Ui;

				Li = Perm[Ui];
				
			}
			if( PathEnd != ENUMUNDEF )
			{
				Ui = ColumnBegin;
				while(Ui != EOL)
				{
					U[Ui] += Dist[Ui] - AugmentPath;
					if( Perm[Ui] != ENUMUNDEF ) V[Perm[Ui]] = C(Perm[Ui],ColumnPosition[Perm[Ui]]) - U[Ui];
					Dist[Ui] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					Li = ColumnList[Ui];
					ColumnList[Ui] = ENUMUNDEF;
					Ui = Li;
				}

				Ui = PathEnd;
				while(Trace != ENUMUNDEF)
				{
					IPermPrev = IPerm[Trace];
					Perm[Ui] = Trace;
					IPerm[Trace] = Ui;

					ColumnPosition[Trace] = Ui;
					V[Trace] = C(Trace,ColumnPosition[Trace]) - U[Ui];

					Ui = IPermPrev;
					Trace = Parent[Trace];

				}
				Heap.Clear();
			}
		}
		
		if( SL || SR )
		{
			for (int k = 0; k < std::min(n,m); ++k)
			{
				l = (V[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() ? 1 : exp(V[k]));
				u = (U[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || Cmax[k] == 0 ? 1 : exp(U[k])/Cmax[k]);
				
				if( l*get_value(A(k,Perm[k]))*u < 0.0 ) l *= -1;
				
				if( SL ) SL[Perm[k]] = u;
				if( SR ) SR[k] = l;
			}
			if( SR ) for(int k = std::min(n,m); k < n; ++k) SR[k] = 1;
			if( SL ) for(int k = std::min(n,m); k < m; ++k) SL[k] = 1;
		}
		//check that there are no gaps in Perm
		std::fill(IPerm.begin(),IPerm.end(),ENUMUNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> gaps;
		for(int k = 0; k < m; ++k)
		{
			if( Perm[k] != ENUMUNDEF )
				IPerm[Perm[k]] = 0;
		}
		for(int k = 0; k < m; ++k)
		{
			if( IPerm[k] == ENUMUNDEF )
				gaps.push_back(k);
		}
		for(int k = 0; k < m; ++k)
		{
			if( Perm[k] == ENUMUNDEF )
			{
				Perm[k] = gaps.back();
				gaps.pop_back();
			}
		}
	}
	
	/// shortcut for matrix of integer values.
	typedef Matrix<INMOST_DATA_INTEGER_TYPE> iMatrix;
	/// shortcut for matrix of real values.
	typedef Matrix<INMOST_DATA_REAL_TYPE> rMatrix;
	/// shortcut for matrix of integer values in stack-preallocated and dynamically reallocated array.
	typedef Matrix<INMOST_DATA_INTEGER_TYPE,dynarray<INMOST_DATA_INTEGER_TYPE,128> > idMatrix;
	/// shortcut for matrix of real values in stack-preallocated and dynamically reallocated array.
	typedef Matrix<INMOST_DATA_REAL_TYPE, dynarray<INMOST_DATA_REAL_TYPE,128> > rdMatrix;
	/// shortcut for matrix of integer values in pool-allocated array (beaware of deallocation order issue).
	typedef Matrix<INMOST_DATA_INTEGER_TYPE,pool_array_t<INMOST_DATA_INTEGER_TYPE> > ipMatrix;
	/// shortcut for matrix of real values in pool-allocated array (beaware of deallocation order issue).
	typedef Matrix<INMOST_DATA_REAL_TYPE, pool_array_t<INMOST_DATA_REAL_TYPE> > rpMatrix;
	/// shortcut for matrix of integer values in existing array.
	typedef Matrix<INMOST_DATA_INTEGER_TYPE,shell<INMOST_DATA_INTEGER_TYPE> > iaMatrix;
	/// shortcut for matrix of real values in existing array.
	typedef Matrix<INMOST_DATA_REAL_TYPE,shell<INMOST_DATA_REAL_TYPE> > raMatrix;
	/// shortcut for symmetric matrix of integer values.
	typedef SymmetricMatrix<INMOST_DATA_INTEGER_TYPE> iSymmetricMatrix;
	/// shortcut for symmetric matrix of real values.
	typedef SymmetricMatrix<INMOST_DATA_REAL_TYPE> rSymmetricMatrix;
	/// shortcut for symmetric matrix of integer values in stack-preallocated and dynamically reallocated array.
	typedef SymmetricMatrix<INMOST_DATA_INTEGER_TYPE,dynarray<INMOST_DATA_INTEGER_TYPE,128> > idSymmetricMatrix;
	/// shortcut for symmetric matrix of real values in stack-preallocated and dynamically reallocated array.
	typedef SymmetricMatrix<INMOST_DATA_REAL_TYPE,dynarray<INMOST_DATA_REAL_TYPE,128> > rdSymmetricMatrix;
	/// shortcut for symmetric matrix of integer values in pool-allocated array (beaware of deallocation order issue).
	typedef SymmetricMatrix<INMOST_DATA_INTEGER_TYPE,pool_array_t<INMOST_DATA_INTEGER_TYPE> > ipSymmetricMatrix;
	/// shortcut for symmetric matrix of real values in pool-allocated array (beaware of deallocation order issue).
	typedef SymmetricMatrix<INMOST_DATA_REAL_TYPE,pool_array_t<INMOST_DATA_REAL_TYPE> > rpSymmetricMatrix;
	/// shortcut for matrix of integer values in existing array.
	typedef SymmetricMatrix<INMOST_DATA_INTEGER_TYPE,shell<INMOST_DATA_INTEGER_TYPE> > iaSymmetricMatrix;
	/// shortcut for matrix of real values in existing array.
	typedef SymmetricMatrix<INMOST_DATA_REAL_TYPE,shell<INMOST_DATA_REAL_TYPE> > raSymmetricMatrix;
	/// return a matrix
	static iaMatrix iaMatrixMake(INMOST_DATA_INTEGER_TYPE * p, iaMatrix::enumerator n, iaMatrix::enumerator m) {return iaMatrix(shell<INMOST_DATA_INTEGER_TYPE>(p,n*m),n,m);}
	static raMatrix raMatrixMake(INMOST_DATA_REAL_TYPE * p, raMatrix::enumerator n, raMatrix::enumerator m) {return raMatrix(shell<INMOST_DATA_REAL_TYPE>(p,n*m),n,m);}
	static iaSymmetricMatrix iaSymmetricMatrixMake(INMOST_DATA_INTEGER_TYPE * p, iaSymmetricMatrix::enumerator n) {return iaSymmetricMatrix(shell<INMOST_DATA_INTEGER_TYPE>(p,n*(n+1)/2),n);}
	static raSymmetricMatrix raSymmetricMatrixMake(INMOST_DATA_REAL_TYPE * p, raSymmetricMatrix::enumerator n) {return raSymmetricMatrix(shell<INMOST_DATA_REAL_TYPE>(p,n*(n+1)/2),n);}
#if defined(USE_AUTODIFF)
	/// shortcut for matrix of variables with single unit entry of first order derivative.
	typedef Matrix<unknown> uMatrix;
	/// shortcut for matrix of variables with first order derivatives.
	typedef Matrix<variable> vMatrix;
	//< shortcut for matrix of variables with first and second order derivatives.
	typedef Matrix<hessian_variable> hMatrix;
	/// shortcut for matrix of variables with single unit entry of first order derivative in stack-preallocated and dynamically reallocated array.
	typedef Matrix<unknown, dynarray<unknown,128> > udMatrix;
	/// shortcut for matrix of variables with first order derivatives in stack-preallocated and dynamically reallocated array.
	typedef Matrix<variable, dynarray<variable,128> > vdMatrix;
	//< shortcut for matrix of variables with first and second order derivatives in stack-preallocated and dynamically reallocated array.
	typedef Matrix<hessian_variable, dynarray<hessian_variable,128> > hdMatrix;
	/// shortcut for matrix of variables with single unit entry of first order derivative in pool-allocated array (beaware of deallocation order issue).
	typedef Matrix<unknown, pool_array_t<unknown> > upMatrix;
	/// shortcut for matrix of variables with first order derivatives in pool-allocated array (beaware of deallocation order issue).
	typedef Matrix<variable, pool_array_t<variable> > vpMatrix;
	//< shortcut for matrix of variables with first and second order derivatives in pool-allocated array (beaware of deallocation order issue).
	typedef Matrix<hessian_variable, pool_array_t<hessian_variable> > hpMatrix;
	/// shortcut for matrix of unknowns in existing array.
	typedef Matrix<unknown,shell<unknown> > uaMatrix;
	/// shortcut for matrix of variables in existing array.
	typedef Matrix<variable,shell<variable> > vaMatrix;
	/// shortcut for matrix of variables in existing array.
	typedef Matrix<hessian_variable,shell<hessian_variable> > haMatrix;
	/// shortcut for matrix of unknowns in existing array.
	typedef SymmetricMatrix<unknown,shell<unknown> > uaSymmetricMatrix;
	/// shortcut for matrix of variables in existing array.
	typedef SymmetricMatrix<variable,shell<variable> > vaSymmetricMatrix;
	/// shortcut for matrix of variables in existing array.
	typedef SymmetricMatrix<hessian_variable,shell<hessian_variable> > haSymmetricMatrix;
	/// shortcut for symmetric matrix of variables with single unit entry of first order derivative in stack-preallocated and dynamically reallocated array.
	typedef SymmetricMatrix<unknown, dynarray<unknown,128> > udSymmetricMatrix;
	/// shortcut for symmetric matrix of variables with first order derivatives in stack-preallocated and dynamically reallocated array.
	typedef SymmetricMatrix<variable, dynarray<variable,128> > vdSymmetricMatrix;
	//< shortcut for symmetric matrix of variables with first and second order derivatives in stack-preallocated and dynamically reallocated array.
	typedef SymmetricMatrix<hessian_variable, dynarray<hessian_variable,128> > hdSymmetricMatrix;
	/// shortcut for symmetric matrix of variables with single unit entry of first order derivative in pool-allocated array (beaware of deallocation order issue).
	typedef SymmetricMatrix<unknown, pool_array_t<unknown> > upSymmetricMatrix;
	/// shortcut for symmetric matrix of variables with first order derivatives in pool-allocated array (beaware of deallocation order issue).
	typedef SymmetricMatrix<variable, pool_array_t<variable> > vpSymmetricMatrix;
	//< shortcut for symmetric matrix of variables with first and second order derivatives in pool-allocated array (beaware of deallocation order issue).
	typedef SymmetricMatrix<hessian_variable, pool_array_t<hessian_variable> > hpSymmetricMatrix;
	
	
	
	static uaMatrix uaMatrixMake(unknown * p, uaMatrix::enumerator n, uaMatrix::enumerator m) {return uaMatrix(shell<unknown>(p,n*m),n,m);}
	static vaMatrix vaMatrixMake(variable * p, vaMatrix::enumerator n, vaMatrix::enumerator m) {return vaMatrix(shell<variable>(p,n*m),n,m);}
	static haMatrix vaMatrixMake(hessian_variable * p, haMatrix::enumerator n, haMatrix::enumerator m) {return haMatrix(shell<hessian_variable>(p,n*m),n,m);}
	static uaSymmetricMatrix uaSymmetricMatrixMake(unknown * p, uaSymmetricMatrix::enumerator n) {return uaSymmetricMatrix(shell<unknown>(p,n*(n+1)/2),n);}
	static vaSymmetricMatrix vaSymmetricMatrixMake(variable * p, vaSymmetricMatrix::enumerator n) {return vaSymmetricMatrix(shell<variable>(p,n*(n+1)/2),n);}
	static haSymmetricMatrix vaSymmetricMatrixMake(hessian_variable * p, haSymmetricMatrix::enumerator n) {return haSymmetricMatrix(shell<hessian_variable>(p,n*(n+1)/2),n);}
#endif
}
/// Multiplication of matrix by constant from left.
/// @param coef Constant coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a constant.
template<typename typeB>
INMOST::Matrix<typename INMOST::Promote<INMOST_DATA_REAL_TYPE,typeB>::type, INMOST::pool_array_t<typename INMOST::Promote<INMOST_DATA_REAL_TYPE,typeB>::type> >
operator *(INMOST_DATA_REAL_TYPE coef, const INMOST::AbstractMatrix<typeB> & other)
{return other*coef;}
#if defined(USE_AUTODIFF)
/// Multiplication of matrix by a unknown from left.
/// Takes account for derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB>
INMOST::Matrix<typename INMOST::Promote<INMOST::unknown,typeB>::type, INMOST::pool_array_t<typename INMOST::Promote<INMOST::unknown,typeB>::type> >
operator *(const INMOST::unknown & coef, const INMOST::AbstractMatrix<typeB> & other)
{return other*coef;}
/// Multiplication of matrix by a variable from left.
/// Takes account for derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB>
INMOST::Matrix<typename INMOST::Promote<INMOST::variable,typeB>::type, INMOST::pool_array_t<typename INMOST::Promote<INMOST::variable,typeB>::type> >
operator *(const INMOST::variable & coef, const INMOST::AbstractMatrix<typeB> & other)
{return other*coef;}
/// Multiplication of matrix by a variable with first and
/// second order derivatives from left.
/// Takes account for first and second order derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB>
INMOST::Matrix<typename INMOST::Promote<INMOST::hessian_variable,typeB>::type, INMOST::pool_array_t<typename INMOST::Promote<INMOST::hessian_variable,typeB>::type> >
operator *(const INMOST::hessian_variable & coef, const INMOST::AbstractMatrix<typeB> & other)
{return other*coef;}
#endif

template<typename T>
__INLINE bool check_nans(const INMOST::AbstractMatrix<T> & A) {return A.CheckNans();}
template<typename T>
__INLINE bool check_infs(const INMOST::AbstractMatrix<T> & A) {return A.CheckInfs();}
template<typename T>
__INLINE bool check_nans_infs(const INMOST::AbstractMatrix<T> & A) {return A.CheckNansInfs();}

#endif //INMOST_DENSE_INCLUDED
