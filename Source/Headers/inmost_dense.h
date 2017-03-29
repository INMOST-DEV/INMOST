#ifndef INMOST_DENSE_INCLUDED
#define INMOST_DENSE_INCLUDED
#include "inmost_common.h"
#include "inmost_expression.h"
#include <iomanip>


namespace INMOST
{
	
	/// Structure that selects desired class, depending on the operation.
	template<class A, class B> struct Promote;
	template<> struct Promote<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
#if defined(USE_AUTODIFF)
	template<> struct Promote<INMOST_DATA_REAL_TYPE, variable>  {typedef variable type;};
	template<> struct Promote<variable, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
	template<> struct Promote<variable, variable> {typedef variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, hessian_variable>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, INMOST_DATA_REAL_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<variable, hessian_variable>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, variable>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, hessian_variable> {typedef hessian_variable type;};
#endif
	
	template<typename Var, typename Storage = array<Var> >
	class SubMatrix;
	
	template<typename Var, typename Storage = array<Var> >
	class Matrix;
	
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
		bool CheckNans()
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					if( check_nans((*this)(i,j)) ) return true;
			return false;
		}
		/// Singular value decomposition.
		/// Reconstruct matrix: A = U*Sigma*V.Transpose().
		/// Source http://www.public.iastate.edu/~dicook/JSS/paper/code/svd.c
		/// With few corrections for robustness.
		/// @param U Left unitary matrix, U^T U = I.
		/// @param Sigma Diagonal matrix with singular values.
		/// @param V Right unitary matrix, not transposed.
		/// @param order_singular_values
		/// \warning Somehow result differ in auto-differential mode.
		/// \todo Test implementation for auto-differentiation.
		bool SVD(AbstractMatrix<Var> & U, AbstractMatrix<Var> & Sigma, AbstractMatrix<Var> & V, bool order_singular_values = true) const
		{
			int flag, i, its, j, jj, k, l, nm;
			int n = Rows();
			int m = Cols();
			Var c, f, h, s, x, y, z;
			Var g = 0.0, scale = 0.0;
			INMOST_DATA_REAL_TYPE anorm = 0.0;
			if (n < m)
			{
				bool success = Transpose().SVD(U,Sigma,V);
				if( success )
				{
					U.Swap(V);
					U = U.Transpose();
					V = V.Transpose();
					return true;
				}
				else return false;
			} //m <= n
			array<Var> _rv1(m);
			shell<Var> rv1(_rv1);
			U = (*this);
			Sigma.Resize(m,m);
			Sigma.Zero();
			V.Resize(m,m);
			std::swap(n,m); //this how original algorithm takes it
			// Householder reduction to bidiagonal form
			for (i = 0; i < (int)n; i++)
			{
				// left-hand reduction
				l = i + 1;
				rv1[i] = scale * g;
				g = s = scale = 0.0;
				if (i < (int)m)
				{
					for (k = i; k < (int)m; k++) scale += fabs(U(k,i));
					if (get_value(scale))
					{
						for (k = i; k < (int)m; k++)
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
							for (j = l; j < (int)n; j++)
							{
								for (s = 0.0, k = i; k < (int)m; k++) s += U(k,i) * U(k,j);
								f = s / h;
								for (k = i; k < (int)m; k++) U(k,j) += f * U(k,i);
							}
						}
						for (k = i; k < (int)m; k++) U(k,i) *= scale;
					}
				}
				Sigma(i,i) = scale * g;
				// right-hand reduction
				g = s = scale = 0.0;
				if (i < (int)m && i != n - 1)
				{
					for (k = l; k < (int)n; k++) scale += fabs(U(i,k));
					if (get_value(scale))
					{
						for (k = l; k < (int)n; k++)
						{
							U(i,k) = U(i,k)/scale;
							s += U(i,k) * U(i,k);
						}
						f = U(i,l);
						g = -sign_func(sqrt(s), f);
						h = f * g - s;
						U(i,l) = f - g;
						for (k = l; k < (int)n; k++) rv1[k] = U(i,k) / h;
						if (i != m - 1)
						{
							for (j = l; j < (int)m; j++)
							{
								for (s = 0.0, k = l; k < (int)n; k++) s += U(j,k) * U(i,k);
								for (k = l; k < (int)n; k++) U(j,k) += s * rv1[k];
							}
						}
						for (k = l; k < (int)n; k++) U(i,k) *= scale;
					}
				}
				anorm = max_func(anorm,fabs(get_value(Sigma(i,i))) + fabs(get_value(rv1[i])));
			}
			
			// accumulate the right-hand transformation
			for (i = n - 1; i >= 0; i--)
			{
				if (i < (int)(n - 1))
				{
					if (get_value(g))
					{
						for (j = l; j < (int)n; j++) V(j,i) = ((U(i,j) / U(i,l)) / g);
						// double division to avoid underflow
						for (j = l; j < (int)n; j++)
						{
							for (s = 0.0, k = l; k < (int)n; k++) s += U(i,k) * V(k,j);
							for (k = l; k < (int)n; k++) V(k,j) += s * V(k,i);
						}
					}
					for (j = l; j < (int)n; j++) V(i,j) = V(j,i) = 0.0;
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
				if (i < (int)(n - 1))
					for (j = l; j < (int)n; j++)
						U(i,j) = 0.0;
				if (get_value(g))
				{
					g = 1.0 / g;
					if (i != n - 1)
					{
						for (j = l; j < (int)n; j++)
						{
							for (s = 0.0, k = l; k < (int)m; k++) s += (U(k,i) * U(k,j));
							f = (s / U(i,i)) * g;
							for (k = i; k < (int)m; k++) U(k,j) += f * U(k,i);
						}
					}
					for (j = i; j < (int)m; j++) U(j,i) = U(j,i)*g;
				}
				else for (j = i; j < (int)m; j++) U(j,i) = 0.0;
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
								for (j = 0; j < (int)m; j++)
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
						if (z < 0.0)
						{// make singular value nonnegative
							Sigma(k,k) = -z;
							for (j = 0; j < (int)n; j++) V(j,k) = -V(j,k);
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
						for (jj = 0; jj < (int)n; jj++)
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
						for (jj = 0; jj < (int)m; jj++)
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
				for(i = 0; i < (int)n; i++)
				{
					k = i;
					for(j = i+1; j < (int)n; ++j)
						if( Sigma(k,k) < Sigma(j,j) ) k = j;
					Var temp;
					if( Sigma(k,k) > Sigma(i,i) )
					{
						temp       = Sigma(k,k);
						Sigma(k,k) = Sigma(i,i);
						Sigma(i,i) = temp;
						// U is m by n
						for(int j = 0; j < (int)m; ++j)
						{
							temp   = U(j,k);
							U(j,k) = U(j,i);
							U(j,i) = temp;
						}
						// V is n by n
						for(int j = 0; j < (int)n; ++j)
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
		Matrix<Var> Transpose() const;
		///Exchange contents of two matrices.
		virtual void Swap(AbstractMatrix<Var> & b);
		/// Unary minus. Change sign of each element of the matrix.
		Matrix<Var> operator-() const;
		/// Subtract a matrix.
		/// @param other The matrix to be subtracted.
		/// @return Difference of current matrix with another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type> operator-(const AbstractMatrix<typeB> & other) const;
		/// Subtract a matrix and store result in the current.
		/// @param other The matrix to be subtracted.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator-=(const AbstractMatrix<typeB> & other);
		/// Add a matrix.
		/// @param other The matrix to be added.
		/// @return Sum of current matrix with another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type> operator+(const AbstractMatrix<typeB> & other) const;
		/// Add a matrix and store result in the current.
		/// @param other The matrix to be added.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator+=(const AbstractMatrix<typeB> & other);
		/// Multiply the matrix by another matrix.
		/// @param other Matrix to be multiplied from right.
		/// @return Matrix multipled by another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type> operator*(const AbstractMatrix<typeB> & other) const;
		/// Multiply matrix with another matrix in-place.
		/// @param B Another matrix to the right in multiplication.
		/// @return Reference to current matrix.
		template<typename typeB>
		AbstractMatrix & operator*=(const AbstractMatrix<typeB> & B);
		/// Multiply the matrix by a coefficient.
		/// @param coef Coefficient.
		/// @return Matrix multiplied by the coefficient.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type> operator*(typeB coef) const;
		/// Multiply the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator*=(typeB coef);
		/// Divide the matrix by a coefficient of a different type.
		/// @param coef Coefficient.
		/// @return Matrix divided by the coefficient.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type>
		operator/(typeB coef) const;
		/// Divide the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix &
		operator/=(typeB coef);
		/// Performs B^{-1}*A, multiplication by inverse matrix from left.
		/// @param other Matrix to be inverted and multiplied from left.
		/// @return Multiplication of current matrix by inverse of another
		/// matrix from left and boolean.
		/// If boolean is true, then the matrix was inverted successfully.
		/// \todo (test) Use Solve here.
		template<typename typeB>
		std::pair<Matrix<typename Promote<Var,typeB>::type>,bool>
		operator/(const AbstractMatrix<typeB> & other) const
		{
			return other.Solve(*this);
		}
		/// Kronecker product, latex symbol \otimes.
		/// @param other Matrix on the right of the Kronecker product.
		/// @return Result of the Kronecker product of current and another matrix.
		template<typename typeB>
		Matrix<typename Promote<Var,typeB>::type>
		Kronecker(const AbstractMatrix<typeB> & other) const;
		/// Inverts matrix using Crout-LU decomposition with full pivoting for
		/// maximum element. Works well for non-singular matrices, for singular
		/// matrices look see Matrix::PseudoInvert.
		/// @param print_fail Verbose output for singular matrices.
		/// @return A pair of inverse matrix and boolean. If boolean is true,
		/// then the matrix was inverted successfully.
		/// @see Matrix::PseudoInvert.
		/// \todo (test) Activate and test implementation with Solve.
		std::pair<Matrix<Var>,bool> Invert(bool print_fail = false) const;
		/// Finds X in A*X=B, where A and B are matrices.
		/// Converts system into A^T*A*X=A^T*B.
		/// Inverts matrix A^T*A using Crout-LU decomposition with full pivoting for
		/// maximum element. Works well for non-singular matrices, for singular
		/// matrices see Matrix::PseudoInvert.
		/// @param print_fail Verbose output for singular matrices.
		/// @return A pair of inverse matrix and boolean. If boolean is true,
		/// then the matrix was inverted successfully.
		/// @see Matrix::PseudoInvert.
		/// \warning Number of rows in matrices A and B should match.
		/// \todo
		/// 1. Test implementation.
		/// 2. Maximum product transversal + block pivoting instead of pivoting
		///    by maximum element.
		template<typename typeB>
		std::pair<Matrix<typename Promote<Var,typeB>::type>,bool>
		Solve(const AbstractMatrix<typeB> & B, bool print_fail = false) const;
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
					if( fabs(get_value((*this)(k,l))) > threshold )
						std::cout << std::setw(10) << get_value((*this)(k,l));
					else
						std::cout << std::setw(10) << 0;
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
			for(enumerator i = 0; i < Cols(); ++i)
				for(enumerator j = 0; j < Rows(); ++j)
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
		///                      then threshold are considered to be zero.
		/// @param print_fail Verbose output if singular value decomposition
		///                   algorithm has failed.
		/// @return A pair of pseudo-inverse matrix and boolean. If boolean is true,
		///         then the matrix was inverted successfully.
		std::pair<Matrix<Var>,bool> PseudoInvert(INMOST_DATA_REAL_TYPE tol = 0, bool print_fail = false) const;
		/// Solves the system of equations of the form A*X=B, with A and B matrices.
		/// Uses Moore-Penrose pseudo-inverse of the matrix A and calculates X = A^+*B.
		/// @param B Matrix at the right hand side.
		/// @param tol Thershold for singular values. Singular values smaller
		///                      then threshold are considered to be zero.
		/// @param print_fail Verbose output if singular value decomposition
		///                   algorithm has failed.
		/// @return A pair of the solution matrix X and boolean. If boolean is true,
		///         then the matrix was inverted successfully.
		template<typename typeB>
		std::pair<Matrix<typename Promote<Var,typeB>::type>,bool>
		PseudoSolve(const AbstractMatrix<typeB> & B, INMOST_DATA_REAL_TYPE tol = 0, bool print_fail = false) const;
		/// Extract submatrix of a matrix.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param ibeg Starting row in the original matrix.
		/// @param iend Last row (excluded) in the original matrix.
		/// @param jbeg Starting column in the original matrix.
		/// @param jend Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		Matrix<Var> ExtractSubMatrix(enumerator ibeg, enumerator iend, enumerator jbeg, enumerator jend) const;
		/// Change representation of the matrix into matrix of another size.
		/// Useful to change representation from matrix into vector and back.
		/// Replaces original number of columns and rows with a new one.
		/// @return Matrix with same entries and provided number of rows and columns.
		Matrix<Var> Repack(enumerator rows, enumerator cols) const;
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
			if( dynamic_cast<Matrix<Var,storage_type> *>(&b) != NULL )
			{
				Matrix<Var,storage_type> & bb = dynamic_cast<Matrix<Var,storage_type> &>(b);
				space.swap(bb.space);
				std::swap(n,bb.n);
				std::swap(m,bb.m);
			}
			else AbstractMatrix<Var>::Swap(b);
		}
		/// Construct empty matrix.
		Matrix() : space(), n(0), m(0) {}
		/// Construct the matrix from provided array and sizes.
		/// @param pspace Array of elements of the matrix, stored in row-wise format.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		Matrix(const Var * pspace, enumerator pn, enumerator pm) : space(pspace,pspace+pn*pm), n(pn), m(pm) {}
		/// Construct the matrix with the provided storage.
		/// Could be used to wrap existing array.
		/// @param pspace Storage of elements of the matrix, stored in row-wise format.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// \todo Do we need reference for pspace or just pspace?
		Matrix(storage_type pspace, enumerator pn, enumerator pm) : space(pspace), n(pn), m(pm) {}
		/// Construct a matrix with provided sizes.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// \warning The matrix does not necessery have zero entries.
		Matrix(enumerator pn, enumerator pm) : space(pn*pm), n(pn), m(pm) {}
		/// Copy matrix.
		/// @param other Another matrix of the same type.
		Matrix(const Matrix & other) : space(other.n*other.m), n(other.n), m(other.m)
		{
			for(enumerator i = 0; i < n*m; ++i)
				space[i] = other.space[i];
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
			if( n*m != other.n*other.m )
				space.resize(other.n*other.m);
			for(enumerator i = 0; i < other.n*other.m; ++i)
				space[i] = other.space[i];
			n = other.n;
			m = other.m;
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
		/// Convert values in array into square matrix.
		/// Supports the following representation, depending on the size
		/// of input array and size of side of final tensors' matrix:
		///
		/// representation | (array size, tensor size)
		///
		/// scalar         | (1,1), (1,2), (1,3), (1,6)
		///
		/// diagonal       | (2,2), (2,2), (3,3), (6,6)
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
		static Matrix<Var> FromTensor(Var * K, enumerator size, enumerator matsize = 3)
		{
			Matrix<Var> Kc(matsize,matsize);
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
		static Matrix FromVector(Var * r, enumerator size)
		{
			return Matrix(r,size,1);
		}
		/// Create diagonal matrix from array
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by array, other elements are zero.
		static Matrix FromDiagonal(Var * r, enumerator size)
		{
			Matrix ret(size,size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
			return ret;
		}
		/// Create diagonal matrix from array of values that have to be inversed.
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by inverse of array elements.
		static Matrix FromDiagonalInverse(Var * r, enumerator size)
		{
			Matrix ret(size,size);
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
		static Matrix CrossProduct(Var vec[3])
		{
			// |  0  -z   y |
			// |  z   0  -x |
			// | -y   x   0 |
			Matrix ret(3,3);
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
		/// and fills the diagonal with ones.
		/// @param pn Number of rows and columns in the matrix.
		/// @return Returns a unit matrix.
		static Matrix Unit(enumerator pn)
		{
			Matrix ret(pn,pn);
			ret.Zero();
			for(enumerator i = 0; i < pn; ++i) ret(i,i) = 1.0;
			return ret;
		}
		/// Concatenate B matrix as columns of current matrix.
		/// Assumes that number of rows of current matrix is
		/// equal to number of rows of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatRows
		Matrix ConcatCols(const Matrix & B)
		{
			assert(Rows() == B.Rows());
			Matrix ret(Rows(),Cols()+B.Cols());
			Matrix & A = *this;
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
		Matrix ConcatRows(const Matrix & B)
		{
			assert(Cols() == B.Cols());
			Matrix ret(Rows()+B.Rows(),Cols());
			Matrix & A = *this;
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
		Matrix JointDiagonalization(INMOST_DATA_REAL_TYPE threshold = 1.0e-7)
		{
			enumerator N = Rows();
			enumerator M = Cols() / Rows();
			Matrix V = Matrix::Unit(m);
			Matrix R(2,M);
			Matrix G(2,2);
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
		SubMatrix<Var,storage_type> MakeSubMatrix(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);
	};
	/// This class allows for in-place operations on submatrix of the matrix elements.
	template<typename Var, typename Storage>
	class SubMatrix : public AbstractMatrix<Var>
	{
	public:
		typedef unsigned enumerator; //< Integer type for indexes.
	private:
		Matrix<Var,Storage> * M;
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
		SubMatrix(Matrix<Var,Storage> & rM, enumerator first_row, enumerator last_row, enumerator first_column, enumerator last_column) : M(&rM), brow(first_row), erow(last_row), bcol(first_column), ecol(last_column)
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
		template<typename typeB,typename storageB>
		SubMatrix & operator =(SubMatrix<typeB,storageB> const & other)
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
		Matrix<Var> MakeMatrix()
		{
			Matrix<Var> ret(Rows(),Cols());
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
		}
	};
	
	template<typename Var>
	Matrix<Var>
	AbstractMatrix<Var>::Transpose() const
	{
		Matrix<Var> ret(Cols(),Rows());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(j,i) = (*this)(i,j);
		return ret;
	}
	
	template<typename Var>
	void
	AbstractMatrix<Var>::Swap(AbstractMatrix<Var> & b)
	{
		Matrix<Var> tmp = b;
		b = (*this);
		(*this) = tmp;
	}
	
	template<typename Var>
	Matrix<Var>
	AbstractMatrix<Var>::operator-() const
	{
		Matrix<Var> ret(Rows(),Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = -(*this)(i,j);
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrix<Var>::operator-(const AbstractMatrix<typeB> & other) const
	{
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols()); //check RVO
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
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrix<Var>::operator+(const AbstractMatrix<typeB> & other) const
	{
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols()); //check RVO
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
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrix<Var>::operator*(const AbstractMatrix<typeB> & other) const
	{
		assert(Cols() == other.Rows());
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),other.Cols()); //check RVO
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
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrix<Var>::operator*(typeB coef) const
	{
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols()); //check RVO
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
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrix<Var>::operator/(typeB coef) const
	{
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols());
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
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrix<Var>::Kronecker(const AbstractMatrix<typeB> & other) const
	{
		Matrix<typename Promote<Var,typeB>::type> ret(Rows()*other.Rows(),Cols()*other.Cols());
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
	std::pair<Matrix<Var>,bool>
	AbstractMatrix<Var>::Invert(bool print_fail) const
	{
		return Solve(Matrix<Var>::Unit(Rows()),print_fail);
	}
	
	template<typename Var>
	template<typename typeB>
	std::pair<Matrix<typename Promote<Var,typeB>::type>,bool>
	AbstractMatrix<Var>::Solve(const AbstractMatrix<typeB> & B, bool print_fail) const
	{
		// for A^T B
		assert(Rows() == B.Rows());
		Matrix<Var> At = this->Transpose(); //m by n matrix
		Matrix<typename Promote<Var,typeB>::type> AtB = At*B; //m by l matrix
		Matrix<Var> AtA = At*(*this); //m by m matrix
		enumerator l = AtB.Cols();
		enumerator n = Rows();
		enumerator m = Cols();
		enumerator * order = new enumerator [m];
		std::pair<Matrix<typename Promote<Var,typeB>::type>,bool>
		ret = std::make_pair(Matrix<typename Promote<Var,typeB>::type>(m,l),true);
		for(enumerator i = 0; i < m; ++i) order[i] = i;
		for(enumerator i = 0; i < m; i++)
		{
			enumerator maxk = i, maxq = i, temp2;
			Var max, temp;
			max = fabs(AtA(maxk,maxq));
			//Find best pivot
			for(enumerator k = i; k < m; k++) // over rows
			{
				for(enumerator q = i; q < m; q++) // over columns
				{
					if( fabs(AtA(k,q)) > max )
					{
						max = fabs(AtA(k,q));
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
					temp = AtA(maxk,q);
					AtA(maxk,q) = AtA(i,q);
					AtA(i,q) = temp;
				}
				//exchange rhs
				for(enumerator q = 0; q < l; q++) // over columns of B
				{
					temp = AtB(maxk,q);
					AtB(maxk,q) = AtB(i,q);
					AtB(i,q) = temp;
				}
			}
			//Exchange columns
			if( maxq != i )
			{
				for(enumerator k = 0; k < m; k++) //over rows
				{
					temp = AtA(k,maxq);
					AtA(k,maxq) = AtA(k,i);
					AtA(k,i) = temp;
				}
				//remember order in sol
				{
					temp2 = order[maxq];
					order[maxq] = order[i];
					order[i] = temp2;
				}
			}
			//Check small entry
			if( fabs(AtA(i,i)) < 1.0e-54 )
			{
				bool ok = true;
				for(enumerator k = 0; k < l; k++) // over columns of B
				{
					if( fabs(AtB(i,k)/1.0e-54) > 1 )
					{
						ok = false;
						break;
					}
				}
				if( ok ) AtA(i,i) = AtA(i,i) < 0.0 ? - 1.0e-12 : 1.0e-12;
				else
				{
					if( print_fail ) std::cout << "Failed to invert matrix" << std::endl;
					ret.second = false;
					delete [] order;
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
				ret.first(order[i],k) = AtB(i,k);
		}
		delete [] order;
		return ret;
	}
	
	template<typename Var>
	Matrix<Var>
	AbstractMatrix<Var>::ExtractSubMatrix(enumerator ibeg, enumerator iend, enumerator jbeg, enumerator jend) const
	{
		assert(ibeg < Rows());
		assert(iend < Rows());
		assert(jbeg < Cols());
		assert(jend < Cols());
		Matrix<Var> ret(iend-ibeg,jend-jbeg);
		for(enumerator i = ibeg; i < iend; ++i)
		{
			for(enumerator j = jbeg; j < jend; ++j)
				ret(i-ibeg,j-jbeg) = (*this)(i,j);
		}
		return ret;
	}
	
	template<typename Var>
	Matrix<Var>
	AbstractMatrix<Var>::Repack(enumerator rows, enumerator cols) const
	{
		assert(Cols()*Rows()==rows*cols);
		Matrix<Var> ret(*this);
		ret.Rows() = rows;
		ret.Cols() = cols;
		return ret;
	}
	
	template<typename Var>
	std::pair<Matrix<Var>,bool>
	AbstractMatrix<Var>::PseudoInvert(INMOST_DATA_REAL_TYPE tol, bool print_fail) const
	{
		std::pair<Matrix<Var>,bool> ret = std::make_pair(Matrix<Var>(Cols(),Rows()),true);
		Matrix<Var> U,S,V;
		ret.second = SVD(U,S,V);
		if( print_fail && !ret.second )
			std::cout << "Failed to compute Moore-Penrose inverse of the matrix" << std::endl;
		for(int k = 0; k < S.Cols(); ++k)
		{
			if( S(k,k) > tol )
				S(k,k) = 1.0/S(k,k);
			else
				S(k,k) = 0.0;
		}
		ret.first = V*S*U.Transpose();
		return ret;
	}
	
	template<typename Var>
	template<typename typeB>
	std::pair<Matrix<typename Promote<Var,typeB>::type>,bool>
	AbstractMatrix<Var>::PseudoSolve(const AbstractMatrix<typeB> & B, INMOST_DATA_REAL_TYPE tol, bool print_fail) const
	{
		std::pair<Matrix<typename Promote<Var,typeB>::type>,bool>
		ret = std::make_pair(Matrix<typename Promote<Var,typeB>::type>(Cols(),B.Cols()),true);
		Matrix<Var> U,S,V;
		ret.second = SVD(U,S,V);
		if( print_fail && !ret.second )
			std::cout << "Failed to compute Moore-Penrose inverse of the matrix" << std::endl;
		for(int k = 0; k < S.Cols(); ++k)
		{
			if( S(k,k) > tol )
				S(k,k) = 1.0/S(k,k);
			else
				S(k,k) = 0.0;
		}
		ret.first = V*S*U.Transpose()*B;
		return ret;
	}
	
	template<typename Var,typename storage_type>
	SubMatrix<Var,storage_type> Matrix<Var,storage_type>::MakeSubMatrix(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col)
	{
		return SubMatrix<Var,storage_type>(*this,first_row,last_row,first_col,last_col);
	}
	/// shortcut for matrix of real values.
	typedef Matrix<INMOST_DATA_REAL_TYPE> rMatrix;
#if defined(USE_AUTODIFF)
	/// shortcut for matrix of variables with first order derivatives.
	typedef Matrix<variable> vMatrix;
	//< shortcut for matrix of variables with first and second order derivatives.
	typedef Matrix<hessian_variable> hMatrix;
#endif
}
/// Multiplication of matrix by constant from left.
/// @param coef Constant coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a constant.
template<typename typeB,typename storageB>
INMOST::Matrix<typename INMOST::Promote<INMOST_DATA_REAL_TYPE,typeB>::type> operator *(INMOST_DATA_REAL_TYPE coef, const INMOST::Matrix<typeB,storageB> & other)
{return other*coef;}
#if defined(USE_AUTODIFF)
/// Multiplication of matrix by a variable from left.
/// Takes account for derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB,typename storageB>
INMOST::Matrix<typename INMOST::Promote<INMOST::variable,typeB>::type> operator *(const INMOST::variable & coef, const INMOST::Matrix<typeB,storageB> & other)
{return other*coef;}
/// Multiplication of matrix by a variable with first and
/// second order derivatives from left.
/// Takes account for first and second order derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB,typename storageB>
INMOST::Matrix<typename INMOST::Promote<INMOST::hessian_variable,typeB>::type> operator *(const INMOST::hessian_variable & coef, const INMOST::Matrix<typeB,storageB> & other)
{return other*coef;}
#endif


#endif //INMOST_DENSE_INCLUDED
