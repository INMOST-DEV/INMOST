#ifndef INMOST_DENSE_INCLUDED
#define INMOST_DENSE_INCLUDED
#include "inmost_common.h"
#include "inmost_expression.h"
#include <iomanip>
#include <stdarg.h>
#include <type_traits>

namespace INMOST
{
	
	/// Structure that selects desired class, depending on the operation.
	template<class A> struct SelfPromote;
	template<class A> struct ComplexType;
	template<class A, class B> struct Promote;
	template<> struct ComplexType<INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<> struct ComplexType<INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct ComplexType<INMOST_DATA_CPLX_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<typename T> inline typename ComplexType<T>::type real_part(T const& val) { return std::real(val); }
	//template<typename T> typename ComplexType<T>::type imag_part(T const& val) { return std::imag(val); }
	//template<> typename ComplexType<INMOST_DATA_INTEGER_TYPE>::type real_part(INMOST_DATA_INTEGER_TYPE const& a) { return a; }
	//template<> typename ComplexType<INMOST_DATA_REAL_TYPE>::type    real_part(INMOST_DATA_REAL_TYPE const& a) { return a; }
	//template<> typename ComplexType<INMOST_DATA_CPLX_TYPE>::type    real_part(INMOST_DATA_CPLX_TYPE const& a) { return a.real(); }
	//template<> typename ComplexType<INMOST_DATA_INTEGER_TYPE>::type  imag_part(INMOST_DATA_INTEGER_TYPE const& a) { return 0.0; }
	//template<> typename ComplexType<INMOST_DATA_REAL_TYPE>::type     imag_part(INMOST_DATA_REAL_TYPE const& a) { return 0.0; }
	//template<> typename ComplexType<INMOST_DATA_CPLX_TYPE>::type     imag_part(INMOST_DATA_CPLX_TYPE const& a) { return a.imag(); }
	template<> struct SelfPromote<INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<> struct SelfPromote<INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct SelfPromote<INMOST_DATA_CPLX_TYPE> { typedef INMOST_DATA_CPLX_TYPE type; };
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_INTEGER_TYPE> {typedef INMOST_DATA_INTEGER_TYPE type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE>    {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_CPLX_TYPE>    {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, INMOST_DATA_INTEGER_TYPE>    {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>       {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, INMOST_DATA_CPLX_TYPE>       {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, INMOST_DATA_INTEGER_TYPE>    {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, INMOST_DATA_REAL_TYPE>       {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, INMOST_DATA_CPLX_TYPE>       {typedef INMOST_DATA_CPLX_TYPE type;};
#if defined(USE_FP64)
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, float> {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<float, INMOST_DATA_INTEGER_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, float>    {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<float, INMOST_DATA_REAL_TYPE>    {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, float>    {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<float, INMOST_DATA_CPLX_TYPE>    {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<float, float> { typedef float type; };
#else
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, double> {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<double, INMOST_DATA_INTEGER_TYPE> {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, double>    {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<double, INMOST_DATA_REAL_TYPE>    {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, double>    {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<double, INMOST_DATA_CPLX_TYPE>    {typedef INMOST_DATA_CPLX_TYPE type;};
	template<> struct Promote<double, double> { typedef double type; };
#endif
#if defined(USE_AUTODIFF)
	template<> struct SelfPromote<unknown> { typedef variable type; };
	template<> struct SelfPromote<variable> { typedef variable type; };
	template<> struct SelfPromote<value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct SelfPromote<multivar_expression_reference> { typedef variable type; };
	template<> struct SelfPromote<hessian_multivar_expression_reference> { typedef hessian_variable type; };
	template<> struct SelfPromote<hessian_variable> { typedef hessian_variable type; };
	template<> struct ComplexType<unknown> { typedef unknown type; };
	template<> struct ComplexType<variable> { typedef variable type; };
	template<> struct ComplexType<value_reference> { typedef value_reference type; };
	template<> struct ComplexType<multivar_expression_reference> { typedef variable type; };
	template<> struct ComplexType<hessian_multivar_expression_reference> { typedef hessian_variable type; };
	template<> struct ComplexType<hessian_variable> { typedef hessian_variable type; };
	template<> struct ComplexType< std::complex<unknown> > { typedef unknown type; };
	template<> struct ComplexType< std::complex<variable> > { typedef variable type; };
	template<> struct ComplexType< std::complex<value_reference> > { typedef value_reference type; };
	template<> struct ComplexType< std::complex<multivar_expression_reference> > { typedef multivar_expression_reference type; };
	template<> struct ComplexType< std::complex<hessian_multivar_expression_reference> > { typedef hessian_multivar_expression_reference type; };
	template<> struct ComplexType< std::complex<hessian_variable> > { typedef hessian_variable type; };
	template<> inline typename ComplexType<unknown>::type                                                real_part(unknown const& a) { return a; }
	template<> inline typename ComplexType<variable>::type                                               real_part(variable const& a) { return a; }
	template<> inline typename ComplexType<value_reference>::type                                        real_part(value_reference const& a) { return a; }
	template<> inline typename ComplexType<multivar_expression_reference>::type                          real_part(multivar_expression_reference const& a) { return a; }
	template<> inline typename ComplexType<hessian_multivar_expression_reference>::type                  real_part(hessian_multivar_expression_reference const& a) { return a; }
	template<> inline typename ComplexType<hessian_variable>::type                                       real_part(hessian_variable const& a) { return a; }
	template<> inline typename ComplexType< std::complex<unknown> >::type                                real_part(std::complex<unknown> const& a) { return a.real(); }
	template<> inline typename ComplexType< std::complex<variable> >::type                               real_part(std::complex<variable> const& a) { return a.real(); }
	template<> inline typename ComplexType< std::complex<value_reference> >::type                        real_part(std::complex<value_reference> const& a) { return a.real(); }
	template<> inline typename ComplexType< std::complex<multivar_expression_reference> >::type          real_part(std::complex<multivar_expression_reference> const& a) { return a.real(); }
	template<> inline typename ComplexType< std::complex<hessian_multivar_expression_reference> >::type  real_part(std::complex<hessian_multivar_expression_reference> const& a) { return a.real(); }
	template<> inline typename ComplexType< std::complex<hessian_variable> >::type                       real_part(std::complex<hessian_variable> const& a) { return a.real(); }
	//template<> typename ComplexType<unknown>::type                                                imag_part(unknown const& a) { return 0.0; }
	//template<> typename ComplexType<variable>::type                                               imag_part(variable const& a) { return 0.0; }
	//template<> typename ComplexType<value_reference>::type                                        imag_part(value_reference const& a) { return 0.0; }
	//template<> typename ComplexType<multivar_expression_reference>::type	                        imag_part(multivar_expression_reference const& a) { return 0.0; }
	//template<> typename ComplexType<hessian_multivar_expression_reference>::type                  imag_part(hessian_multivar_expression_reference const& a) { return 0.0; }
	//template<> typename ComplexType<hessian_variable>::type                                       imag_part(hessian_variable const& a) { return 0.0; }
	//template<> typename ComplexType< std::complex<unknown> >::type                                imag_part(std::complex<unknown> const& a) { return a.imag(); }
	//template<> typename ComplexType< std::complex<variable> >::type                               imag_part(std::complex<variable> const& a) { return a.imag(); }
	//template<> typename ComplexType< std::complex<value_reference> >::type                        imag_part(std::complex<value_reference> const& a) { return a.imag(); }
	//template<> typename ComplexType< std::complex<multivar_expression_reference> >::type          imag_part(std::complex<multivar_expression_reference> const& a) { return a.imag(); }
	//template<> typename ComplexType< std::complex<hessian_multivar_expression_reference> >::type  imag_part(std::complex<hessian_multivar_expression_reference> const& a) { return a.imag(); }
	//template<> typename ComplexType< std::complex<hessian_variable> >::type                       imag_part(std::complex<hessian_variable> const& a) { return a.imag(); }

	//For INMOST_DATA_INTEGER_TYPE
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, unknown>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, variable>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, multivar_expression_reference>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<INMOST_DATA_INTEGER_TYPE, hessian_variable>  {typedef hessian_variable type;};
	//For INMOST_DATA_REAL_TYPE
	template<> struct Promote<INMOST_DATA_REAL_TYPE, unknown>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, variable>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct Promote<INMOST_DATA_REAL_TYPE, multivar_expression_reference>  {typedef variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<INMOST_DATA_REAL_TYPE, hessian_variable>  {typedef hessian_variable type;};
	//For INMOST_DATA_CPLX_TYPE
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, unknown> { typedef std::complex<variable> type; };
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, variable> { typedef std::complex<variable> type; };
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, value_reference> { typedef INMOST_DATA_CPLX_TYPE type; };
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, multivar_expression_reference> { typedef std::complex<variable> type; };
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, hessian_multivar_expression_reference> { typedef std::complex<hessian_variable> type; };
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, hessian_variable> { typedef std::complex<hessian_variable> type; };
	//for unknown
	template<> struct Promote<unknown, INMOST_DATA_INTEGER_TYPE>  {typedef variable type;};
	template<> struct Promote<unknown, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
	template<> struct Promote<unknown, unknown>  {typedef variable type;};
	template<> struct Promote<unknown, variable>  {typedef variable type;};
	template<> struct Promote<unknown, value_reference> { typedef variable type; };
	template<> struct Promote<unknown, multivar_expression_reference>  {typedef variable type;};
	template<> struct Promote<unknown, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<unknown, hessian_variable>  {typedef hessian_variable type;};
	//for value_reference
	template<> struct Promote<value_reference, INMOST_DATA_INTEGER_TYPE>  {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<value_reference, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct Promote<value_reference, unknown> { typedef variable type; };
	template<> struct Promote<value_reference, variable> { typedef variable type; };
	template<> struct Promote<value_reference, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct Promote<value_reference, multivar_expression_reference> { typedef variable type; };
	template<> struct Promote<value_reference, hessian_multivar_expression_reference> { typedef hessian_variable type; };
	template<> struct Promote<value_reference, hessian_variable> { typedef hessian_variable type; };
	//For multivar_expression_reference
	template<> struct Promote<multivar_expression_reference, INMOST_DATA_INTEGER_TYPE>  {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, unknown> {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, variable> {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, value_reference> { typedef variable type; };
	template<> struct Promote<multivar_expression_reference, multivar_expression_reference> {typedef variable type;};
	template<> struct Promote<multivar_expression_reference, hessian_multivar_expression_reference> {typedef hessian_variable type;};
	template<> struct Promote<multivar_expression_reference, hessian_variable> {typedef hessian_variable type;};
	//For variable
	template<> struct Promote<variable, INMOST_DATA_INTEGER_TYPE>  {typedef variable type;};
	template<> struct Promote<variable, INMOST_DATA_REAL_TYPE>  {typedef variable type;};
	template<> struct Promote<variable, unknown> {typedef variable type;};
	template<> struct Promote<variable, variable> {typedef variable type;};
	template<> struct Promote<variable, value_reference> { typedef variable type; };
	template<> struct Promote<variable, multivar_expression_reference> {typedef variable type;};
	template<> struct Promote<variable, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<variable, hessian_variable>  {typedef hessian_variable type;};
	//For hessian_multivar_expression_reference
	template<> struct Promote<hessian_multivar_expression_reference, INMOST_DATA_INTEGER_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, INMOST_DATA_REAL_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, unknown>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, variable>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, value_reference> { typedef hessian_variable type; };
	template<> struct Promote<hessian_multivar_expression_reference, multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, hessian_multivar_expression_reference> {typedef hessian_variable type;};
	template<> struct Promote<hessian_multivar_expression_reference, hessian_variable> {typedef hessian_variable type;};
	//For hessian_variable
	template<> struct Promote<hessian_variable, INMOST_DATA_INTEGER_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, INMOST_DATA_REAL_TYPE>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, unknown>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, variable>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, value_reference> { typedef hessian_variable type; };
	template<> struct Promote<hessian_variable, multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, hessian_multivar_expression_reference> {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, hessian_variable> {typedef hessian_variable type;};
#if defined(USE_FP64)
	template<> struct Promote<unknown, float>  {typedef variable type;};
	template<> struct Promote<float, unknown> { typedef variable type; };
	template<> struct Promote<value_reference, float>  {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<float, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct Promote<multivar_expression_reference, float> { typedef variable type; };
	template<> struct Promote<float, multivar_expression_reference> { typedef variable type; };
	template<> struct Promote<variable, float>  {typedef variable type;};
	template<> struct Promote<float, variable> { typedef variable type; };
	template<> struct Promote<hessian_multivar_expression_reference, float>  {typedef hessian_variable type;};
	template<> struct Promote<float, hessian_multivar_expression_reference> { typedef hessian_variable type; };
	template<> struct Promote<hessian_variable, float>  {typedef hessian_variable type;};
	template<> struct Promote<float, hessian_variable> { typedef hessian_variable type; };
#else
	template<> struct Promote<unknown, double>  {typedef variable type;};
	template<> struct Promote<double, unknown> { typedef variable type; };
	template<> struct Promote<value_reference, double>  {typedef INMOST_DATA_REAL_TYPE type;};
	template<> struct Promote<double, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct Promote<multivar_expression_reference, double> { typedef variable type; };
	template<> struct Promote<double, multivar_expression_reference> { typedef variable type; };
	template<> struct Promote<variable, double>  {typedef variable type;};
	template<> struct Promote<double, variable> { typedef variable type; };
	template<> struct Promote<hessian_multivar_expression_reference, double> { typedef hessian_variable type; };
	template<> struct Promote<double, hessian_multivar_expression_reference>  {typedef hessian_variable type;};
	template<> struct Promote<hessian_variable, double>  {typedef hessian_variable type;};
	template<> struct Promote<double, hessian_variable> { typedef hessian_variable type; };
#endif
#endif

	class AbstractMatrixBase
	{
	protected:
		static thread_private<Sparse::RowMerger> merger;
	};


	template<typename Var>
	class AbstractMatrixReadOnly : public AbstractMatrixBase
	{
	private:
		/// Sign function for SVD algorithm.
		static Var sign_func(const Var& a, const Var& b) { return (real_part(get_value(b)) >= 0.0 ? fabs(a) : -fabs(a)); }
		/// Max function for SVD algorithm.
		static INMOST_DATA_REAL_TYPE max_func(INMOST_DATA_REAL_TYPE x, INMOST_DATA_REAL_TYPE y) { return x > y ? x : y; }
		/// Function for QR rotation in SVD algorithm.
		static Var pythag(const Var& a, const Var& b)
		{
			Var at = fabs(a), bt = fabs(b), ct;
			if (real_part(get_value(at)) > real_part(get_value(bt))) { ct = bt / at; return at * sqrt(1.0 + ct * ct); }
			else if (real_part(get_value(bt)) > 0.0) { ct = at / bt; return bt * sqrt(1.0 + ct * ct); }
			else return 0.0;
			//return sqrt(a * conj(a) + b * conj(b));
		}
	public:
		typedef unsigned enumerator;
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		virtual bool TrivialArguments() const = 0;
		/// Obtain number of rows.
		/// @return Number of rows.
		virtual enumerator Rows() const = 0;
		/// Obtain number of columns.
		/// @return Number of columns.
		virtual enumerator Cols() const = 0;
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		virtual Var operator () (enumerator i, enumerator j) const = 0;
		/// Check all matrix entries for nans.
		/// Also checks derivatives for matrices of variables.
		bool CheckNans() const
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					if (check_nans((*this)(i, j))) return true;
			return false;
		}
		bool CheckInfs() const
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					if (check_infs((*this)(i, j))) return true;
			return false;
		}
		bool CheckNansInfs() const
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					if (check_nans_infs((*this)(i, j))) return true;
			return false;
		}
		/// Matrix determinant
		Var Det() const
		{
			assert(Rows() == Cols());
			const double eps = 1.0e-13;
			const enumerator n = Rows();
			Matrix<Var> A(*this);
			for (enumerator d = 0; d < n; ++d)
				for (enumerator i = d + 1; i < n; ++i)
				{
					if (fabs(get_value(A(d, d))) < eps)
						A(d, d) = eps;
					for (enumerator j = 0; j < n; ++j)
						A(i, j) = A(i, j) - A(i, d) * A(d, j) / A(d, d);
				}
			Var ret = 1.0;
			for (enumerator d = 0; d < n; ++d)
				ret *= A(d, d);
			return ret;
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
		void MPT(INMOST_DATA_ENUM_TYPE* Perm, INMOST_DATA_REAL_TYPE* SL = NULL, INMOST_DATA_REAL_TYPE* SR = NULL) const;
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
		bool SVD(AbstractMatrix<Var>& U, AbstractMatrix<Var>& Sigma, AbstractMatrix<Var>& V, bool order_singular_values = true, bool nonnegative = true) const;

		bool cSVD(AbstractMatrix<Var>& U, AbstractMatrix<Var>& Sigma, AbstractMatrix<Var>& V) const;

        /// Eigen decomposition.
        /// Reconstruct matrix: A = V*Lambda*V.Transpose().
        /// @param V Eigenvectors matrix.
        /// @param Lambda Diagonal matrix with eigen values.
        /// @param order_eigen_values Return eigen values in ascending order of their real part.
        bool EigenDecomposition(AbstractMatrix<INMOST_DATA_CPLX_TYPE> & V, AbstractMatrix<INMOST_DATA_CPLX_TYPE> & Lambda, bool order_eigen_values = true,
                double abs_eps=1e-8, double rel_eps=1e-8, int max_iters=256) const;

        /// Householder reflector.
        /// Computes householder reflector for a vector.
        /// @param first_n_untouched Number of components left intact by the reflector.
        template<typename RealType = Var>
		typename std::enable_if<!is_complex<RealType>::value, Matrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type>>::type HouseholderReflector(enumerator first_n_untouched = 0) const
        {
            size_t n_rows = Rows();
            assert(Cols() == 1);
            assert(first_n_untouched < n_rows);

            // Проекция вектора на первые first_n_untouched осей.
            Matrix<Var> projected_vector(n_rows,1);
            projected_vector.Zero();
            for (enumerator index = first_n_untouched; index < n_rows; ++index)
            {
                projected_vector(index, 0) = (*this)(index, 0);
            }

            // Единичный вектор first_n_untouched-ой оси.
            Matrix<Var> axis(n_rows,1);
            axis(first_n_untouched, 0) = sign_func(1.0, (*this)(first_n_untouched, 0)) * projected_vector.FrobeniusNorm();

            // Вектор рефлектора.
            Matrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type> reflector = projected_vector - axis;
            auto reflector_norm = reflector.FrobeniusNorm();
            if (reflector_norm > 0.0)
            { reflector /= reflector_norm; }

            return reflector;
        }
        template<typename ComplexType = Var>
		typename std::enable_if<is_complex<ComplexType>::value, Matrix<INMOST_DATA_CPLX_TYPE>>::type HouseholderReflector(enumerator first_n_untouched = 0) const
        {
            size_t n_rows = Rows();
            assert(Cols() == 1);
            assert(first_n_untouched < n_rows);

            // Проекция вектора на первые first_n_untouched осей.
            Matrix<Var> projected_vector(n_rows,1);
            projected_vector.Zero();
            for (enumerator index = first_n_untouched; index < n_rows; ++index)
            {
                projected_vector(index, 0) = (*this)(index, 0);
            }

            // Фазовый множитель.
            auto phase_mult = (*this)(first_n_untouched, 0);
            if (std::norm(phase_mult) > 0.0)
            { phase_mult /= std::abs(phase_mult); }

            // Единичный вектор first_n_untouched-ой оси.
            Matrix<Var> axis(n_rows,1);
            axis(first_n_untouched, 0) = phase_mult * projected_vector.FrobeniusNorm();

            // Вектор рефлектора.
            Matrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type> reflector = projected_vector - axis;
            auto reflector_norm = reflector.FrobeniusNorm();
            if (reflector_norm > 0.0)
            { reflector /= reflector_norm; }

            return reflector;
        }


        /// Hessenberg form.
        /// Reconstruct matrix: A = Q*H*Q.T
        /// @param Q Unitary matrix.
        /// @param H Hessenberg matrix.
        bool Hessenberg(AbstractMatrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type> & Q,
                        AbstractMatrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type> & H) const
        {
            enumerator dim = Rows();
            if (dim != Cols()) return false;

            // Унитарная матрица.
            Q = Matrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type>::Unit(dim);

            // Верхнехессенбергова матрица.
            H = (*this);

            for (enumerator step = 0; step + 2 < dim; ++step)
            {
                // Получение отражающего преобразования.
                Matrix<Var> reflector = H.ExtractSubMatrix(0, dim, step, step + 1).HouseholderReflector(step + 1);
                Matrix<Var> outer_product = reflector * reflector.ConjugateTranspose();

                // Применение отражающего преобразования.
                H = H - 2.0 * H * outer_product;
                H = H - 2.0 * outer_product * H;
                Q = Q * (Matrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type>::Unit(dim) - 2.0 * outer_product);
            }

            return true;
        }


        /// QR-decomposition.
        /// Reconstruct matrix: A = Q*R.
        /// @param Q Unitary matrix.
        /// @param R Upper triangular matrix.
        bool QRDecompoition(AbstractMatrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type> & Q,
                            AbstractMatrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type> & R) const
        {
            enumerator n_rows = Rows();

            // Унитарная матрица.
            Q = Matrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type>::Unit(n_rows);

            // Верхнетреугольная матрица.
            R = (*this);

            for (enumerator step = 0; step < n_rows; ++step)
            {
                // Получение отражающего преобразования.
                Matrix<Var> reflector = R.ExtractSubMatrix(0, n_rows, step, step + 1).HouseholderReflector(step);
                Matrix<Var> outer_product = reflector * reflector.ConjugateTranspose();

                // Применение отражающего преобразования.
                R = R - 2.0 * outer_product * R;
                Q = Q * (Matrix<typename Promote<Var, INMOST_DATA_REAL_TYPE>::type>::Unit(n_rows) - 2.0 * outer_product);
            }

            return true;
        }



		/// Transpose current matrix.
		/// @return Transposed matrix.
		//Matrix<Var> Transpose() const;
		ConstMatrixTranspose<Var> Transpose() const { return ConstMatrixTranspose<Var>(*this); }
		/// Transpose and conjugate current matrix.
		/// @return Transposed and conjugated matrix.
		ConstMatrixConjugateTranspose<Var> ConjugateTranspose() const { return ConstMatrixConjugateTranspose<Var>(*this); }
		/// Conjugate  current matrix.
		/// @return Conjugated matrix.
		ConstMatrixConjugate<Var> Conjugate() const { return ConstMatrixConjugate<Var>(*this); }
		/// Cross-product operation for a vector.
		/// Both right hand side and left hand side should be a vector
		/// @param other The right hand side of cross product.
		/// @return The cross product of current and right hand side vector.
		/// \warning Works for 3x1 vector and 3xm m-vectors as right hand side.
		template<typename typeB>
		Matrix<typename Promote<Var, typeB>::type>
			CrossProduct(const AbstractMatrixReadOnly<typeB>& other) const;
		/// Transformation matrix from current vector to provided vector using shortest arc rotation.
		/// @param other Vector to transform to.
		/// @return A sqaure (rotation) matrix that transforms current vector into right hand side vector.
		/// \warning Works only for 3x1 vectors.
		template<typename typeB>
		Matrix<typename Promote<Var, typeB>::type>
			Transform(const AbstractMatrixReadOnly<typeB>& other) const;
		/// Subtract a matrix.
		/// @param other The matrix to be subtracted.
		/// @return Difference of current matrix with another matrix.
		template<typename typeB>
		//Matrix<typename Promote<Var, typeB>::type>
		MatrixDifference<Var,typeB>
			operator-(const AbstractMatrixReadOnly<typeB>& other) const;
		/// Add a matrix.
		/// @param other The matrix to be added.
		/// @return Sum of current matrix with another matrix.
		template<typename typeB>
		//Matrix<typename Promote<Var, typeB>::type>
		MatrixSum<Var,typeB>
			operator+(const AbstractMatrixReadOnly<typeB>& other) const;
		/// Multiply the matrix by another matrix.
		/// @param other Matrix to be multiplied from right.
		/// @return Matrix multipled by another matrix.
		template<typename typeB>
		//Matrix<typename Promote<Var, typeB>::type>
		MatrixMul<Var,typeB, typename Promote<Var,typeB>::type>
			operator*(const AbstractMatrixReadOnly<typeB>& other) const;
		/// Unary minus. Change sign of each element of the matrix.
		MatrixUnaryMinus<Var> operator-() const { return MatrixUnaryMinus<Var>(*this); }
		/// Kronecker product, latex symbol \otimes.
		/// @param other Matrix on the right of the Kronecker product.
		/// @return Result of the Kronecker product of current and another matrix.
		template<typename typeB>
		//Matrix<typename Promote<Var, typeB>::type>
		KroneckerProduct<Var, typeB>
			Kronecker(const AbstractMatrixReadOnly<typeB>& other) const;
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
		Matrix<Var> Invert(int* ierr = NULL) const;
		/// Inverts symmetric positive-definite matrix using Cholesky decomposition.
		Matrix<Var> CholeskyInvert(int* ierr = NULL) const;
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
		Matrix<typename Promote<Var, typeB>::type>
			Solve(const AbstractMatrixReadOnly<typeB>& B, int* ierr = NULL) const;
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
		Matrix<typename Promote<Var, typeB>::type>
			CholeskySolve(const AbstractMatrixReadOnly<typeB>& B, int* ierr = NULL) const;
		/// Calculate sum of the diagonal elements of the matrix.
		/// @return Trace of the matrix.
		Var Trace() const
		{
			assert(Cols() == Rows());
			Var ret = 0.0;
			for (enumerator i = 0; i < Rows(); ++i) ret += (*this)(i, i);
			return ret;
		}
		/// Output matrix to screen.
		/// Does not output derivatices.
		/// @param threshold Elements smaller then the number are considered to be zero.
		void Print(INMOST_DATA_REAL_TYPE threshold = 1.0e-10, std::ostream & sout = std::cout) const
		{
			for (enumerator k = 0; k < Rows(); ++k)
			{
				for (enumerator l = 0; l < Cols(); ++l)
				{
					//if (__isinf__(real_part(get_value((*this)(k, l)))))
					//	sout << std::setw(12) << "inf";
					//else if (std::isnan(get_value((*this)(k, l))))
					//	sout << std::setw(12) << "nan";
					//else
					if (fabs(get_value((*this)(k, l))) > threshold)
						sout << std::setw(12) << get_value((*this)(k, l));
					else
						sout << std::setw(12) << 0;
					sout << " ";
				}
				sout << std::endl;
			}
		}
		/// Check if the matrix is symmetric.
		/// @return Returns true if the matrix is symmetric, otherwise false.
		bool isSymmetric(double eps = 1.0e-7) const
		{
			if (Rows() != Cols()) return false;
			for (enumerator k = 0; k < Rows(); ++k)
			{
				for (enumerator l = k + 1; l < Rows(); ++l)
					if (fabs((*this)(k, l) - (*this)(l, k)) > eps)
						return false;
			}
			return true;
		}
		/// Computes dot product by summing up multiplication of entries with the
		/// same indices in the current and the provided matrix.
		/// @param other Provided matrix.
		/// @return Dot product of two matrices.
		template<typename typeB>
		typename Promote<Var, typeB>::type
			DotProduct(const AbstractMatrixReadOnly<typeB>& other) const
		{
			assert(Cols() == other.Cols());
			assert(Rows() == other.Rows());
			typename Promote<Var, typeB>::type ret = 0.0;
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					ret += ((*this)(i, j)) * other(i, j);
			return ret;
		}
		/// Computes dot product by summing up multiplication of entries with the
		/// same indices in the current and the provided matrix.
		/// @param other Provided matrix.
		/// @return Dot product of two matrices.
		template<typename typeB>
		typename Promote<Var, typeB>::type operator ^(const AbstractMatrixReadOnly<typeB>& other) const
		{
			return DotProduct(other);
		}
		/// Computes frobenious norm of the matrix.
		/// @return Frobenius norm of the matrix.
        template<typename RealType = Var>
		typename std::enable_if<!is_complex<RealType>::value, typename SelfPromote<RealType>::type>::type FrobeniusNorm() const
		{
            typename SelfPromote<RealType>::type ret = 0;
            for (enumerator i = 0; i < Rows(); ++i)
                for (enumerator j = 0; j < Cols(); ++j)
                    ret += (*this)(i, j) * (*this)(i, j);
            return sqrt(ret);
		}
        template<typename ComplexType = Var>
		typename std::enable_if<is_complex<ComplexType>::value, typename SelfPromote<typename ComplexType::value_type>::type>::type FrobeniusNorm() const
		{
            typename SelfPromote<typename ComplexType::value_type>::type ret = 0;
            for (enumerator i = 0; i < Rows(); ++i)
                for (enumerator j = 0; j < Cols(); ++j)
                    ret += std::real((*this)(i, j) * std::conj((*this)(i, j)));
            return sqrt(fabs(ret));
		}



		/// Computes maximum absolute value of the matrix.
		/// @return Maximum norm of the matrix.
		Var MaxNorm() const
		{
			Var ret = 0;
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
				{
					if (ret < fabs((*this)(i, j)))
						ret = fabs((*this)(i, j));
				}
			return ret;
		}
		/// Calculates Moore-Penrose pseudo-inverse of the matrix.
		/// @param tol Thershold for singular values. Singular values smaller
		///            then threshold are considered to be zero.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr = 1, in case of no failure *ierr = 0.
		/// @return A pseudo-inverse of the matrix.
		Matrix<Var> PseudoInvert(INMOST_DATA_REAL_TYPE tol = 0, int* ierr = NULL) const;
		/// Calcuate A^n, where n is some real value.
		/// \warning Currently works for symmetric matrices only!
		/// @param n Real value.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr = 1, in case of no failure *ierr = 0.
		///             The error may be caused by error in SVD algorithm.
		/// @return Matrix in power of n.
		Matrix<Var> Power(INMOST_DATA_REAL_TYPE n, int * ierr = NULL) const;
		/// Calculate square root of A matrix by Babylonian method.
		/// @param iter Number of iterations.
		/// @param tol  Convergence tolerance.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr = 1, in case of no failure *ierr = 0.
		/// @return Square root of a matrix.
		Matrix<Var> Root(INMOST_DATA_ENUM_TYPE iter = 25, INMOST_DATA_REAL_TYPE tol = 1.0e-7, int* ierr = NULL) const;
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
		Matrix<typename Promote<Var, typeB>::type>
			PseudoSolve(const AbstractMatrixReadOnly<typeB>& B, INMOST_DATA_REAL_TYPE tol = 0, int* ierr = NULL) const;
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
		ConstMatrixRepack<Var> Repack(enumerator rows, enumerator cols) const {return ConstMatrixRepack<Var>(*this, rows, cols);}
		/// Extract submatrix of a matrix for in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param first_row Starting row in the original matrix.
		/// @param last_row Last row (excluded) in the original matrix.
		/// @param first_col Starting column in the original matrix.
		/// @param last_col Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		ConstSubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const;
		/// Define matrix as a part of a matrix of larger size with in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, if i in [offset_row,offset_row+n),
		/// and j in [offset_col,offset_col+m) and B = {0} otherwise.
		/// @param nrows Number of rows in larger matrix.
		/// @param ncols Number of columns in larger matrix.
		/// @param offset_row Offset for row number.
		/// @param offset_col Offset for column number.
		/// @return Submatrix of the original matrix.
		ConstBlockOfMatrix<Var> BlockOf(enumerator nrows, enumerator ncols, enumerator offset_row, enumerator offset_col) const;
		/// Multiply the matrix by a coefficient.
		/// @param coef Coefficient.
		/// @return Matrix multiplied by the coefficient.
		template<typename typeB>
		//Matrix<typename Promote<Var, typeB>::type>
		MatrixMulCoef<Var, typeB, typename Promote<Var, typeB>::type>
			operator*(const typeB& coef) const;
#if defined(USE_AUTODIFF)
		/// Multiply the matrix by a coefficient.
		/// @param coef Coefficient.
		/// @return Matrix multiplied by the coefficient.
		template<class A>
		//Matrix<typename Promote<Var, variable>::type>
		MatrixMulShellCoef<Var, shell_expression<A>, typename Promote<Var, variable>::type>
			operator*(shell_expression<A> const& coef) const;// { return operator*(variable(coef)); }
#endif //USE_AUTODIFF
		/// Divide the matrix by a coefficient of a different type.
		/// @param coef Coefficient.
		/// @return Matrix divided by the coefficient.
		template<typename typeB>
		//Matrix<typename Promote<Var, typeB>::type>
		MatrixDivCoef<Var, typeB, typename Promote<Var, typeB>::type>
			operator/(const typeB & coef) const;
#if defined(USE_AUTODIFF)
		/// Divide the matrix by a coefficient of a different type.
		/// @param coef Coefficient.
		/// @return Matrix divided by the coefficient.
		template<class A>
		//Matrix<typename Promote<Var, variable>::type>
		MatrixDivShellCoef<Var, shell_expression<A>, typename Promote<Var, variable>::type>
			operator/(shell_expression<A> const& coef) const;// { return operator/(variable(coef)); }
#endif //USE_AUTODIFF
		/// Performs B^{-1}*A, multiplication by inverse matrix from left.
		/// Throws exception if matrix is not invertable. See Mesh::PseudoSolve for
		/// singular matrices.
		/// @param other Matrix to be inverted and multiplied from left.
		/// @return Multiplication of current matrix by inverse of another
		/// @see Matrix::PseudoInvert.
		template<typename typeB>
		Matrix<typename Promote<Var, typeB>::type>
			operator/(const AbstractMatrixReadOnly<typeB>& other) const
		{
			Matrix<typename Promote<Var, typeB>::type> ret(other.Cols(), Cols());
			ret = other.Solve(*this);
			return ret;
		}
		/// Concatenate B matrix as columns of current matrix.
		/// Assumes that number of rows of current matrix is
		/// equal to number of rows of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatRows
		ConstMatrixConcatCols<Var> ConcatCols(const AbstractMatrixReadOnly<Var>& B) const
		{
			return ConstMatrixConcatCols<Var>(*this, B);
		}
		/// Concatenate B matrix as columns of current matrix.
		/// Assumes that number of rows of current matrix is
		/// equal to number of rows of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatRows
		template<typename VarB>
		ConstMatrixConcatCols2<Var, VarB, typename Promote<Var, VarB>::type>
			ConcatCols(const AbstractMatrixReadOnly<VarB>& B) const
		{
			return ConstMatrixConcatCols2<Var, VarB, typename Promote<Var, VarB>::type>(*this, B);
		}
		/*
		{

			assert(Rows() == B.Rows());
			Matrix<Var> ret(Rows(),Cols()+B.Cols());
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
		*/
		/// Concatenate B matrix as rows of current matrix.
		/// Assumes that number of colums of current matrix is
		/// equal to number of columns of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatCols
		ConstMatrixConcatRows<Var> ConcatRows(const AbstractMatrixReadOnly<Var>& B) const
		{
			return ConstMatrixConcatRows<Var>(*this, B);
		}
		/// Concatenate B matrix as rows of current matrix.
		/// Assumes that number of colums of current matrix is
		/// equal to number of columns of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatCols
		template<typename VarB>
		ConstMatrixConcatRows2<Var,VarB,typename Promote<Var,VarB>::type>
			ConcatRows(const AbstractMatrixReadOnly<VarB>& B) const
		{
			return ConstMatrixConcatRows2<Var, VarB, typename Promote<Var, VarB>::type>(*this, B);
		}
		/*
		{
			assert(Cols() == B.Cols());
			Matrix<Var> ret(Rows()+B.Rows(),Cols());
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
		*/
		/// Destructor
		virtual ~AbstractMatrixReadOnly() {}
		/// Retrive interval of indices of derivatives.
		__INLINE virtual void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			for (INMOST_DATA_ENUM_TYPE i = 0; i < Rows(); ++i)
				for (INMOST_DATA_ENUM_TYPE j = 0; j < Cols(); ++j)
					GetInterval((*this)(i, j), beg, end, cnt);
		}

	};
	/// Abstract class for a matrix,
	/// used to abstract away all the data storage and access
	/// and provide common implementation of the algorithms.
	template<typename Var>
	class AbstractMatrix : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		using AbstractMatrixReadOnly<Var>::Rows;
		using AbstractMatrixReadOnly<Var>::Cols;
		using AbstractMatrixReadOnly<Var>::Transpose;
		using AbstractMatrixReadOnly<Var>::BlockOf;
		using AbstractMatrixReadOnly<Var>::Repack;
		using AbstractMatrixReadOnly<Var>::ConcatRows;
		using AbstractMatrixReadOnly<Var>::ConcatCols;

		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator;
		/// Construct empty matrix.
		AbstractMatrix() {}
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; }
		/*
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
		*/
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
		AbstractMatrix & operator =(AbstractMatrixReadOnly<typeB> const & other)
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
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		virtual Var & operator () (enumerator i, enumerator j) = 0;
		/// Resize the matrix into different size.
		/// @param nrows New number of rows.
		/// @param ncols New number of columns.
		virtual void Resize(enumerator rows, enumerator cols) = 0;
		/// Set all the elements of the matrix to zero.
		void Zero()
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					(*this)(i,j) = 0.0;
		}
		///Exchange contents of two matrices.
		virtual void Swap(AbstractMatrix<Var> & b);
		/// Subtract a matrix and store result in the current.
		/// @param other The matrix to be subtracted.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator-=(const AbstractMatrixReadOnly<typeB> & other);
		/// Add a matrix and store result in the current.
		/// @param other The matrix to be added.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator+=(const AbstractMatrixReadOnly<typeB> & other);
		/// Multiply matrix with another matrix in-place.
		/// @param B Another matrix to the right in multiplication.
		/// @return Reference to current matrix.
		template<typename typeB>
		AbstractMatrix & operator*=(const AbstractMatrixReadOnly<typeB> & B);
		/// Multiply the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator*=(typeB coef);
		/// Divide the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator/=(typeB coef);
		/// Extract submatrix of a matrix for in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param first_row Starting row in the original matrix.
		/// @param last_row Last row (excluded) in the original matrix.
		/// @param first_col Starting column in the original matrix.
		/// @param last_col Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		SubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);
		/// Define matrix as a part of a matrix of larger size with in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, if i in [offset_row,offset_row+n),
		/// and j in [offset_col,offset_col+m) and B = {0} otherwise.
		/// @param nrows Number of rows in larger matrix.
		/// @param ncols Number of columns in larger matrix.
		/// @param offset_row Offset for row number.
		/// @param offset_col Offset for column number.
		/// @return Submatrix of the original matrix.
		BlockOfMatrix<Var> BlockOf(enumerator nrows, enumerator ncols, enumerator offset_row, enumerator offset_col);
		/// Transpose the current matrix with access to elements.
		/// @return Transposed matrix.
		//Matrix<Var> Transpose() const;
		MatrixTranspose<Var> Transpose() { return MatrixTranspose<Var>(*this); }
		/// Change representation of the matrix into matrix of another size.
		/// Useful to change representation from matrix into vector and back.
		/// Replaces original number of columns and rows with a new one.
		/// @return Matrix with same entries and provided number of rows and columns.
		MatrixRepack<Var> Repack(enumerator rows, enumerator cols) { return MatrixRepack<Var>(*this, rows, cols); }
		/// Concatenate B matrix as columns of current matrix.
		/// Assumes that number of rows of current matrix is
		/// equal to number of rows of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatRows
		MatrixConcatCols<Var> ConcatCols(AbstractMatrix<Var>& B) { return MatrixConcatCols<Var>(*this, B); }
		/// Concatenate B matrix as rows of current matrix.
		/// Assumes that number of colums of current matrix is
		/// equal to number of columns of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatCols
		MatrixConcatRows<Var> ConcatRows(AbstractMatrix<Var>& B) { return MatrixConcatRows<Var>(*this, B); }
		/// Destructor
		virtual ~AbstractMatrix() {};
		/// Retrive interval of indices of derivatives.
		//virtual void GetMatrixInterval(INMOST_DATA_ENUM_TYPE & beg, INMOST_DATA_ENUM_TYPE & end, INMOST_DATA_ENUM_TYPE & cnt) const
		//{
		//	for (INMOST_DATA_ENUM_TYPE i = 0; i < Rows(); ++i)
		//		for (INMOST_DATA_ENUM_TYPE j = 0; j < Cols(); ++j)
		//			GetInterval((*this)(i, j), beg, end, cnt);
		//}

	};


	template<typename Var, typename storage_type>
	class SymmetricMatrix : public AbstractMatrix<Var>
	{
	public:
		using AbstractMatrix<Var>::operator();
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
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
		/// Constract a matrix with provided elements.
		/// The elements are ordered row-wise starting from the diagonal element.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		static SymmetricMatrix Make(enumerator pn, ...)
		{
			SymmetricMatrix A(pn);
			va_list argptr;
			va_start(argptr, pn);
			for (enumerator j = 0; j < pn; ++j)
				for (enumerator i = j; i < pn; ++i)
				{
					Var val = va_arg(argptr, Var);
					A(i, j) = val;
				}
			va_end(argptr);
			return A;
		}
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
		SymmetricMatrix(const AbstractMatrixReadOnly<typeB> & other) : space(other.Rows()*(other.Rows()+1)/2), n(other.Rows())
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
			(void)ncols;
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
		__INLINE Var & operator()(enumerator i, enumerator j)
		{
			assert(i < n);
			assert(j < n);
			if( i > j ) std::swap(i,j);
			return space[j+n*i-i*(i+1)/2];
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			assert(i < n);
			assert(j < n);
			if( i > j ) std::swap(i,j);
			return space[j+n*i-i*(i+1)/2];
		}

		/// Return raw pointer to matrix data, stored in row-wise format.
		/// @return Pointer to data.
		__INLINE Var * data() {return space.data();}
		/// Return raw pointer to matrix data without right of change,
		/// stored in row-wise format.
		/// @return Pointer to constant data.
		__INLINE const Var * data() const {return space.data();}
		/// Obtain number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const {return n;}
		/// Obtain number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const {return n;}
		/// Obtain number of rows.
		/// @return Reference to number of rows.
		__INLINE enumerator & Rows() {return n;}
		/// Obtain number of rows.
		/// @return Reference to number of columns.
		__INLINE enumerator & Cols() {return n;}
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
		static SymmetricMatrix<Var> FromTensor(const Var * K, enumerator size, enumerator matsize = 3)
		{
			SymmetricMatrix<Var> Kc(matsize,matsize);
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
		static SymmetricMatrix<Var> FromDiagonal(const Var * r, enumerator size)
		{
			SymmetricMatrix<Var> ret(size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
			return ret;
		}
		/// Create diagonal matrix from array of values that have to be inversed.
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by inverse of array elements.
		static SymmetricMatrix<Var> FromDiagonalInverse(const Var * r, enumerator size)
		{
			SymmetricMatrix<Var> ret(size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = 1.0/r[k];
			return ret;
		}
		/// Unit matrix. Creates a square matrix of size pn by pn
		/// and fills the diagonal with c.
		/// @param pn Number of rows and columns in the matrix.
		/// @param c Value to put onto diagonal.
		/// @return Returns a unit matrix.
		static SymmetricMatrix<Var> Unit(enumerator pn, const Var & c = 1.0)
		{
			SymmetricMatrix<Var> ret(pn,0.0);
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
		//::INMOST::SubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);

		//::INMOST::ConstSubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const;
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
		using AbstractMatrix<Var>::operator();
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
	protected:
		storage_type space; //< Array of row-wise stored elements.
		enumerator n; //< Number of rows.
		enumerator m; //< Number of columns.
	public:
		Var& operator [](enumerator k) { return space[k]; }
		const Var& operator[](enumerator k) const { return space[k]; }
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
		/// Constract a matrix with provided elements.
		/// The elements are ordered row-wise.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		static Matrix Make(enumerator pn, enumerator pm, ...)
		{
			Matrix A(pn, pm);
			va_list argptr;
			va_start(argptr, pm);
			for (enumerator i = 0; i < pn; ++i)
				for (enumerator j = 0; j < pm; ++j)
				{
					Var val = va_arg(argptr, Var);
					A(i, j) = val;
				}
			va_end(argptr);
			return A;
		}
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
		Matrix(const AbstractMatrixReadOnly<typeB> & other) : space(other.Cols()*other.Rows()), n(other.Rows()), m(other.Cols())
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
		Matrix & operator =(AbstractMatrixReadOnly<typeB> const & other)
		{
			if( Cols()*Rows() != other.Cols()*other.Rows() )
				space.resize(other.Cols()*other.Rows());
			n = other.Rows();
			m = other.Cols();
			for(enumerator i = 0; i < other.Rows(); ++i)
				for(enumerator j = 0; j < other.Cols(); ++j)
					assign((*this)(i,j),other(i,j));
			return *this;
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		__INLINE Var & operator()(enumerator i, enumerator j)
		{
			assert(i < n);
			assert(j < m);
			assert(i*m+j < n*m); //overflow check?
			return space[i*m+j];
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			assert(i < n);
			assert(j < m);
			assert(i*m+j < n*m); //overflow check?
			return space[i*m+j];
		}

		/// Return raw pointer to matrix data, stored in row-wise format.
		/// @return Pointer to data.
		__INLINE Var * data() {return space.data();}
		/// Return raw pointer to matrix data without right of change,
		/// stored in row-wise format.
		/// @return Pointer to constant data.
		__INLINE const Var * data() const {return space.data();}
		/// Obtain number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const {return n;}
		/// Obtain number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const {return m;}
		/// Obtain number of rows.
		/// @return Reference to number of rows.
		__INLINE enumerator & Rows() {return n;}
		/// Obtain number of rows.
		/// @return Reference to number of columns.
		__INLINE enumerator & Cols() {return m;}
		/// Construct row permutation matrix from array of new positions for rows.
		/// Row permutation matrix multiplies matrix from left.
		/// Column permutation matrix is obtained by transposition and is multiplied
		/// from the right.
		/// Argument Perm is filled with the new position for rows, i.e.
		/// i-th row takes new position Perm[i]
		/// @param Perm Array with new positions for rows.
		/// @param size Size of the array and the resulting matrix.
		/// @return Permutation matrix.
		static Matrix<Var> Permutation(const INMOST_DATA_ENUM_TYPE * Perm, enumerator size)
		{
			Matrix<Var> Ret(size,size);
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
		static Matrix<Var> FromTensor(const Var * K, enumerator size, enumerator matsize = 3)
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
		static Matrix<Var> FromVector(const Var * r, enumerator size)
		{
			return Matrix<Var>(r,size,1);
		}
		/// Create diagonal matrix from array
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by array, other elements are zero.
		static Matrix<Var> FromDiagonal(const Var * r, enumerator size)
		{
			Matrix<Var> ret(size,size);
			ret.Zero();
			for(enumerator k = 0; k < size; ++k) ret(k,k) = r[k];
			return ret;
		}
		/// Create diagonal matrix from array of values that have to be inversed.
		/// @param r Array of diagonal elements.
		/// @param size Size of the matrix.
		/// @return Matrix with diagonal defined by inverse of array elements.
		static Matrix<Var> FromDiagonalInverse(const Var * r, enumerator size)
		{
			Matrix<Var> ret(size,size);
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
		static Matrix<Var> CrossProductMatrix(const Var vec[3])
		{
			// |  0  -z   y |
			// |  z   0  -x |
			// | -y   x   0 |
			Matrix<Var> ret(3,3);
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
		static MatrixUnit<Var> Unit(enumerator pn, const Var & c = 1.0)
		{
			return MatrixUnit<Var>(pn, c);
			/*
			Matrix<Var> ret(pn,pn);
			ret.Zero();
			for(enumerator i = 0; i < pn; ++i) ret(i,i) = c;
			return ret;
			*/
		}
		/// Matix with 1 row, Create a matrix of size 1 by pn and
		/// fills it with c.
		/// @param pn Number of columns.
		/// @param c Value to fill the matrix.
		/// @return Returns a matrix with 1 row.
		static MatrixRow<Var> Row(enumerator pn, const Var & c = 1.0)
		{ return MatrixRow<Var>(pn,c); }
		/// Matix with 1 column, Create a matrix of size pn by 1 and
		/// fills it with c.
		/// @param pn Number of rows.
		/// @param c Value to fill the matrix.
		/// @return Returns a matrix with 1 column.
		static MatrixCol<Var> Col(enumerator pn, const Var & c = 1.0)
		{ return MatrixCol<Var>(pn, c); }
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
		Matrix<Var> JointDiagonalization(INMOST_DATA_REAL_TYPE threshold = 1.0e-7)
		{
			enumerator N = Rows();
			enumerator M = Cols() / Rows();
			Matrix<Var> V((MatrixUnit<Var>(m)));
			Matrix<Var> R(2,M);
			Matrix<Var> G(2,2);
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
		//::INMOST::SubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col);

		//::INMOST::ConstSubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const;
	};
	/// This class allows for in-place operations on submatrix of the matrix elements.
	template<typename Var>
	class SubMatrix : public AbstractMatrix<Var>
	{
	public:
		using AbstractMatrix<Var>::operator();
		using AbstractMatrix<Var>::operator =;
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		AbstractMatrix<Var> * M;
		enumerator brow; //< First row in matrix M.
		enumerator erow; //< Last row in matrix M.
		enumerator bcol; //< First column in matrix M.
		enumerator ecol; //< Last column in matrix M.
	public:
		/// Number of rows in submatrix.
		/// @return Number of rows.
		__INLINE enumerator Rows() const {return erow-brow;}
		/// Number of columns in submatrix.
		/// @return Number of columns.
		__INLINE enumerator Cols() const {return ecol-bcol;}
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
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element.
		__INLINE Var & operator()(enumerator i, enumerator j)
		{
			assert(i < Rows());
			assert(j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return (*M)(i+brow,j+bcol);
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			assert(i < Rows());
			assert(j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return (*M)(i+brow,j+bcol);
		}
		/// Convert submatrix into matrix.
		/// Note, that modifying returned matrix does
		/// not affect elements of the submatrix or original matrix
		/// used to create submatrix.
		/// @return Matrix with same entries as submatrix.
		::INMOST::Matrix<Var> MakeMatrix()
		{
			::INMOST::Matrix<Var> ret(Rows(),Cols());
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
	class ConstSubMatrix : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var> * M;
		enumerator brow; //< First row in matrix M.
		enumerator erow; //< Last row in matrix M.
		enumerator bcol; //< First column in matrix M.
		enumerator ecol; //< Last column in matrix M.
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const {return M->TrivialArguments();};
		/// Number of rows in submatrix.
		/// @return Number of rows.
		__INLINE enumerator Rows() const {return erow-brow;}
		/// Number of columns in submatrix.
		/// @return Number of columns.
		__INLINE enumerator Cols() const {return ecol-bcol;}
		/// Create submatrix for a matrix.
		/// @param rM Reference to the matrix that stores elements.
		/// @param first_row First row in the matrix.
		/// @param last_row Last row in the matrix.
		/// @param first_column First column in the matrix.
		/// @param last_column Last column in the matrix.
		ConstSubMatrix(const AbstractMatrixReadOnly<Var> & rM, enumerator first_row, enumerator last_row, enumerator first_column, enumerator last_column) : M(&rM), brow(first_row), erow(last_row), bcol(first_column), ecol(last_column)
		{}
		ConstSubMatrix(const ConstSubMatrix & b) : M(b.M), brow(b.brow), erow(b.erow), bcol(b.bcol), ecol(b.ecol) {}

		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			assert(i < Rows());
			assert(j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return (*M)(i+brow,j+bcol);
		}
		/// Convert submatrix into matrix.
		/// Note, that modifying returned matrix does
		/// not affect elements of the submatrix or original matrix
		/// used to create submatrix.
		/// @return Matrix with same entries as submatrix.
		::INMOST::Matrix<Var> MakeMatrix()
		{
			::INMOST::Matrix<Var> ret(Rows(),Cols());
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i,j) = (*this)(i,j);
			return ret;
		}
	};

	/// This class allows to address a matrix as a block of an empty matrix of larger size.
	template<typename Var>
	class BlockOfMatrix : public AbstractMatrix<Var>
	{
	public:
		using AbstractMatrix<Var>::operator();
		using AbstractMatrix<Var>::operator =;
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		AbstractMatrix<Var> * M;
		enumerator nrows; //< Number of rows in larger matrix.
		enumerator ncols; //< Number of colums in larger matrix.
		enumerator orow; //< Row offset in larger matrix.
		enumerator ocol; //< Column offset in larger matrix.
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows in submatrix.
		/// @return Number of rows.
		__INLINE enumerator Rows() const {return nrows;}
		/// Number of columns in submatrix.
		/// @return Number of columns.
		__INLINE enumerator Cols() const {return ncols;}
		/// Create submatrix for a matrix.
		/// @param rM Reference to the matrix that stores elements.
		/// @param num_rows Number of rows in the larger matrix.
		/// @param num_cols Number of columns in the larger matrix.
		/// @param first_row Offset for row index in the larger matrix.
		/// @param first_column Offset for column index in the larger matrix.
		BlockOfMatrix(AbstractMatrix<Var> & rM, enumerator num_rows, enumerator num_cols, enumerator offset_row, enumerator offset_col) : M(&rM), nrows(num_rows), ncols(num_cols), orow(offset_row), ocol(offset_col)
		{}
		BlockOfMatrix(const BlockOfMatrix & b) : M(b.M), nrows(b.nrows), ncols(b.ncols), orow(b.orow), ocol(b.ocol) {}
		/// Assign matrix of another type to the block of matrix.
		/// \warning Only part of the matrix related to non-empty block is copied, other entries are ignored.
		/// @param other Another matrix of different type.
		/// @return Reference to current submatrix.
		template<typename typeB>
		BlockOfMatrix & operator =(AbstractMatrix<typeB> const & other)
		{
			assert( Cols() == other.Cols() );
			assert( Rows() == other.Rows() );
			for(enumerator i = orow; i < orow+M->Rows(); ++i)
				for(enumerator j = ocol; j < ocol+M->Cols(); ++j)
					assign((*M)(i-orow,j-ocol),other(i,j));
			return *this;
		}
		/// Assign matrix of another type to the block of matrix.
		/// \warning Only part of the matrix related to non-empty block is copied, other entries are ignored.
		/// @return Reference to current submatrix.
		template<typename typeB>
		BlockOfMatrix & operator =(BlockOfMatrix<typeB> const & other)
		{
			assert( Cols() == other.Cols() );
			assert( Rows() == other.Rows() );
			for(enumerator i = orow; i < orow+M->Rows(); ++i)
				for(enumerator j = ocol; j < ocol+M->Cols(); ++j)
					assign((*M)(i-orow,j-ocol),other(i,j));
			return *this;
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		__INLINE Var & operator()(enumerator i, enumerator j)
		{
			static Var zero = 0;
			assert(i < Rows());
			assert(j < Cols());
			if( i < orow || i >= orow+M->Rows() || j < ocol || j >= ocol+M->Cols() )
				return zero;
			else
				return (*M)(i-orow,j-ocol);
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			assert(i < Rows());
			assert(j < Cols());
			if( i < orow || i >= orow+M->Rows() || j < ocol || j >= ocol+M->Cols() )
				return Var(0.0);
			else
				return (*M)(i-orow,j-ocol);
		}
		/// Convert block of matrix into matrix.
		/// Note, that modifying returned matrix does
		/// not affect elements of the matrix or original matrix
		/// used to create block of matrix.
		/// @return Matrix with same entries as block of matrix including empty elements.
		::INMOST::Matrix<Var> MakeMatrix()
		{
			::INMOST::Matrix<Var> ret(Rows(),Cols());
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i,j) = (*this)(i,j);
			return ret;
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. BlockOfMatrix cannot change it's size,
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
            (void)cols; (void)rows;
		}
	};

	/// This class allows to address a matrix as a block of an empty matrix of larger size.
	template<typename Var>
	class ConstBlockOfMatrix : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var> * M;
		enumerator nrows; //< Number of rows in larger matrix.
		enumerator ncols; //< Number of colums in larger matrix.
		enumerator orow; //< Row offset in larger matrix.
		enumerator ocol; //< Column offset in larger matrix.
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return M->TrivialArguments(); };
		/// Number of rows in submatrix.
		/// @return Number of rows.
		__INLINE enumerator Rows() const {return nrows;}
		/// Number of columns in submatrix.
		/// @return Number of columns.
		__INLINE enumerator Cols() const {return ncols;}
		/// Create submatrix for a matrix.
		/// @param rM Reference to the matrix that stores elements.
		/// @param num_rows Number of rows in the larger matrix.
		/// @param num_cols Number of columns in the larger matrix.
		/// @param first_row Offset for row index in the larger matrix.
		/// @param first_column Offset for column index in the larger matrix.
		ConstBlockOfMatrix(const AbstractMatrixReadOnly<Var> & rM, enumerator num_rows, enumerator num_cols, enumerator offset_row, enumerator offset_col) : M(&rM), nrows(num_rows), ncols(num_cols), orow(offset_row), ocol(offset_col)
		{}
		ConstBlockOfMatrix(const ConstBlockOfMatrix & b) : M(b.M), nrows(b.nrows), ncols(b.ncols), orow(b.orow), ocol(b.ocol) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			assert(i < Rows());
			assert(j < Cols());
			if( i < orow || i >= orow+M->Rows() || j < ocol || j >= ocol+M->Cols() )
				return Var(0.0);
			else
				return (*M)(i-orow,j-ocol);
		}
		/// Convert block of matrix into matrix.
		/// Note, that modifying returned matrix does
		/// not affect elements of the matrix or original matrix
		/// used to create block of matrix.
		/// @return Matrix with same entries as block of matrix including empty elements.
		::INMOST::Matrix<Var> MakeMatrix()
		{
			::INMOST::Matrix<Var> ret(Rows(),Cols());
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret(i,j) = (*this)(i,j);
			return ret;
		}

	};

	template<typename Var>
	class MatrixUnit : public AbstractMatrixReadOnly< Var >
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		enumerator n;
		Var c;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		MatrixUnit(enumerator n, const Var & c = Var(1.0)) :n(n), c(c) {}
		MatrixUnit(const MatrixUnit& b) : n(b.n), c(b.c) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return n; }
		__INLINE Var operator()(enumerator i, enumerator j) const { return i == j ? c : Var(0.0); }
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			GetInterval(c, beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixRow : public AbstractMatrixReadOnly< Var >
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		enumerator n;
		Var c;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		MatrixRow(enumerator n, const Var& c = Var(1.0)) :n(n), c(c) {}
		MatrixRow(const MatrixRow& b) : n(b.n), c(b.c) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return 1; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return n; }
		__INLINE Var operator()(enumerator i, enumerator j) const { return c; }
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			GetInterval(c, beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixCol : public AbstractMatrixReadOnly< Var >
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		enumerator n;
		Var c;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		MatrixCol(enumerator n, const Var& c = Var(1.0)) :n(n), c(c) {}
		MatrixCol(const MatrixCol& b) : n(b.n), c(b.c) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return 1; }
		__INLINE Var operator()(enumerator i, enumerator j) const { return c; }
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			GetInterval(c, beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixDiag : public AbstractMatrixReadOnly< Var >
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		enumerator n;
		Var* diag;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		MatrixDiag(Var* diag, enumerator n) : diag(diag), n(n) {}
		MatrixDiag(const MatrixDiag& b) : diag(b.diag), n(b.n) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return n; }
		__INLINE Var operator()(enumerator i, enumerator j) const { return i == j ? diag[i] : Var(0.0); }
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			for(enumerator k = 0; k < n; ++k)
				GetInterval(diag[k], beg, end, cnt);
		}
	};

	template<typename VarA, typename VarB>
	class MatrixSum : public AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type >
	{
	public:
		using AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type >::operator();
		typedef typename AbstractMatrixReadOnly<typename Promote<VarA,VarB>::type >::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const AbstractMatrixReadOnly<VarB>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return false; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixSum(const AbstractMatrixReadOnly<VarA>& rA, const AbstractMatrixReadOnly<VarB>& rB)
			: A(&rA), B(&rB)
		{
			assert(A->Rows() == B->Rows());
			assert(A->Cols() == B->Cols());
		}
		MatrixSum(const MatrixSum& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE typename Promote<VarA, VarB>::type operator()(enumerator i, enumerator j) const
		{
			return (*A)(i,j) + (*B)(i,j);
		}

		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename VarA, typename VarB>
	class MatrixDifference : public AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type >
	{
	public:
		using AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type >::operator();
		typedef typename AbstractMatrixReadOnly<typename Promote<VarA, VarB>::type>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const AbstractMatrixReadOnly<VarB>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return false; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixDifference(const AbstractMatrixReadOnly<VarA>& rA, const AbstractMatrixReadOnly<VarB>& rB)
			: A(&rA), B(&rB)
		{
			assert(A->Rows() == B->Rows());
			assert(A->Cols() == B->Cols());
		}
		MatrixDifference(const MatrixDifference& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE typename Promote<VarA, VarB>::type operator()(enumerator i, enumerator j) const
		{
			return (*A)(i, j) - (*B)(i, j);
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class ConstMatrixTranspose : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var>* A;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Cols(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Rows(); }
		ConstMatrixTranspose(const AbstractMatrixReadOnly<Var>& rA)
			: A(&rA) {}
		ConstMatrixTranspose(const ConstMatrixTranspose& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const { return (*A)(j, i); }

		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class ConstMatrixConjugateTranspose : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var>* A;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Cols(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Rows(); }
		ConstMatrixConjugateTranspose(const AbstractMatrixReadOnly<Var>& rA)
			: A(&rA) {}
		ConstMatrixConjugateTranspose(const ConstMatrixConjugateTranspose& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const { return conj((*A)(j, i)); }
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class ConstMatrixConjugate : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var>* A;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		ConstMatrixConjugate(const AbstractMatrixReadOnly<Var>& rA)
			: A(&rA) {}
		ConstMatrixConjugate(const ConstMatrixConjugate& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const { return conj((*A)(i, j)); }
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixUnaryMinus : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var>* A;
	public:
		/// Consider unary minus is trivial.
		bool TrivialArguments() const { return A->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixUnaryMinus(const AbstractMatrixReadOnly<Var>& rA)
			: A(&rA) {}
		MatrixUnaryMinus(const MatrixUnaryMinus& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const { return -(*A)(i, j); }
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixTranspose : public AbstractMatrix<Var>
	{
	public:
		using AbstractMatrix<Var>::operator();
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		AbstractMatrix<Var>* A;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Cols(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Rows(); }
		MatrixTranspose(AbstractMatrix<Var>& rA)
			: A(&rA) {}
		MatrixTranspose(const MatrixTranspose& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const { return (*A)(j, i); }
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var& operator()(enumerator i, enumerator j) { return (*A)(j, i); }
		/// This is a stub function to fulfill abstract
		/// inheritance. BlockOfMatrix cannot change it's size,
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class ConstMatrixConcatRows : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var>* A;
		const AbstractMatrixReadOnly<Var>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments() && B->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows() + B->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		ConstMatrixConcatRows(const AbstractMatrixReadOnly<Var>& rA, const AbstractMatrixReadOnly<Var>& rB)
			: A(&rA), B(&rB) {
			assert(A->Cols() == B->Cols());
		}
		ConstMatrixConcatRows(const ConstMatrixConcatRows& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			return i < A->Rows() ? (*A)(i, j) : (*B)(i - A->Rows(), j);
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename VarA, typename VarB, typename VarR>
	class ConstMatrixConcatRows2 : public AbstractMatrixReadOnly<VarR>
	{
	public:
		using AbstractMatrixReadOnly<VarR>::operator();
		typedef typename AbstractMatrixReadOnly<VarR>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const AbstractMatrixReadOnly<VarB>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments() && B->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows() + B->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		ConstMatrixConcatRows2(const AbstractMatrixReadOnly<VarA>& rA, const AbstractMatrixReadOnly<VarB>& rB)
			: A(&rA), B(&rB) {
			assert(A->Cols() == B->Cols());
		}
		ConstMatrixConcatRows2(const ConstMatrixConcatRows2& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE VarR operator()(enumerator i, enumerator j) const
		{
			if (i < A->Rows())
				return (*A)(i, j);
			else
				return (*B)(i - A->Rows(), j);
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixConcatRows : public AbstractMatrix<Var>
	{
	public:
		using AbstractMatrix<Var>::operator();
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		AbstractMatrix<Var>* A;
		AbstractMatrix<Var>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows() + B->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixConcatRows(AbstractMatrix<Var>& rA, AbstractMatrix<Var>& rB)
			: A(&rA), B(&rB) { assert(A->Cols() == B->Cols());	}
		MatrixConcatRows(const MatrixConcatRows& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			return i < A->Rows() ? (*A)(i, j) : (*B)(i - A->Rows(), j);
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var& operator()(enumerator i, enumerator j)
		{
			return i < A->Rows() ? (*A)(i, j) : (*B)(i - A->Rows(), j);
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. BlockOfMatrix cannot change it's size,
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class ConstMatrixConcatCols : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var>* A;
		const AbstractMatrixReadOnly<Var>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments() && B->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols() + B->Cols(); }
		ConstMatrixConcatCols(const AbstractMatrixReadOnly<Var>& rA, const AbstractMatrixReadOnly<Var>& rB)
			: A(&rA), B(&rB) {
			assert(A->Rows() == B->Rows());
		}
		ConstMatrixConcatCols(const ConstMatrixConcatCols& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			return j < A->Cols() ? (*A)(i, j) : (*B)(i, j - A->Cols());
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename VarA, typename VarB, typename VarR>
	class ConstMatrixConcatCols2 : public AbstractMatrixReadOnly<VarR>
	{
	public:
		using AbstractMatrixReadOnly<VarR>::operator();
		typedef typename AbstractMatrixReadOnly<VarR>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const AbstractMatrixReadOnly<VarB>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments() && B->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols() + B->Cols(); }
		ConstMatrixConcatCols2(const AbstractMatrixReadOnly<VarA>& rA, const AbstractMatrixReadOnly<VarB>& rB)
			: A(&rA), B(&rB) {
			assert(A->Rows() == B->Rows());
		}
		ConstMatrixConcatCols2(const ConstMatrixConcatCols2& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE VarR operator()(enumerator i, enumerator j) const
		{
			if (j < A->Cols())
				return (*A)(i, j);
			else
				return (*B)(i, j - A->Cols());
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixConcatCols : public AbstractMatrix<Var>
	{
	public:
		using AbstractMatrix<Var>::operator();
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		AbstractMatrix<Var>* A;
		AbstractMatrix<Var>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols() + B->Cols(); }
		MatrixConcatCols(AbstractMatrix<Var>& rA, AbstractMatrix<Var>& rB)
			: A(&rA), B(&rB) { assert(A->Rows() == B->Rows()); }
		MatrixConcatCols(const MatrixConcatCols& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			return j < A->Cols() ? (*A)(i, j) : (*B)(i, j - A->Cols());
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var& operator()(enumerator i, enumerator j)
		{
			return j < A->Cols() ? (*A)(i, j) : (*B)(i, j - A->Cols());
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. BlockOfMatrix cannot change it's size,
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class ConstMatrixRepack : public AbstractMatrixReadOnly<Var>
	{
	public:
		using AbstractMatrixReadOnly<Var>::operator();
		typedef typename AbstractMatrixReadOnly<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<Var>* A;
		enumerator n, m;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return A->TrivialArguments(); };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return m; }
		ConstMatrixRepack(const AbstractMatrixReadOnly<Var>& rA, enumerator n, enumerator m)
			: A(&rA), n(n), m(m) { assert(A->Rows() * A->Cols() == n * m); }
		ConstMatrixRepack(const ConstMatrixRepack& b) : A(b.A), n(b.n), m(b.m) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			enumerator ind = i * m + j;
			return (*A)(ind / A->Cols(), ind % A->Cols());
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename Var>
	class MatrixRepack : public AbstractMatrix<Var>
	{
	public:
		using AbstractMatrix<Var>::operator();
		typedef typename AbstractMatrix<Var>::enumerator enumerator; //< Integer type for indexes.
	private:
		AbstractMatrix<Var>* A;
		enumerator n, m;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return m; }
		MatrixRepack(AbstractMatrix<Var>& rA, enumerator n, enumerator m)
			: A(&rA), n(n), m(m) {
			assert(A->Rows() * A->Cols() == n * m);
		}
		MatrixRepack(const MatrixRepack& b) : A(b.A), n(b.n), m(b.m) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var operator()(enumerator i, enumerator j) const
		{
			enumerator p = i * m + j;
			return (*A)(p / A->Cols(), p % A->Cols());
			//return A->data()[i * m + j];
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Var & operator()(enumerator i, enumerator j)
		{
			enumerator p = i * m + j;
			return (*A)(p / A->Cols(), p % A->Cols());
			//return A->data()[i * m + j];
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. BlockOfMatrix cannot change it's size,
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
		}
	};

	template<typename VarA, typename VarB, typename VarR>
	class MatrixMul : public AbstractMatrixReadOnly< VarR >
	{
	public:
		using AbstractMatrixReadOnly< VarR >::operator();
		typedef typename AbstractMatrixReadOnly< VarR >::enumerator enumerator; //< Integer type for indexes.
	private:
		Matrix<VarR> M;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return M.Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return M.Cols(); }
		MatrixMul(const AbstractMatrixReadOnly<VarA>& rA, const AbstractMatrixReadOnly<VarB>& rB)
		{
			assert(rA.Cols() == rB.Rows());
			const AbstractMatrixReadOnly<VarA>* pA = &rA;
			const AbstractMatrixReadOnly<VarB>* pB = &rB;
			if (!pA->TrivialArguments())
			{
				static thread_private< Matrix<VarA> > tmpA;
				*tmpA = rA;
				pA = &(*tmpA);
			}
			if (!pB->TrivialArguments())
			{
				static thread_private< Matrix<VarB> > tmpB;
				*tmpB = rB;
				pB = &(*tmpB);
			}
			M.Resize(rA.Rows(), rB.Cols());
			for (enumerator i = 0; i < pA->Rows(); ++i)
			{
				for (enumerator  k = 0; k < pB->Cols(); ++k)
					M(i, k) = (*pA)(i, 0) * (*pB)(0, k);
				for (enumerator j = 1; j < pA->Cols(); ++j)
					for (enumerator k = 0; k < pB->Cols(); ++k)
						M(i, k) += (*pA)(i, j) * (*pB)(j, k);
			}
		}
		MatrixMul(const MatrixMul& b) : M(b.M) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE VarR operator()(enumerator i, enumerator j) const
		{
			return M(i, j);
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			M.GetMatrixInterval(beg, end, cnt);
		}
	};
#if defined(USE_AUTODIFF)
	template<>
	class MatrixMul<INMOST_DATA_REAL_TYPE, variable, Promote<INMOST_DATA_REAL_TYPE,variable>::type > : public AbstractMatrixReadOnly< Promote<INMOST_DATA_REAL_TYPE, variable>::type >
	{
	public:
		using AbstractMatrixReadOnly< Promote<INMOST_DATA_REAL_TYPE, variable>::type >::operator();
		typedef typename AbstractMatrixReadOnly< Promote<INMOST_DATA_REAL_TYPE, variable>::type >::enumerator enumerator; //< Integer type for indexes.
	private:
		Matrix<Promote<INMOST_DATA_REAL_TYPE, variable>::type> M;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return M.Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return M.Cols(); }
		MatrixMul(const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE>& rA, const AbstractMatrixReadOnly<variable>& rB)
		{
			assert(rA.Cols() == rB.Rows());
			const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE>* pA = &rA;
			const AbstractMatrixReadOnly<variable>* pB = &rB;
			if (!pA->TrivialArguments())
			{
				static thread_private< Matrix<INMOST_DATA_REAL_TYPE> > tmpA;
				*tmpA = rA;
				pA = &(*tmpA);
			}
			if (!pB->TrivialArguments())
			{
				static thread_private< Matrix<variable> > tmpB;
				*tmpB = rB;
				pB = &(*tmpB);
			}
			M.Resize(rA.Rows(), rB.Cols());
			INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
			pB->GetMatrixInterval(beg, end, cnt);
			if (cnt >= CNT_USE_MERGER)
			//if(USE_MERGER)
			{
				//Sparse::RowMerger merger;// (beg, end);
				merger->Resize(beg, end);
				for (enumerator i = 0; i < pA->Rows(); ++i)
				{
					for (enumerator j = 0; j < pB->Cols(); ++j)
					{
						INMOST_DATA_REAL_TYPE value = 0.0;
						for (enumerator k = 0; k < pA->Cols(); ++k)
						{
							value += (*pA)(i, k) * (*pB)(k, j).GetValue();
							merger->AddRow((*pA)(i, k), (*pB)(k, j).GetRow());
						}
						M(i, j).SetValue(value);
						merger->RetriveRow(M(i, j).GetRow());
						merger->Clear();
					}
				}
			}
			else
			{
				for (enumerator i = 0; i < pA->Rows(); ++i)
				{
					for (enumerator k = 0; k < pB->Cols(); ++k)
						M(i, k) = (*pA)(i, 0) * (*pB)(0, k);
					for (enumerator j = 1; j < pA->Cols(); ++j)
						for (enumerator k = 0; k < pB->Cols(); ++k)
							M(i, k) += (*pA)(i, j) * (*pB)(j, k);
				}
			}
		}
		MatrixMul(const MatrixMul& b) : M(b.M) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Promote<INMOST_DATA_REAL_TYPE, variable>::type operator()(enumerator i, enumerator j) const
		{
			return M(i, j);
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			M.GetMatrixInterval(beg, end, cnt);
		}
	};

	template<>
	class MatrixMul<variable, INMOST_DATA_REAL_TYPE, Promote<variable, INMOST_DATA_REAL_TYPE>::type > : public AbstractMatrixReadOnly< Promote<variable, INMOST_DATA_REAL_TYPE>::type >
	{
	public:
		using AbstractMatrixReadOnly< Promote<variable, INMOST_DATA_REAL_TYPE>::type >::operator();
		typedef typename AbstractMatrixReadOnly< Promote<variable, INMOST_DATA_REAL_TYPE>::type >::enumerator enumerator; //< Integer type for indexes.
	private:
		Matrix<Promote<variable, INMOST_DATA_REAL_TYPE>::type> M;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return M.Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return M.Cols(); }
		MatrixMul(const AbstractMatrixReadOnly<variable>& rA, const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE>& rB)
		{
			assert(rA.Cols() == rB.Rows());
			const AbstractMatrixReadOnly<variable>* pA = &rA;
			const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE>* pB = &rB;
			if (!pA->TrivialArguments())
			{
				static thread_private< Matrix<variable> > tmpA;
				*tmpA = rA;
				pA = &(*tmpA);
			}
			if (!pB->TrivialArguments())
			{
				static thread_private< Matrix<INMOST_DATA_REAL_TYPE> > tmpB;
				*tmpB = rB;
				pB = &(*tmpB);
			}
			M.Resize(rA.Rows(), rB.Cols());
			INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
			pA->GetMatrixInterval(beg, end, cnt);
			if (cnt >= CNT_USE_MERGER)
			//if(USE_MERGER)
			{
				//Sparse::RowMerger merger;// (beg, end);
				merger->Resize(beg, end);
				for (enumerator i = 0; i < pA->Rows(); ++i)
				{
					for (enumerator j = 0; j < pB->Cols(); ++j)
					{
						INMOST_DATA_REAL_TYPE value = 0.0;
						for (enumerator k = 0; k < pA->Cols(); ++k)
						{
							value += (*pA)(i, k).GetValue() * (*pB)(k, j);
							merger->AddRow((*pB)(k, j), (*pA)(i, k).GetRow());
						}
						M(i, j).SetValue(value);
						merger->RetriveRow(M(i, j).GetRow());
						merger->Clear();
					}
				}
			}
			else
			{
				for (enumerator i = 0; i < pA->Rows(); ++i)
				{
					for (enumerator k = 0; k < pB->Cols(); ++k)
						M(i, k) = (*pA)(i, 0) * (*pB)(0, k);
					for (enumerator j = 1; j < pA->Cols(); ++j)
						for (enumerator k = 0; k < pB->Cols(); ++k)
							M(i, k) += (*pA)(i, j) * (*pB)(j, k);
				}
			}
		}
		MatrixMul(const MatrixMul& b) : M(b.M) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Promote<variable, INMOST_DATA_REAL_TYPE>::type operator()(enumerator i, enumerator j) const
		{
			return M(i, j);
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			M.GetMatrixInterval(beg, end, cnt);
		}
	};

	template<>
	class MatrixMul<variable, variable, Promote<variable, variable>::type > : public AbstractMatrixReadOnly< Promote<variable, variable>::type >
	{
	public:
		using AbstractMatrixReadOnly< Promote<variable, variable>::type >::operator();
		typedef typename AbstractMatrixReadOnly< Promote<variable, variable>::type >::enumerator enumerator; //< Integer type for indexes.
	private:
		Matrix<Promote<variable, variable>::type> M;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return true; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return M.Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return M.Cols(); }
		MatrixMul(const AbstractMatrixReadOnly<variable>& rA, const AbstractMatrixReadOnly<variable>& rB)
		{
			assert(rA.Cols() == rB.Rows());
			const AbstractMatrixReadOnly<variable>* pA = &rA;
			const AbstractMatrixReadOnly<variable>* pB = &rB;
			if (!pA->TrivialArguments())
			{
				static thread_private< Matrix<variable> > tmpA;
				*tmpA = rA;
				pA = &(*tmpA);
			}
			if (!pB->TrivialArguments())
			{
				static thread_private< Matrix<variable> > tmpB;
				*tmpB = rB;
				pB = &(*tmpB);
			}
			M.Resize(rA.Rows(), rB.Cols());
			INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
			pA->GetMatrixInterval(beg, end, cnt);
			pB->GetMatrixInterval(beg, end, cnt);
			if (cnt >= CNT_USE_MERGER)
			//if(USE_MERGER)
			{
				//Sparse::RowMerger merger;// (beg, end);
				merger->Resize(beg, end);
				for (enumerator i = 0; i < pA->Rows(); ++i)
				{
					for (enumerator j = 0; j < pB->Cols(); ++j)
					{
						INMOST_DATA_REAL_TYPE value = 0.0;
						for (enumerator k = 0; k < pA->Cols(); ++k)
						{
							value += (*pA)(i, k).GetValue() * (*pB)(k, j).GetValue();
							merger->AddRow((*pA)(i, k).GetValue(), (*pB)(k, j).GetRow());
							merger->AddRow((*pB)(k, j).GetValue(), (*pA)(i, k).GetRow());
						}
						M(i, j).SetValue(value);
						merger->RetriveRow(M(i, j).GetRow());
						merger->Clear();
					}
				}
			}
			else
			{
				for (enumerator i = 0; i < pA->Rows(); ++i)
				{
					for (enumerator k = 0; k < pB->Cols(); ++k)
						M(i, k) = (*pA)(i, 0) * (*pB)(0, k);
					for (enumerator j = 1; j < pA->Cols(); ++j)
						for (enumerator k = 0; k < pB->Cols(); ++k)
							M(i, k) += (*pA)(i, j) * (*pB)(j, k);
				}
			}
		}
		MatrixMul(const MatrixMul& b) : M(b.M) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Promote<variable, variable>::type operator()(enumerator i, enumerator j) const
		{
			return M(i, j);
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			M.GetMatrixInterval(beg, end, cnt);
		}
	};
#endif //USE_AUTODIFF
	//template<typename VarA, typename VarB, typename VarR>
	//thread_private< Matrix<VarR> > MatrixMul<VarA, VarB, VarR>::M;

	template<typename VarA, typename VarB, typename VarR>
	class MatrixMulCoef : public AbstractMatrixReadOnly< VarR >
	{
	public:
		using AbstractMatrixReadOnly< VarR >::operator();
		typedef typename AbstractMatrixReadOnly< VarR >::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const VarB* coef;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return false; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixMulCoef(const AbstractMatrixReadOnly<VarA>& rA, const VarB& rcoef)
			: A(&rA), coef(&rcoef) {}
		MatrixMulCoef(const MatrixMulCoef& b) : A(b.A), coef(b.coef) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE VarR operator()(enumerator i, enumerator j) const
		{
			VarR ret;
			assign(ret, (*A)(i, j) * (*coef));
			return ret;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			GetInterval(*coef, beg, end, cnt);
		}
	};

	template<typename VarA, typename VarB, typename VarR>
	class MatrixDivCoef : public AbstractMatrixReadOnly< VarR >
	{
	public:
		using AbstractMatrixReadOnly< VarR >::operator();
		typedef typename AbstractMatrixReadOnly< VarR >::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const VarB* coef;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return false; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixDivCoef(const AbstractMatrixReadOnly<VarA>& rA, const VarB& rcoef)
			: A(&rA), coef(&rcoef) {}
		MatrixDivCoef(const MatrixDivCoef& b)
			: A(b.A), coef(b.coef) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE VarR operator()(enumerator i, enumerator j) const
		{
			VarR ret;
			assign(ret, (*A)(i, j) / (*coef));
			return ret;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			GetInterval(*coef, beg, end, cnt);
		}
	};
#if defined(USE_AUTODIFF)
	template<typename VarA, typename VarB, typename VarR>
	class MatrixMulShellCoef : public AbstractMatrixReadOnly< VarR >
	{
	public:
		using AbstractMatrixReadOnly< VarR >::operator();
		typedef typename AbstractMatrixReadOnly< VarR >::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const variable coef;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return false; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixMulShellCoef(const AbstractMatrixReadOnly<VarA>& rA, const VarB& rcoef)
			: A(&rA), coef(rcoef) {}
		MatrixMulShellCoef(const MatrixMulShellCoef& b) : A(b.A), coef(b.coef) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE VarR operator()(enumerator i, enumerator j) const
		{
			VarR ret;
			assign(ret, (*A)(i, j) * coef);
			return ret;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			GetInterval(coef, beg, end, cnt);
		}
	};

	template<typename VarA, typename VarB, typename VarR>
	class MatrixDivShellCoef : public AbstractMatrixReadOnly< VarR >
	{
	public:
		using AbstractMatrixReadOnly< VarR >::operator();
		typedef typename AbstractMatrixReadOnly< VarR >::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const variable coef;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return false; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixDivShellCoef(const AbstractMatrixReadOnly<VarA>& rA, const VarB& rcoef)
			: A(&rA), coef(rcoef) {}
		MatrixDivShellCoef(const MatrixDivShellCoef& b)
			: A(b.A), coef(b.coef) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE VarR operator()(enumerator i, enumerator j) const
		{
			VarR ret;
			assign(ret, (*A)(i, j) / coef);
			return ret;
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
			A->GetMatrixInterval(beg, end, cnt);
			GetInterval(coef, beg, end, cnt);
		}
	};
#endif //USE_AUTODIFF
	template<typename VarA, typename VarB>
	class KroneckerProduct : public AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type >
	{
	public:
		using AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type >::operator();
		typedef typename AbstractMatrixReadOnly<typename Promote<VarA, VarB>::type>::enumerator enumerator; //< Integer type for indexes.
	private:
		const AbstractMatrixReadOnly<VarA>* A;
		const AbstractMatrixReadOnly<VarB>* B;
	public:
		/// Check that this matrix does not require calculations for element access.
		/// The function is used during multiplication.
		bool TrivialArguments() const { return false; };
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows() * B->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols() * B->Cols(); }
		KroneckerProduct(const AbstractMatrixReadOnly<VarA>& rA, const AbstractMatrixReadOnly<VarB>& rB)
			: A(&rA), B(&rB) {}
		KroneckerProduct(const KroneckerProduct& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE typename Promote<VarA, VarB>::type operator()(enumerator i, enumerator j) const
		{
			return (*A)(i / B->Rows(), j / B->Cols()) * (*B)(i % B->Rows(), j % B->Cols());
		}
		__INLINE void GetMatrixInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end, INMOST_DATA_ENUM_TYPE& cnt) const
		{
            A->GetMatrixInterval(beg, end, cnt);
			B->GetMatrixInterval(beg, end, cnt);
		}
	};

    /*
    template<>
    typename SelfPromote<INMOST_DATA_CPLX_TYPE>::type AbstractMatrixReadOnly<INMOST_DATA_CPLX_TYPE>::FrobeniusNorm() const
    {
        typename SelfPromote<INMOST_DATA_CPLX_TYPE>::type ret = 0;
        for (enumerator i = 0; i < Rows(); ++i)
            for (enumerator j = 0; j < Cols(); ++j)
                ret += (*this)(i, j) * std::conj((*this)(i, j));
        return sqrt(fabs(ret));
    }
    */


	/*
	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::Transpose() const
	{
		Matrix<Var> ret(Cols(),Rows());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(j,i) = (*this)(i,j);
		return ret;
	}
	*/

	template<typename Var>
	void
	AbstractMatrix<Var>::Swap(AbstractMatrix<Var> & b)
	{
		Matrix<Var> tmp = b;
		b = (*this);
		(*this) = tmp;
	}

	/*
	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::operator-() const
	{
		Matrix<Var> ret(Rows(),Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = -(*this)(i,j);
		return ret;
	}
	*/
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var>::CrossProduct(const AbstractMatrixReadOnly<typeB> & B) const
	{
		const AbstractMatrixReadOnly<Var> & A = *this;
		assert(A.Rows() == 3);
		assert(A.Cols() == 1);
		assert(B.Rows() == 3);
		Matrix<typename Promote<Var,typeB>::type> ret(3,B.Cols()); //check RVO
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
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var>::Transform(const AbstractMatrixReadOnly<typeB> & B) const
	{
		const AbstractMatrixReadOnly<Var> & A = *this;
		assert(A.Rows() == 3);
		assert(A.Cols() == 1);
		assert(B.Rows() == 3);
		assert(B.Cols() == 1);
		typedef typename Promote<Var,typeB>::type var_t;
		Matrix<var_t> Q(3,3);
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
	//Matrix<typename Promote<Var,typeB>::type>
	MatrixDifference<Var,typeB>
	AbstractMatrixReadOnly<Var>::operator-(const AbstractMatrixReadOnly<typeB> & other) const
	{
		return MatrixDifference<Var, typeB>(*this, other);
		/*
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = (*this)(i,j)-other(i,j);
		return ret;
		*/
	}

	template<typename Var>
	template<typename typeB>
	AbstractMatrix<Var> &
	AbstractMatrix<Var>::operator-=(const AbstractMatrixReadOnly<typeB> & other)
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
	//Matrix<typename Promote<Var,typeB>::type>
	MatrixSum<Var,typeB>
	AbstractMatrixReadOnly<Var>::operator+(const AbstractMatrixReadOnly<typeB> & other) const
	{
		return MatrixSum<Var, typeB>(*this, other);
		/*
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				ret(i,j) = (*this)(i,j)+other(i,j);
		return ret;
		*/
	}

	template<typename Var>
	template<typename typeB>
	AbstractMatrix<Var> &
	AbstractMatrix<Var>::operator+=(const AbstractMatrixReadOnly<typeB> & other)
	{
		assert(Rows() == other.Rows());
		assert(Cols() == other.Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				assign((*this)(i,j),(*this)(i,j)+other(i,j));
		return *this;
	}

	/*
#if defined(USE_AUTODIFF)
	template<>
	template<>
	__INLINE Matrix<Promote<INMOST_DATA_REAL_TYPE,variable>::type>
	AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE>::operator*<variable>(const AbstractMatrixReadOnly<variable> & other) const
	{
		//~ std::cout << __FUNCTION__ << std::endl;
		assert(Cols() == other.Rows());
		Matrix<Promote<INMOST_DATA_REAL_TYPE,variable>::type> ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				if( CheckCurrentAutomatizator() )
				{
					Sparse::RowMerger & merger = GetCurrentMerger();
					double value = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
					{
						value += (*this)(i,k)*other(k,j).GetValue();
						merger.AddRow((*this)(i,k),other(k,j).GetRow());
					}
					ret(i,j).SetValue(value);
					merger.RetriveRow(ret(i,j).GetRow());
					merger.Clear();
				}
				else
				{
					//~ typename Promote<INMOST_DATA_REAL_TYPE,variable>::type tmp = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
						ret(i,j) += (*this)(i,k)*other(k,j);
					//~ ret(i,j) = tmp;
				}
			}
		}
		return ret;
	}


	template<>
	template<>
	__INLINE Matrix<Promote<variable,INMOST_DATA_REAL_TYPE>::type>
	AbstractMatrixReadOnly<variable>::operator*<INMOST_DATA_REAL_TYPE>(const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE> & other) const
	{
		//~ std::cout << __FUNCTION__ << std::endl;
		assert(Cols() == other.Rows());
		Matrix<Promote<variable,INMOST_DATA_REAL_TYPE>::type> ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				if( CheckCurrentAutomatizator() )
				{
					Sparse::RowMerger & merger = GetCurrentMerger();
					double value = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
					{
						value += (*this)(i,k).GetValue()*other(k,j);
						merger.AddRow(other(k,j),(*this)(i,k).GetRow());
					}
					ret(i,j).SetValue(value);
					merger.RetriveRow(ret(i,j).GetRow());
					merger.Clear();
				}
				else
				{
					//~ typename Promote<INMOST_DATA_REAL_TYPE,variable>::type tmp = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
						ret(i,j) += (*this)(i,k)*other(k,j);
					//~ ret(i,j) = tmp;
				}
			}
		}
		return ret;
	}

	template<>
	template<>
	__INLINE Matrix<Promote<variable,variable>::type>
	AbstractMatrixReadOnly<variable>::operator*<variable>(const AbstractMatrixReadOnly<variable> & other) const
	{
		//~ std::cout << __FUNCTION__ << std::endl;
		assert(Cols() == other.Rows());
		Matrix<Promote<variable,variable>::type> ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				if( CheckCurrentAutomatizator() )
				{
					Sparse::RowMerger & merger = GetCurrentMerger();
					double value = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
					{
						value += (*this)(i,k).GetValue()*other(k,j).GetValue();
						merger.AddRow(other(k,j).GetValue(),(*this)(i,k).GetRow());
						merger.AddRow((*this)(i,k).GetValue(),other(k,j).GetRow());
					}
					ret(i,j).SetValue(value);
					merger.RetriveRow(ret(i,j).GetRow());
					merger.Clear();
				}
				else
				{
					//~ typename Promote<INMOST_DATA_REAL_TYPE,variable>::type tmp = 0.0;
					for(enumerator k = 0; k < Cols(); ++k)
						ret(i,j) += (*this)(i,k)*other(k,j);
					//~ ret(i,j) = tmp;
				}
			}
		}
		return ret;
	}
#endif //USE_AUTODIFF
	*/
	template<typename Var>
	template<typename typeB>
	//Matrix<typename Promote<Var,typeB>::type>
	MatrixMul<Var, typeB, typename Promote<Var, typeB>::type>
	AbstractMatrixReadOnly<Var>::operator*(const AbstractMatrixReadOnly<typeB> & other) const
	{
		return MatrixMul<Var, typeB, typename Promote<Var, typeB>::type>(*this, other);
		/*
		assert(Cols() == other.Rows());
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),other.Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Cols(); ++j) //loop columns
			{
				//~ typename Promote<Var,typeB>::type tmp = 0.0;
				for(enumerator k = 0; k < Cols(); ++k)
					ret(i,j) += (*this)(i,k)*other(k,j);
				//~ ret(i,j) = tmp;
			}
		}
		return ret;
		*/
	}

#if defined(USE_AUTODIFF)
	template<>
	template<>
	__INLINE Promote<INMOST_DATA_REAL_TYPE,variable>::type
	AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE>::DotProduct<variable>(const AbstractMatrixReadOnly<variable> & other) const
	{
		assert(Cols() == other.Cols());
		assert(Rows() == other.Rows());
		Promote<INMOST_DATA_REAL_TYPE,variable>::type ret = 0.0;
		INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
		other.GetMatrixInterval(beg, end, cnt);
		if( cnt >= CNT_USE_MERGER )
		//if(USE_MERGER)
		{
			//Sparse::RowMerger merger;// (beg, end);
			merger->Resize(beg, end);
			double value = 0.0;
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
				{
					value += (*this)(i,j)*other(i,j).GetValue();
					merger->AddRow((*this)(i,j),other(i,j).GetRow());
				}
			ret.SetValue(value);
			merger->RetriveRow(ret.GetRow());
			merger->Clear();
		}
		else
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret += ((*this)(i,j))*other(i,j);
		}
		return ret;
	}

	template<>
	template<>
	__INLINE Promote<variable,INMOST_DATA_REAL_TYPE>::type
	AbstractMatrixReadOnly<variable>::DotProduct<INMOST_DATA_REAL_TYPE>(const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE> & other) const
	{
		assert(Cols() == other.Cols());
		assert(Rows() == other.Rows());
		Promote<variable,INMOST_DATA_REAL_TYPE>::type ret = 0.0;
		INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
		GetMatrixInterval(beg, end, cnt);
		if( cnt >= CNT_USE_MERGER )
		//if(USE_MERGER)
		{
			//Sparse::RowMerger merger;// (beg, end);
			merger->Resize(beg, end);
			double value = 0.0;
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
				{
					value += (*this)(i,j).GetValue()*other(i,j);
					merger->AddRow(other(i,j),(*this)(i,j).GetRow());
				}
			ret.SetValue(value);
			merger->RetriveRow(ret.GetRow());
			merger->Clear();
		}
		else
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret += ((*this)(i,j))*other(i,j);
		}
		return ret;
	}

	template<>
	template<>
	__INLINE Promote<variable,variable>::type
	AbstractMatrixReadOnly<variable>::DotProduct<variable>(const AbstractMatrixReadOnly<variable> & other) const
	{
		assert(Cols() == other.Cols());
		assert(Rows() == other.Rows());
		Promote<variable,variable>::type ret = 0.0;
		INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
		GetMatrixInterval(beg, end, cnt);
		other.GetMatrixInterval(beg, end, cnt);
		if( cnt >= CNT_USE_MERGER )
		//if(USE_MERGER)
		{
			//Sparse::RowMerger merger;// (beg, end);
			merger->Resize(beg, end);
			double value = 0.0;
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
				{
					value += (*this)(i,j).GetValue()*other(i,j).GetValue();
					merger->AddRow(other(i,j).GetValue(),(*this)(i,j).GetRow());
					merger->AddRow((*this)(i,j).GetValue(),other(i,j).GetRow());
				}
			ret.SetValue(value);
			merger->RetriveRow(ret.GetRow());
			merger->Clear();
		}
		else
		{
			for(enumerator i = 0; i < Rows(); ++i)
				for(enumerator j = 0; j < Cols(); ++j)
					ret += ((*this)(i,j))*other(i,j);
		}
		return ret;
	}
#endif //USE_AUTODIFF



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
	AbstractMatrix<Var>::operator*=(const AbstractMatrixReadOnly<typeB> & B)
	{
		(*this) = (*this)*B;
		return *this;
	}

	template<typename Var>
	template<typename typeB>
	//Matrix<typename Promote<Var,typeB>::type>
	MatrixMulCoef<Var, typeB, typename Promote<Var, typeB>::type>
	AbstractMatrixReadOnly<Var>::operator*(const typeB& coef) const
	{
		return MatrixMulCoef<Var, typeB, typename Promote<Var, typeB>::type>(*this, coef);
		/*
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols()); //check RVO
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				assign(ret(i,j),(*this)(i,j)*coef);
		return ret;
		*/
	}

#if defined(USE_AUTODIFF)
	template<typename Var>
	template<typename A>
	MatrixMulShellCoef<Var, shell_expression<A>, typename Promote<Var, variable>::type>
		AbstractMatrixReadOnly<Var>::operator*(shell_expression<A> const & coef) const
	{
		return MatrixMulShellCoef<Var, shell_expression<A>, typename Promote<Var, variable>::type>(*this, coef);
	}
#endif //USE_AUTODIFF

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
#if defined(USE_AUTODIFF)
	template<typename Var>
	template<typename A>
	MatrixDivShellCoef<Var, shell_expression<A>, typename Promote<Var, variable>::type>
	AbstractMatrixReadOnly<Var>::operator/(shell_expression<A> const & coef) const
	{
		return MatrixDivShellCoef<Var, shell_expression<A>, typename Promote<Var, variable>::type>(*this, coef);
	}
#endif //USE_AUTODIFF

	template<typename Var>
	template<typename typeB>
	//Matrix<typename Promote<Var,typeB>::type>
	MatrixDivCoef<Var, typeB, typename Promote<Var, typeB>::type>
		AbstractMatrixReadOnly<Var>::operator/(const typeB& coef) const
	{
		return MatrixDivCoef<Var, typeB, typename Promote<Var, typeB>::type>(*this, coef);
		/*
		Matrix<typename Promote<Var,typeB>::type> ret(Rows(),Cols());
		for(enumerator i = 0; i < Rows(); ++i)
			for(enumerator j = 0; j < Cols(); ++j)
				assign(ret(i,j),(*this)(i,j)/coef);
		return ret;
		*/
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
	//Matrix<typename Promote<Var,typeB>::type>
	KroneckerProduct<Var,typeB>
	AbstractMatrixReadOnly<Var>::Kronecker(const AbstractMatrixReadOnly<typeB> & other) const
	{
		return KroneckerProduct<Var, typeB>(*this, other);
		/*
		Matrix<typename Promote<Var,typeB>::type> ret(Rows()*other.Rows(),Cols()*other.Cols());
		for(enumerator i = 0; i < Rows(); ++i) //loop rows
		{
			for(enumerator j = 0; j < other.Rows(); ++j) //loop columns
			{
				for(enumerator k = 0; k < Cols(); ++k)
				{
					for(enumerator l = 0; l < other.Cols(); ++l)
					{
						assign(ret(i*other.Rows()+j,k*other.Cols()+l),other(j,l)*(*this)(i,k));
					}
				}
			}
		}
		return ret;
		*/
	}

	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::Invert(int * ierr) const
	{
		Matrix<Var> ret(Cols(),Rows());
		ret = Solve(MatrixUnit<Var>(Rows()),ierr);
		return ret;
	}

	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::CholeskyInvert(int * ierr) const
	{
		Matrix<Var> ret(Rows(),Rows());
		ret = CholeskySolve(MatrixUnit<Var>(Rows()),ierr);
		return ret;
	}


	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var>::CholeskySolve(const AbstractMatrixReadOnly<typeB> & B, int * ierr) const
	{
		const AbstractMatrixReadOnly<Var> & A = *this;
		assert(A.Rows() == A.Cols());
		assert(A.Rows() == B.Rows());
		enumerator n = A.Rows();
		enumerator l = B.Cols();
		Matrix<typename Promote<Var,typeB>::type> ret(B);
		static thread_private< SymmetricMatrix<Var> > L;// (A);
		*L = A;

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
			if( real_part((*L)(k,k)) < 0.0 )
			{
				if( ierr )
				{
					if( *ierr == -1 ) std::cout << "Negative diagonal pivot " << get_value((*L)(k,k)) << " row " << k << std::endl;
					*ierr = k+1;
				}
				else throw MatrixCholeskySolveFail;
				return ret;
			}

			(*L)(k,k) = sqrt((*L)(k,k));

			if( fabs((*L)(k,k)) < 1.0e-24 )
			{
				if( ierr )
				{
					if( *ierr == -1 ) std::cout << "Diagonal pivot is too small " << get_value((*L)(k,k)) << " row " << k << std::endl;
					*ierr = k+1;
				}
				else throw MatrixCholeskySolveFail;
				return ret;
			}

			for(enumerator i = k+1; i < n; ++i)
				(*L)(i,k) /= (*L)(k,k);

			for(enumerator j = k+1; j < n; ++j)
			{
				for(enumerator i = j; i < n; ++i)
					(*L)(i,j) -= (*L)(i,k)*(*L)(j,k);
			}
		}
		// LY=B
		Matrix<typename Promote<Var,typeB>::type> & Y = ret;
		for(enumerator i = 0; i < n; ++i)
		{
			for(enumerator k = 0; k < l; ++k)
			{
				for(enumerator j = 0; j < i; ++j)
					Y(i,k) -= Y(j,k)*(*L)(j,i);
				Y(i,k) /= (*L)(i,i);
			}
		}
		// L^TX = Y
		Matrix<typename Promote<Var,typeB>::type> & X = ret;
		for(enumerator it = n; it > 0; --it)
		{
			enumerator i = it-1;
			for(enumerator k = 0; k < l; ++k)
			{
				for(enumerator jt = n; jt > it; --jt)
				{
					enumerator j = jt-1;
					X(i,k) -= X(j,k)*(*L)(i,j);
				}
				X(i,k) /= (*L)(i,i);
			}
		}
		if( ierr ) *ierr = 0;
		return ret;
	}
#if defined(USE_AUTODIFF)
	template<>
	template<>
	__INLINE Matrix<Promote<variable,variable>::type>
	AbstractMatrixReadOnly<variable>::CholeskySolve(const AbstractMatrixReadOnly<variable> & B, int * ierr) const
	{
		const AbstractMatrixReadOnly<variable> & A = *this;
		assert(A.Rows() == A.Cols());
		assert(A.Rows() == B.Rows());
		enumerator n = A.Rows();
		enumerator l = B.Cols();
		Matrix<Promote<variable,variable>::type> ret(B);
		SymmetricMatrix<variable> L(A);
		INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
		A.GetMatrixInterval(beg, end, cnt);
		B.GetMatrixInterval(beg, end, cnt);
		Sparse::RowMerger* pmerger = NULL;
		if (cnt >= CNT_USE_MERGER)
		{
			merger->Resize(beg, end);
			pmerger = &(*merger);
		}
		//Outer product
		for(enumerator k = 0; k < n; ++k)
		{
			if( L(k,k).GetValue() < 0.0 )
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

			if( fabs(L(k,k).GetValue()) < 1.0e-24 )
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
				L(i,k) /= L(k,k);

			for(enumerator j = k+1; j < n; ++j)
			{
				for(enumerator i = j; i < n; ++i)
					L(i,j) -= L(i,k)*L(j,k);
			}
		}
		// LY=B
		Matrix<Promote<variable,variable>::type> & Y = ret;
		for(enumerator i = 0; i < n; ++i)
		{
			for(enumerator k = 0; k < l; ++k)
			{
				if( pmerger )
				{
					double value = Y(i,k).GetValue();
					pmerger->PushRow(1.0,Y(i,k).GetRow());
					for(enumerator j = 0; j < i; ++j)
					{
						value -= Y(j,k).GetValue()*L(j,i).GetValue();
						pmerger->AddRow(-Y(j,k).GetValue(),L(j,i).GetRow());
						pmerger->AddRow(-L(j,i).GetValue(),Y(j,k).GetRow());
					}
					Y(i,k).SetValue(value);
					pmerger->RetriveRow(Y(i,k).GetRow());
					pmerger->Clear();
				}
				else
				{
					for(enumerator j = 0; j < i; ++j)
						Y(i,k) -= Y(j,k)*L(j,i);
				}
				Y(i,k) /= L(i,i);
			}
		}
		// L^TX = Y
		Matrix<Promote<variable,variable>::type> & X = ret;
		for(enumerator it = n; it > 0; --it)
		{
			enumerator i = it-1;
			for(enumerator k = 0; k < l; ++k)
			{
				if( pmerger )
				{
					double value = X(i,k).GetValue();
					pmerger->PushRow(1.0,X(i,k).GetRow());
					for(enumerator jt = n; jt > it; --jt)
					{
						enumerator j = jt-1;
						value -= X(j,k).GetValue()*L(i,j).GetValue();
						pmerger->AddRow(-X(j,k).GetValue(),L(i,j).GetRow());
						pmerger->AddRow(-L(i,j).GetValue(),X(j,k).GetRow());
					}
					X(i,k).SetValue(value);
					pmerger->RetriveRow(X(i,k).GetRow());
					pmerger->Clear();
				}
				else
				{
					for(enumerator jt = n; jt > it; --jt)
					{
						enumerator j = jt-1;
						X(i,k) -= X(j,k)*L(i,j);
					}
				}
				X(i,k) /= L(i,i);
			}
		}
		if( ierr ) *ierr = 0;
		return ret;
	}


	template<>
	template<>
	__INLINE Matrix<Promote<INMOST_DATA_REAL_TYPE,variable>::type>
	AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE>::CholeskySolve(const AbstractMatrixReadOnly<variable> & B, int * ierr) const
	{
		const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE> & A = *this;
		assert(A.Rows() == A.Cols());
		assert(A.Rows() == B.Rows());
		enumerator n = A.Rows();
		enumerator l = B.Cols();
		Matrix<Promote<INMOST_DATA_REAL_TYPE,variable>::type> ret(B);
		SymmetricMatrix<INMOST_DATA_REAL_TYPE> L(A);
		INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
		B.GetMatrixInterval(beg, end, cnt);
		Sparse::RowMerger* pmerger = NULL;
		if (cnt >= CNT_USE_MERGER)
		{
			merger->Resize(beg, end);
			pmerger = &(*merger);
		}
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
				L(i,k) /= L(k,k);

			for(enumerator j = k+1; j < n; ++j)
			{
				for(enumerator i = j; i < n; ++i)
					L(i,j) -= L(i,k)*L(j,k);
			}
		}
		// LY=B
		Matrix<Promote<variable,variable>::type> & Y = ret;
		for(enumerator i = 0; i < n; ++i)
		{
			for(enumerator k = 0; k < l; ++k)
			{
				if( pmerger )
				{
					double value = Y(i,k).GetValue();
					pmerger->PushRow(1.0,Y(i,k).GetRow());
					for(enumerator j = 0; j < i; ++j)
					{
						value -= Y(j,k).GetValue()*L(j,i);
						pmerger->AddRow(-L(j,i),Y(j,k).GetRow());
					}
					Y(i,k).SetValue(value);
					pmerger->RetriveRow(Y(i,k).GetRow());
					pmerger->Clear();
				}
				else
				{
					for(enumerator j = 0; j < i; ++j)
						Y(i,k) -= Y(j,k)*L(j,i);
				}
				Y(i,k) /= L(i,i);
			}
		}
		// L^TX = Y
		Matrix<Promote<variable,variable>::type> & X = ret;
		for(enumerator it = n; it > 0; --it)
		{
			enumerator i = it-1;
			for(enumerator k = 0; k < l; ++k)
			{
				if( pmerger )
				{
					double value = X(i,k).GetValue();
					pmerger->PushRow(1.0,X(i,k).GetRow());
					for(enumerator jt = n; jt > it; --jt)
					{
						enumerator j = jt-1;
						value -= X(j,k).GetValue()*L(i,j);
						pmerger->AddRow(-L(i,j),X(j,k).GetRow());
					}
					X(i,k).SetValue(value);
					pmerger->RetriveRow(X(i,k).GetRow());
					pmerger->Clear();
				}
				else
				{
					for(enumerator jt = n; jt > it; --jt)
					{
						enumerator j = jt-1;
						X(i,k) -= X(j,k)*L(i,j);
					}
				}
				X(i,k) /= L(i,i);
			}
		}
		if( ierr ) *ierr = 0;
		return ret;
	}

	template<>
	template<>
	__INLINE Matrix<Promote<variable,INMOST_DATA_REAL_TYPE>::type>
	AbstractMatrixReadOnly<variable>::CholeskySolve(const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE> & B, int * ierr) const
	{
		const AbstractMatrixReadOnly<variable> & A = *this;
		assert(A.Rows() == A.Cols());
		assert(A.Rows() == B.Rows());
		enumerator n = A.Rows();
		enumerator l = B.Cols();
		Matrix<Promote<variable,INMOST_DATA_REAL_TYPE>::type> ret(B);
		SymmetricMatrix<variable> L(A);
		INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end = 0, cnt = 0;
		A.GetMatrixInterval(beg, end, cnt);
		Sparse::RowMerger* pmerger = NULL;
		if (cnt >= CNT_USE_MERGER)
		{
			merger->Resize(beg, end);
			pmerger = &(*merger);
		}
		//Outer product
		for(enumerator k = 0; k < n; ++k)
		{
			if( L(k,k).GetValue() < 0.0 )
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

			if( fabs(L(k,k).GetValue()) < 1.0e-24 )
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
				L(i,k) /= L(k,k);

			for(enumerator j = k+1; j < n; ++j)
			{
				for(enumerator i = j; i < n; ++i)
					L(i,j) -= L(i,k)*L(j,k);
			}
		}
		// LY=B
		Matrix<Promote<variable,INMOST_DATA_REAL_TYPE>::type> & Y = ret;
		for(enumerator i = 0; i < n; ++i)
		{
			for(enumerator k = 0; k < l; ++k)
			{
				if( pmerger )
				{
					double value = Y(i,k).GetValue();
					pmerger->PushRow(1.0,Y(i,k).GetRow());
					for(enumerator j = 0; j < i; ++j)
					{
						value -= Y(j,k).GetValue()*L(j,i).GetValue();
						pmerger->AddRow(-Y(j,k).GetValue(),L(j,i).GetRow());
						pmerger->AddRow(-L(j,i).GetValue(),Y(j,k).GetRow());
					}
					Y(i,k).SetValue(value);
					pmerger->RetriveRow(Y(i,k).GetRow());
					pmerger->Clear();
				}
				else
				{
					for(enumerator j = 0; j < i; ++j)
						Y(i,k) -= Y(j,k)*L(j,i);
				}
				Y(i,k) /= L(i,i);
			}
		}
		// L^TX = Y
		Matrix<Promote<variable,INMOST_DATA_REAL_TYPE>::type> & X = ret;
		for(enumerator it = n; it > 0; --it)
		{
			enumerator i = it-1;
			for(enumerator k = 0; k < l; ++k)
			{
				if( pmerger )
				{
					double value = X(i,k).GetValue();
					pmerger->PushRow(1.0,X(i,k).GetRow());
					for(enumerator jt = n; jt > it; --jt)
					{
						enumerator j = jt-1;
						value -= X(j,k).GetValue()*L(i,j).GetValue();
						pmerger->AddRow(-X(j,k).GetValue(),L(i,j).GetRow());
						pmerger->AddRow(-L(i,j).GetValue(),X(j,k).GetRow());
					}
					X(i,k).SetValue(value);
					pmerger->RetriveRow(X(i,k).GetRow());
					pmerger->Clear();
				}
				else
				{
					for(enumerator jt = n; jt > it; --jt)
					{
						enumerator j = jt-1;
						X(i,k) -= X(j,k)*L(i,j);
					}
				}
				X(i,k) /= L(i,i);
			}
		}
		if( ierr ) *ierr = 0;
		return ret;
	}
#endif //USE_AUTODIFF
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var>::Solve(const AbstractMatrixReadOnly<typeB> & B, int * ierr) const
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
		Matrix<typename Promote<Var,typeB>::type> ret(m,l);
		ConstMatrixTranspose<Var> At = this->Transpose(); //m by n matrix
		Matrix<typename Promote<Var,typeB>::type> AtB = At*B; //m by l matrix
		Matrix<Var> AtA = At*(*this); //m by m matrix
		assert(l == AtB.Cols());
		//enumerator l = AtB.Cols();
        //enumerator n = Rows();

		std::vector<enumerator> order(m);

		Var temp;
		INMOST_DATA_REAL_TYPE max,v;
		typename Promote<Var,typeB>::type tempb;
		for(enumerator i = 0; i < m; ++i) order[i] = i;
		for(enumerator i = 0; i < m; i++)
		{
			enumerator maxk = i, maxq = i, temp2;
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
						//std::swap(AtA(maxk,q),AtA(i,q));
						temp = AtA(maxk,q);
						AtA(maxk,q) = AtA(i,q);
						AtA(i,q) = temp;
					}
					//exchange rhs
					for(enumerator q = 0; q < l; q++) // over columns of B
					{
						//std::swap(AtB(maxk,q),AtB(i,q));
						tempb = AtB(maxk,q);
						AtB(maxk,q) = AtB(i,q);
						AtB(i,q) = tempb;
					}
				}
				//Exchange columns
				if( maxq != i )
				{
					for(enumerator k = 0; k < m; k++) //over rows
					{
						//std::swap(AtA(k,maxq),AtA(k,i));
						temp = AtA(k,maxq);
						AtA(k,maxq) = AtA(k,i);
						AtA(k,i) = temp;
					}
					//remember order in sol
					{
						//std::swap(order[maxq],order[i]);
						temp2 = order[maxq];
						order[maxq] = order[i];
						order[i] = temp2;
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
				if( ok ) AtA(i,i) = real_part(AtA(i,i)) < 0.0 ? - 1.0e-12 : 1.0e-12;
				else
				{
					if( ierr )
					{
						if( *ierr == -1 )
						{
							std::cout << "Failed to invert matrix diag " << get_value(AtA(i,i)) << std::endl;
							std::cout << "rhs:";
							for(enumerator k = 0; k < l; k++)
								std::cout << " " << get_value(AtB(i,k));
							std::cout << std::endl;
						}
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
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::ExtractSubMatrix(enumerator ibeg, enumerator iend, enumerator jbeg, enumerator jend) const
	{
        // TODO по-моему, здесь ошибка.
		//assert(ibeg < Rows());
		//assert(iend < Rows());
		//assert(jbeg < Cols());
		//assert(jend < Cols());

        assert(ibeg < Rows());
		assert(iend <= Rows());
		assert(jbeg < Cols());
		assert(jend <= Cols());

		Matrix<Var> ret(iend-ibeg,jend-jbeg);
		for(enumerator i = ibeg; i < iend; ++i)
		{
			for(enumerator j = jbeg; j < jend; ++j)
				ret(i-ibeg,j-jbeg) = (*this)(i,j);
		}
		return ret;
	}

	/*
	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::Repack(enumerator rows, enumerator cols) const
	{
		assert(Cols()*Rows()==rows*cols);
		Matrix<Var> ret(*this);
		ret.Rows() = rows;
		ret.Cols() = cols;
		return ret;
	}
	*/

	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::PseudoInvert(INMOST_DATA_REAL_TYPE tol, int * ierr) const
	{
		Matrix<Var> ret(Cols(),Rows());
		Matrix<Var> U,S,V;
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
			if (get_value(fabs(S(k, k))) > tol)
				S(k, k) = 1.0 / S(k, k);
			else
				S(k, k) = 0.0;
		}
		//ret = V*S*U.Transpose();
		ret = V.ConjugateTranspose().Transpose() * S * U.ConjugateTranspose();
		if( ierr ) *ierr = 0;
		return ret;
	}
	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::Root(INMOST_DATA_ENUM_TYPE iter, INMOST_DATA_REAL_TYPE tol, int *ierr) const
	{
		assert(Rows() == Cols());
		Matrix<Var> ret(Cols(),Cols());
		Matrix<Var> ret0(Cols(),Cols());
		ret.Zero();
		ret0.Zero();
		int k = 0;
		for(k = 0; k < Cols(); ++k) ret(k,k) = ret0(k,k) = 1;
		while(k < iter)
		{
			ret0 = ret;
			ret = 0.5*(ret + (*this)*ret.Invert());
			if( (ret - ret0).FrobeniusNorm() < tol ) return ret;
		}
		if( ierr )
		{
			if( *ierr == -1 ) std::cout << "Failed to find square root of matrix by Babylonian method" << std::endl;
			*ierr = 1;
			return ret;
		}
		return ret;
	}
	template<typename Var>
	Matrix<Var>
	AbstractMatrixReadOnly<Var>::Power(INMOST_DATA_REAL_TYPE n, int * ierr) const
	{
		Matrix<Var> ret(Cols(),Rows());
		Matrix<Var> U, S, V;
		bool success = SVD(U,S,V);
		if( !success )
		{
			if( ierr )
			{
				if( *ierr == -1 )
					std::cout << "Failed to compute singular value decomposition of the matrix" << std::endl;
				*ierr = 1;
				return ret;
			}
			else throw MatrixPseudoSolveFail;
		}
        for(INMOST_DATA_ENUM_TYPE k = 0; k < S.Cols(); ++k)
			if( get_value(S(k,k)) ) S(k,k) = pow(S(k,k),n);
        if( n >= 0 )
			ret = U*S*V.Transpose();
		else
			ret = V*S*U.Transpose();
		if( ierr ) *ierr = 0;
		return ret;
	}
	template<typename Var>
	template<typename typeB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var>::PseudoSolve(const AbstractMatrixReadOnly<typeB> & B, INMOST_DATA_REAL_TYPE tol, int * ierr) const
	{
		Matrix<typename Promote<Var,typeB>::type> ret(Cols(),B.Cols());
		Matrix<Var> U,S,V;
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
		for(int k = 0; k < (int)S.Cols(); ++k)
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
	/*
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
	*/
	template<typename Var>
	SubMatrix<Var> AbstractMatrix<Var>::operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col)
	{
		return SubMatrix<Var>(*this,first_row,last_row,first_col,last_col);
	}
	template<typename Var>
	ConstSubMatrix<Var> AbstractMatrixReadOnly<Var>::operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const
	{
		return ConstSubMatrix<Var>(*this, first_row, last_row, first_col, last_col);
	}

	template<typename Var>
	BlockOfMatrix<Var> AbstractMatrix<Var>::BlockOf(enumerator nrows, enumerator ncols, enumerator offset_row, enumerator offset_col)
	{
		return BlockOfMatrix<Var>(*this,nrows,ncols,offset_row,offset_col);
	}
	template<typename Var>
	ConstBlockOfMatrix<Var> AbstractMatrixReadOnly<Var>::BlockOf(enumerator nrows, enumerator ncols, enumerator offset_row, enumerator offset_col) const
	{
		return ConstBlockOfMatrix<Var>(*this,nrows,ncols,offset_row,offset_col);
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
			size_t j;
			while(2*k <= size)
			{
				j = 2*static_cast<size_t>(k);
				if(j < size && keys[heap[j]] > keys[heap[j+1]])
					j++;
				if(keys[heap[k]] <= keys[heap[j]])
					break;
				swap(k, static_cast<INMOST_DATA_ENUM_TYPE>(j));
				k = static_cast<INMOST_DATA_ENUM_TYPE>(j);
			}
		}
	public:
		BinaryHeapDense(INMOST_DATA_REAL_TYPE * pkeys, INMOST_DATA_ENUM_TYPE len)
		{
			size_max = len;
			keys = pkeys;
			size = 0;
			heap.resize(static_cast<size_t>(size_max)+1);
			index.resize(static_cast<size_t>(size_max)+1);
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
			heap[static_cast<size_t>(size)+1] = ENUMUNDEF;
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
	void AbstractMatrixReadOnly<Var>::MPT(INMOST_DATA_ENUM_TYPE * Perm, INMOST_DATA_REAL_TYPE * SL, INMOST_DATA_REAL_TYPE * SR) const
	{
		const INMOST_DATA_ENUM_TYPE EOL = ENUMUNDEF-1;
		int n = Rows();
		int m = Cols();
		INMOST_DATA_REAL_TYPE u, l;
		std::vector<INMOST_DATA_REAL_TYPE> Cmax(m,0.0);
		std::vector<INMOST_DATA_REAL_TYPE> U(m,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		std::vector<INMOST_DATA_REAL_TYPE> V(n,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		std::fill(Perm,Perm+m,ENUMUNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> IPerm(std::max(n,m),ENUMUNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> ColumnList(m,ENUMUNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> ColumnPosition(n,ENUMUNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> Parent(n,ENUMUNDEF);
		std::vector<INMOST_DATA_REAL_TYPE> Dist(m,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
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

	template<typename Var>
	bool AbstractMatrixReadOnly<Var>::cSVD(AbstractMatrix<Var>& U, AbstractMatrix<Var>& Sigma, AbstractMatrix<Var>& V) const
	{
		Var cs, eta, f, g, h, q, r, sn, w, x, y, z;
		int i, j, k, l, m, n, p = 0;
		std::vector<Var> t, b, c, s;
		Matrix<Var> A = *this;
		m = Rows();
		n = Cols();
		//data eta / 1.1920929e-07 /
		//data tol / 1.5e-31 /
		if (m > n)
		{
			bool success = ConjugateTranspose().cSVD(V, Sigma, U);
			if (success)
				return true;
			else return false;
		}
		//  Householder reduction.
		c.resize(n);
		t.resize(n);
		s.resize(n);
		b.resize(n);
		c[0] = 0.0;
		k = 0;

		while (true)
		{
			// Elimination of A(I, K), I = K + 1, ..., M.
			z = 0.0;
			for (i = k; i < m; ++i)
				z += A(i, k) * conj(A(i, k));
			b[k] = 0.0;

			if (fabs(get_value(z)))
			{
				z = sqrt(z);
				b[k] = z;
				w = fabs(A(k, k));
				q = 1.0;
				if (fabs(get_value(w))) q = A(k, k) / w;
				A(k, k) = q * (z + w);

				if (k != n + p - 1)
				{
					for (j = k + 1; j < n + p; ++j)
					{
						q = 0.0;
						for (i = k; i < m; ++i)
							q += conj(A(i, k)) * A(i, j);
						q = q / (z * (z + w));
						for (i = k; i < m; ++i)
							A(i, j) -= q * A(i, k);
					}
					//  Phase transformation.
					q = -conj(A(k, k)) / fabs(A(k, k));
					for (j = k + 1; j < n + p; ++j)
						A(k, j) *= q;
				}
			}
			//  Elimination of A(K, J), J = K + 2, ..., N
			if (k == n - 1) break;
			z = 0.0;
			for (j = k + 1; j < n; ++j)
				z += conj(A(k, j)) * A(k, j);
			c[k + 1] = 0.0;

			if (fabs(get_value(z)))
			{
				z = sqrt(z);
				c[k + 1] = z;
				w = fabs(A(k, k + 1));
				q = 1.0;
				if (fabs(get_value(w))) q = A(k, k + 1) / w;
				A(k, k + 1) = q * (z + w);
				for (i = k + 1; i < m; ++i)
				{
					q = 0.0;
					for (j = k + 1; j < n; ++j)
						q += conj(A(k, j)) * A(i, j);
					q = q / (z * (z + w));

					for (j = k + 1; j < n; ++j)
						A(i, j) -= q * A(k, j);
				}
				//  Phase transformation.
				q = -conj(A(k, k + 1)) / fabs(A(k, k + 1));
				for (i = k + 1; i < m; ++i)
					A(i, k + 1) = A(i, k + 1) * q;
			}
			k++;
		}
		for (k = 0; k < n; k++)
		{
			s[k] = b[k];
			t[k] = c[k];
		}
		//  Initialization of U and V.
		U.Resize(m,m);
		U.Zero();
		for(j = 0; j < m; ++j) U(j,j) = 1.0;
		V.Resize(n,n);
		V.Zero();
		for(j = 0; j < n; ++j) V(j,j) = 1.0;
		//  QR diagonalization.
		for (k = n - 1; k >= 0; k--)
		{
			// Test for split.
			while (true)
			{
				bool skip = false;
				for (l = k; l >= 0; l--)
				{
					if (!fabs(get_value(t[l])))
					{
						skip = true;
						break;
					}
					if (!fabs(get_value(s[l - 1]))) break;
				}
				//  Cancellation of E(L).
				if (!skip)
				{
					cs = 0.0;
					sn = 1.0;
					for (i = l; i <= k; ++i)
					{
						f = sn * t[i];
						t[i] = cs * t[i];
						if (!fabs(get_value(f))) break;
						h = s[i];
						w = sqrt(f * f + h * h);
						s[i] = w;
						cs = h / w;
						sn = -f / w;
						for (j = 0; j < n; ++j)
						{
							x = real_part(U(j, l - 1));
							y = real_part(U(j, i));
							U(j, l - 1) = x * cs + y * sn;
							U(j, i) = y * cs - x * sn;
						}

						if (p == 0) continue;

						for (j = n; j < n + p; ++j)
						{
							q = A(l - 1, j);
							r = A(i, j);
							A(l - 1, j) = q * cs + r * sn;
							A(i, j) = r * cs - q * sn;
						}
					}
				}
				//  Test for convergence.
				w = s[k];
				if (l == k) break;
				//  Origin shift.
				x = s[l];
				y = s[k - 1];
				g = t[k - 1];
				h = t[k];
				f = ((y - w) * (y + w) + (g - h) * (g + h)) / (2.0 * h * y);
				g = sqrt(f * f + 1.0);
				if (real_part(get_value(f)) < 0.0) g = -g;
				f = ((x - w) * (x + w) + (y / (f + g) - h) * h) / x;
				//  QR step.
				cs = 1.0;
				sn = 1.0;
				for (i = l+1; i <= k; ++i)
				{
					g = t[i];
					y = s[i];
					h = sn * g;
					g = cs * g;
					w = sqrt(h * h + f * f);
					t[i - 1] = w;
					cs = f / w;
					sn = h / w;
					f = x * cs + g * sn;
					g = g * cs - x * sn;
					h = y * sn;
					y = y * cs;
					for (j = 0; j < n; ++j)
					{
						x = real_part(V(j, i - 1));
						w = real_part(V(j, i));
						V(j, i - 1) = x * cs + w * sn;
						V(j, i) = w * cs - x * sn;
					}
					w = sqrt(h * h + f * f);
					s[i - 1] = w;
					cs = f / w;
					sn = h / w;
					f = cs * g + sn * y;
					x = cs * y - sn * g;
					for (j = 0; j < n; ++j)
					{
						y = real_part(U(j, i - 1));
						w = real_part(U(j, i));
						U(j, i - 1) = y * cs + w * sn;
						U(j, i) = w * cs - y * sn;
					}
					if (p == 0) break;
					for (j = n; j < n + p; ++j)
					{
						q = A(i - 1, j);
						r = A(i, j);
						A(i - 1, j) = q * cs + r * sn;
						A(i, j) = r * cs - q * sn;
					}
				}
				t[l] = 0.0;
				t[k] = f;
				s[k] = x;
			}
			//  Convergence.
			if (real_part(get_value(w)) >= 0.0) continue;
			s[k] = -w;
			for (j = 0; j < n; ++j) V(j, k) = -V(j, k);
		}
		// Sort the singular values.
		for (k = 0; k < n; ++k)
		{
			g = -1.0;
			j = k;
			for (i = k; i < n; ++i)
			{
				if (real_part(get_value(g)) < real_part(get_value(s[i])))
				{
					g = s[i];
					j = i;
				}
			}
			if (j == k) continue;

			s[j] = s[k];
			s[k] = g;
			//  Interchange V(1:N, J) and V(1:N, K).
			for (i = 0; i < n; ++i)
			{
				q = V(i, j);
				V(i, j) = V(i, k);
				V(i, k) = q;
			}
			//  Interchange U(1:N, J) and U(1:N, K).
			for (i = 0; i < n; ++i)
			{
				q = U(i, j);
				U(i, j) = U(i, k);
				U(i, k) = q;
			}
			//  Interchange A(J, N1:NP) and A(K, N1:NP).
			if (p == 0) continue;
			for (i = n; i < n + p; ++i)
			{
				q = A(j, i);
				A(j, i) = A(k, i);
				A(k, i) = q;
			}
		}
		//  Back transformation.
		for (k = n - 1; k >= 0; k--)
		{
			if (b[k] == 0.0) continue;
			q = -A(k, k) / fabs(A(k, k));
			for (j = 0; j < m; ++j)
				U(k, j) *= q;
			for (j = 0; j < m; ++j)
			{
				q = 0.0;
				for (i = k; i < m; ++i)
					q += conj(A(i, k)) * U(i, j);
				q = q / (fabs(A(k, k)) * b[k]);
				for (i = k; i < m; ++i)
					U(i, j) -= q * A(i, k);
			}
		}
		if (n > 1)
		{
			for (k = n - 2; k >= 0; k--)
			{
				if (c[k+1] == 0.0) continue;
				q = -conj(A(k, k+1)) / fabs(A(k, k+1));
				for (j = 0; j < n; ++j)
					V(k+1, j) *= q;

				for (j = 0; j < n; ++j)
				{
					q = 0.0;
					for (i = k+1; i < n; ++i)
						q += A(k, i) * V(i, j);
					q = q / (fabs(A(k, k+1)) * c[k+1]);
					for (i = k+1; i < n; ++i)
						V(i, j) -= q * conj(A(k, i));
				}
			}
		}
		Sigma.Resize(m, n);
		Sigma.Zero();
		for (i = 0; i < n; ++i)
			Sigma(i, i) = s[i];
		return true;
	}

	template<typename Var>
	bool AbstractMatrixReadOnly<Var>::SVD(AbstractMatrix<Var>& U, AbstractMatrix<Var>& Sigma, AbstractMatrix<Var>& V, bool order_singular_values, bool nonnegative) const
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
			Sigma.Resize(m, m);
			Sigma.Zero();
			V.Resize(m, m);
		}
		else // n < m
		{
			U.Resize(n, n);
			Sigma.Resize(n, n);
			Sigma.Zero();
			V.Resize(m, n);
		}
		if (n < m)
		{
			bool success = Transpose().SVD(V, Sigma, U);
			if (success)
			{
				//U.Swap(V);
				//U = U.Transpose();
				//V = V.Transpose();
				return true;
			}
			else return false;
		} //m <= n
		std::vector<Var> rv1(m);
		std::swap(n, m); //this how original algorithm takes it
		// Householder reduction to bidiagonal form
		for (i = 0; i < n; i++)
		{
			// left-hand reduction
			l = i + 1;
			rv1[i] = scale * g;
			g = s = scale = 0.0;
			if (i < m)
			{
				for (k = i; k < m; k++) scale += fabs(U(k, i));
				if (fabs(get_value(scale)))
				{
					for (k = i; k < m; k++)
					{
						U(k, i) /= scale;
						//s += U(k, i) * U(k, i); //original
						s += U(k, i) * conj(U(k, i)); //for complex
					}
					f = U(i, i);
					g = -sign_func(sqrt(s), f);
					h = f * g - s;
					U(i, i) = f - g;
					if (i != n - 1 && fabs(get_value(h)))
					{
						for (j = l; j < n; j++)
						{
							//for (s = 0.0, k = i; k < m; k++) s += U(k, i) * U(k, j);
							for (s = 0.0, k = i; k < m; k++) s += conj(U(k, i)) * U(k, j);
							f = s / h;
							for (k = i; k < m; k++) U(k, j) += f * U(k, i);
						}
					}
					for (k = i; k < m; k++) U(k, i) *= scale;
				}
			}
			Sigma(i, i) = scale * g;
			// right-hand reduction
			g = s = scale = 0.0;
			if (i < m && i != n - 1)
			{
				for (k = l; k < n; k++) scale += fabs(U(i, k));
				if (fabs(get_value(scale)))
				{
					for (k = l; k < n; k++)
					{
						U(i, k) = U(i, k) / scale;
						//s += U(i, k) * U(i, k); //original
						s += U(i, k) * conj(U(i, k));
					}
					f = U(i, l);
					g = -sign_func(sqrt(s), f);
					h = f * g - s;
					U(i, l) = f - g;
					for (k = l; k < n; k++) rv1[k] = U(i, k) / h;
					if (i != m - 1)
					{
						for (j = l; j < m; j++)
						{
							for (s = 0.0, k = l; k < n; k++) s += conj(U(i, k)) * U(j, k);
							for (k = l; k < n; k++) U(j, k) += s * rv1[k];
						}
					}
					for (k = l; k < n; k++) U(i, k) *= scale;
				}
			}
			anorm = max_func(anorm, fabs(get_value(Sigma(i, i))) + fabs(get_value(rv1[i])));
		}

		// accumulate the right-hand transformation
		for (i = n - 1; i >= 0; i--)
		{
			if (i < (n - 1))
			{
				if (fabs(get_value(g)))
				{
					for (j = l; j < n; j++) V(j, i) = ((U(i, j) / U(i, l)) / g);
					// double division to avoid underflow
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = l; k < n; k++) s += U(i, k) * V(k, j);
						for (k = l; k < n; k++) V(k, j) += s * V(k, i);
					}
				}
				for (j = l; j < n; j++) V(i, j) = V(j, i) = 0.0;
			}
			V(i, i) = 1.0;
			g = rv1[i];
			l = i;
		}

		// accumulate the left-hand transformation
		for (i = n - 1; i >= 0; i--)
		{
			l = i + 1;
			g = Sigma(i, i);
			if (i < (n - 1))
				for (j = l; j < n; j++)
					U(i, j) = 0.0;
			if (fabs(get_value(g)))
			{
				g = 1.0 / g;
				if (i != n - 1)
				{
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = l; k < m; k++) s += conj(U(k, i)) * U(k, j);
						f = (s / U(i, i)) * g;
						for (k = i; k < m; k++) U(k, j) += f * U(k, i);
					}
				}
				for (j = i; j < m; j++) U(j, i) *= g;
			}
			else for (j = i; j < m; j++) U(j, i) = 0.0;
			U(i, i) += 1;
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
					if (l == 0) break;
					if (fabs(get_value(Sigma(nm, nm))) + anorm == anorm)
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
							g = Sigma(i, i);
							h = pythag(f, g);
							Sigma(i, i) = h;
							h = 1.0 / h;
							c = g * h;
							s = (-f * h);
							for (j = 0; j < m; j++)
							{
								y = U(j, nm);
								z = U(j, i);
								U(j, nm) = (y * c + z * s);
								U(j, i) = (z * c - y * s);
							}
						}
					}
				}
				z = Sigma(k, k);
				if (l == k)
				{// convergence
					if (real_part(get_value(z)) < 0.0 && nonnegative)
					{// make singular value nonnegative
						Sigma(k, k) = -z;
						for (j = 0; j < n; j++) V(j, k) = -V(j, k);
					}
					break;
				}
				if (its >= 30)
				{
					std::cout << "No convergence after " << its << " iterations" << std::endl;
					std::swap(n, m);
					return false;
				}
				// shift from bottom 2 x 2 minor
				x = Sigma(l, l);
				nm = k - 1;
				y = Sigma(nm, nm);
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
					y = Sigma(i, i);
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
						x = V(jj, j);
						z = V(jj, i);
						V(jj, j) = (x * c + z * s);
						V(jj, i) = (z * c - x * s);
					}
					z = pythag(f, h);
					Sigma(j, j) = z;
					if (fabs(get_value(z)))
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = (c * g) + (s * y);
					x = (c * y) - (s * g);
					for (jj = 0; jj < m; jj++)
					{
						y = U(jj, j);
						z = U(jj, i);
						U(jj, j) = (y * c + z * s);
						U(jj, i) = (z * c - y * s);
					}
				}
				rv1[l] = 0.0;
				rv1[k] = f;
				Sigma(k, k) = x;
			}
		}
		//CHECK THIS!
		if (order_singular_values)
		{
			Var temp;
			for (i = 0; i < n; i++)
			{
				k = i;
				for (j = i + 1; j < n; ++j)
					if (real_part(get_value(Sigma(k, k))) < real_part(get_value(Sigma(j, j)))) k = j;
				if (real_part(get_value(Sigma(k, k))) > real_part(get_value(Sigma(i, i))))
				{
					temp = Sigma(k, k);
					Sigma(k, k) = Sigma(i, i);
					Sigma(i, i) = temp;
					// U is m by n
					for (int j = 0; j < m; ++j)
					{
						temp = U(j, k);
						U(j, k) = U(j, i);
						U(j, i) = temp;
					}
					// V is n by n
					for (int j = 0; j < n; ++j)
					{
						temp = V(j, k);
						V(j, k) = V(j, i);
						V(j, i) = temp;
					}
				}
			}
		}

		std::swap(n, m);
		return true;
	}


	template<typename Var>
	bool AbstractMatrixReadOnly<Var>::EigenDecomposition(AbstractMatrix<INMOST_DATA_CPLX_TYPE> & V, AbstractMatrix<INMOST_DATA_CPLX_TYPE> & Lambda, bool order_eigen_values,
            double abs_eps, double rel_eps, int max_iters) const
    {
        // Вспомогательная структура.
        struct FrancisQR
        {
            // Получение собственных чисел матрицы 2 на 2.
            static std::pair<std::complex<double>, std::complex<double>> get_2x2_eigenvalues(double a, double b, double c, double d)
            {
                double det = a * d - c * b;
                double trace_2 = (a + d) / 2.0;

                // Два случая: либо корни действительные, либо комплексные
                // (и, как следствие, они сопряжены).
                double gap = trace_2 * trace_2 - det;
                if (gap >= 0.0)
                {
                    gap = std::sqrt(gap);
                    return std::make_pair(trace_2 + gap, trace_2 - gap);
                }
                else
                {
                    gap = std::sqrt(-gap);
                    std::complex<double> eigenvalue = {trace_2, gap};
                    return std::make_pair(eigenvalue, std::conj(eigenvalue));
                }
            }

            static std::pair<size_t, size_t> deflate_hessenberg(const AbstractMatrix<Var>& matrix, std::vector<std::complex<double>>& eigenvalues,
                    double abs_eps, double rel_eps)
            {
                size_t dim = matrix.Rows();
                size_t lower = 0;
                size_t upper = dim;

                double eps = abs_eps + matrix.FrobeniusNorm() * rel_eps;

                // Верхняя граница.
                while (upper > 2)
                {
                    if (std::abs(matrix(upper-1, upper-2)) < eps)
                    {
                        eigenvalues.push_back(matrix(upper-1, upper-1));
                        --upper;
                    }
                    else if (std::abs(matrix(upper-2, upper-3)) < eps)
                    {
                        auto evs = get_2x2_eigenvalues(matrix(upper - 2, upper - 2),
                                matrix(upper - 2, upper - 1),
                                matrix(upper - 1, upper - 2),
                                matrix(upper - 1, upper - 1));

                        eigenvalues.push_back(evs.first);
                        eigenvalues.push_back(evs.second);

                        upper -= 2;
                    }
                    else break;
                }

                // Нижняя граница.
                while (lower < upper)
                {
                    if (std::abs(matrix(lower+1, lower)) < eps)
                    {
                        eigenvalues.push_back(matrix(lower, lower));
                        ++lower;
                    }
                    else if (std::abs(matrix(lower+2, lower+1)) < eps)
                    {
                        auto evs = get_2x2_eigenvalues(matrix(lower, lower),
                                matrix(lower, lower + 1),
                                matrix(lower + 1, lower),
                                matrix(lower + 1, lower + 1));

                        eigenvalues.push_back(evs.first);
                        eigenvalues.push_back(evs.second);

                        lower += 2;
                    }
                    else break;
                }

                return std::make_pair(lower, upper);
            }
        };

        // Проверка квадратности матрицы.
        size_t n_rows    = Rows();
        size_t n_columns = Cols();
        if (n_rows != n_columns)
        {
            throw "Matrix has to be square.";
        }

        // Подготовка матриц для вывода.
        V.Resize(n_rows, n_rows);
        V.Zero();
        Lambda.Resize(n_rows, n_rows);
        Lambda.Zero();

        // Массив для собственных значений.
        std::vector<std::complex<double>> eigenvalues;

        // ПОИСК СОБСТВЕННЫХ ЗНАЧЕНИЙ.
        {
            size_t dim = Rows();

            // Приводим матрицу к верхнехессенберговой форме.
            //auto pair_Q_H = FrancisQR::to_hessenberg(Matrix<Var>(*this));
            //Matrix<Var> H = pair_Q_H.second;
            Matrix<Var> Q;
            Matrix<Var> H;
            bool success = Hessenberg(Q, H);

            // Индексы, выделяющие текущую подматрицу.
            size_t lower = 0;
            size_t upper = dim;

            int step = 0;
            for (; step < max_iters; ++step)
            {
                //H.Print();

                // Дефляция.
                auto pair_lower_upper = FrancisQR::deflate_hessenberg(H, eigenvalues, abs_eps, rel_eps);
                lower = pair_lower_upper.first;
                upper = pair_lower_upper.second;

                // Алгоритм сошёлся, если подматрица пустая.
                if (lower == upper) { break; }

                // Выделяем подматрицу.
                H = H.ExtractSubMatrix(lower, upper, lower, upper);
                dim = upper - lower;

                // Единичный вектор вдоль первой оси (в текущем подпространстве).
                Matrix<Var> e_1(dim,1);
                e_1.Zero();
                e_1(0,0) = 1.0;

                // Вычисление сдвигов по правой нижней подматрице.
                double trace = H(dim-2, dim-2) + H(dim-1, dim-1);
                double det   = H(dim-2, dim-2) * H(dim-1, dim-1) - H(dim-2, dim-1) * H(dim-1, dim-2);

                // Вычисление первого столбца сдвинутой матрицы.
                Matrix<Var> vector = H * (H * e_1 - trace * e_1) + det * e_1;

                // Возмущение рефлектором.
                Matrix<Var> reflector = vector.HouseholderReflector();
                Matrix<Var> outer_product = Matrix<Var>(reflector * reflector.Transpose());
                H = H - 2.0 * H * outer_product;
                H = H - 2.0 * outer_product * H;

                // Возвращение к хессенберговой форме.
                //pair_Q_H = FrancisQR::to_hessenberg(H);
                //H = pair_Q_H.second;
                Matrix<Var> _H;
                bool success = H.Hessenberg(Q, _H);
                H = _H;
            }

            if (step == max_iters)
            {
                return false;
            }
        }

        // Поиск близких собственных значений.
        std::vector<std::complex<Var>> processed_eigenvalues;
        processed_eigenvalues.reserve(n_rows);
        while (processed_eigenvalues.size() < n_rows)
        {
            size_t best_index = 0;
            std::vector<size_t> best_indices;
            for (size_t index = 0; index < eigenvalues.size(); ++index)
            {
                // Радиус окружности.
                double radius = abs_eps + rel_eps * std::abs(eigenvalues[index]);

                std::vector<size_t> indices;
                for (size_t jndex = 0; jndex < eigenvalues.size(); ++jndex)
                {
                    // Достаточно ли близко два выбранных собственных числа.
                    double delta = std::abs(eigenvalues[index] - eigenvalues[jndex]);
                    if (delta < radius)
                    { indices.push_back(jndex); }
                }

                if (indices.size() > best_indices.size())
                {
                    best_index = index;
                    best_indices = indices;
                }
            }

            // Если все собственные значения достаточно хорошо различимы, ничего делать не надо.
            if (best_indices.size() <= 1)
            {
                processed_eigenvalues.insert(processed_eigenvalues.begin(), eigenvalues.begin(), eigenvalues.end());
            }
            else
            {
                // Иначе заменяем часть собственных значений средним.
                std::complex<Var> average = 0.0;
                for (ssize_t subindex = best_indices.size() - 1; subindex >= 0; --subindex) // Здесь используется тот факт, что значения в best_indices идут по возрастанию.
                {
                    size_t index = best_indices[subindex];
                    average += eigenvalues[index];

                    // Удалям обработанные собственные значения.
                    eigenvalues[index] = eigenvalues.back();
                    eigenvalues.pop_back();
                }
                average /= best_indices.size();

                processed_eigenvalues.insert(processed_eigenvalues.end(), best_indices.size(), average);
            }
        }
        eigenvalues = processed_eigenvalues;

        // Сортировка собственных значений по действительной части, если требуется.
        if (order_eigen_values)
        {
            struct ComplexComparator
            {
                static bool compare(const std::complex<Var>& left, const std::complex<Var>& right)
                { return std::real(left) == std::real(right) ? std::imag(left) > std::imag(right) : std::real(left) > std::real(right); };
            };

            std::sort(eigenvalues.begin(), eigenvalues.end(), ComplexComparator::compare);
        }

        // Запись собственных значений в диагональную матрицу.
        for (size_t index = 0; index < n_rows; ++index)
        { Lambda(index,index) = eigenvalues[index]; }

        // Поиск собственных векторов.
        {
            // Единичная матрица.
            Matrix<std::complex<Var>> I = Matrix<std::complex<Var>>::Unit(n_rows);

            // Итерация по найденным собственным значениям.
            for (size_t index = 0; index < eigenvalues.size();) //; ++index)
            {
                // Поиск одинаковых СЗ.
                size_t jndex = index + (eigenvalues[index].imag() > 0 ? 2 : 1);
                for (; jndex < eigenvalues.size(); ++jndex)
                {
                    if (eigenvalues[index] != eigenvalues[jndex])
                    { break; }
                }

                // Сдвинутая матрица.
                auto eigenvalue = eigenvalues[index];
                Matrix<std::complex<Var>> A = Matrix<std::complex<Var>>(*this) - I * eigenvalue;

                // Используем SVD-разложение для поиска ядра.
                // Собственный вектор - последний столбец V
                // (соответствующий наименьшему сингулярному значению).
                // (в случае большой кратности СЗ берется всё ядро).
                Matrix<INMOST_DATA_CPLX_TYPE> _U, _Sigma, _V;
                bool success = A.SVD(_U, _Sigma, _V);
                if (!success) { return false; }

                V(0, n_rows, index, jndex) = _V(0, n_rows, n_rows - (jndex - index), n_rows);

                index = jndex;
            }
        }

        return true;
    }



	/// shortcut for matrix of integer values.
	typedef Matrix<INMOST_DATA_INTEGER_TYPE> iMatrix;
	/// shortcut for matrix of real values.
	typedef Matrix<INMOST_DATA_REAL_TYPE> rMatrix;
	/// shortcut for matrix of complex values.
	typedef Matrix<INMOST_DATA_CPLX_TYPE> cMatrix;
	/// shortcut for matrix of integer values in existing array.
	typedef Matrix<INMOST_DATA_INTEGER_TYPE,shell<INMOST_DATA_INTEGER_TYPE> > iaMatrix;
	/// shortcut for matrix of real values in existing array.
	typedef Matrix<INMOST_DATA_REAL_TYPE,shell<INMOST_DATA_REAL_TYPE> > raMatrix;
	/// shortcut for matrix of complex values in existing array.
	typedef Matrix<INMOST_DATA_CPLX_TYPE, shell<INMOST_DATA_CPLX_TYPE> > caMatrix;
	/// shortcut for symmetric matrix of integer values.
	typedef SymmetricMatrix<INMOST_DATA_INTEGER_TYPE> iSymmetricMatrix;
	/// shortcut for symmetric matrix of real values.
	typedef SymmetricMatrix<INMOST_DATA_REAL_TYPE> rSymmetricMatrix;
	/// shortcut for symmetric matrix of complex values.
	typedef SymmetricMatrix<INMOST_DATA_CPLX_TYPE> cSymmetricMatrix;
	/// shortcut for matrix of integer values in existing array.
	typedef SymmetricMatrix<INMOST_DATA_INTEGER_TYPE,shell<INMOST_DATA_INTEGER_TYPE> > iaSymmetricMatrix;
	/// shortcut for matrix of real values in existing array.
	typedef SymmetricMatrix<INMOST_DATA_REAL_TYPE,shell<INMOST_DATA_REAL_TYPE> > raSymmetricMatrix;
	/// shortcut for matrix of real values in existing array.
	typedef SymmetricMatrix<INMOST_DATA_CPLX_TYPE, shell<INMOST_DATA_CPLX_TYPE> > caSymmetricMatrix;
	/// return a matrix
	__INLINE iaMatrix iaMatrixMake(INMOST_DATA_INTEGER_TYPE * p, iaMatrix::enumerator n, iaMatrix::enumerator m) {return iaMatrix(shell<INMOST_DATA_INTEGER_TYPE>(p,n*m),n,m);}
	__INLINE raMatrix raMatrixMake(INMOST_DATA_REAL_TYPE * p, raMatrix::enumerator n, raMatrix::enumerator m) {return raMatrix(shell<INMOST_DATA_REAL_TYPE>(p,n*m),n,m);}
	__INLINE caMatrix caMatrixMake(INMOST_DATA_CPLX_TYPE* p, raMatrix::enumerator n, caMatrix::enumerator m) { return caMatrix(shell<INMOST_DATA_CPLX_TYPE>(p, n * m), n, m); }
	__INLINE iaSymmetricMatrix iaSymmetricMatrixMake(INMOST_DATA_INTEGER_TYPE * p, iaSymmetricMatrix::enumerator n) {return iaSymmetricMatrix(shell<INMOST_DATA_INTEGER_TYPE>(p,n*(n+1)/2),n);}
	__INLINE raSymmetricMatrix raSymmetricMatrixMake(INMOST_DATA_REAL_TYPE * p, raSymmetricMatrix::enumerator n) {return raSymmetricMatrix(shell<INMOST_DATA_REAL_TYPE>(p,n*(n+1)/2),n);}
	__INLINE caSymmetricMatrix caSymmetricMatrixMake(INMOST_DATA_CPLX_TYPE* p, caSymmetricMatrix::enumerator n) { return caSymmetricMatrix(shell<INMOST_DATA_CPLX_TYPE>(p, n * (n + 1) / 2), n); }
#if defined(USE_AUTODIFF)
	/// shortcut for matrix of variables with single unit entry of first order derivative.
	typedef Matrix<unknown> uMatrix;
	/// shortcut for matrix of variables with first order derivatives.
	typedef Matrix<variable> vMatrix;
	//< shortcut for matrix of variables with first and second order derivatives.
	typedef Matrix<hessian_variable> hMatrix;
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



	__INLINE uaMatrix uaMatrixMake(unknown * p, uaMatrix::enumerator n, uaMatrix::enumerator m) {return uaMatrix(shell<unknown>(p,n*m),n,m);}
	__INLINE vaMatrix vaMatrixMake(variable * p, vaMatrix::enumerator n, vaMatrix::enumerator m) {return vaMatrix(shell<variable>(p,n*m),n,m);}
	__INLINE haMatrix vaMatrixMake(hessian_variable * p, haMatrix::enumerator n, haMatrix::enumerator m) {return haMatrix(shell<hessian_variable>(p,n*m),n,m);}
	__INLINE uaSymmetricMatrix uaSymmetricMatrixMake(unknown * p, uaSymmetricMatrix::enumerator n) {return uaSymmetricMatrix(shell<unknown>(p,n*(n+1)/2),n);}
	__INLINE vaSymmetricMatrix vaSymmetricMatrixMake(variable * p, vaSymmetricMatrix::enumerator n) {return vaSymmetricMatrix(shell<variable>(p,n*(n+1)/2),n);}
	__INLINE haSymmetricMatrix vaSymmetricMatrixMake(hessian_variable * p, haSymmetricMatrix::enumerator n) {return haSymmetricMatrix(shell<hessian_variable>(p,n*(n+1)/2),n);}
#endif
}
/// Multiplication of matrix by constant from left.
/// @param coef Constant coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a constant.
template<typename typeB>
//INMOST::Matrix<typename INMOST::Promote<INMOST_DATA_REAL_TYPE,typeB>::type>
INMOST::MatrixMulCoef<typeB, INMOST_DATA_REAL_TYPE, typename INMOST::Promote<INMOST_DATA_REAL_TYPE, typeB>::type>
operator *(const INMOST_DATA_REAL_TYPE& coef, const INMOST::AbstractMatrixReadOnly<typeB>& other)
//{return other*coef;}
{return INMOST::MatrixMulCoef<typeB, INMOST_DATA_REAL_TYPE, typename INMOST::Promote<INMOST_DATA_REAL_TYPE, typeB>::type>(other, coef);}
#if defined(USE_AUTODIFF)
/// Multiplication of matrix by a unknown from left.
/// Takes account for derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB>
//INMOST::Matrix<typename INMOST::Promote<INMOST::unknown,typeB>::type>
INMOST::MatrixMulCoef<typeB, INMOST::unknown, typename INMOST::Promote<INMOST::unknown, typeB>::type>
operator *(const INMOST::unknown& coef, const INMOST::AbstractMatrixReadOnly<typeB>& other)
//{return other*coef;}
{return INMOST::MatrixMulCoef<typeB, INMOST::unknown, typename INMOST::Promote<INMOST::unknown, typeB>::type>(other, coef);}
/// Multiplication of matrix by a variable from left.
/// Takes account for derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB>
//INMOST::Matrix<typename INMOST::Promote<INMOST::variable,typeB>::type>
INMOST::MatrixMulCoef<typeB, INMOST::variable, typename INMOST::Promote<INMOST::variable, typeB>::type>
operator *(const INMOST::variable& coef, const INMOST::AbstractMatrixReadOnly<typeB>& other)
//{return other*coef;}
{return INMOST::MatrixMulCoef<typeB, INMOST::variable, typename INMOST::Promote<INMOST::variable, typeB>::type>(other, coef);}
/// Multiplication of matrix by a variable with first and
/// second order derivatives from left.
/// Takes account for first and second order derivatives of variable.
/// @param coef Variable coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a variable.
template<typename typeB>
//INMOST::Matrix<typename INMOST::Promote<INMOST::hessian_variable,typeB>::type>
INMOST::MatrixMulCoef<typeB, INMOST::hessian_variable, typename INMOST::Promote<INMOST::hessian_variable, typeB>::type>
operator *(const INMOST::hessian_variable& coef, const INMOST::AbstractMatrixReadOnly<typeB>& other)
//{return other*coef;}
{return INMOST::MatrixMulCoef<typeB, INMOST::hessian_variable, typename INMOST::Promote<INMOST::hessian_variable, typeB>::type>(other, coef);}
/// Multiplication of matrix by constant from left.
/// @param coef Constant coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a constant.
template<class A, typename typeB>
//INMOST::Matrix<typename INMOST::Promote<INMOST::variable,typeB>::type>
INMOST::MatrixMulShellCoef<typeB, INMOST::variable, typename INMOST::Promote<INMOST::variable, typeB>::type>
operator *(INMOST::shell_expression<A> const& coef, const INMOST::AbstractMatrixReadOnly<typeB>& other)
//{return other*INMOST::variable(coef);}
{return INMOST::MatrixMulShellCoef<typeB, INMOST::variable, typename INMOST::Promote<INMOST::variable, typeB>::type>(other, INMOST::variable(coef));}
#endif

template<typename T>
__INLINE bool check_nans(const INMOST::AbstractMatrixReadOnly<T> & A) {return A.CheckNans();}
template<typename T>
__INLINE bool check_infs(const INMOST::AbstractMatrixReadOnly<T> & A) {return A.CheckInfs();}
template<typename T>
__INLINE bool check_nans_infs(const INMOST::AbstractMatrixReadOnly<T> & A) {return A.CheckNansInfs();}

#endif //INMOST_DENSE_INCLUDED
