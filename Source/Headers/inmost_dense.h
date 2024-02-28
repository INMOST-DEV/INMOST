#ifndef INMOST_DENSE_INCLUDED
#define INMOST_DENSE_INCLUDED
#include "inmost_common.h"
#include "inmost_expression.h"
//#include "inmost_variable.h"
#include <iomanip>
#include <stdarg.h>

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
	template<> struct Promote<INMOST_DATA_CPLX_TYPE, value_reference> { typedef std::complex<INMOST_DATA_REAL_TYPE> type; };
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
	template<> struct Promote<value_reference, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
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
#else //USE_FP64
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
#endif //USE_FP64
#endif //USE_AUTODIFF

	// begin for matrices
	// all these pool expressions are needed to avoid dangling
	template<class A, class B> __INLINE         binary_pool_expression<branch_expression<A, B>, A, B> branch(bool cond, shell_expression<A> const& Left, shell_expression<B> const& Right) 
	{ 
		binary_pool<branch_expression<A, B>, A, B> pool(cond, Left, Right);
		return binary_pool_expression<branch_expression<A, B>, A, B>(pool); 
	}
	template<class A>          __INLINE          unary_pool_expression<const_branch_expression<A>, A> branch(bool cond, shell_expression<A> const& Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		unary_pool<const_branch_expression<A>, A> pool(cond, Left, Right);
		return unary_pool_expression<const_branch_expression<A>, A>(pool); 
	}
	template<class B>          __INLINE          unary_pool_expression<const_branch_expression<B>, B> branch(bool cond, INMOST_DATA_REAL_TYPE Left, shell_expression<B> const& Right) 
	{ 
		unary_pool<const_branch_expression<B>, B> pool(!cond, Right, Left);
		return unary_pool_expression<const_branch_expression<B>, B>(pool); 
	}
	__INLINE                                                                    INMOST_DATA_REAL_TYPE branch(bool cond, INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		return cond ? Left : Right; 
	}
	template<class A, class B> __INLINE       binary_pool_expression<addition_expression<A, B>, A, B> add(shell_expression<A> const& Left, shell_expression<B> const& Right) 
	{ 
		binary_pool< addition_expression<A, B>, A, B> pool(Left, Right);
		return binary_pool_expression< addition_expression<A, B>, A, B>(pool); 
	}
	template<class A>          __INLINE        unary_pool_expression<const_addition_expression<A>, A> add(shell_expression<A> const& Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		unary_pool<const_addition_expression<A>, A> pool(Left, Right);
		return unary_pool_expression<const_addition_expression<A>, A>(pool); 
	}
	template<class B>          __INLINE        unary_pool_expression<const_addition_expression<B>, B> add(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const& Right) 
	{ 
		unary_pool<const_addition_expression<B>, B> pool(Right, Left);
		return unary_pool_expression<const_addition_expression<B>, B>(pool);
	}
	__INLINE                                                                    INMOST_DATA_REAL_TYPE add(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) { return Left + Right; }
	template<class A, class B> __INLINE    binary_pool_expression<subtraction_expression<A, B>, A, B> sub(shell_expression<A> const& Left, shell_expression<B> const& Right) 
	{ 
		binary_pool< subtraction_expression<A, B>, A, B> pool(Left, Right);
		return binary_pool_expression< subtraction_expression<A, B>, A, B>(pool); 
	}
	template<class A>          __INLINE        unary_pool_expression<const_addition_expression<A>, A> sub(shell_expression<A> const& Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		unary_pool<const_addition_expression<A>, A> pool(Left, -Right);
		return unary_pool_expression<const_addition_expression<A>, A>(pool); 
	}
	template<class B>          __INLINE     unary_pool_expression<const_subtraction_expression<B>, B> sub(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const& Right) 
	{ 
		unary_pool<const_subtraction_expression<B>, B> pool(Right, Left);
		return unary_pool_expression<const_subtraction_expression<B>, B>(pool); 
	}
	__INLINE                                                                    INMOST_DATA_REAL_TYPE sub(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		return Left - Right; 
	}
	template<class A, class B> __INLINE binary_pool_expression<multiplication_expression<A, B>, A, B> mul(shell_expression<A> const& Left, shell_expression<B> const& Right) 
	{ 
		binary_pool<multiplication_expression<A, B>, A, B> pool(Left, Right);
		return binary_pool_expression<multiplication_expression<A, B>, A, B>(pool); 
	}
	template<class A>          __INLINE  unary_pool_expression<const_multiplication_expression<A>, A> mul(shell_expression<A> const& Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		unary_pool<const_multiplication_expression<A>, A> pool(Left, Right);
		return unary_pool_expression<const_multiplication_expression<A>, A>(pool); 
	}
	template<class B>          __INLINE  unary_pool_expression<const_multiplication_expression<B>, B> mul(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const& Right) 
	{ 
		unary_pool<const_multiplication_expression<B>, B> pool(Right, Left);
		return unary_pool_expression<const_multiplication_expression<B>, B>(pool); 
	}
	__INLINE                                                                    INMOST_DATA_REAL_TYPE mul(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		return Left * Right; 
	}
	template<class A, class B> __INLINE       binary_pool_expression<division_expression<A, B>, A, B> div(shell_expression<A> const& Left, shell_expression<B> const& Right) 
	{ 
		binary_pool<division_expression<A, B>, A, B> pool(Left, Right);
		return binary_pool_expression<division_expression<A, B>, A, B>(pool); 
	}
	template<class A>          __INLINE        unary_pool_expression<const_division_expression<A>, A> div(shell_expression<A> const& Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		unary_pool<const_division_expression<A>, A> pool(Left, Right);
		return unary_pool_expression<const_division_expression<A>, A>(pool); 
	}
	template<class B>          __INLINE            unary_pool_expression<reciprocal_expression<B>, B> div(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const& Right) 
	{ 
		unary_pool<reciprocal_expression<B>, B> pool(Right, Left);
		return unary_pool_expression<reciprocal_expression<B>, B>(pool); 
	}
	__INLINE                                                                    INMOST_DATA_REAL_TYPE div(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) 
	{ 
		return Left / Right; 
	}
	// end for matrices

	

	template<class A> struct OpRef { typedef A type; };
	template<> struct OpRef<variable> { typedef const_multivar_expression_reference type; };
	template<> struct OpRef<hessian_variable> { typedef const_hessian_multivar_expression_reference type; };

	template<class A, template<class> class Op> struct Op1 { typedef unary_pool_expression<Op<A>, A> type; };
	template<template<class> class Op> struct Op1<INMOST_DATA_REAL_TYPE, Op> { typedef INMOST_DATA_REAL_TYPE type; };
	template<template<class> class Op> struct Op1<INMOST_DATA_INTEGER_TYPE, Op> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<template<class> class Op> struct Op1<value_reference, Op> { typedef INMOST_DATA_REAL_TYPE type; };

	template<class A, class B> struct OpAdd { typedef binary_pool_expression<addition_expression<A, B>, A, B> type; };
	template<class A> struct OpAdd<A, INMOST_DATA_REAL_TYPE> { typedef unary_pool_expression<const_addition_expression<A>, A> type; };
	template<class A> struct OpAdd<A, INMOST_DATA_INTEGER_TYPE> { typedef unary_pool_expression<const_addition_expression<A>, A> type; };
	template<class A> struct OpAdd<A, value_reference> { typedef unary_pool_expression<const_addition_expression<A>, A> type; };
	template<class B> struct OpAdd<INMOST_DATA_REAL_TYPE, B> { typedef unary_pool_expression<const_addition_expression<B>, B> type; };
	template<class B> struct OpAdd<INMOST_DATA_INTEGER_TYPE, B> { typedef unary_pool_expression<const_addition_expression<B>, B> type; };
	template<class B> struct OpAdd<value_reference, B> { typedef unary_pool_expression<const_addition_expression<B>, B> type; };
	template<> struct OpAdd<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpAdd<INMOST_DATA_REAL_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpAdd<INMOST_DATA_REAL_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpAdd<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpAdd<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<> struct OpAdd<INMOST_DATA_INTEGER_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpAdd<value_reference, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpAdd<value_reference, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpAdd<value_reference, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };

	template<class A, class B> struct OpSub { typedef binary_pool_expression<subtraction_expression<A, B>, A, B> type; };
	template<class A> struct OpSub<A, INMOST_DATA_REAL_TYPE> { typedef unary_pool_expression<const_addition_expression<A>, A> type; };
	template<class A> struct OpSub<A, INMOST_DATA_INTEGER_TYPE> { typedef unary_pool_expression<const_addition_expression<A>, A> type; };
	template<class A> struct OpSub<A, value_reference> { typedef unary_pool_expression<const_addition_expression<A>, A> type; };
	template<class B> struct OpSub<INMOST_DATA_REAL_TYPE, B> { typedef unary_pool_expression<const_subtraction_expression<B>, B> type; };
	template<class B> struct OpSub<INMOST_DATA_INTEGER_TYPE, B> { typedef unary_pool_expression<const_subtraction_expression<B>, B> type; };
	template<class B> struct OpSub<value_reference, B> { typedef unary_pool_expression<const_subtraction_expression<B>, B> type; };
	template<> struct OpSub<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpSub<INMOST_DATA_REAL_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpSub<INMOST_DATA_REAL_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpSub<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpSub<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<> struct OpSub<INMOST_DATA_INTEGER_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpSub<value_reference, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpSub<value_reference, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpSub<value_reference, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };

	template<class A, class B> struct OpMul { typedef binary_pool_expression<multiplication_expression<A, B>, A, B> type; };
	template<class A> struct OpMul<A, INMOST_DATA_REAL_TYPE> { typedef unary_pool_expression<const_multiplication_expression<A>, A> type; };
	template<class A> struct OpMul<A, INMOST_DATA_INTEGER_TYPE> { typedef unary_pool_expression<const_multiplication_expression<A>, A> type; };
	template<class A> struct OpMul<A, value_reference> { typedef unary_pool_expression<const_multiplication_expression<A>, A> type; };
	template<class B> struct OpMul<INMOST_DATA_REAL_TYPE, B> { typedef unary_pool_expression<const_multiplication_expression<B>, B> type; };
	template<class B> struct OpMul<INMOST_DATA_INTEGER_TYPE, B> { typedef unary_pool_expression<const_multiplication_expression<B>, B> type; };
	template<class B> struct OpMul<value_reference, B> { typedef unary_pool_expression<const_multiplication_expression<B>, B> type; };
	template<> struct OpMul<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpMul<INMOST_DATA_REAL_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpMul<INMOST_DATA_REAL_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpMul<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpMul<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<> struct OpMul<INMOST_DATA_INTEGER_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpMul<value_reference, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpMul<value_reference, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpMul<value_reference, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	
	template<class A, class B> struct OpDiv { typedef binary_pool_expression<division_expression<A, B>, A, B> type; };
	template<class A> struct OpDiv<A, INMOST_DATA_REAL_TYPE> { typedef unary_pool_expression<const_division_expression<A>, A> type; };
	template<class A> struct OpDiv<A, INMOST_DATA_INTEGER_TYPE> { typedef unary_pool_expression<const_division_expression<A>, A> type; };
	template<class A> struct OpDiv<A, value_reference> { typedef unary_pool_expression<const_division_expression<A>, A> type; };
	template<class B> struct OpDiv<INMOST_DATA_REAL_TYPE, B> { typedef unary_pool_expression<reciprocal_expression<B>, B> type; };
	template<class B> struct OpDiv<INMOST_DATA_INTEGER_TYPE, B> { typedef unary_pool_expression<reciprocal_expression<B>, B> type; };
	template<class B> struct OpDiv<value_reference, B> { typedef unary_pool_expression<reciprocal_expression<B>, B> type; };
	template<> struct OpDiv<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpDiv<INMOST_DATA_REAL_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpDiv<INMOST_DATA_REAL_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpDiv<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpDiv<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<> struct OpDiv<INMOST_DATA_INTEGER_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpDiv<value_reference, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpDiv<value_reference, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpDiv<value_reference, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };

	template<class A, class B> struct OpCond { typedef binary_pool_expression<branch_expression<A, B>, A, B> type; };
	template<class A> struct OpCond<A, INMOST_DATA_REAL_TYPE> { typedef unary_pool_expression<const_branch_expression<A>, A> type; };
	template<class A> struct OpCond<A, INMOST_DATA_INTEGER_TYPE> { typedef unary_pool_expression<const_branch_expression<A>, A> type; };
	template<class A> struct OpCond<A, value_reference> { typedef unary_pool_expression<const_branch_expression<A>, A> type; };
	template<class B> struct OpCond<INMOST_DATA_REAL_TYPE, B> { typedef unary_pool_expression<const_branch_expression<B>, B> type; };
	template<class B> struct OpCond<INMOST_DATA_INTEGER_TYPE, B> { typedef unary_pool_expression<const_branch_expression<B>, B> type; };
	template<class B> struct OpCond<value_reference, B> { typedef unary_pool_expression<const_branch_expression<B>, B> type; };
	template<> struct OpCond<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpCond<INMOST_DATA_REAL_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpCond<INMOST_DATA_REAL_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpCond<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpCond<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_INTEGER_TYPE type; };
	template<> struct OpCond<INMOST_DATA_INTEGER_TYPE, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpCond<value_reference, INMOST_DATA_REAL_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpCond<value_reference, INMOST_DATA_INTEGER_TYPE> { typedef INMOST_DATA_REAL_TYPE type; };
	template<> struct OpCond<value_reference, value_reference> { typedef INMOST_DATA_REAL_TYPE type; };

	class AbstractMatrixBase
	{
	public:
		typedef INMOST_DATA_ENUM_TYPE enumerator; 
		typedef basic_expression::merger_type merger_type;
		static merger_type& GetMerger() { return *merger; }
	protected:
		static thread_private<merger_type> merger;
	};

	template<typename VarA, typename RetA, typename VarB, typename RetB>
	__INLINE typename Promote<VarA, VarB>::type DotProduct(const AbstractMatrixReadOnly<VarA, RetA>& A, const AbstractMatrixReadOnly<VarB, RetB>& B);
#if defined(USE_AUTODIFF)
	template<typename RetB>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>& A, const AbstractMatrixReadOnly<variable, RetB>& B);
	template<typename RetA>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<variable, RetA>& A, const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>& B);
	template<typename RetA, typename RetB>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<variable, RetA>& A, const AbstractMatrixReadOnly<variable, RetB>& B);	
#endif

	template<typename Var, typename RetType>
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
		/// Obtain number of rows.
		/// @return Number of rows.
		virtual enumerator Rows() const = 0;
		/// Obtain number of columns.
		/// @return Number of columns.
		virtual enumerator Cols() const = 0;
		/// Compute element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		virtual RetType compute(enumerator i, enumerator j) const = 0;
		//TODO: comment this and use compute in internal functions
		RetType operator()(enumerator i, enumerator j) const {return compute(i,j);}
		/// Check all matrix entries for not a number.
		/// Also checks derivatives for matrices of variables.
		bool CheckNans() const 
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					if (check_nans(compute(i, j))) return true;
			return false;
		}
		/// Check all matrix entries for infinity.
		/// Also checks derivatives for matrices of variables.
		bool CheckInfs() const 
		{ 
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					if (check_infs(compute(i, j))) return true;
			return false;
		}
		/// Check all matrix entries for not a number and infinity.
		/// Also checks derivatives for matrices of variables.
		bool CheckNansInfs() const 
		{ 
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					if (check_nans_infs(compute(i, j))) return true;
			return false;
		}
		/// Matrix determinant
		Var Det() const
		{
			assert(Rows() == Cols());
			const INMOST_DATA_REAL_TYPE eps = 1.0e-13;
			const enumerator n = Rows();
			Matrix<Var> A(*this);
			Matrix<Var> row_visited(n,1,-1);
			Var ret = 1.0, coef;
			INMOST_DATA_REAL_TYPE sign;
			for (enumerator d = 0; d < n; ++d)
			{
				enumerator r = 0;
				sign = 1;
				// find unvisited row with non-zero value in column d
				while(row_visited(r,0) != -1 || fabs(get_value(A(r, d))) < eps)
				{
					if(r == n-1) return Var(0);
					if(row_visited(r,0) == -1) sign *= -1;
					++r;
				}
				row_visited(r,0) = d;
				ret *= sign * A(r, d);
				for (enumerator i = 0; i < n; ++i) if(row_visited(i,0) == -1)
				{
					coef = A(i, d) / A(r, d);
					if(fabs(coef) > eps)
						for (enumerator j = 0; j < n; ++j)
							A(i, j) = A(i, j) - coef * A(r, j);
				}
			}
			return ret;
		}
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
		/// Singular value decomposition.
		/// Specific for complex matrices.
		bool cSVD(AbstractMatrix<Var>& U, AbstractMatrix<Var>& Sigma, AbstractMatrix<Var>& V) const;
		/// Transpose current matrix.
		/// @return Transposed matrix.
		//Matrix<Var> Transpose() const;
		AbstractMatrixTranspose<Var, RetType> Transpose() const { return AbstractMatrixTranspose<Var, RetType>(*this); }
		/// Transpose and conjugate current matrix.
		/// @return Transposed and conjugated matrix.
		AbstractMatrixConjugateTranspose<Var, RetType> ConjugateTranspose() const { return AbstractMatrixConjugateTranspose<Var, RetType>(*this); }
		/// Conjugate  current matrix.
		/// @return Conjugated matrix.
		AbstractMatrixConjugate<Var, RetType> Conjugate() const { return AbstractMatrixConjugate<Var, RetType>(*this); }
		/// Cross-product operation for a vector.
		/// Both right hand side and left hand side should be a vector
		/// @param other The right hand side of cross product.
		/// @return The cross product of current and right hand side vector.
		/// \warning Works for 3x1 vector and 3xm m-vectors as right hand side.
		template<typename typeB, typename retB>
		Matrix<typename Promote<Var, typeB>::type>
			CrossProduct(const AbstractMatrixReadOnly<typeB, retB>& other) const;
		/// Transformation matrix from current vector to provided vector using shortest arc rotation.
		/// @param other Vector to transform to.
		/// @return A sqaure (rotation) matrix that transforms current vector into right hand side vector.
		/// \warning Works only for 3x1 vectors.
		template<typename typeB, typename retB>
		Matrix<typename Promote<Var, typeB>::type>
			Transform(const AbstractMatrixReadOnly<typeB, retB>& other) const;
		/// Subtract a matrix.
		/// @param other The matrix to be subtracted.
		/// @return Difference of current matrix with another matrix.
		template<typename typeB, typename retB>
		MatrixDifference<Var, typeB, RetType, retB>
			operator-(const AbstractMatrixReadOnly<typeB, retB>& other) const
		{ return MatrixDifference<Var, typeB, RetType, retB>(*this, other);	}
		/// Add a matrix.
		/// @param other The matrix to be added.
		/// @return Sum of current matrix with another matrix.
		template<typename typeB, typename retB>
		MatrixSum<Var, typeB, RetType, retB>
			operator+(const AbstractMatrixReadOnly<typeB, retB>& other) const
		{ return MatrixSum<Var, typeB, RetType, retB>(*this, other); }
		/// Multiply the matrix by another matrix.
		/// @param other Matrix to be multiplied from right.
		/// @return Matrix multipled by another matrix.
		template<typename typeB, typename retB>
		MatrixMul<Var, typeB, RetType, retB>
			operator*(const AbstractMatrixReadOnly<typeB, retB>& other) const
		{ return MatrixMul<Var, typeB, RetType, retB>(*this, other); }
		/// Unary minus. Change sign of each element of the matrix.
		MatrixUnaryMinus<Var, RetType> operator-() const { return MatrixUnaryMinus<Var, RetType>(*this); }
		/// Kronecker product, latex symbol \otimes.
		/// @param other Matrix on the right of the Kronecker product.
		/// @return Result of the Kronecker product of current and another matrix.
		template<typename typeB, typename retB>
		KroneckerProduct<Var, typeB, RetType, retB>
			Kronecker(const AbstractMatrixReadOnly<typeB, retB>& other) const
		{ return KroneckerProduct<Var, typeB, RetType, retB>(*this, other); }
		/// Inverts matrix using Crout-LU decomposition with full pivoting for
		/// maximum element. Works well for non-singular matrices, for singular
		/// matrices look see Matrix::PseudoInvert.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr equal to a positive integer that represents the row
		///             on which the failure occurred (starting from 1), in case of no failure *ierr = 0.
		/// @return Inverse matrix.
		/// @see Matrix::PseudoInvert.
		/// \todo (test) Activate and test implementation with Solve.
		Matrix<Var> Invert(int* ierr = NULL) const
		{
			Matrix<Var> ret(Cols(), Rows());
			ret = Solve(MatrixUnit<Var>(Rows()), ierr);
			return ret;
		}
		/// Inverts symmetric positive-definite matrix using Cholesky decomposition.
		Matrix<Var> CholeskyInvert(int* ierr = NULL) const
		{
			Matrix<Var> ret(Rows(), Rows());
			ret = CholeskySolve(MatrixUnit<Var>(Rows()), ierr);
			return ret;
		}
		/// Finds X in A*X=B, where A and B are general matrices.
		/// Converts system into A^T*A*X=A^T*B.
		/// Inverts matrix A^T*A using Crout-LU decomposition with full pivoting for
		/// maximum element. Works well for non-singular matrices, for singular
		/// matrices see Matrix::PseudoInvert.
		/// @param B Right hand side matrix.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr equal to a positive integer that represents the row
		///             on which the failure occurred (starting from 1), in case of no failure *ierr = 0.
		/// @return Inverse matrix,
		/// @see Matrix::PseudoInvert.
		/// \warning Number of rows in matrices A and B should match.
		/// \todo
		/// 1. Test implementation.
		/// 2. Maximum product transversal + block pivoting instead of pivoting
		///    by maximum element.
		template<typename typeB, typename retB>
		Matrix<typename Promote<Var, typeB>::type>
			Solve(const AbstractMatrixReadOnly<typeB, retB>& B, int* ierr = NULL) const;
		/// Finds X in A*X=B, where A is a square symmetric positive definite matrix.
		/// Uses Cholesky decomposition algorithm.
		/// @param B Right hand side matrix.
		/// @param ierr Returns error on fail. If ierr is NULL, then throws an exception.
		///             If *ierr == -1 on input, then prints out information in case of failure.
		///             In case of failure *ierr equal to a positive integer that represents the row
		///             on which the failure occurred (starting from 1), in case of no failure *ierr = 0.
		/// @see Matrix::PseudoInvert.
		/// @return Inverse matrix,
		template<typename typeB, typename retB>
		Matrix<typename Promote<Var, typeB>::type>
			CholeskySolve(const AbstractMatrixReadOnly<typeB, retB>& B, int* ierr = NULL) const;
		/// Calculate sum of the diagonal elements of the matrix.
		/// @return Trace of the matrix.
		Var Trace() const
		{
			assert(Cols() == Rows());
			Var ret = 0.0;
			for (enumerator i = 0; i < Rows(); ++i) ret += compute(i, i);
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
					//if (__isinf__(real_part(get_value(compute(k, l)))))
					//	sout << std::setw(12) << "inf";
					//else if (std::isnan(get_value(compute(k, l))))
					//	sout << std::setw(12) << "nan";
					//else 
					if (fabs(get_value(compute(k, l))) > threshold)
						sout << std::setw(12) << get_value(compute(k, l));
					else
						sout << std::setw(12) << 0;
					sout << " ";
				}
				sout << std::endl;
			}
		}
		/// Check if the matrix is symmetric.
		/// @return Returns true if the matrix is symmetric, otherwise false.
		bool isSymmetric(INMOST_DATA_REAL_TYPE eps = 1.0e-7) const
		{
			if (Rows() != Cols()) return false;
			for (enumerator k = 0; k < Rows(); ++k)
			{
				for (enumerator l = k + 1; l < Rows(); ++l)
					if (fabs(compute(k, l) - compute(l, k)) > eps)
						return false;
			}
			return true;
		}
		/// Computes dot product by summing up multiplication of entries with the
		/// same indices in the current and the provided matrix.
		/// @param other Provided matrix.
		/// @return Dot product of two matrices.
		template<typename typeB, typename retB>
		typename Promote<Var, typeB>::type
			DotProduct(const AbstractMatrixReadOnly<typeB, retB>& other) const
		{
			/*
			assert(Cols() == other.Cols());
			assert(Rows() == other.Rows());
			typename Promote<Var, typeB>::type ret = 0.0;
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					ret += compute(i, j) * other.compute(i, j);
			return ret;
			*/
			return INMOST::DotProduct(*this, other);
		}
		/// Computes dot product by summing up multiplication of entries with the
		/// same indices in the current and the provided matrix.
		/// @param other Provided matrix.
		/// @return Dot product of two matrices.
		template<typename typeB, typename retB>
		typename Promote<Var, typeB>::type operator ^(const AbstractMatrixReadOnly<typeB, retB>& other) const
		{
			return DotProduct(other);
		}
		/// Computes frobenius norm of the matrix.
		/// @return Frobenius norm of the matrix.
		typename SelfPromote<Var>::type FrobeniusNorm() const
		{
			typename SelfPromote<Var>::type ret = 0;
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					ret += pow(compute(i, j), 2);
			return sqrt(ret);
		}
		/// Computes maximum absolute value of the matrix.
		/// @return Maximum norm of the matrix.
		Var MaxNorm() const
		{
			Var ret = 0;
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					ret = std::max<Var>(ret, compute(i, j));
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
		template<typename typeB, typename retB>
		Matrix<typename Promote<Var, typeB>::type>
			PseudoSolve(const AbstractMatrixReadOnly<typeB, retB>& B, INMOST_DATA_REAL_TYPE tol = 0, int* ierr = NULL) const;
		/// Extract submatrix of a matrix.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param ibeg Starting row in the original matrix.
		/// @param iend Last row (excluded) in the original matrix.
		/// @param jbeg Starting column in the original matrix.
		/// @param jend Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		Matrix<Var> ExtractSubMatrix(enumerator ibeg, enumerator iend, enumerator jbeg, enumerator jend) const
		{
			assert(ibeg < Rows());
			assert(iend < Rows());
			assert(jbeg < Cols());
			assert(jend < Cols());
			Matrix<Var> ret(iend - ibeg, jend - jbeg);
			for (enumerator i = ibeg; i < iend; ++i)
			{
				for (enumerator j = jbeg; j < jend; ++j)
					ret(i - ibeg, j - jbeg) = compute(i, j);
			}
			return ret;
		}
		/// Change representation of the matrix into matrix of another size.
		/// Useful to change representation from matrix into vector and back.
		/// Replaces original number of columns and rows with a new one.
		/// @return Matrix with same entries and provided number of rows and columns.
		AbstractMatrixRepack<Var, RetType> Repack(enumerator rows, enumerator cols) const {return AbstractMatrixRepack<Var, RetType>(*this, rows, cols);}
		/// Extract submatrix of a matrix for in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, i in [ibeg,iend),
		/// and j in [jbeg,jend).
		/// @param first_row Starting row in the original matrix.
		/// @param last_row Last row (excluded) in the original matrix.
		/// @param first_col Starting column in the original matrix.
		/// @param last_col Last column (excluded) in the original matrix.
		/// @return Submatrix of the original matrix.
		AbstractSubMatrix<Var, RetType> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col) const
		{ return AbstractSubMatrix<Var, RetType>(*this, first_row, last_row, first_col, last_col); }
		/// Define matrix as a part of a matrix of larger size with in-place manipulation of elements.
		/// Let A = {a_ij}, i in [0,n), j in [0,m) be original matrix.
		/// Then the method returns B = {a_ij}, if i in [offset_row,offset_row+n),
		/// and j in [offset_col,offset_col+m) and B = {0} otherwise.
		/// @param nrows Number of rows in larger matrix.
		/// @param ncols Number of columns in larger matrix.
		/// @param offset_row Offset for row number.
		/// @param offset_col Offset for column number.
		/// @return Submatrix of the original matrix.
		AbstractBlockOfMatrix<Var, RetType> BlockOf(enumerator nrows, enumerator ncols, enumerator offset_row, enumerator offset_col) const
		{ return AbstractBlockOfMatrix<Var, RetType>(*this, nrows, ncols, offset_row, offset_col); }
		/// Multiply the matrix by a coefficient.
		/// @param coef Coefficient.
		/// @return Matrix multiplied by the coefficient.
		template<typename typeB>
		MatrixMulCoef<Var, typeB, RetType, typename Promote<Var, typeB>::type>
			operator*(const typeB& coef) const
		{ return MatrixMulCoef<Var, typeB, RetType, typename Promote<Var, typeB>::type>(*this, coef); }
		
#if defined(USE_AUTODIFF)
		/// Multiply the matrix by a coefficient.
		/// @param coef Coefficient.
		/// @return Matrix multiplied by the coefficient.
		template<class A>
		MatrixMulCoef<Var, variable, RetType, typename Promote<Var, variable>::type>
			operator*(shell_expression<A> const& coef) const { return operator*(variable(coef)); }
#endif //USE_AUTODIFF
		
		
#if defined(USE_AUTODIFF)
		/// Divide the matrix by a coefficient of a different type.
		/// @param coef Coefficient.
		/// @return Matrix divided by the coefficient.
		template<class A>
		MatrixDivCoef<Var, variable, RetType, typename Promote<Var, variable>::type>
			operator/(shell_expression<A> const& coef) const { return operator/(variable(coef)); }
#endif //USE_AUTODIFF
		
		/// Performs B^{-1}*A, multiplication by inverse matrix from left.
		/// Throws exception if matrix is not invertible. See Mesh::PseudoSolve for
		/// singular matrices.
		/// @param other Matrix to be inverted and multiplied from left.
		/// @return Multiplication of current matrix by inverse of another
		/// @see Matrix::PseudoInvert.
		template<typename typeB, typename retB>
		Matrix<typename Promote<Var, typeB>::type>
			operator/(const AbstractMatrixReadOnly<typeB, retB>& other) const
		{
			Matrix<typename Promote<Var, typeB>::type> ret(other.Cols(), Cols());
			ret = other.Solve(*this);
			return ret;
		}
		/// Divide the matrix by a coefficient of a different type.
		/// @param coef Coefficient.
		/// @return Matrix divided by the coefficient.
		template<typename typeB>
		MatrixDivCoef<Var, typeB, RetType, typename Promote<Var, typeB>::type>
			operator/(const typeB& coef) const
		{
			return MatrixDivCoef<Var, typeB, RetType, typename Promote<Var, typeB>::type>(*this, coef);
		}
		/// Concatenate B matrix as columns of current matrix.
		/// Assumes that number of rows of current matrix is
		/// equal to number of rows of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatRows
		template<typename VarB, typename RetB>
		AbstractMatrixConcatCols<Var, VarB, RetType, RetB>
			ConcatCols(const AbstractMatrixReadOnly<VarB, RetB>& B) const
		{
			return AbstractMatrixConcatCols<Var, VarB, RetType, RetB>(*this, B);
		}
		/// Concatenate B matrix as rows of current matrix.
		/// Assumes that number of colums of current matrix is
		/// equal to number of columns of B matrix.
		/// @param B Matrix to be concatenated to current matrix.
		/// @return Result of concatenation of current matrix and parameter.
		/// @see Matrix::ConcatCols
		template<typename VarB, typename RetB>
		AbstractMatrixConcatRows<Var,VarB,RetType,RetB> 
			ConcatRows(const AbstractMatrixReadOnly<VarB, RetB>& B) const
		{
			return AbstractMatrixConcatRows<Var, VarB, RetType, RetB>(*this, B);
		}
		/// Convert into matrix.
		/// @return Matrix with same entries.
		Matrix<Var> MakeMatrix()
		{
			Matrix<Var> ret(Rows(), Cols());
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					ret(i, j) = compute(i, j);
			return ret;
		}
		/// Destructor
		virtual ~AbstractMatrixReadOnly() {}
	};
	/// Abstract class for a matrix,
	/// used to abstract away all the data storage and access
	/// and provide common implementation of the algorithms.
	template<typename Var>
	class AbstractMatrix : public AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;	
	 	using AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>::operator ();
		using AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>::Transpose;
		using AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>::Repack;
		using AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>::ConcatRows;
		using AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>::ConcatCols;
		using AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>::Rows;
		using AbstractMatrixReadOnly<Var, typename OpRef<Var>::type>::Cols;
		/// Construct empty matrix.
		AbstractMatrix() {}
		/// Assign matrix of the same type.
		/// @param other Another matrix of the same type.
		/// @return Reference to matrix.
		AbstractMatrix & operator =(AbstractMatrix const & other)
		{
			Resize(other.Rows(),other.Cols());
			for (enumerator i = 0; i < other.Rows(); ++i)
				for (enumerator j = 0; j < other.Cols(); ++j)
					assign(get(i, j), other.get(i, j));
			return *this;
		}
		/// Assign matrix of another type.
		/// @param other Another matrix of different type.
		/// @return Reference to matrix.
		template<typename typeB, typename retB>
		AbstractMatrix & operator =(AbstractMatrixReadOnly<typeB, retB> const & other)
		{
			Resize(other.Rows(),other.Cols());
			for (enumerator i = 0; i < other.Rows(); ++i)
				for (enumerator j = 0; j < other.Cols(); ++j)
					assign(get(i, j), other.compute(i, j));
			return *this;
		}
		/// Assign value to all entries of the matrix.
		/// @param b Assigned value.
		/// @return Reference to matrix.
		AbstractMatrix & operator =(Var const & b)
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					assign(get(i, j), b);
			return *this;
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to element.
		__INLINE Var& operator () (enumerator i, enumerator j) { return get(i, j); }
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		virtual Var& get(enumerator i, enumerator j) = 0;
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		virtual const Var& get(enumerator i, enumerator j) const = 0;
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE typename OpRef<Var>::type compute(enumerator i, enumerator j) const { return typename OpRef<Var>::type(get(i, j)); }
		/// Resize the matrix into different size.
		/// @param nrows New number of rows.
		/// @param ncols New number of columns.
		virtual void Resize(enumerator rows, enumerator cols) = 0;
		/// Set all the elements of the matrix to zero.
		virtual void Zero()
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					assign(get(i, j), 0.0);
		}
		///Exchange contents of two matrices.
		virtual void Swap(AbstractMatrix<Var> & b);
		/// Subtract a matrix and store result in the current.
		/// @param other The matrix to be subtracted.
		/// @return Reference to the current matrix.
		template<typename typeB, typename retB>
		AbstractMatrix & operator-=(const AbstractMatrixReadOnly<typeB, retB> & other)
		{
			assert(Rows() == other.Rows());
			assert(Cols() == other.Cols());
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					assign(get(i, j), sub(compute(i, j), other.compute(i, j)));
			return *this;
		}
		/// Add a matrix and store result in the current.
		/// @param other The matrix to be added.
		/// @return Reference to the current matrix.
		template<typename typeB, typename retB>
		AbstractMatrix & operator+=(const AbstractMatrixReadOnly<typeB, retB> & other)
		{
			assert(Rows() == other.Rows());
			assert(Cols() == other.Cols());
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					assign(get(i, j), add(compute(i, j), other.compute(i, j)));
			return *this;
		}
		/// Multiply matrix with another matrix in-place.
		/// @param B Another matrix to the right in multiplication.
		/// @return Reference to current matrix.
		template<typename typeB, typename retB>
		AbstractMatrix & operator*=(const AbstractMatrixReadOnly<typeB, retB> & B)
		{
			(*this) = (*this) * B;
			return *this;
		}
		/// Multiply the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator*=(typeB coef)
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					assign(get(i, j), mul(compute(i, j), coef));
			return *this;
		}
		/// Divide the matrix by the coefficient of the same type and store the result.
		/// @param coef Coefficient.
		/// @return Reference to the current matrix.
		template<typename typeB>
		AbstractMatrix & operator/=(typeB coef)
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = 0; j < Cols(); ++j)
					assign(get(i, j), div(compute(i, j), coef));
			return *this;
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
		SubMatrix<Var> operator()(enumerator first_row, enumerator last_row, enumerator first_col, enumerator last_col)
		{ return SubMatrix<Var>(*this, first_row, last_row, first_col, last_col); }
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
		/// Destructor
		virtual ~AbstractMatrix() {};
	};
	
	
	template<typename Var, typename storage_type>
	class SymmetricMatrix : public AbstractMatrix<Var>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
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
		/// \warning The matrix does not necessary have zero entries.
		SymmetricMatrix(enumerator pn) : space(pn*(pn+1)/2), n(pn) {}
		/// Construct a matrix with provided elements.
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
		/// @param pn Number of rows and columns.
		/// @param c Value to fill the matrix.
		SymmetricMatrix(enumerator pn, const Var & c) : space(pn*(pn+1)/2,c), n(pn) {}
		/// Copy matrix.
		/// @param other Another matrix of the same type.
		SymmetricMatrix(const SymmetricMatrix & other) : space(other.space), n(other.n) {}
		/// Construct matrix from matrix of different type.
		/// Function assumes that the other matrix is square and symmetric.
		/// Copies only top-right triangular part.
		/// Uses assign function declared in inmost_expression.h.
		/// Copies derivative information if possible.
		/// @param other Another matrix of different type.
		template<typename typeB, typename retB>
		SymmetricMatrix(const AbstractMatrixReadOnly<typeB, retB> & other) : space(other.Rows()*(other.Rows()+1)/2), n(other.Rows())
		{
			assert(other.Rows() == other.Cols());
			for (enumerator i = 0; i < n; ++i)
				for (enumerator j = i; j < n; ++j)
					assign(get(i, j), other.compute(i, j));
		}
		/// Delete matrix.
		~SymmetricMatrix() {}
		/// Resize the matrix into different size.
		/// Number of rows must match the number of columns for symmetric matrix.
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
		/// Set all the elements of the matrix to zero.
		void Zero()
		{
			for (enumerator i = 0; i < Rows(); ++i)
				for (enumerator j = i; j < Cols(); ++j)
					assign(get(i, j), 0.0);
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
		SymmetricMatrix& operator =(SymmetricMatrix&& other)
		{
			if (this != &other)
			{
				space = std::move(other.space);
				n = other.n;
			}
			return *this;
		}
		/// Assign matrix of another type.
		/// Function assumes that the other matrix is square and symmetric.
		/// Copies only top-right triangular part.
		/// @param other Another matrix of different type.
		/// @return Reference to matrix.
		template<typename typeB, typename retB>
		SymmetricMatrix & operator =(AbstractMatrixReadOnly<typeB, retB> const & other)
		{
			assert(other.Rows() == other.Cols());
			if (Rows() != other.Rows())
				space.resize((other.Rows() + 1) * other.Rows() / 2);
			n = other.Rows();
			for (enumerator i = 0; i < other.Rows(); ++i)
				for (enumerator j = i; j < other.Cols(); ++j)
					assign(get(i, j), other.compute(i, j));
			return *this;
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to const element.
		__INLINE Var& get(enumerator i, enumerator j)
		{
			assert(i < n);
			assert(j < n);
			if (i > j) std::swap(i, j);
			return space[j + n * i - i * (i + 1) / 2];
		}
		/// Access element of the matrix by row and column indices.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to const element.
		__INLINE const Var& get(enumerator i, enumerator j) const
		{
			assert(i < n);
			assert(j < n);
			if (i > j) std::swap(i, j);
			return space[j + n * i - i * (i + 1) / 2];
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
			SymmetricMatrix<Var> Kc(matsize);
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
		typedef typename AbstractMatrixBase::enumerator enumerator;
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
					get(k-1,l) = get(k,l);
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
					get(k-shift,l) = get(k,l);
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
					tmp(k,l) = get(k,l);
				for(enumerator l = col+1; l < m; ++l)
					tmp(k,l-1) = get(k,l);
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
					tmp(k,l) = get(k,l);
				for(enumerator l = last; l < m; ++l)
					tmp(k,l-shift) = get(k,l);
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
					tmp(k,l) = get(k,l);
				for(enumerator l = lastcol; l < m; ++l)
					tmp(k,l-shiftcol) = get(k,l);
			}
			for(enumerator k = lastrow; k < n; ++k)
			{
				for(enumerator l = 0; l < firstcol; ++l)
					tmp(k-shiftrow,l) = get(k,l);
				for(enumerator l = lastcol; l < m; ++l)
					tmp(k-shiftrow,l-shiftcol) = get(k,l);
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
		/// \warning The matrix does not necessary have zero entries.
		Matrix(enumerator pn, enumerator pm) : space(pn*pm), n(pn), m(pm) {}
		/// Construct a matrix with provided sizes and fills with value.
		/// @param pn Number of rows.
		/// @param pm Number of columns.
		/// @param c Value to fill the matrix.
		Matrix(enumerator pn, enumerator pm, const Var & c) : space(pn*pm,c), n(pn), m(pm) {}
		/// Construct a matrix with provided elements.
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
		}
		Matrix(Matrix && other) : space(std::move(other.space)), n(other.n), m(other.m)
		{
		}
		/// Construct matrix from matrix of different type.
		/// Uses assign function declared in inmost_expression.h.
		/// Copies derivative information if possible.
		/// @param other Another matrix of different type.
		template<typename typeB, typename retB>
		Matrix(const AbstractMatrixReadOnly<typeB, retB> & other) : space(other.Cols()*other.Rows()), n(other.Rows()), m(other.Cols())
		{
			for(enumerator i = 0; i < n; ++i)
				for(enumerator j = 0; j < m; ++j)
					assign(get(i,j),other.compute(i,j));
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
				if (n * m != other.n * other.m)
					space.resize(other.n * other.m);
				for (enumerator i = 0; i < other.n * other.m; ++i)
					space[i] = other.space[i];
				n = other.n;
				m = other.m;
			}
			return *this;
		}
		/*
		Matrix& operator =(Matrix && other)
		{
			if (this != &other)
			{
				space = std::move(other.space);
				n = other.n;
				m = other.m;
			}
			return *this;
		}
		*/
		/// Assign matrix of another type.
		/// @param other Another matrix of different type.
		/// @return Reference to matrix.
		template<typename typeB>
		Matrix& operator =(AbstractMatrix<typeB> const& other)
		{
			if (Cols() * Rows() != other.Cols() * other.Rows())
				space.resize(other.Cols() * other.Rows());
			n = other.Rows();
			m = other.Cols();
			for (enumerator i = 0; i < n; ++i)
				for (enumerator j = 0; j < m; ++j)
					assign(get(i, j), other.get(i, j));
			return *this;
		}
		/// Assign matrix of another type.
		/// @param other Another matrix of different type.
		/// @return Reference to matrix.
		template<typename typeB, typename retB>
		Matrix & operator =(AbstractMatrixReadOnly<typeB, retB> const & other)
		{
			if (Cols() * Rows() != other.Cols() * other.Rows())
				space.resize(other.Cols() * other.Rows());
			n = other.Rows();
			m = other.Cols();
			for (enumerator i = 0; i < n; ++i)
				for (enumerator j = 0; j < m; ++j)
					assign(get(i, j), other.compute(i, j));
			return *this;
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE Var& get(enumerator i, enumerator j)
		{
			assert(i < n);
			assert(j < m);
			assert(i * m + j < n * m); //overflow check?
			return space[i * m + j];
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE const Var & get(enumerator i, enumerator j) const
		{
			assert(i < n);
			assert(j < m);
			assert(i * m + j < n* m); //overflow check?
			return space[i * m + j];
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
		{ return MatrixCol<Var>(pn, c); 	}
		/// Joint diagonalization algorithm by Cardoso.
		/// Source http://perso.telecom-paristech.fr/~cardoso/Algo/Joint_Diag/joint_diag_r.m
		/// Current matrix should have size n by n*m
		/// And represent concatenation of m n by n matrices.
		/// Current matrix is replaced by diagonalized matrices.
		/// For correct result it is required that input matrices are
		/// exactly diagonalizable, otherwise the result may be approximate.
		/// @param threshold Optional small number.
		/// @return A unitary n by n matrix V used to diagonalize array of
		/// initial matrices. Current matrix is replaced by concatenation of
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
	};
	/// This class allows for in-place operations on submatrix of the matrix elements.
	template<typename Var>
	class SubMatrix : public AbstractMatrix<Var>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
		using AbstractMatrix<Var>::operator =;
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
		SubMatrix(AbstractMatrix<Var> & rM, enumerator first_row, enumerator last_row, enumerator first_column, enumerator last_column) : M(&rM), brow(first_row), erow(last_row), bcol(first_column), ecol(last_column) {}
		SubMatrix(const SubMatrix & b) : M(b.M), brow(b.brow), erow(b.erow), bcol(b.bcol), ecol(b.ecol) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to const element.
		__INLINE const Var & get(enumerator i, enumerator j) const
		{
			assert(i < Rows());
			assert(j < Cols());
			assert(i * Cols() + j < Rows()* Cols()); //overflow check?
			return M->get(i + brow, j + bcol);
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to const element.
		__INLINE Var& get(enumerator i, enumerator j)
		{
			assert(i < Rows());
			assert(j < Cols());
			assert(i * Cols() + j < Rows() * Cols()); //overflow check?
			return M->get(i + brow, j + bcol);
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. SubMatrix cannot change it's size,
		/// since it just points to a part of the original matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
            (void)cols; (void)rows;
			if( Cols() != cols || Rows() != rows) throw Impossible;
		}
	};

	/// This class allows for in-place operations on submatrix of the matrix elements.
	template<typename Var, typename Ret>
	class AbstractSubMatrix : public AbstractMatrixReadOnly<Var, Ret>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<Var, Ret> * M;
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
		AbstractSubMatrix(const AbstractMatrixReadOnly<Var, Ret> & rM, enumerator first_row, enumerator last_row, enumerator first_column, enumerator last_column) : M(&rM), brow(first_row), erow(last_row), bcol(first_column), ecol(last_column)
		{}
		AbstractSubMatrix(const AbstractSubMatrix & b) : M(b.M), brow(b.brow), erow(b.erow), bcol(b.bcol), ecol(b.ecol) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE Ret compute(enumerator i, enumerator j) const
		{
			assert(i < Rows());
			assert(j < Cols());
			assert(i*Cols()+j < Rows()*Cols()); //overflow check?
			return M->compute(i+brow,j+bcol);
		}
	};
	
	/// This class allows to address a matrix as a block of an empty matrix of larger size.
	template<typename Var, typename Ret>
	class AbstractBlockOfMatrix : public AbstractMatrixReadOnly<Var, typename Op1<Ret, const_branch_expression>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<Var,Ret> * M;
		enumerator nrows; //< Number of rows in larger matrix.
		enumerator ncols; //< Number of colums in larger matrix.
		enumerator orow; //< Row offset in larger matrix.
		enumerator ocol; //< Column offset in larger matrix.
	public:
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
		AbstractBlockOfMatrix(const AbstractMatrixReadOnly<Var,Ret> & rM, enumerator num_rows, enumerator num_cols, enumerator offset_row, enumerator offset_col) : M(&rM), nrows(num_rows), ncols(num_cols), orow(offset_row), ocol(offset_col)
		{}
		AbstractBlockOfMatrix(const AbstractBlockOfMatrix & b) : M(b.M), nrows(b.nrows), ncols(b.ncols), orow(b.orow), ocol(b.ocol) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Column index.
		/// @param j Row index.
		/// @return Reference to constant element.
		__INLINE typename Op1<Ret, const_branch_expression>::type compute(enumerator i, enumerator j) const
		{
			assert(i < Rows());
			assert(j < Cols());
			return branch(!(i < orow || i >= orow + M->Rows() || j < ocol || j >= ocol + M->Cols()), M->compute(i - orow, j - ocol), 0.0);
		}
	};

	template<typename Var>
	class MatrixUnit : public AbstractMatrixReadOnly< Var, typename Op1<Var, const_branch_expression>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		enumerator n;
		Var c;
	public:
		MatrixUnit(enumerator n, const Var & c = Var(1.0)) :n(n), c(c) {}
		MatrixUnit(const MatrixUnit& b) : n(b.n), c(b.c) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return n; }
		__INLINE typename Op1<Var, const_branch_expression>::type compute(enumerator i, enumerator j) const
		{
			return branch(i == j, c, 0.0);
		}
	};

	template<typename Var>
	class MatrixRow : public AbstractMatrixReadOnly< Var, typename OpRef<Var>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		enumerator n;
		Var c;
	public:
		MatrixRow(enumerator n, const Var& c = Var(1.0)) :n(n), c(c) {}
		MatrixRow(const MatrixRow& b) : n(b.n), c(b.c) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return 1; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return n; }
		__INLINE typename OpRef<Var>::type compute(enumerator i, enumerator j) const { return c; }
	};

	template<typename Var>
	class MatrixCol : public AbstractMatrixReadOnly< Var, typename OpRef<Var>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		enumerator n;
		Var c;
	public:
		MatrixCol(enumerator n, const Var& c = Var(1.0)) :n(n), c(c) {}
		MatrixCol(const MatrixCol& b) : n(b.n), c(b.c) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return 1; }
		__INLINE  typename OpRef<Var>::type compute(enumerator i, enumerator j) const { return c; }
	};

	template<typename Var>
	class MatrixDiag : public AbstractMatrixReadOnly< Var, typename Op1<Var, const_branch_expression>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		enumerator n;
		Var* diag;
	public:
		MatrixDiag(Var* diag, enumerator n) : diag(diag), n(n) {}
		MatrixDiag(const MatrixDiag& b) : diag(b.diag), n(b.n) {}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return n; }
		__INLINE typename Op1<Var, const_branch_expression>::type compute(enumerator i, enumerator j) const
		{
			return branch(i == j, diag[i], 0.0);
			//return i == j ? diag[i] : Var(0.0); 
		}
	};

	template<typename VarA, typename VarB, typename RetA, typename RetB>
	class MatrixSum : public AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type, typename OpAdd<RetA, RetB>::type>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<VarA, RetA>* A;
		const AbstractMatrixReadOnly<VarB, RetB>* B;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixSum(const AbstractMatrixReadOnly<VarA, RetA>& rA, const AbstractMatrixReadOnly<VarB, RetB>& rB)
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
		__INLINE typename OpAdd<RetA, RetB>::type compute(enumerator i, enumerator j) const
		{
			return add(A->compute(i, j), B->compute(i, j));
		}
	};
	
	template<typename VarA, typename VarB, typename RetA, typename RetB>
	class MatrixDifference : public AbstractMatrixReadOnly< typename Promote<VarA, VarB>::type, typename OpSub<RetA, RetB>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<VarA, RetA>* A;
		const AbstractMatrixReadOnly<VarB, RetB>* B;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixDifference(const AbstractMatrixReadOnly<VarA, RetA>& rA, const AbstractMatrixReadOnly<VarB, RetB>& rB)
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
		__INLINE typename OpSub<RetA, RetB>::type compute(enumerator i, enumerator j) const
		{
			return sub(A->compute(i, j), B->compute(i, j));
		}
	};

	template<typename Var, typename Ret>
	class AbstractMatrixTranspose : public AbstractMatrixReadOnly<Var, Ret>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<Var, Ret>* A;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Cols(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Rows(); }
		AbstractMatrixTranspose(const AbstractMatrixReadOnly<Var, Ret>& rA) : A(&rA) {}
		AbstractMatrixTranspose(const AbstractMatrixTranspose& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Ret compute(enumerator i, enumerator j) const { return A->compute(j, i); }
	};

	template<typename Var, typename Ret>
	class AbstractMatrixConjugateTranspose : public AbstractMatrixReadOnly<Var, Ret>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<Var,Ret>* A;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Cols(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Rows(); }
		AbstractMatrixConjugateTranspose(const AbstractMatrixReadOnly<Var, Ret>& rA)
			: A(&rA) {}
		AbstractMatrixConjugateTranspose(const AbstractMatrixConjugateTranspose& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Ret compute(enumerator i, enumerator j) const { return conj(A->compute(j, i)); }
	};
	
	template<typename Var, typename Ret>
	class AbstractMatrixConjugate : public AbstractMatrixReadOnly<Var,Ret>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<Var,Ret>* A;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		AbstractMatrixConjugate(const AbstractMatrixReadOnly<Var,Ret>& rA)
			: A(&rA) {}
		AbstractMatrixConjugate(const AbstractMatrixConjugate& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Ret compute(enumerator i, enumerator j) const { return conj(A->compute(i, j)); }
	};

	template<typename Var, typename Ret>
	class MatrixUnaryMinus : public AbstractMatrixReadOnly<Var, typename Op1<Ret, unary_minus_expression>::type>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<Var, Ret>* A;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixUnaryMinus(const AbstractMatrixReadOnly<Var, Ret>& rA)
			: A(&rA) {}
		MatrixUnaryMinus(const MatrixUnaryMinus& b) : A(b.A) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE typename Op1<Ret, unary_minus_expression>::type compute(enumerator i, enumerator j) const { return -A->compute(i, j); }
	};
	
	template<typename Var>
	class MatrixTranspose : public AbstractMatrix<Var>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		AbstractMatrix<Var>* A;
	public:
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
		__INLINE const Var & get(enumerator i, enumerator j) const { return A->get(j, i); }
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element.
		__INLINE Var& get(enumerator i, enumerator j) { return A->get(j, i); }
		/// This is a stub function to fulfill abstract
		/// inheritance. 
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
			if (Cols() != cols || Rows() != rows) throw Impossible;
		}
	};

	template<typename VarA, typename VarB, typename RetA, typename RetB>
	class AbstractMatrixConcatRows : public AbstractMatrixReadOnly<typename Promote<VarA,VarB>::type, typename OpCond<RetA, RetB>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<VarA, RetA>* A;
		const AbstractMatrixReadOnly<VarB, RetB>* B;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows() + B->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		AbstractMatrixConcatRows(const AbstractMatrixReadOnly<VarA, RetA>& rA, const AbstractMatrixReadOnly<VarB, RetB>& rB)
			: A(&rA), B(&rB) 
		{
			assert(A->Cols() == B->Cols());
		}
		AbstractMatrixConcatRows(const AbstractMatrixConcatRows& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE typename OpCond<RetA, RetB>::type compute(enumerator i, enumerator j) const
		{
			return branch(i < A->Rows(), A->compute(i, j), B->compute(i - A->Rows(), j));
		}
	};
	
	template<typename Var>
	class MatrixConcatRows : public AbstractMatrix<Var>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		AbstractMatrix<Var>* A;
		AbstractMatrix<Var>* B;
	public:
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
		/// @return Reference to element.
		__INLINE Var& get(enumerator i, enumerator j)  
		{ 
			return i < A->Rows() ? A->get(i, j) : B->get(i - A->Rows(), j);
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE const Var& get(enumerator i, enumerator j) const
		{
			return i < A->Rows() ? A->get(i, j) : B->get(i - A->Rows(), j);
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. 
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
			if (Cols() != cols || Rows() != rows) throw Impossible;
		}
	};

	
	template<typename VarA, typename VarB, typename RetA, typename RetB>
	class AbstractMatrixConcatCols : public AbstractMatrixReadOnly<typename Promote<VarA, VarB>::type, typename OpCond<RetA, RetB>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<VarA, RetA>* A;
		const AbstractMatrixReadOnly<VarB, RetB>* B;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols() + B->Cols(); }
		AbstractMatrixConcatCols(const AbstractMatrixReadOnly<VarA, RetA>& rA, const AbstractMatrixReadOnly<VarB, RetB>& rB)
			: A(&rA), B(&rB) {
			assert(A->Rows() == B->Rows());
		}
		AbstractMatrixConcatCols(const AbstractMatrixConcatCols& b) : A(b.A), B(b.B) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE typename OpCond<RetA, RetB>::type compute(enumerator i, enumerator j) const
		{
			return branch(j < A->Cols(), A->compute(i, j), B->compute(i, j - A->Cols()));
		}
	};
	
	template<typename Var>
	class MatrixConcatCols : public AbstractMatrix<Var>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		AbstractMatrix<Var>* A;
		AbstractMatrix<Var>* B;
	public:
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
		/// @return Reference to element.
		__INLINE Var& get(enumerator i, enumerator j) 
		{
			return j < A->Cols() ? A->get(i, j) : B->get(i, j - A->Cols());
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE const Var& get(enumerator i, enumerator j) const
		{
			return j < A->Cols() ? A->get(i, j) : B->get(i, j - A->Cols());
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. 
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
			if (Cols() != cols || Rows() != rows) throw Impossible;
		}
	};

	template<typename Var, typename Ret>
	class AbstractMatrixRepack : public AbstractMatrixReadOnly<Var, Ret>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<Var, Ret>* A;
		enumerator n, m;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return n; }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return m; }
		AbstractMatrixRepack(const AbstractMatrixReadOnly<Var, Ret>& rA, enumerator n, enumerator m)
			: A(&rA), n(n), m(m) { assert(A->Rows() * A->Cols() == n * m); }
		AbstractMatrixRepack(const AbstractMatrixRepack& b) : A(b.A), n(b.n), m(b.m) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE Ret compute(enumerator i, enumerator j) const
		{
			enumerator ind = i * m + j;
			return A->compute(ind / A->Cols(), ind % A->Cols());
		}
	};

	template<typename Var>
	class MatrixRepack : public AbstractMatrix<Var>
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		AbstractMatrix<Var>* A;
		enumerator n, m;
	public:
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
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element.
		__INLINE Var & get(enumerator i, enumerator j) 
		{ 
			enumerator p = i * m + j;
			return A->get(p / A->Cols(), p % A->Cols());
			//return A->data()[i * m + j]; 
		}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE const Var& get(enumerator i, enumerator j) const
		{
			enumerator p = i * m + j;
			return A->get(p / A->Cols(), p % A->Cols());
			//return A->data()[i * m + j]; 
		}
		/// This is a stub function to fulfill abstract
		/// inheritance. 
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
			if (Cols() != cols || Rows() != rows) throw Impossible;
		}
	};

	template<typename VarA, typename VarB, typename RetA, typename RetB>
	class MatrixMul : public AbstractMatrix< typename Promote<VarA, VarB>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		Matrix<typename Promote<VarA, VarB>::type> M;
	public:
		/// This is a stub function to fulfill abstract
		/// inheritance. 
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
			if (Cols() != cols || Rows() != rows) throw Impossible;
		}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return M.Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return M.Cols(); }
		MatrixMul(const AbstractMatrixReadOnly<VarA,RetA>& rA, const AbstractMatrixReadOnly<VarB,RetB>& rB)
		{
			assert(rA.Cols() == rB.Rows());
			const AbstractMatrix<VarA>* pA = dynamic_cast<const AbstractMatrix<VarA> *>(&rA);
			const AbstractMatrix<VarB>* pB = dynamic_cast<const AbstractMatrix<VarB> *>(&rB);
			if (!pA)
			{
				static thread_private< Matrix<VarA> > tmpA;
				*tmpA = rA;
				pA = &(*tmpA);
			}
			if (!pB)
			{
				static thread_private< Matrix<VarB> > tmpB;
				*tmpB = rB;
				pB = &(*tmpB);
			}
			M.Resize(rA.Rows(), rB.Cols());
			enumerator Arows = pA->Rows(), Acols = pA->Cols(), Bcols = pB->Cols(), Brows = pB->Rows();
			for (enumerator i = 0; i < Arows; ++i)
				for (enumerator j = 0; j < Bcols; ++j)
					M(i, j) = (*pA)(i, i + 1, 0, Acols).DotProduct((*pB)(0, Brows, j, j + 1).Transpose());
		}
		MatrixMul(const MatrixMul& b) : M(b.M) {}
		/// Access element of the matrix by row and column indices.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element.
		__INLINE typename Promote<VarA, VarB>::type& get(enumerator i, enumerator j) { return M.get(i, j); }
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE const typename Promote<VarA, VarB>::type& get(enumerator i, enumerator j) const { return M.get(i, j); }
	};

	template<typename VarA, typename VarB, typename RetA, typename VarAB>
	class MatrixMulCoef : public AbstractMatrixReadOnly< VarAB, typename OpMul<RetA, VarB>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<VarA, RetA>* A;
		VarB coef;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixMulCoef(const AbstractMatrixReadOnly<VarA, RetA>& rA, const VarB& rcoef)
			: A(&rA), coef(rcoef) {}
		MatrixMulCoef(const MatrixMulCoef& b) : A(b.A), coef(b.coef) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE typename OpMul<RetA, VarB>::type compute(enumerator i, enumerator j) const { return mul(A->compute(i, j), coef); }
	};
	
	template<typename VarA, typename VarB, typename RetA, typename VarAB>
	class MatrixDivCoef : public AbstractMatrixReadOnly< VarAB, typename OpDiv<RetA, VarB>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		const AbstractMatrixReadOnly<VarA, RetA>* A;
		VarB coef;
	public:
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return A->Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return A->Cols(); }
		MatrixDivCoef(const AbstractMatrixReadOnly<VarA, RetA>& rA, const VarB& rcoef)
			: A(&rA), coef(rcoef) {}
		MatrixDivCoef(const MatrixDivCoef& b)
			: A(b.A), coef(b.coef) {}
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE typename OpDiv<RetA, VarB>::type compute(enumerator i, enumerator j) const { return div(A->compute(i, j), coef); }
	};

	template<typename VarA, typename VarB, typename RetA, typename RetB>
	class KroneckerProduct : public AbstractMatrix< typename Promote<VarA, VarB>::type >
	{
	public:
		typedef typename AbstractMatrixBase::enumerator enumerator;
	private:
		Matrix<typename Promote<VarA, VarB>::type> M;
	public:
		/// This is a stub function to fulfill abstract
		/// inheritance. 
		/// since it just points to a part of the larger empty matrix.
		void Resize(enumerator rows, enumerator cols)
		{
			assert(Cols() == cols);
			assert(Rows() == rows);
			(void)cols; (void)rows;
			if (Cols() != cols || Rows() != rows) throw Impossible;
		}
		/// Number of rows.
		/// @return Number of rows.
		__INLINE enumerator Rows() const { return M.Rows(); }
		/// Number of columns.
		/// @return Number of columns.
		__INLINE enumerator Cols() const { return M.Cols(); }
		KroneckerProduct(const AbstractMatrixReadOnly<VarA, RetA>& rA, const AbstractMatrixReadOnly<VarB, RetB>& rB)
		{
			const AbstractMatrix<VarA>* pA = dynamic_cast<const AbstractMatrix<VarA> *>(&rA);
			const AbstractMatrix<VarB>* pB = dynamic_cast<const AbstractMatrix<VarB> *>(&rB);
			if (!pA)
			{
				static thread_private< Matrix<VarA> > tmpA;
				*tmpA = rA;
				pA = &(*tmpA);
			}
			if (!pB)
			{
				static thread_private< Matrix<VarB> > tmpB;
				*tmpB = rB;
				pB = &(*tmpB);
			}
			M.Resize(rA.Rows() * rB.Rows(), rA.Cols() * rB.Cols());
			enumerator Mrows = M.Rows(), Mcols = M.Cols();
			enumerator Brows = pB->Rows(), Bcols = pB->Cols();
			for (enumerator i = 0; i < Mrows; ++i)
				for (enumerator j = 0; j < Mcols; ++j)
					M(i, j) = pA->get(i / Brows, j / Bcols) * pB->get(i % Brows, j % Bcols);
		}
		KroneckerProduct(const KroneckerProduct& b) : M(b.M) {}
		/// Access element of the matrix by row and column indices.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to element.
		__INLINE typename Promote<VarA, VarB>::type& get(enumerator i, enumerator j) { return M.get(i, j); }
		/// Access element of the matrix by row and column indices
		/// without right to change the element.
		/// @param i Row index.
		/// @param j Column index.
		/// @return Reference to constant element.
		__INLINE const typename Promote<VarA, VarB>::type& get(enumerator i, enumerator j) const { return M.get(i, j); }
	};
	
	template<typename Var>
	void AbstractMatrix<Var>::Swap(AbstractMatrix<Var> & b)
	{
		Matrix<Var> tmp = b;
		b = (*this);
		(*this) = tmp;
	}
	
	template<typename Var, typename RetType>
	template<typename typeB, typename retB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var, RetType>::CrossProduct(const AbstractMatrixReadOnly<typeB,retB> & B) const
	{
		const AbstractMatrixReadOnly<Var,RetType> & A = *this;
		assert(A.Rows() == 3);
		assert(A.Cols() == 1);
		assert(B.Rows() == 3);
		Matrix<typename Promote<Var,typeB>::type> ret(3,B.Cols()); //check RVO
		for(enumerator k = 0; k < B.Cols(); ++k)
		{
			ret(0,k) = (A.compute(1,0)*B.compute(2,k) - A.compute(2,0)*B.compute(1,k));
			ret(1,k) = (A.compute(2,0)*B.compute(0,k) - A.compute(0,0)*B.compute(2,k));
			ret(2,k) = (A.compute(0,0)*B.compute(1,k) - A.compute(1,0)*B.compute(0,k));
		}
		return ret;
	}
		
	template<typename Var, typename RetType>
	template<typename typeB, typename retB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var, RetType>::Transform(const AbstractMatrixReadOnly<typeB, retB> & B) const
	{
		const AbstractMatrixReadOnly<Var, RetType> & A = *this;
		assert(A.Rows() == 3);
		assert(A.Cols() == 1);
		assert(B.Rows() == 3);
		assert(B.Cols() == 1);
		typedef typename Promote<Var,typeB>::type var_t;
		Matrix<var_t> Q(3,3);
		var_t x,y,z,w,n,l;
		
		x = -(B.compute(1,0)*A.compute(2,0) - B.compute(2,0)*A.compute(1,0));
		y = -(B.compute(2,0)*A.compute(0,0) - B.compute(0,0)*A.compute(2,0));
		z = -(B.compute(0,0)*A.compute(1,0) - B.compute(1,0)*A.compute(0,0));
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
	
	template<typename Var, typename RetType>
	template<typename typeB, typename retB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var,RetType>::CholeskySolve(const AbstractMatrixReadOnly<typeB,retB> & B, int * ierr) const
	{
		const AbstractMatrixReadOnly<Var, RetType> & A = *this;
		assert(A.Rows() == A.Cols());
		assert(A.Rows() == B.Rows());
		enumerator n = A.Rows();
		enumerator l = B.Cols();
		Matrix<typename Promote<Var,typeB>::type> ret(B);
		static thread_private< SymmetricMatrix<Var> > L;// (A);
		*L = A;
		
		//SAXPY
#if 1
		for(enumerator i = 0; i < n; ++i)
		{
			for (enumerator j = i; j < n; ++j)
				(*L)(i, j) -= (*L)(i, i + 1, 0, i).DotProduct((*L)(j, j + 1, 0, i));
				//for (enumerator k = 0; k < i; ++k)
				//	(*L)(i, j) -= (*L)(i, k) * (*L)(j, k);
			if (real_part((*L)(i, i)) < 0.0)
			{
				if (ierr)
				{
					if (*ierr == -1) std::cout << "Negative diagonal pivot " << get_value((*L)(i, i)) << " row " << i << std::endl;
					*ierr = i + 1;
				}
				else throw MatrixCholeskySolveFail;
				return ret;
			}
			(*L)(i, i) = sqrt((*L)(i, i));
			if (fabs((*L)(i, i)) < 1.0e-24)
			{
				if (ierr)
				{
					if (*ierr == -1) std::cout << "Diagonal pivot is too small " << get_value((*L)(i, i)) << " row " << i << std::endl;
					*ierr = i + 1;
				}
				else throw MatrixCholeskySolveFail;
				return ret;
			}
			for (enumerator j = i + 1; j < n; ++j)
				(*L)(i, j) = (*L)(i, j) / (*L)(i, i);
		}
#else
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
#endif
		// LY=B
		Matrix<typename Promote<Var,typeB>::type> & Y = ret;
		for(enumerator i = 0; i < n; ++i)
		{
			for(enumerator k = 0; k < l; ++k)
			{
				//for (enumerator j = 0; j < i; ++j)
				//	Y(i, k) -= Y(j, k) * (*L)(j, i);
				Y(i, k) -= Y(0, i, k, k + 1).DotProduct((*L)(0, i, i, i + 1));
				Y(i, k) /= (*L)(i, i);
			}
		}
		// L^TX = Y
		Matrix<typename Promote<Var,typeB>::type> & X = ret;
		for(enumerator it = n; it > 0; --it)
		{
			enumerator i = it-1;
			for(enumerator k = 0; k < l; ++k)
			{
				//for(enumerator jt = n; jt > it; --jt)
				//{
				//	enumerator j = jt-1;
				//	X(i,k) -= X(j,k)*(*L)(i,j);
				//}
				X(i, k) -= X(i + 1, n, k, k + 1).DotProduct((*L)(i, i + 1, i + 1, n).Transpose());
				X(i, k) /= (*L)(i, i);
			}
		}
		if( ierr ) *ierr = 0;
		return ret;
	}

	template<typename Var, typename RetType>
	template<typename typeB, typename retB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var, RetType>::Solve(const AbstractMatrixReadOnly<typeB, retB> & B, int * ierr) const
	{
		// for A^T B
		assert(Rows() == B.Rows());
		
		if( Rows() != Cols() )
		{
			AbstractMatrixTranspose<Var,RetType> At = this->Transpose(); //m by n matrix
			return (At * (*this)).Solve(At * B, ierr);
		}
		Matrix<typename Promote<Var,typeB>::type> AtB = B; //m by l matrix
		Matrix<Var> AtA = (*this); //m by m matrix
		enumerator l = B.Cols();
		enumerator m = Cols();
		Matrix<typename Promote<Var,typeB>::type> ret(m,l);
		//ConstMatrixTranspose<Var> At = this->Transpose(); //m by n matrix
		//Matrix<typename Promote<Var,typeB>::type> AtB = At*B; //m by l matrix
		//Matrix<Var> AtA = At*(*this); //m by m matrix
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
	
	template<typename Var, typename RetType>
	Matrix<Var>
	AbstractMatrixReadOnly<Var, RetType>::PseudoInvert(INMOST_DATA_REAL_TYPE tol, int * ierr) const
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
		ret = V * S * U.Transpose();
		//ret = V.ConjugateTranspose().Transpose() * S * U.ConjugateTranspose();
		if( ierr ) *ierr = 0;
		return ret;
	}

	template<typename Var, typename RetType>
	Matrix<Var>
	AbstractMatrixReadOnly<Var, RetType>::Root(INMOST_DATA_ENUM_TYPE iter, INMOST_DATA_REAL_TYPE tol, int *ierr) const
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
	template<typename Var, typename RetType>
	Matrix<Var>
	AbstractMatrixReadOnly<Var, RetType>::Power(INMOST_DATA_REAL_TYPE n, int * ierr) const
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
	template<typename Var, typename RetType>
	template<typename typeB, typename retB>
	Matrix<typename Promote<Var,typeB>::type>
	AbstractMatrixReadOnly<Var, RetType>::PseudoSolve(const AbstractMatrixReadOnly<typeB, retB> & B, INMOST_DATA_REAL_TYPE tol, int * ierr) const
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
	void AbstractMatrix<Var>::MPT(INMOST_DATA_ENUM_TYPE * Perm, INMOST_DATA_REAL_TYPE * SL, INMOST_DATA_REAL_TYPE * SR) const
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
	
	template<typename Var, typename RetType>
	bool AbstractMatrixReadOnly<Var, RetType>::cSVD(AbstractMatrix<Var>& U, AbstractMatrix<Var>& Sigma, AbstractMatrix<Var>& V) const
	{
		int m = Rows();
		int n = Cols();
		//data eta / 1.1920929e-07 /
		//data tol / 1.5e-31 /
		if (m > n)
		{
			bool success = ConjugateTranspose().cSVD(V, Sigma, U);
			if (success)
				return true;
			else return false;
		}
		Var cs, eta, f, g, h, q, r, sn, w, x, y, z;
		int i, j, k, l, p = 0;
		std::vector<Var> t, b, c, s;
		Matrix<Var> A = *this;
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

	template<typename Var, typename RetType>
	bool AbstractMatrixReadOnly<Var, RetType>::SVD(AbstractMatrix<Var>& U, AbstractMatrix<Var>& Sigma, AbstractMatrix<Var>& V, bool order_singular_values, bool nonnegative) const
	{
		int flag, i, its, j, jj, k, l, nm;
		int n = Rows();
		int m = Cols();
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
		Var c, f, h, s, x, y, z;
		Var g = 0.0, scale = 0.0;
		INMOST_DATA_REAL_TYPE anorm = 0.0;
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


	template<typename VarA, typename RetA, typename VarB, typename RetB>
	__INLINE typename Promote<VarA, VarB>::type DotProduct(const AbstractMatrixReadOnly<VarA, RetA>& A, const AbstractMatrixReadOnly<VarB, RetB>& B)
	{
		assert(A.Cols() == B.Cols() && A.Rows() == B.Rows());
		typename Promote<VarA, VarB>::type ret = 0.0;
		unsigned Arows = A.Rows(), Acols = A.Cols();
		for (unsigned i = 0; i < Arows; ++i)
			for (unsigned j = 0; j < Acols; ++j)
				ret += A.compute(i, j) * B.compute(i, j);
		return ret;
	}


#if defined(USE_AUTODIFF)
	template<typename RetB>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>& A, const AbstractMatrixReadOnly<variable, RetB>& B, Sparse::RowMerger & m)
	{
		assert(A.Cols() == B.Cols() && A.Rows() == B.Rows());
		variable ret = 0.0;
		INMOST_DATA_REAL_TYPE value = 0.0;
		unsigned Arows = A.Rows(), Acols = A.Cols();
		m.Clear();
		for (unsigned i = 0; i < Arows; ++i)
			for (unsigned j = 0; j < Acols; ++j)
			{
				value += A.compute(i, j) * B.compute(i, j).GetValue();
				B.compute(i, j).GetJacobian(A.compute(i, j), m);
			}
		ret.SetValue(value);
		m.RetrieveRow(ret.GetRow());
		return ret;
	}

	template<typename RetA, typename RetB>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<variable, RetA>& A, const AbstractMatrixReadOnly<variable, RetB>& B, Sparse::RowMerger & m)
	{
		assert(A.Cols() == B.Cols() && A.Rows() == B.Rows());
		variable ret = 0.0;
		INMOST_DATA_REAL_TYPE value = 0.0;
		unsigned Arows = A.Rows(), Acols = A.Cols();
		m.Clear();
		for (unsigned i = 0; i < Arows; ++i)
			for (unsigned j = 0; j < Acols; ++j)
			{
				value += A.compute(i, j).GetValue() * B.compute(i, j).GetValue();
				A.compute(i, j).GetJacobian(B.compute(i, j).GetValue(), m);
				B.compute(i, j).GetJacobian(A.compute(i, j).GetValue(), m);
			}
		ret.SetValue(value);
		m.RetrieveRow(ret.GetRow());
		return ret;
	}
	

	template<typename RetB>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>& A, const AbstractMatrixReadOnly<variable, RetB>& B)
	{
		return DotProduct(A, B, AbstractMatrixBase::GetMerger());
	}

	template<typename RetA>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<variable, RetA>& A, const AbstractMatrixReadOnly<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE>& B)
	{
		return DotProduct(B, A, AbstractMatrixBase::GetMerger());
	}

	template<typename RetA, typename RetB>
	__INLINE variable DotProduct(const AbstractMatrixReadOnly<variable, RetA>& A, const AbstractMatrixReadOnly<variable, RetB>& B)
	{
		return DotProduct(A, B, AbstractMatrixBase::GetMerger());
	}


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
#endif //USE_AUTODIFF
}
/// Multiplication of matrix by constant from left.
/// @param coef Constant coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a constant.
template<typename typeA, typename typeB, typename retB>
INMOST::MatrixMulCoef<typeB, typeA, retB, typename INMOST::Promote<typeB, typeA>::type>
operator *(const typeA & coef, const INMOST::AbstractMatrixReadOnly<typeB, retB>& other)
{return INMOST::MatrixMulCoef<typeB, typeA, retB, typename INMOST::Promote<typeB, typeA>::type>(other, coef);}

#if defined(USE_AUTODIFF)
/// Multiplication of matrix by constant from left.
/// @param coef Constant coefficient multiplying matrix.
/// @param other Matrix to be multiplied.
/// @return Matrix, each entry multiplied by a constant.
template<class A, typename typeB, typename retB>
INMOST::MatrixMulCoef<typeB, INMOST::variable, retB, typename INMOST::Promote<typeB, INMOST::variable>::type>
operator *(INMOST::shell_expression<A> const& coef, const INMOST::AbstractMatrixReadOnly<typeB, retB>& other)
{return INMOST::MatrixMulCoef<typeB, INMOST::variable, retB, typename INMOST::Promote<typeB, INMOST::variable>::type>(other, INMOST::variable(coef));}
#endif //USE_AUTODIFF

template<typename T, typename R>
__INLINE bool check_nans(const INMOST::AbstractMatrixReadOnly<T,R> & A) {return A.CheckNans();}
template<typename T, typename R>
__INLINE bool check_infs(const INMOST::AbstractMatrixReadOnly<T,R> & A) {return A.CheckInfs();}
template<typename T, typename R>
__INLINE bool check_nans_infs(const INMOST::AbstractMatrixReadOnly<T,R> & A) {return A.CheckNansInfs();}

#endif //INMOST_DENSE_INCLUDED
