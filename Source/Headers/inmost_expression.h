#ifndef INMOST_AUTODIFF_ETEXPR_H_INCLUDED
#define INMOST_AUTODIFF_ETEXPR_H_INCLUDED
#include "inmost_common.h"
#include "inmost_sparse.h"
#include <sstream> //for debug
#include <new>

#if defined(USE_AUTODIFF) && !defined(USE_SOLVER)
#warning "USE_AUTODIFF require USE_SOLVER"
#undef USE_AUTODIFF
#endif


//TODO:
// 1. Incorporate tables
// 2. (ok, test) implement condition
// 3. (ok, test) implement stencil
// 4. (???) copying of basic_dynamic_variable
// 5. Consider optimization by checking zero variation multipliers, check that assembly do not degrade.
// 6. floor, ceil, atan, acos, asin, max, min functions
// 7. choice of directional derivatives at discontinuities for abs, pow, max, min (see ADOL-C)


#if defined(USE_AUTODIFF)
namespace INMOST
{





	class basic_expression
	{
	public:
    basic_expression() {}//std::cout << this << " Created" << std::endl;}
    basic_expression(const basic_expression & other) {};//std::cout << this << " Created from " << &other << std::endl;}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue() const = 0;
		__INLINE virtual void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const = 0;
    __INLINE virtual void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const = 0;
    //virtual ~basic_expression() {std::cout << this << " Destroied" << std::endl;}
	};

	template<class Derived>
	class shell_expression : virtual public basic_expression
	{
	public:
		shell_expression() {}//std::cout << this << " Shell Created for " << dynamic_cast<basic_expression *>(this) << std::endl;}
		shell_expression(const shell_expression & other) {}//std::cout << this << " Shell Created from " << &other << std::endl;}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue() const {return static_cast<const Derived *>(this)->GetValue(); }
    __INLINE virtual void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const { return static_cast<const Derived *>(this)->GetDerivative(mult,r); }
    __INLINE virtual void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const { return static_cast<const Derived *>(this)->GetDerivative(mult,r); }
    operator Derived & () {return *static_cast<Derived *>(this);}
    operator const Derived & () const {return *static_cast<const Derived *>(this);}
    //~shell_expression() {std::cout << this << " Shell Destroied for " << dynamic_cast<basic_expression *>(this) << std::endl;}
    //Derived * GetDerived() { return dynamic_cast<Derived *>(this); }
	};

  
 
  
  

  class var_expression : public shell_expression<var_expression>
  {
    INMOST_DATA_REAL_TYPE value;
    INMOST_DATA_ENUM_TYPE index;
  public:
    var_expression(const var_expression & other) :value(other.value), index(other.index) {}
    var_expression(INMOST_DATA_REAL_TYPE pvalue, INMOST_DATA_ENUM_TYPE pindex) : value(pvalue), index(pindex) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {if( index != ENUMUNDEF ) r[index] += mult;}
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {if( index != ENUMUNDEF ) r[index] += mult;}
    __INLINE var_expression & operator =(var_expression const & other)
    {
      value = other.value;
      index = other.index;
      return *this;
    }
    bool check_nans() const
    {
      return value != value;
    }
  };


  class multivar_expression : public shell_expression<multivar_expression>
  {
    INMOST_DATA_REAL_TYPE value;
    Sparse::Row entries;
  public:
    multivar_expression() :value(0) {}
    multivar_expression(INMOST_DATA_REAL_TYPE pvalue) : value(pvalue) {}
    multivar_expression(const multivar_expression & other) : value(other.value), entries(other.entries) {}
    multivar_expression(INMOST_DATA_REAL_TYPE pvalue, Sparse::Row & pentries)
     : value(pvalue), entries(pentries) {}
    multivar_expression(INMOST_DATA_REAL_TYPE pvalue, INMOST_DATA_ENUM_TYPE pindex, INMOST_DATA_REAL_TYPE pdmult = 1.0)
     : value(pvalue)
    {
      entries.Push(pindex,pdmult);
    }
    multivar_expression(const basic_expression & expr)
    {
      expr.GetDerivative(1.0,entries);
      value = expr.GetValue();
    }
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
        r[it->first] += it->second*mult;
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
        r[it->first] += it->second*mult;
    }
    __INLINE multivar_expression & operator = (INMOST_DATA_REAL_TYPE pvalue)
    {
      value = pvalue;
      entries.Clear();
      return *this;
    }
    __INLINE multivar_expression & operator = (basic_expression const & expr)
    {
      value = expr.GetValue();
      Sparse::Row tmp;
      expr.GetDerivative(1.0,tmp);
      entries.Swap(tmp);
      return *this;
    }
    __INLINE multivar_expression & operator = (multivar_expression const & other)
    {
      value = other.value;
      entries = other.entries;
      return *this;
    }
    __INLINE Sparse::Row & GetRow() {return entries;}
    __INLINE const Sparse::Row & GetRow() const {return entries;}
    __INLINE multivar_expression & operator +=(basic_expression const & expr)
    {
      value += expr.GetValue();
      Sparse::Row tmp(entries);
      expr.GetDerivative(1.0,tmp);
      entries.Swap(tmp);
      return *this;
    }
    __INLINE multivar_expression & operator -=(basic_expression const & expr)
    {
      value -= expr.GetValue();
      Sparse::Row tmp(entries);
      expr.GetDerivative(-1.0,tmp);
      entries.Swap(tmp);
      return *this;
    }
    __INLINE multivar_expression & operator *=(basic_expression const & expr)
    {
      INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
      Sparse::Row tmp(entries);
      for(Sparse::Row::iterator it = tmp.Begin(); it != tmp.End(); ++it) it->second *= rval;
      expr.GetDerivative(lval,tmp);
      entries.Swap(tmp);
      value *= rval;
      return *this;
    }
    __INLINE multivar_expression & operator /=(basic_expression const & expr)
    {
      INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
      INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0/rval;
      value *= reciprocial_rval;
      Sparse::Row tmp(entries);
      for(Sparse::Row::iterator it = tmp.Begin(); it != tmp.End(); ++it) it->second *= reciprocial_rval;
      expr.GetDerivative(-value*reciprocial_rval,tmp); 
      entries.Swap(tmp);
      return *this;
    }
    __INLINE multivar_expression & operator +=(INMOST_DATA_REAL_TYPE right)
    {
      value += right;
      return *this;
    }
    __INLINE multivar_expression & operator -=(INMOST_DATA_REAL_TYPE right)
    {
      value -= right;
      return *this;
    }
    __INLINE multivar_expression & operator *=(INMOST_DATA_REAL_TYPE right)
    {
      value *= right;
      for(Sparse::Row::iterator it = entries.Begin(); it != entries.End(); ++it) it->second *= right;
      return *this;
    }
    __INLINE multivar_expression & operator /=(INMOST_DATA_REAL_TYPE right)
    {
      value /= right;
      for(Sparse::Row::iterator it = entries.Begin(); it != entries.End(); ++it) it->second /= right;
      return *this;
    }
    bool check_nans() const
    {
      if( value != value ) return true;
      for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
        if( it->second != it->second ) return true;
      return false;
    }
  };

  template<class A>
  class const_multiplication_expression : public shell_expression<const_multiplication_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
  public:
    const_multiplication_expression(const shell_expression<A> & parg,INMOST_DATA_REAL_TYPE pdmult) : arg(parg), dmult(pdmult)
    {
      value = arg.GetValue()*dmult;
    }
    const_multiplication_expression(const const_multiplication_expression & other) : arg(other.arg), value(other.value), dmult(other.dmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetDerivative(mult*dmult,r);
    }
  };

  template<class A>
  class variation_multiplication_expression : public shell_expression<variation_multiplication_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
  public:
    variation_multiplication_expression(const shell_expression<A> & parg,INMOST_DATA_REAL_TYPE pdmult) : arg(parg), dmult(pdmult)
    {
      value = arg.GetValue();
    }
    variation_multiplication_expression(const variation_multiplication_expression & other) : arg(other.arg), value(other.value), dmult(other.dmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetDerivative(mult*dmult,r);
    }
  };


  template<class A>
  class const_division_expression : public shell_expression<const_division_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
  public:
    const_division_expression(const shell_expression<A> & parg,INMOST_DATA_REAL_TYPE pdmult) : arg(parg), dmult(pdmult)
    {
      dmult = 1.0/dmult;
      value = arg.GetValue()*dmult;
    }
    const_division_expression(const const_division_expression & other) : arg(other.arg), value(other.value), dmult(other.dmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetDerivative(mult*dmult,r);
    }
  };
  
  template<class A>
  class const_addition_expression : public shell_expression<const_addition_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value;
  public:
    const_addition_expression(const shell_expression<A> & parg,INMOST_DATA_REAL_TYPE padd) : arg(parg)
    {
      value = arg.GetValue()+padd;
    }
    const_addition_expression(const const_addition_expression & other) : arg(other.arg), value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetDerivative(mult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetDerivative(mult,r);
    }
  };

  template<class A>
  class const_subtraction_expression : public shell_expression<const_subtraction_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value;
  public:
    const_subtraction_expression(const shell_expression<A> & parg,INMOST_DATA_REAL_TYPE pleft) : arg(parg)
    {
      value = pleft-arg.GetValue();
    }
    const_subtraction_expression(const const_subtraction_expression & other) : arg(other.arg), value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetDerivative(-mult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetDerivative(-mult,r);
    }
  };

  template<class A>
  class reciprocal_expression : public shell_expression<reciprocal_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
  public:
    reciprocal_expression(const shell_expression<A> & parg,INMOST_DATA_REAL_TYPE pdmult) : arg(parg), dmult(pdmult)
    {
      value = dmult/arg.GetValue();
    }
    reciprocal_expression(const reciprocal_expression & other) : arg(other.arg), value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetDerivative(-mult*dmult*value*value,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetDerivative(-mult*dmult*value*value,r);
    }
  };

  template<class A>
	class unary_minus_expression : public shell_expression<unary_minus_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value;
	public:
    unary_minus_expression(const shell_expression<A> & parg) : arg(parg) {value = -arg.GetValue();}
    unary_minus_expression(const unary_minus_expression & b) : arg(b.arg) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative(-mult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative(-mult,r);
    }
	};

  template<class A>
  class abs_expression : public shell_expression<abs_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
	public:
    abs_expression(const shell_expression<A> & parg) : arg(parg) 
    {
      value = arg.GetValue();
      dmult = value < 0.0 ? -1.0 : 1.0;
      value = ::fabs(value);
    }
    abs_expression(const abs_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative( mult * dmult, r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative( mult * dmult, r);
    }
	};

  template<class A>
	class exp_expression : public shell_expression<exp_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value;
	public:
    exp_expression(const shell_expression<A> & parg) : arg(parg) 
    {
      value = arg.GetValue();
      value = ::exp(value);
    }
		exp_expression(const exp_expression & b) : unary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative( mult * value, r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative( mult * value, r);
    }
  };

  template<class A>
	class log_expression : public shell_expression<log_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
	public:
    log_expression(const shell_expression<A> & parg) : arg(parg) 
    {
      value = arg.GetValue();
      dmult = 1.0/value;
      value = ::log(value);
    }
    log_expression(const log_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
	};


  template<class A>
	class sin_expression : public shell_expression<sin_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
	public:
    sin_expression(const shell_expression<A> & parg) : arg(parg) 
    {
      value = arg.GetValue();
      dmult = ::cos(value);
      value = ::sin(value);
    }
    sin_expression(const sin_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
	};

  template<class A>
	class cos_expression : public shell_expression<cos_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
	public:
    cos_expression(const shell_expression<A> & parg) : arg(parg) 
    {
      value = arg.GetValue();
      dmult = -(::sin(value));
      value = ::cos(value);
    }
    cos_expression(const cos_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
	};

  template<class A>
	class sqrt_expression : public shell_expression<sqrt_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value;
	public:
    sqrt_expression(const shell_expression<A> & parg) : arg(parg) 
    {
      value = ::sqrt(arg.GetValue());
    }
    sqrt_expression(const sqrt_expression & b) : arg(b.arg), value(b.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative(0.5*mult/value,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative(0.5*mult/value,r);
    }
	};


  template<class A>
	class soft_abs_expression : public shell_expression<soft_abs_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
	public:
    soft_abs_expression(const shell_expression<A> & parg, INMOST_DATA_REAL_TYPE tol) : arg(parg) 
    {
      INMOST_DATA_REAL_TYPE lval = arg.GetValue();
      value = ::sqrt(lval*lval+tol*tol);
      dmult = lval/value;
    }
    soft_abs_expression(const soft_abs_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
	};

  template<class A>
	class soft_sign_expression : public shell_expression<soft_sign_expression<A> >
	{
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult;
	public:
    soft_sign_expression(const shell_expression<A> & parg, INMOST_DATA_REAL_TYPE tol) : arg(parg) 
    {
      INMOST_DATA_REAL_TYPE lval = arg.GetValue(), lval2 = lval*lval;
      INMOST_DATA_REAL_TYPE div = lval2+tol*tol;
      INMOST_DATA_REAL_TYPE sdiv = ::sqrt(div);
      value = lval/sdiv;
      dmult = (1.0 - lval2/div)/sdiv;
    }
    soft_sign_expression(const soft_sign_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetDerivative(mult*dmult,r);
    }
	};

  template<class A, class B>
  class soft_max_expression : public shell_expression<soft_max_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value, ldmult, rdmult;
  public:
    soft_max_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright, INMOST_DATA_REAL_TYPE tol) : left(pleft), right(pright)
    {
      INMOST_DATA_REAL_TYPE lval = left.GetValue(), rval = right.GetValue();
      INMOST_DATA_REAL_TYPE diff = lval-rval, root = ::sqrt(diff*diff+tol*tol);
      value = 0.5*(lval+rval + root);
      ldmult = 0.5*(1 + diff/root);
      rdmult = 0.5*(1 - diff/root);
    }
    soft_max_expression(const soft_max_expression & other) 
      : left(other.left), right(other.right), value(other.value), ldmult(other.ldmult), rdmult(other.rdmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
  };

  template<class A, class B>
  class soft_min_expression : public shell_expression<soft_min_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value, ldmult, rdmult;
  public:
    soft_min_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright, INMOST_DATA_REAL_TYPE tol) : left(pleft), right(pright)
    {
      INMOST_DATA_REAL_TYPE lval = left.GetValue(), rval = right.GetValue();
      INMOST_DATA_REAL_TYPE diff = lval-rval, root = ::sqrt(diff*diff+tol*tol);
      value = 0.5*(lval+rval - root);
      ldmult = 0.5*(1 - diff/root);
      rdmult = 0.5*(1 + diff/root);
    }
    soft_min_expression(const soft_min_expression & other) 
      : left(other.left), right(other.right), value(other.value), ldmult(other.ldmult), rdmult(other.rdmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
  };


  template<class A, class B>
  class multiplication_expression : public shell_expression<multiplication_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value, ldmult, rdmult;
  public:
    multiplication_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright)
    {
      rdmult = left.GetValue();
      ldmult = right.GetValue();
      value = rdmult*ldmult;
    }
    multiplication_expression(const multiplication_expression & other) 
      : left(other.left), right(other.right), value(other.value), ldmult(other.ldmult), rdmult(other.rdmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
  };

  template<class A, class B>
  class division_expression : public shell_expression<division_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value, reciprocal_rval;
  public:
    division_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright)
    {
      INMOST_DATA_REAL_TYPE lval = left.GetValue();
      INMOST_DATA_REAL_TYPE rval = right.GetValue();
      reciprocal_rval = 1.0 / rval;
      value = lval * reciprocal_rval;
    }
    division_expression(const division_expression & other) : left(other.left), right(other.right), value(other.value), reciprocal_rval(other.reciprocal_rval) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult * reciprocal_rval,r);
      right.GetDerivative(- mult * value * reciprocal_rval,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult * reciprocal_rval,r);
      right.GetDerivative(- mult * value * reciprocal_rval,r);
    }
  };

  template<class A, class B>
  class addition_expression : public shell_expression<addition_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value;
  public:
    addition_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright) 
    {
      value = left.GetValue() + right.GetValue();
    }
    addition_expression(const addition_expression & other) 
      : left(other.left), right(other.right), value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult,r);
      right.GetDerivative(mult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult,r);
      right.GetDerivative(mult,r);
    }
  };

  template<class A, class B>
  class subtraction_expression : public shell_expression<subtraction_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value;
  public:
    subtraction_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright) 
    {
      value = left.GetValue() - right.GetValue();
    }
    subtraction_expression(const subtraction_expression & other) 
      : left(other.left), right(other.right),value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult,r);
      right.GetDerivative(-mult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult,r);
      right.GetDerivative(-mult,r);
    }
  };

  template<class A, class B>
  class pow_expression : public shell_expression<pow_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value, ldmult, rdmult;
  public:
    pow_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright) 
    {
      INMOST_DATA_REAL_TYPE lval = left.GetValue();
      INMOST_DATA_REAL_TYPE rval = right.GetValue();
      value = ::pow(lval,rval);
      if( lval != 0 )
      {
        ldmult = value * rval / lval;
        rdmult = value * ::log(lval);
      }
      else
      {
        ldmult = 0;
        rdmult = 0;
      }
    }
    pow_expression(const pow_expression & other)
      :left(other.left), right(other.right), value(other.value), 
       ldmult(other.ldmult), rdmult(other.rdmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult*ldmult,r);
      right.GetDerivative(mult*rdmult,r);
    }
  };

  template<class A>
  class pow_const_expression : public shell_expression<pow_const_expression<A> >
  {
    const A & left;
    INMOST_DATA_REAL_TYPE value, ldmult;
  public:
    pow_const_expression(const shell_expression<A> & pleft, INMOST_DATA_REAL_TYPE pright) : left(pleft)
    {
      INMOST_DATA_REAL_TYPE lval = left.GetValue();
      value = ::pow(lval,pright);
      if( lval != 0 )
        ldmult = value * pright / lval;
      else
        ldmult = 0;
    }
    pow_const_expression(const pow_const_expression & other)
      :left(other.left), value(other.value), ldmult(other.ldmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetDerivative(mult*ldmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetDerivative(mult*ldmult,r);
    }
  };

  template<class A>
  class const_pow_expression : public shell_expression<const_pow_expression<A> >
  {
    const A & right;
    INMOST_DATA_REAL_TYPE value, rdmult;
  public:
    const_pow_expression(INMOST_DATA_REAL_TYPE pleft, const shell_expression<A> & pright) : right(pright) 
    {
      INMOST_DATA_REAL_TYPE rval = right.GetValue();
      value = ::pow(pleft,rval);
      if( lval != 0 )
        rdmult = value * ::log(lval);
      else
        rdmult = 0;
    }
    const_pow_expression(const const_pow_expression & other)
      :right(other.right), value(other.value), rdmult(other.rdmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      right.GetDerivative(mult*rdmult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      right.GetDerivative(mult*rdmult,r);
    }
  };

  template<class A, class B, class C>
  class condition_expression : public shell_expression<condition_expression<A,B,C> >
  {
    const A & cond;
    const B & left;
    const C & right;
    INMOST_DATA_REAL_TYPE value, cond_value;
  public:
    condition_expression(const shell_expression<A> & pcond, const shell_expression<B> & pleft, const shell_expression<C> & pright) : cond(pcond), left(pleft), right(pright) 
    {
      cond_value = cond.GetValue();
      value = cond_value >= 0.0 ? left.GetValue() : right.GetValue();
    }
    condition_expression(const condition_expression & other)
      :cond(other.cond), left(other.left), right(other.right), 
      value(other.value), cond_value(other.cond_value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      if( cond_value >= 0.0 )
        left.GetDerivative(mult,r);
      else
        right.GetDerivative(mult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      if( cond_value >= 0.0 )
        left.GetDerivative(mult,r);
      else
        right.GetDerivative(mult,r);
    }
  };

  template<class A> 
  class stencil_expression : public shell_expression<stencil_expression<A> >
  {
    dynarray< std::pair<INMOST_DATA_REAL_TYPE, A >, 64 > arg;
    INMOST_DATA_REAL_TYPE value;
  public:
    stencil_expression(const dynarray< std::pair<INMOST_DATA_REAL_TYPE, A >, 64 > & parg) : arg(parg) 
    {
      value = 0.0;
      for(dynarray< std::pair<INMOST_DATA_REAL_TYPE, A >, 64 >::iterator it = arg.begin(); it != arg.end(); ++it)
        value += it->first * it->second.GetValue();
    }
    stencil_expression(const stencil_expression & other) : arg(other.arg), value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      for(dynarray< std::pair<INMOST_DATA_REAL_TYPE, A >, 64 >::iterator it = arg.begin(); it != arg.end(); ++it)
        it->second.GetDerivative(it->first*mult,r);
    }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      for(dynarray< std::pair<INMOST_DATA_REAL_TYPE, A >, 64 >::iterator it = arg.begin(); it != arg.end(); ++it)
        it->second.GetDerivative(it->first*mult,r);
    }
  };


  

  template<class A, class B, class C> __INLINE   condition_expression<A,B,C> condition(shell_expression<A> const & control, shell_expression<B> const & if_ge_zero, shell_expression<C> const & if_lt_zero) { return condition_expression<A,B,C>(control,if_ge_zero,if_lt_zero); }
                             __INLINE                  INMOST_DATA_REAL_TYPE condition(INMOST_DATA_REAL_TYPE control, INMOST_DATA_REAL_TYPE if_ge_zero, INMOST_DATA_REAL_TYPE if_lt_zero) {return control >= 0.0 ? if_ge_zero : if_lt_zero;}
  template<class A>          __INLINE              unary_minus_expression<A> operator-(shell_expression<A> const & Arg) { return unary_minus_expression<A>(Arg); }
	template<class A>          __INLINE                      abs_expression<A>       abs(shell_expression<A> const & Arg) { return abs_expression<A>(Arg); }
	template<class A>          __INLINE                      exp_expression<A>       exp(shell_expression<A> const & Arg) { return exp_expression<A> (Arg); }
	template<class A>          __INLINE                      log_expression<A>       log(shell_expression<A> const & Arg) { return log_expression<A> (Arg); }
	template<class A>          __INLINE                      sin_expression<A>       sin(shell_expression<A> const & Arg) { return sin_expression<A> (Arg ); }
	template<class A>          __INLINE                      cos_expression<A>       cos(shell_expression<A> const & Arg) { return cos_expression<A> (Arg); }
  template<class A>          __INLINE                     sqrt_expression<A>      sqrt(shell_expression<A> const & Arg) { return sqrt_expression<A> (Arg); }
  template<class A>          __INLINE variation_multiplication_expression<A> variation(shell_expression<A> const & Arg,INMOST_DATA_REAL_TYPE Mult) {return variation_multiplication_expression<A>(Arg,Mult);}
                             __INLINE                  INMOST_DATA_REAL_TYPE variation(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE) {return Arg;}
  template<class A>          __INLINE                  INMOST_DATA_REAL_TYPE get_value(shell_expression<A> const & Arg) { return Arg.GetValue(); }
                             __INLINE                  INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE Arg) {return Arg;}
  template<class A>          __INLINE                 soft_abs_expression<A>  soft_abs(shell_expression<A> const & Arg, INMOST_DATA_REAL_TYPE tol) { return soft_abs_expression<A>(Arg,tol); }
                             __INLINE                  INMOST_DATA_REAL_TYPE  soft_abs(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol) {return ::sqrt(Arg*Arg+tol*tol);}
  template<class A>          __INLINE                soft_sign_expression<A> soft_sign(shell_expression<A> const & Arg, INMOST_DATA_REAL_TYPE tol) { return soft_sign_expression<A>(Arg,tol); }
                             __INLINE                  INMOST_DATA_REAL_TYPE soft_sign(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol) {return Arg/::sqrt(Arg*Arg+tol*tol);}
	
  template<class A, class B> __INLINE        multiplication_expression<A, B> operator*(shell_expression<A> const & Left, shell_expression<B> const & Right) { return multiplication_expression<A, B> (Left, Right); }
	template<class A, class B> __INLINE              division_expression<A, B> operator/(shell_expression<A> const & Left, shell_expression<B> const & Right) { return division_expression<A, B> (Left, Right); }
	template<class A, class B> __INLINE              addition_expression<A, B> operator+(shell_expression<A> const & Left, shell_expression<B> const & Right) { return addition_expression<A, B> (Left, Right); }
	template<class A, class B> __INLINE           subtraction_expression<A, B> operator-(shell_expression<A> const & Left, shell_expression<B> const & Right) { return subtraction_expression<A, B> (Left, Right); }
	template<class A, class B> __INLINE                   pow_expression<A, B>       pow(shell_expression<A> const & Left, shell_expression<B> const & Right) { return pow_expression<A, B> (Left, Right); }
  template<class A, class B> __INLINE              soft_max_expression<A, B>  soft_max(shell_expression<A> const & Left, shell_expression<B> const & Right ,INMOST_DATA_REAL_TYPE tol) { return soft_max_expression<A, B> (Left, Right,tol); }
                             __INLINE                  INMOST_DATA_REAL_TYPE  soft_max(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol) {return 0.5*(Left+Right+::sqrt((Left-Right)*(Left-Right)+tol*tol));}
  template<class A, class B> __INLINE              soft_min_expression<A, B>  soft_min(shell_expression<A> const & Left, shell_expression<B> const & Right ,INMOST_DATA_REAL_TYPE tol) { return soft_min_expression<A, B> (Left, Right,tol); }
                             __INLINE                  INMOST_DATA_REAL_TYPE  soft_min(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol) {return 0.5*(Left+Right-::sqrt((Left-Right)*(Left-Right)+tol*tol));}

  template<class B>          __INLINE                const_pow_expression<B>       pow(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return const_pow_expression<B> (Left, Right); }
  template<class A>          __INLINE                pow_const_expression<A>       pow(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return pow_const_expression<A> (Left, Right); }
  template<class B>          __INLINE     const_multiplication_expression<B> operator*(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return const_multiplication_expression<B>(Right,Left); }
	template<class A>          __INLINE     const_multiplication_expression<A> operator*(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return const_multiplication_expression<A>(Left,Right); }
	template<class B>          __INLINE               reciprocal_expression<B> operator/(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return reciprocal_expression<B>(Right,Left); }
	template<class A>          __INLINE           const_division_expression<A> operator/(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return const_division_expression<A>(Left, Right); }
	template<class B>          __INLINE           const_addition_expression<B> operator+(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return const_addition_expression<B>(Right,Left); }
	template<class A>          __INLINE           const_addition_expression<A> operator+(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return const_addition_expression<A>(Left,Right); }
	template<class B>          __INLINE        const_subtraction_expression<B> operator-(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return const_subtraction_expression<B>(Right, Left); }
	template<class A>          __INLINE           const_addition_expression<A> operator-(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return const_addition_expression<A>(Left, -Right); }


  template<class A, class B> __INLINE bool operator == (shell_expression<A> const & Left, shell_expression<B> const & Right) {return Left.GetValue() == Right.GetValue();}
  template<class A, class B> __INLINE bool operator != (shell_expression<A> const & Left, shell_expression<B> const & Right) {return Left.GetValue() != Right.GetValue();}
  template<class A, class B> __INLINE bool operator <  (shell_expression<A> const & Left, shell_expression<B> const & Right) {return Left.GetValue() <  Right.GetValue();}
  template<class A, class B> __INLINE bool operator >  (shell_expression<A> const & Left, shell_expression<B> const & Right) {return Left.GetValue() >  Right.GetValue();}
  template<class A, class B> __INLINE bool operator <= (shell_expression<A> const & Left, shell_expression<B> const & Right) {return Left.GetValue() <= Right.GetValue();}
  template<class A, class B> __INLINE bool operator >= (shell_expression<A> const & Left, shell_expression<B> const & Right) {return Left.GetValue() >= Right.GetValue();}

  template<class A>          __INLINE bool operator == (shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() == Right;}
  template<class A>          __INLINE bool operator != (shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() != Right;}
  template<class A>          __INLINE bool operator <  (shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() <  Right;}
  template<class A>          __INLINE bool operator >  (shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() >  Right;}
  template<class A>          __INLINE bool operator <= (shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() <= Right;}
  template<class A>          __INLINE bool operator >= (shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() >= Right;}


  template<class B>          __INLINE bool operator == (INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) {return Left == Right.GetValue();}
  template<class B>          __INLINE bool operator != (INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) {return Left != Right.GetValue();}
  template<class B>          __INLINE bool operator <  (INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) {return Left <  Right.GetValue();}
  template<class B>          __INLINE bool operator >  (INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) {return Left >  Right.GetValue();}
  template<class B>          __INLINE bool operator <= (INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) {return Left <= Right.GetValue();}
  template<class B>          __INLINE bool operator >= (INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) {return Left >= Right.GetValue();}

  /*
  template<typename T>
  class get_value
  {
    operator INMOST_DATA_REAL_TYPE() = 0;
  };

  template<>
  class get_value<INMOST_DATA_REAL_TYPE>
  {
    INMOST_DATA_REAL_TYPE val;
  public:
    get_value(const INMOST_DATA_REAL_TYPE & v) :val(v) {}
    operator INMOST_DATA_REAL_TYPE() {return val;}
  };

  template<>
  class get_value<var_expression>
  {
    INMOST_DATA_REAL_TYPE val;
  public:
    get_value(const var_expression & v) :val(v.GetValue()) {}
    operator INMOST_DATA_REAL_TYPE() {return val;}
  };

  template<>
  class get_value<multivar_expression>
  {
    INMOST_DATA_REAL_TYPE val;
  public:
    get_value(const multivar_expression & v) :val(v.GetValue()) {}
    operator INMOST_DATA_REAL_TYPE() {return val;}
  };
  */

  __INLINE bool check_nans(INMOST_DATA_REAL_TYPE val) {return val != val;}
  __INLINE bool check_nans(var_expression const & e) {return e.check_nans();}
  __INLINE bool check_nans(multivar_expression const & e) {return e.check_nans();}


  typedef multivar_expression variable;
  typedef var_expression unknown;  
};

#endif
#endif
