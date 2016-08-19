#ifndef INMOST_AUTODIFF_ETEXPR_H_INCLUDED
#define INMOST_AUTODIFF_ETEXPR_H_INCLUDED
#include "inmost_common.h"
#include "inmost_sparse.h"
#include <sstream> //for debug
#include <new>


//TODO:
// 1. Incorporate tables
// 2. (ok, test) implement condition
// 3. (ok, test) implement stencil
// 4. (???) copying of basic_dynamic_variable
// 5. Consider optimization by checking zero variation multipliers, check that assembly do not degrade.
// 6. floor, ceil, atan, acos, asin, max, min functions
// 7. choice of directional derivatives at discontinuities for abs, pow, max, min (see ADOL-C)
// 8. replace stencil with foreach for provided iterators
// 9. enclose in namespace
//10. CheckCurrentAutomatizator -> CheckCurrentRowMerger
//10.0 structure/service to handle multiple RowMerger objects in openmp environment
//10.1 user should be able to provide RowMerger when Automatizator is not compiled
//10.2 Automatizator may provide internal structure for RowMerger

#ifdef _MSC_VER
#pragma warning(disable : 4503)
#endif

#if defined(USE_AUTODIFF)
namespace INMOST
{
	

	class basic_expression
	{
	public:
		basic_expression() {}//if( GetAutodiffPrint() ) std::cout << this << " Created" << std::endl;}
		basic_expression(const basic_expression & other) {};//std::cout << this << " Created from " << &other << std::endl;}
		virtual INMOST_DATA_REAL_TYPE GetValue() const = 0;
		virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const = 0;
		virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const = 0;
		virtual void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const = 0;
		virtual ~basic_expression() {}//if( GetAutodiffPrint() ) std::cout << this << " Destroied" << std::endl;}
	};

	bool CheckCurrentAutomatizator();
	void FromBasicExpression(Sparse::Row & entries, const basic_expression & expr);
	void AddBasicExpression(Sparse::Row & entries, INMOST_DATA_REAL_TYPE multme, INMOST_DATA_REAL_TYPE multit, const basic_expression & expr);
	void FromGetJacobian(const basic_expression & expr, INMOST_DATA_REAL_TYPE mult, Sparse::Row & r);
	//bool GetAutodiffPrint();
	//void SetAutodiffPrint(bool set);

	template<class Derived>
	class shell_expression : virtual public basic_expression
	{
	public:
		shell_expression() {}// if( GetAutodiffPrint() ) std::cout << this << " Shell Created for " << dynamic_cast<basic_expression *>(this) << std::endl;}
		shell_expression(const shell_expression & other) {}//std::cout << this << " Shell Created from " << &other << std::endl;}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue() const {return static_cast<const Derived *>(this)->GetValue(); }
		__INLINE virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const { return static_cast<const Derived *>(this)->GetJacobian(mult,r); }
		__INLINE virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const { return static_cast<const Derived *>(this)->GetJacobian(mult,r); }
		__INLINE virtual void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {return static_cast<const Derived *>(this)->GetHessian(multJ,J,multH,H); }
		operator Derived & () {return *static_cast<Derived *>(this);}
		operator const Derived & () const {return *static_cast<const Derived *>(this);}
		~shell_expression() {}// if( GetAutodiffPrint() ) std::cout << this << " Shell Destroied for " << dynamic_cast<basic_expression *>(this) << std::endl;}
		//Derived * GetDerived() { return dynamic_cast<Derived *>(this); }
	};

  
 
	class const_expression : public shell_expression<const_expression>
	{
		INMOST_DATA_REAL_TYPE value;
	public:
		const_expression(const const_expression & other) :value(other.value) {}
		const_expression(INMOST_DATA_REAL_TYPE pvalue) : value(pvalue) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {}
		__INLINE const_expression & operator =(const_expression const & other)
		{
			value = other.value;
			return *this;
		}
		operator INMOST_DATA_REAL_TYPE () {return value;}
	};
  

	class var_expression : public shell_expression<var_expression>
	{
		INMOST_DATA_REAL_TYPE value;
		INMOST_DATA_ENUM_TYPE index;
	public:
		var_expression(const var_expression & other) :value(other.value), index(other.index) {}
		var_expression(INMOST_DATA_REAL_TYPE pvalue, INMOST_DATA_ENUM_TYPE pindex) : value(pvalue), index(pindex) {}
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {if( index != ENUMUNDEF ) r[index] += mult;}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {if( index != ENUMUNDEF ) r[index] += mult;}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {if( index != ENUMUNDEF ) J.Push(index,multJ);}
		__INLINE var_expression & operator =(var_expression const & other)
		{
			value = other.value;
			index = other.index;
			return *this;
		}
		__INLINE var_expression & operator =(INMOST_DATA_REAL_TYPE other)
		{
			value = other;
			index = ENUMUNDEF;
			return *this;
		}
		bool check_nans() const
		{
			return value != value;
		}
	};
	

#if defined(PACK_ARRAY)
#pragma pack(push,r1,4)
#endif


	class multivar_expression : public shell_expression<multivar_expression>
	{
		INMOST_DATA_REAL_TYPE value;
		Sparse::Row entries;
		Sparse::HessianRow hessian_entries;
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
			value = expr.GetValue();
			if( CheckCurrentAutomatizator() )
				FromBasicExpression(entries,expr); //Optimized version
			else expr.GetJacobian(1.0,entries);
		}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
		{
			for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
				r[it->first] += it->second*mult;
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
		{
			if( CheckCurrentAutomatizator() )
				FromGetJacobian(*this,mult,r);
			else
			{
				for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
					r[it->first] += it->second*mult;
			}
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			J = entries;
			if( !J.isSorted() ) std::sort(J.Begin(),J.End());
			for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second *= multJ;
			H = hessian_entries;
			for(Sparse::HessianRow::iterator it = H.Begin(); it != H.End(); ++it) it->second *= multH;
		}
		__INLINE multivar_expression & operator = (INMOST_DATA_REAL_TYPE pvalue)
		{
			value = pvalue;
			entries.Clear();
			hessian_entries.Clear();
			return *this;
		}
		__INLINE multivar_expression & operator = (basic_expression const & expr)
		{
			value = expr.GetValue();
			if( CheckCurrentAutomatizator() )
				FromBasicExpression(entries,expr);
			else
			{
				Sparse::Row tmp;
				expr.GetJacobian(1.0,tmp);
				entries.Swap(tmp);
			}
			return *this;
		}
		__INLINE multivar_expression & operator = (multivar_expression const & other)
		{
			value = other.value;
			entries = other.entries;
			hessian_entries = other.hessian_entries;
			return *this;
		}
		__INLINE Sparse::Row & GetRow() {return entries;}
		__INLINE Sparse::HessianRow & GetHessianRow() {return hessian_entries;}
		__INLINE const Sparse::Row & GetRow() const {return entries;}
		__INLINE const Sparse::HessianRow & GetHessianRow() const {return hessian_entries;}
		__INLINE multivar_expression & operator +=(basic_expression const & expr)
		{
			value += expr.GetValue();
			if( CheckCurrentAutomatizator() )
				AddBasicExpression(entries,1.0,1.0,expr);
			else
			{ 
				Sparse::Row tmp(entries);
				expr.GetJacobian(1.0,tmp);
				entries.Swap(tmp);
			}
			return *this;
		}
		__INLINE multivar_expression & operator -=(basic_expression const & expr)
		{
			value -= expr.GetValue();
			if( CheckCurrentAutomatizator() )
				AddBasicExpression(entries,1.0,-1.0,expr);
			else
			{
				Sparse::Row tmp(entries);
				expr.GetJacobian(-1.0,tmp);
				entries.Swap(tmp);
			}
			return *this;
		}
		__INLINE multivar_expression & operator *=(basic_expression const & expr)
		{
			INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
			if( CheckCurrentAutomatizator() )
				AddBasicExpression(entries,rval,lval,expr);
			else
			{
				Sparse::Row tmp(entries);
				for(Sparse::Row::iterator it = tmp.Begin(); it != tmp.End(); ++it) it->second *= rval;
				expr.GetJacobian(lval,tmp);
				entries.Swap(tmp);
			}
			value *= rval;
			return *this;
		}
		__INLINE multivar_expression & operator /=(basic_expression const & expr)
		{
			INMOST_DATA_REAL_TYPE rval = expr.GetValue();
			INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0/rval;
			value *= reciprocial_rval;
			if( CheckCurrentAutomatizator() )
				AddBasicExpression(entries,reciprocial_rval,-value*reciprocial_rval,expr);
			else
			{
				Sparse::Row tmp(entries);
				for(Sparse::Row::iterator it = tmp.Begin(); it != tmp.End(); ++it) it->second *= reciprocial_rval;
				expr.GetJacobian(-value*reciprocial_rval,tmp); 
				entries.Swap(tmp);
			}
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
		/// Write variable into array of entries.
		/// Size of array can be determined via RecordSize.
		/// Used internally by Mesh::GetData.
		/// @param v Array of entries that will store data of the variable.
		/// @return Number of entries used.
		INMOST_DATA_ENUM_TYPE Record(Sparse::Row::entry * v) const
		{
			INMOST_DATA_ENUM_TYPE k = 0;
			v[k].first = entries.Size();
			v[k].second = value;
			k++;
			for(INMOST_DATA_ENUM_TYPE r = 0; r < entries.Size(); ++r)
			{
				v[k].first = entries.GetIndex(r);
				v[k].second = entries.GetValue(r);
				k++;
			}
			return k;
		}
		/// Number of entries required to record the variable.
		INMOST_DATA_ENUM_TYPE RecordSize() const
		{
			return 1 + entries.Size();
		}
		/// Retrive variable from array of entries.
		/// Size of array without retrival can be determined via RetriveSize.
		/// @param v Array of entries that will store data of the variable.
		/// @return Number of entries red.
		INMOST_DATA_ENUM_TYPE Retrive(const Sparse::Row::entry * v)
		{
			int k = 0;
			value = v[k].second;
			entries.Resize(v[k].first);
			k++;
			for(int r = 0; r < (int)entries.Size(); ++r)
			{
				entries.GetIndex(r) = v[k].first;
				entries.GetValue(r) = v[k].second;
				k++;
			}
			return k;
		}
		/// Number of entries used.
		static INMOST_DATA_ENUM_TYPE RetriveSize(const Sparse::Row::entry * v)
		{
			return 1 + v[0].first;
		}
		
		friend class multivar_expression_reference;
	};

#if defined(PACK_ARRAY)
#pragma pack(pop,r1)
#endif


  
  class multivar_expression_reference : public shell_expression<multivar_expression_reference>
  {
    INMOST_DATA_REAL_TYPE & value;
    Sparse::Row * entries;
    Sparse::HessianRow * hessian_entries;
  public:
    /// Constructor, set links to the provided value and entries
    multivar_expression_reference(INMOST_DATA_REAL_TYPE & _value, Sparse::Row * _entries, Sparse::HessianRow * _hentries = NULL) 
      : value(_value), entries(_entries), hessian_entries(_hentries) {}
    /// Copy constructor, sets links to the same reference of value and entries
    multivar_expression_reference(const multivar_expression_reference & other) 
      : value(other.value), entries(other.entries) {}
	  /// Copy constructor from multivar_expression, sets links to the same reference of value and entries
	  multivar_expression_reference(multivar_expression & other)
	  : value(other.value), entries(&other.entries) {}
    /// Retrive value
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    /// Set value without changing derivatives
    __INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
    /// Retrive derivatives with multiplier into Sparse::RowMerger structure.
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it)
        r[it->first] += it->second*mult;
    }
    /// Retrive derivatives with multiplier into Sparse::Row structure.
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      if( CheckCurrentAutomatizator() )
        FromGetJacobian(*this,mult,r);
      else
      {
        for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it)
          r[it->first] += it->second*mult;
      }
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J,INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        J = *entries;
        if( !J.isSorted() ) std::sort(J.Begin(),J.End());
        for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second *= multJ;
        if( hessian_entries )
        {
            H = *hessian_entries;
            for(Sparse::HessianRow::iterator it = H.Begin(); it != H.End(); ++it) it->second *= multH;
        }
    }
    __INLINE multivar_expression_reference & operator = (INMOST_DATA_REAL_TYPE pvalue)
    {
      value = pvalue;
      entries->Clear();
      return *this;
    }
    __INLINE multivar_expression_reference & operator = (basic_expression const & expr)
    {
      value = expr.GetValue();
      if( CheckCurrentAutomatizator() )
        FromBasicExpression(*entries,expr);
      else
      {
        Sparse::Row tmp;
        expr.GetJacobian(1.0,tmp);
        entries->Swap(tmp);
      }
      return *this;
    }
    __INLINE multivar_expression_reference & operator = (multivar_expression_reference const & other)
    {
      value = other.GetValue();
      *entries = other.GetRow();
      return *this;
    }
    __INLINE multivar_expression_reference & operator = (multivar_expression const & other)
    {
      value = other.GetValue();
      *entries = other.GetRow();
      return *this;
    }
    __INLINE Sparse::Row & GetRow() {return *entries;}
    __INLINE const Sparse::Row & GetRow() const {return *entries;}
    __INLINE Sparse::HessianRow & GetHessianRow() {return *hessian_entries;}
    __INLINE const Sparse::HessianRow & GetHessianRow() const {return *hessian_entries;}
    __INLINE multivar_expression_reference & operator +=(basic_expression const & expr)
    {
      value += expr.GetValue();
      if( CheckCurrentAutomatizator() )
        AddBasicExpression(*entries,1.0,1.0,expr);
      else
      { 
        Sparse::Row tmp(*entries);
        expr.GetJacobian(1.0,tmp);
        entries->Swap(tmp);
      }
      return *this;
    }
    __INLINE multivar_expression_reference & operator -=(basic_expression const & expr)
    {
      value -= expr.GetValue();
      if( CheckCurrentAutomatizator() )
        AddBasicExpression(*entries,1.0,-1.0,expr);
      else
      {
        Sparse::Row tmp(*entries);
        expr.GetJacobian(-1.0,tmp);
        entries->Swap(tmp);
      }
      return *this;
    }
    __INLINE multivar_expression_reference & operator *=(basic_expression const & expr)
    {
      INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
      if( CheckCurrentAutomatizator() )
        AddBasicExpression(*entries,rval,lval,expr);
      else
      {
        Sparse::Row tmp(*entries);
        for(Sparse::Row::iterator it = tmp.Begin(); it != tmp.End(); ++it) it->second *= rval;
        expr.GetJacobian(lval,tmp);
        entries->Swap(tmp);
      }
      value *= rval;
      return *this;
    }
    __INLINE multivar_expression_reference & operator /=(basic_expression const & expr)
    {
      INMOST_DATA_REAL_TYPE rval = expr.GetValue();
      INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0/rval;
      value *= reciprocial_rval;
      if( CheckCurrentAutomatizator() )
        AddBasicExpression(*entries,reciprocial_rval,-value*reciprocial_rval,expr);
      else
      {
        Sparse::Row tmp(*entries);
        for(Sparse::Row::iterator it = tmp.Begin(); it != tmp.End(); ++it) it->second *= reciprocial_rval;
        expr.GetJacobian(-value*reciprocial_rval,tmp); 
        entries->Swap(tmp);
      }
      return *this;
    }
    __INLINE multivar_expression_reference & operator +=(INMOST_DATA_REAL_TYPE right)
    {
      value += right;
      return *this;
    }
    __INLINE multivar_expression_reference & operator -=(INMOST_DATA_REAL_TYPE right)
    {
      value -= right;
      return *this;
    }
    __INLINE multivar_expression_reference & operator *=(INMOST_DATA_REAL_TYPE right)
    {
      value *= right;
      for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it) it->second *= right;
      return *this;
    }
    __INLINE multivar_expression_reference & operator /=(INMOST_DATA_REAL_TYPE right)
    {
      value /= right;
      for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it) it->second /= right;
      return *this;
    }
    bool check_nans() const
    {
      if( value != value ) return true;
      for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it)
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        arg.GetHessian(multJ*dmult,J,multH*dmult,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
      arg.GetHessian(multJ*dmult,J,multH*dmult,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        arg.GetHessian(multJ*dmult,J,multH*dmult,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetJacobian(mult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetJacobian(mult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        arg.GetHessian(multJ,J,multH,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetJacobian(-mult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetJacobian(-mult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        arg.GetHessian(-multJ,J,-multH,H);
    }
  };

  /// (c/x)' = -c dx / (x*x)
  /// (c/x)'' = 2 c dx dx / (x*x*x)
  template<class A>
  class reciprocal_expression : public shell_expression<reciprocal_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value, reciprocial_val;
  public:
    reciprocal_expression(const shell_expression<A> & parg,INMOST_DATA_REAL_TYPE coef) : arg(parg)
    {
      reciprocial_val = 1.0/arg.GetValue();
      value = coef*reciprocial_val;
    }
    reciprocal_expression(const reciprocal_expression & other) 
      : arg(other.arg), value(other.value), reciprocial_val(other.reciprocial_val) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetJacobian(-mult*value*reciprocial_val,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetJacobian(-mult*value*reciprocial_val,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        Sparse::HessianRow ArgH;
        double coefJ = -multJ*value*reciprocial_val;
        double signJ = coefJ < 0 ? -1 : 1;
        arg.GetHessian(signJ,J,-2*multH*value*reciprocial_val,ArgH);
        Sparse::HessianRow::MergeJacobianHessian(2*value*reciprocial_val*reciprocial_val*signJ,J,J,1.0,ArgH,H);
        for(INMOST_DATA_ENUM_TYPE k = 0; k < J.Size(); ++k) J.GetValue(k) *= coefJ*signJ;
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
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
      {
          arg.GetJacobian(-mult,r);
      }
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
      {
          arg.GetJacobian(-mult,r);
      }
      __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
      {
          arg.GetHessian(-multJ,J,-multH,H);
      }
  };

  template<class A>
  class unary_plus_expression : public shell_expression<unary_plus_expression<A> >
  {
      const A & arg;
      INMOST_DATA_REAL_TYPE value;
  public:
      unary_plus_expression(const shell_expression<A> & parg) : arg(parg) {value = arg.GetValue();}
      unary_plus_expression(const unary_plus_expression & b) : arg(b.arg) {}
      __INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
      {
          arg.GetJacobian(mult,r);
      }
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
      {
          arg.GetJacobian(mult,r);
      }
      __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
      {
          arg.GetHessian(multJ,J,multH,H);
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
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
      {
          arg.GetJacobian( (value == 0 ? (mult < 0.0 ? -1 : 1) : 1) * mult * dmult, r);
      }
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
      {
          arg.GetJacobian( (value == 0 ? (mult < 0.0 ? -1 : 1) : 1) * mult * dmult, r);
      }
      __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ,Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
      {
          double a = (value == 0 ? (multJ < 0.0 ? -1 : 1) : 1);
          double b = (value == 0 ? (multH < 0.0 ? -1 : 1) : 1);
          arg.GetHessian( a * multJ * dmult,J, b*multH * dmult,H);
      }
  };

  // ex
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
      exp_expression(const exp_expression & b) : arg(b.arg), value(b.value) {}
      __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
      {
          arg.GetJacobian( mult * value, r);
      }
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
      {
          arg.GetJacobian( mult * value, r);
      }
      __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
      {
          Sparse::HessianRow ArgH;
          double coefJ = multJ*value;
          double signJ = coefJ < 0.0 ? -1 : 1;
          arg.GetHessian(signJ, J, multH*value, ArgH); //check
          Sparse::HessianRow::MergeJacobianHessian(coefJ*signJ,J,J,1.0,ArgH,H);
          for(INMOST_DATA_ENUM_TYPE k = 0; k < J.Size(); ++k) J.GetValue(k) *= coefJ*signJ;
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
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
      {
          arg.GetJacobian(mult*dmult,r);
      }
      __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
      {
          arg.GetJacobian(mult*dmult,r);
      }
      __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
      {
          Sparse::HessianRow ArgH;
          double coefJ = multJ*dmult;
          double signJ = coefJ < 0.0 ? -1 : 1;
          arg.GetHessian(signJ, J, 2*multH*dmult, ArgH); //check
          Sparse::HessianRow::MergeJacobianHessian(-coefJ*signJ*dmult,J,J,1.0,ArgH,H);
          for(INMOST_DATA_ENUM_TYPE k = 0; k < J.Size(); ++k) J.GetValue(k) *= coefJ*signJ;
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        arg.GetHessian(multJ*dmult,J,-multH*value,H);
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        arg.GetHessian(multJ*dmult,J,-multH*value,H);
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetJacobian(0.5*mult/value,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetJacobian(0.5*mult/value,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        arg.GetHessian(0.5*multJ/value,J,-0.25*multH/::pow(value,3),H);
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const 
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult*ldmult,r);
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult*ldmult,r);
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult*ldmult,r);
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult*ldmult,r);
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
    }
  };


  template<class A, class B>
  class multiplication_expression : public shell_expression<multiplication_expression<A,B> >
  {
    const A & left;
    const B & right;
    INMOST_DATA_REAL_TYPE value;
  public:
    multiplication_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright)
    {
      value = left.GetValue()*right.GetValue();
    }
    multiplication_expression(const multiplication_expression & other) 
      : left(other.left), right(other.right), value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult*right.GetValue(),r);
      right.GetJacobian(mult*left.GetValue(),r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult*right.GetValue(),r);
      right.GetJacobian(mult*left.GetValue(),r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
        Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
        left.GetHessian(multJ,JL,multH,HL); //retrive jacobian row and hessian matrix of the left expression
        right.GetHessian(multJ,JR,multH,HR); //retrive jacobian row and hessian matrix of the right expression
        //assume rows are sorted (this is to be ensured by corresponding GetHessian functions)
        //preallocate J to JL.Size+JR.Size
        //perform merging of two sorted arrays
        //resize to correct size
        Sparse::Row::MergeSortedRows(right.GetValue(),JL,left.GetValue(),JR,J);
        //preallocate H to HL.Size+HR.Size+JL.Size*JR.Size
        //merge sorted
        Sparse::HessianRow::MergeJacobianHessian(2.0,JL,JR,right.GetValue(),HL,left.GetValue(),HR,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult * reciprocal_rval,r);
      right.GetJacobian(- mult * value * reciprocal_rval,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult * reciprocal_rval,r);
      right.GetJacobian(- mult * value * reciprocal_rval,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
        Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
        Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
        left.GetHessian(multJ,JL,multH,HL); //retrive jacobian row and hessian matrix of the left expression
        right.GetHessian(multJ,JR,multH,HR); //retrive jacobian row and hessian matrix of the right expression
        //assume rows are sorted (this is to be ensured by corresponding GetHessian functions)
        //preallocate J to JL.Size+JR.Size
        //perform merging of two sorted arrays
        //resize to correct size
        Sparse::Row::MergeSortedRows(reciprocal_rval,JL,-value * reciprocal_rval,JR,J);
        //preallocate H to HL.Size+HR.Size+JL.Size*JR.Size
        //merge sorted
        Sparse::HessianRow::MergeJacobianHessian(2.0,JL,JR,reciprocal_rval,HL,2*left.GetValue()*multH,HR,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult,r);
      right.GetJacobian(mult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult,r);
      right.GetJacobian(mult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
        Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
        left.GetHessian(multJ,JL,multH,HL); //retrive jacobian row and hessian matrix of the left expression
        right.GetHessian(multJ,JR,multH,HR); //retrive jacobian row and hessian matrix of the right expression
        Sparse::Row::MergeSortedRows(1.0,JL,1.0,JR,J);
        Sparse::HessianRow::MergeSortedRows(1.0,HL,1.0,HR,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult,r);
      right.GetJacobian(-mult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult,r);
      right.GetJacobian(-mult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
        Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
        left.GetHessian(multJ,JL,multH,HL); //retrive jacobian row and hessian matrix of the left expression
        right.GetHessian(multJ,JR,multH,HR); //retrive jacobian row and hessian matrix of the right expression
        Sparse::Row::MergeSortedRows(1.0,JL,-1.0,JR,J);
        Sparse::HessianRow::MergeSortedRows(1.0,HL,-1.0,HR,H);
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult*ldmult,r);
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult*ldmult,r);
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      left.GetJacobian(mult*ldmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      left.GetJacobian(mult*ldmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
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
      if( pleft != 0 )
        rdmult = value * ::log(pleft);
      else
        rdmult = 0;
    }
    const_pow_expression(const const_pow_expression & other)
      :right(other.right), value(other.value), rdmult(other.rdmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      right.GetJacobian(mult*rdmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        throw NotImplemented;
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
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      if( cond_value >= 0.0 )
        left.GetJacobian(mult,r);
      else
        right.GetJacobian(mult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      if( cond_value >= 0.0 )
        left.GetJacobian(mult,r);
      else
        right.GetJacobian(mult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
      if( cond_value >= 0.0 )
        left.GetHessian(multJ,J,multH,H);
      else
        right.GetHessian(multJ,J,multH,H);
    }
  };

  template<class A, class B>
  class branch_expression : public shell_expression<branch_expression<A,B> >
  {
	bool cond;
    const A & left;
    const B & right;
	INMOST_DATA_REAL_TYPE value;
  public:
    branch_expression(bool pcond, const shell_expression<A> & pleft, const shell_expression<B> & pright) : cond(pcond), left(pleft), right(pright) 
    {
      value = cond ? left.GetValue() : right.GetValue();
    }
    branch_expression(const condition_expression & other)
      :cond(other.cond), left(other.left), right(other.right), 
      value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      if( cond )
        left.GetJacobian(mult,r);
      else
        right.GetJacobian(mult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      if( cond )
        left.GetJacobian(mult,r);
      else
        right.GetJacobian(mult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
      if( cond )
        left.GetHessian(multJ,J,multH,H);
      else
        right.GetHessian(multJ,J,multH,H);
    }
  };

  template<class A> 
  class stencil_expression : public shell_expression<stencil_expression<A> >
  {
    dynarray< const_multiplication_expression<A>, 64 > arg;
    INMOST_DATA_REAL_TYPE value;
  public:
    stencil_expression(const dynarray< const_multiplication_expression<A>, 64 > & parg) : arg(parg) 
    {
      value = 0.0;
      for(typename dynarray< const_multiplication_expression<A>, 64 >::iterator it = arg.begin(); it != arg.end(); ++it)
        value += it->GetValue();
    }
    stencil_expression(const stencil_expression & other) : arg(other.arg), value(other.value) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      for(typename dynarray< const_multiplication_expression<A>, 64 >::iterator it = arg.begin(); it != arg.end(); ++it)
        it->GetJacobian(mult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      for(typename dynarray< const_multiplication_expression<A>, 64 >::iterator it = arg.begin(); it != arg.end(); ++it)
        it->GetJacobian(mult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
        Sparse::Row tmpJ, curJ;
        Sparse::HessianRow tmpH, curH;
        for(typename dynarray< const_multiplication_expression<A>, 64 >::iterator it = arg.begin(); it != arg.end(); ++it)
        {
            curJ.Clear();
            curH.Clear();
            it->GetHessian(multJ,curJ,multH,curH);
            Sparse::Row::MergeSortedRows(1.0,curJ,1.0,J,tmpJ);
            Sparse::HessianRow::MergeSortedRows(1.0,curH,1.0,H,tmpH);
            J.Swap(tmpJ);
            H.Swap(tmpH);
        }
    }
	void SetCondition(bool _cond) { cond = _cond; }
  };


  template<class A>
  class function_expression : public shell_expression< function_expression<A> >
  {
    const A & arg;
    INMOST_DATA_REAL_TYPE value, dmult, ddmult;
  public:
    function_expression(const shell_expression<A> & _arg)
      :arg(_arg), value(1), dmult(0) {}
    function_expression(const shell_expression<A> & _arg, INMOST_DATA_REAL_TYPE pvalue, INMOST_DATA_REAL_TYPE pdmult, INMOST_DATA_REAL_TYPE pddmult = 0)
      :arg(_arg), value(pvalue), dmult(pdmult), ddmult(pddmult) {}
    function_expression(const function_expression & other) 
      : arg(other.arg), value(other.value), dmult(other.dmult), ddmult(other.ddmult) {}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
    {
      arg.GetJacobian(mult*dmult,r);
    }
    __INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
    {
      arg.GetHessian(multJ*dmult,J,multH*ddmult,H);
    }
    void SetFunctionValue(INMOST_DATA_REAL_TYPE val) {value = val;}
    void SetFunctionDerivative(INMOST_DATA_REAL_TYPE val) {dmult = val;}
  };

  
  class keyval_table
  {
    std::string name;
    INMOST_DATA_REAL_TYPE * vals;
    INMOST_DATA_REAL_TYPE * args;
    INMOST_DATA_ENUM_TYPE size;
    INMOST_DATA_ENUM_TYPE binary_search(INMOST_DATA_REAL_TYPE arg) const
	  {
		  int l = 0, r = static_cast<int>(size)-1, mid = 0;
		  while (r >= l)
		  {
			  mid = (l + r) / 2;
			  if (args[mid] > arg) r = mid - 1;
			  else if (args[mid] < arg) l = mid + 1;
			  else return mid;
		  }
		  mid = (l + r) / 2;
		  if (mid > static_cast<int>(size - 2)) mid = static_cast<int>(size - 2);
		  return static_cast<INMOST_DATA_ENUM_TYPE>(mid);
	  }
  public:
    INMOST_DATA_REAL_TYPE GetValue(INMOST_DATA_REAL_TYPE arg) const
	  {
		  if (arg < args[0]) return vals[0];
		  INMOST_DATA_ENUM_TYPE i = binary_search(arg);
		  return vals[i] + (vals[i + 1] - vals[i]) * (arg - args[i]) / (args[i + 1] - args[i]);
	  }
    INMOST_DATA_REAL_TYPE GetDerivative(INMOST_DATA_REAL_TYPE arg) const
	  {
		  if (arg < args[0]) return 0.0;
		  INMOST_DATA_ENUM_TYPE i = binary_search(arg);
		  return (vals[i + 1] - vals[i]) / (args[i + 1] - args[i]);
	  }
    std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> GetBoth(INMOST_DATA_REAL_TYPE arg) const
	  {
		  if (arg < args[0]) return std::make_pair(vals[0], 0.0);
		  INMOST_DATA_ENUM_TYPE i = binary_search(arg);
		  INMOST_DATA_REAL_TYPE der = (vals[i + 1] - vals[i]) / (args[i + 1] - args[i]);
		  return std::make_pair(vals[i] + der * (arg - args[i]), der);
	  }
    keyval_table() :name(""), vals(NULL), args(NULL), size(0) {}
    keyval_table(std::string _name, INMOST_DATA_REAL_TYPE * _args, INMOST_DATA_REAL_TYPE * _vals, INMOST_DATA_ENUM_TYPE _size)
    {
      name = _name;
      size = _size;
      vals = new INMOST_DATA_REAL_TYPE[size];
      memcpy(vals,_vals,sizeof(INMOST_DATA_REAL_TYPE)*size);
      //std::copy(_vals,_vals+size,vals);
      args = new INMOST_DATA_REAL_TYPE[size];
      memcpy(args,_args,sizeof(INMOST_DATA_REAL_TYPE)*size);
      //std::copy(_args,_args+size,args);
    }
    keyval_table(const keyval_table & other)
    {
      name = other.name;
      size = other.size;
      vals = new INMOST_DATA_REAL_TYPE[size];
      memcpy(vals,other.vals,sizeof(INMOST_DATA_REAL_TYPE)*size);
      //std::copy(other.vals,other.vals+size,vals);
      args = new INMOST_DATA_REAL_TYPE[size];
      memcpy(args,other.args,sizeof(INMOST_DATA_REAL_TYPE)*size);
      //std::copy(other.args,other.args+size,args);
    }
    keyval_table & operator = (keyval_table const & other)
    {
      Clear();
      name = other.name;
      size = other.size;
      vals = new INMOST_DATA_REAL_TYPE[size];
      memcpy(vals,other.vals,sizeof(INMOST_DATA_REAL_TYPE)*size);
      //std::copy(other.vals,other.vals+size,vals);
      args = new INMOST_DATA_REAL_TYPE[size];
      memcpy(args,other.args,sizeof(INMOST_DATA_REAL_TYPE)*size);
      //std::copy(other.args,other.args+size,args);
      return * this;
    }
    ~keyval_table()
    {
      Clear();
    }
    void Clear()
    {
      name = "";
      if( args ) {delete [] args; args = NULL;}
      if( vals ) {delete [] vals; vals = NULL;}
      size = 0;
    }
    bool Empty() const
    {
      return size == 0;
    }
    const INMOST_DATA_REAL_TYPE * GetTableArguments() const {return args;}
    const INMOST_DATA_REAL_TYPE * GetTableValues() const {return vals;}
    INMOST_DATA_REAL_TYPE * GetTableArguments() {return args;}
    INMOST_DATA_REAL_TYPE * GetTableValues() {return vals;}
    INMOST_DATA_ENUM_TYPE GetSize() const {return size;}
  };

  typedef multivar_expression variable;
  typedef var_expression unknown;
}


__INLINE bool check_nans(INMOST_DATA_REAL_TYPE val) {return val != val;}
__INLINE bool check_nans(INMOST::var_expression const & e) {return e.check_nans();}
__INLINE bool check_nans(INMOST::multivar_expression const & e) {return e.check_nans();}
__INLINE bool check_nans(INMOST::multivar_expression_reference const & e) {return e.check_nans();}

template<class A, class B, class C> __INLINE   INMOST::condition_expression<A,B,C> condition(INMOST::shell_expression<A> const & control, INMOST::shell_expression<B> const & if_ge_zero, INMOST::shell_expression<C> const & if_lt_zero) { return INMOST::condition_expression<A,B,C>(control,if_ge_zero,if_lt_zero); }
                                    __INLINE                 INMOST_DATA_REAL_TYPE condition(INMOST_DATA_REAL_TYPE control, INMOST_DATA_REAL_TYPE if_ge_zero, INMOST_DATA_REAL_TYPE if_lt_zero) {return control >= 0.0 ? if_ge_zero : if_lt_zero;}
template<class A>          __INLINE              INMOST::unary_minus_expression<A> operator-(INMOST::shell_expression<A> const & Arg) { return INMOST::unary_minus_expression<A>(Arg); }
template<class A>          __INLINE               INMOST::unary_plus_expression<A> operator+(INMOST::shell_expression<A> const & Arg) { return INMOST::unary_plus_expression<A>(Arg); }
template<class A>          __INLINE                      INMOST::abs_expression<A>      fabs(INMOST::shell_expression<A> const & Arg) { return INMOST::abs_expression<A>(Arg); }
template<class A>          __INLINE                      INMOST::exp_expression<A>       exp(INMOST::shell_expression<A> const & Arg) { return INMOST::exp_expression<A> (Arg); }
template<class A>          __INLINE                      INMOST::log_expression<A>       log(INMOST::shell_expression<A> const & Arg) { return INMOST::log_expression<A> (Arg); }
template<class A>          __INLINE                      INMOST::sin_expression<A>       sin(INMOST::shell_expression<A> const & Arg) { return INMOST::sin_expression<A> (Arg ); }
template<class A>          __INLINE                      INMOST::cos_expression<A>       cos(INMOST::shell_expression<A> const & Arg) { return INMOST::cos_expression<A> (Arg); }
template<class A>          __INLINE                     INMOST::sqrt_expression<A>      sqrt(INMOST::shell_expression<A> const & Arg) { return INMOST::sqrt_expression<A> (Arg); }
template<class A>          __INLINE INMOST::variation_multiplication_expression<A> variation(INMOST::shell_expression<A> const & Arg,INMOST_DATA_REAL_TYPE Mult) {return INMOST::variation_multiplication_expression<A>(Arg,Mult);}
                           __INLINE                          INMOST_DATA_REAL_TYPE variation(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE) {return Arg;}
template<class A>          __INLINE                          INMOST_DATA_REAL_TYPE get_value(INMOST::shell_expression<A> const & Arg) { return Arg.GetValue(); }
                           __INLINE                          INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE Arg) {return Arg;}
                           __INLINE                                           void set_value(INMOST::var_expression & Arg, INMOST_DATA_REAL_TYPE Val) {Arg.SetValue(Val); }
                           __INLINE                                           void set_value(INMOST::multivar_expression & Arg, INMOST_DATA_REAL_TYPE Val) {Arg.SetValue(Val); }
                           __INLINE                                           void set_value(INMOST::multivar_expression_reference & Arg, INMOST_DATA_REAL_TYPE Val) {Arg.SetValue(Val); }
                           __INLINE                                           void set_value(INMOST_DATA_REAL_TYPE & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val;}
						   __INLINE                                           void set_value(INMOST_DATA_REAL_TYPE & Arg, const INMOST::var_expression & Val) {Arg = Val.GetValue();}
						   __INLINE                                           void set_value(INMOST_DATA_REAL_TYPE & Arg, const INMOST::multivar_expression & Val) {Arg = Val.GetValue();}
						   __INLINE                                           void set_value(INMOST_DATA_REAL_TYPE & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val.GetValue();}
                           __INLINE                                           void set_value(INMOST::multivar_expression & Arg, const INMOST::var_expression & Val) {Arg.SetValue(Val.GetValue()); }
						   __INLINE                                           void set_value(INMOST::multivar_expression & Arg, const INMOST::multivar_expression & Val) {Arg.SetValue(Val.GetValue()); }
						   __INLINE                                           void set_value(INMOST::multivar_expression & Arg, const INMOST::multivar_expression_reference & Val) {Arg.SetValue(Val.GetValue()); }
						   __INLINE                                           void set_value(INMOST::multivar_expression_reference & Arg, const INMOST::var_expression & Val) {Arg.SetValue(Val.GetValue()); }
						   __INLINE                                           void set_value(INMOST::multivar_expression_reference & Arg, const INMOST::multivar_expression & Val) {Arg.SetValue(Val.GetValue()); }
						   __INLINE                                           void set_value(INMOST::multivar_expression_reference & Arg, const INMOST::multivar_expression_reference & Val) {Arg.SetValue(Val.GetValue()); }
                           __INLINE                                           void    assign(INMOST::var_expression & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val;}
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val;}
						   __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::var_expression & Val) {Arg = Val.GetValue();}
						   __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::multivar_expression & Val) {Arg = Val.GetValue();}
						   __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val.GetValue();}
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::var_expression & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::multivar_expression & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::var_expression & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::multivar_expression & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val; }
template<class A>          __INLINE                 INMOST::soft_abs_expression<A> soft_fabs(INMOST::shell_expression<A> const & Arg, INMOST_DATA_REAL_TYPE tol) { return INMOST::soft_abs_expression<A>(Arg,tol); }
                           __INLINE                          INMOST_DATA_REAL_TYPE soft_fabs(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol) {return ::sqrt(Arg*Arg+tol*tol);}
template<class A>          __INLINE                INMOST::soft_sign_expression<A> soft_sign(INMOST::shell_expression<A> const & Arg, INMOST_DATA_REAL_TYPE tol) { return INMOST::soft_sign_expression<A>(Arg,tol); }
                           __INLINE                          INMOST_DATA_REAL_TYPE soft_sign(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol) {return Arg/::sqrt(Arg*Arg+tol*tol);}
template<class A, class B> __INLINE        INMOST::multiplication_expression<A, B> operator*(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::multiplication_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE              INMOST::division_expression<A, B> operator/(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::division_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE              INMOST::addition_expression<A, B> operator+(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::addition_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE           INMOST::subtraction_expression<A, B> operator-(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::subtraction_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE                   INMOST::pow_expression<A, B>       pow(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::pow_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE              INMOST::soft_max_expression<A, B>  soft_max(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right ,INMOST_DATA_REAL_TYPE tol) { return INMOST::soft_max_expression<A, B> (Left, Right,tol); }
                           __INLINE                          INMOST_DATA_REAL_TYPE  soft_max(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol) {return 0.5*(Left+Right+::sqrt((Left-Right)*(Left-Right)+tol*tol));}
template<class A, class B> __INLINE              INMOST::soft_min_expression<A, B>  soft_min(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right ,INMOST_DATA_REAL_TYPE tol) { return INMOST::soft_min_expression<A, B> (Left, Right,tol); }
                           __INLINE                          INMOST_DATA_REAL_TYPE  soft_min(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol) {return 0.5*(Left+Right-::sqrt((Left-Right)*(Left-Right)+tol*tol));}
template<class B>          __INLINE                INMOST::const_pow_expression<B>       pow(INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) { return INMOST::const_pow_expression<B> (Left, Right); }
template<class A>          __INLINE                INMOST::pow_const_expression<A>       pow(INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::pow_const_expression<A> (Left, Right); }
template<class B>          __INLINE     INMOST::const_multiplication_expression<B> operator*(INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) { return INMOST::const_multiplication_expression<B>(Right,Left); }
template<class A>          __INLINE     INMOST::const_multiplication_expression<A> operator*(INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::const_multiplication_expression<A>(Left,Right); }
template<class B>          __INLINE               INMOST::reciprocal_expression<B> operator/(INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) { return INMOST::reciprocal_expression<B>(Right,Left); }
template<class A>          __INLINE           INMOST::const_division_expression<A> operator/(INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::const_division_expression<A>(Left, Right); }
template<class B>          __INLINE           INMOST::const_addition_expression<B> operator+(INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) { return INMOST::const_addition_expression<B>(Right,Left); }
template<class A>          __INLINE           INMOST::const_addition_expression<A> operator+(INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::const_addition_expression<A>(Left,Right); }
template<class B>          __INLINE        INMOST::const_subtraction_expression<B> operator-(INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) { return INMOST::const_subtraction_expression<B>(Right, Left); }
template<class A>          __INLINE           INMOST::const_addition_expression<A> operator-(INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::const_addition_expression<A>(Left, -Right); }
template<class A, class B> __INLINE                                        bool operator == (INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) {return Left.GetValue() == Right.GetValue();}
template<class A, class B> __INLINE                                        bool operator != (INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) {return Left.GetValue() != Right.GetValue();}
template<class A, class B> __INLINE                                        bool operator <  (INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) {return Left.GetValue() <  Right.GetValue();}
template<class A, class B> __INLINE                                        bool operator >  (INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) {return Left.GetValue() >  Right.GetValue();}
template<class A, class B> __INLINE                                        bool operator <= (INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) {return Left.GetValue() <= Right.GetValue();}
template<class A, class B> __INLINE                                        bool operator >= (INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) {return Left.GetValue() >= Right.GetValue();}
template<class A>          __INLINE                                        bool operator == (INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() == Right;}
template<class A>          __INLINE                                        bool operator != (INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() != Right;}
template<class A>          __INLINE                                        bool operator <  (INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() <  Right;}
template<class A>          __INLINE                                        bool operator >  (INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() >  Right;}
template<class A>          __INLINE                                        bool operator <= (INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() <= Right;}
template<class A>          __INLINE                                        bool operator >= (INMOST::shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) {return Left.GetValue() >= Right;}
template<class B>          __INLINE                                        bool operator == (INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) {return Left == Right.GetValue();}
template<class B>          __INLINE                                        bool operator != (INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) {return Left != Right.GetValue();}
template<class B>          __INLINE                                        bool operator <  (INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) {return Left <  Right.GetValue();}
template<class B>          __INLINE                                        bool operator >  (INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) {return Left >  Right.GetValue();}
template<class B>          __INLINE                                        bool operator <= (INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) {return Left <= Right.GetValue();}
template<class B>          __INLINE                                        bool operator >= (INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const & Right) {return Left >= Right.GetValue();}
template<class A>          __INLINE                 INMOST::function_expression<A> get_table(INMOST::shell_expression<A> const & Arg, const INMOST::keyval_table & Table)
{
  std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> both = Table.GetBoth(Arg.GetValue());
  return INMOST::function_expression<A>(Arg,both.first,both.second);
}
                           __INLINE                          INMOST_DATA_REAL_TYPE get_table(INMOST_DATA_REAL_TYPE Arg, const INMOST::keyval_table & Table) {return Table.GetValue(Arg);}



#endif //USE_AUTODIFF
#endif //INMOST_AUTODIFF_ETEXPR_H_INCLUDED
