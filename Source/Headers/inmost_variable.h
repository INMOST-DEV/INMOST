
#ifndef INMOST_AUTODIFF_ETVAR_H_INCLUDED
#define INMOST_AUTODIFF_ETVAR_H_INCLUDED
#include "inmost_common.h"
#include "inmost_expression.h"
#include "inmost_mesh.h"
#include "inmost_autodiff.h"
#include "inmost_solver.h"
#include <sstream> //for debug
#include <new>

#if defined(USE_AUTODIFF) && (!defined(USE_MESH) || !defined(USE_SOLVER))
#warning "USE_AUTODIFF require USE_MESH"
#undef USE_AUTODIFF
#endif


//TODO:
// 1. Incorporate tables
// 2. (ok, test) implement condition
// 3. (ok, test) implement stencil
// 4. (???) copying of basic_dynamic_variable
// 5. Consider optimization by checking zero variation multipliers, check that assembly do not degrade.
// 6. Document everything



#if defined(USE_AUTODIFF)
namespace INMOST
{
  
  template<class Op, class A>
  class unary_pool
  {
    A arg;
    Op operand;
  public:
    unary_pool(const A & parg) : arg(parg), operand(arg) {}
    unary_pool(const unary_pool & other) : arg(other.arg), operand(arg) {}
    unary_pool & operator = (unary_pool const & other) {arg = other.arg; operand = Op(arg); return * this;}
    const shell_expression<A> & get_arg() {return arg;}
    Op & get_op() {return operand;}
    const Op & get_op() const {return operand;}
  };


  template<class Op, class A>
  class unary_const_pool
  {
    A left;
    INMOST_DATA_REAL_TYPE right;
    Op operand;
  public:
    unary_const_pool(const A & pleft, INMOST_DATA_REAL_TYPE pright) : left(pleft), right(pright), operand(left,right) {}
    unary_const_pool(const unary_const_pool & other) : left(other.left), right(other.right), operand(left,right) {}
    unary_const_pool & operator = (unary_const_pool const & other) {left = other.left; right = other.right; operand = Op(left,right); return * this;}
    const shell_expression<A> & get_arg() {return left;}
    Op & get_op() {return operand;}
    const Op & get_op() const {return operand;}
  };

  template<class Op, class A, class B>
  class binary_pool
  {
    
    A left;
    B right;
    Op operand;
  public:
    binary_pool(const A & pleft, const B & pright) : left(pleft), right(pright), operand(left,right) {}
    binary_pool(const binary_pool & other) : left(other.left), right(other.right), operand(left,right) {}
    binary_pool & operator =(binary_pool const & other) {left = other.left; right = other.right; operand = Op(left,right); return * this;}
    const shell_expression<A> & get_left() {return left;}
    const shell_expression<B> & get_right() {return right;}
    Op & get_op() {return operand;}
    const Op & get_op() const {return operand;}
    ~binary_pool() {}
  };

  template<class Op, class A, class B, class C>
  class ternary_pool
  {
    A cond;
    B left;
    C right;
    Op operand;
  public:
    ternary_pool(const A & pcond, const B & pleft, const C & pright) : cond(pcond), left(pleft), right(pright), operand(cond,left,right) {}
    ternary_pool(const ternary_pool & other) : cond(other.cond), left(other.left), right(other.right), operand(cond,left,right) {}
    ternary_pool & operator =(ternary_pool const & other) {cond = other.cond; left = other.left; right = other.right; operand = Op(cond,left,right); return * this;}
    const shell_expression<A> & get_cond() {return cond;}
    const shell_expression<B> & get_left() {return left;}
    const shell_expression<C> & get_right() {return right;}
    Op & get_op() {return operand;}
    const Op & get_op() const {return operand;}
    ~ternary_pool() {}
  };

  template<class A, class ArgA>
  class unary_pool_expression : public shell_expression<unary_pool_expression<A,ArgA> >
  {
    unary_pool<A,ArgA> pool;
   public:
    unary_pool_expression(const unary_pool<A,ArgA> & ppool) : pool(ppool) {}
    unary_pool_expression(const unary_pool_expression & other) : pool(other.pool) {}
    unary_pool_expression & operator = (unary_pool_expression const & other) {pool = other.pool; return * this;}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return pool.get_op().GetValue(); }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {pool.get_op().GetDerivative(mult,r);}
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {pool.get_op().GetDerivative(mult,r);}
  };

  template<class A, class ArgA>
  class unary_const_pool_expression : public shell_expression<unary_const_pool_expression<A,ArgA> >
  {
    unary_const_pool<A,ArgA> pool;
   public:
    unary_const_pool_expression(const unary_const_pool<A,ArgA> & ppool) : pool(ppool) {}
    unary_const_pool_expression(const unary_const_pool_expression & other) : pool(other.pool) {}
    unary_const_pool_expression & operator = (unary_const_pool_expression const & other) {pool = other.pool; return * this;}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return pool.get_op().GetValue(); }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {pool.get_op().GetDerivative(mult,r);}
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {pool.get_op().GetDerivative(mult,r);}
  };

  template<class A, class ArgA, class ArgB>
  class binary_pool_expression : public shell_expression<binary_pool_expression<A,ArgA,ArgB> >
  {
    binary_pool<A,ArgA,ArgB> pool;
   public:
    binary_pool_expression(const binary_pool<A,ArgA,ArgB> & ppool) : pool(ppool) {}
    binary_pool_expression(const binary_pool_expression & other) : pool(other.pool) {}
    binary_pool_expression & operator = (binary_pool_expression const & other) {pool = other.pool; return * this;}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return pool.get_op().GetValue(); }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {pool.get_op().GetDerivative(mult,r);}
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {pool.get_op().GetDerivative(mult,r);}
  };

  template<class A, class ArgA, class ArgB, class ArgC>
  class ternary_pool_expression : public shell_expression<ternary_pool_expression<A,ArgA,ArgB,ArgC> >
  {
    ternary_pool<A,ArgA,ArgB,ArgC> pool;
   public:
    ternary_pool_expression(const ternary_pool<A,ArgA,ArgB,ArgC> & ppool) : pool(ppool) {}
    ternary_pool_expression(const ternary_pool_expression & other) : pool(other.pool) {}
    ternary_pool_expression & operator = (ternary_pool_expression const & other) {pool = other.pool; return * this;}
    __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return pool.get_op().GetValue(); }
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {pool.get_op().GetDerivative(mult,r);}
    __INLINE void GetDerivative(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {pool.get_op().GetDerivative(mult,r);}
  };

  class abstract_dynamic_variable
  {
  public:
    virtual INMOST_DATA_REAL_TYPE Value (const Storage & e) const = 0;
    virtual multivar_expression Variable(const Storage & e) const = 0;
    virtual void GetVariation(const Storage & e, Sparse::Row & r) const = 0;
    virtual void GetVariation(const Storage & e, Sparse::RowMerger & r) const = 0;
    virtual abstract_dynamic_variable * Copy() const = 0;
  };

  template<typename RetType>
  class get_variable
  {
  public:
    virtual RetType operator()(const Storage & e) const = 0;
  };

  template<>
  class get_variable<multivar_expression>
  {
    const abstract_dynamic_variable & var;
  public:
    typedef multivar_expression type;
    get_variable(const abstract_dynamic_variable & var) : var(var) {}
    multivar_expression operator()(const Storage & e) const {return var.Variable(e);}
  };

  template<>
  class get_variable<INMOST_DATA_REAL_TYPE>
  {
    const abstract_dynamic_variable & var;
  public:
    typedef INMOST_DATA_REAL_TYPE type;
    get_variable(const abstract_dynamic_variable & var) : var(var) {}
    INMOST_DATA_REAL_TYPE operator()(const Storage & e) const {return var.Value(e);}
  };

  class store_variable
  {
    abstract_dynamic_variable * var;
  public:
    store_variable(const abstract_dynamic_variable & var) : var(var.Copy()) {}
    ~store_variable() {delete var;}
    template<typename T>
    get_variable<T> get_variable() {return get_variable<T>(*var);}
    abstract_dynamic_variable & retrive() {return *var;}
    const abstract_dynamic_variable & retrive() const {return *var;}
  };



  template<class VariableType>
  class basic_dynamic_variable : public abstract_dynamic_variable
  {
  public:
    typedef VariableType Var;
    virtual INMOST_DATA_REAL_TYPE Value(const Storage & e) const = 0;
    virtual multivar_expression Variable(const Storage & e) const = 0;
    virtual VariableType operator[](const Storage & e) const = 0;
    virtual void GetVariation(const Storage & e, Sparse::Row & r) const = 0;
    virtual void GetVariation(const Storage & e, Sparse::RowMerger & r) const = 0;
    virtual abstract_dynamic_variable * Copy() const = 0;
  };

  template<class VariableType, class Derived>
  class shell_dynamic_variable : virtual public basic_dynamic_variable<VariableType>
  {
  public:
    typedef VariableType Var;
    virtual INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return static_cast<const Derived *>(this)->Value(e);}
    virtual multivar_expression operator ()(const Storage & e) const {return static_cast<const Derived *>(this)->Variable(e);}
    virtual VariableType operator[](const Storage & e) const {return (*static_cast<const Derived *>(this))[e];}
    virtual void GetVariation(const Storage & e, Sparse::Row & r) const {static_cast<const Derived *>(this)->GetVariation(e,r);}
    virtual void GetVariation(const Storage & e, Sparse::RowMerger & r) const {static_cast<const Derived *>(this)->GetVariation(e,r);}
    operator Derived & () {return *static_cast<Derived *>(this);}
    operator const Derived & () const {return *static_cast<const Derived *>(this);}
    virtual abstract_dynamic_variable * Copy() const { return static_cast<const Derived *>(this)->Copy(); }
  };
	
  
  class dynamic_variable : public shell_dynamic_variable<var_expression,dynamic_variable >
  {
  private:
    Automatizator & aut;
    Tag index_tag, value_tag;
    MarkerType mask;
    INMOST_DATA_ENUM_TYPE comp;
  public:
    dynamic_variable(Automatizator & paut, INMOST_DATA_ENUM_TYPE pregval, INMOST_DATA_ENUM_TYPE pcomp = 0) : aut(paut), comp(pcomp) 
    {
      if( pregval != ENUMUNDEF )
      {
        mask = aut.GetDynamicMask(pregval);
        if( pregval < AD_CTAG )
        {
          value_tag = aut.GetDynamicValueTag(pregval);
          index_tag = aut.GetDynamicIndexTag(pregval);
        }
        else
        {
          value_tag = aut.GetStaticValueTag(pregval);
          index_tag = Tag();
        }
      }
    }
    dynamic_variable(const dynamic_variable & other) : aut(other.aut), index_tag(other.index_tag), value_tag(other.value_tag), mask(other.mask), comp(other.comp) {}
    dynamic_variable & operator =(const dynamic_variable & other) 
    {
      aut = other.aut; 
      index_tag = other.index_tag; 
      value_tag = other.value_tag;
      mask = other.mask;
      comp = other.comp;
      return * this;
    }
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return e->RealArray(value_tag)[comp];}
    INMOST_DATA_ENUM_TYPE Index(const Storage & e) const {return (!mask || e->GetMarker(mask))?e->IntegerArray(index_tag)[comp]:ENUMUNDEF;}
    multivar_expression Variable(const Storage & e) const 
    {
      if( !mask || e->GetMarker(mask) )
        return multivar_expression(e->RealArray(value_tag)[comp],e->IntegerArray(index_tag)[comp]);
      else
        return multivar_expression(e->Real(value_tag));
    }
    var_expression operator [](const Storage & e) const {return var_expression(e->Real(value_tag),(!mask || e->GetMarker(mask))?e->Integer(index_tag):ENUMUNDEF);}
    Tag IndexTag() {return index_tag;}
    Tag ValueTag() {return value_tag;}
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    bool isUnknown(const Storage & e) const {return (!mask || e->GetMarker(mask))?true:false;}
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new dynamic_variable(*this));}
  };

  
  class static_variable : public shell_dynamic_variable<const_expression,static_variable>
  {
  private:
    Tag value_tag;
    INMOST_DATA_ENUM_TYPE comp;
  public:
    static_variable(Tag t, INMOST_DATA_ENUM_TYPE pcomp = 0) : value_tag(t), comp(pcomp)  {}
    static_variable(const static_variable & other) : value_tag(other.value_tag), comp(other.comp) {}
    static_variable & operator =(const static_variable & other) 
    {
      value_tag = other.value_tag;
      comp = other.comp;
      return * this;
    }
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return e->RealArray(value_tag)[comp];}
    multivar_expression Variable(const Storage & e) const 
    {
      return multivar_expression(e->RealArray(value_tag)[comp]);
    }
    const_expression operator [](const Storage & e) const {return const_expression(e->RealArray(value_tag)[comp]);}
    Tag ValueTag() {return value_tag;}
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    bool isUnknown(const Storage & e) const {return false;}
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new static_variable(*this));}
  };

  class stored_variable : public shell_dynamic_variable<multivar_expression,stored_variable>
  {
  private:
    Tag variable_tag;
    INMOST_DATA_ENUM_TYPE comp;
  public:
    stored_variable(Tag t, INMOST_DATA_ENUM_TYPE pcomp = 0) : variable_tag(t), comp(pcomp)  {}
    stored_variable(const stored_variable & other) : variable_tag(other.variable_tag), comp(other.comp) {}
    stored_variable & operator =(const stored_variable & other) 
    {
      variable_tag = other.variable_tag;
      comp = other.comp;
      return * this;
    }
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return e->VariableArray(variable_tag)[comp].GetValue();}
    multivar_expression Variable(const Storage & e) const 
    {
      return e->VariableArray(variable_tag)[comp];
    }
    multivar_expression operator [](const Storage & e) const {return e->VariableArray(variable_tag)[comp];}
    Tag VariableTag() {return variable_tag;}
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    bool isUnknown(const Storage & e) const {return false;}
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new stored_variable(*this));}
  };

  template<class A>
  class stencil_variable : public shell_dynamic_variable< stencil_expression<typename A::Var>, stencil_variable<A> >
  {
  private:
    Automatizator & aut;
    INMOST_DATA_ENUM_TYPE stnclind;
    const A & Arg;
    void * user_data;
  public:
    stencil_variable(Automatizator & paut, INMOST_DATA_ENUM_TYPE pstnclind, const shell_dynamic_variable<typename A::Var,A> & parg, void * puser_data = NULL) : aut(paut), stnclind(pstnclind), Arg(parg), user_data(puser_data) {}
    stencil_variable(const stencil_variable & other) : aut(other.aut), stnclind(other.stnclind), Arg(other.Arg), user_data(other.user_data) {}
    stencil_variable & operator =(const stencil_variable & other) {aut = other.aut; stnclind = other.stnclind; Arg = other.Arg; user_data = other.user_data; return * this;}
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
    multivar_expression Variable(const Storage & e) const
    {
      multivar_expression ret = (*this)[e];
      return ret;
    }
    stencil_expression<typename A::Var> operator [](const Storage & e) const
    {
      dynarray<std::pair<INMOST_DATA_REAL_TYPE, typename A::Var>, 64> tmp;
      Automatizator::stencil_pairs stncl;
			aut.GetStencil(stnclind, e, user_data, stncl);
      tmp.resize(stncl.size());
      for(INMOST_DATA_ENUM_TYPE k = 0; k < stncl.size(); ++k)
        tmp[k] = std::make_pair(stncl[k].first, Arg[stncl[k].second]);
      return stencil_expression<typename A::Var>(tmp);
    }
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new stencil_variable(*this));}
  };

  template<class A>
  class table_variable : public shell_dynamic_variable< unary_pool_expression< function_expression<typename A::Var>,typename A::Var > , table_variable<A> >
  {
    const A & Arg;
    const keyval_table & Table;
  public:
    table_variable(const shell_dynamic_variable<typename A::Var,A> & parg, const keyval_table & ptable) : Arg(parg), Table(ptable) {}
    table_variable(const table_variable & other) : Arg(other.Arg), Table(other.table) {}
    table_variable & operator = (table_variable const & other) {Arg = other.Arg; Table = other.Table; return * this;}
    multivar_expression Variable(const Storage & e) const
    {
      multivar_expression ret = (*this)[e];
      return ret;
    }
    unary_pool_expression< function_expression<typename A::Var> ,typename A::Var > operator [](const Storage & e) const
    {
      typename A::Var arg = Arg[e];
      unary_pool< function_expression<typename A::Var>, typename A::Var> pool(arg);
      std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> both = Table.GetBoth(arg.GetValue());
      pool.get_op().SetFunctionValue(both.first);
      pool.get_op().SetFunctionDerivative(both.second);
      return unary_pool_expression< function_expression<typename A::Var>, typename A::Var >(pool);
    }
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new table_variable(*this));}
  };

  template<class Expr, class A>
  class unary_custom_variable : public shell_dynamic_variable< unary_pool_expression<Expr, typename A::Var >,unary_custom_variable<Expr,A> >
  {
  private:
    const A & Arg;
  public:
    unary_custom_variable(const shell_dynamic_variable<typename A::Var,A> & parg) : Arg(parg) {}
    unary_custom_variable(const unary_custom_variable & other) : Arg(other.Arg) {}
    unary_custom_variable & operator =(unary_custom_variable const & other) {Arg = other.Arg; return * this;}
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
    multivar_expression Variable(const Storage & e) const
    {
      multivar_expression ret = (*this)[e];
      return ret;
    }
    unary_pool_expression<Expr, typename A::Var > operator [](const Storage & e) const 
    {
      unary_pool<Expr,typename A::Var> pool(Arg[e]);
      return unary_pool_expression<Expr, typename A::Var >(pool);
    }
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new unary_custom_variable(*this));}
  };


  template<class Expr, class A, class B>
  class binary_custom_variable : public shell_dynamic_variable< binary_pool_expression<Expr, typename A::Var, typename B::Var >,binary_custom_variable<Expr,A,B> >
  {
  private:
    const A & Left;
    const B & Right;
  public:
    binary_custom_variable(const shell_dynamic_variable<typename A::Var,A> & pleft, const shell_dynamic_variable<typename B::Var,B> & pright) 
    : Left(pleft), Right(pright) {}
    binary_custom_variable(const binary_custom_variable & other) : Left(other.Left), Right(other.Right) {}
    binary_custom_variable & operator =(binary_custom_variable const & other) {Left = other.Left; Right = other.Right; return * this;}
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
    multivar_expression Variable(const Storage & e) const
    {
      multivar_expression ret = (*this)[e];
      return ret;
    }
    binary_pool_expression<Expr, typename A::Var, typename B::Var > operator [](const Storage & e) const 
    {
      binary_pool<Expr,typename A::Var,typename B::Var> pool(Left[e],Right[e]);
      return binary_pool_expression<Expr, typename A::Var, typename B::Var >(pool);
    }
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new binary_custom_variable(*this));}
  };

  template<class Expr, class A>
  class unary_const_custom_variable : public shell_dynamic_variable< unary_const_pool_expression<Expr, typename A::Var >,unary_const_custom_variable<Expr,A> >
  {
  private:
    const A & Left;
    INMOST_DATA_REAL_TYPE Right;
  public:
    unary_const_custom_variable(const shell_dynamic_variable<typename A::Var,A> & pleft, INMOST_DATA_REAL_TYPE pright) 
    : Left(pleft), Right(pright) {}
    unary_const_custom_variable(const unary_const_custom_variable & other) : Left(other.Left), Right(other.Right) {}
    unary_const_custom_variable & operator =(unary_const_custom_variable const & other) {Left = other.Left; Right = other.Right; return * this;}
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
    multivar_expression Variable(const Storage & e) const
    {
      multivar_expression ret = (*this)[e];
      return ret;
    }
    unary_const_pool_expression<Expr, typename A::Var > operator [](const Storage & e) const 
    {
      unary_const_pool<Expr,typename A::Var> pool(Left[e],Right);
      return unary_const_pool_expression<Expr, typename A::Var >(pool);
    }
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new unary_const_custom_variable(*this));}
  };

  template<class Expr, class A, class B, class C>
  class ternary_custom_variable : public shell_dynamic_variable< ternary_pool_expression<Expr, typename A::Var, typename B::Var, typename C::Var >,ternary_custom_variable<Expr,A,B,C> >
  {
  private:
    const A & Cond;
    const B & Left;
    const C & Right;
  public:
    ternary_custom_variable(const shell_dynamic_variable<typename A::Var,A> & pcond, const shell_dynamic_variable<typename B::Var,B> & pleft, const shell_dynamic_variable<typename C::Var,C> & pright) 
    : Cond(pcond), Left(pleft), Right(pright) {}
    ternary_custom_variable(const ternary_custom_variable & other) : Left(other.Left), Right(other.Right) {}
    ternary_custom_variable & operator =(ternary_custom_variable const & other) {Left = other.Left; Right = other.Right; return * this;}
    INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
    multivar_expression Variable(const Storage & e) const
    {
      multivar_expression ret = (*this)[e];
      return ret;
    }
    ternary_pool_expression<Expr, typename A::Var, typename B::Var, typename C::Var > operator [](const Storage & e) const 
    {
      ternary_pool<Expr,typename A::Var,typename B::Var, typename C::Var> pool(Cond[e],Left[e],Right[e]);
      return ternary_pool_expression<Expr, typename A::Var, typename B::Var, typename C::Var>(pool);
    }
    void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetDerivative(1.0,r); }
    void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetDerivative(1.0,r); }
    abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new ternary_custom_variable(*this));}
  };
  typedef abstract_dynamic_variable abstract_variable;
}

template<class A, class B, class C> __INLINE INMOST::ternary_custom_variable<INMOST::condition_expression<typename A::Var, typename B::Var, typename C::Var>,A,B,C> condition(INMOST::shell_dynamic_variable<typename A::Var, A> const & control, INMOST::shell_dynamic_variable<typename B::Var, B> const & if_ge_zero, INMOST::shell_dynamic_variable<typename C::Var, C> const & if_lt_zero) { return INMOST::ternary_custom_variable<INMOST::condition_expression<typename A::Var, typename B::Var, typename C::Var>,A,B,C>(control,if_ge_zero,if_lt_zero); }
template<class A>          __INLINE                                 INMOST::unary_custom_variable<INMOST::unary_minus_expression<typename A::Var>,A> operator-(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::unary_minus_expression<typename A::Var>,A>(Arg); }
template<class A>          __INLINE                                         INMOST::unary_custom_variable<INMOST::abs_expression<typename A::Var>,A>      fabs(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::abs_expression<typename A::Var>,A>(Arg); }
template<class A>          __INLINE                                         INMOST::unary_custom_variable<INMOST::exp_expression<typename A::Var>,A>       exp(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::exp_expression<typename A::Var>,A>(Arg); }
template<class A>          __INLINE                                         INMOST::unary_custom_variable<INMOST::log_expression<typename A::Var>,A>       log(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::log_expression<typename A::Var>,A>(Arg); }
template<class A>          __INLINE                                         INMOST::unary_custom_variable<INMOST::sin_expression<typename A::Var>,A>       sin(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::sin_expression<typename A::Var>,A>(Arg ); }
template<class A>          __INLINE                                         INMOST::unary_custom_variable<INMOST::cos_expression<typename A::Var>,A>       cos(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::cos_expression<typename A::Var>,A>(Arg); }
template<class A>          __INLINE                                        INMOST::unary_custom_variable<INMOST::sqrt_expression<typename A::Var>,A>      sqrt(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::sqrt_expression<typename A::Var>,A>(Arg); }
template<class A>          __INLINE              INMOST::unary_const_custom_variable<INMOST::variation_multiplication_expression<typename A::Var>,A> variation(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg, INMOST_DATA_REAL_TYPE Mult) {return INMOST::unary_const_custom_variable<INMOST::variation_multiplication_expression<typename A::Var>,A>(Arg,Mult);}
template<class A, class B> __INLINE                INMOST::binary_custom_variable<INMOST::addition_expression<typename A::Var,typename B::Var>,A, B> operator+(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::binary_custom_variable<INMOST::addition_expression<typename A::Var,typename B::Var>,A, B> (Left, Right); }
template<class A, class B> __INLINE             INMOST::binary_custom_variable<INMOST::subtraction_expression<typename A::Var,typename B::Var>,A, B> operator-(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::binary_custom_variable<INMOST::subtraction_expression<typename A::Var,typename B::Var>, A, B> (Left, Right); }
template<class A, class B> __INLINE          INMOST::binary_custom_variable<INMOST::multiplication_expression<typename A::Var,typename B::Var>,A, B> operator*(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::binary_custom_variable<INMOST::multiplication_expression<typename A::Var,typename B::Var>, A, B> (Left, Right); }
template<class A, class B> __INLINE                INMOST::binary_custom_variable<INMOST::division_expression<typename A::Var,typename B::Var>,A, B> operator/(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::binary_custom_variable<INMOST::division_expression<typename A::Var,typename B::Var>, A, B> (Left, Right); }
template<class A, class B> __INLINE                     INMOST::binary_custom_variable<INMOST::pow_expression<typename A::Var,typename B::Var>,A, B>       pow(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::binary_custom_variable<INMOST::pow_expression<typename A::Var,typename B::Var>,A, B>(Left, Right); }
template<class B>          __INLINE                             INMOST::unary_const_custom_variable<INMOST::const_pow_expression<typename B::Var>,B>       pow(INMOST_DATA_REAL_TYPE Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::unary_const_custom_variable<INMOST::const_pow_expression<typename B::Var>,B>(Left, Right); }
template<class A>          __INLINE                             INMOST::unary_const_custom_variable<INMOST::pow_const_expression<typename A::Var>,A>       pow(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::unary_const_custom_variable<INMOST::pow_const_expression<typename A::Var>,A>(Left, Right); }
template<class B>          __INLINE                  INMOST::unary_const_custom_variable<INMOST::const_multiplication_expression<typename B::Var>,B> operator*(INMOST_DATA_REAL_TYPE Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::unary_const_custom_variable<INMOST::const_multiplication_expression<typename B::Var>,B>(Right,Left); }
template<class A>          __INLINE                  INMOST::unary_const_custom_variable<INMOST::const_multiplication_expression<typename A::Var>,A> operator*(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::unary_const_custom_variable<INMOST::const_multiplication_expression<typename A::Var>,A>(Left,Right); }
template<class B>          __INLINE                            INMOST::unary_const_custom_variable<INMOST::reciprocal_expression<typename B::Var>,B> operator/(INMOST_DATA_REAL_TYPE Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::unary_const_custom_variable<INMOST::reciprocal_expression<typename B::Var>,B>(Right,Left); }
template<class A>          __INLINE                        INMOST::unary_const_custom_variable<INMOST::const_division_expression<typename A::Var>,A> operator/(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::unary_const_custom_variable<INMOST::const_division_expression<typename A::Var>,A>(Left, Right); }
template<class B>          __INLINE                        INMOST::unary_const_custom_variable<INMOST::const_addition_expression<typename B::Var>,B> operator+(INMOST_DATA_REAL_TYPE Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::unary_const_custom_variable<INMOST::const_addition_expression<typename B::Var>,B>(Right,Left); }
template<class A>          __INLINE                        INMOST::unary_const_custom_variable<INMOST::const_addition_expression<typename A::Var>,A> operator+(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::unary_const_custom_variable<INMOST::const_addition_expression<typename A::Var>,A>(Left,Right); }
template<class B>          __INLINE                     INMOST::unary_const_custom_variable<INMOST::const_subtraction_expression<typename B::Var>,B> operator-(INMOST_DATA_REAL_TYPE Left, INMOST::shell_dynamic_variable<typename B::Var,B> const & Right) { return INMOST::unary_const_custom_variable<INMOST::const_subtraction_expression<typename B::Var>,B>(Right, Left); }
template<class A>          __INLINE                        INMOST::unary_const_custom_variable<INMOST::const_addition_expression<typename A::Var>,A> operator-(INMOST::shell_dynamic_variable<typename A::Var,A> const & Left, INMOST_DATA_REAL_TYPE Right) { return INMOST::unary_const_custom_variable<INMOST::const_addition_expression<typename A::Var>,A>(Left, -Right); }
template<class A>          __INLINE                                                                                      INMOST::stencil_variable<A>   stencil(INMOST::Automatizator & aut, INMOST_DATA_ENUM_TYPE stncl, INMOST::shell_dynamic_variable<typename A::Var,A> const & Arg, void * user_data = NULL) { return INMOST::stencil_variable<A>(aut,stncl,Arg,user_data); }
template<class A>          __INLINE                                                                                        INMOST::table_variable<A> get_table(INMOST::shell_dynamic_variable<typename A::Var,A> const & Arg, const INMOST::keyval_table & Table) {return INMOST::table_variable<A>(Table,Arg);}

  
  


#endif
#endif
