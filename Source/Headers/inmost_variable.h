
#ifndef INMOST_AUTODIFF_ETVAR_H_INCLUDED
#define INMOST_AUTODIFF_ETVAR_H_INCLUDED
#include "inmost_common.h"
#include "inmost_expression.h"
#include "inmost_mesh.h"
#include "inmost_autodiff.h"
#include "inmost_solver.h"
#include <sstream> //for debug
#include <new>

#if defined(USE_AUTODIFF) && defined(USE_MESH)

//TODO:
// 1. Incorporate tables
// 2. (ok, test) implement condition
// 3. (ok, test) implement stencil
// 4. (???) copying of basic_dynamic_variable
// 5. Consider optimization by checking zero variation multipliers, check that assembly do not degrade.
// 6. Document everything
// 7. change stencil_variable with foreach_variable and introduce function foreach(iterator beg, iterator end, arg)
// 8. enclose in namespace
// 9. maybe should not use const A & in classes, since some class may be destroyed prior use - investigate

//This should stop Visual Studio from complaining of very long auto-generated class types
#ifdef _MSC_VER
#pragma warning(disable : 4503)
#endif



namespace INMOST
{
	
	template<class Op, class A>
	class unary_pool
	{
		A arg;
		Op operand;
		
		unary_pool & operator = (unary_pool const & other) {arg = other.arg; operand.assign(other.operand,arg); return * this;}
	public:
		unary_pool(const A & parg) : arg(parg), operand(arg) {}
		unary_pool(const A & parg, INMOST_DATA_REAL_TYPE pmult) : arg(parg), operand(arg,pmult) {}		
		unary_pool(const unary_pool & other) : arg(other.arg), operand(other.operand,arg) {}
		
		const shell_expression<A> & get_arg() {return arg;}
		Op & get_op() {return operand;}
		const Op & get_op() const {return operand;}
	};
	

	
	template<class Op, class A, class B>
	class binary_pool
	{
		
		A left;
		B right;
		Op operand;
		
		binary_pool & operator = (binary_pool const & other) {left = other.left; right = other.right; operand.assign(other.operand,left,right); return * this;}
	public:
		binary_pool(const A & pleft, const B & pright) : left(pleft), right(pright), operand(left,right) {}
		binary_pool(const binary_pool & other) : left(other.left), right(other.right), operand(other.operand,left,right) {}
		
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
		
		ternary_pool & operator =(ternary_pool const & other) {cond = other.cond; left = other.left; right = other.right; operand.assign(other.operand,cond,left,right); return * this;}
	public:
		ternary_pool(const A & pcond, const B & pleft, const C & pright) : cond(pcond), left(pleft), right(pright), operand(cond,left,right) {}
		ternary_pool(const ternary_pool & other) : cond(other.cond), left(other.left), right(other.right), operand(other.operand,cond,left,right) {}
		
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {pool.get_op().GetJacobian(mult,r);}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {pool.get_op().GetJacobian(mult,r);}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {pool.get_op().GetHessian(multJ,J,multH,H);}
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {pool.get_op().GetJacobian(mult,r);}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {pool.get_op().GetJacobian(mult,r);}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {pool.get_op().GetHessian(multJ,J,multH,H);}
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
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {pool.get_op().GetJacobian(mult,r);}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const {pool.get_op().GetJacobian(mult,r);}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {pool.get_op().GetHessian(multJ,J,multH,H);}
	};
	
	class abstract_dynamic_variable
	{
	public:
		virtual INMOST_DATA_REAL_TYPE Value (const Storage & e) const = 0;
		virtual multivar_expression Variable(const Storage & e) const = 0;
		virtual void GetVariation(const Storage & e, Sparse::Row & r) const = 0;
		virtual void GetVariation(const Storage & e, Sparse::RowMerger & r) const = 0;
		virtual abstract_dynamic_variable * Copy() const = 0;
		virtual ~abstract_dynamic_variable() {}
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
		virtual ~basic_dynamic_variable() {}
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
	
	class stored_variable_expression : public shell_dynamic_variable<multivar_expression,stored_variable_expression>
	{
		abstract_dynamic_variable * var;
	public:
		stored_variable_expression() : var(NULL) {}
		stored_variable_expression(const abstract_dynamic_variable & pvar) : var(pvar.Copy()) {}
		stored_variable_expression(const stored_variable_expression & other) : var(other.var->Copy()) {}
		~stored_variable_expression() {delete var; var = NULL;}
		stored_variable_expression operator =(stored_variable_expression const & other) {var = other.var->Copy(); return *this;}
		stored_variable_expression operator =(const abstract_dynamic_variable & pvar) {var = pvar.Copy(); return *this;}
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return var->Value(e);}
		multivar_expression Variable(const Storage & e) const {return var->Variable(e);}
		multivar_expression operator [](const Storage & e) const {return var->Variable(e);}
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		
		template<typename T>
		get_variable<T> get_variable() {return get_variable<T>(*var);}
		abstract_dynamic_variable & retrive_expression() {return *var;}
		const abstract_dynamic_variable & retrive_expression() const {return *var;}
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new stored_variable_expression(*this));}
	};
	
	
	class dynamic_variable : public shell_dynamic_variable<var_expression,dynamic_variable >
	{
	private:
		const AbstractEntry * entry;
		INMOST_DATA_ENUM_TYPE comp;
	public:
		dynamic_variable() :entry(NULL), comp(ENUMUNDEF) {}
		dynamic_variable(Automatizator & aut, INMOST_DATA_ENUM_TYPE reg_index, INMOST_DATA_ENUM_TYPE comp = 0) : entry(reg_index==ENUMUNDEF?NULL:&aut.GetEntry(reg_index)), comp(comp) {}
		dynamic_variable(const AbstractEntry * re, INMOST_DATA_ENUM_TYPE comp = 0) : entry(re), comp(comp) {}
		dynamic_variable(const dynamic_variable & other) : entry(other.entry), comp(other.comp) {}
		dynamic_variable & operator =(const dynamic_variable & other)
		{
			entry = other.entry;
			comp = other.comp;
			return * this;
		}
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return entry->Value(e,comp);}
		INMOST_DATA_ENUM_TYPE Index(const Storage & e) const {return entry->isValid(e) ? entry->Index(e,comp):ENUMUNDEF;}
		multivar_expression Variable(const Storage & e) const
		{
			if( entry->isValid(e) )
				return entry->Unknown(e,comp);
			else
				return entry->Value(e,comp);
		}
		var_expression operator [](const Storage & e) const {return var_expression(entry->Value(e,comp),entry->isValid(e)?entry->Index(e,comp):ENUMUNDEF);}
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		bool isUnknown(const Storage & e) const {return entry->isValid(e)?true:false;}
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new dynamic_variable(*this));}
	};
	
	class const_variable : public shell_dynamic_variable<const_expression,const_variable>
	{
	private:
		INMOST_DATA_REAL_TYPE value;
	public:
		const_variable(INMOST_DATA_REAL_TYPE _value) : value(_value)  {}
		const_variable(const const_variable & other) : value(other.value) {}
		const_variable & operator =(const const_variable & other)
		{
			value = other.value;
			return * this;
		}
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return value;}
		multivar_expression Variable(const Storage & e) const
		{
			return multivar_expression(value);
		}
		const_expression operator [](const Storage & e) const {return const_expression(value);}
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new const_variable(*this));}
	};
	
	class const_link_variable : public shell_dynamic_variable<const_expression,const_link_variable>
	{
	private:
		const INMOST_DATA_REAL_TYPE * value;
	public:
		const_link_variable(const INMOST_DATA_REAL_TYPE * _value) : value(_value)  {}
		const_link_variable(const const_link_variable & other) : value(other.value) {}
		const_link_variable & operator =(const const_link_variable & other)
		{
			value = other.value;
			return * this;
		}
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return *value;}
		multivar_expression Variable(const Storage & e) const
		{
			return multivar_expression(*value);
		}
		const_expression operator [](const Storage & e) const {return const_expression(*value);}
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new const_link_variable(*this));}
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
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		bool isUnknown(const Storage & e) const {(void)e; return false;}
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new static_variable(*this));}
	};
	
	class stored_variable : public shell_dynamic_variable<multivar_expression_reference,stored_variable>
	{
	private:
		Tag variable_tag;
		INMOST_DATA_ENUM_TYPE comp;
	public:
		stored_variable() : variable_tag(), comp(ENUMUNDEF) {}
		stored_variable(Tag t, INMOST_DATA_ENUM_TYPE pcomp = 0) : variable_tag(t), comp(pcomp)
		{
			assert(t.GetDataType() == DATA_REAL || t.GetDataType() == DATA_VARIABLE);
		}
		stored_variable(const stored_variable & other) : variable_tag(other.variable_tag), comp(other.comp) {}
		stored_variable & operator =(const stored_variable & other)
		{
			variable_tag = other.variable_tag;
			comp = other.comp;
			return * this;
		}
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const
		{
			if( variable_tag.GetDataType() == DATA_VARIABLE )
				return e->VariableArray(variable_tag)[comp].GetValue();
			else if( variable_tag.GetDataType() == DATA_REAL )
				return e->RealArray(variable_tag)[comp];
			else throw NotImplemented;
		}
		multivar_expression Variable(const Storage & e) const
		{
			if( variable_tag.GetDataType() == DATA_VARIABLE )
				return e->VariableArray(variable_tag)[comp];
			else if( variable_tag.GetDataType() == DATA_REAL )
				return variable(e->RealArray(variable_tag)[comp]);
			else throw NotImplemented;
		}
		multivar_expression_reference operator [](const Storage & e) const
		{
			if( variable_tag.GetDataType() == DATA_VARIABLE )
				return e->VariableArray(variable_tag)[comp];
			else if( variable_tag.GetDataType() == DATA_REAL )
				return multivar_expression_reference(e->RealArray(variable_tag)[comp],NULL);
			else throw NotImplemented;
		}
		Tag VariableTag() {return variable_tag;}
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		bool isUnknown(const Storage & e) const {(void)e; return false;}
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new stored_variable(*this));}
	};
	
	template<class A>
	class stencil_expression : public shell_expression<stencil_expression<A> >
	{
	public:
		typedef const_multiplication_expression<A> argument;
		typedef unary_pool_expression< argument, A> pool;
		typedef dynarray< argument, 64 > container;
	private:
		container arg;
		INMOST_DATA_REAL_TYPE value;
	public:
		stencil_expression(const container & parg) : arg(parg)
		{
			value = 0.0;
			for(typename container::iterator it = arg.begin(); it != arg.end(); ++it)
				value += it->GetValue();
		}
		stencil_expression(const stencil_expression & other) : arg(other.arg), value(other.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			for(typename container::iterator it = arg.begin(); it != arg.end(); ++it)
				it->GetJacobian(mult,r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::Row & r) const
		{
			for(typename container::iterator it = arg.begin(); it != arg.end(); ++it)
				it->GetJacobian(mult,r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			Sparse::Row tmpJ, curJ;
			Sparse::HessianRow tmpH, curH;
			for(typename container::iterator it = arg.begin(); it != arg.end(); ++it)
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
	};
	
	template<class A>
	class stencil_variable : public shell_dynamic_variable< stencil_expression<typename A::Var>, stencil_variable<A> >
	{
	private:
		TagReferenceArray tag_elems;
		TagRealArray      tag_coefs;
		A Arg;
	public:
		stencil_variable(Tag tag_elems, Tag tag_coefs, const shell_dynamic_variable<typename A::Var,A> & parg) : tag_elems(tag_elems), tag_coefs(tag_coefs), Arg(parg) {}
		stencil_variable(const stencil_variable & other) : tag_elems(other.tag_elems), tag_coefs(other.tag_coefs), Arg(other.Arg) {}
		stencil_variable & operator =(const stencil_variable & other) {tag_elems = other.tag_elems; tag_coefs = other.tag_coefs; Arg = other.Arg; return * this;}
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
		multivar_expression Variable(const Storage & e) const
		{
			multivar_expression ret = (*this)[e];
			return ret;
		}
		stencil_expression<typename A::Var> operator [](const Storage & e) const
		{
			typename stencil_expression<typename A::Var>::container tmp;
			Storage::real_array      coefs = tag_coefs[e];
			Storage::reference_array elems = tag_elems[e];
			assert(coefs.size() == elems.size());
			tmp.resize(elems.size());
			for(INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k)
			{
				typename stencil_expression<typename A::Var>::argument arg(Arg[elems[k]],coefs[k]);
				typename stencil_expression<typename A::Var>::pool pool(arg);
				tmp[k] = pool;
			}
			return stencil_expression<typename A::Var>(tmp);
		}
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new stencil_variable(*this));}
	};
	
	template<class A>
	class table_variable : public shell_dynamic_variable< unary_pool_expression< function_expression<typename A::Var>,typename A::Var > , table_variable<A> >
	{
		A Arg;
		const keyval_table & Table;
	public:
		table_variable(const shell_dynamic_variable<typename A::Var,A> & parg,  const keyval_table  & ptable) : Arg(parg), Table(ptable) {}
		table_variable(const table_variable & other) : Arg(other.Arg), Table(other.Table) {}
		table_variable & operator = (table_variable const & other) {Arg = other.Arg; Table = other.Table; return * this;}
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
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
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new table_variable(*this));}
	};
	
	/// This class makes possible to evaluate different expressions on different element types.
	/// See etype_branch function.
	template<class A, class B>
	class etype_branch_variable : public shell_dynamic_variable< multivar_expression, etype_branch_variable<A,B> >
	{
	private:
		A ArgA; //< Variable expression to be evaluated when type of provided element matches selected types.
		B ArgB; //< Variable expression to be evaluated when type of provided element does not match selected types.
		ElementType types_true; //< Selected types of elements.
	public:
		/// Constructor. Used by etype_branch function.
		etype_branch_variable(ElementType _types_true, const A & _ArgA, const B & _ArgB) : types_true(_types_true), ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		etype_branch_variable(const etype_branch_variable & other) : types_true(other.types_true), ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		etype_branch_variable & operator =(etype_branch_variable const & other)
		{
			types_true = other.types_true;
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		multivar_expression Variable(const Storage & e) const { return (*this)[e]; }
		/// Build an expression associated with variable expression on provided element e.
		multivar_expression operator [](const Storage & e) const
		{
			if( e->GetElementType() & types_true )
				return ArgA[e];
			else return ArgB[e];
		}
		/// Retrive first derivatives of variable expression on provided element e, default approach.
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		/// Retrive first derivatives of variable expression on provided element e, with supplimentary structure Sparse::RowMerger.
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new etype_branch_variable(*this));}
	};
	
	/// This class makes possible to evaluate different expressions depending on the markers.
	/// Works similarly for shared and private markers.
	/// See marker_branch function.
	template<class A, class B>
	class marker_branch_variable : public shell_dynamic_variable< multivar_expression, marker_branch_variable<A,B> >
	{
	private:
		A ArgA; //< Variable expression to be evaluated when marker is set on the element.
		B ArgB; //< Variable expression to be evaluated when marker is not set on the element.
		MarkerType marker; //< Marker.
	public:
		/// Constructor. Used by marker_branch function.
		marker_branch_variable(MarkerType _marker, const A & _ArgA, const B & _ArgB) : marker(_marker), ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		marker_branch_variable(const marker_branch_variable & other) : marker(other.marker), ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		marker_branch_variable & operator =(marker_branch_variable const & other)
		{
			marker = other.marker;
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		INMOST_DATA_REAL_TYPE Value(const Storage & e) const {return (*this)[e].GetValue();}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		multivar_expression Variable(const Storage & e) const { return (*this)[e]; }
		/// Build an expression associated with variable expression on provided element e.
		multivar_expression operator [](const Storage & e) const
		{
			if( isPrivate(marker) ? e->GetPrivateMarker(marker) : e->GetMarker(marker) )
				return ArgA[e];
			else return ArgB[e];
		}
		/// Retrive first derivatives of variable expression on provided element e, default approach.
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		/// Retrive first derivatives of variable expression on provided element e, with supplimentary structure Sparse::RowMerger.
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new marker_branch_variable(*this));}
	};
	
	
	
	template<class Expr, class A>
	class unary_custom_variable : public shell_dynamic_variable< unary_pool_expression<Expr, typename A::Var >,unary_custom_variable<Expr,A> >
	{
	private:
		A Arg;
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
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new unary_custom_variable(*this));}
	};
	
	template<class Expr, class A>		
 	class unary_const_custom_variable : public shell_dynamic_variable< unary_pool_expression<Expr, typename A::Var >,unary_const_custom_variable<Expr,A> >		
 	{		
 	private:		
 		A Left;		
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
 		unary_pool_expression<Expr, typename A::Var > operator [](const Storage & e) const		
 		{		
 			unary_pool<Expr,typename A::Var> pool(Left[e],Right);		
 			return unary_pool_expression<Expr, typename A::Var >(pool);		
 		}		
 		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }		
 		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }		
 		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new unary_const_custom_variable(*this));}		
 	};
	
	template<class Expr, class A, class B>
	class binary_custom_variable : public shell_dynamic_variable< binary_pool_expression<Expr, typename A::Var, typename B::Var >,binary_custom_variable<Expr,A,B> >
	{
	private:
		A Left;
		B Right;
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
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new binary_custom_variable(*this));}
	};

	template<class Expr, class A, class B, class C>
	class ternary_custom_variable : public shell_dynamic_variable< ternary_pool_expression<Expr, typename A::Var, typename B::Var, typename C::Var >,ternary_custom_variable<Expr,A,B,C> >
	{
	private:
		A Cond;
		B Left;
		C Right;
	public:
		ternary_custom_variable(const shell_dynamic_variable<typename A::Var,A> & pcond, const shell_dynamic_variable<typename B::Var,B> & pleft, const shell_dynamic_variable<typename C::Var,C> & pright)
		: Cond(pcond), Left(pleft), Right(pright) {}
		ternary_custom_variable(const ternary_custom_variable & other) : Cond(other.Cond), Left(other.Left), Right(other.Right) {}
		ternary_custom_variable & operator =(ternary_custom_variable const & other) {Cond = other.Cond; Left = other.Left; Right = other.Right; return * this;}
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
		void GetVariation(const Storage & e, Sparse::Row & r) const { (*this)[e].GetJacobian(1.0,r); }
		void GetVariation(const Storage & e, Sparse::RowMerger & r) const { (*this)[e].GetJacobian(1.0,r); }
		abstract_dynamic_variable * Copy() const {return static_cast<abstract_dynamic_variable *>(new ternary_custom_variable(*this));}
	};
	typedef abstract_dynamic_variable abstract_variable;
}

template<class A, class B, class C> 
__INLINE 
INMOST::ternary_custom_variable<INMOST::condition_expression<typename A::Var, typename B::Var, typename C::Var>,A,B,C> condition(INMOST::shell_dynamic_variable<typename A::Var, A> const & control, INMOST::shell_dynamic_variable<typename B::Var, B> const & if_ge_zero, INMOST::shell_dynamic_variable<typename C::Var, C> const & if_lt_zero) { return INMOST::ternary_custom_variable<INMOST::condition_expression<typename A::Var, typename B::Var, typename C::Var>,A,B,C>(control,if_ge_zero,if_lt_zero); }
template<class A>          __INLINE                                 INMOST::unary_custom_variable<INMOST::unary_minus_expression<typename A::Var>,A> operator-(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::unary_minus_expression<typename A::Var>,A>(Arg); }
template<class A>          __INLINE                                  INMOST::unary_custom_variable<INMOST::unary_plus_expression<typename A::Var>,A> operator+(INMOST::shell_dynamic_variable<typename A::Var, A> const & Arg) { return INMOST::unary_custom_variable<INMOST::unary_plus_expression<typename A::Var>,A>(Arg); }
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
template<class A>          __INLINE                                                                                      INMOST::stencil_variable<A>   stencil(INMOST::Tag tag_elems, INMOST::Tag tag_coefs, INMOST::shell_dynamic_variable<typename A::Var,A> const & Arg) { return INMOST::stencil_variable<A>(tag_elems,tag_coefs,Arg); }
template<class A>          __INLINE                                                                                        INMOST::table_variable<A> get_table(INMOST::shell_dynamic_variable<typename A::Var,A> const & Arg, const INMOST::keyval_table & Table) {return INMOST::table_variable<A>(Arg,Table);}
template<class A>          __INLINE                                                                                    INMOST::stencil_expression<A>   stencil(INMOST::HandleType * elems, INMOST_DATA_REAL_TYPE * coefs, INMOST_DATA_ENUM_TYPE num, INMOST::shell_dynamic_variable<typename A::Var,A> const & Arg)
{
	INMOST::dynarray< INMOST::const_multiplication_expression<typename A::Var>, 64> tmp;
	for( INMOST_DATA_ENUM_TYPE k = 0; k < num; ++k) if( elems[k] != 0 )
		tmp.push_back(INMOST::const_multiplication_expression<typename A::Var>(Arg[elems[k]],coefs[k]));
	return INMOST::stencil_expression<typename A::Var>(tmp);
}
template<class A, class B> __INLINE INMOST::etype_branch_variable<A,B> etype_branch(INMOST::ElementType true_type, INMOST::shell_dynamic_variable<typename A::Var,A> const & iftrue, INMOST::shell_dynamic_variable<typename B::Var,B> const & iffalse) {return INMOST::etype_branch_variable<A,B>(true_type,iftrue,iffalse);}
template<class A, class B> __INLINE INMOST::marker_branch_variable<A,B> marker_branch(INMOST::MarkerType marker, INMOST::shell_dynamic_variable<typename A::Var,A> const & iftrue, INMOST::shell_dynamic_variable<typename B::Var,B> const & iffalse) {return INMOST::marker_branch_variable<A,B>(marker,iftrue,iffalse);}
__INLINE INMOST::const_link_variable extval(const INMOST_DATA_REAL_TYPE & pvar) {return INMOST::const_link_variable(&pvar);}

#endif //defined(USE_AUTODIFF) && defined(USE_MESH)




#endif //INMOST_AUTODIFF_ETVAR_H_INCLUDED

