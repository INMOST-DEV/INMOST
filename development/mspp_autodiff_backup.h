#pragma once
#ifndef INMOST_AUTODIFF_H_INCLUDED
#define INMOST_AUTODIFF_H_INCLUDED
#include "inmost_common.h"
#include "inmost_mesh.h"
#include <sstream> //for debug


#if defined(USE_AUTODIFF) && (!defined(USE_MESH))
#warning "USE_AUTODIFF require USE_MESH"
#undef USE_AUTODIFF
#endif

//#define DPRNT


//#define EXPRESSION_TEMPLATES

//TODO:


// Recheck computation of derivatives!!!!
// Detect similar parts of the tree in Automatizator
// Restructure tree expressions
// Intorduce multivariate tables
// Generate opencl code

// Make so that stencil may be represented by tags, by set or by callback_function

//RegisterTable and table implementation should account for boundary conditions beyond data

#if defined(USE_AUTODIFF)
#include <math.h>

#if defined(USE_AUTODIFF_ASMJIT)
#include "asmjit.h"
#endif

#define FILTER_EPS 1e-12

#define AD_SPACE 16384 //typical number of registered entries

#define AD_NONE  0 //this should generate an error
//binary operations below
#define AD_PLUS  1 //sum x+y
#define AD_MINUS 2 //substraction x-y
#define AD_MULT  3 //multiplication x*y
#define AD_DIV   4 //devision x/y
#define AD_POW   6 //power operation x^y
#define AD_COND  7 //condition
#define AD_ALTR  8 //alternatives for condition

#define AD_PRECOMP 15 // this expression points to precomputed replacement, expr * left is position in array of precomputed values, expr * right points to original expression

//unary operations below
#define AD_COS   20 //cos(x)
#define AD_INV   21 // 1/x
#define AD_ABS   22 // |x|
#define AD_EXP   23 // exp(x)
#define AD_LOG   24 // log(x)
#define AD_SIN   25 // sin(x)

#define AD_CONST 50 // expr * left represents const value
#define AD_MES   51 // volume(e) for cell, area(e) for face, length(e) for edge

#define AD_TAG   100               //tag with dynamic data
#define AD_CTAG  (AD_TAG+AD_SPACE)   //tag with constant data
#define AD_STNCL (AD_CTAG+AD_SPACE)  //stencil that sets current elements and their coefs
#define AD_TABLE (AD_STNCL+AD_SPACE) //table of values
#define AD_FUNC  (AD_TABLE+AD_SPACE) //register function that calculates some value from element

namespace INMOST
{
	class Automatizator; //forward declaration

	typedef dynarray<INMOST_DATA_REAL_TYPE, 2048> precomp_values_t;

#if defined(USE_AUTODIFF_EXPRESSION_TEMPLATES)
	typedef small_hash< INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE, 128> DerivativeArray;

	static void AddArray(DerivativeArray & target, DerivativeArray & op) { for (DerivativeArray::iterator k = op.begin(); k != op.end(); ++k) target[k->first] += k->second; }
	static void MultArray(DerivativeArray & target, INMOST_DATA_REAL_TYPE mult) {for (DerivativeArray::iterator k = target.begin(); k != target.end(); ++k) k->second *= mult;}


	class basic_expression
	{
	public:
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const = 0;
		__INLINE virtual INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const = 0;
		__INLINE virtual INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const = 0;
		__INLINE virtual void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const = 0;
		virtual ~basic_expression() {}
		virtual basic_expression * Copy() const = 0;
	};

	template<class Derived>
	class shell_expression : virtual public basic_expression
	{
	public:
		shell_expression() {}
		shell_expression(const shell_expression & other) {}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const {return static_cast<const Derived *>(this)->GetValue(elem, aut, user_data); }
		__INLINE virtual INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const { return static_cast<const Derived *>(this)->GetDerivative(elem, aut, user_data, out); }
		__INLINE virtual INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const { return static_cast<const Derived *>(this)->DerivativePrecompute(elem, aut, user_data, out); };
		__INLINE virtual void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{ static_cast<const Derived *>(this)->DerivativeFill(elem, aut, user_data, values, entries, multval); }
		virtual ~shell_expression() {}
		virtual basic_expression * Copy() const { return static_cast<const Derived*>(this)->Copy(); };
		operator Derived &() { return static_cast<Derived &>(*this); }
	};


	typedef basic_expression * basic_expression_ptr;
#endif
	class expr;
	
	

	class Automatizator
	{
	public:
		typedef INMOST_DATA_REAL_TYPE(*func_callback)(Storage * current_element, void * user_data);
		typedef std::pair<Storage *, INMOST_DATA_REAL_TYPE> stencil_pair;
		typedef dynarray<stencil_pair, 64> stencil_pairs;
		typedef void(*stencil_callback)(Storage * current_element, stencil_pairs & out_stencil, void * user_data);
	private:
		typedef struct{ Tag t; MIDType domain_mask; } tagdomain;
		typedef small_hash<INMOST_DATA_ENUM_TYPE, tagdomain, 128> const_tag_type;
		typedef struct{ tagdomain d; Tag indices; } tagpair;
		typedef small_hash<INMOST_DATA_ENUM_TYPE, tagpair, 128> tagpairs_type;
		typedef std::vector<tagpair> index_enum;
		typedef struct { Tag elements, coefs; } stencil_tag;
		typedef struct { std::string name; INMOST_DATA_ENUM_TYPE kind; MIDType domainmask; void * link; } stencil_kind_domain;
		typedef small_hash<INMOST_DATA_ENUM_TYPE, stencil_kind_domain, 128> stencil_type;
		typedef struct func_name_callback_t { std::string name; func_callback func; } func_name_callback;
		typedef small_hash<INMOST_DATA_ENUM_TYPE, func_name_callback, 128> func_type;
	public:
		typedef struct
		{
			std::string name;
			INMOST_DATA_REAL_TYPE * args;
			INMOST_DATA_REAL_TYPE * vals;
			INMOST_DATA_ENUM_TYPE size;
			__INLINE INMOST_DATA_ENUM_TYPE binary_search(INMOST_DATA_REAL_TYPE arg)
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
			__INLINE INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE arg)
			{
				if (arg < args[0]) return vals[0];
				INMOST_DATA_ENUM_TYPE i = binary_search(arg);
				return vals[i] + (vals[i + 1] - vals[i]) * (arg - args[i]) / (args[i + 1] - args[i]);
			}
			__INLINE INMOST_DATA_REAL_TYPE get_derivative(INMOST_DATA_REAL_TYPE arg)
			{
				if (arg < args[0]) return 0.0;
				INMOST_DATA_ENUM_TYPE i = binary_search(arg);
				return (vals[i + 1] - vals[i]) / (args[i + 1] - args[i]);
			}
			__INLINE std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> get_both(INMOST_DATA_REAL_TYPE arg)
			{
				if (arg < args[0]) return std::make_pair(vals[0], 0.0);
				INMOST_DATA_ENUM_TYPE i = binary_search(arg);
				INMOST_DATA_REAL_TYPE der = (vals[i + 1] - vals[i]) / (args[i + 1] - args[i]);
				return std::make_pair(vals[i] + der * (arg - args[i]), der);
			}
		} table;
		typedef table * table_ptr;
	private:
		typedef small_hash<INMOST_DATA_ENUM_TYPE, table_ptr, 32> table_type;
		const_tag_type reg_ctags;
		index_enum index_tags;
		tagpairs_type reg_tags;
		table_type reg_tables;
		stencil_type reg_stencils;
		func_type reg_funcs;
		INMOST_DATA_ENUM_TYPE first_num;
		INMOST_DATA_ENUM_TYPE last_num;
		Mesh * m;
		INMOST_DATA_REAL_TYPE DerivativePrecompute(const expr & var, Storage * e, precomp_values_t & values, void * user_data);
		void DerivativeFill(const expr & var, Storage * e, Solver::Row & entries, precomp_values_t & values, INMOST_DATA_REAL_TYPE multval, void * user_data);
	public:
		Automatizator(Mesh * m);
		~Automatizator();
		__INLINE INMOST_DATA_ENUM_TYPE GetFirstIndex() { return first_num; }
		__INLINE INMOST_DATA_ENUM_TYPE GetLastIndex() { return last_num; }
		INMOST_DATA_ENUM_TYPE RegisterFunc(std::string name, func_callback func);
		INMOST_DATA_ENUM_TYPE RegisterStencil(std::string name, Tag elements_tag, Tag coefs_tag, MIDType domain_mask = 0);
		INMOST_DATA_ENUM_TYPE RegisterStencil(std::string name, stencil_callback func, MIDType domain_mask = 0);
		INMOST_DATA_ENUM_TYPE RegisterTable(std::string name, INMOST_DATA_REAL_TYPE * Arguments, INMOST_DATA_REAL_TYPE * Values, INMOST_DATA_ENUM_TYPE size);
		INMOST_DATA_ENUM_TYPE RegisterDynamicTag(Tag t, ElementType typemask, MIDType domain_mask = 0);
		INMOST_DATA_ENUM_TYPE RegisterStaticTag(Tag t, MIDType domain_mask = 0);
		void EnumerateDynamicTags();
		__INLINE Tag                 GetDynamicValueTag(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind].d.t; }
		__INLINE Tag                 GetDynamicIndexTag(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind].indices; }
		__INLINE MIDType             GetDynamicMask(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind].d.domain_mask; }
		__INLINE Tag                 GetStaticValueTag(INMOST_DATA_ENUM_TYPE ind) { return reg_ctags[ind].t; }
		__INLINE MIDType             GetStaticMask(INMOST_DATA_ENUM_TYPE ind) { return reg_ctags[ind].domain_mask; }
		__INLINE INMOST_DATA_REAL_TYPE GetDynamicValue(Storage * e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->RealArray(GetDynamicValueTag(ind))[comp]; }
		__INLINE INMOST_DATA_ENUM_TYPE GetDynamicIndex(Storage * e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->IntegerArray(GetDynamicIndexTag(ind))[comp]; }
		__INLINE bool                isDynamicValid(Storage * e, INMOST_DATA_ENUM_TYPE ind) { MIDType mask = GetDynamicMask(ind); return mask == 0 || e->GetMarker(mask); }
		__INLINE INMOST_DATA_REAL_TYPE GetStaticValue(Storage * e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->RealArray(GetStaticValueTag(ind))[comp]; }
		__INLINE bool                isStaticValid(Storage * e, INMOST_DATA_ENUM_TYPE ind) { MIDType mask = GetStaticMask(ind); return mask == 0 || e->GetMarker(mask); }
#if defined(USE_AUTODIFF_EXPRESSION_TEMPLATES)
		__INLINE INMOST_DATA_REAL_TYPE Evaluate(basic_expression_ptr  var, Storage * e, void * user_data) { return var->GetValue(e, *this, user_data); }
		__INLINE INMOST_DATA_REAL_TYPE Derivative(basic_expression_ptr  var, Storage * e, Solver::Row & out, void * user_data) 
		{ 
			//DerivativeArray  inout;
			//INMOST_DATA_REAL_TYPE ret =  var->GetDerivative(e, *this, user_data, inout); 
			//for (DerivativeArray::iterator k = inout.begin(); k != inout.end(); ++k) out[k->first] += k->second;
			precomp_values_t vals;
			INMOST_DATA_REAL_TYPE ret = var->DerivativePrecompute(e, *this, user_data, vals);
			var->DerivativeFill(e,*this,user_data,vals,out,1.0);
			return ret;
		}
		static table_ptr GetInnerTable(Automatizator & aut, INMOST_DATA_ENUM_TYPE tableind) { return aut.reg_tables[tableind]; }
		static func_callback GetInnerFunction(Automatizator & aut, INMOST_DATA_ENUM_TYPE funcind) { return aut.reg_funcs[funcind].func; }
#endif
		INMOST_DATA_REAL_TYPE Evaluate(const expr & var, Storage * e, void * user_data);
		INMOST_DATA_REAL_TYPE Derivative(const expr & var, Storage * e, Solver::Row & out, void * user_data);
		__INLINE INMOST_DATA_REAL_TYPE                                 GetIndex(Storage * e, INMOST_DATA_ENUM_TYPE tagind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->IntegerArray(GetDynamicIndexTag(tagind))[comp]; }
		__INLINE INMOST_DATA_ENUM_TYPE                                 GetComponents(Storage *e, INMOST_DATA_ENUM_TYPE tagind) { return e->IntegerArray(GetDynamicIndexTag(tagind)).size(); }
		__INLINE Mesh *                                              GetMesh() { return m; }
		__INLINE INMOST_DATA_REAL_TYPE                                 GetTableValue(INMOST_DATA_ENUM_TYPE tableind, INMOST_DATA_REAL_TYPE arg) { return reg_tables[tableind]->get_value(arg); }
		__INLINE INMOST_DATA_REAL_TYPE                                 GetTableDerivative(INMOST_DATA_ENUM_TYPE tableind, INMOST_DATA_REAL_TYPE arg) { return reg_tables[tableind]->get_derivative(arg); }
		__INLINE std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> GetTableBoth(INMOST_DATA_ENUM_TYPE tableind, INMOST_DATA_REAL_TYPE arg) { return reg_tables[tableind]->get_both(arg); }
		__INLINE INMOST_DATA_ENUM_TYPE                                 GetStencil(INMOST_DATA_ENUM_TYPE stnclind, Storage * elem, void * user_data, stencil_pairs & ret)
		{
			stencil_kind_domain st = reg_stencils[stnclind];
			assert(st.domainmask == 0 || elem->GetMarker(st.domainmask));
			if (st.kind == 0)
			{
				Storage::reference_array elems = elem->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
				Storage::real_array coefs = elem->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
				assert(elems.size() == coefs.size());
				for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k) ret.push_back(std::make_pair(elems[k], coefs[k]));
			}
			else if (st.kind == 1) reinterpret_cast<stencil_callback>(st.link)(elem, ret, user_data);
			return ret.size();
		}
		__INLINE INMOST_DATA_REAL_TYPE GetFunction(INMOST_DATA_ENUM_TYPE funcid, Storage * elem, void * user_data) { return reg_funcs[funcid].func(elem, user_data); }
#if defined(USE_AUTODIFF_ASMJIT)
	private:
#if defined(USE_OMP)
#define MAX_NESTED 32
#define MAX_THREADS 64
#define THREAD_NUM omp_get_thread_num()
#else
#define MAX_NESTED 32
#define MAX_THREADS 1
#define THREAD_NUM 0
#endif
		struct local_storage_pairs
		{
			stencil_pairs pairs;
			char padding[4096 - sizeof(stencil_pairs)];
		} stencil_storage[MAX_THREADS][MAX_NESTED];
		precomp_values_t stacks[MAX_THREADS];
		asmjit::JitRuntime runtime;
		typedef dynarray<asmjit::host::XmmVar, 128> precomp_vars_t;
		asmjit::host::XmmVar          GenerateEvaluateSubASMJIT(const expr & var, asmjit::host::Compiler & comp, asmjit::host::GpVar & elem, asmjit::host::GpVar & user_data, INMOST_DATA_ENUM_TYPE nested);
		asmjit::host::XmmVar          GenerateDerivativePrecomputeSubASMJIT(const expr & var, asmjit::host::Compiler & comp, asmjit::host::GpVar & elem, asmjit::host::GpVar & user_data, INMOST_DATA_ENUM_TYPE nested);
		void                          GenerateDerivativeFillSubASMJIT(const expr & var, asmjit::host::Compiler & comp, asmjit::host::GpVar & elem, asmjit::host::GpVar & user_data, asmjit::host::GpVar & row, asmjit::host::XmmVar & multval, INMOST_DATA_ENUM_TYPE nested);
		static INMOST_DATA_REAL_TYPE    GetGeometricDataJIT(Mesh * m, Storage * e);
		static INMOST_DATA_REAL_TYPE    GetTableValueJIT(INMOST_DATA_REAL_TYPE arg, table_ptr tab);
		static INMOST_DATA_REAL_TYPE    GetTableDerivativeJIT(INMOST_DATA_REAL_TYPE arg, table_ptr tab);
		static INMOST_DATA_REAL_TYPE    GetValueJIT(Storage * e, Tag * t, INMOST_DATA_ENUM_TYPE comp);
		static INMOST_DATA_ENUM_TYPE    GetIndexJIT(Storage * e, Tag * t, INMOST_DATA_ENUM_TYPE comp);
		static INMOST_DATA_ENUM_TYPE    GetStencilJIT(Automatizator * aut, Storage * elem, void * user_data, INMOST_DATA_ENUM_TYPE id, INMOST_DATA_ENUM_TYPE nested);
		static INMOST_DATA_REAL_TYPE    GetStencilCoefJIT(Automatizator * aut, INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE nested);
		static Storage *              GetStencilElemJIT(Automatizator * aut, INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE nested);
		static void                   ClearStencilJIT(Automatizator * aut, INMOST_DATA_ENUM_TYPE id, INMOST_DATA_ENUM_TYPE nested);
		static void                   PushJIT(INMOST_DATA_REAL_TYPE val, Automatizator * aut);
		static INMOST_DATA_REAL_TYPE    PopJIT(Automatizator * aut);
		static INMOST_DATA_REAL_TYPE    InsertPairJIT(INMOST_DATA_REAL_TYPE val, INMOST_DATA_ENUM_TYPE ind, Solver::Row * row);
		void PushCode(asmjit::host::Compiler & comp, asmjit::host::XmmVar & var);
		void PopCode(asmjit::host::Compiler & comp, asmjit::host::XmmVar & var);
		void * GenerateEvaluateASMJIT(const expr & var);
		void * GenerateDerivativeASMJIT(const expr & var);
	public:
#if defined(DPRNT)
		static std::string opname(Automatizator * aut, INMOST_DATA_ENUM_TYPE op)
		{
			std::stringstream str;
			switch (op)
			{
			case AD_NONE: return "none";
			case AD_PLUS: return "add";
			case AD_MINUS: return "sub";
			case AD_MULT: return "mul";
			case AD_DIV: return "div";
			case AD_INV: return "rec";
			case AD_COND: return "cond";
			case AD_POW: return "pow";
			case AD_ABS: return "abs";
			case AD_EXP: return "exp";
			case AD_LOG: return "log";
			case AD_SIN: return "sin";
			case AD_COS: return "cos";
			case AD_CONST: return "const";
			case AD_MES: return "mes";
			}
			if (op >= AD_FUNC) { str << "func(" << aut->reg_funcs[op].name << ")"; return str.str(); }
			else if (op >= AD_TABLE) { str << "tab(" << aut->reg_tables[op]->name << ")"; return str.str(); }
			else if (op >= AD_STNCL) { str << "stncl(" << aut->reg_stencils[op].name << ")"; return str.str(); }
			else if (op >= AD_CTAG) { str << "const_tag(" << aut->GetStaticValueTag(op).GetTagName() << ")"; return str.str(); }
			else if (op >= AD_TAG) { str << "tag(" << aut->GetDynamicValueTag(op).GetTagName() << ")"; return str.str(); }
			assert(false);
			return "unknown";
		}
		static void PrintOp0(Automatizator * aut, INMOST_DATA_ENUM_TYPE op) { std::cout << opname(aut, op) << std::endl; }
		static INMOST_DATA_REAL_TYPE PrintOp1(INMOST_DATA_REAL_TYPE arg, Automatizator * aut, INMOST_DATA_ENUM_TYPE op) { std::cout << opname(aut, op) << " param: " << arg << std::endl; return arg; }
		static INMOST_DATA_REAL_TYPE PrintConst(INMOST_DATA_REAL_TYPE coef) { std::cout << " mult: " << coef << std::endl; return coef; }
		static INMOST_DATA_REAL_TYPE PrintRet(INMOST_DATA_REAL_TYPE ret, Automatizator * aut, INMOST_DATA_ENUM_TYPE op) { std::cout << opname(aut, op) << " return: " << ret << std::endl; return ret; }
#endif
		void * GenerateCode(const expr & var);
		INMOST_DATA_REAL_TYPE Evaluate(void * prog, Storage * elem, void * user_data);
		INMOST_DATA_REAL_TYPE Derivative(void * prog, Storage * elem, Solver::Row & out, void * user_data);
#endif
	};


#if defined(USE_AUTODIFF_EXPRESSION_TEMPLATES)
	
	

	class const_expression : public shell_expression<const_expression>
	{
		INMOST_DATA_REAL_TYPE value;
	public:
		const_expression(INMOST_DATA_REAL_TYPE value) : value(value) {}
		const_expression(const const_expression & b) : value(b.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return value; }
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const { return value; }
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const { return value; }
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const {}

		basic_expression * Copy() const { return static_cast<basic_expression *>(new const_expression(*this)); }
	};

	
	class measure_expression : public shell_expression<measure_expression>
	{
	public:
		measure_expression() {}
		measure_expression(const measure_expression & b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { Storage::real measure; aut.GetMesh()->GetGeometricData(static_cast<Element *>(elem), MEASURE, &measure); return measure; }
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const { Storage::real measure; aut.GetMesh()->GetGeometricData(static_cast<Element *>(elem), MEASURE, &measure); return measure; }
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const { Storage::real measure; aut.GetMesh()->GetGeometricData(static_cast<Element *>(elem), MEASURE, &measure); return measure; }
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const {}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new measure_expression(*this)); }
	};

	
	class function_expression : public shell_expression<function_expression>
	{
		Automatizator::func_callback func;
		INMOST_DATA_ENUM_TYPE funcid;
	public:
		function_expression(Automatizator & aut, INMOST_DATA_ENUM_TYPE funcid) : func(Automatizator::GetInnerFunction(aut,funcid)), funcid(funcid) {}
		function_expression(const function_expression & b) :func(b.func) , funcid(b.funcid) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return func(elem, user_data); }
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const { return func(elem,user_data); }
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const { return func(elem, user_data); }
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const {}

		basic_expression * Copy() const { return static_cast<basic_expression *>(new function_expression(*this)); }
	};
	
	class tag_expression : public shell_expression<tag_expression>
	{
		INMOST_DATA_ENUM_TYPE tagop;
		INMOST_DATA_ENUM_TYPE comp;
		Tag value;
		Tag index;
		bool is_dynamic;
	public:
		tag_expression(Automatizator & aut, INMOST_DATA_ENUM_TYPE tagop, INMOST_DATA_ENUM_TYPE comp) : tagop(tagop), comp(comp) 
		{
			is_dynamic = tagop < AD_CTAG;
			if (is_dynamic)
			{
				value = aut.GetDynamicValueTag(tagop);
				index = aut.GetDynamicIndexTag(tagop);
			}
			else value = aut.GetStaticValueTag(tagop);
		}
		tag_expression(const tag_expression & b) : tagop(b.tagop), comp(b.comp), value(b.value), index(b.index), is_dynamic(b.is_dynamic) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return elem->RealArray(value)[comp]; }
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{ 
			if( is_dynamic ) out[elem->IntegerArray(index)[comp]] += 1.0; 
			return elem->RealArray(value)[comp];
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{ 
			return elem->RealArray(value)[comp];
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{ 
			if (is_dynamic) entries[elem->IntegerArray(index)[comp]] += multval;
		}
		virtual basic_expression * Copy() const { return static_cast<basic_expression *>(new tag_expression(*this)); }
	};
	
	template <class A>
	class unary_expression : virtual public basic_expression //: public shell_expression<unary_expression<A> >
	{
	protected:
		A * Arg;
	public:
		unary_expression(const A & Arg) : Arg(new A(Arg)) { }
		unary_expression(const unary_expression & B) : Arg(new A(*B.Arg)) {}
		virtual ~unary_expression() {}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const = 0;
		__INLINE virtual INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const = 0;
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const = 0;
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const = 0;
		virtual basic_expression * Copy() const = 0;
	};
	
	template <class A, class B>
	class binary_expression : virtual public basic_expression //: public shell_expression<binary_expression<A, B> >
	{
	protected:
		A * Left;
		B * Right;
	public:
		binary_expression(const A & Left, const B & Right) : Left(new A(Left)), Right(new B(Right)) { }
		binary_expression(const binary_expression & Q) : Left(new A(*Q.Left)), Right(new B(*Q.Right)) {}
		virtual ~binary_expression() {}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const = 0;
		__INLINE virtual INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const = 0;
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const = 0;
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const = 0;
		virtual basic_expression * Copy() const = 0;
	};
	
	template <class A, class B, class C>
	class ternary_expression : virtual public basic_expression //: public shell_expression< ternary_expression<A, B, C> >
	{
	protected:
		A * A1;
		B * A2;
		C * A3;
	public:
		ternary_expression(const A & A1, const B & A2, const C & A3) : A1(new A(A1)), A2(new B(A2)), A3(new C(A3)) { }
		ternary_expression(const ternary_expression & Q) : A1(new A(*Q.A1)), A2(new B(*Q.A2)), A3(new C(*Q.A3)) {}
		virtual ~ternary_expression() {}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const = 0;
		__INLINE virtual INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const = 0;
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const = 0;
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const = 0;
		virtual basic_expression * Copy() const = 0;
	};
	
	template<class A>
	class unary_minus_expression : public shell_expression<unary_minus_expression<A> >, protected unary_expression<A>
	{
	public:
		unary_minus_expression(const A & Arg) : unary_expression(Arg) {}
		unary_minus_expression(const unary_minus_expression & b) : unary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return -Arg->GetValue(elem, aut, user_data); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE ret = Arg->GetDerivative(elem, aut, user_data, out);
			MultArray(out, -1.0);
			return -ret;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const { return -Arg->DerivativePrecompute(elem, aut, user_data, out); }
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			Arg->DerivativeFill(elem, aut, user_data, values, entries, multval * -1.0);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new unary_minus_expression<A>(*this)); }
	};
	
	template<class A>
	class abs_expression : public shell_expression<abs_expression<A> >, protected unary_expression<A>
	{
	public:
		abs_expression(const A & Arg) : unary_expression(Arg) {}
		abs_expression(const abs_expression & b) : unary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return ::fabs(Arg->GetValue(elem, aut, user_data)); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->GetDerivative(elem, aut, user_data, out);
			MultArray(out, (arg > 0.0 ? 1.0 : -1.0))
			return ::fabs(arg);
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const 
		{
			INMOST_DATA_REAL_TYPE arg = Arg->DerivativePrecompute(elem, aut, user_data, out); 
			out.push_back(arg);
			return ::fabs(arg);
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE arg = values.back();  values.pop_back();
			Arg->DerivativeFill(elem, aut, user_data, values, entries, multval * (arg > 0.0 ? 1.0 : -1.0));
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new abs_expression<A>(*this)); }
	};

	template<class A>
	class exp_expression : public shell_expression<exp_expression<A> >, protected unary_expression<A>
	{
	public:
		exp_expression(const A & Arg) : unary_expression(Arg) {}
		exp_expression(const exp_expression & b) : unary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return ::exp(Arg->GetValue(elem, aut, user_data)); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->GetDerivative(elem, aut, user_data, out), ret = ::exp(arg);
			MultArray(out, ret);
			return ret;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->DerivativePrecompute(elem, aut, user_data, out), ret = ::exp(arg);
			out.push_back(ret);
			return ret;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE ret = values.back();  values.pop_back();
			Arg->DerivativeFill(elem, aut, user_data, values, entries, multval * ret);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new exp_expression<A>(*this)); }
	};

	template<class A>
	class log_expression : public shell_expression<log_expression<A> >, protected unary_expression<A>
	{
	public:
		log_expression(const A & Arg) : unary_expression(Arg) {}
		log_expression(const log_expression & b) : unary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return ::log(Arg->GetValue(elem, aut, user_data)); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->GetDerivative(elem, aut, user_data, out);
			MultArray(out, 1.0 / arg);
			return log(arg);
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->DerivativePrecompute(elem, aut, user_data, out), ret = ::log(arg);
			out.push_back(arg);
			return ret;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE arg = values.back();  values.pop_back();
			Arg->DerivativeFill(elem, aut, user_data, values, entries, multval / arg);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new log_expression<A>(*this)); }
	};

	template<class A>
	class sin_expression : public shell_expression<sin_expression<A> >, protected unary_expression<A>
	{
	public:
		sin_expression(const A & Arg) : unary_expression(Arg) {}
		sin_expression(const sin_expression & b) : unary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return ::sin(Arg->GetValue(elem, aut, user_data)); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->GetDerivative(elem, aut, out);
			MultArray(out, ::cos(arg));
			return ::sin(arg);
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->DerivativePrecompute(elem, aut, user_data, out), ret = ::sin(arg);
			out.push_back(arg);
			return ret;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE arg = values.back();  values.pop_back();
			Arg->DerivativeFill(elem, aut, user_data, values, entries, multval * :: cos(arg));
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new sin_expression<A>(*this)); }
	};

	template<class A>
	class cos_expression : public shell_expression<cos_expression<A> >, protected unary_expression<A>
	{
	public:
		cos_expression(const A & Arg) : unary_expression(Arg) {}
		cos_expression(const cos_expression & b) : unary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return ::cos(Arg->GetValue(elem, aut, user_data)); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->GetDerivative(elem, aut, out);
			MultArray(out, -::sin(arg));
			return ::cos(arg);
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->DerivativePrecompute(elem, aut, user_data, out), ret = ::cos(arg);
			out.push_back(arg);
			return ret;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE arg = values.back();  values.pop_back();
			Arg->DerivativeFill(elem, aut, user_data, values, entries, multval * -::sin(arg));
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new cos_expression<A>(*this)); }
	};

	template<class A>
	class table_expression : public shell_expression<table_expression<A> >, protected unary_expression<A>
	{
		Automatizator::table_ptr table;
		INMOST_DATA_ENUM_TYPE tableind;
	public:
		table_expression(Automatizator & aut, INMOST_DATA_ENUM_TYPE tableind, const A & Arg) : unary_expression(Arg), table(Automatizator::GetInnerTable(aut,tableind)), tableind(tableind) {}
		table_expression(const table_expression & b) : unary_expression(b), table(b.table), tableind(b.tableind) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return table->get_value(Arg->GetValue(elem, aut, user_data)); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->GetDerivative(elem, aut, user_data, out);
			std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> vd = table->get_both(arg);
			MultArray(out, vd.second);
			return vd.first;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE arg = Arg->DerivativePrecompute(elem, aut, user_data, out);
			std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> vd = table->get_both(arg);
			out.push_back(vd.second);
			return vd.first;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE der = values.back();  values.pop_back();
			Arg->DerivativeFill(elem, aut, user_data, values, entries, multval * der);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new table_expression<A>(*this)); }
	};

	template<class A>
	class stencil_expression : public shell_expression<stencil_expression<A> >, protected unary_expression<A>
	{
		INMOST_DATA_ENUM_TYPE stnclind;
	public:
		stencil_expression(INMOST_DATA_ENUM_TYPE stnclind, const A & Arg) : stnclind(stnclind), unary_expression(Arg) {}
		stencil_expression(const stencil_expression & b) : unary_expression(b), stnclind(b.stnclind) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const
		{
			INMOST_DATA_REAL_TYPE ret = 0.0;
			Automatizator::stencil_pairs stncl;
			aut.GetStencil(stnclind, elem, user_data, stncl);
			for (INMOST_DATA_ENUM_TYPE k = 0; k < stncl.size(); ++k) ret += Arg->GetValue(stncl[k].first, aut, user_data)*stncl[k].second;
			return ret;
		};
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			INMOST_DATA_REAL_TYPE ret = 0.0;
			Automatizator::stencil_pairs stncl;
			aut.GetStencil(stnclind, elem, user_data, stncl);
			for (INMOST_DATA_ENUM_TYPE k = 0; k < stncl.size(); ++k)
			{
				DerivativeArray temp;
				ret += Arg->GetDerivative(stncl[k].first, aut, user_data, temp)*stncl[k].second;
				MultArray(temp, stncl[k].second);
				AddArray(out, temp);
			}
			return ret;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE ret = 0.0;
			Automatizator::stencil_pairs stncl;
			aut.GetStencil(stnclind, elem, user_data, stncl);
			for (INMOST_DATA_ENUM_TYPE k = 0; k < stncl.size(); ++k) ret += Arg->DerivativePrecompute(stncl[k].first, aut, user_data,out)*stncl[k].second;
			return ret;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE ret = 0.0;
			Automatizator::stencil_pairs stncl;
			aut.GetStencil(stnclind, elem, user_data, stncl);
			for (INMOST_DATA_ENUM_TYPE k = stncl.size(); k > 0; --k) Arg->DerivativeFill(stncl[k - 1].first, aut, user_data, values, entries, multval*stncl[k-1].second);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new stencil_expression<A>(*this)); }
	};


	template<class A, class B>
	class multiply_expression : public shell_expression<multiply_expression<A, B> >, protected binary_expression<A, B>
	{
	public:
		multiply_expression(const A & Left, const B & Right) : binary_expression(Left, Right) {}
		multiply_expression(const multiply_expression & b) : binary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return Left->GetValue(elem, aut, user_data)*Right->GetValue(elem, aut, user_data); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			DerivativeArray temp;
			INMOST_DATA_REAL_TYPE lval = Left->GetDerivative(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->GetDerivative(elem, aut, user_data, temp);
			MultArray(out, rval);
			MultArray(temp, lval);
			AddArray(out, temp);
			return lval*rval;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE lval = Left->DerivativePrecompute(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->DerivativePrecompute(elem, aut, user_data, out);
			out.push_back(lval);
			out.push_back(rval);
			return lval*rval;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE rval = values.back(); values.pop_back();
			INMOST_DATA_REAL_TYPE lval = values.back(); values.pop_back();
			Right->DerivativeFill(elem, aut, user_data, values, entries, multval*lval);
			Left->DerivativeFill(elem, aut, user_data, values, entries, multval*rval);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new multiply_expression<A,B>(*this)); }
	};

	template<class A, class B>
	class division_expression : public shell_expression<division_expression<A, B> >, protected binary_expression<A, B>
	{
	public:
		division_expression(const A & Left, const B & Right) : binary_expression(Left, Right) {}
		division_expression(const division_expression & b) : binary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return Left->GetValue(elem, aut, user_data) / Right->GetValue(elem, aut, user_data); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			DerivativeArray temp;
			INMOST_DATA_REAL_TYPE lval = Left->GetDerivative(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->GetDerivative(elem, aut, user_data, temp);
			MultArray(temp, -lval / (rval*rval));
			MultArray(out, 1.0 / rval);
			AddArray(out, temp);
			return lval / rval;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE lval = Left->DerivativePrecompute(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->DerivativePrecompute(elem, aut, user_data, out);
			out.push_back(lval);
			out.push_back(rval);
			return lval/rval;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE rval = values.back(); values.pop_back();
			INMOST_DATA_REAL_TYPE lval = values.back(); values.pop_back();
			Right->DerivativeFill(elem, aut, user_data, values, entries, -multval*lval/(rval*rval));
			Left->DerivativeFill(elem, aut, user_data, values, entries, multval/rval);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new division_expression<A, B>(*this)); }
	};

	template<class A, class B>
	class addition_expression : public shell_expression<addition_expression<A, B> >, protected binary_expression<A, B>
	{
	public:
		addition_expression(const A & Left, const B & Right) : binary_expression(Left, Right) {}
		addition_expression(const addition_expression & b) : binary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return Left->GetValue(elem, aut, user_data) + Right->GetValue(elem, aut, user_data); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			DerivativeArray temp;
			INMOST_DATA_REAL_TYPE lval = Left->GetDerivative(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->GetDerivative(elem, aut, user_data, temp);
			AddArray(out, temp);
			return lval + rval;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE lval = Left->DerivativePrecompute(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->DerivativePrecompute(elem, aut, user_data, out);
			return lval + rval;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			Right->DerivativeFill(elem, aut, user_data, values, entries, multval);
			Left->DerivativeFill(elem, aut, user_data, values, entries, multval);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new addition_expression<A, B>(*this)); }
	};

	template<class A, class B>
	class substraction_expression : public shell_expression<substraction_expression<A, B> >, protected binary_expression<A, B>
	{
	public:
		substraction_expression(const A & Left, const B & Right) : binary_expression(Left, Right) {}
		substraction_expression(const substraction_expression & b) : binary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return Left->GetValue(elem, aut, user_data) - Right->GetValue(elem, aut, user_data); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			DerivativeArray temp;
			INMOST_DATA_REAL_TYPE lval = Left->GetDerivative(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->GetDerivative(elem, aut, user_data, temp);
			MultArray(temp, -1.0);
			AddArray(out, temp);
			return lval - rval;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE lval = Left->DerivativePrecompute(elem, aut, user_data, out);
			INMOST_DATA_REAL_TYPE rval = Right->DerivativePrecompute(elem, aut, user_data, out);
			return lval - rval;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			Right->DerivativeFill(elem, aut, user_data, values, entries, -multval);
			Left->DerivativeFill(elem, aut, user_data, values, entries, multval);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new substraction_expression<A, B>(*this)); }
	};

	template<class A, class B>
	class pow_expression : public shell_expression<pow_expression<A, B> >, protected binary_expression<A, B>
	{
	public:
		pow_expression(const A & Left, const B & Right) : binary_expression(Left, Right) {}
		pow_expression(const pow_expression & b) : binary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return pow(Left->GetValue(elem, aut, user_data), Right->GetValue(elem, aut, user_data)); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			DerivativeArray temp;
			INMOST_DATA_REAL_TYPE lval = Left->GetDerivative(elem, aut, user_data, out), rval = Right->GetDerivative(elem, aut, user_data, temp), ret = ::pow(lval,rval);
			MultArray(temp, ret * ::log(lval));
			MultArray(out,ret * rval / lval);
			AddArray(out, temp);
			return ret;
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE lval = Left->DerivativePrecompute(elem, aut, user_data, out), rval = Right->DerivativePrecompute(elem, aut, user_data, out), ret = ::pow(lval,rval);
			out.push_back(lval);
			out.push_back(rval);
			out.push_back(ret);
			return ret;
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE ret = values.back(); values.pop_back();
			INMOST_DATA_REAL_TYPE rval = values.back(); values.pop_back();
			INMOST_DATA_REAL_TYPE lval = values.back(); values.pop_back();
			Right->DerivativeFill(elem, aut, user_data, values, entries, multval*ret* ::log(lval));
			Left->DerivativeFill(elem, aut, user_data, values, entries, multval*ret * rval / lval);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new pow_expression<A, B>(*this)); }
	};
	
	template<class A, class B, class C>
	class condition_expression : public shell_expression<condition_expression<A, B, C> >, protected ternary_expression<A, B, C>
	{
	public:
		condition_expression(const A & A1, const B & A2, const C & A3) : ternary_expression(A1, A2, A3) {}
		condition_expression(const condition_expression & b) : ternary_expression(b) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue(Storage * elem, Automatizator & aut, void * user_data) const { return A1->GetValue(elem, aut, user_data) > 0.0 ? A2->GetValue(elem, aut, user_data) : A3->GetValue(elem, aut, user_data); };
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(Storage * elem, Automatizator & aut, void * user_data, DerivativeArray & out) const
		{
			return A1->GetValue(elem, aut, user_data) > 0.0 ? A2->GetDerivative(elem, aut,user_data,out) : A3->GetDerivative(elem, aut, user_data,out);
		}
		__INLINE INMOST_DATA_REAL_TYPE DerivativePrecompute(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & out) const
		{
			INMOST_DATA_REAL_TYPE cond = A1->GetValue(elem, aut, user_data);
			out.push_back(cond);
			return cond > 0.0 ? A2->DerivativePrecompute(elem, aut, user_data, out) : A3->DerivativePrecompute(elem, aut, user_data, out);
		}
		__INLINE void                DerivativeFill(Storage * elem, Automatizator & aut, void * user_data, precomp_values_t & values, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval) const
		{
			INMOST_DATA_REAL_TYPE cond = values.back(); values.pop_back();
			cond > 0.0 ? A2->DerivativeFill(elem, aut, user_data, values, entries, multval) : A3->DerivativeFill(elem, aut, user_data, values, entries, multval);
		}
		basic_expression * Copy() const { return static_cast<basic_expression *>(new condition_expression<A, B, C>(*this)); }
	};

	
										__INLINE         const_expression              cval(INMOST_DATA_REAL_TYPE val) { return const_expression(val); }
	                                    __INLINE       measure_expression         measuret() { return measure_expression(); }
										__INLINE           tag_expression          tagvalt(Automatizator & aut, INMOST_DATA_ENUM_TYPE reg_tag, INMOST_DATA_ENUM_TYPE comp = 0) { assert(reg_tag >= AD_TAG && reg_tag < AD_STNCL); return tag_expression(aut, reg_tag, comp); }
										__INLINE      function_expression         funcvalt(Automatizator & aut, INMOST_DATA_ENUM_TYPE funcid) { return function_expression(aut, funcid); }

	//unary operators
	template<class A>                   __INLINE   unary_minus_expression<A>      operator-(shell_expression<A> const & Arg) { return unary_minus_expression<A>(static_cast<A const &>(Arg)); }
	template<class A>                   __INLINE           abs_expression<A>            abs(shell_expression<A> const & Arg) { return abs_expression<A>(static_cast<A const &>(Arg)); }
	template<class A>                   __INLINE           exp_expression<A>            exp(shell_expression<A> const & Arg) { return exp_expression<A> (static_cast<A const &>(Arg)); }
	template<class A>                   __INLINE           log_expression<A>            log(shell_expression<A> const & Arg) { return log_expression<A> (static_cast<A const &>(Arg)); }
	template<class A>                   __INLINE           sin_expression<A>            sin(shell_expression<A> const & Arg) { return sin_expression<A> (static_cast<A const &>(Arg)); }
	template<class A>                   __INLINE           cos_expression<A>            cos(shell_expression<A> const & Arg) { return cos_expression<A> (static_cast<A const &>(Arg)); }
	template<class A>                   __INLINE       stencil_expression<A>        stencil(INMOST_DATA_ENUM_TYPE stnclid, shell_expression<A> const & Arg) { return stencil_expression<A> (stnclid, static_cast<A const &>(Arg)); }
	template<class A>                   __INLINE         table_expression<A>         tabval(Automatizator & aut, INMOST_DATA_ENUM_TYPE tableid, shell_expression<A> const & Arg) { return table_expression<A> (aut, tableid, static_cast<A const &>(Arg)); }
	//binary operators
	template<class A, class B>          __INLINE      multiply_expression<A, B>   operator*(shell_expression<A> const & Left, shell_expression<B> const & Right) { return multiply_expression<A, B> (static_cast<A const &>(Left), static_cast<B const &>(Right)); }
	template<class A, class B>          __INLINE      division_expression<A, B>   operator/(shell_expression<A> const & Left, shell_expression<B> const & Right) { return division_expression<A, B> (static_cast<A const &>(Left), static_cast<B const &>(Right)); }
	template<class A, class B>          __INLINE      addition_expression<A, B>   operator+(shell_expression<A> const & Left, shell_expression<B> const & Right) { return addition_expression<A, B> (static_cast<A const &>(Left), static_cast<B const &>(Right)); }
	template<class A, class B>          __INLINE  substraction_expression<A, B>   operator-(shell_expression<A> const & Left, shell_expression<B> const & Right) { return substraction_expression<A, B> (static_cast<A const &>(Left), static_cast<B const &>(Right)); }
	template<class A, class B>          __INLINE           pow_expression<A, B>         pow(shell_expression<A> const & Left, shell_expression<B> const & Right) { return pow_expression<A, B> (static_cast<A const &>(Left), static_cast<B const &>(Right)); }

	template<class A, class B, class C> __INLINE    condition_expression<A, B, C> condition(shell_expression<A> const & A1, shell_expression<B> const & A2, shell_expression<C> const & A3) { return condition_expression<A, B, C> (static_cast<A const &>(A1), static_cast<B const &>(A2), static_cast<C const &>(A3)); }

	/*
	template<class B>          __INLINE      multiply_expression<const_expression, B>    operator*(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return multiply_expression<const_expression, B>(const_expression(Left), Right); }
	template<class A>          __INLINE      multiply_expression<A, const_expression>    operator*(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return multiply_expression<A, const_expression>(Left, const_expression(Right)); }
	template<class B>          __INLINE      division_expression<const_expression, B>    operator/(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return division_expression<const_expression, B>(const_expression(Left), Right); }
	template<class A>          __INLINE      division_expression<A, const_expression>    operator/(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return division_expression<A, const_expression>(Left, const_expression(Right)); }
	template<class B>          __INLINE      addition_expression<const_expression, B>    operator+(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return addition_expression<const_expression, B>(const_expression(Left), Right); }
	template<class A>          __INLINE      addition_expression<A, const_expression>    operator+(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return addition_expression<A, const_expression>(Left, const_expression(Right)); }
	template<class B>          __INLINE  substraction_expression<const_expression, B>    operator-(INMOST_DATA_REAL_TYPE Left, shell_expression<B> const & Right) { return substraction_expression<const_expression, B>(const_expression(Left), Right); }
	template<class A>          __INLINE  substraction_expression<A, const_expression>    operator-(shell_expression<A> const & Left, INMOST_DATA_REAL_TYPE Right) { return substraction_expression<A, const_expression>(Left, const_expression(Right)); }
	*/

	
#endif


	class expr
	{
		INMOST_DATA_ENUM_TYPE op;
		INMOST_DATA_REAL_TYPE coef;
		expr * left;
		expr * right;
	public:
		expr() : op(AD_NONE), coef(1), left(NULL), right(NULL) { }
		expr(const expr & other) : op(other.op), coef(other.coef)
		{
			if (other.op >= AD_TAG && other.op < AD_STNCL) left = other.left; //copy component number
			else if (other.left != NULL) left = new expr(*other.left); else left = NULL;
			if (other.right != NULL) right = new expr(*other.right); else right = NULL;
		}
		expr(INMOST_DATA_REAL_TYPE val) : op(AD_CONST), coef(val), left(NULL), right(NULL) {}
		expr(INMOST_DATA_ENUM_TYPE new_op, expr * l, expr * r) : op(new_op), coef(1), left(l), right(r) {}
		expr(INMOST_DATA_ENUM_TYPE tag_op, INMOST_DATA_ENUM_TYPE comp = 0) :op(tag_op), coef(1), left(reinterpret_cast<expr *>(comp)), right(NULL) {}
		~expr()
		{
			if (op < AD_COS)
			{
				delete left;
				delete right;
			}
			else if (op < AD_CONST)
				delete left;
			else if (op >= AD_STNCL && op < AD_TABLE)
				delete left;
		}
		expr & operator =(expr const & other)
		{
			op = other.op; coef = other.coef;
			if (other.op >= AD_TAG && other.op < AD_STNCL)
			{
				left = other.left; //copy component number
				right = other.right;
			}
			else if (other.left != NULL) left = new expr(*other.left); else left = NULL;
			if (other.right != NULL) right = new expr(*other.right); else right = NULL;
			return *this;
		}
		expr operator +() { return expr(*this); }
		expr operator -() { expr var(*this); var.coef *= -1.0; return var; }
		expr operator +(const expr & other) const { return expr(AD_PLUS, new expr(*this), new expr(other)); }
		expr operator -(const expr & other) const { return expr(AD_MINUS, new expr(*this), new expr(other)); }
		expr operator *(const expr & other) const { return expr(AD_MULT, new expr(*this), new expr(other)); }
		expr operator /(const expr & other) const { return expr(AD_DIV, new expr(*this), new expr(other)); }
		expr operator +(const INMOST_DATA_REAL_TYPE & other) const { return expr(AD_PLUS, new expr(*this), new expr(other)); }
		expr operator -(const INMOST_DATA_REAL_TYPE & other) const { return expr(AD_MINUS, new expr(*this), new expr(other)); }
		expr operator *(const INMOST_DATA_REAL_TYPE & other) const { expr var(*this); var.coef *= other; return var; }
		expr operator /(const INMOST_DATA_REAL_TYPE & other) const { expr var(*this); var.coef /= other; return var; }
		__INLINE friend expr operator +(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		__INLINE friend expr operator -(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		__INLINE friend expr operator *(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		__INLINE friend expr operator /(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		bool operator ==(const expr & other) {return this == &other || (op == other.op && left == other.left && right == other.right);}
		bool is_endpoint() { return op >= AD_TAG && op < AD_STNCL; }
		friend class Automatizator;
	};
	

	__INLINE expr operator +(const INMOST_DATA_REAL_TYPE & left, const expr & right) { return expr(AD_PLUS, new expr(left), new expr(right)); }
	__INLINE expr operator -(const INMOST_DATA_REAL_TYPE & left, const expr & right) { return expr(AD_MINUS, new expr(left), new expr(right)); }
	__INLINE expr operator *(const INMOST_DATA_REAL_TYPE & left, const expr & right) { expr var(right); var.coef *= left; return var; }
	__INLINE expr operator /(const INMOST_DATA_REAL_TYPE & left, const expr & right) { expr var; var.op = AD_INV; var.coef = left; var.left = new expr(right); var.right = NULL; return var; }
	__INLINE expr pow(const expr & v, const expr n) { return expr(AD_POW, new expr(v), new expr(n)); }
	__INLINE expr abs(const expr & v) { return expr(AD_ABS, new expr(v), NULL); }
	__INLINE expr exp(const expr & v) { return expr(AD_EXP, new expr(v), NULL); }
	__INLINE expr log(const expr & v) { return expr(AD_LOG, new expr(v), NULL); }
	__INLINE expr sin(const expr & v) { return expr(AD_SIN, new expr(v), NULL); }
	__INLINE expr cos(const expr & v) { return expr(AD_COS, new expr(v), NULL); }
	__INLINE expr measure() { return expr(AD_MES, NULL, NULL); }
	__INLINE expr condition(const expr & cond, const expr & if_true, const expr & if_false) { return expr(AD_COND, new expr(cond), new expr(AD_ALTR, new expr(if_true), new expr(if_false))); }
	__INLINE expr stencil(INMOST_DATA_ENUM_TYPE stncl, const expr & v) { assert(stncl >= AD_STNCL && stncl < AD_TABLE); return expr(stncl, new expr(v), NULL); }
	__INLINE expr tabval(INMOST_DATA_ENUM_TYPE tabl, const expr & v) { assert(tabl >= AD_TABLE && tabl < AD_FUNC); return expr(tabl, new expr(v), NULL); }
	__INLINE expr tagval(INMOST_DATA_ENUM_TYPE reg_tag, INMOST_DATA_ENUM_TYPE comp = 0) { assert(reg_tag >= AD_TAG && reg_tag < AD_STNCL); return expr(reg_tag, comp); }
	__INLINE expr funcval(INMOST_DATA_ENUM_TYPE reg_func) { assert(reg_func >= AD_FUNC); return expr(reg_func, NULL, NULL); }
}

#endif
#endif
