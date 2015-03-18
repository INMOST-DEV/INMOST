
#ifndef INMOST_AUTODIFF_H_INCLUDED
#define INMOST_AUTODIFF_H_INCLUDED
#include "inmost_common.h"
#include "inmost_mesh.h"
#include <sstream> //for debug

//#define NEW_VERSION

#if defined(USE_AUTODIFF) && (!defined(USE_MESH))
#warning "USE_AUTODIFF require USE_MESH"
#undef USE_AUTODIFF
#endif

//#define DPRNT



#if defined(USE_AUTODIFF)
#include <math.h>


#define FILTER_EPS 1e-12

#define AD_SPACE 16384 //typical number of registered entries

#define AD_NONE  0 //this should generate an error
//binary operations below
#define AD_PLUS  1 //sum x+y
#define AD_MINUS 2 //substraction x-y
#define AD_MULT  3 //multiplication x*y
#define AD_DIV   4 //devision x/y
#define AD_INV   5 //inversion 1/x, depricated
#define AD_POW   6 //power operation x^y
#define AD_SQRT  7 //square root operation


#define AD_PRECOMP 15 // this expression points to precomputed replacement, expr * left is position in array of precomputed values, expr * right points to original expression

#define AD_VAL   19 // calculate only value of current expression

//unary operations below
#define AD_COS   20 //cos(x)

#define AD_ABS   22 // |x|
#define AD_EXP   23 // exp(x)
#define AD_LOG   24 // log(x)
#define AD_SIN   25 // sin(x)

#define AD_CONST 50 // expr * left represents const value
#define AD_MES   51 // volume(e) for cell, area(e) for face, length(e) for edge

#define AD_COND  60 //condition
#define AD_COND_TYPE 61 //condition on element type
#define AD_COND_MARK 62 //condition on marker
#define AD_ALTR  65 //alternatives for condition



#define AD_EXT   99 // pointer to external expression data

#define AD_TAG   100               //tag with dynamic data
#define AD_CTAG  (AD_TAG+AD_SPACE)   //tag with constant data
#define AD_STNCL (AD_CTAG+AD_SPACE)  //stencil that sets current elements and their coefs
#define AD_TABLE (AD_STNCL+AD_SPACE) //table of values
#define AD_FUNC  (AD_TABLE+AD_SPACE) //register function that calculates some value from element


//types of operands

#define AD_TYPE_INVALID   ENUMUNDEF
#define AD_TYPE_ENDPOINT  0
#define AD_TYPE_UNARY     1
#define AD_TYPE_BINARY    2
#define AD_TYPE_TERNARY   3
#define AD_TYPE_VALUE     4

namespace INMOST
{
	class Automatizator; //forward declaration
	typedef dynarray<INMOST_DATA_REAL_TYPE, 2048> precomp_values_t;
	class expr;
	
	

	class Automatizator
	{
	public:
		typedef INMOST_DATA_REAL_TYPE(*func_callback)(const Storage & current_element, void * user_data);
		typedef std::pair<HandleType, INMOST_DATA_REAL_TYPE> stencil_pair;
		typedef dynarray<stencil_pair, 64> stencil_pairs;
		typedef void(*stencil_callback)(const Storage & current_element, stencil_pairs & out_stencil, void * user_data);
	private:
		typedef struct{ Tag t; MarkerType domain_mask; } tagdomain;
		typedef dynarray<tagdomain, 128> const_tag_type;
		typedef struct{ tagdomain d; Tag indices; } tagpair;
		typedef dynarray<tagpair, 128> tagpairs_type;
		typedef std::vector<tagpair> index_enum;
		typedef struct { Tag elements, coefs; } stencil_tag;
		typedef struct { std::string name; INMOST_DATA_ENUM_TYPE kind; MarkerType domainmask; void * link; } stencil_kind_domain;
		typedef dynarray<stencil_kind_domain, 128> stencil_type;
		typedef struct func_name_callback_t { std::string name; func_callback func; } func_name_callback;
		typedef dynarray<func_name_callback, 128> func_type;
	public:
		typedef struct
		{
			std::string name;
			INMOST_DATA_REAL_TYPE * args;
			INMOST_DATA_REAL_TYPE * vals;
			INMOST_DATA_ENUM_TYPE size;
			INMOST_DATA_ENUM_TYPE binary_search(INMOST_DATA_REAL_TYPE arg);
			INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE arg);
			INMOST_DATA_REAL_TYPE get_derivative(INMOST_DATA_REAL_TYPE arg);
			std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> get_both(INMOST_DATA_REAL_TYPE arg);
		} table;
		typedef table * table_ptr;
	private:
		typedef dynarray<INMOST_DATA_REAL_TYPE, 128> values_container;
		typedef dynarray<table_ptr, 32>              table_type;
		const_tag_type        reg_ctags;
		index_enum            index_tags;
		tagpairs_type         reg_tags;
		table_type            reg_tables;
		stencil_type          reg_stencils;
		func_type             reg_funcs;
		INMOST_DATA_ENUM_TYPE first_num;
		INMOST_DATA_ENUM_TYPE last_num;
		Mesh * m;
#if defined(NEW_VERSION)
		INMOST_DATA_REAL_TYPE                                            EvaluateSub(expr & var,INMOST_DATA_ENUM_TYPE element,INMOST_DATA_ENUM_TYPE parent, void * user_data);
		void                                                             DerivativeFill(expr & var, INMOST_DATA_ENUM_TYPE element, INMOST_DATA_ENUM_TYPE parent, Solver::Row & entries, INMOST_DATA_REAL_TYPE multval, void * user_data);
#else
		INMOST_DATA_REAL_TYPE                                            DerivativePrecompute(const expr & var, const Storage & e, precomp_values_t & values, void * user_data);
		void                                                             DerivativeFill(const expr & var, const Storage & e, Solver::Row & entries, precomp_values_t & values, INMOST_DATA_REAL_TYPE multval, void * user_data);
#endif
	public:
		Automatizator(Mesh * m);
		~Automatizator();
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetFirstIndex() { return first_num; }
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetLastIndex() { return last_num; }
		INMOST_DATA_ENUM_TYPE                                            RegisterFunc(std::string name, func_callback func);
		INMOST_DATA_ENUM_TYPE                                            RegisterStencil(std::string name, Tag elements_tag, Tag coefs_tag, MarkerType domain_mask = 0);
		INMOST_DATA_ENUM_TYPE                                            RegisterStencil(std::string name, stencil_callback func, MarkerType domain_mask = 0);
		INMOST_DATA_ENUM_TYPE                                            RegisterTable(std::string name, INMOST_DATA_REAL_TYPE * Arguments, INMOST_DATA_REAL_TYPE * Values, INMOST_DATA_ENUM_TYPE size);
		INMOST_DATA_ENUM_TYPE                                            RegisterDynamicTag(Tag t, ElementType typemask, MarkerType domain_mask = 0);
		INMOST_DATA_ENUM_TYPE                                            RegisterStaticTag(Tag t, MarkerType domain_mask = 0);
		void                                                             EnumerateDynamicTags();
		__INLINE Tag                                                     GetDynamicValueTag(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind-AD_TAG].d.t; }
		__INLINE Tag                                                     GetDynamicIndexTag(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind-AD_TAG].indices; }
		__INLINE MarkerType                                              GetDynamicMask(INMOST_DATA_ENUM_TYPE ind) { return reg_tags[ind-AD_TAG].d.domain_mask; }
		__INLINE Tag                                                     GetStaticValueTag(INMOST_DATA_ENUM_TYPE ind) { return reg_ctags[ind-AD_CTAG].t; }
		__INLINE MarkerType                                              GetStaticMask(INMOST_DATA_ENUM_TYPE ind) { return reg_ctags[ind-AD_CTAG].domain_mask; }
		__INLINE INMOST_DATA_REAL_TYPE                                   GetDynamicValue(const Storage & e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->RealArray(GetDynamicValueTag(ind))[comp]; }
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetDynamicIndex(const Storage & e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->IntegerArray(GetDynamicIndexTag(ind))[comp]; }
		__INLINE bool                                                    isDynamicValid(const Storage & e, INMOST_DATA_ENUM_TYPE ind) { MarkerType mask = GetDynamicMask(ind); return mask == 0 || e->GetMarker(mask); }
		__INLINE INMOST_DATA_REAL_TYPE                                   GetStaticValue(const Storage & e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->RealArray(GetStaticValueTag(ind))[comp]; }
		__INLINE bool                                                    isStaticValid(const Storage & e, INMOST_DATA_ENUM_TYPE ind) { MarkerType mask = GetStaticMask(ind); return mask == 0 || e->GetMarker(mask); }
#if defined(NEW_VERSION)
		INMOST_DATA_REAL_TYPE                                            Evaluate(expr & var, const Storage & e, void * user_data);
		INMOST_DATA_REAL_TYPE                                            Derivative(expr & var, const Storage & e, Solver::Row & out, Storage::real multiply, void * user_data);
#else
		INMOST_DATA_REAL_TYPE                                            Evaluate(const expr & var, const Storage & e, void * user_data);
		INMOST_DATA_REAL_TYPE                                            Derivative(const expr & var, const Storage & e, Solver::Row & out, Storage::real multiply,  void * user_data);
#endif
		__INLINE INMOST_DATA_REAL_TYPE                                   GetIndex(const Storage & e, INMOST_DATA_ENUM_TYPE tagind, INMOST_DATA_ENUM_TYPE comp = 0) { return e->IntegerArray(GetDynamicIndexTag(tagind))[comp]; }
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetComponents(const Storage & e, INMOST_DATA_ENUM_TYPE tagind) { return static_cast<INMOST_DATA_ENUM_TYPE>(e->IntegerArray(GetDynamicIndexTag(tagind)).size()); }
		__INLINE Mesh *                                                  GetMesh() { return m; }
		__INLINE INMOST_DATA_REAL_TYPE *                                 GetTableArguments(INMOST_DATA_ENUM_TYPE tableind) {return reg_tables[tableind-AD_TABLE]->args;}
		__INLINE INMOST_DATA_REAL_TYPE *                                 GetTableValues(INMOST_DATA_ENUM_TYPE tableind) {return reg_tables[tableind-AD_TABLE]->vals;}
		__INLINE INMOST_DATA_ENUM_TYPE                                   GetTableSize(INMOST_DATA_ENUM_TYPE tableind) {return reg_tables[tableind-AD_TABLE]->size;}
		__INLINE INMOST_DATA_REAL_TYPE                                   GetTableValue(INMOST_DATA_ENUM_TYPE tableind, INMOST_DATA_REAL_TYPE arg) { return reg_tables[tableind-AD_TABLE]->get_value(arg); }
		__INLINE INMOST_DATA_REAL_TYPE                                   GetTableDerivative(INMOST_DATA_ENUM_TYPE tableind, INMOST_DATA_REAL_TYPE arg) { return reg_tables[tableind-AD_TABLE]->get_derivative(arg); }
		__INLINE std::pair<INMOST_DATA_REAL_TYPE, INMOST_DATA_REAL_TYPE> GetTableBoth(INMOST_DATA_ENUM_TYPE tableind, INMOST_DATA_REAL_TYPE arg) { return reg_tables[tableind-AD_TABLE]->get_both(arg); }
		INMOST_DATA_ENUM_TYPE                                            GetStencil(INMOST_DATA_ENUM_TYPE stnclind, const Storage & elem, void * user_data, stencil_pairs & ret);
		__INLINE INMOST_DATA_REAL_TYPE                                   GetFunction(INMOST_DATA_ENUM_TYPE funcid, const Storage & elem, void * user_data) { return reg_funcs[funcid-AD_FUNC].func(elem, user_data); }
	};

#if defined(NEW_VERSION)
	class expr
	{
		
		typedef small_hash<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE, 1024> replace_type;
		
		class expr_data
		{
			typedef union { INMOST_DATA_REAL_TYPE r; INMOST_DATA_ENUM_TYPE i; expr * e; expr_data * q; } operand;
			INMOST_DATA_ENUM_TYPE op, priority;
			operand left, right;
		public:
			expr_data() : op(AD_NONE) { memset(&left, 0, sizeof(operand)); memset(&right, 0, sizeof(operand)); }
			expr_data(const expr_data & other) : op(other.op) 
			{
				if (other.op == AD_STNCL)
					left.e = new expr(*other.left.e);
				else if (other.op == AD_ALTR)
				{
					left.e = new expr(*other.left.e);
					right.e = new expr(*other.right.e);
				}
				else if (other.op == AD_COND)
				{
					left = other.left;
					right.q = new expr_data(*other.right.q);
				}
				else
				{
					left = other.left; 
					right = other.right;
				}
			}
			expr_data & operator =(expr_data const & other) 
			{ 
				op = other.op; 
				if (other.op == AD_STNCL)
					left.e = new expr(*other.left.e);
				else if (other.op == AD_ALTR)
				{
					left.e = new expr(*other.left.e);
					right.e = new expr(*other.right.e);
				}
				else if (other.op == AD_COND)
				{
					left = other.left;
					right.q = new expr_data(*other.right.q);
				}
				else
				{
					left = other.left;
					right = other.right;
				}
				return *this;
			}
			expr_data(INMOST_DATA_REAL_TYPE val) : op(AD_CONST) { left.r = val; }
			expr_data(INMOST_DATA_ENUM_TYPE new_op, INMOST_DATA_ENUM_TYPE l, INMOST_DATA_ENUM_TYPE r) : op(new_op) { left.i = l; right.i = r; }
			expr_data(INMOST_DATA_ENUM_TYPE new_op, INMOST_DATA_ENUM_TYPE l, INMOST_DATA_REAL_TYPE r) : op(new_op) { left.i = l; right.r = r; }
			expr_data(INMOST_DATA_ENUM_TYPE new_op, INMOST_DATA_ENUM_TYPE l, const expr & e) : op(new_op) { left.i = l; right.e = new expr(e); }
			expr_data(INMOST_DATA_ENUM_TYPE l, const expr_data & e) : op(AD_COND) { left.i = l; right.q = new expr_data(e); }
			expr_data(INMOST_DATA_ENUM_TYPE op, const expr * e) : op(op) { left.e = new expr(*e); }
			expr_data(const expr * a, const expr * b) : op(AD_ALTR) { left.e = new expr(*a); right.e = new expr(*b); }
			~expr_data()
			{
				if (op == AD_STNCL)
				{
					delete left.e;
					op = AD_NONE;
				}
				else if (op == AD_COND)
				{
					delete right.q;
					op = AD_NONE;
				}
				else if (op == AD_ALTR)
				{
					delete left.e;
					delete right.e;
					op = AD_NONE;
				}
			}
			bool operator ==(const expr_data & other) const
			{
				assert(op != AD_NONE);
				if (op == AD_EXT) return left.e->data[right.i] == other;
				if (other.op == AD_EXT) return *this == other.left.e->data[other.right.i];
				if (op != other.op) return false;
				switch (op)
				{
				case AD_MES: return true;
				case AD_CONST: return fabs(left.r - other.left.r) < 1e-9;
				case AD_PLUS:
				case AD_MINUS:
				case AD_MULT:
				case AD_DIV:
				case AD_POW:
					return left.i == other.left.i && right.i == other.right.i;
				case AD_COS:
				case AD_ABS:
				case AD_EXP:
				case AD_LOG:
				case AD_SIN:
				case AD_SQRT:
					return left.i == other.left.i;
				case AD_VAL:
					return left.i == other.right.i && (right.e == other.right.e || *right.e == *other.right.e);
				case AD_ALTR:
					if (left.e == other.left.e && right.e == other.right.e) return true;
					return *left.e == *other.left.e && *right.e == *other.right.e;
				case AD_COND:
					return left.i == other.left.i && (right.q == other.right.q || *right.q == *other.right.q);
				case AD_COND_TYPE:
				case AD_COND_MARK:
					return left.i == other.left.i;
				}
				if (op >= AD_FUNC) return true;
				if (op >= AD_TABLE) return left.i == other.left.i;
				if (op >= AD_STNCL) return (left.e == other.left.e) || (*left.e == *other.left.e);
				if (op >= AD_CTAG) return left.i == other.left.i;
				if (op >= AD_TAG) return left.i == other.left.i;
				assert(false); //cannot reach here
				return false;
			}
			bool deep_cmp(const expr_data & other, const expr & my_set, const expr & other_set ) const
			{
				assert(op != AD_NONE);
				
				if (op == AD_EXT) return left.e->data[right.i].deep_cmp(other,*left.e,other_set);
				if (other.op == AD_EXT) return deep_cmp(other.left.e->data[other.right.i],my_set, *other.left.e);
				if (op != other.op) return false;
				switch (op)
				{
				case AD_MES: return true;
				case AD_CONST: return fabs(left.r - other.left.r) < 1e-9;
				case AD_PLUS:
				case AD_MINUS:
				case AD_MULT:
				case AD_DIV:
				case AD_POW:
					return my_set.data[left.i].deep_cmp(other_set.data[other.left.i],my_set,other_set) 
						&& my_set.data[right.i].deep_cmp(other_set.data[other.right.i],my_set,other_set);
				case AD_COS:
				case AD_ABS:
				case AD_EXP:
				case AD_LOG:
				case AD_SIN:
				case AD_SQRT:
					return my_set.data[left.i].deep_cmp(other_set.data[other.left.i], my_set, other_set);
				case AD_VAL:
					return (right.e == other.right.e || *right.e == *other.right.e) && my_set.data[left.i].deep_cmp(other_set.data[other.left.i], my_set, other_set);
				case AD_ALTR:
					if (left.e == other.left.e && right.e == other.right.e) return true;
					return *left.e == *other.left.e && *right.e == *other.right.e;
				case AD_COND:
					return my_set.data[left.i].deep_cmp(other_set.data[other.left.i], my_set, other_set) && (right.q == other.right.q || *right.q == *other.right.q);
				case AD_COND_TYPE:
				case AD_COND_MARK:
					return left.i == other.left.i;
				}
				if (op >= AD_FUNC) return true;
				if (op >= AD_TABLE) return my_set.data[left.i].deep_cmp(other_set.data[other.left.i], my_set, other_set);
				if (op >= AD_STNCL) return (left.e == other.left.e) || (*left.e == *other.left.e);
				if (op >= AD_CTAG) return left.i == other.left.i;
				if (op >= AD_TAG) return left.i == other.left.i;
				assert(false); //cannot reach here
				return false;
			}
			bool cmp(const expr_data & other,replace_type & r) const
			{
				assert(op != AD_NONE);
				if (op == AD_EXT) return left.e->data[right.i] == other;
				if (other.op == AD_EXT) return *this == other.left.e->data[other.right.i];
				if (op != other.op) return false;
				switch (op)
				{
				case AD_MES: return true;
				case AD_CONST: return fabs(left.r - other.left.r) < 1e-9;
				case AD_PLUS:
				case AD_MINUS:
				case AD_MULT:
				case AD_DIV:
				case AD_POW:
					return left.i == r[other.left.i] && right.i == r[other.right.i];
				case AD_COS:
				case AD_ABS:
				case AD_EXP:
				case AD_LOG:
				case AD_SIN:
				case AD_SQRT:
					return left.i == r[other.left.i];
				case AD_VAL:
					return left.i == r[other.left.i] && (right.e == other.right.e || *right.e == *other.right.e);
				case AD_ALTR:
					if (left.e == other.left.e && right.e == other.right.e) return true;
					return *left.e == *other.left.e && *right.e == *other.right.e;
				case AD_COND:
					return left.i == r[other.left.i] && (right.q == other.right.q || *right.q == *other.right.q);
				case AD_COND_TYPE:
				case AD_COND_MARK:
					return left.i == other.left.i;
				}
				if (op >= AD_FUNC) return true;
				if (op >= AD_TABLE) return left.i == r[other.left.i];
				if (op >= AD_STNCL) return (left.e == other.left.e) || (*left.e == *other.left.e);
				if (op >= AD_CTAG) return left.i == other.left.i;
				if (op >= AD_TAG) return left.i == other.left.i;
				assert(false); //cannot reach here
				return false;
			}
			INMOST_DATA_ENUM_TYPE type() const
			{ 
				assert(op != AD_NONE);
				if (op == AD_EXT) return left.e->data[right.i].type();
				switch (op)
				{
				case AD_COND_MARK:
				case AD_COND_TYPE:
				case AD_CONST:
				case AD_MES: return AD_TYPE_ENDPOINT;
				case AD_PLUS:
				case AD_MINUS:
				case AD_MULT:
				case AD_DIV:
				case AD_POW:
				case AD_ALTR: return AD_TYPE_BINARY;
				case AD_COS:
				case AD_ABS:
				case AD_EXP:
				case AD_LOG:
				case AD_SQRT:
				case AD_SIN: return AD_TYPE_UNARY;
				case AD_COND: return AD_TYPE_TERNARY;
				case AD_VAL: return AD_TYPE_VALUE;
				}
				if (op >= AD_FUNC) return AD_TYPE_ENDPOINT;
				if (op >= AD_TABLE) return AD_TYPE_UNARY;
				if (op >= AD_STNCL) return AD_TYPE_UNARY;
				if (op >= AD_CTAG) return AD_TYPE_ENDPOINT;
				if (op >= AD_TAG) return AD_TYPE_ENDPOINT;
				assert(false); //cannot reach here
				return AD_TYPE_INVALID;
			}
			bool is_func() { return op >= AD_FUNC; }
			bool is_table() { return op >= AD_TABLE && op < AD_FUNC; }
			bool is_stncl() { return op >= AD_STNCL && op < AD_TABLE; }
			bool is_ctag() { return op >= AD_CTAG && op < AD_STNCL; }
			bool is_tag() { return op >= AD_TAG && op < AD_CTAG; }
			friend class expr;
			friend class Automatizator;
		};
		friend class expr_data;
		//typedef std::vector<expr_data> data_type;
		typedef std::vector<expr_data> data_type;
		data_type data;
		std::vector<INMOST_DATA_REAL_TYPE> values;
		Automatizator::stencil_pairs current_stencil;
		__INLINE INMOST_DATA_ENUM_TYPE values_offset(INMOST_DATA_ENUM_TYPE element) {return element * static_cast<INMOST_DATA_ENUM_TYPE>(data.size());}
		__INLINE INMOST_DATA_ENUM_TYPE derivatives_offset(INMOST_DATA_ENUM_TYPE element) {return element * static_cast<INMOST_DATA_ENUM_TYPE>(data.size()) + static_cast<INMOST_DATA_ENUM_TYPE>(data.size()*current_stencil.size());}
		void resize_for_stencil()
		{
			values.resize(current_stencil.size()*data.size()*2);
			memset(&values[current_stencil.size()*data.size()],0,sizeof(INMOST_DATA_REAL_TYPE)*current_stencil.size()*data.size());
		}

		void relink_data()
		{
			for (data_type::iterator it = data.begin(); it != data.end(); ++it)
				if (it->op == AD_COND )
				{
					for (data_type::iterator jt = it->right.q->left.e->data.begin(); jt != it->right.q->left.e->data.end(); ++jt)
					{
						if (jt->op == AD_EXT) jt->left.e = this;
					}
					it->right.q->left.e->relink_data();
					for (data_type::iterator jt = it->right.q->right.e->data.begin(); jt != it->right.q->right.e->data.end(); ++jt) 
					{
						if (jt->op == AD_EXT) jt->left.e = this;
					}
					it->right.q->left.e->relink_data();
				}
		}
		void link_data(data_type & other, replace_type & r)
		{
			INMOST_DATA_ENUM_TYPE i = 0, j;
			for (data_type::iterator it = other.begin(); it != other.end(); ++it)
			{
				j = 0;
				if (it->op == AD_EXT)
				{
					it->left.e = this;
					it->right.i = r[it->right.i];
				}
				else for (data_type::iterator jt = data.begin(); jt != data.end(); ++jt)
				{
					if (*it == *jt)
					{
						it->op = AD_EXT;
						it->left.e = this;
						it->right.i = j;
						break;
					}
					j++;
				}
				i++;
			}
		}
		void replace_operands(expr_data & e, replace_type & r)
		{
			if (e.op == AD_EXT) replace_operands(e.left.e->data[e.right.i],r);
			else if (e.op == AD_ALTR)
			{
				link_data(e.left.e->data,r);
				link_data(e.right.e->data,r);
			}
			else if (e.is_stncl()) {}
			else switch (e.type())
			{
			case AD_TYPE_UNARY:
				assert(r.is_present(e.left.i));
				e.left.i = r[e.left.i];
				break;
			case AD_TYPE_BINARY:
				assert(r.is_present(e.left.i));
				assert(r.is_present(e.right.i));
				e.left.i = r[e.left.i];
				e.right.i = r[e.right.i];
				break;
			case AD_TYPE_TERNARY:
				assert(r.is_present(e.left.i));
				e.left.i = r[e.left.i];
				replace_operands(*e.right.q, r);
				break;
			case AD_TYPE_VALUE:
				e.left.i = r[e.left.i];
				link_data(e.right.e->data,r);
				break;
			}
		}
		INMOST_DATA_ENUM_TYPE merge_data(const data_type & other)
		{
			replace_type from_to;
			INMOST_DATA_ENUM_TYPE i = 0, j;
			for (data_type::const_iterator it = other.begin(); it != other.end(); ++it)
			{
				j = 0;
				for (data_type::iterator jt = data.begin(); jt != data.end(); ++jt)
				{
					if (jt->cmp(*it,from_to))
					{
						from_to[i] = j;
						break;
					}
					j++;
				}
				if (!from_to.is_present(i))
				{
					from_to[i] = static_cast<INMOST_DATA_ENUM_TYPE>(data.size());
					data.push_back(*it);
					replace_operands(data.back(),from_to);
				}
				i++;
			}
			return from_to[static_cast<int>(other.size()) - 1];
		}
		
	public:
		expr() : data() { }
		expr(const expr & other) : data(other.data) { relink_data(); }
		expr(INMOST_DATA_ENUM_TYPE op, const expr * sub) :data() { data.push_back(expr_data(op,sub)); }
		expr(const expr * l, const expr * r) :data() { data.push_back(expr_data(l,r)); }
		expr(INMOST_DATA_REAL_TYPE val) : data() { data.push_back(expr_data(val)); }
		expr(INMOST_DATA_ENUM_TYPE op, INMOST_DATA_ENUM_TYPE comp) : data() { data.push_back(expr_data(op, comp, ENUMUNDEF)); }
		expr(INMOST_DATA_ENUM_TYPE op, const expr & operand) : data(operand.data) { relink_data();  data.push_back(expr_data(op, static_cast<INMOST_DATA_ENUM_TYPE>(data.size() - 1), ENUMUNDEF)); }
		expr(const expr & operand, const expr & multiplyer) : data(operand.data) { relink_data();  data.push_back(expr_data(AD_VAL, static_cast<INMOST_DATA_ENUM_TYPE>(data.size() - 1), multiplyer)); }
		expr(const expr & cond, const expr & if_true, const expr & if_false) :data(cond.data)
		{ 
			relink_data();  
			data.push_back(expr_data(static_cast<INMOST_DATA_ENUM_TYPE>(data.size() - 1), expr_data(&if_true, &if_false))); 
		}
		expr(INMOST_DATA_ENUM_TYPE op, const expr & l, const expr & r) :data(l.data) 
		{
			relink_data();
			INMOST_DATA_ENUM_TYPE lp = static_cast<INMOST_DATA_ENUM_TYPE>(data.size() - 1);
			INMOST_DATA_ENUM_TYPE rp = merge_data(r.data); 
			data.push_back(expr_data(op, lp, rp));
		}
		~expr() {}
		expr & operator =(expr const & other) {data = other.data; return *this;}
		expr operator +() const { return expr(*this); }
		expr operator -() const { expr zero(0.0); return zero - *this; }
		expr operator +(const expr & other) const {return expr(AD_PLUS,*this,other);}
		expr operator -(const expr & other) const { return expr(AD_MINUS, *this, other); }
		expr operator *(const expr & other) const { return expr(AD_MULT, *this, other); }
		expr operator /(const expr & other) const { return expr(AD_DIV, *this, other); }
		expr operator +(INMOST_DATA_REAL_TYPE other) const { return expr(AD_PLUS, *this, expr(other)); }
		expr operator -(INMOST_DATA_REAL_TYPE other) const { return expr(AD_MINUS, *this, expr(other)); }
		expr operator *(INMOST_DATA_REAL_TYPE other) const { return expr(AD_MULT, *this, expr(other)); }
		expr operator /(INMOST_DATA_REAL_TYPE other) const { return expr(AD_DIV, *this, expr(other)); }
		bool operator ==(const expr & other) //this should account for possible reorder! here x*y is not the same as y*x
		{
			if (data.size() != other.data.size()) return false;
			data_type::iterator it = data.begin();
			data_type::const_iterator jt = other.data.begin();
			while (it != data.end() && jt != other.data.end())
			{
				if (!it->deep_cmp(*jt, *this, other)) return false;
				++it; ++jt;
			}
			return true;
		}
		
		friend class Automatizator;
	};
	

	__INLINE expr operator +(INMOST_DATA_REAL_TYPE left, const expr & right) { return expr(left) + right; }
	__INLINE expr operator -(INMOST_DATA_REAL_TYPE left, const expr & right) { return expr(left) - right; }
	__INLINE expr operator *(INMOST_DATA_REAL_TYPE left, const expr & right) { return expr(left) * right; }
	__INLINE expr operator /(INMOST_DATA_REAL_TYPE left, const expr & right) { return expr(left) / right; }
	__INLINE expr ad_pow(const expr & v, const expr & n) { return expr(AD_POW, v, n); }
	//__INLINE expr ad_pow(const expr & v, INMOST_DATA_REAL_TYPE n) { return expr(AD_POW, v, expr(n)); }
	//__INLINE expr ad_pow(const expr & v, INMOST_DATA_ENUM_TYPE n) { return expr(AD_POW, v, expr(static_cast<INMOST_DATA_REAL_TYPE>(n))); }
	__INLINE expr ad_abs(const expr & v) { return expr(AD_ABS, v); }
	__INLINE expr ad_exp(const expr & v) { return expr(AD_EXP, v); }
	__INLINE expr ad_log(const expr & v) { return expr(AD_LOG, v); }
	__INLINE expr ad_sin(const expr & v) { return expr(AD_SIN, v); }
	__INLINE expr ad_cos(const expr & v) { return expr(AD_COS, v); }
	__INLINE expr ad_sqrt(const expr & v) {return expr(AD_SQRT, v);}
	__INLINE expr ad_val(const expr & v, const expr & multiplyer = expr(0.0)) {return expr(v, multiplyer);}
	__INLINE expr measure() { return expr(AD_MES,ENUMUNDEF); }
	__INLINE expr condition(const expr & cond, const expr & if_gt_zero, const expr & if_le_zero) { return expr(cond, if_gt_zero, if_le_zero); }
	__INLINE expr condition_etype(ElementType etypes, const expr & if_true, const expr & if_false) { return expr(expr(AD_COND_TYPE, etypes), if_true, if_false); }
	__INLINE expr condition_marker(MarkerType marker, const expr & if_true, const expr & if_false) { return expr(expr(AD_COND_MARK, marker), if_true, if_false); }
	__INLINE expr stencil(INMOST_DATA_ENUM_TYPE stncl, const expr & v) { assert(stncl >= AD_STNCL && stncl < AD_TABLE); return expr(stncl, &v); }
	__INLINE expr tabval(INMOST_DATA_ENUM_TYPE tabl, const expr & v) { assert(tabl >= AD_TABLE && tabl < AD_FUNC); return expr(tabl, v); }
	__INLINE expr tagval(INMOST_DATA_ENUM_TYPE reg_tag, INMOST_DATA_ENUM_TYPE comp = 0) { assert(reg_tag >= AD_TAG && reg_tag < AD_STNCL); return expr(reg_tag, comp); }
	__INLINE expr funcval(INMOST_DATA_ENUM_TYPE reg_func) { assert(reg_func >= AD_FUNC); return expr(reg_func, ENUMUNDEF); }
#else

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
			if ((other.op >= AD_TAG && other.op < AD_STNCL) || other.op == AD_COND_TYPE || other.op == AD_COND_MARK) left = other.left; //copy component number
			else if (other.left != NULL) left = new expr(*other.left); else left = NULL;
			if (other.right != NULL) right = new expr(*other.right); else right = NULL;
		}
		expr(INMOST_DATA_REAL_TYPE val) : op(AD_CONST), coef(val), left(NULL), right(NULL) {}
		//expr(INMOST_DATA_ENUM_TYPE val) : op(AD_CONST), coef(val), left(NULL), right(NULL) {}
		expr(INMOST_DATA_ENUM_TYPE new_op, expr * l, expr * r) : op(new_op), coef(1), left(l), right(r) {}
		expr(INMOST_DATA_ENUM_TYPE tag_op, INMOST_DATA_ENUM_TYPE comp) :op(tag_op), coef(1), right(NULL) 
		{
			*(INMOST_DATA_ENUM_TYPE *)(&left) = comp;
		}
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
	__INLINE expr ad_pow(const expr & v, const expr n) { return expr(AD_POW, new expr(v), new expr(n)); }
	__INLINE expr ad_abs(const expr & v) { return expr(AD_ABS, new expr(v), NULL); }
	__INLINE expr ad_exp(const expr & v) { return expr(AD_EXP, new expr(v), NULL); }
	__INLINE expr ad_log(const expr & v) { return expr(AD_LOG, new expr(v), NULL); }
	__INLINE expr ad_sin(const expr & v) { return expr(AD_SIN, new expr(v), NULL); }
	__INLINE expr ad_cos(const expr & v) { return expr(AD_COS, new expr(v), NULL); }
	__INLINE expr ad_sqrt(const expr & v) {return expr(AD_SQRT, new expr(v), NULL);}
	__INLINE expr ad_val(const expr & v, const expr & multiplyer = expr(0.0)) {return expr(AD_VAL,new expr(v), new expr(multiplyer));}
	__INLINE expr measure() { return expr(AD_MES, NULL, NULL); }
	__INLINE expr condition_etype(ElementType etype, const expr & if_true, const expr & if_false) { return expr(AD_COND, new expr(AD_COND_TYPE,etype), new expr(AD_ALTR, new expr(if_true), new expr(if_false))); }
	__INLINE expr condition_marker(MarkerType marker, const expr & if_true, const expr & if_false) { return expr(AD_COND, new expr(AD_COND_MARK,marker), new expr(AD_ALTR, new expr(if_true), new expr(if_false))); }
	__INLINE expr condition(const expr & cond, const expr & if_gt_zero, const expr & if_le_zero) { return expr(AD_COND, new expr(cond), new expr(AD_ALTR, new expr(if_gt_zero), new expr(if_le_zero))); }
	__INLINE expr stencil(INMOST_DATA_ENUM_TYPE stncl, const expr & v) { assert(stncl >= AD_STNCL && stncl < AD_TABLE); return expr(stncl, new expr(v), NULL); }
	__INLINE expr tabval(INMOST_DATA_ENUM_TYPE tabl, const expr & v) { assert(tabl >= AD_TABLE && tabl < AD_FUNC); return expr(tabl, new expr(v), NULL); }
	__INLINE expr tagval(INMOST_DATA_ENUM_TYPE reg_tag, INMOST_DATA_ENUM_TYPE comp = 0) { assert(reg_tag >= AD_TAG && reg_tag < AD_STNCL); return expr(reg_tag, comp); }
	__INLINE expr funcval(INMOST_DATA_ENUM_TYPE reg_func) { assert(reg_func >= AD_FUNC); return expr(reg_func, NULL, NULL); }
#endif
}

#endif
#endif
