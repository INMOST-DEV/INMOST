#pragma once
#ifndef INMOST_AUTODIFF_H_INCLUDED
#define INMOST_AUTODIFF_H_INCLUDED
#include "inmost_common.h"
#include "inmost_mesh.h"
#include <sstream> //for debug
#include <CL/cl.h>


#if defined(USE_AUTODIFF) && (!defined(USE_MESH))
#warning "USE_AUTODIFF require USE_MESH"
#undef USE_AUTODIFF
#endif


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

	
		
	class expr
	{
		INMOST_DATA_ENUM_TYPE op;
		INMOST_DATA_REAL_TYPE coef;
		expr * left;
		expr * right;
	public:
		expr() : op(AD_NONE), coef(1), left(NULL), right(NULL) {}
		expr(const expr & other) : op(other.op), coef(other.coef) 
		{
			if (other.op >= AD_TAG && other.op < AD_STNCL) left = other.left; //copy component number
			else if( other.left != NULL ) left = new expr(*other.left); else left = NULL; 
			if( other.right != NULL ) right = new expr(*other.right); else right = NULL;
		}
		expr(INMOST_DATA_REAL_TYPE val) : op(AD_CONST), coef(val), left(NULL) , right(NULL) {}
		expr(INMOST_DATA_ENUM_TYPE new_op, expr * l, expr * r) : op(new_op), coef(1), left(l), right(r) {}
		expr(INMOST_DATA_ENUM_TYPE tag_op, INMOST_DATA_ENUM_TYPE comp = 0) :op(tag_op), coef(1), left(reinterpret_cast<expr *>(comp)), right(NULL) {}
		~expr() 
		{
			if( op < AD_COS ) 
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
			if (other.op >= AD_TAG && other.op < AD_STNCL) left = other.left; //copy component number
			else if( other.left != NULL ) left = new expr(*other.left); else left = NULL; 
			if( other.right != NULL ) right = new expr(*other.right); else right = NULL;
			return *this;
		}
		expr operator +() {return expr(*this);}
		expr operator -() {expr var(*this); var.coef *= -1.0; return var;}
		expr operator +(const expr & other) const {return expr(AD_PLUS ,new expr(*this),new expr(other));}
		expr operator -(const expr & other) const {return expr(AD_MINUS,new expr(*this),new expr(other));}
		expr operator *(const expr & other) const {return expr(AD_MULT ,new expr(*this),new expr(other));}
		expr operator /(const expr & other) const {return expr(AD_DIV  ,new expr(*this),new expr(other));}
		expr operator +(const INMOST_DATA_REAL_TYPE & other) const {return expr(AD_PLUS ,new expr(*this),new expr(other));}
		expr operator -(const INMOST_DATA_REAL_TYPE & other) const {return expr(AD_MINUS,new expr(*this),new expr(other));}
		expr operator *(const INMOST_DATA_REAL_TYPE & other) const {expr var(*this); var.coef *= other; return var;}
		expr operator /(const INMOST_DATA_REAL_TYPE & other) const {expr var(*this); var.coef /= other; return var;}
		__INLINE friend expr operator +(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		__INLINE friend expr operator -(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		__INLINE friend expr operator *(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		__INLINE friend expr operator /(const INMOST_DATA_REAL_TYPE & left, const expr & right);
		bool operator ==(const expr & other) {return op == other.op && left == other.left && right == other.right;}
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

	class Automatizator
	{

		static void __stdcall pfn_notify(const char * errinfo, const void * private_info, size_t cb, void * user_data)
		{
			std::cout << errinfo << std::endl;
		}
		typedef struct skope_t
		{
			int unknowns;
			INMOST_DATA_ENUM_TYPE stencil;
			bool measure_data;
			dynarray<INMOST_DATA_ENUM_TYPE, 32> function_data;
			dynarray<std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>, 32> static_data;
			dynarray<std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>, 32> dynamic_data;
			dynarray<skope_t *, 32> skopes;
		} skope;
		void AnalyzeSkopes(const expr & e, skope & data)
		{
			int ret = 0;
			switch (e.op)
			{
			case AD_PLUS:
				AnalyzeSkopes(*e.left, data);
				AnalyzeSkopes(*e.right, data);
				return;
			case AD_MINUS:
				AnalyzeSkopes(*e.left, data);
				AnalyzeSkopes(*e.right, data);
				return;
			case AD_MULT:
				AnalyzeSkopes(*e.left, data);
				AnalyzeSkopes(*e.right, data);
				return;
			case AD_DIV:
				AnalyzeSkopes(*e.left, data);
				AnalyzeSkopes(*e.right, data);
				return;
			case AD_POW:
				AnalyzeSkopes(*e.left, data);
				AnalyzeSkopes(*e.right, data);
				return;
			case AD_INV:
			case AD_ABS:
			case AD_EXP:
			case AD_LOG:
			case AD_SIN:
			case AD_COS:
				AnalyzeSkopes(*e.left, data);
				return;
			case AD_CONST:
			case AD_MES:
				return;
			case AD_COND:
				AnalyzeSkopes(*e.left, data);
				//skope * left = new skope, *right = new skope;
				//left->stencil = right->stencil = AD_COND;
				//left->unknowns = right->unknowns = 0;
				//data.skopes.push_back(left);
				//data.skopes.push_back(right);
				AnalyzeSkopes(*e.right->left, data);// *left);
				AnalyzeSkopes(*e.right->right, data);// *right);
				return;
			}
			if (e.op >= AD_FUNC)
			{
				for (dynarray<INMOST_DATA_ENUM_TYPE, 32>::iterator it = data.function_data.begin(); it != data.function_data.end(); ++it)
				{
					if (*it == e.op)
						return;
				}
				data.function_data.push_back(e.op);
				return;
			}
			if (e.op >= AD_TABLE)
			{
				AnalyzeSkopes(*e.left, data);
				return;
			}
			if (e.op >= AD_STNCL)
			{
				skope * s = new skope;
				s->stencil = e.op;
				s->unknowns = 0;
				data.skopes.push_back(s);
				AnalyzeSkopes(*e.left, *s);
				return;
			}
			if (e.op >= AD_CTAG)
			{
				std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> p(e.op, *(INMOST_DATA_ENUM_TYPE *)&e.left);
				for (dynarray<std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>, 32>::iterator it = data.static_data.begin(); it != data.static_data.end(); ++it)
				{
					if (*it == p)
						return;
				}
				data.static_data.push_back(p);
				return;
			}
			if (e.op >= AD_TAG)
			{
				data.unknowns++;
				std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> p(e.op, *(INMOST_DATA_ENUM_TYPE *)&e.left);
				for (dynarray<std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>, 32>::iterator it = data.dynamic_data.begin(); it != data.dynamic_data.end(); ++it)
				{
					if (*it == p)
						return;
				}
				data.dynamic_data.push_back(p);
				return;
			}
		}
		typedef struct execution_data_t
		{
			cl_context ctx;
			cl_uint nkernels;
			cl_kernel * kernels;
			cl_command_queue queue;
			cl_uint alloc;
			cl_mem buffers[2];
			char * local_buffers[2];
			std::vector<std::pair<expr, std::string> > expressions;
		} execution_data;
		std::vector<execution_data> exec_data;
		//structures to support static data types
		typedef struct{Tag t; MIDType domain_mask;} tagdomain;
		typedef small_hash<INMOST_DATA_ENUM_TYPE,tagdomain,128> const_tag_type;
		//data to support static data types
		const_tag_type reg_ctags;

		//structures to support dynamic data tags
		typedef struct{tagdomain d; Tag indices;} tagpair;
		typedef small_hash<INMOST_DATA_ENUM_TYPE,tagpair,128> tagpairs_type;
		typedef std::vector<tagpair> index_enum;
		//data to support dynamic data tags
		index_enum index_tags;
		tagpairs_type reg_tags;
		INMOST_DATA_ENUM_TYPE first_num;
		INMOST_DATA_ENUM_TYPE last_num;
		
		

		//structures to support tables
		typedef struct
		{
			std::string name;
			INMOST_DATA_REAL_TYPE * args;
			INMOST_DATA_REAL_TYPE * vals;
			INMOST_DATA_ENUM_TYPE size;
			INMOST_DATA_ENUM_TYPE binary_search(INMOST_DATA_REAL_TYPE arg)
			{
				int l = 0, r = static_cast<int>(size) - 1, mid = 0;
				while( r>=l )
				{
					mid = (l+r)/2;
					if( args[mid] > arg ) r = mid-1;
					else if( args[mid] < arg ) l = mid+1;
					else return mid;
				}
				mid = (l + r) / 2;
				if( mid > static_cast<int>(size-2) ) mid = static_cast<int>(size-2);
				return static_cast<INMOST_DATA_ENUM_TYPE>(mid);
			}
			INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE arg)
			{
				if (arg < args[0]) return vals[0];
				INMOST_DATA_ENUM_TYPE i = binary_search(arg);
				return vals[i] + (vals[i + 1] - vals[i]) * (arg - args[i]) / (args[i + 1] - args[i]);
			}
			INMOST_DATA_REAL_TYPE get_derivative(INMOST_DATA_REAL_TYPE arg)
			{
				if (arg < args[0]) return 0.0;
				INMOST_DATA_ENUM_TYPE i = binary_search(arg);
				return (vals[i + 1] - vals[i]) / (args[i + 1] - args[i]);
			}
		} table;
		typedef table * table_ptr;
		typedef small_hash<INMOST_DATA_ENUM_TYPE,table_ptr,32> table_type;
		//data to support tables
		table_type reg_tables;

		//structures to support stencil
	public:
		typedef std::pair<Storage *, INMOST_DATA_REAL_TYPE> stencil_pair;
		typedef dynarray<stencil_pair,32> stencil_pairs;
		typedef void (*stencil_callback)(Storage * current_element, stencil_pairs & out_stencil, void * user_data);
	private:
		typedef struct { Tag elements, coefs; } stencil_tag;
		typedef struct { std::string name; INMOST_DATA_ENUM_TYPE kind; MIDType domainmask; void * link; } stencil_kind_domain;
		typedef small_hash<INMOST_DATA_ENUM_TYPE,stencil_kind_domain,128> stencil_type;
		stencil_type reg_stencils;

	public:
		typedef INMOST_DATA_REAL_TYPE (*func_callback)(Storage * current_element, void * user_data);
	private:
		typedef struct func_name_callback_t { std::string name; func_callback func; } func_name_callback;
		typedef small_hash<INMOST_DATA_ENUM_TYPE, func_name_callback, 128> func_type;
		func_type reg_funcs;

		Mesh * m;
//	public:
//		typedef tiny_map<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE, 128> row_entries;
//	private:
		typedef dynarray<INMOST_DATA_REAL_TYPE,1024> precomp_values_t;
		std::string opname(const expr & var)
		{
			std::stringstream str;
			switch (var.op)
			{
			case AD_NONE: return "none";
			case AD_PLUS: return "+";
			case AD_MINUS: return "-";
			case AD_MULT: return "*";
			case AD_DIV: return "/";
			case AD_INV: return "1/";
			case AD_COND: return "?";
			case AD_POW: return "^";
			case AD_ABS: return "|.|";
			case AD_EXP: return "exp";
			case AD_LOG: return "log";
			case AD_SIN: return "sin";
			case AD_COS: return "cos";
			case AD_CONST: return "const";
			case AD_MES: return "measure";
			}
			if (var.op >= AD_FUNC)
			{
				str << "function(" << reg_funcs[var.op].name << ")";
				return str.str();
			}
			else if (var.op >= AD_TABLE)
			{
				str << "table(" << reg_tables[var.op]->name << ")";
				return str.str();
			}
			else if (var.op >= AD_STNCL)
			{
				str << "stencil(" << reg_stencils[var.op].name << ")";
				return str.str();
			}
			else if (var.op >= AD_CTAG)
			{
				str << "const_tag(" << GetStaticValueTag(var.op).GetTagName() << ")";
				return str.str();
			}
			else if (var.op >= AD_TAG)
			{
				str << "tag(" << GetDynamicValueTag(var.op).GetTagName() << ")";
				return str.str();
			}
			assert(false);
			return "unknown";
		}
		typedef dynarray<int, 1024> savevars_t;
		void GenerateTables(std::ostream & ostream)
		{
			int size = 0;
			for (table_type::iterator it = reg_tables.begin(); it != reg_tables.end(); ++it)
				size += static_cast<int>(it->second->size);
			std::stringstream values, arguments, offsets;
			values << "__constant real_t table_values[" << size << "]=" << std::endl << "{" << std::endl << "\t";
			arguments << "__constant real_t table_arguments[" << size << "]=" << std::endl << "{" << std::endl << "\t";
			offsets << "__constant int table_offset[" << reg_tables.size() + 1 << "]=" << std::endl << "{" << std::endl << "\t";
			bool first = true;
			size = 0;
			for (table_type::iterator it = reg_tables.begin(); it != reg_tables.end(); ++it)
			{
				if (first)
					first = false;
				else
				{
					values << ", " << std::endl << "\t";
					arguments << ", " << std::endl << "\t";
					offsets << ", " << std::endl << "\t";
				}
				values << "// values of " << it->second->name << " array entries " << size << ".." << size + it->second->size << std::endl << "\t";
				arguments << "// derivatives of " << it->second->name << " array entries " << size << ".." << size + it->second->size << std::endl << "\t";
				offsets << "// offset of " << it->second->name << std::endl << "\t";
				for (INMOST_DATA_ENUM_TYPE k = 0; k < it->second->size - 1; ++k)
				{
					values << std::scientific << it->second->vals[k] << ", ";
					arguments << std::scientific << it->second->args[k] << ", ";
					if ((k + 1) % 6 == 0)
					{
						values << std::endl << "\t";
						arguments << std::endl << "\t";
					}
				}
				values << std::scientific << it->second->vals[it->second->size - 1];
				arguments << std::scientific << it->second->args[it->second->size - 1];
				offsets << size;
				size += it->second->size;
			}
			values << std::endl << "};" << std::endl;
			arguments << std::endl << "};" << std::endl;
			offsets << "," << std::endl << "\t" << size;
			offsets << std::endl << "};" << std::endl;
			ostream << values.rdbuf();
			ostream << arguments.rdbuf();
			ostream << offsets.rdbuf();
			ostream << "int binary_search(int table_number, real_t argument)" << std::endl;
			ostream << "{" << std::endl;
			ostream << "\tint l = 0, r = table_offset[table_number+1]-table_offset[table_number] - 1, mid = 0;" << std::endl;
			ostream << "\twhile(r>=l)" << std::endl;
			ostream << "\t{" << std::endl;
			ostream << "\t\tmid=(l+r)/2;" << std::endl;
			ostream << "\t\tif(table_arguments[table_offset[table_number]+mid]>argument) r = mid-1;" << std::endl;
			ostream << "\t\telse if(table_arguments[table_offset[table_number]+mid]<argument) l = mid+1;" << std::endl;
			ostream << "\t\telse return mid;" << std::endl;
			ostream << "\t}" << std::endl;
			ostream << "\tif(mid>table_offset[table_number+1]-table_offset[table_number]-2) mid=table_offset[table_number+1]-table_offset[table_number]-2;" << std::endl;
			ostream << "\treturn mid;" << std::endl;
			ostream << "}" << std::endl;
			ostream << "real_t get_table_value(int table_number, real_t argument)" << std::endl;
			ostream << "{" << std::endl;
			ostream << "\tif(argument<table_arguments[table_offset[table_number]]) return table_values[table_offset[table_number]];" << std::endl;
			ostream << "\tint i = binary_search(table_number,argument);" << std::endl;
			ostream << "\treturn table_values[table_offset[table_number]+i]" << std::endl;
			ostream << "\t\t+(table_values[table_offset[table_number]+i+1]-table_values[table_offset[table_number]+i])" << std::endl;
			ostream << "\t\t*(argument - table_arguments[table_offset[table_number] + i])" << std::endl;
			ostream << "\t\t/(table_arguments[table_offset[table_number]+i+1]-table_arguments[table_offset[table_number]+i]);" << std::endl;
			ostream << "}" << std::endl;
			ostream << "real_t get_table_derivative(int table_number, real_t argument)" << std::endl;
			ostream << "{" << std::endl;
			ostream << "\tif(argument<table_arguments[table_offset[table_number]]) return 0.0;" << std::endl;
			ostream << "\tint i = binary_search(table_number,argument);" << std::endl;
			ostream << "\treturn (table_values[table_offset[table_number]+i+1]-table_values[table_offset[table_number]+i])" << std::endl;
			ostream << "\t\t/(table_arguments[table_offset[table_number]+i+1]-table_arguments[table_offset[table_number]+i]);" << std::endl;
			ostream << "}" << std::endl;
		}
		std::ostream & WriteTabs(std::ostream & out, int tabs)
		{
			for (int i = 0; i < tabs; i++) out << "\t";
			return out;
		}
		
		int GenerateEvaluateCode(const expr & var, std::ostream & ostream, int & varnum, int loopnum, savevars_t & savevars, skope & sdata, int tabs)
		{
			assert(var.op != AD_NONE);
			int v1, v2, v3, vret;
			savevars_t savestub;
			switch (var.op)
			{
			case AD_PLUS:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				v2 = GenerateEvaluateCode(*var.right, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=(" << "var" << v1 << "+" << "var" << v2 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				return vret;
			case AD_MINUS:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				v2 = GenerateEvaluateCode(*var.right, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=(" << "var" << v1 << "-" << "var" << v2 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				return vret;
			case AD_MULT:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				v2 = GenerateEvaluateCode(*var.right, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=(" << "var" << v1 << "*" << "var" << v2 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(v1);
				savevars.push_back(v2);
				//ostream << opname(var) << ":" << v1 << " " << v2 << " " << savevars.size() << std::endl;;
				return vret;
			case AD_DIV:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				v2 = GenerateEvaluateCode(*var.right, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=(" << "var" << v1 << "/" << "var" << v2 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(v1);
				savevars.push_back(v2);
				//ostream << opname(var) << ":" << v1 << " " << v2 << " " << savevars.size() << std::endl;;
				return vret;
			case AD_INV:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=(" << var.coef << "/" << "var" << v1 << ");" << std::endl;
				savevars.push_back(v1);
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;;
				return vret;
			case AD_POW:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				v2 = GenerateEvaluateCode(*var.right, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=pow(" << "var" << v1 << "," << "var" << v2 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(v1);
				savevars.push_back(v2);
				savevars.push_back(vret);
				//ostream << opname(var) << ":" << v1 << " " << v2 << " " << vret << " " << savevars.size() << std::endl;;
				return vret;
			case AD_ABS:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=abs(" << "var" << v1 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(v1);
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;;
				return vret;
			case AD_EXP:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=exp(" << "var" << v1 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(vret);
				//ostream << opname(var) << ":" << vret << " " << savevars.size() << std::endl;;
				return vret;
			case AD_LOG:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=log(" << "var" << v1 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(v1);
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;;
				return vret;
			case AD_SIN:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=sin(" << "var" << v1 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(v1);
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;;
				return vret;
			case AD_COS:
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=cos(" << "var" << v1 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				savevars.push_back(v1);
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;;
				return vret;
			case AD_CONST:
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=" << var.coef << "; // constant value" << std::endl;
				return vret;
			case AD_MES:
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=function_data[elem_" << loopnum << "]";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				return vret;
			case AD_COND:
				v3 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savestub, sdata, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "; // condition statement variable" << std::endl;
				WriteTabs(ostream, tabs) << "if(var" << v3 << ">0.0)" << std::endl;
				WriteTabs(ostream, tabs) << "{" << std::endl;
				v1 = GenerateEvaluateCode(*var.right->left, ostream, varnum, loopnum, savevars, sdata, tabs + 1);
				WriteTabs(ostream, tabs + 1) << "var" << vret << "=var" << v1;
				if (fabs(var.coef - 1.0) > 1e-13)  ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				WriteTabs(ostream, tabs) << "}" << std::endl;
				WriteTabs(ostream, tabs) << "else" << std::endl;
				WriteTabs(ostream, tabs) << "{" << std::endl;
				v2 = GenerateEvaluateCode(*var.right->right, ostream, varnum, loopnum, savevars, sdata, tabs + 1);
				WriteTabs(ostream, tabs + 1) << "var" << vret << "=var" << v2;
				if (fabs(var.coef - 1.0) > 1e-13)  ostream << "*" << var.coef;
				ostream << ";" << std::endl;
				WriteTabs(ostream, tabs) << "}" << std::endl;
				savevars.push_back(v3);
				//ostream << opname(var) << ":" << v3 << " " << savevars.size() << std::endl;;
				return vret;
			}
			if (var.op >= AD_FUNC)
			{
				vret = varnum;
				WriteTabs(ostream, tabs) << "real_t var" << varnum++ << "=function_data[" << (var.op - AD_FUNC) + 1 << "*nelems+elem_" << loopnum << "]";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";";
				ostream << "// function " << reg_funcs[var.op].name << std::endl;
				return vret;
			}
			if (var.op >= AD_TABLE)
			{
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum, savevars, tabs);
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=get_table_value(" << (var.op - AD_TABLE) << ",var" << v1 << ")";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << ";";
				ostream << "// table " << reg_tables[var.op]->name << std::endl;
				savevars.push_back(v1);
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;
				return vret;
			}
			if (var.op >= AD_STNCL)
			{
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=0.0; // initial variable for loop" << std::endl;
				WriteTabs(ostream, tabs) << "for(int iter=stncl_offset[" << (var.op - AD_STNCL) << "]; iter < stncl_offset[" << (var.op - AD_STNCL) + 1 << "]; iter++) // stencil " << reg_stencils[var.op].name << std::endl;
				WriteTabs(ostream, tabs) << "{" << std::endl;
				WriteTabs(ostream, tabs + 1) << "int elem_" << loopnum + 1 << "=stncl_elem[iter];" << std::endl;
				v1 = GenerateEvaluateCode(*var.left, ostream, varnum, loopnum + 1, savevars, tabs + 1);
				WriteTabs(ostream, tabs + 1) << "var" << vret << "=var" << vret << "+var" << v1 << "*stncl_coef[iter];" << std::endl;
				WriteTabs(ostream, tabs) << "}" << std::endl;
				if (fabs(var.coef - 1) > 1.0e-13)
				{
					WriteTabs(ostream, tabs) << "var" << vret << "=var" << vret << "*" << var.coef << ";" << std::endl;
				}
				savevars.push_back(vret);
				return vret;
			}
			if (var.op >= AD_CTAG)
			{
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=static_data[" << (var.op - AD_CTAG) << "*nelems+elem_" << loopnum << "]";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << "; // static tag " << GetStaticValueTag(var.op).GetTagName() << std::endl;
				return vret;
			}
			if (var.op >= AD_TAG)
			{
				vret = varnum++;
				WriteTabs(ostream, tabs) << "real_t var" << vret << "=dynamic_data[" << (var.op - AD_TAG) << "*nelems+elem_" << loopnum << "]";
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << "; // dynamic tag " << GetDynamicValueTag(var.op).GetTagName() << std::endl;
				return vret;
			}
			assert(false);
			return vret;
		}
		int GenerateDerivativeCode(const expr & var, std::ostream & ostream, int loopnum, savevars_t & savevars, std::vector<std::string> txtexpr, int tabs)
		{
			assert(var.op != AD_NONE);
			int v1, v2, v3, vret, vars = 0;
			std::stringstream str, temp, templ, tempr;
			switch (var.op)
			{
			case AD_PLUS:
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
				}
				vars += GenerateDerivativeCode(*var.right, ostream, loopnum, savevars, txtexpr, tabs);
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_MINUS:
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
				}
				txtexpr.push_back("*(-1)");
				vars += GenerateDerivativeCode(*var.right, ostream, loopnum, savevars, txtexpr, tabs);
				txtexpr.pop_back();
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_MULT:
				v2 = savevars.back(); savevars.pop_back();
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << v2 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*var" << v1;
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.right, ostream, loopnum, savevars, txtexpr, tabs);
				txtexpr.pop_back();
				str.clear();
				str << "*var" << v2;
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_DIV:
				v2 = savevars.back(); savevars.pop_back();
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << v2 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*var" << v1 << "/(var" << v2 << "*var" << v2 << ")";
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.right, ostream, loopnum, savevars, txtexpr, tabs);
				txtexpr.pop_back();
				str.clear();
				str << "/var" << v2;
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_INV:
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;
				str << "*" << var.coef << "/(var" << v1 << "*var" << v1 << ")";
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_POW:
				vret = savevars.back(); savevars.pop_back();
				v2 = savevars.back(); savevars.pop_back();
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << v2 << " " << vret << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*var" << vret;
				txtexpr.push_back(str.str());
				str.clear();
				str << "*var" << v2 << "/var" << v1;
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.right, ostream, loopnum, savevars, txtexpr, tabs);
				txtexpr.pop_back();
				str.clear();
				str << "*log(var" << v1 << ")";
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_ABS:
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*(var" << v1 << " > 0.0? 1.0 : -1.0)";
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_EXP:
				vret = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << vret << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*var" << vret;
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_LOG:
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "/var" << v1;
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_SIN:
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*cos(var" << v1 << ")";
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_COS:
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*(-sin(var" << v1 << "))";
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			case AD_CONST:
				return 0;
			case AD_MES:
				return 0;
			case AD_COND:
				savevars_t stub;
				v3 = savevars.back(); savevars.pop_back();
				vret = v3 + 1;
				GenerateEvaluateCode(*var.right->left, templ, vret, loopnum, stub, tabs + 1);
				GenerateEvaluateCode(*var.right->right, tempr, vret, loopnum, stub, tabs + 1);
				//ostream << opname(var) << ":" << v3 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
				}
				WriteTabs(temp, tabs) << "if(var" << v3 << ">0.0)" << std::endl;
				WriteTabs(temp, tabs) << "{" << std::endl;
				temp << tempr.rdbuf();
				vars += GenerateDerivativeCode(*var.right->right, temp, loopnum, savevars, txtexpr, tabs + 1);
				WriteTabs(temp, tabs) << "}" << std::endl;
				WriteTabs(temp, tabs) << "else" << std::endl;
				WriteTabs(temp, tabs) << "{" << std::endl;
				temp << templ.rdbuf();
				vars += GenerateDerivativeCode(*var.right->left, temp, loopnum, savevars, txtexpr, tabs + 1);
				WriteTabs(temp, tabs) << "}" << std::endl;
				if (vars > 0) ostream << temp.rdbuf();
				return vars;
			}
			if (var.op >= AD_FUNC)
			{
				return 0;
			}
			if (var.op >= AD_TABLE)
			{
				v1 = savevars.back(); savevars.pop_back();
				//ostream << opname(var) << ":" << v1 << " " << savevars.size() << std::endl;
				if (fabs(var.coef - 1.0) > 1e-13)
				{
					str << "*" << var.coef;
					txtexpr.push_back(str.str());
					str.clear();
				}
				str << "*get_table_derivative(" << (var.op - AD_TABLE) << ",var" << v1 << ")";
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, ostream, loopnum, savevars, txtexpr, tabs);
				return vars;
			}
			if (var.op >= AD_STNCL)
			{
				savevars_t stub;
				vret = savevars.back(); savevars.pop_back();
				vret = vret + 1;
				WriteTabs(ostream, tabs) << "for(int iter=stncl_offset[" << (var.op - AD_STNCL) << "]; iter < stncl_offset[" << (var.op - AD_STNCL) + 1 << "]; iter++) // stencil " << reg_stencils[var.op].name << std::endl;
				WriteTabs(temp, tabs) << "{" << std::endl;
				WriteTabs(temp, tabs + 1) << "int elem_" << loopnum + 1 << "=stncl_elem[iter];" << std::endl;
				GenerateEvaluateCode(*var.left, temp, vret, loopnum, stub, tabs + 1);
				str << "*stncl_coef[iter]";
				txtexpr.push_back(str.str());
				vars += GenerateDerivativeCode(*var.left, temp, loopnum + 1, savevars, txtexpr, tabs + 1);
				WriteTabs(temp, tabs) << "}" << std::endl;
				if (vars > 0) ostream << temp.rdbuf();
				return vars;
			}
			if (var.op >= AD_CTAG)
			{
				return 0;
			}
			if (var.op >= AD_TAG)
			{
				WriteTabs(ostream, tabs) << "derivative[" << (var.op - AD_TAG) << "*nelems+elem_" << loopnum << "] += 1";
				for (size_t k = 0; k < txtexpr.size(); k++) ostream << txtexpr[k];
				if (fabs(var.coef - 1.0) > 1e-13) ostream << "*" << var.coef;
				ostream << "; // derivative of " << GetDynamicValueTag(var.op).GetTagName() << std::endl;
				return 1;
			}
			assert(false);
			return 0;
		}
		//! returns value of computed expression, fills array of intermediate values
		INMOST_DATA_REAL_TYPE DerivativePrecompute(const expr & var, Storage * e, precomp_values_t & values, void * user_data, int print)
		{
			int newprint = print ? print + 1 : 0;
			assert(var.op != AD_NONE);
			INMOST_DATA_REAL_TYPE lval,rval,ret = 0.0;
			switch(var.op)
			{
			case AD_COND:  
				lval = Evaluate(*var.left,e,user_data); 
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if( print ) std::cout << " ( " << lval << " > 0.0 ? ";
				rval = DerivativePrecompute(*(lval > 0.0 ? var.right->left : var.right->right), e, values, user_data, newprint);
				values.push_back(lval); 
				if (print)  std::cout << " ) ";
				return rval*var.coef;
			case AD_PLUS:  
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " ( ";
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				if (print) std::cout << " + ";
				rval = DerivativePrecompute(*var.right, e, values, user_data, newprint);
				if (print) std::cout << " ) ";
				return (lval + rval)*var.coef;
			case AD_MINUS: 
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " ( ";
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				if (print) std::cout << " - ";
				rval = DerivativePrecompute(*var.right, e, values, user_data, newprint);
				if (print) std::cout << " ) ";
				return (lval - rval)*var.coef;
			case AD_MULT:
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " ( ";
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				if (print) std::cout << " * ";
				rval = DerivativePrecompute(*var.right, e, values, user_data, newprint);
				if (print) std::cout << " ) ";
				values.push_back(lval); 
				values.push_back(rval); 
				return (lval * rval)*var.coef;
			case AD_DIV:
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " ( ";
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				if (print) std::cout << " / ";
				rval = DerivativePrecompute(*var.right, e, values, user_data, newprint);
				if (print) std::cout << " ) ";
				values.push_back(lval); 
				values.push_back(rval); 
				return (lval / rval)*var.coef;
			case AD_POW:  
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				rval = DerivativePrecompute(*var.right, e, values, user_data, newprint);
				values.push_back(lval); 
				values.push_back(rval); 
				ret = ::pow(lval,rval); 
				values.push_back(ret);
				return ret*var.coef;
			case AD_INV:
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " ( 1.0 / ";
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				if (print) std::cout << " ) ";
				values.push_back(lval); 
				return var.coef / lval;
			case AD_ABS:   
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				values.push_back(lval); 
				return ::fabs(lval)*var.coef;
			case AD_EXP:   
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				ret = ::exp(lval);
				values.push_back(ret); 
				return ret*var.coef;
			case AD_LOG:   
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				values.push_back(lval); 
				return ::log(lval)*var.coef;
			case AD_SIN:   
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				values.push_back(lval); 
				return ::sin(lval)*var.coef;
			case AD_COS:   
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				values.push_back(lval); 
				return ::cos(lval)*var.coef;
			case AD_CONST: 
				if (print) std::cout << "const = " << var.coef;
				return var.coef;
			case AD_MES: 
				assert(!(e->GetElementType() & (ESET | MESH))); 
				m->GetGeometricData(static_cast<Element *>(e), MEASURE, &ret); 
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " vol = " << ret;
				return ret*var.coef;
			}
			if(var.op >= AD_FUNC )
			{
				ret =  reg_funcs[var.op].func(e,user_data);
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << reg_funcs[var.op].name << " = " << ret;
				return ret*var.coef;
			}
			else if(var.op >= AD_TABLE ) 
			{
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << reg_tables[var.op]->name << " [ ";
				lval = DerivativePrecompute(*var.left, e, values, user_data, newprint);
				values.push_back(lval); 
				ret = reg_tables[var.op]->get_value(lval);
				if (print) std::cout << " ] = " << ret;
				return ret*var.coef;
			}
			else if(var.op >= AD_STNCL ) 
			{
				stencil_kind_domain st = reg_stencils[var.op];
				assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " (" << st.name << ":";
				if( st.kind == 0 )
				{
					Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
					Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
					assert(elems.size() == coefs.size());
					for(INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k)
					{
						if (print) std::cout << " [ " << coefs[k] << " * ";
						lval = DerivativePrecompute(*var.left, elems[k], values, user_data, newprint);
						ret += lval * coefs[k];
						if (print) std::cout << " ] ";
					}
				}
				else if( st.kind == 1 )
				{
					stencil_pairs get_st;
					reinterpret_cast<stencil_callback>(st.link)(e,get_st,user_data);
					for(INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k)
					{
						if (print) std::cout << " [ " << get_st[k].second << " * ";
						lval = DerivativePrecompute(*var.left, get_st[k].first, values, user_data, newprint);
						ret += lval * get_st[k].second;
						if (print) std::cout << " ] ";
					}
				}
				if (print) std::cout << " ) ";
				return ret*var.coef;
			}
			else if (var.op >= AD_CTAG)
			{
				INMOST_DATA_ENUM_TYPE comp = *(INMOST_DATA_ENUM_TYPE *)(&var.left);
				ret = GetStaticValue(e, var.op, comp);
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " " << GetStaticValueTag(var.op).GetTagName() << "[" << comp << "] = " << ret;
				return ret * var.coef;
			}
			else if (var.op >= AD_TAG)
			{
				INMOST_DATA_ENUM_TYPE comp = *(INMOST_DATA_ENUM_TYPE *)(&var.left);
				ret = GetDynamicValue(e, var.op, comp);
				if (print) if (fabs(var.coef - 1.0) > 1e-3) std::cout << var.coef << " *";
				if (print) std::cout << " " << GetDynamicValueTag(var.op).GetTagName() << "[" << comp << "] = " << ret;
				return ret*var.coef;
			}
			assert(false);
			return 0.0;
		}
		//! returns offset from the end of precomputed values
		void DerivativeFill(const expr & var, Storage * e, Solver::Row & entries, precomp_values_t & values, INMOST_DATA_REAL_TYPE multval, void * user_data, int print)
		{
			int newprint = print ? print + 1 : 0;
			if (print)
			{
				for (int i = 1; i < print; i++) std::cout << "  ";
				std::cout << "entered " << opname(var) << " coef " << multval << std::endl;
			}
			assert(var.op != AD_NONE);
			INMOST_DATA_REAL_TYPE lval, rval, ret;
			switch(var.op)
			{
			case AD_COND:  
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*(lval > 0.0 ? var.right->left : var.right->right),e,entries,values,multval*var.coef,user_data, newprint); 
				goto exit;
			case AD_PLUS: 
				DerivativeFill(*var.right, e, entries, values, multval*var.coef, user_data,newprint);
				DerivativeFill(*var.left,e,entries,values,multval*var.coef,user_data,newprint);  
				goto exit;
			case AD_MINUS: 
				DerivativeFill(*var.right, e, entries, values, -multval*var.coef, user_data,newprint);
				DerivativeFill(*var.left, e, entries, values, multval*var.coef, user_data, newprint);
				goto exit;
			case AD_MULT:  
				rval = values.back(); values.pop_back(); 
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.right, e, entries, values, lval*multval*var.coef, user_data, newprint);
				DerivativeFill(*var.left, e, entries, values, rval*multval*var.coef, user_data, newprint);
				goto exit;
			case AD_DIV:   
				rval = values.back(); values.pop_back(); 
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.right, e, entries, values, -multval * lval / (rval*rval) * var.coef, user_data, newprint);
				DerivativeFill(*var.left, e, entries, values, multval / rval*var.coef, user_data, newprint);
				goto exit;
			case AD_POW:
				ret = values.back(); values.pop_back();
				rval = values.back(); values.pop_back(); 
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.right, e, entries, values, multval * ret * ::log(lval) * var.coef, user_data, newprint);
				DerivativeFill(*var.left, e, entries, values, multval * ret * rval / lval * var.coef, user_data, newprint);
				goto exit;
			case AD_INV:
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.left, e, entries, values, -multval / (lval * lval) * var.coef, user_data, newprint);
				goto exit;
			case AD_ABS:   
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.left, e, entries, values, multval * var.coef * (lval > 0 ? 1.0 : -1.0), user_data, newprint);
				goto exit;
			case AD_EXP:   
				ret = values.back(); values.pop_back(); 
				DerivativeFill(*var.left, e, entries, values, multval * ret * var.coef, user_data, newprint);
				goto exit;
			case AD_LOG:   
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.left, e, entries, values, multval * var.coef / lval, user_data, newprint);
				goto exit;
			case AD_SIN:   
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.left, e, entries, values, multval * var.coef * ::cos(lval), user_data, newprint);
				goto exit;
			case AD_COS:   
				lval = values.back(); values.pop_back(); 
				DerivativeFill(*var.left, e, entries, values, -multval * var.coef * ::sin(lval), user_data, newprint);
				goto exit;
			case AD_CONST: goto exit;
			case AD_MES: goto exit;
			}
			if(var.op >= AD_FUNC )
			{
				goto exit;
			}
			else if(var.op >= AD_TABLE ) 
			{
				lval = values.back(); values.pop_back();
				DerivativeFill(*var.left, e, entries, values, multval * var.coef * reg_tables[var.op]->get_derivative(lval), user_data, newprint);
				goto exit;
			}
			else if(var.op >= AD_STNCL ) 
			{
				stencil_kind_domain st = reg_stencils[var.op];
				assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
				if( st.kind == 0 )
				{
					Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
					Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
					assert(elems.size() == coefs.size());
					for (INMOST_DATA_ENUM_TYPE k = elems.size(); k > 0; --k)
						DerivativeFill(*var.left, elems[k - 1], entries, values, var.coef * coefs[k - 1] * multval, user_data, newprint);
				}
				else if( st.kind == 1 )
				{
					stencil_pairs get_st;
					reinterpret_cast<stencil_callback>(st.link)(e,get_st,user_data);
					for (INMOST_DATA_ENUM_TYPE k = get_st.size(); k > 0; --k)
						DerivativeFill(*var.left, get_st[k - 1].first, entries, values, var.coef * get_st[k - 1].second*multval, user_data, newprint);
				}
				goto exit;
			}
			else if (var.op >= AD_CTAG) goto exit;
			else if(var.op >= AD_TAG ) 
			{
				if (isDynamicValid(e, var.op))
				{
					INMOST_DATA_ENUM_TYPE ind = GetDynamicIndex(e, var.op, *(INMOST_DATA_ENUM_TYPE *)(&var.left));
					if (print)
					{
						for (int i = 1; i < print+1; i++) std::cout << "  ";
						std::cout << GetDynamicValueTag(var.op).GetTagName() << "[" << ind << "] +=" << multval * var.coef << std::endl;
					}
					entries[ind] += multval * var.coef;
				}
				goto exit;
			}
			assert(false);
		exit:
			if (print)
			{
				for (int i = 1; i < print; i++) std::cout << "  ";
				std::cout << "exited " << opname(var) << std::endl;
			}
			return;
		}
	public:
		
		INMOST_DATA_ENUM_TYPE GetFirstIndex() {return first_num;}
		INMOST_DATA_ENUM_TYPE GetLastIndex() {return last_num;}
		Automatizator(Mesh * m) :first_num(0),last_num(0),m(m) {}
		~Automatizator() 
		{
			for(unsigned k = 0; k < index_tags.size(); k++)
				index_tags[k].indices = m->DeleteTag(index_tags[k].indices);
			for(table_type::iterator it = reg_tables.begin(); it != reg_tables.end(); ++it)
			{
				delete [] it->second->args;
				delete [] it->second->vals;
				delete it->second;
			}
			for(stencil_type::iterator it = reg_stencils.begin(); it != reg_stencils.end(); ++it)
				if( it->second.kind == 0 )
					delete static_cast<stencil_tag *>(it->second.link);
		}
		// register function that calculates value for current element,
		// value should not depend on dynamic tags
		INMOST_DATA_ENUM_TYPE RegisterFunc(std::string name, func_callback func)
		{
			INMOST_DATA_ENUM_TYPE ret = reg_funcs.size() + AD_FUNC;
			func_name_callback v;
			v.name = name;
			v.func = func;
			reg_funcs[ret] = v;
			return ret;
		}
		//register stencil that can be got from tags
		INMOST_DATA_ENUM_TYPE RegisterStencil(std::string name, Tag elements_tag, Tag coefs_tag, MIDType domain_mask = 0)
		{
			INMOST_DATA_ENUM_TYPE ret = reg_stencils.size() + AD_STNCL;
			stencil_kind_domain st;
			stencil_tag * save = new stencil_tag;
			st.name = name;
			save->coefs = coefs_tag;
			save->elements = elements_tag;
			st.kind = 0;
			st.link = static_cast<void *>(save);
			st.domainmask = domain_mask;
			reg_stencils[ret] = st;
			return ret;
		}
		//register stencil that can be got from function
		INMOST_DATA_ENUM_TYPE RegisterStencil(std::string name, stencil_callback func, MIDType domain_mask = 0)
		{
			INMOST_DATA_ENUM_TYPE ret = reg_stencils.size() + AD_STNCL;
			stencil_kind_domain st;
			st.name = name;
			st.kind = 1;
			st.link = reinterpret_cast<void *>(func);
			st.domainmask = domain_mask;
			reg_stencils[ret] = st;
			return ret;
		}
		INMOST_DATA_ENUM_TYPE RegisterTable(std::string name, INMOST_DATA_REAL_TYPE * Arguments, INMOST_DATA_REAL_TYPE * Values, INMOST_DATA_ENUM_TYPE size)
		{
			INMOST_DATA_ENUM_TYPE ret = reg_tables.size() + AD_TABLE;
			table_ptr t = new table;
			t->name = name;
			t->args = new INMOST_DATA_REAL_TYPE [size];
			memcpy(t->args,Arguments,sizeof(INMOST_DATA_REAL_TYPE)*size);
			t->vals = new INMOST_DATA_REAL_TYPE [size];
			memcpy(t->vals,Values,sizeof(INMOST_DATA_REAL_TYPE)*size);
			t->size = size;
			reg_tables[ret] = t;
			return ret;
		}
		/// set data of tag t defined on domain_mask to be dynamic data
		/// don't register tag twice
		INMOST_DATA_ENUM_TYPE RegisterDynamicTag(Tag t, ElementType typemask, MIDType domain_mask = 0)
		{
			tagpair p;
			p.d.domain_mask = domain_mask;
			p.d.t = t;
			ElementType def = NONE, sparse = NONE;
			for(ElementType q = NODE; q <= MESH; q = q << 1) if( q & typemask )
			{
				if( t.isDefined(q) ) def |= q;
				if( t.isSparse(q) ) sparse |= q;
			}
			p.indices = m->CreateTag(t.GetTagName()+"_index",DATA_INTEGER,def,sparse,t.GetSize());
			INMOST_DATA_ENUM_TYPE ret = reg_tags.size() + AD_TAG;
			reg_tags[ret] = p;
			index_tags.push_back(p);
			return ret;
		}
		/// set index for every data entry of dynamic tag
		void EnumerateDynamicTags()
		{
			first_num = last_num = 0;
			const ElementType paralleltypes = NODE | EDGE | FACE | CELL;
			for(index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
				for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
					if( it->indices.isDefined(etype) && it->indices.isSparse(etype) )
						for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
							jt->DelData(it->indices);
					

			for(index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
			{
				for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
					if( it->indices.isDefined(etype) )
					{
						if(it->indices.GetSize() == ENUMUNDEF )
						{
							if( !it->indices.isSparse(etype) )
							{
								for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
									if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)) )
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										indarr.resize(jt->RealArray(it->d.t).size());
										for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt = last_num++;
									}
							}
							else
							{
								for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
									if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))) )
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										indarr.resize(jt->RealArray(it->d.t).size());
										for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt = last_num++;
									}
							}
						}
						else
						{
							if( !it->indices.isSparse(etype) )
							{
								for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
									if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)) )
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt = last_num++;
									}
							}
							else
							{
								for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
									if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))) )
									{
										Storage::integer_array indarr = jt->IntegerArray(it->indices);
										for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
											*qt = last_num++;
									}
							}
						}
					}
			}
#if defined(USE_MPI)
			if( m->GetProcessorsNumber() > 1 )
			{
				MPI_Scan(&last_num,&first_num,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,m->GetCommunicator());
				first_num -= last_num;
				ElementType exch_mask = NONE;
				if( first_num > 0 ) for(index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it)
				{
					for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
						if( it->indices.isDefined(etype) )
						{
							exch_mask |= etype;
							if(it->indices.GetSize() == ENUMUNDEF )
							{
								if( !it->indices.isSparse(etype) )
								{
									for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
										if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)) )
										{
											Storage::integer_array indarr = jt->IntegerArray(it->indices);
											for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
												*qt += first_num;
										}
								}
								else
								{
									for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
										if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))) )
										{
											Storage::integer_array indarr = jt->IntegerArray(it->indices);
											indarr.resize(jt->RealArray(it->d.t).size());
											for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
												*qt += first_num;
										}
								}
							}
							else
							{
								if( !it->indices.isSparse(etype) )
								{
									for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
										if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask)) )
										{
											Storage::integer_array indarr = jt->IntegerArray(it->indices);
											for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
												*qt += first_num;
										}
								}
								else
								{
									for(Mesh::base_iterator jt = m->Begin(etype); jt != m->End(); ++jt)
										if( ((etype & paralleltypes) && static_cast<Element *>(&*jt)->GetStatus() != Element::Ghost) && (jt->HaveData(it->d.t) && (it->d.domain_mask == 0 || jt->GetMarker(it->d.domain_mask))) )
										{
											Storage::integer_array indarr = jt->IntegerArray(it->indices);
											for(Storage::integer_array::iterator qt = indarr.begin(); qt != indarr.end(); ++qt)
												*qt += first_num;
										}
								}
							}
						}
				}
				last_num+=first_num;
				{
					std::vector<Tag> exch_tags;
					for(index_enum::iterator it = index_tags.begin(); it != index_tags.end(); ++it) exch_tags.push_back(it->indices);
					m->ExchangeData(exch_tags,exch_mask);
				}
			}
#endif

		}
		/// register tag, data for which don't change through iterations
		/// don't register tag twice
		INMOST_DATA_ENUM_TYPE RegisterStaticTag(Tag t, MIDType domain_mask = 0)
		{
			INMOST_DATA_ENUM_TYPE ret = reg_ctags.size() + AD_CTAG;
			tagdomain d;
			d.t = t;
			d.domain_mask = domain_mask;
			reg_ctags[ret] = d;
			return ret;
		}
		Tag					GetDynamicValueTag(INMOST_DATA_ENUM_TYPE ind) {return reg_tags[ind].d.t;}
		Tag					GetDynamicIndexTag(INMOST_DATA_ENUM_TYPE ind) {return reg_tags[ind].indices;}
		MIDType             GetDynamicMask(INMOST_DATA_ENUM_TYPE ind) {return reg_tags[ind].d.domain_mask;}
		Tag					GetStaticValueTag(INMOST_DATA_ENUM_TYPE ind) {return reg_ctags[ind].t;}
		MIDType             GetStaticMask(INMOST_DATA_ENUM_TYPE ind) {return reg_ctags[ind].domain_mask;}

		INMOST_DATA_REAL_TYPE GetDynamicValue(Storage * e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) {return e->RealArray(GetDynamicValueTag(ind))[comp];}
		INMOST_DATA_ENUM_TYPE GetDynamicIndex(Storage * e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) {return e->IntegerArray(GetDynamicIndexTag(ind))[comp];}
		bool                isDynamicValid(Storage * e, INMOST_DATA_ENUM_TYPE ind) {MIDType mask = GetDynamicMask(ind); return mask == 0 || e->GetMarker(mask);} 

		INMOST_DATA_REAL_TYPE GetStaticValue(Storage * e, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE comp = 0) {return e->RealArray(GetStaticValueTag(ind))[comp];}
		bool                isStaticValid(Storage * e, INMOST_DATA_ENUM_TYPE ind) {MIDType mask = GetStaticMask(ind); return mask == 0 || e->GetMarker(mask);} 
		//! Return value of given tree of expr
		INMOST_DATA_REAL_TYPE Evaluate(const expr & var, Storage * e, void * user_data)
		{
			assert(var.op != AD_NONE);
			switch(var.op)
			{
			case AD_COND:  return Evaluate(*(Evaluate(*var.left,e,user_data) > 0.0 ? var.right->left : var.right->right),e,user_data)*var.coef;
			case AD_PLUS:  return (Evaluate(*var.left,e,user_data) + Evaluate(*var.right,e,user_data))*var.coef;
			case AD_MINUS: return (Evaluate(*var.left,e,user_data) - Evaluate(*var.right,e,user_data))*var.coef;
			case AD_MULT:  return (Evaluate(*var.left,e,user_data) * Evaluate(*var.right,e,user_data))*var.coef;
			case AD_DIV:   return (Evaluate(*var.left,e,user_data) / Evaluate(*var.right,e,user_data))*var.coef;
			case AD_INV:   return var.coef / Evaluate(*var.left,e,user_data);
			case AD_POW:   return ::pow(Evaluate(*var.left,e,user_data),Evaluate(*var.right,e,user_data))*var.coef;
			case AD_ABS:   return ::fabs(Evaluate(*var.left,e,user_data))*var.coef;
			case AD_EXP:   return ::exp(Evaluate(*var.left,e,user_data))*var.coef;
			case AD_LOG:   return ::log(Evaluate(*var.left,e,user_data))*var.coef;
			case AD_SIN:   return ::sin(Evaluate(*var.left,e,user_data))*var.coef;
			case AD_COS:   return ::cos(Evaluate(*var.left,e,user_data))*var.coef;
			case AD_CONST: return var.coef;
			case AD_MES: assert(!(e->GetElementType() & (ESET | MESH))); Storage::real ret; m->GetGeometricData(static_cast<Element *>(e),MEASURE,&ret); return ret*var.coef;
			}
			if(var.op >= AD_FUNC ) return reg_funcs[var.op].func(e,user_data);
			if(var.op >= AD_TABLE ) return reg_tables[var.op]->get_value(Evaluate(*var.left,e,user_data))*var.coef;
			if(var.op >= AD_STNCL ) 
			{
				INMOST_DATA_REAL_TYPE ret = 0.0;
				stencil_kind_domain st = reg_stencils[var.op];
				assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
				if( st.kind == 0 )
				{
					Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
					Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
					assert(elems.size() == coefs.size());
					for(INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k)
						ret += var.coef * Evaluate(*var.left,elems[k],user_data) * coefs[k];
				}
				else if( st.kind == 1 )
				{
					stencil_pairs get_st;
					reinterpret_cast<stencil_callback>(st.link)(e,get_st, user_data);
					for(INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k)
						ret += var.coef * Evaluate(*var.left,get_st[k].first,user_data) * get_st[k].second;
				}
				return ret;
			}
			if(var.op >= AD_CTAG ) return GetStaticValue(e,var.op,*(INMOST_DATA_ENUM_TYPE *)(&var.left))*var.coef;
			if(var.op >= AD_TAG ) return  GetDynamicValue(e,var.op,*(INMOST_DATA_ENUM_TYPE *)(&var.left))*var.coef;
			assert(false);
			return 0.0;
		}
		//! Fill the row with expr derivatives and return value of expr
		INMOST_DATA_REAL_TYPE Derivative(const expr & var, Storage * e, Solver::Row & out, void * user_data, bool print = false)
		{
			INMOST_DATA_REAL_TYPE ret;
			precomp_values_t values;
			ret = DerivativePrecompute(var,e,values,user_data,print);
			if( print ) std::cout << std::endl;
			DerivativeFill(var,e,out,values,1.0,user_data,print);
			return ret;
		}
		//! Return the solution for given tag index of given element
		INMOST_DATA_REAL_TYPE GetIndex(Storage * e, INMOST_DATA_ENUM_TYPE tagind, INMOST_DATA_ENUM_TYPE comp = 0) {return e->IntegerArray(GetDynamicIndexTag(tagind))[comp];}
		//! Retrive number of components on given element for given tag
		INMOST_DATA_ENUM_TYPE GetComponents(Storage *e, INMOST_DATA_ENUM_TYPE tagind) {return e->IntegerArray(GetDynamicIndexTag(tagind)).size();}

		enum Precision { Double, Single };
		std::vector<skope *> GenerateCode(std::vector< std::pair<expr, std::string> > & expressions ,std::ostream & ostream, Precision prec)
		{
			std::stringstream header, body;
			int vnum = 0, vret;
			std::vector<std::string> txtexpr;
			std::vector<skope *> ret(expressions.size());
			savevars_t vars;
			if (prec == Single) header << "typedef float real_t;" << std::endl;
			else if (prec == Double)
			{
				header << "#pragma OPENCL EXTENSION cl_khr_fp64: enable" << std::endl;
				header << "typedef double real_t;" << std::endl;
			}
			else throw - 1;
			GenerateTables(header);
			for (size_t k = 0; k < expressions.size(); k++)
			{
				AnalyzeSkopes(expressions[k].first,ret[k])
				vnum = 0;
				header.clear();
				body.clear();
				header << "__kernel void " << expressions[k].second << std::endl;
				header << "\t\t\t\t(" << std::endl;
				header << "\t\t\t\t";
				header << "__global const int * offsets," << std::endl;
				header << "\t\t\t\t";
				header << "__global char * global_work_data," << std::endl;
				header << "\t\t\t\t)" << std::endl;
				body << "{" << std::endl;
				body << "\tint gid = get_global_id(0);" << std::endl;
				body << "\tchar * work_data = global_work_data + offsets[gid];" << std::endl;
				vret = GenerateEvaluateCode(expressions[k].first, body, vnum, 0, vars, 1);
				GenerateDerivativeCode(expressions[k].first, body, 0, vars, txtexpr, 1);
				body << "\t((real_t *)work_data)=var" << vret << ";" << std::endl;
				body << "}" << std::endl;
				ostream << header.rdbuf();
				ostream << body.rdbuf();
			}
			return ret;
		}
		template<typename real_t, typename int_t>
		int PrepareDataInner(Storage * elem, void * user_data, skope & sdata, char * buf, long pos, long off, long end)
		{
			if (sdata.measure_data)
			{
				off += sizeof(real_t); //remember measure data
				if (pos + off >= end) return 0;
				Storage::real measure;
				m->GetGeometricData(elem, MEASURE, &measure);
				((real_t)buf + pos + off - sizeof(real_t))[0] = static_cast<real_t>(measure);
			}
			for (dynarray<INMOST_DATA_ENUM_TYPE, 32>::iterator it = sdata.function_data.begin(); it != sdata.function_data.end(); ++it)
			{
				off += sizeof(real_t); //remember function data
				if (pos + off >= end) return 0;
				((real_t)buf + pos + off - sizeof(real_t))[0] = static_cast<real_t>(reg_funcs[*it].func(elem,user_data));
			}
			for (dynarray<std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>, 32>::iterator it = sdata.static_data.begin(); it != sdata.static_data.end(); ++it)
			{
				off += sizeof(real_t); //remember static data
				if (pos + off >= end) return 0;
				((real_t)buf + pos + off - sizeof(real_t))[0] = static_cast<real_t>(GetStaticValue(elem,it->first,it->second));
			}
			for (dynarray<std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>, 32>::iterator it = sdata.dynamic_data.begin(); it != sdata.dynamic_data.end(); ++it)
			{
				off += sizeof(real_t); //remember dynamic data
				if (pos + off >= end) return 0;
				((real_t)buf + pos + off - sizeof(real_t))[0] = static_cast<real_t>(GetDynamicValue(elem, it->first, it->second));
				off += sizeof(real_t); //remember dynamic data derivative
				if (pos + off >= end) return 0;
				((real_t)buf + pos + off - sizeof(real_t))[0] = static_cast<real_t>(0.0);
				off += sizeof(int_t); //remember dynamic data index
				if (pos + off >= end) return 0;
				((int_t)buf + pos + off - sizeof(int_t))[0] = static_cast<int_t>(GetDynamicIndex(elem, it->first, it->second));
			}
			for (dynarray<skope *, 32>::iterator it = sdata.stencil.begin(); it != sdata.stencil.end(); ++it)
			{
				stencil_kind_domain st = reg_stencils[(*it)->stencil];
				assert(st.domainmask == 0 || e->GetMarker(st.domainmask));
				stencil_pairs get_st;
				if (st.kind == 0)
				{
					Storage::reference_array elems = e->ReferenceArray(static_cast<stencil_tag *>(st.link)->elements);
					Storage::real_array coefs = e->RealArray(static_cast<stencil_tag *>(st.link)->coefs);
					assert(elems.size() == coefs.size());
					for (INMOST_DATA_ENUM_TYPE k = 0; k < elems.size(); ++k) get_st.push_back(std::make_pair(elems[k], coefs[k]));
				}
				else if (st.kind == 1) reinterpret_cast<stencil_callback>(st.link)(e, get_st, user_data);

				off += sizeof(int_t); //remember size of loop
				if (pos + off >= end) return 0;
				((int_t)buf + pos + off - sizeof(int_t))[0] = static_cast<int_t>(get_st.size());
				for (INMOST_DATA_ENUM_TYPE k = 0; k < get_st.size(); ++k)
				{
					off += sizeof(real_t);
					if (pos + off >= end) return 0;
					((real_t)buf + pos + off - sizeof(real_t))[0] = static_cast<real_t>(get_st[k].second); //remember coefficient
					long ret = PrepareDataInner(get_st[k].first, user_data, *it, buf, pos, off, end, elemkey);
					if (ret == 0) return 0;
					off += ret;
				}
			}
		}
		void CompileCode(std::vector<std::pair<expr, std::string> > & expressions, Precision prec)
		{
			char * program_source;
			{
				std::fstream test("out.cl", std::ios::out);
				std::stringstream source;
				GenerateCode(expressions, source, prec);
				
				std::string source_str = source.str();
				program_source = new char[source_str.size()+1];
				memcpy(program_source, source_str.c_str(), source_str.size());
				program_source[source_str.size()] = '\0';
				test << program_source;
			}
			cl_int status;
			cl_uint numPlatforms = 0;
			status = clGetPlatformIDs(0, NULL, &numPlatforms);
			if (status != CL_SUCCESS) { std::cout << "OpenCL error" << std::endl; throw(-1); }
			//std::cout << "Number of platforms: " << numPlatforms << std::endl;
			cl_platform_id * platforms = new cl_platform_id[numPlatforms];
			status = clGetPlatformIDs(numPlatforms, platforms, NULL);
			if (status != CL_SUCCESS) { std::cout << "OpenCL error" << std::endl; throw(-1); }
			for (cl_uint k = 0; k < numPlatforms; k++)
			{

				cl_uint numDevices = 0;
				status = clGetDeviceIDs(platforms[k], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
				//std::cout << "Number of devices on platform " << k << ": " << numDevices << std::endl;
				cl_device_id * devices = new cl_device_id[numDevices];
				status = clGetDeviceIDs(platforms[k], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
				if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }

				for (cl_uint q = 0; q < numDevices; q++)
				{
					char devname[4096];
					size_t inp = 4096, outp;
					status = clGetDeviceInfo(devices[q], CL_DEVICE_NAME, inp, devname, &outp);
					if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
					//std::cout << "Device " << q << " name: " << devname << std::endl;
					if (prec == Double)
					{
						cl_uint fp64;
						status = clGetDeviceInfo(devices[q], CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(cl_uint), (void *)&fp64, NULL);
						if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
						if (fp64 == 0)
						{
							std::cout << "Device " << devname << " don't support double precision " << std::endl;
							continue;
						}
					}
					cl_ulong alloc;
					status = clGetDeviceInfo(devices[q], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), (void *)&alloc, NULL);
					if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
					std::cout << "Can allocate " << alloc << " bytes on device " << devname << std::endl;
					alloc = alloc;
					cl_context context = clCreateContext(NULL, 1, &devices[q], pfn_notify, NULL, &status);
					if (status == CL_INVALID_PLATFORM) std::cout << "bad platform" << std::endl;
					else if (status == CL_INVALID_VALUE) std::cout << "invalid value" << std::endl;
					else if (status == CL_INVALID_DEVICE) std::cout << "invalid device" << std::endl;
					else if (status == CL_DEVICE_NOT_AVAILABLE) std::cout << "device not available" << std::endl;
					else if (status == CL_OUT_OF_HOST_MEMORY) std::cout << "no more memory" << std::endl;
					if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
					cl_program program = clCreateProgramWithSource(context, 1, (const char **)&program_source, NULL, &status);
					if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
					double t = Timer();
					status = clBuildProgram(program, 1, &devices[q], NULL, NULL, NULL);
					std::cout << "Build program time: " << Timer()-t << std::endl;
					if (status != CL_SUCCESS) 
					{ 
						char str[24576];
						size_t size_ret;
						status = clGetProgramBuildInfo(program, devices[q], CL_PROGRAM_BUILD_LOG, 24576, str, &size_ret);
						if (size_ret > 24576) std::cout << "not enough space for log" << std::endl;
						else std::cout << str << std::endl;
						std::cout << "Build failed for " << devname << std::endl;
					}
					else
					{
						std::cout << "Build succeded for " << devname << std::endl;
						
						execution_data r;
						r.alloc = alloc;
						r.ctx = context;
						r.queue = clCreateCommandQueue(context, devices[q], 0, &status);
						if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
						r.buffers[0] = clCreateBuffer(context, CL_MEM_READ_ONLY, alloc / 32, NULL, &status);
						r.local_buffers[0] = (char *)malloc(alloc / 32);
						if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
						r.buffers[1] = clCreateBuffer(context, CL_MEM_READ_WRITE, alloc, NULL, &status);
						if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
						r.local_buffers[0] = (char *)malloc(alloc);
						r.nkernels = expressions.size();
						r.kernels = static_cast<cl_kernel *>(malloc(sizeof(cl_kernel)*r.nkernels));
						r.expressions = expressions;
						t = Timer();
						for (size_t k = 0; k < expressions.size(); ++k)
						{
							r.kernels[k] = clCreateKernel(program, &expressions[k].second[0], &status);
							if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
							status = clSetKernelArg(r.kernels[k], 0, sizeof(cl_mem), r.buffers + 0);
							if (status != CL_SUCCESS) { std::cout << __FILE__ << ":" << __LINE__ << " " << "OpenCL error" << std::endl; throw(-1); }
							status = clSetKernelArg(r.kernels[k], 1, sizeof(cl_mem), r.buffers + 1);
						}
						std::cout << "Create kernels: " << Timer() - t << std::endl;
						exec_data.push_back(r);
					}
				}
				delete[] devices;
			}
			delete[] platforms;
			delete[] program_source;
		}
		
	};
}


#endif
#endif
