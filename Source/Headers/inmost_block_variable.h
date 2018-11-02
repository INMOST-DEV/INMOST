
#ifndef INMOST_AUTODIFF_ETBVAR_H_INCLUDED
#define INMOST_AUTODIFF_ETBVAR_H_INCLUDED
#include "inmost_common.h"
#include "inmost_expression.h"
#include "inmost_mesh.h"
#include "inmost_autodiff.h"
#include "inmost_solver.h"
#include "inmost_variable.h"
#include <sstream> //for debug
#include <new>

#if defined(USE_AUTODIFF) && defined(USE_MESH)

//TODO:
// 1. Add operations ConcatRows, ConcatCols that help assemble matrix from pieces


//This should stop Visual Studio from complaining of very long auto-generated class types
#ifdef _MSC_VER
#pragma warning(disable : 4503)
#endif



namespace INMOST
{
	
	class abstract_dynamic_block_variable
	{
	public:
		virtual rMatrix Value (const Storage & e) const = 0;
		virtual vMatrix Variable(const Storage & e) const = 0;
		virtual abstract_dynamic_block_variable * Copy() const = 0;
		virtual ~abstract_dynamic_block_variable() {}
	};
	
	template<typename RetType>
	class get_block_variable
	{
	public:
		virtual RetType operator()(const Storage & e) const = 0;
	};
	
	template<>
	class get_block_variable<vMatrix>
	{
		const abstract_dynamic_block_variable & var;
	public:
		typedef vMatrix type;
		get_block_variable(const abstract_dynamic_block_variable & var) : var(var) {}
		vMatrix operator()(const Storage & e) const {return var.Variable(e);}
	};
	
	template<>
	class get_block_variable<INMOST_DATA_REAL_TYPE>
	{
		const abstract_dynamic_block_variable & var;
	public:
		typedef rMatrix type;
		get_block_variable(const abstract_dynamic_block_variable & var) : var(var) {}
		rMatrix operator()(const Storage & e) const {return var.Value(e);}
	};
	
	
	class stored_block_variable_expression : public abstract_dynamic_block_variable
	{
		abstract_dynamic_block_variable * var;
	public:
		stored_block_variable_expression() : var(NULL) {}
		stored_block_variable_expression(const abstract_dynamic_block_variable & pvar) : var(pvar.Copy()) {}
		stored_block_variable_expression(const stored_block_variable_expression & other) : var(other.var->Copy()) {}
		~stored_block_variable_expression() {delete var; var = NULL;}
		stored_block_variable_expression operator =(stored_block_variable_expression const & other) {var = other.var->Copy(); return *this;}
		stored_block_variable_expression operator =(const abstract_dynamic_block_variable & pvar) {var = pvar.Copy(); return *this;}
		rMatrix Value(const Storage & e) const {return var->Value(e);}
		vMatrix Variable(const Storage & e) const {return var->Variable(e);}
		
		template<typename T>
		get_variable<T> get_variable() {return get_variable<T>(*var);}
		abstract_dynamic_block_variable & retrive_expression() {return *var;}
		const abstract_dynamic_block_variable & retrive_expression() const {return *var;}
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new stored_block_variable_expression(*this));}
		/// Checks that the stored expresison was defined.
		bool isDefined() const {return var != NULL;}
	};
	
	
	class dynamic_block_variable : public abstract_dynamic_block_variable
	{
	private:
		const AbstractEntry * entry;
	public:
		dynamic_block_variable() :entry(NULL) {}
		dynamic_block_variable(Automatizator & aut, INMOST_DATA_ENUM_TYPE reg_index) : entry(reg_index==ENUMUNDEF?NULL:&aut.GetEntry(reg_index)) {}
		dynamic_block_variable(const AbstractEntry * re) : entry(re) {}
		dynamic_block_variable(const dynamic_block_variable & other) : entry(other.entry) {}
		dynamic_block_variable & operator =(const dynamic_block_variable & other)
		{
			entry = other.entry;
			return * this;
		}
		rMatrix Value(const Storage & e) const {return entry->Value(e);}
		//iMatrix Index(const Storage & e) const {return entry->isValid(e) ? entry->Index(e):iMatrix(entry->MatrixSize(e),1,ENUMUNDEF);}
		vMatrix Variable(const Storage & e) const {return entry->Unknown(e);}
		bool isUnknown(const Storage & e) const {return entry->isValid(e)?true:false;}
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new dynamic_block_variable(*this));}
	};
	
	class const_block_variable : public abstract_dynamic_block_variable
	{
	private:
		rMatrix value;
	public:
                const_block_variable(const rMatrix & _value) : value(_value)  {}
		const_block_variable(const const_block_variable & other) : value(other.value) {}
		const_block_variable & operator =(const const_block_variable & other)
		{
			value = other.value;
			return * this;
		}
                rMatrix Value(const Storage & e) const {(void)e; return value;}
                vMatrix Variable(const Storage & e) const {(void)e; return value;}
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new const_block_variable(*this));}
	};
	
	class static_block_variable : public abstract_dynamic_block_variable
	{
	private:
		TagRealArray value_tag;
		INMOST_DATA_ENUM_TYPE n,m;
	public:
		static_block_variable(TagRealArray t)
		: value_tag(t), n(t.GetSize()), m(1) {assert(t.GetDataType() == DATA_REAL);}
		static_block_variable(TagRealArray t, INMOST_DATA_ENUM_TYPE pn, INMOST_DATA_ENUM_TYPE pm)
		: value_tag(t), n(pn), m(pm) {assert(t.GetDataType() == DATA_REAL);}
		static_block_variable(const static_block_variable & other)
		: value_tag(other.value_tag), n(other.n), m(other.m) {}
		static_block_variable & operator =(const static_block_variable & other)
		{
			value_tag = other.value_tag;
			n = other.n;
			m = other.m;
			return * this;
		}
		rMatrix Value(const Storage & e) const {return value_tag(e,n,m);}
		vMatrix Variable(const Storage & e) const {return value_tag(e,n,m);}
		TagRealArray ValueTag() {return value_tag;}
		bool isUnknown(const Storage & e) const {(void)e; return false;}
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new static_block_variable(*this));}
	};
	
	class stored_block_variable : public abstract_dynamic_block_variable
	{
	private:
		TagVariableArray variable_tag;
		INMOST_DATA_ENUM_TYPE n,m;
	public:
		stored_block_variable(TagVariableArray t)
		: variable_tag(t), n(t.GetSize()), m(1) {}
		stored_block_variable(TagVariableArray t, INMOST_DATA_ENUM_TYPE pn, INMOST_DATA_ENUM_TYPE pm)
		: variable_tag(t), n(pn), m(pm) {assert(t.GetDataType() == DATA_VARIABLE);}
		stored_block_variable(const stored_block_variable & other)
		: variable_tag(other.variable_tag), n(other.n), m(other.m) {}
		stored_block_variable & operator =(const stored_block_variable & other)
		{
			variable_tag = other.variable_tag;
			n = other.n;
			m = other.m;
			return * this;
		}
		rMatrix Value(const Storage & e) const { return variable_tag(e,n,m);}
		vMatrix Variable(const Storage & e) const { return variable_tag(e,n,m);}
		TagVariableArray VariableTag() {return variable_tag;}
		bool isUnknown(const Storage & e) const {(void)e; return false;}
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new stored_block_variable(*this));}
	};
	
	
	/// \todo coefficients could be matrices here, introduce another class?
	class stencil_block_variable : public abstract_dynamic_block_variable
	{
	private:
		TagReferenceArray tag_elems;
		TagRealArray      tag_coefs;
		stored_block_variable_expression Arg;
	public:
		stencil_block_variable(Tag tag_elems, Tag tag_coefs, const abstract_dynamic_block_variable & parg)
		: tag_elems(tag_elems), tag_coefs(tag_coefs), Arg(parg) {}
		stencil_block_variable(const stencil_block_variable & other)
		: tag_elems(other.tag_elems), tag_coefs(other.tag_coefs), Arg(other.Arg) {}
		stencil_block_variable & operator =(const stencil_block_variable & other)
		{
			tag_elems = other.tag_elems;
			tag_coefs = other.tag_coefs;
			Arg = other.Arg;
			return * this;
		}
		rMatrix Value(const Storage & e) const
		{
			Storage::real_array      coefs = tag_coefs[e];
			Storage::reference_array elems = tag_elems[e];
			assert(coefs.size() == elems.size());
			rMatrix ret = coefs[0]*Arg.Value(elems[0]);
			for(INMOST_DATA_ENUM_TYPE k = 1; k < elems.size(); ++k)
				ret += coefs[k]*Arg.Value(elems[k]);
			return ret;
		}
		vMatrix Variable(const Storage & e) const
		{
			Storage::real_array      coefs = tag_coefs[e];
			Storage::reference_array elems = tag_elems[e];
			assert(coefs.size() == elems.size());
			vMatrix ret = coefs[0]*Arg.Variable(elems[0]);
			for(INMOST_DATA_ENUM_TYPE k = 1; k < elems.size(); ++k)
				ret += coefs[k]*Arg.Variable(elems[k]);
			return ret;
		}
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new stencil_block_variable(*this));}
	};
	
	///Apply table component-wise on argument matrix.
	class table_block_variable : public abstract_dynamic_block_variable
	{
		stored_block_variable_expression Arg;
		keyval_table Table;
	public:
		table_block_variable(const abstract_dynamic_block_variable & parg, const keyval_table & ptable) : Arg(parg), Table(ptable) {}
		table_block_variable(const table_block_variable & other) : Arg(other.Arg), Table(other.Table) {}
		table_block_variable & operator = (table_block_variable const & other) {Arg = other.Arg; Table = other.Table; return * this;}
		rMatrix Value(const Storage & e) const
		{
			rMatrix ret = Arg.Value(e);
			for(INMOST_DATA_ENUM_TYPE k = 0; k < ret.Rows(); ++k)
				for(INMOST_DATA_ENUM_TYPE l = 0; l < ret.Cols(); ++l)
					ret(k,l) = get_table(ret(k,l),Table);
			return ret;
		}
		vMatrix Variable(const Storage & e) const
		{
			vMatrix ret = Arg.Variable(e);
			for(INMOST_DATA_ENUM_TYPE k = 0; k < ret.Rows(); ++k)
				for(INMOST_DATA_ENUM_TYPE l = 0; l < ret.Cols(); ++l)
					ret(k,l) = get_table(ret(k,l),Table);
			return ret;
		}
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new table_block_variable(*this));}
	};
	
	/// This class makes possible to evaluate different expressions on different element types.
	/// See etype_branch function.
	class etype_branch_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Variable expression to be evaluated when type of provided element matches selected types.
		stored_block_variable_expression ArgB; //< Variable expression to be evaluated when type of provided element does not match selected types.
		ElementType types_true; //< Selected types of elements.
	public:
		/// Constructor. Used by etype_branch function.
		etype_branch_block_variable(ElementType _types_true, const abstract_dynamic_block_variable & _ArgA, const abstract_dynamic_block_variable & _ArgB)
		: types_true(_types_true), ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		etype_branch_block_variable(const etype_branch_block_variable & other)
		: types_true(other.types_true), ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		etype_branch_block_variable & operator =(etype_branch_block_variable const & other)
		{
			types_true = other.types_true;
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return (e->GetElementType() & types_true) ? ArgA.Value(e) : ArgB.Value(e);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return (e->GetElementType() & types_true) ? ArgA.Variable(e) : ArgB.Variable(e);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new etype_branch_block_variable(*this));}
	};
	
	/// This class makes possible to evaluate different expressions depending on the markers.
	/// Works similarly for shared and private markers.
	/// See marker_branch function.
	class marker_branch_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Variable expression to be evaluated when marker is set on the element.
		stored_block_variable_expression ArgB; //< Variable expression to be evaluated when marker is not set on the element.
		MarkerType marker; //< Marker.
	public:
		/// Constructor. Used by marker_branch function.
		marker_branch_block_variable(MarkerType _marker, const abstract_dynamic_block_variable & _ArgA, const abstract_dynamic_block_variable & _ArgB) : marker(_marker), ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		marker_branch_block_variable(const marker_branch_block_variable & other) : marker(other.marker), ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		marker_branch_block_variable & operator =(marker_branch_block_variable const & other)
		{
			marker = other.marker;
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ( isPrivate(marker) ? e->GetPrivateMarker(marker) : e->GetMarker(marker) ) ? ArgA.Value(e) : ArgB.Value(e);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ( isPrivate(marker) ? e->GetPrivateMarker(marker) : e->GetMarker(marker) ) ? ArgA.Variable(e) : ArgB.Variable(e);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new marker_branch_block_variable(*this));}
	};
	
	/// This class represents addition of two matrices.
	class addition_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression on the left.
		stored_block_variable_expression ArgB; //< Block variable expression on the right.
	public:
		/// Constructor.
		addition_block_variable(const abstract_dynamic_block_variable & _ArgA,
								const abstract_dynamic_block_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		addition_block_variable(const addition_block_variable & other)
		: ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		addition_block_variable & operator =(addition_block_variable const & other)
		{
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e) + ArgB.Value(e);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e) + ArgB.Variable(e);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new addition_block_variable(*this));}
	};
	
	/// This class represents subtraction of two matrices.
	class subtraction_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression on the left.
		stored_block_variable_expression ArgB; //< Block variable expression on the right.
	public:
		/// Constructor.
		subtraction_block_variable(const abstract_dynamic_block_variable & _ArgA,
								const abstract_dynamic_block_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		subtraction_block_variable(const subtraction_block_variable & other)
		: ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		subtraction_block_variable & operator =(subtraction_block_variable const & other)
		{
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e) - ArgB.Value(e);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e) - ArgB.Variable(e);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new subtraction_block_variable(*this));}
	};
	
	/// This class represents multiplication of two matrices.
	class multiplication_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression on the left.
		stored_block_variable_expression ArgB; //< Block variable expression on the right.
	public:
		/// Constructor.
		multiplication_block_variable(const abstract_dynamic_block_variable & _ArgA,
								   const abstract_dynamic_block_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		multiplication_block_variable(const multiplication_block_variable & other)
		: ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		multiplication_block_variable & operator =(multiplication_block_variable const & other)
		{
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e)*ArgB.Value(e);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e)*ArgB.Variable(e);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new multiplication_block_variable(*this));}
	};
	
	/// This class represents division of two matrices, this is technically B^{-1}A.
	class division_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression on the left.
		stored_block_variable_expression ArgB; //< Block variable expression on the right.
	public:
		/// Constructor.
		division_block_variable(const abstract_dynamic_block_variable & _ArgA,
								const abstract_dynamic_block_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		division_block_variable(const division_block_variable & other)
		: ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		division_block_variable & operator =(division_block_variable const & other)
		{
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return (ArgA.Value(e)/ArgB.Value(e));}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return (ArgA.Variable(e)/ArgB.Variable(e));}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new division_block_variable(*this));}
	};
	
	/// This class represents division of two matrices, this is technically B^{-1}A.
	class concat_cols_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression on the left.
		stored_block_variable_expression ArgB; //< Block variable expression on the right.
	public:
		/// Constructor.
		concat_cols_block_variable(const abstract_dynamic_block_variable & _ArgA,
								   const abstract_dynamic_block_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		concat_cols_block_variable(const concat_cols_block_variable & other)
		: ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		concat_cols_block_variable & operator =(concat_cols_block_variable const & other)
		{
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e).ConcatCols(ArgB.Value(e));}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e).ConcatCols(ArgB.Variable(e));}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new concat_cols_block_variable(*this));}
	};
	
	/// This class represents division of two matrices, this is technically B^{-1}A.
	class concat_rows_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression on the left.
		stored_block_variable_expression ArgB; //< Block variable expression on the right.
	public:
		/// Constructor.
		concat_rows_block_variable(const abstract_dynamic_block_variable & _ArgA,
								   const abstract_dynamic_block_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		concat_rows_block_variable(const concat_rows_block_variable & other)
		: ArgA(other.ArgA), ArgB(other.ArgB) {}
		/// Assignment operator.
		concat_rows_block_variable & operator =(concat_rows_block_variable const & other)
		{
			ArgA = other.ArgA;
			ArgB = other.ArgB;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e).ConcatRows(ArgB.Value(e));}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e).ConcatRows(ArgB.Variable(e));}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new concat_rows_block_variable(*this));}
	};
	
	
	
	/// This class represents transposition of the matrix, this is A^T.
	class transpose_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression
	public:
		/// Constructor.
		transpose_block_variable(const abstract_dynamic_block_variable & _ArgA)
		: ArgA(_ArgA) {}
		/// Copy constructor.
		transpose_block_variable(const transpose_block_variable & other)
		: ArgA(other.ArgA) {}
		/// Assignment operator.
		transpose_block_variable & operator =(transpose_block_variable const & other)
		{
			ArgA = other.ArgA;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e).Transpose();}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e).Transpose();}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new transpose_block_variable(*this));}
	};
	
	/// This class represents submatrix of the matrix.
	class submatrix_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression
		INMOST_DATA_ENUM_TYPE row1,row2,col1,col2;
	public:
		/// Constructor.
		submatrix_block_variable(const abstract_dynamic_block_variable & _ArgA,
								 INMOST_DATA_ENUM_TYPE row1,
								 INMOST_DATA_ENUM_TYPE row2,
								 INMOST_DATA_ENUM_TYPE col1,
								 INMOST_DATA_ENUM_TYPE col2)
		: ArgA(_ArgA), row1(row1), row2(row2), col1(col1), col2(col2) {}
		/// Copy constructor.
		submatrix_block_variable(const submatrix_block_variable & b)
		: ArgA(b.ArgA), row1(b.row1), row2(b.row2), col1(b.col1), col2(b.col2) {}
		/// Assignment operator.
		submatrix_block_variable & operator =(submatrix_block_variable const & b)
		{
			row1 = b.row1;
			row2 = b.row2;
			col1 = b.col1;
			col2 = b.col2;
			ArgA = b.ArgA;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e)(row1,row2,col1,col2);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e)(row1,row2,col1,col2);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new submatrix_block_variable(*this));}
	};
	
	/// This class represents multiplication of the matrix by the constant.
	class multiply_const_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression
		INMOST_DATA_REAL_TYPE value;
	public:
		/// Constructor.
		multiply_const_block_variable(const abstract_dynamic_block_variable & _ArgA,
								 INMOST_DATA_REAL_TYPE val)
		: ArgA(_ArgA), value(val) {}
		/// Copy constructor.
		multiply_const_block_variable(const multiply_const_block_variable & b)
		: ArgA(b.ArgA), value(b.value) {}
		/// Assignment operator.
		multiply_const_block_variable & operator =(multiply_const_block_variable const & b)
		{
			value = b.value;
			ArgA = b.ArgA;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e)*value;}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e)*value;}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new multiply_const_block_variable(*this));}
	};
	
	/// This class represents multiplication of the matrix by the variable.
	class condition_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression
		stored_block_variable_expression ArgB; //< Block variable expression
		stored_variable_expression ArgC;
	public:
		/// Constructor.
		condition_block_variable(const abstract_dynamic_variable & _ArgC,
								 const abstract_dynamic_block_variable & _ArgA,
								 const abstract_dynamic_block_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB), ArgC(_ArgC) {}
		/// Copy constructor.
		condition_block_variable(const condition_block_variable & b)
		: ArgA(b.ArgA), ArgB(b.ArgB), ArgC(b.ArgC) {}
		/// Assignment operator.
		condition_block_variable & operator =(condition_block_variable const & b)
		{
			ArgB = b.ArgB;
			ArgA = b.ArgA;
			ArgC = b.ArgC;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgC.Value(e) > 0.0 ? ArgA.Value(e) : ArgB.Value(e);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgC.Value(e) > 0.0 ? ArgA.Variable(e) : ArgB.Variable(e);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new condition_block_variable(*this));}
	};
	
	/// This class represents multiplication of the matrix by the variable.
	class multiply_variable_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression
		stored_variable_expression ArgB;
	public:
		/// Constructor.
		multiply_variable_block_variable(const abstract_dynamic_block_variable & _ArgA,
										 const abstract_dynamic_variable & _ArgB)
		: ArgA(_ArgA), ArgB(_ArgB) {}
		/// Copy constructor.
		multiply_variable_block_variable(const multiply_variable_block_variable & b)
		: ArgA(b.ArgA), ArgB(b.ArgB) {}
		/// Assignment operator.
		multiply_variable_block_variable & operator =(multiply_variable_block_variable const & b)
		{
			ArgB = b.ArgB;
			ArgA = b.ArgA;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e)*ArgB.Value(e);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e)*ArgB.Variable(e);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new multiply_variable_block_variable(*this));}
	};
	
	/// This class represents inverse of the matrix, this is A^{-1}.
	class inverse_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression
	public:
		/// Constructor.
		inverse_block_variable(const abstract_dynamic_block_variable & _ArgA)
		: ArgA(_ArgA) {}
		/// Copy constructor.
		inverse_block_variable(const inverse_block_variable & other)
		: ArgA(other.ArgA) {}
		/// Assignment operator.
		inverse_block_variable & operator =(inverse_block_variable const & other)
		{
			ArgA = other.ArgA;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e).Invert();}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e).Invert();}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new inverse_block_variable(*this));}
	};
	
	/// This class represents pseudo-inverse of the matrix, this is A^{+}.
	class pseudo_inverse_block_variable : public abstract_dynamic_block_variable
	{
	private:
		stored_block_variable_expression ArgA; //< Block variable expression
		double eps;
	public:
		/// Constructor.
		pseudo_inverse_block_variable(const abstract_dynamic_block_variable & _ArgA, double _eps = 1.0e-13)
		: ArgA(_ArgA), eps(_eps) {}
		/// Copy constructor.
		pseudo_inverse_block_variable(const pseudo_inverse_block_variable & other)
		: ArgA(other.ArgA), eps(other.eps) {}
		/// Assignment operator.
		pseudo_inverse_block_variable & operator =(pseudo_inverse_block_variable const & other)
		{
			ArgA = other.ArgA;
			eps = other.eps;
			return *this;
		}
		/// Get value of variable expression on provided element e.
		rMatrix Value(const Storage & e) const
		{return ArgA.Value(e).PseudoInvert(eps);}
		/// Get value with derivatives of variable expression on provided element e.
		/// This function collapses associated expression tree into multivar_expression.
		vMatrix Variable(const Storage & e) const
		{return ArgA.Variable(e).PseudoInvert(eps);}
		/// Make a copy of this class, used to reproduce and store a tree of variable expressions.
		abstract_dynamic_block_variable * Copy() const {return static_cast<abstract_dynamic_block_variable *>(new pseudo_inverse_block_variable(*this));}
	};
	
	
	
	
	typedef abstract_dynamic_block_variable abstract_block_variable;
}

/// If control evaluates to non-negative number, returns A, otherwise B
__INLINE 
INMOST::condition_block_variable
condition(INMOST::abstract_dynamic_variable const & control,
		  INMOST::abstract_dynamic_block_variable const & if_ge_zero,
		  INMOST::abstract_dynamic_block_variable const & if_lt_zero)
{ return INMOST::condition_block_variable(control,if_ge_zero,if_lt_zero); }
/// Matrix operation A+B
__INLINE
INMOST::addition_block_variable
operator+(INMOST::abstract_block_variable const & Left,
		  INMOST::abstract_block_variable const & Right)
{ return INMOST::addition_block_variable(Left, Right); }
/// Matrix operation A-B
__INLINE
INMOST::subtraction_block_variable
operator-(INMOST::abstract_block_variable const & Left,
		  INMOST::abstract_block_variable const & Right)
{ return INMOST::subtraction_block_variable(Left, Right); }
/// Matrix operation A*B
__INLINE
INMOST::multiplication_block_variable
operator*(INMOST::abstract_block_variable const & Left,
		  INMOST::abstract_block_variable const & Right)
{ return INMOST::multiplication_block_variable(Left, Right); }
/// Matrix operation B^{-1}*A
__INLINE
INMOST::division_block_variable
operator/(INMOST::abstract_block_variable const & Left,
		  INMOST::abstract_block_variable const & Right)
{ return INMOST::division_block_variable(Left, Right); }
/// Attach two matrices by rows, i.e. C = [A,B]
__INLINE
INMOST::concat_rows_block_variable
concat_rows(INMOST::abstract_block_variable const & Left,
			INMOST::abstract_block_variable const & Right)
{ return INMOST::concat_rows_block_variable(Left, Right); }
/// Attach two matrices by columns, i.e. C = [A^T,B^T]^T
__INLINE
INMOST::concat_cols_block_variable
concat_cols(INMOST::abstract_block_variable const & Left,
			INMOST::abstract_block_variable const & Right)
{ return INMOST::concat_cols_block_variable(Left, Right); }
/// Matrix operation A^T
__INLINE
INMOST::transpose_block_variable
transpose(INMOST::abstract_block_variable const & Left)
{ return INMOST::transpose_block_variable(Left); }
/// Matrix operation A^{-1}
__INLINE
INMOST::inverse_block_variable
inv(INMOST::abstract_block_variable const & Left)
{ return INMOST::inverse_block_variable(Left); }
/// Matrix operation A^{+}
__INLINE
INMOST::pseudo_inverse_block_variable
pinv(INMOST::abstract_block_variable const & Left)
{ return INMOST::pseudo_inverse_block_variable(Left); }
/// Submatrix of a matrix
__INLINE
INMOST::submatrix_block_variable
submatrix(INMOST::abstract_block_variable const & Left,
	 INMOST_DATA_ENUM_TYPE row1,
	 INMOST_DATA_ENUM_TYPE row2,
	 INMOST_DATA_ENUM_TYPE col1,
	 INMOST_DATA_ENUM_TYPE col2)
{ return INMOST::submatrix_block_variable(Left,row1,row2,col1,col2); }
/// Matrix multiplication by a constant A*a
__INLINE
INMOST::multiply_const_block_variable
operator *(INMOST::abstract_block_variable const & Left,
		   INMOST_DATA_REAL_TYPE Right)
{ return INMOST::multiply_const_block_variable(Left,Right); }
/// Matrix multiplication by a constant a*A
__INLINE
INMOST::multiply_const_block_variable
operator *(INMOST_DATA_REAL_TYPE Left,
		   INMOST::abstract_block_variable const & Right)
{ return INMOST::multiply_const_block_variable(Right,Left); }
/// Matrix multiplication by inverse of a constant A/a
__INLINE
INMOST::multiply_const_block_variable
operator /(INMOST::abstract_block_variable const & Left,
		   INMOST_DATA_REAL_TYPE Right)
{ return INMOST::multiply_const_block_variable(Left,1.0/Right); }
/// Matrix multiplication by a scalar expression A*a
__INLINE
INMOST::multiply_variable_block_variable
operator *(INMOST::abstract_block_variable const & Left,
		   INMOST::abstract_variable const & Right)
{ return INMOST::multiply_variable_block_variable(Left,Right); }
/// Matrix multiplication by a scalar expression a*A
__INLINE
INMOST::multiply_variable_block_variable
operator *(INMOST::abstract_variable const & Left,
		   INMOST::abstract_block_variable const & Right)
{ return INMOST::multiply_variable_block_variable(Right,Left); }
/// Matrix multiplication by inverse of a scalar expression A/a
__INLINE
INMOST::multiply_variable_block_variable
operator /(INMOST::abstract_block_variable const & Left,
		   INMOST::abstract_variable const & Right)
{ return INMOST::multiply_variable_block_variable(Left,1.0/INMOST::stored_variable_expression(Right)); }
/// Calculation of matrix convex combination on stencil
__INLINE
INMOST::stencil_block_variable
stencil(INMOST::Tag tag_elems,
		INMOST::Tag tag_coefs,
		INMOST::abstract_dynamic_block_variable const & Arg)
{ return INMOST::stencil_block_variable(tag_elems,tag_coefs,Arg); }
/// Operation for key-value table
__INLINE
INMOST::table_block_variable
get_table(INMOST::abstract_dynamic_block_variable const & Arg,
		  INMOST::keyval_table const & Table)
{return INMOST::table_block_variable(Arg,Table);}
/// Branching expression by element type
__INLINE
INMOST::etype_branch_block_variable
etype_branch(INMOST::ElementType true_type,
			 INMOST::abstract_dynamic_block_variable const & iftrue,
			 INMOST::abstract_dynamic_block_variable const & iffalse)
{return INMOST::etype_branch_block_variable(true_type,iftrue,iffalse);}
/// Branching expression by marker
__INLINE
INMOST::marker_branch_block_variable
marker_branch(INMOST::MarkerType marker,
			  INMOST::abstract_dynamic_block_variable const & iftrue,
			  INMOST::abstract_dynamic_block_variable const & iffalse)
{return INMOST::marker_branch_block_variable(marker,iftrue,iffalse);}

#endif //defined(USE_AUTODIFF) && defined(USE_MESH)




#endif //INMOST_AUTODIFF_ETBVAR_H_INCLUDED

