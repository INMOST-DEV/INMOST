#ifndef INMOST_AUTODIFF_ETEXPR_H_INCLUDED
#define INMOST_AUTODIFF_ETEXPR_H_INCLUDED
#include "inmost_common.h"
#include "inmost_sparse.h"
#include <sstream> //for debug
#include <new>

//TODO:
// 1. Incorporate tables
// 2. Consider optimization by checking zero variation multipliers, check that assembly do not degrade.
// 3. floor, ceil, atan, acos, asin, max, min functions
// 4. choice of directional derivatives at discontinuities for abs, pow, max, min (see ADOL-C)


#ifdef _MSC_VER
#pragma warning(disable : 4503)
#endif


#if defined(USE_AUTODIFF)
namespace INMOST
{
	

	class basic_expression
	{
	public:
		typedef Sparse::RowMerger merger_type;
		static merger_type& GetMerger() { return *merger; }
	protected:
		static thread_private<merger_type> merger;
	public:
		basic_expression() {}
		virtual INMOST_DATA_REAL_TYPE GetValue() const = 0;
		virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const = 0;
		virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const = 0;
		virtual void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const = 0;
		virtual ~basic_expression() {}
		void GetDerivatives(INMOST_DATA_REAL_TYPE mult, Sparse::Row& r) const 
		{ 
			Sparse::RowMerger& m = GetMerger();
			m.Clear();
			m.AddRow(1.0, r);
			GetJacobian(mult, m);
			m.RetrieveRow(r);
		}
		__INLINE static void UseMerger(INMOST_DATA_REAL_TYPE coefa, const Sparse::Row& rowa, INMOST_DATA_REAL_TYPE coefb, Sparse::Row& rowb)
		{
			Sparse::RowMerger& m = GetMerger();
			m.Clear();
			m.AddRow(coefb, rowb);
			m.AddRow(coefa, rowa);
			m.RetrieveRow(rowb);
		}
		__INLINE static void UseMerger(INMOST_DATA_REAL_TYPE coefa, const basic_expression& expra, INMOST_DATA_REAL_TYPE coefb, Sparse::Row& rowb)
		{
			Sparse::RowMerger& m = GetMerger();
			m.Clear();
			m.AddRow(coefb, rowb);
			expra.GetJacobian(coefa, m);
			m.RetrieveRow(rowb);
		}
	};

	template<class Derived>
	class shell_expression : virtual public basic_expression
	{
	public:
		shell_expression() {}
		__INLINE virtual INMOST_DATA_REAL_TYPE GetValue() const {return static_cast<const Derived *>(this)->GetValue(); }
		__INLINE virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const { if( mult ) return static_cast<const Derived *>(this)->GetJacobian(mult,r); }
		__INLINE virtual void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const { if (mult) return static_cast<const Derived*>(this)->GetJacobian(mult, r); }
		__INLINE virtual void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {return static_cast<const Derived *>(this)->GetHessian(multJ,J,multH,H); }
		operator Derived & () {return *static_cast<Derived *>(this);}
		operator const Derived & () const {return *static_cast<const Derived *>(this);}
		~shell_expression() {}
	};
	
	
	class const_expression : public shell_expression<const_expression>
	{
		INMOST_DATA_REAL_TYPE value;
	public:
		const_expression(const const_expression & other) :value(other.value) {}
		const_expression(INMOST_DATA_REAL_TYPE pvalue) : value(pvalue) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {(void)mult; (void)r;}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const { (void)mult; (void)r; }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {(void)multJ; (void)J; (void)multH; (void)H;}
		__INLINE const_expression & operator =(const_expression const & other)
		{
			value = other.value;
			return *this;
		}
	};
	
	
	class var_expression : public shell_expression<var_expression>
	{
		INMOST_DATA_REAL_TYPE value;
		INMOST_DATA_ENUM_TYPE index;
	public:
		var_expression() : value(0), index(ENUMUNDEF) {}
		var_expression(const var_expression & other) :value(other.value), index(other.index) {}
		var_expression(INMOST_DATA_REAL_TYPE pvalue, INMOST_DATA_ENUM_TYPE pindex) : value(pvalue), index(pindex) {}
		var_expression(INMOST_DATA_REAL_TYPE pvalue) : value(pvalue), index(ENUMUNDEF) {}
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE INMOST_DATA_ENUM_TYPE GetIndex() const { return index; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const {if( mult && index != ENUMUNDEF ) r.AddValue(mult,index);}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const { if (mult && index != ENUMUNDEF) r[index] += mult; }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const {if( index != ENUMUNDEF ) J.Push(index,multJ);  (void)multH; (void)H;}
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(INMOST_DATA_ENUM_TYPE i) const {return index == i? 1.0 : 0.0;}
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
		bool check_infs() const
		{
			return __isinf__(value);
		}
	};
	
#if defined(PACK_ARRAY)
#pragma pack(push,r1,4)
#endif
	/// A class that represents a variable with multiple
	/// first order variations.
	/// Short type name is variable.
	class multivar_expression : public shell_expression<multivar_expression>
	{
		INMOST_DATA_REAL_TYPE value; //< Value of the variable.
		Sparse::Row entries; //< Sparse vector of variations.
	public:
		multivar_expression() :value(0) {}
		multivar_expression(INMOST_DATA_REAL_TYPE pvalue) : value(pvalue) {}
		multivar_expression(const multivar_expression & other) : value(other.value), entries(other.entries) {}
		multivar_expression(multivar_expression&& other) : value(other.value), entries(std::move(other.entries)) {}
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
			UseMerger(1.0, expr, 0.0, entries);
		}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			r.AddRow(mult, entries);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			for (Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
				r[it->first] += it->second * mult;
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			J = entries;
			if( !J.isSorted() ) std::sort(J.Begin(),J.End());
			for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second *= multJ;
			H.Clear();
            (void)multH;
		}
		__INLINE multivar_expression & operator = (INMOST_DATA_REAL_TYPE pvalue)
		{
			value = pvalue;
			entries.Clear();
			return *this;
		}
		__INLINE multivar_expression& operator = (multivar_expression && other)
		{
			value = other.value;
			entries = std::move(other.entries);
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
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(INMOST_DATA_ENUM_TYPE index) const {return GetRow().get_safe(index);}
		__INLINE bool operator < (INMOST_DATA_REAL_TYPE right) const {return value < right;}
		__INLINE bool operator > (INMOST_DATA_REAL_TYPE right) const {return value > right;}
		__INLINE bool operator <=(INMOST_DATA_REAL_TYPE right) const {return value <= right;}
		__INLINE bool operator >=(INMOST_DATA_REAL_TYPE right) const {return value >= right;}
		__INLINE bool operator ==(INMOST_DATA_REAL_TYPE right) const {return value == right;}
		__INLINE bool operator !=(INMOST_DATA_REAL_TYPE right) const {return value != right;}
		__INLINE bool operator < (basic_expression const & expr) const {return value < expr.GetValue();}
		__INLINE bool operator > (basic_expression const & expr) const {return value > expr.GetValue();}
		__INLINE bool operator <=(basic_expression const & expr) const {return value <= expr.GetValue();}
		__INLINE bool operator >=(basic_expression const & expr) const {return value >= expr.GetValue();}
		__INLINE bool operator ==(basic_expression const & expr) const {return value == expr.GetValue();}
		__INLINE bool operator !=(basic_expression const & expr) const {return value != expr.GetValue();}
		__INLINE bool operator < (multivar_expression const & expr) const {return value < expr.GetValue();}
		__INLINE bool operator > (multivar_expression const & expr) const {return value > expr.GetValue();}
		__INLINE bool operator <=(multivar_expression const & expr) const {return value <= expr.GetValue();}
		__INLINE bool operator >=(multivar_expression const & expr) const {return value >= expr.GetValue();}
		__INLINE bool operator ==(multivar_expression const & expr) const {return value == expr.GetValue();}
		__INLINE bool operator !=(multivar_expression const & expr) const {return value != expr.GetValue();}
		__INLINE multivar_expression & operator +=(basic_expression const & expr)
		{
			value += expr.GetValue();
			UseMerger(1.0, expr, 1.0, entries);
			return *this;
		}
		__INLINE multivar_expression & operator -=(basic_expression const & expr)
		{
			value -= expr.GetValue();
			UseMerger(-1.0, expr, 1.0, entries);
			return *this;
		}
		__INLINE multivar_expression & operator *=(basic_expression const & expr)
		{
			INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
			UseMerger(lval, expr, rval, entries);
			value *= rval;
			return *this;
		}
		__INLINE multivar_expression & operator /=(basic_expression const & expr)
		{
			INMOST_DATA_REAL_TYPE rval = expr.GetValue();
			INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0/rval;
			value *= reciprocial_rval;
			UseMerger(-value * reciprocial_rval, expr, reciprocial_rval, entries);
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
		bool check_infs() const
		{
			if( __isinf__(value) ) return true;
			for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
				if( __isinf__(it->second) ) return true;
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
		/// Retrieve variable from array of entries.
		/// Size of array without retrieval can be determined via RetrieveSize.
		/// @param v Array of entries that will store data of the variable.
		/// @return Number of entries red.
		INMOST_DATA_ENUM_TYPE Retrieve(const Sparse::Row::entry * v)
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
		static INMOST_DATA_ENUM_TYPE RetrieveSize(const Sparse::Row::entry * v)
		{
			return 1 + v[0].first;
		}
		void Print(double eps = -1, std::ostream & sout = std::cout) const
		{
			sout << value << " ";
			entries.Print(eps, sout);
			if (entries.Empty()) std::cout << std::endl;
		}
		void swap(multivar_expression & b)
		{
			std::swap(value,b.value);
			entries.Swap(b.entries);
		}
		friend class multivar_expression_reference;
	};

	class multivar_dense_expression : public shell_expression<multivar_dense_expression>
	{
		INMOST_DATA_REAL_TYPE value; //< Value of the variable.
		INMOST_DATA_REAL_TYPE * entries; //< Dense vector of variations.
		INMOST_DATA_ENUM_TYPE nentries; //
	public:
		multivar_dense_expression() : value(0), entries(NULL), nentries(0) {}
		multivar_dense_expression(INMOST_DATA_REAL_TYPE _value, INMOST_DATA_REAL_TYPE * _entries, INMOST_DATA_ENUM_TYPE _nentries) 
			: value(_value), entries(_entries), nentries(_nentries) {}
		multivar_dense_expression(const multivar_dense_expression& other) 
			: value(other.value), entries(other.entries), nentries(other.nentries) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const
		{
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k)
				r.AddValue(entries[k] * mult, k);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const
		{
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k)
				r[k] += entries[k] * mult;
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const
		{
			throw NotImplemented;
		}
		__INLINE multivar_dense_expression& operator = (INMOST_DATA_REAL_TYPE pvalue)
		{
			value = pvalue;
			memset(entries, 0, sizeof(INMOST_DATA_REAL_TYPE) * nentries);
			return *this;
		}
		__INLINE multivar_dense_expression& operator = (multivar_dense_expression const& other)
		{
			value = other.value;
			memcpy(entries, other.entries, sizeof(INMOST_DATA_REAL_TYPE) * nentries);
			return *this;
		}
		__INLINE multivar_dense_expression& operator = (basic_expression const& expr)
		{
			value = expr.GetValue();
			memset(entries, 0, sizeof(INMOST_DATA_REAL_TYPE) * nentries);
			expr.GetJacobian(1.0, entries);
			return *this;
		}
		
		__INLINE INMOST_DATA_REAL_TYPE * GetRow() { return entries; }
		__INLINE const INMOST_DATA_REAL_TYPE * GetRow() const { return entries; }
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(INMOST_DATA_ENUM_TYPE index) const { return GetRow()[index]; }
		__INLINE bool operator < (INMOST_DATA_REAL_TYPE right) const { return value < right; }
		__INLINE bool operator > (INMOST_DATA_REAL_TYPE right) const { return value > right; }
		__INLINE bool operator <=(INMOST_DATA_REAL_TYPE right) const { return value <= right; }
		__INLINE bool operator >=(INMOST_DATA_REAL_TYPE right) const { return value >= right; }
		__INLINE bool operator ==(INMOST_DATA_REAL_TYPE right) const { return value == right; }
		__INLINE bool operator !=(INMOST_DATA_REAL_TYPE right) const { return value != right; }
		__INLINE bool operator < (basic_expression const& expr) const { return value < expr.GetValue(); }
		__INLINE bool operator > (basic_expression const& expr) const { return value > expr.GetValue(); }
		__INLINE bool operator <=(basic_expression const& expr) const { return value <= expr.GetValue(); }
		__INLINE bool operator >=(basic_expression const& expr) const { return value >= expr.GetValue(); }
		__INLINE bool operator ==(basic_expression const& expr) const { return value == expr.GetValue(); }
		__INLINE bool operator !=(basic_expression const& expr) const { return value != expr.GetValue(); }
		__INLINE bool operator < (multivar_dense_expression const& expr) const { return value < expr.GetValue(); }
		__INLINE bool operator > (multivar_dense_expression const& expr) const { return value > expr.GetValue(); }
		__INLINE bool operator <=(multivar_dense_expression const& expr) const { return value <= expr.GetValue(); }
		__INLINE bool operator >=(multivar_dense_expression const& expr) const { return value >= expr.GetValue(); }
		__INLINE bool operator ==(multivar_dense_expression const& expr) const { return value == expr.GetValue(); }
		__INLINE bool operator !=(multivar_dense_expression const& expr) const { return value != expr.GetValue(); }
		__INLINE multivar_dense_expression& operator +=(basic_expression const& expr)
		{
			value += expr.GetValue();
			expr.GetJacobian(1.0, entries);
			return *this;
		}
		__INLINE multivar_dense_expression& operator -=(basic_expression const& expr)
		{
			value -= expr.GetValue();
			expr.GetJacobian(-1.0, entries);
			return *this;
		}
		__INLINE multivar_dense_expression& operator *=(basic_expression const& expr)
		{
			INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k) 
					entries[k] *= rval;
			expr.GetJacobian(lval, entries);
			value *= rval;
			return *this;
		}
		__INLINE multivar_dense_expression& operator /=(basic_expression const& expr)
		{
			INMOST_DATA_REAL_TYPE rval = expr.GetValue();
			INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0 / rval;
			value *= reciprocial_rval;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k) 
				entries[k] *= reciprocial_rval;
			expr.GetJacobian(-value * reciprocial_rval, entries);
			return *this;
		}
		__INLINE multivar_dense_expression& operator +=(INMOST_DATA_REAL_TYPE right)
		{
			value += right;
			return *this;
		}
		__INLINE multivar_dense_expression& operator -=(INMOST_DATA_REAL_TYPE right)
		{
			value -= right;
			return *this;
		}
		__INLINE multivar_dense_expression& operator *=(INMOST_DATA_REAL_TYPE right)
		{
			value *= right;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k) 
				entries[k] *= right;
			return *this;
		}
		__INLINE multivar_dense_expression& operator /=(INMOST_DATA_REAL_TYPE right)
		{
			value /= right;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k) 
				entries[k] /= right;
			return *this;
		}
		bool check_nans() const
		{
			if (value != value) return true;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k)
				if (entries[k] != entries[k]) return true;
			return false;
		}
		bool check_infs() const
		{
			if (__isinf__(value)) return true;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nentries; ++k)
				if (__isinf__(entries[k])) return true;
			return false;
		}
		void swap(multivar_dense_expression& b)
		{
			std::swap(value, b.value);
			std::swap(entries, b.entries);
			std::swap(nentries, b.nentries);
		}
	};
	
	/// A class that represents a variable with multiple
	/// first order and second order variations.
	/// Short type name is hessian_variable.
	class hessian_multivar_expression : public shell_expression<hessian_multivar_expression>
	{
		INMOST_DATA_REAL_TYPE value;
		Sparse::Row entries;
		Sparse::HessianRow hessian_entries;
	public:
		void swap(hessian_multivar_expression & b)
		{
			std::swap(value,b.value);
			entries.Swap(b.entries);
			hessian_entries.Swap(b.hessian_entries);
		}
		/// Sets zero value and no first or second order variations.
		hessian_multivar_expression() :value(0) {}
		/// Sets value and no first or second order variations.
		hessian_multivar_expression(INMOST_DATA_REAL_TYPE pvalue) : value(pvalue) {}
		/// Copy value and all first and second order variations from hessian_variable.
		hessian_multivar_expression(const hessian_multivar_expression & other) : value(other.value), entries(other.entries), hessian_entries(other.hessian_entries) {}
		/// Copy value and all first variations from variable, no second order variations.
		//hessian_multivar_expression(const multivar_expression & other) : value(other.GetValue()), entries(other.GetRow()) {}
		/// Sets value and all first variations, no second order variations.
		hessian_multivar_expression(INMOST_DATA_REAL_TYPE pvalue, Sparse::Row & pentries)
		: value(pvalue), entries(pentries) {}
		/// Sets value and all first and second order variations.
		hessian_multivar_expression(INMOST_DATA_REAL_TYPE pvalue, Sparse::Row & pentries, Sparse::HessianRow & phentries)
		: value(pvalue), entries(pentries), hessian_entries(phentries) {}
		/// Sets value and it's variation index.
		hessian_multivar_expression(INMOST_DATA_REAL_TYPE pvalue, INMOST_DATA_ENUM_TYPE pindex, INMOST_DATA_REAL_TYPE pdmult = 1.0)
		: value(pvalue)
		{
			entries.Push(pindex,pdmult);
		}
		/// Evaluates argument expression to get the value and all variations of variable.
		hessian_multivar_expression(const basic_expression & expr)
		{
			value = expr.GetValue();
			expr.GetHessian(1.0,entries,1.0,hessian_entries);
		}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			r.AddRow(mult, entries);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			for (Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
				r[it->first] += it->second * mult;
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			J = entries;
			if( !J.isSorted() ) std::sort(J.Begin(),J.End());
			for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second *= multJ;
			H = hessian_entries;
			for(Sparse::HessianRow::iterator it = H.Begin(); it != H.End(); ++it) it->second *= multH;
		}
		__INLINE multivar_expression GetVariable(INMOST_DATA_ENUM_TYPE index)
		{
			multivar_expression ret(0);
			for(int k = 0; k < (int)entries.Size(); ++k)
			{
				if( entries.GetIndex(k) == index )
				{
					ret.SetValue(entries.GetValue(k));
					break;
				}
				
			}
			Sparse::Row& r = ret.GetRow();
			for (int q = 0; q < (int)hessian_entries.Size(); ++q)
			{
				Sparse::HessianRow::index& i = hessian_entries.GetIndex(q);
				if (i.first == index)
					r.Push(i.second, hessian_entries.GetValue(q) * (i.first == i.second ? 1.0 : 0.5));
				else if (i.second == index)
					r.Push(i.first, hessian_entries.GetValue(q) * (i.first == i.second ? 1.0 : 0.5));
			}
			return ret;
		}
		__INLINE hessian_multivar_expression & operator = (INMOST_DATA_REAL_TYPE pvalue)
		{
			value = pvalue;
			entries.Clear();
			hessian_entries.Clear();
			return *this;
		}
		/*
		__INLINE hessian_multivar_expression & operator = (basic_expression const & expr)
		{
			value = expr.GetValue();
			Sparse::Row tmp;
			Sparse::HessianRow htmp;
			expr.GetHessian(1.0,tmp,1.0,htmp);
			entries.Swap(tmp);
			hessian_entries.Swap(htmp);
			return *this;
		}
		*/
		/*
		__INLINE hessian_multivar_expression & operator = (multivar_expression const & other)
		{
			value = other.GetValue();
			entries = other.GetRow();
			hessian_entries.Clear();
			return *this;
		}
		*/
		__INLINE hessian_multivar_expression & operator = (hessian_multivar_expression const & other)
		{
			value = other.value;
			entries = other.entries;
			hessian_entries = other.hessian_entries;
			return *this;
		}
		__INLINE Sparse::Row & GetRow() {return entries;}
		__INLINE Sparse::HessianRow & GetHessianRow() {return hessian_entries;}
		__INLINE const Sparse::Row & GetRow() const {return entries;}
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(INMOST_DATA_ENUM_TYPE index) const {return GetRow().get_safe(index);}
		__INLINE const Sparse::HessianRow & GetHessianRow() const {return hessian_entries;}
		__INLINE hessian_multivar_expression & operator +=(basic_expression const & expr)
		{
			value += expr.GetValue();
			Sparse::Row tmpr, tmp;
			Sparse::HessianRow htmpr, htmp;
			expr.GetHessian(1.0,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(1.0,entries,1.0,tmpr,tmp);
			entries.Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(1.0,hessian_entries,1.0,htmpr,htmp);
			hessian_entries.Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression & operator -=(basic_expression const & expr)
		{
			value -= expr.GetValue();
			Sparse::Row tmpr, tmp;
			Sparse::HessianRow htmpr, htmp;
			expr.GetHessian(1.0,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(1.0,entries,-1.0,tmpr,tmp);
			entries.Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(1.0,hessian_entries,-1.0,htmpr,htmp);
			hessian_entries.Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression & operator *=(basic_expression const & expr)
		{
			throw NotImplemented; //check code below is correct
			INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
			Sparse::Row tmp, tmpr;
			Sparse::HessianRow htmp, htmpr;
			expr.GetHessian(1.0,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(rval,entries,lval,tmpr,tmp);
			entries.Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(lval,hessian_entries,rval,htmpr,htmp);
			hessian_entries.Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression & operator /=(basic_expression const & expr)
		{
			throw NotImplemented; //check code below is correct
			INMOST_DATA_REAL_TYPE rval = expr.GetValue();
			INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0/rval;
			value *= reciprocial_rval;
			Sparse::Row tmp, tmpr;
			Sparse::HessianRow htmp, htmpr;
			expr.GetHessian(-value*reciprocial_rval,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(reciprocial_rval,entries,-value*reciprocial_rval,tmpr,tmp);
			entries.Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(reciprocial_rval,hessian_entries,-value*reciprocial_rval,htmpr,htmp);
			hessian_entries.Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression & operator +=(INMOST_DATA_REAL_TYPE right)
		{
			value += right;
			return *this;
		}
		__INLINE hessian_multivar_expression & operator -=(INMOST_DATA_REAL_TYPE right)
		{
			value -= right;
			return *this;
		}
		__INLINE hessian_multivar_expression & operator *=(INMOST_DATA_REAL_TYPE right)
		{
			value *= right;
			for(Sparse::Row::iterator it = entries.Begin(); it != entries.End(); ++it) it->second *= right;
			for(Sparse::HessianRow::iterator it = hessian_entries.Begin(); it != hessian_entries.End(); ++it) it->second *= right;
			return *this;
		}
		__INLINE hessian_multivar_expression & operator /=(INMOST_DATA_REAL_TYPE right)
		{
			value /= right;
			for(Sparse::Row::iterator it = entries.Begin(); it != entries.End(); ++it) it->second /= right;
			for(Sparse::HessianRow::iterator it = hessian_entries.Begin(); it != hessian_entries.End(); ++it) it->second /= right;
			return *this;
		}
		bool check_nans() const
		{
			if( value != value ) return true;
			for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
				if( it->second != it->second ) return true;
			for(Sparse::HessianRow::const_iterator it = hessian_entries.Begin(); it != hessian_entries.End(); ++it)
				if( it->second != it->second ) return true;
			return false;
		}
		bool check_infs() const
		{
			if( __isinf__(value)) return true;
			for(Sparse::Row::const_iterator it = entries.Begin(); it != entries.End(); ++it)
				if( __isinf__(it->second) ) return true;
			for(Sparse::HessianRow::const_iterator it = hessian_entries.Begin(); it != hessian_entries.End(); ++it)
				if( __isinf__(it->second) ) return true;
			return false;
		}
		friend class hessian_multivar_expression_reference;
	};
#if defined(PACK_ARRAY)
#pragma pack(pop,r1)
#endif
	
	static INMOST_DATA_REAL_TYPE stub_multivar_expression_reference_value; //for default constructor in multivar_expression_reference

	
	class multivar_expression_reference : public shell_expression<multivar_expression_reference>
	{
		INMOST_DATA_REAL_TYPE & value;
		Sparse::Row * entries;
	public:
		/// Default constructor
		multivar_expression_reference() : value(stub_multivar_expression_reference_value), entries(NULL) {}
		/// Constructor, set links to the provided value and entries
		multivar_expression_reference(INMOST_DATA_REAL_TYPE & _value, Sparse::Row * _entries)
		: value(_value), entries(_entries) {}
		/// Copy constructor, sets links to the same reference of value and entries
		multivar_expression_reference(const multivar_expression_reference & other)
		: value(other.value), entries(other.entries) {}
		/// Copy constructor from multivar_expression, sets links to the same reference of value and entries
		multivar_expression_reference(multivar_expression & other)
		: value(other.value), entries(&other.entries) {}
		/// Retrieve value
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		/// Set value without changing derivatives
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		/// Retrieve derivatives with multiplier into Sparse::RowMerger structure.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			if( entries )
				r.AddRow(mult,*entries);
		}
		/// Retrieve derivatives with multiplier into array.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			if (entries)
			{
				for (Sparse::Row::const_iterator it = entries->Begin(); it != entries->End(); ++it)
					r[it->first] += it->second * mult;
			}
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J,INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			if( entries )
			{
				J = *entries;
				if( !J.isSorted() ) std::sort(J.Begin(),J.End());
				for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second *= multJ;
				H.Clear();
			}
			(void)multH;
		}
		multivar_expression_reference & operator = (INMOST_DATA_REAL_TYPE pvalue)
		{
			value = pvalue;
			entries->Clear();
			return *this;
		}
		__INLINE multivar_expression_reference & operator = (basic_expression const & expr)
		{
			value = expr.GetValue();
			if (entries) UseMerger(1.0, expr, 0.0, *entries);
			return *this;
		}
		__INLINE multivar_expression_reference & operator = (multivar_expression_reference const & other)
		{
			value = other.GetValue();
			*entries = other.GetRow();
			return *this;
		}
		__INLINE Sparse::Row & GetRow() {return *entries;}
		__INLINE const Sparse::Row & GetRow() const {return *entries;}
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(INMOST_DATA_ENUM_TYPE index) const {return GetRow().get_safe(index);}
		__INLINE multivar_expression_reference & operator +=(basic_expression const & expr)
		{
			value += expr.GetValue();
			if (entries) UseMerger(1.0, expr, 1.0, *entries);
			return *this;
		}
		__INLINE multivar_expression_reference & operator -=(basic_expression const & expr)
		{
			value -= expr.GetValue();
			if (entries) UseMerger(-1.0, expr, 1.0, *entries);
			return *this;
		}
		__INLINE multivar_expression_reference & operator *=(basic_expression const & expr)
		{
			INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
			if (entries) UseMerger(lval, expr, rval, *entries);
			value *= rval;
			return *this;
		}
		__INLINE multivar_expression_reference & operator /=(basic_expression const & expr)
		{
			INMOST_DATA_REAL_TYPE rval = expr.GetValue();
			INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0/rval;
			value *= reciprocial_rval;
			if (entries) UseMerger(-value * reciprocial_rval, expr, reciprocial_rval, *entries);
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
		bool check_infs() const
		{
			if( __isinf__(value) ) return true;
			for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it)
				if( __isinf__(it->second) ) return true;
			return false;
		}
	};

	class value_reference 
	{
		INMOST_DATA_REAL_TYPE& value;
	public:
		/// Default constructor
		value_reference() : value(stub_multivar_expression_reference_value) {}
		/// Constructor, set links to the provided value and entries
		value_reference(INMOST_DATA_REAL_TYPE& _value)
			: value(_value) {}
		/// Copy constructor, sets links to the same reference of value and entries
		value_reference(const value_reference& other)
			: value(other.value) {}
		__INLINE value_reference& operator = (INMOST_DATA_REAL_TYPE pvalue) 
		{ 
			value = pvalue; 
			return *this; 
		}
		__INLINE value_reference& operator = (value_reference const& other)
		{
			value = other.value;
			return *this;
		}
		__INLINE value_reference& operator +=(INMOST_DATA_REAL_TYPE right)
		{
			value += right;
			return *this;
		}
		__INLINE value_reference& operator -=(INMOST_DATA_REAL_TYPE right)
		{
			value -= right;
			return *this;
		}
		__INLINE value_reference& operator *=(INMOST_DATA_REAL_TYPE right)
		{
			value *= right;
			return *this;
		}
		__INLINE value_reference& operator /=(INMOST_DATA_REAL_TYPE right)
		{
			value /= right;
			return *this;
		}
		bool check_nans() const
		{
			if (value != value) return true;
			return false;
		}
		bool check_infs() const
		{
			if (__isinf__(value)) return true;
			return false;
		}
		operator INMOST_DATA_REAL_TYPE() const { return value; }
		operator INMOST_DATA_REAL_TYPE& () { return value; }
	};
	
	
	
	class hessian_multivar_expression_reference : public shell_expression<hessian_multivar_expression_reference>
	{
		INMOST_DATA_REAL_TYPE & value;
		Sparse::Row * entries;
		Sparse::HessianRow * hentries;
	public:
		/// Constructor, set links to the provided value and entries
		hessian_multivar_expression_reference(INMOST_DATA_REAL_TYPE & _value, Sparse::Row * _entries, Sparse::HessianRow * _hentries)
		: value(_value), entries(_entries), hentries(_hentries) {}
		/// Copy constructor, sets links to the same reference of value and entries
		hessian_multivar_expression_reference(const hessian_multivar_expression_reference & other)
		: value(other.value), entries(other.entries), hentries(other.hentries) {}
		/// Copy constructor from hessian_multivar_expression, sets links to the same reference of value and entries
		hessian_multivar_expression_reference(hessian_multivar_expression & other)
		: value(other.value), entries(&other.entries), hentries(&other.hessian_entries) {}
		/// Retrieve value
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		/// Set value without changing derivatives
		__INLINE void SetValue(INMOST_DATA_REAL_TYPE val) { value = val; }
		/// Retrieve derivatives with multiplier into Sparse::RowMerger structure.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			r.AddRow(mult,*entries);
		}
		/// Retrieve derivatives with multiplier into array.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			for (Sparse::Row::const_iterator it = entries->Begin(); it != entries->End(); ++it)
				r[it->first] += it->second * mult;
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J,INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			J = *entries;
			if( !J.isSorted() ) std::sort(J.Begin(),J.End());
			for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second *= multJ;
			H = *hentries;
			for(Sparse::HessianRow::iterator it = H.Begin(); it != H.End(); ++it) it->second *= multH;
		}
		__INLINE INMOST_DATA_REAL_TYPE GetDerivative(INMOST_DATA_ENUM_TYPE index) const {return GetRow().get_safe(index);}
		__INLINE multivar_expression GetVariable(INMOST_DATA_ENUM_TYPE index)
		{
			multivar_expression ret(0);
			for(int k = 0; k < (int)entries->Size(); ++k)
			{
				if( entries->GetIndex(k) == index )
				{
					ret.SetValue(entries->GetValue(k));
					Sparse::Row & r = ret.GetRow();
					for(int q = 0; q < (int)hentries->Size(); ++q)
					{
						Sparse::HessianRow::index & i = hentries->GetIndex(q);
						if( i.first == index )
							r.Push(i.second,hentries->GetValue(q)*(i.first == i.second ? 1.0 : 0.5));
						else if( i.second == index )
							r.Push(i.first,hentries->GetValue(q)*(i.first == i.second ? 1.0 : 0.5));
					}
					break;
				}
			}
			return ret;
		}
		__INLINE hessian_multivar_expression_reference & operator = (INMOST_DATA_REAL_TYPE pvalue)
		{
			value = pvalue;
			entries->Clear();
			hentries->Clear();
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator = (basic_expression const & expr)
		{
			value = expr.GetValue();
			Sparse::Row tmp;
			Sparse::HessianRow htmp;
			expr.GetHessian(1.0,tmp,1.0,htmp);
			entries->Swap(tmp);
			hentries->Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator = (multivar_expression_reference const & other)
		{
			value = other.GetValue();
			*entries = other.GetRow();
			hentries->Clear();
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator = (multivar_expression const & other)
		{
			value = other.GetValue();
			*entries = other.GetRow();
			hentries->Clear();
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator = (hessian_multivar_expression_reference const & other)
		{
			value = other.GetValue();
			*entries = other.GetRow();
			*hentries = other.GetHessianRow();
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator = (hessian_multivar_expression const & other)
		{
			value = other.GetValue();
			*entries = other.GetRow();
			*hentries = other.GetHessianRow();
			return *this;
		}
		__INLINE Sparse::Row & GetRow() {return *entries;}
		__INLINE Sparse::HessianRow & GetHessianRow() {return *hentries;}
		__INLINE const Sparse::Row & GetRow() const {return *entries;}
		__INLINE const Sparse::HessianRow & GetHessianRow() const {return *hentries;}
		__INLINE hessian_multivar_expression_reference & operator +=(basic_expression const & expr)
		{
			value += expr.GetValue();
			Sparse::Row tmpr, tmp;
			Sparse::HessianRow htmpr, htmp;
			expr.GetHessian(1.0,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(1.0,*entries,1.0,tmpr,tmp);
			entries->Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(1.0,*hentries,1.0,htmpr,htmp);
			hentries->Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator -=(basic_expression const & expr)
		{
			value -= expr.GetValue();
			Sparse::Row tmpr, tmp;
			Sparse::HessianRow htmpr, htmp;
			expr.GetHessian(1.0,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(1.0,*entries,-1.0,tmpr,tmp);
			entries->Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(1.0,*hentries,-1.0,htmpr,htmp);
			hentries->Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator *=(basic_expression const & expr)
		{
			throw NotImplemented; //check code below is correct
			INMOST_DATA_REAL_TYPE lval = value, rval = expr.GetValue();
			Sparse::Row tmp, tmpr;
			Sparse::HessianRow htmp, htmpr;
			expr.GetHessian(1.0,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(rval,*entries,lval,tmpr,tmp);
			entries->Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(lval,*hentries,rval,htmpr,htmp);
			hentries->Swap(htmp);
			value *= rval;
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator /=(basic_expression const & expr)
		{
			throw NotImplemented; //check code below is correct
			INMOST_DATA_REAL_TYPE rval = expr.GetValue();
			INMOST_DATA_REAL_TYPE reciprocial_rval = 1.0/rval;
			value *= reciprocial_rval;
			Sparse::Row tmp, tmpr;
			Sparse::HessianRow htmp, htmpr;
			expr.GetHessian(-value*reciprocial_rval,tmpr,1.0,htmpr);
			Sparse::Row::MergeSortedRows(reciprocial_rval,*entries,-value*reciprocial_rval,tmpr,tmp);
			entries->Swap(tmp);
			Sparse::HessianRow::MergeSortedRows(reciprocial_rval,*hentries,-value*reciprocial_rval,htmpr,htmp);
			hentries->Swap(htmp);
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator +=(INMOST_DATA_REAL_TYPE right)
		{
			value += right;
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator -=(INMOST_DATA_REAL_TYPE right)
		{
			value -= right;
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator *=(INMOST_DATA_REAL_TYPE right)
		{
			value *= right;
			for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it) it->second *= right;
			for(Sparse::HessianRow::iterator it = hentries->Begin(); it != hentries->End(); ++it)
				it->second *= right;
			return *this;
		}
		__INLINE hessian_multivar_expression_reference & operator /=(INMOST_DATA_REAL_TYPE right)
		{
			value /= right;
			for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it) it->second /= right;
			for(Sparse::HessianRow::iterator it = hentries->Begin(); it != hentries->End(); ++it) it->second /= right;
			return *this;
		}
		bool check_nans() const
		{
			if( value != value ) return true;
			for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it)
				if( it->second != it->second ) return true;
			for(Sparse::HessianRow::iterator it = hentries->Begin(); it != hentries->End(); ++it)
				if( it->second != it->second ) return true;
			return false;
		}
		bool check_infs() const
		{
			if( __isinf__(value) ) return true;
			for(Sparse::Row::iterator it = entries->Begin(); it != entries->End(); ++it)
				if( __isinf__(it->second) ) return true;
			for(Sparse::HessianRow::iterator it = hentries->Begin(); it != hentries->End(); ++it)
				if( __isinf__(it->second) ) return true;
			return false;
		}
	};
	

	class const_multivar_expression_reference : public shell_expression<const_multivar_expression_reference>
	{
		const multivar_expression& arg;
	public:
		/// Copy constructor, sets links to the same reference of value and entries
		const_multivar_expression_reference(const const_multivar_expression_reference& other)
			: arg(other.arg) {}
		/// Copy constructor from multivar_expression, sets links to the same reference of value and entries
		const_multivar_expression_reference(const multivar_expression& other)
			: arg(other) {}
		/// Retrieve value
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return arg.GetValue(); }
		/// Retrieve derivatives with multiplier into Sparse::RowMerger structure.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const {	arg.GetJacobian(mult, r); }
		/// Retrieve derivatives with multiplier into array.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const {	arg.GetJacobian(mult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const
		{
			arg.GetHessian(multJ, J, multH, H);
		}
	};

	class const_hessian_multivar_expression_reference : public shell_expression<const_hessian_multivar_expression_reference>
	{
		const hessian_multivar_expression& arg;
	public:
		/// Copy constructor, sets links to the same reference of value and entries
		const_hessian_multivar_expression_reference(const const_hessian_multivar_expression_reference& other)
			: arg(other.arg) {}
		/// Copy constructor from multivar_expression, sets links to the same reference of value and entries
		const_hessian_multivar_expression_reference(const hessian_multivar_expression& other)
			: arg(other) {}
		/// Retrieve value
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return arg.GetValue(); }
		/// Retrieve derivatives with multiplier into Sparse::RowMerger structure.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const
		{
			arg.GetJacobian(mult, r);
		}
		/// Retrieve derivatives with multiplier into array.
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const
		{
			arg.GetJacobian(mult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const
		{
			arg.GetHessian(multJ, J, multH, H);
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
			value = arg.GetValue() * dmult;
		}
		const_multiplication_expression(const const_multiplication_expression & other) : arg(other.arg), value(other.value), dmult(other.dmult) {}
        const_multiplication_expression(const const_multiplication_expression & other, const A & parg) : arg(parg), value(other.value), dmult(other.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
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
        variation_multiplication_expression(const variation_multiplication_expression & other, const A & parg) : arg(parg), value(other.value), dmult(other.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			arg.GetHessian(multJ * dmult, J, multH * dmult, H);
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
			dmult = 1.0 / dmult;
			value = arg.GetValue() * dmult;
		}
		const_division_expression(const const_division_expression & other) : arg(other.arg), value(other.value), dmult(other.dmult) {}
        const_division_expression(const const_division_expression & other, const A & parg) : arg(parg), value(other.value), dmult(other.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
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
			value = arg.GetValue() + padd;
		}
		const_addition_expression(const const_addition_expression & other) : arg(other.arg), value(other.value) {}
        const_addition_expression(const const_addition_expression & other, const A & parg) : arg(parg), value(other.value) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult, r);
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
        const_subtraction_expression(const const_subtraction_expression & other, const A & parg) : arg(parg), value(other.value) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(-mult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(-mult, r);
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
        reciprocal_expression(const reciprocal_expression & other, const A & parg)
                : arg(parg), value(other.value), reciprocial_val(other.reciprocial_val) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(-mult * value * reciprocial_val, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(-mult * value * reciprocial_val, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			Sparse::HessianRow ArgH;
			INMOST_DATA_REAL_TYPE coefJ = -multJ*value*reciprocial_val;
			INMOST_DATA_REAL_TYPE signJ = coefJ < 0 ? -1 : 1;
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
        unary_minus_expression(const unary_minus_expression & b, const A & parg) : arg(parg) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(-mult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(-mult, r);
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
        unary_plus_expression(const unary_plus_expression & b, const A & parg) : arg(parg) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult, r);
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
        abs_expression(const abs_expression & b, const A & parg) : arg(parg), value(b.value), dmult(b.dmult) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian( (value == 0 ? (mult < 0.0 ? -1 : 1) : 1) * mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian((value == 0 ? (mult < 0.0 ? -1 : 1) : 1) * mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ,Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			INMOST_DATA_REAL_TYPE a = (value == 0 ? (multJ < 0.0 ? -1 : 1) : 1);
			INMOST_DATA_REAL_TYPE b = (value == 0 ? (multH < 0.0 ? -1 : 1) : 1);
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
        exp_expression(const exp_expression & b, const A & parg) : arg(parg), value(b.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * value, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * value, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			Sparse::HessianRow ArgH;
			INMOST_DATA_REAL_TYPE coefJ = multJ*value;
			INMOST_DATA_REAL_TYPE signJ = coefJ < 0.0 ? -1 : 1;
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
        log_expression(const log_expression & b, const A & parg) : arg(parg), value(b.value), dmult(b.dmult) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			Sparse::HessianRow ArgH;
			INMOST_DATA_REAL_TYPE coefJ = multJ*dmult;
			INMOST_DATA_REAL_TYPE signJ = coefJ < 0.0 ? -1 : 1;
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
        sin_expression(const sin_expression & b, const A & parg) : arg(parg), value(b.value), dmult(b.dmult) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			Sparse::HessianRow htmp;
			arg.GetHessian(1,J,1,htmp);
			assert(J.isSorted());
			assert(htmp.isSorted());
			Sparse::HessianRow::MergeJacobianHessian(-value*multH,J,J,dmult*multH,htmp,H);
			for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second*=dmult*multJ;
			assert(H.isSorted());
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
        cos_expression(const cos_expression & b, const A & parg) : arg(parg), value(b.value), dmult(b.dmult) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			//arg.GetHessian(multJ*dmult,J,-multH*value,H);
			Sparse::HessianRow htmp;
			arg.GetHessian(1,J,1,htmp);
			Sparse::HessianRow::MergeJacobianHessian(-value*multH,J,J,dmult*multH,htmp,H);
			for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second*=dmult*multJ;
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
        sqrt_expression(const sqrt_expression & b, const A & parg) : arg(parg), value(b.value) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			if (value) arg.GetJacobian(0.5 * mult / value, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			if (value) arg.GetJacobian(0.5 * mult / value, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			//general formula:
			// (F(G))'' = F'(G) G'' + F''(G) G'.G'
			if( value )
			{
				Sparse::HessianRow htmp;
				arg.GetHessian(1,J,1,htmp);
				Sparse::HessianRow::MergeJacobianHessian(-0.25/::pow(value,3.0)*multH,J,J,0.5/value*multH,htmp,H);
				for(Sparse::Row::iterator it = J.Begin(); it != J.End(); ++it) it->second *= 0.5/value*multJ;
			}
			//arg.GetHessian(0.5*multJ/value,J,-0.25*multH/::pow(value,3),H);
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
			dmult = value ? lval/value : 0.0;
		}
		soft_abs_expression(const soft_abs_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
        soft_abs_expression(const soft_abs_expression & b, const A & parg) : arg(parg), value(b.value), dmult(b.dmult) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			(void)multJ,(void)J,(void)multH,(void)H;
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
			value = sdiv ? lval/sdiv : 0.0;
			dmult = sdiv ? (1.0 - (div ? lval2/div : 0.0))/sdiv : 0.0;
		}
		soft_sign_expression(const soft_sign_expression & b) : arg(b.arg), value(b.value), dmult(b.dmult) {}
        soft_sign_expression(const soft_sign_expression & b, const A & parg) : arg(parg), value(b.value), dmult(b.dmult) {}
        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; };
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			(void)multJ,(void)J,(void)multH,(void)H;
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
			if (root)
			{
				ldmult = 0.5 * (1 + diff / root);
				rdmult = 0.5 * (1 - diff / root);
			}
			else ldmult = rdmult = 0.5;
		}
		soft_max_expression(const soft_max_expression & other)
		: left(other.left), right(other.right), value(other.value), ldmult(other.ldmult), rdmult(other.rdmult) {}
        soft_max_expression(const soft_max_expression & other, const A & pleft, const B & pright)
                : left(pleft), right(pright), value(other.value), ldmult(other.ldmult), rdmult(other.rdmult) {}

        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult*ldmult,r);
			right.GetJacobian(mult*rdmult,r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult * ldmult, r);
			right.GetJacobian(mult * rdmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			(void)multJ,(void)J,(void)multH,(void)H;
			throw NotImplemented;
		}
	};

	template<class A>
	class soft_max_const_expression : public shell_expression<soft_max_const_expression<A> >
	{
		const A& left;
		INMOST_DATA_REAL_TYPE right;
		INMOST_DATA_REAL_TYPE value, ldmult;
	public:
		soft_max_const_expression(const shell_expression<A>& pleft, INMOST_DATA_REAL_TYPE pright, INMOST_DATA_REAL_TYPE tol) : left(pleft), right(pright)
		{
			INMOST_DATA_REAL_TYPE lval = left.GetValue(), rval = right;
			INMOST_DATA_REAL_TYPE diff = lval - rval, root = ::sqrt(diff * diff + tol * tol);
			value = 0.5 * (lval + rval + root);
			ldmult = 0.5 * (1 + (root ? diff / root : 0.0));
		}
		soft_max_const_expression(const soft_max_const_expression& other)
			: left(other.left), right(other.right), value(other.value), ldmult(other.ldmult) {}
		soft_max_const_expression(const soft_max_const_expression& other, const A& pleft, INMOST_DATA_REAL_TYPE pright)
			: left(pleft), right(pright), value(other.value), ldmult(other.ldmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const {left.GetJacobian(mult * ldmult, r);}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const { left.GetJacobian(mult * ldmult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const
		{
			(void)multJ, (void)J, (void)multH, (void)H;
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
			if (root)
			{
				ldmult = 0.5 * (1 - diff / root);
				rdmult = 0.5 * (1 + diff / root);
			}
			else ldmult = rdmult = 0.5;
		}
		soft_min_expression(const soft_min_expression & other)
		: left(other.left), right(other.right), value(other.value), ldmult(other.ldmult), rdmult(other.rdmult) {}
        soft_min_expression(const soft_min_expression & other, const A & pleft, const B & pright)
                : left(pleft), right(pright), value(other.value), ldmult(other.ldmult), rdmult(other.rdmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult*ldmult,r);
			right.GetJacobian(mult*rdmult,r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult * ldmult, r);
			right.GetJacobian(mult * rdmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			(void)multJ,(void)J,(void)multH,(void)H;
			throw NotImplemented;
		}
	};
	

	template<class A>
	class soft_min_const_expression : public shell_expression<soft_min_const_expression<A> >
	{
		const A& left;
		INMOST_DATA_REAL_TYPE right;
		INMOST_DATA_REAL_TYPE value, ldmult;
	public:
		soft_min_const_expression(const shell_expression<A>& pleft, INMOST_DATA_REAL_TYPE pright, INMOST_DATA_REAL_TYPE tol) : left(pleft), right(pright)
		{
			INMOST_DATA_REAL_TYPE lval = left.GetValue(), rval = right;
			INMOST_DATA_REAL_TYPE diff = lval - rval, root = ::sqrt(diff * diff + tol * tol);
			value = 0.5 * (lval + rval - root);
			ldmult = 0.5 * (1 - (root ? diff / root : 0.0));
		}
		soft_min_const_expression(const soft_min_const_expression& other)
			: left(other.left), right(other.right), value(other.value), ldmult(other.ldmult) {}
		soft_min_const_expression(const soft_min_const_expression& other, const A& pleft, INMOST_DATA_REAL_TYPE pright)
			: left(pleft), right(pright), value(other.value), ldmult(other.ldmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const {left.GetJacobian(mult * ldmult, r);}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const { left.GetJacobian(mult * ldmult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const
		{
			(void)multJ, (void)J, (void)multH, (void)H;
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
        multiplication_expression(const multiplication_expression & other, const A & pleft, const B & pright)
                : left(pleft), right(pright), value(other.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult*right.GetValue(),r);
			right.GetJacobian(mult*left.GetValue(),r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult * right.GetValue(), r);
			right.GetJacobian(mult * left.GetValue(), r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			// (F*G)'' = (F'G+G'F)' = (F''G + F'G' + G''F + G'F') = (F''G + G''F + 2F'G')
			Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
			Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
			left.GetHessian(1,JL,1,HL); //retrieve jacobian row and hessian matrix of the left expression
			assert(JL.isSorted());
			assert(HL.isSorted());
			right.GetHessian(1,JR,1,HR); //retrieve jacobian row and hessian matrix of the right expression
			assert(JR.isSorted());
			assert(HR.isSorted());
			//assume rows are sorted (this is to be ensured by corresponding GetHessian functions)
			//preallocate J to JL.Size+JR.Size
			//perform merging of two sorted arrays
			//resize to correct size
			Sparse::Row::MergeSortedRows(right.GetValue()*multJ,JL,left.GetValue()*multJ,JR,J);
			assert(J.isSorted());
			//preallocate H to HL.Size+HR.Size+JL.Size*JR.Size
			//merge sorted
			Sparse::HessianRow::MergeJacobianHessian(2.0*multH,JL,JR,right.GetValue()*multH,HL,left.GetValue()*multH,HR,H);
			assert(H.isSorted());
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
			if( rval )
			{
				reciprocal_rval = 1.0 / rval;
				value = lval * reciprocal_rval;
			}
			else
			{
				reciprocal_rval = 0;
				value = 0;
			}
		}
		division_expression(const division_expression & other) : left(other.left), right(other.right), value(other.value), reciprocal_rval(other.reciprocal_rval) {}
        division_expression(const division_expression & other, const A & pleft, const B & pright) :
                left(pleft), right(pright), value(other.value), reciprocal_rval(other.reciprocal_rval) {}

        __INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult * reciprocal_rval,r);
			right.GetJacobian(- mult * value * reciprocal_rval,r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult * reciprocal_rval, r);
			right.GetJacobian(-mult * value * reciprocal_rval, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			throw NotImplemented;
			Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
			Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
			left.GetHessian(multJ,JL,multH,HL); //retrieve jacobian row and hessian matrix of the left expression
			right.GetHessian(multJ,JR,multH,HR); //retrieve jacobian row and hessian matrix of the right expression
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
        addition_expression(const addition_expression & other, const A & pleft, const B & pright)
                : left(pleft), right(pright), value(other.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult, r);
			right.GetJacobian(mult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult, r);
			right.GetJacobian(mult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
			Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
			left.GetHessian(1,JL,1,HL); //retrieve jacobian row and hessian matrix of the left expression
			assert(JL.isSorted());
			assert(HL.isSorted());
			right.GetHessian(1,JR,1,HR); //retrieve jacobian row and hessian matrix of the right expression
			assert(JR.isSorted());
			assert(HR.isSorted());
			Sparse::Row::MergeSortedRows(multJ,JL,multJ,JR,J);
			assert(J.isSorted());
			Sparse::HessianRow::MergeSortedRows(multH,HL,multH,HR,H);
			assert(H.isSorted());
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
        subtraction_expression(const subtraction_expression & other, const A & pleft, const B & pright)
                : left(pleft), right(pright),value(other.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult,r);
			right.GetJacobian(-mult,r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult, r);
			right.GetJacobian(-mult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			//(F-G)'' = (F'-G')' = (F''-G'')
			Sparse::Row JL, JR; //temporary jacobian rows from left and right expressions
			Sparse::HessianRow HL, HR; //temporary hessian rows form left and right expressions
			left.GetHessian(1,JL,1,HL); //retrieve jacobian row and hessian matrix of the left expression
			right.GetHessian(1,JR,1,HR); //retrieve jacobian row and hessian matrix of the right expression
			Sparse::Row::MergeSortedRows(multJ,JL,-multJ,JR,J);
			Sparse::HessianRow::MergeSortedRows(multH,HL,-multH,HR,H);
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
        pow_expression(const pow_expression & other, const A & pleft, const B & pright)
                :left(pleft), right(pright), value(other.value),
                 ldmult(other.ldmult), rdmult(other.rdmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult*ldmult,r);
			right.GetJacobian(mult*rdmult,r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult * ldmult, r);
			right.GetJacobian(mult * rdmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			throw NotImplemented;
		}
	};
	
	template<class A>
	class acos_expression : public shell_expression<acos_expression<A> >
	{
		const A& arg;
		INMOST_DATA_REAL_TYPE value, dmult;
	public:
		acos_expression(const shell_expression<A>& parg) : arg(parg)
		{
			INMOST_DATA_REAL_TYPE val = arg.GetValue();
			value = ::acos(val);
			dmult = -1.0 / sqrt(1.0 - val * val);
		}
		acos_expression(const acos_expression& other)
			:arg(other.arg), value(other.value), dmult(other.dmult) {}
		acos_expression(const acos_expression& other, const A& parg)
			:arg(parg), value(other.value), dmult(other.dmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const {	arg.GetJacobian(mult * dmult, r); }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const { arg.GetJacobian(mult * dmult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const	
		{ 
			throw NotImplemented; 
		}
	};
	
	template<class A, class B>
	class atan2_expression : public shell_expression<atan2_expression<A,B> >
	{
		const A & left;
		const B & right;
		INMOST_DATA_REAL_TYPE value, ldmult, rdmult;
	public:
		atan2_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright)
		{
			INMOST_DATA_REAL_TYPE lval = left.GetValue();
			INMOST_DATA_REAL_TYPE rval = right.GetValue();
			value = ::atan2(lval,rval);
			ldmult = rval/(rval*rval+lval*lval);
			rdmult = -lval/(rval*rval+lval*lval);
		}
		atan2_expression(const atan2_expression & other)
		:left(other.left), right(other.right), value(other.value),
		ldmult(other.ldmult), rdmult(other.rdmult) {}
        atan2_expression(const atan2_expression & other, const A & pleft, const B & pright)
                :left(pleft), right(pright), value(other.value),
                 ldmult(other.ldmult), rdmult(other.rdmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult * ldmult, r);
			right.GetJacobian(mult * rdmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult * ldmult, r);
			right.GetJacobian(mult * rdmult, r);
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
		INMOST_DATA_REAL_TYPE value, ldmult, ldmult2;
	public:
		pow_const_expression(const shell_expression<A> & pleft, INMOST_DATA_REAL_TYPE pright) : left(pleft)
		{
			INMOST_DATA_REAL_TYPE lval = left.GetValue();
			value = ::pow(lval,pright);
			if( lval != 0 )
			{
				ldmult = value * pright / lval;
				ldmult2 = ldmult * (pright-1) / lval;
			}
			else
				ldmult = ldmult2 = 0;
		}
		pow_const_expression(const pow_const_expression & other)
		:left(other.left), value(other.value), ldmult(other.ldmult) {}
        pow_const_expression(const pow_const_expression & other, const A & pleft)
                :left(pleft), value(other.value), ldmult(other.ldmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			left.GetJacobian(mult * ldmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			left.GetJacobian(mult * ldmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			//(F(G))'' = (F'(G)G')' = (F''(G)G'+F'(G)G'')
			//(g(x)^n)'' = n g(x)^(n - 2) (g(x) g''(x) + (n - 1) g'(x)^2)
			Sparse::Row JL;
			Sparse::HessianRow HL;
			left.GetHessian(1,JL,1,HL); //retrieve jacobian row and hessian matrix of the left expression
			Sparse::HessianRow::MergeJacobianHessian(multH*ldmult2,JL,JL,multH*ldmult,HL,H);
			//for(Sparse::Row::iterator it = JL.Begin(); it != JL.End(); ++it) it->second *= ldmult*multJ;
			for (Sparse::Row::iterator it = JL.Begin(); it != JL.End(); ++it) J[it->first] += it->second * ldmult * multJ;
            (void)J;
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
        const_pow_expression(const const_pow_expression & other, const A & pright)
                :right(pright), value(other.value), rdmult(other.rdmult) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			right.GetJacobian(mult * rdmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			right.GetJacobian(mult * rdmult, r);
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
        condition_expression(const condition_expression & other, const A & pcond, const B & pleft, const C & pright)
                :cond(pcond), left(pleft), right(pright),
                 value(other.value), cond_value(other.cond_value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			if( cond_value >= 0.0 )
				left.GetJacobian(mult,r);
			else
				right.GetJacobian(mult,r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			if (cond_value >= 0.0)
				left.GetJacobian(mult, r);
			else
				right.GetJacobian(mult, r);
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
		branch_expression(const shell_expression<A> & pleft, const shell_expression<B> & pright) : left(pleft), right(pright)
		{
			value = 0;
		}
		branch_expression(bool pcond, const shell_expression<A> & pleft, const shell_expression<B> & pright) : cond(pcond), left(pleft), right(pright)
		{
			value = cond ? left.GetValue() : right.GetValue();
		}
		branch_expression(const branch_expression & other)
		:cond(other.cond), left(other.left), right(other.right),
		value(other.value) {}
        branch_expression(const branch_expression & other, const A & pleft, const B & pright)
                :cond(other.cond), left(pleft), right(pright),
                 value(other.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			if (cond)
				left.GetJacobian(mult, r);
			else
				right.GetJacobian(mult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			if (cond)
				left.GetJacobian(mult, r);
			else
				right.GetJacobian(mult, r);
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
	class const_branch_expression : public shell_expression<const_branch_expression<A> >
	{
		bool cond;
		const A& left;
		INMOST_DATA_REAL_TYPE right;
		INMOST_DATA_REAL_TYPE value;
	public:
		const_branch_expression(const shell_expression<A>& pleft, INMOST_DATA_REAL_TYPE pright) : left(pleft), right(pright) {value = 0.0;}
		const_branch_expression(bool pcond, const shell_expression<A>& pleft, INMOST_DATA_REAL_TYPE pright) : cond(pcond), left(pleft), right(pright)
		{
			value = cond ? left.GetValue() : right;
		}
		const_branch_expression(const const_branch_expression& other)
			:cond(other.cond), left(other.left), right(other.right), value(other.value) {}
		const_branch_expression(const const_branch_expression& other, const A& pleft)
			:cond(other.cond), left(pleft), right(other.right), value(other.value) {}
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return value; }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const { if (cond) left.GetJacobian(mult, r);}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const { if (cond) left.GetJacobian(mult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const
		{
			if (cond) left.GetHessian(multJ, J, multH, H);
		}
		void SetCondition(bool _cond)
		{
			cond = _cond;
			value = cond ? left.GetValue() : right;
		}
	};
	
		
	template<class A>
	class function_expression : public shell_expression< function_expression<A> > 
	{
        const A &arg;
        INMOST_DATA_REAL_TYPE value, dmult, ddmult;
    public:
        function_expression(const shell_expression<A> &_arg)
                : arg(_arg), value(1), dmult(0) {}

        function_expression(const shell_expression<A> &_arg, INMOST_DATA_REAL_TYPE pvalue, INMOST_DATA_REAL_TYPE pdmult,
                            INMOST_DATA_REAL_TYPE pddmult = 0)
                : arg(_arg), value(pvalue), dmult(pdmult), ddmult(pddmult) {}

        function_expression(const function_expression &other)
                : arg(other.arg), value(other.value), dmult(other.dmult), ddmult(other.ddmult) {}
        function_expression(const function_expression &other, const A & parg)
                : arg(parg), value(other.value), dmult(other.dmult), ddmult(other.ddmult) {}

        function_expression &operator=(function_expression const &b)
        {
            arg = b.arg;
            value = b.value;
            dmult = b.dmult;
            ddmult = b.ddmult;
            return *this;
        }
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const {return value;}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger & r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE * r) const
		{
			arg.GetJacobian(mult * dmult, r);
		}
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row & J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow & H) const
		{
			arg.GetHessian(multJ * dmult, J, multH * ddmult, H);
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
				else break;
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
		keyval_table(std::string _name, const INMOST_DATA_REAL_TYPE * _args, const INMOST_DATA_REAL_TYPE * _vals, INMOST_DATA_ENUM_TYPE _size)
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

	template<class Op, class A>
	class unary_pool
	{
		A arg;
		Op operand;

		unary_pool& operator = (unary_pool const& other) { arg = other.arg; operand.assign(other.operand, arg); return *this; }
	public:
		unary_pool(const A& parg) : arg(parg), operand(arg) {}
		unary_pool(const A& parg, INMOST_DATA_REAL_TYPE pmult) : arg(parg), operand(arg, pmult) {}
		unary_pool(bool cond, const A& parg, INMOST_DATA_REAL_TYPE pright) : arg(parg), operand(cond, arg, pright) {}
		unary_pool(const unary_pool& other) : arg(other.arg), operand(other.operand, arg) {}

		const shell_expression<A>& get_arg() { return arg; }
		Op& get_op() { return operand; }
		const Op& get_op() const { return operand; }
	};



	template<class Op, class A, class B>
	class binary_pool
	{

		A left;
		B right;
		Op operand;

		binary_pool& operator = (binary_pool const& other) { left = other.left; right = other.right; operand.assign(other.operand, left, right); return *this; }
	public:
		binary_pool(const A& pleft, const B& pright) : left(pleft), right(pright), operand(left, right) {}
		binary_pool(bool cond, const A& pleft, const B& pright) : left(pleft), right(pright), operand(cond, left, right) {}
		binary_pool(const binary_pool& other) : left(other.left), right(other.right), operand(other.operand, left, right) {}

		const shell_expression<A>& get_left() { return left; }
		const shell_expression<B>& get_right() { return right; }
		Op& get_op() { return operand; }
		const Op& get_op() const { return operand; }
		~binary_pool() {}
	};

	template<class Op, class A, class B, class C>
	class ternary_pool
	{
		A cond;
		B left;
		C right;
		Op operand;

		ternary_pool& operator =(ternary_pool const& other) { cond = other.cond; left = other.left; right = other.right; operand.assign(other.operand, cond, left, right); return *this; }
	public:
		ternary_pool(const A& pcond, const B& pleft, const C& pright) : cond(pcond), left(pleft), right(pright), operand(cond, left, right) {}
		ternary_pool(const ternary_pool& other) : cond(other.cond), left(other.left), right(other.right), operand(other.operand, cond, left, right) {}

		const shell_expression<A>& get_cond() { return cond; }
		const shell_expression<B>& get_left() { return left; }
		const shell_expression<C>& get_right() { return right; }
		Op& get_op() { return operand; }
		const Op& get_op() const { return operand; }
		~ternary_pool() {}
	};

	template<class A, class ArgA>
	class unary_pool_expression : public shell_expression<unary_pool_expression<A, ArgA> >
	{
		unary_pool<A, ArgA> pool;
	public:
		unary_pool_expression(const unary_pool<A, ArgA>& ppool) : pool(ppool) {}
		unary_pool_expression(const unary_pool_expression& other) : pool(other.pool) {}
		unary_pool_expression& operator = (unary_pool_expression const& other) { pool = other.pool; return *this; }
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return pool.get_op().GetValue(); }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const { pool.get_op().GetJacobian(mult, r); }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const { pool.get_op().GetJacobian(mult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const { pool.get_op().GetHessian(multJ, J, multH, H); }
	};

	template<class A, class ArgA, class ArgB>
	class binary_pool_expression : public shell_expression<binary_pool_expression<A, ArgA, ArgB> >
	{
		binary_pool<A, ArgA, ArgB> pool;
	public:
		binary_pool_expression(const binary_pool<A, ArgA, ArgB>& ppool) : pool(ppool) {}
		binary_pool_expression(const binary_pool_expression& other) : pool(other.pool) {}
		binary_pool_expression& operator = (binary_pool_expression const& other) { pool = other.pool; return *this; }
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return pool.get_op().GetValue(); }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const { pool.get_op().GetJacobian(mult, r); }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const { pool.get_op().GetJacobian(mult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const { pool.get_op().GetHessian(multJ, J, multH, H); }
	};

	template<class A, class ArgA, class ArgB, class ArgC>
	class ternary_pool_expression : public shell_expression<ternary_pool_expression<A, ArgA, ArgB, ArgC> >
	{
		ternary_pool<A, ArgA, ArgB, ArgC> pool;
	public:
		ternary_pool_expression(const ternary_pool<A, ArgA, ArgB, ArgC>& ppool) : pool(ppool) {}
		ternary_pool_expression(const ternary_pool_expression& other) : pool(other.pool) {}
		ternary_pool_expression& operator = (ternary_pool_expression const& other) { pool = other.pool; return *this; }
		__INLINE INMOST_DATA_REAL_TYPE GetValue() const { return pool.get_op().GetValue(); }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, Sparse::RowMerger& r) const { pool.get_op().GetJacobian(mult, r); }
		__INLINE void GetJacobian(INMOST_DATA_REAL_TYPE mult, INMOST_DATA_REAL_TYPE* r) const { pool.get_op().GetJacobian(mult, r); }
		__INLINE void GetHessian(INMOST_DATA_REAL_TYPE multJ, Sparse::Row& J, INMOST_DATA_REAL_TYPE multH, Sparse::HessianRow& H) const { pool.get_op().GetHessian(multJ, J, multH, H); }
	};

	
	typedef multivar_expression variable;
	typedef multivar_dense_expression variable_dense;
	typedef hessian_multivar_expression hessian_variable;
	typedef var_expression unknown;
}


__INLINE bool check_nans(INMOST_DATA_REAL_TYPE val) {return val != val;}
__INLINE bool check_nans(INMOST::var_expression const & e) {return e.check_nans();}
__INLINE bool check_nans(INMOST::multivar_expression const & e) {return e.check_nans();}
__INLINE bool check_nans(INMOST::multivar_expression_reference const & e) {return e.check_nans();}
__INLINE bool check_infs(INMOST_DATA_REAL_TYPE val) {return __isinf__(val);}
__INLINE bool check_infs(INMOST::var_expression const & e) {return e.check_infs();}
__INLINE bool check_infs(INMOST::multivar_expression const & e) {return e.check_infs();}
__INLINE bool check_infs(INMOST::multivar_expression_reference const & e) {return e.check_infs();}
__INLINE bool check_nans_infs(INMOST_DATA_REAL_TYPE val) {return check_nans(val) || check_infs(val);}
__INLINE bool check_nans_infs(INMOST::var_expression const & e) {return e.check_nans() || e.check_infs();}
__INLINE bool check_nans_infs(INMOST::multivar_expression const & e) {return e.check_nans() || e.check_infs();}
__INLINE bool check_nans_infs(INMOST::multivar_expression_reference const & e) {return e.check_nans() || e.check_infs();}


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
template<class A>          __INLINE                     INMOST::acos_expression<A>      acos(INMOST::shell_expression<A> const & Arg) { return INMOST::acos_expression<A>(Arg); }
__INLINE                          INMOST_DATA_REAL_TYPE variation(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE) {return Arg;}
template<class A>          __INLINE                          INMOST_DATA_REAL_TYPE get_value(INMOST::shell_expression<A> const & Arg) { return Arg.GetValue(); }
                           __INLINE                          INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE Arg) {return Arg;}
						   __INLINE                          INMOST_DATA_CPLX_TYPE get_value(INMOST_DATA_CPLX_TYPE Arg) { return Arg; }
						   __INLINE                          INMOST_DATA_CPLX_TYPE get_value(std::complex<INMOST::var_expression> const& Arg) { return INMOST_DATA_CPLX_TYPE(Arg.real().GetValue(), Arg.imag().GetValue()); }
						   __INLINE                          INMOST_DATA_CPLX_TYPE get_value(std::complex<INMOST::multivar_expression> const& Arg) { return INMOST_DATA_CPLX_TYPE(Arg.real().GetValue(), Arg.imag().GetValue()); }
						   __INLINE                          INMOST_DATA_CPLX_TYPE get_value(std::complex<INMOST::hessian_multivar_expression> const& Arg) { return INMOST_DATA_CPLX_TYPE(Arg.real().GetValue(), Arg.imag().GetValue()); }
						   __INLINE                  const INMOST::var_expression& conj(const INMOST::var_expression& Arg) { return Arg; }
						   __INLINE             const INMOST::multivar_expression& conj(const INMOST::multivar_expression& Arg) { return Arg; }
						   __INLINE   const INMOST::multivar_expression_reference& conj(const INMOST::multivar_expression_reference& Arg) { return Arg; }
						   __INLINE     const INMOST::hessian_multivar_expression& conj(const INMOST::hessian_multivar_expression& Arg) { return Arg; }
						   __INLINE     const INMOST::hessian_multivar_expression_reference& conj(const INMOST::hessian_multivar_expression_reference& Arg) { return Arg; }
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
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = Val;}
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val;}
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, const INMOST::var_expression & Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val.GetValue();}
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, const INMOST::multivar_expression & Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val.GetValue();}
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, const INMOST::multivar_expression_reference & Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val.GetValue();}
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, const INMOST::hessian_multivar_expression & Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val.GetValue(); }
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, const INMOST::hessian_multivar_expression_reference & Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val.GetValue(); }
                           __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = (INMOST_DATA_REAL_TYPE)Val;}
                           __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val;}
						   __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::var_expression & Val) {Arg = Val.GetValue();}
						   __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::multivar_expression & Val) {Arg = Val.GetValue();}
						   __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val.GetValue();}
                           __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::hessian_multivar_expression & Val) {Arg = Val.GetValue(); }
                           __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::hessian_multivar_expression_reference & Val) {Arg = Val.GetValue(); }
//						   __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::value_reference& Val) { Arg = Val.GetValue(); }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, const INMOST_DATA_CPLX_TYPE& Val) { Arg = Val; }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, INMOST_DATA_INTEGER_TYPE Val) { Arg = (INMOST_DATA_REAL_TYPE)Val; }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, INMOST_DATA_REAL_TYPE Val) { Arg = Val; }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, const INMOST::var_expression& Val) { Arg = Val.GetValue(); }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, const INMOST::multivar_expression& Val) { Arg = Val.GetValue(); }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, const INMOST::multivar_expression_reference& Val) { Arg = Val.GetValue(); }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, const INMOST::hessian_multivar_expression& Val) { Arg = Val.GetValue(); }
						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, const INMOST::hessian_multivar_expression_reference& Val) { Arg = Val.GetValue(); }
//						   __INLINE                                           void    assign(INMOST_DATA_CPLX_TYPE& Arg, const INMOST::value_reference& Val) { Arg = Val.GetValue(); }
                           __INLINE                                           void    assign(INMOST::var_expression & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::var_expression & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::var_expression & Arg, const INMOST::var_expression & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::var_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::multivar_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::hessian_multivar_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::hessian_multivar_expression_reference & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::var_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::multivar_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::hessian_multivar_expression & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, const INMOST::var_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, const INMOST::multivar_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, const INMOST::hessian_multivar_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, const INMOST::hessian_multivar_expression_reference & Val) {Arg = Val; }
                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, const INMOST::var_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, const INMOST::multivar_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, const INMOST::multivar_expression_reference & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, const INMOST::hessian_multivar_expression & Val) {Arg = Val; }
//                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, const INMOST::hessian_multivar_expression_reference & Val) {Arg = Val; }
template<class A>          __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg, const INMOST::shell_expression<A> & Val) {Arg = (INMOST_DATA_REAL_TYPE)Val.GetValue();}
template<class A>          __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg, const INMOST::shell_expression<A> & Val) {Arg = Val.GetValue();}
template<class A>          __INLINE                                           void    assign(INMOST::multivar_expression & Arg, const INMOST::shell_expression<A> & Val) {Arg = Val;}
//template<class A>          __INLINE                                           void    assign(INMOST::value_reference& Arg, const INMOST::shell_expression<A>& Val) {Arg.SetValue(Val.GetValue());}
template<class A>          __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg, const INMOST::shell_expression<A> & Val) {Arg = Val;}
template<class A>          __INLINE                                           void    assign(INMOST::hessian_multivar_expression & Arg, const INMOST::shell_expression<A> & Val) {Arg = Val;}
template<class A>          __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, const INMOST::shell_expression<A> & Val) {Arg = Val;}
//                           __INLINE                                           void    assign(INMOST::value_reference& Arg, INMOST_DATA_REAL_TYPE Val) { Arg = (INMOST_DATA_REAL_TYPE)Val; }
#if defined(USE_FP64)
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg,                      float Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg,                         float Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
						   __INLINE                                           void    assign(float& Arg,                                          float Val) {Arg = (float)Val; }
                           __INLINE                                           void    assign(INMOST::var_expression & Arg,                        float Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg,                   float Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg,         float Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
//						   __INLINE                                           void    assign(INMOST::value_reference& Arg,                        float Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, float Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
#else //USE_FP64
                           __INLINE                                           void    assign(INMOST_DATA_INTEGER_TYPE & Arg,                      double Val) {Arg = (INMOST_DATA_INTEGER_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST_DATA_REAL_TYPE & Arg,                         double Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
						   __INLINE                                           void    assign(double& Arg,                                         double Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::var_expression & Arg,                        double Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression & Arg,                   double Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::multivar_expression_reference & Arg,         double Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
                           __INLINE                                           void    assign(INMOST::hessian_multivar_expression_reference & Arg, double Val) {Arg = (INMOST_DATA_REAL_TYPE)Val; }
#endif //USE_FP64
template<class A>          __INLINE                 INMOST::soft_abs_expression<A> soft_fabs(INMOST::shell_expression<A> const & Arg, INMOST_DATA_REAL_TYPE tol = 0) { return INMOST::soft_abs_expression<A>(Arg,tol); }
__INLINE                                                     INMOST_DATA_REAL_TYPE soft_fabs(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol = 0) {return ::sqrt(Arg*Arg+tol*tol);}
template<class A>          __INLINE                INMOST::soft_sign_expression<A> soft_sign(INMOST::shell_expression<A> const & Arg, INMOST_DATA_REAL_TYPE tol = 0) { return INMOST::soft_sign_expression<A>(Arg,tol); }
__INLINE                                                     INMOST_DATA_REAL_TYPE soft_sign(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol = 0) {return Arg/::sqrt(Arg*Arg+tol*tol);}
template<class A, class B> __INLINE        INMOST::multiplication_expression<A, B> operator*(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::multiplication_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE              INMOST::division_expression<A, B> operator/(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::division_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE              INMOST::addition_expression<A, B> operator+(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::addition_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE           INMOST::subtraction_expression<A, B> operator-(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::subtraction_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE                   INMOST::pow_expression<A, B>       pow(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::pow_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE                 INMOST::atan2_expression<A, B>     atan2(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right) { return INMOST::atan2_expression<A, B> (Left, Right); }
template<class A, class B> __INLINE              INMOST::soft_max_expression<A, B>  soft_max(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right ,INMOST_DATA_REAL_TYPE tol = 0.0) { return INMOST::soft_max_expression<A, B> (Left, Right,tol); }
template<class A>          __INLINE           INMOST::soft_max_const_expression<A>  soft_max(INMOST::shell_expression<A> const& Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol = 0.0) { return INMOST::soft_max_const_expression<A>(Left, Right, tol); }
template<class B>          __INLINE           INMOST::soft_max_const_expression<B>  soft_max(INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const& Right, INMOST_DATA_REAL_TYPE tol = 0.0) { return INMOST::soft_max_const_expression<B>(Right, Left, tol); }
                           __INLINE                          INMOST_DATA_REAL_TYPE  soft_max(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol = 0.0) {return 0.5*(Left+Right+::sqrt((Left-Right)*(Left-Right)+tol*tol));}
template<class A, class B> __INLINE              INMOST::soft_min_expression<A, B>  soft_min(INMOST::shell_expression<A> const & Left, INMOST::shell_expression<B> const & Right ,INMOST_DATA_REAL_TYPE tol = 0.0) { return INMOST::soft_min_expression<A, B> (Left, Right,tol); }
template<class A>          __INLINE           INMOST::soft_min_const_expression<A>  soft_min(INMOST::shell_expression<A> const& Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol = 0.0) { return INMOST::soft_min_const_expression<A>(Left, Right, tol); }
template<class B>          __INLINE           INMOST::soft_min_const_expression<B>  soft_min(INMOST_DATA_REAL_TYPE Left, INMOST::shell_expression<B> const& Right, INMOST_DATA_REAL_TYPE tol = 0.0) { return INMOST::soft_min_const_expression<B>(Right, Left, tol); }
                           __INLINE                          INMOST_DATA_REAL_TYPE  soft_min(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol = 0.0) {return 0.5*(Left+Right-::sqrt((Left-Right)*(Left-Right)+tol*tol));}
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
#else //USE_AUTODIFF
__INLINE bool check_nans(INMOST_DATA_REAL_TYPE val) {return val != val;}
__INLINE bool check_infs(INMOST_DATA_REAL_TYPE val) {return __isinf__(val);}
__INLINE bool check_nans_infs(INMOST_DATA_REAL_TYPE val) {return check_nans(val) || check_infs(val);}
__INLINE void                     assign(INMOST_DATA_INTEGER_TYPE & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = Val;}
__INLINE void                     assign(INMOST_DATA_INTEGER_TYPE & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = static_cast<INMOST_DATA_INTEGER_TYPE>(Val);}
__INLINE void                     assign(INMOST_DATA_REAL_TYPE & Arg, INMOST_DATA_INTEGER_TYPE Val) {Arg = static_cast<INMOST_DATA_REAL_TYPE>(Val);}
__INLINE void                     assign(INMOST_DATA_REAL_TYPE & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val;}
__INLINE void                     assign(INMOST_DATA_CPLX_TYPE& Arg, INMOST_DATA_REAL_TYPE Val) { Arg = Val; }
__INLINE void                     assign(INMOST_DATA_CPLX_TYPE& Arg, INMOST_DATA_CPLX_TYPE Val) { Arg = Val; }
__INLINE void                  set_value(INMOST_DATA_REAL_TYPE & Arg, INMOST_DATA_REAL_TYPE Val) {Arg = Val;}
__INLINE INMOST_DATA_REAL_TYPE get_value(INMOST_DATA_REAL_TYPE Arg) {return Arg;}
__INLINE INMOST_DATA_CPLX_TYPE get_value(INMOST_DATA_CPLX_TYPE Arg) { return Arg; }
__INLINE INMOST_DATA_REAL_TYPE variation(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE) {return Arg;}
__INLINE INMOST_DATA_REAL_TYPE soft_fabs(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol = 0) {return ::sqrt(Arg*Arg+tol*tol);}
__INLINE INMOST_DATA_REAL_TYPE soft_sign(INMOST_DATA_REAL_TYPE Arg, INMOST_DATA_REAL_TYPE tol = 0) {return Arg/::sqrt(Arg*Arg+tol*tol);}
__INLINE INMOST_DATA_REAL_TYPE  soft_max(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol) {return 0.5*(Left+Right+::sqrt((Left-Right)*(Left-Right)+tol*tol));}
__INLINE INMOST_DATA_REAL_TYPE  soft_min(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right, INMOST_DATA_REAL_TYPE tol) {return 0.5*(Left+Right-::sqrt((Left-Right)*(Left-Right)+tol*tol));}
__INLINE INMOST_DATA_REAL_TYPE cond_zero(bool cond, INMOST_DATA_REAL_TYPE Arg) { return cond ? Arg : 0.0; }
__INLINE INMOST_DATA_REAL_TYPE cond_both(bool cond, INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) { return cond ? Left : Right; }
__INLINE INMOST_DATA_REAL_TYPE       add(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) { return Left + Right; }
__INLINE INMOST_DATA_REAL_TYPE       sub(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) { return Left - Right; }
__INLINE INMOST_DATA_REAL_TYPE       mul(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) { return Left * Right; }
__INLINE INMOST_DATA_REAL_TYPE       div(INMOST_DATA_REAL_TYPE Left, INMOST_DATA_REAL_TYPE Right) { return Left / Right; }
#endif //USE_AUTODIFF
#endif //INMOST_AUTODIFF_ETEXPR_H_INCLUDED
