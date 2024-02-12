#ifndef INMOST_SPARSE_INCLUDED
#define INMOST_SPARSE_INCLUDED


#include "inmost_common.h"
#include <unordered_map>
#include <queue>
#include "robin_hood.h"

#define ASSUME_SORTED
//#define TEST_HASHTABLE
//#define TEST_JUDY1

#if defined(TEST_JUDY1)
#include "judyLArray.h"
#endif
#if defined(TEST_HASHTABLE)
#include "hashtable.h"
#endif

namespace INMOST
{
	namespace Sparse
	{
		typedef bool bit_type; //for RowMerger2

#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		/// Retrieve MPI type for row entry type.
		INMOST_MPI_Type GetRowEntryType();
		/// Create MPI type for row entry type.
		void CreateRowEntryType();
		/// Release MPI type.
		void DestroyRowEntryType();
		/// Check whether MPI type was created.
		bool HaveRowEntryType();
		/// This class can be used to annotate the matrix.
		/// You can optionally provide a pointer to an object of
		/// this class to Sparse::Matrix::Save or Sparse::HessianMatrix::Save.
		/// Then comments would be added to mtx file next to the
		/// non-empty strings.
		class AnnotationService
		{
			interval<INMOST_DATA_ENUM_TYPE,std::string> text; ///< Array of strings corresponding to matrix rows.
		public:
			/// Create a new annotation for local matrix.
			/// @param start First row of the matrix.
			/// @param end Last row of the matrix.
			AnnotationService(INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0) { if( end != start ) SetInterval(start,end); }
			/// Copy annotations for the matrix.
			/// @param other Another class for annotation.
			AnnotationService(const AnnotationService & other) : text(other.text) { }
			/// Assign annotations for the matrix.
			/// @param other Another class for annotation.
			AnnotationService & operator = (AnnotationService const & other) {text = other.text; return *this;}
			/// Destroy annotations for the matrix.
			~AnnotationService() {}
			/// Get the first row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return text.get_interval_beg();}
			/// Get the last row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetLastIndex() const {return text.get_interval_end();}
			/// Check whether the interval of the local matrix was specified and non-empty.
			bool                   Empty() const {return text.empty();}
			/// Specify interval of the local matrix.
			/// @param beg The first row of the local matrix.
			/// @param end The last row of the local matrix.
			void                   SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end) { text.set_interval_beg(beg); text.set_interval_end(end); }
			/// Retrieve interval of the local matrix.
			/// @param start Record the first row of the local matrix into this variable.
			/// @param end Record the last row of the local matrix into this variable.
			void                   GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = text.get_interval_beg(); end = text.get_interval_end();}
			/// Retrieve the text corresponding to a certain row of the matrix.
			/// @param row Row of the matrix.
			std::string &          GetAnnotation(INMOST_DATA_ENUM_TYPE row) {assert(!text.empty()); return text[row];}
			/// Retrieve the text corresponding to a certain row of the matrix without right of modification.
			/// @param row Row of the matrix.
			const std::string &    GetAnnotation(INMOST_DATA_ENUM_TYPE row) const {assert(!text.empty()); return text[row];}
			/// Specify the text to a certain row of the matrix.
			/// @param row Row of the matrix.
			/// @param str Annotation for the row of the matrix.
			void                   SetAnnotation(INMOST_DATA_ENUM_TYPE row, std::string str) {assert(!text.empty()); text[row] = str;}
		};
#endif
#if defined(USE_SOLVER)
		/// Distributed vector class.
		/// This class can be used to store both local and distributed dense data of real type.
		/// For example, to form the right-hand side or initial guess to the solution.
		/// @see Solve
		class Vector
		{
		public:
			typedef interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> Entries; ///< Type for representation of the local array of values.
			typedef Entries::iterator iterator; ///< Type for the iterator over vector of the values.
			typedef Entries::const_iterator const_iterator; ///< Type for the iterator over vector of the constant values.
		private:
			INMOST_MPI_Comm comm; ///< Communicator corresponding to the Vector. Do we need this?
			Entries data; ///< Vector of the vector values on the current processor.
			std::string name; ///< Name of the vector, may be used to specify some parameters.
			bool is_parallel; ///< A flag indicating that the vector has extended range of values via OrderInfo class.
		public:
			/// Main constructor of the Vector class.
			/// @param _name Name of the vector, empty string by default.
			/// @param start Start of the local data interval.
			/// @param end End of the local data interval.
			/// @param _comm Communicator for parallel data exchanges, MPI_COMM_WORLD by default.
			Vector(std::string _name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
			/// Copy constructor.
			/// @param other Another vector.
			Vector(const Vector & other);
			/// Assignment operator.
			/// @param other Another vector.
			Vector & operator =(Vector const & other);
			/// Delete data of the vector.
			~Vector();
			/// Return reference to i-th element of the vector.
			INMOST_DATA_REAL_TYPE & operator [](INMOST_DATA_ENUM_TYPE i) {return data[i];}
			/// Return i-th element of the vector.
			INMOST_DATA_REAL_TYPE   operator [](INMOST_DATA_ENUM_TYPE i) const {return data[i];}
			/// Return a block of elements.
			INMOST::Matrix<INMOST_DATA_REAL_TYPE> operator [](const INMOST::AbstractMatrix<INMOST_DATA_INTEGER_TYPE> & rows) const;
#if defined(USE_AUTODIFF)
			/// Return a block of elements.
			INMOST::Matrix<value_reference> operator [](const INMOST::AbstractMatrix<INMOST_DATA_INTEGER_TYPE>& rows);
#endif //USE_AUTODIFF
			/// Return the global size of the vector.
			INMOST_DATA_ENUM_TYPE  Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			/// Iterator pointing to the first value of the vector.
			iterator             Begin() {return data.begin();}
			/// Iterator pointing to the first constant value of the vector.
			const_iterator       Begin() const {return data.begin();}
			/// Iterator pointing behind the last value of the vector.
			iterator             End() {return data.end();}
			/// Iterator pointing behind the last constant value of the vector.
			const_iterator       End() const {return data.end();}
			/// Test is there any data in the vector.
			bool                 Empty() const {return data.empty();}
			/// Set the start and the end of the distributed vector interval.
			/// @param start The first index of the local part of the vector.
			/// @param end The last index of the local part of the vector.
			void                 SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end)       {assert(start<=end); data.set_interval_beg(start); data.set_interval_end(end);}
			/// Get the start and the end of the distributed vector interval.
			/// @param start Record the first index of the local part of the vector into this variable.
			/// @param end Record the last index of the local part of the vector into this variable.
			void                 GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = data.get_interval_beg(); end = data.get_interval_end();}
			/// Move starting position of local indexes
			/// @param shift Number of positions to shift indices.
			void                 ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
			/// Get the first index of the distributed vector interval.
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
			/// Get the last index of the distributed vector interval.
			INMOST_DATA_ENUM_TYPE  GetLastIndex() const {return data.get_interval_end();}
			/// Get the communicator which the vector is associated with.
			INMOST_MPI_Comm        GetCommunicator() const {return comm;}
			/// Exchange all the data with another vector.
			/// @param other Another vector.
			void Swap(Vector & other) {data.swap(other.data); name.swap(other.name); std::swap(is_parallel,other.is_parallel); std::swap(comm,other.comm);}
			/// Save the distributed vector to a single data file using parallel MPI I/O.
			void                 Save(std::string file);
			/// Load the vector from a single data file using the specified interval.
			/// If interval is not specified, then it will be automatically constructed,
			/// with the about equal block size (the last block may has larger dimension).
			void                 Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg = ENUMUNDEF, INMOST_DATA_ENUM_TYPE mend = ENUMUNDEF, std::string file_ord = "");
			/// Test whether the vector was assigned an extended range of values via OrderInfo class.
			bool                 & isParallel() {return is_parallel;}
			/// Get the vector name specified in the main constructor.
			std::string          GetName() {return name;}
			/// Clear all data of the current vector.
			void Clear() {data.clear();}
		};
		
#endif //defined(USE_SOLVER)
#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		/// Class to store the sparse matrix row.
		/// Represents a sparse vector and facilitates access to individual entries.
		/// This class is used in expressions of automatic differentiation.
		class Row
		{
		public:
			/*
//#pragma pack(push,r1,4)
			/// Entry of the sparse matrix row.
			typedef struct entry_s
			{
				INMOST_DATA_ENUM_TYPE first;  ///< the column number of the row element.
				INMOST_DATA_REAL_TYPE second; ///< the real value of the row element.
				/// Comparison operator that helps sorting entries.
				bool operator < (const entry_s & other) const { return first < other.first || (first == other.first && second < other.second); }
				//entry_s& operator =(entry_s const& b) { first = b.first; second = b.second; return *this; }
				//entry_s(const entry_s& b) : first(b.first), second(b.second) {}
				entry_s() : first(ENUMUNDEF), second(0.0) {}
			} entry;
//#pragma pack(pop,r1)
			*/
			typedef std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> entry;
			/// Assemble an entry of entry_s type.
			/// @param ind Index.
			/// @param val Value.
			__INLINE static entry make_entry(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val)
			{
				//entry ret;
				//ret.first = ind;
				//ret.second = val;
				return std::make_pair(ind,val);
			}
		private:
			typedef std::vector<entry> Entries; ///< Container type for 
		public:
			typedef Entries::iterator iterator; ///< Iterator over pairs of index and value.
			typedef Entries::const_iterator const_iterator; ///< Iterator over constant pairs of index and value.
			typedef Entries::reverse_iterator reverse_iterator; ///< Iterator over pairs of index and value running in backward direction.
			typedef Entries::const_reverse_iterator const_reverse_iterator; ///< Iterator over constant pairs of index and value running in backward direction.
		private:
			Entries data; ///< Array of paris of index and value.
		public:
			/// Construct an empty row.
			Row() : data() {}
			/// Copy all data from another row.
			/// @param other Another row.
			Row(const Row & other) : data(other.data) {}
			/// Construct a row from array of pairs of indices and values.
			/// @param pbegin Pointer to the first position in array.
			/// @param pend Pointer behind the last position of array.
			Row(entry * pbegin, entry * pend)  :data(pbegin, pend) {}
			/// Release all data.
			~Row() {}
			/// Copy all data from another row.
			/// @param other Another row.
			Row & operator = (Row const & other) { data = other.data; return *this; }
			Row& operator = (Row && other) { data = std::move(other.data); return *this; }
			/// Finds and returns value with specified index. Adds a new entry if index was not found.
			/// \warning
			/// The operator [] should be used to fill the sparse matrix row, but not to routinely access individual elements of the row.
			/// You can use Sparse::Row::GetIndex and Sparse::Row::GetValue for rapid access to individual elements.
			/// @param i Index.
			/// @return Value corresponding to specified index.
			INMOST_DATA_REAL_TYPE & operator [](INMOST_DATA_ENUM_TYPE i)
			{
				for(Entries::size_type it = 0; it < data.size(); ++it)
					if( data[it].first == i ) return data[it].second;
				data.push_back(make_entry(i,0));
				return data.back().second;
			}
			/// Finds and returns value with specified index. Rises exception on debug and returns extra large value on release if
			/// index was not found.
			/// \warning
			/// The operator [] should be used to fill the sparse matrix row, but not to routinely access individual elements of the row.
			/// You can use Sparse::Row::GetIndex and Sparse::Row::GetValue for rapid access to individual elements.
			/// @param i Index.
			/// @return Value corresponding to specified index.
			INMOST_DATA_REAL_TYPE operator[](INMOST_DATA_ENUM_TYPE i) const
			{
				for (Entries::size_type it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;
				//you should not come here
				assert(false);
				return 1.0e20;
			}
			/// Finds and returns value with specified index.
			/// Returns zero if no entry was found.
			/// @param i Index.
			/// @return Value corresponding to specified index.
			INMOST_DATA_REAL_TYPE get_safe(INMOST_DATA_ENUM_TYPE i) const
			{
				for (Entries::size_type it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;
				return 0.0;
			}
			/// Clear all data of the current row.
			void                    Clear() { data.clear(); }
			void                    Reserve(INMOST_DATA_ENUM_TYPE size) { data.reserve(size); }
			/// Exchange all the data with another row.
			/// @param other Another row.
			void                    Swap(Row & other) { data.swap(other.data); }
			/// The size of the sparse row, i.e. the total number of nonzero elements.
			INMOST_DATA_ENUM_TYPE   Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			/// Checks are there any nonzero entries in the row.
			bool                    Empty() const { return data.empty(); }
			/// Retrieve an index corresponding to certain position in the array of pairs of index and value.
			/// @param k Position in the array of pairs of index and value.
			/// @return Index corresponding to the position in the array.
			INMOST_DATA_ENUM_TYPE & GetIndex(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->first;}
			/// Retrieve a value corresponding to certain position in the array of pairs of index and value.
			/// @param k Position in the array of pairs of index and value.
			/// @return Value corresponding to the position in the array.
			INMOST_DATA_REAL_TYPE & GetValue(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->second;}
			/// Retrieve an index corresponding to certain position in the array of pairs of index and value.
			/// @param k Position in the array of pairs of index and value.
			/// @return Index corresponding to the position in the array.
			INMOST_DATA_ENUM_TYPE   GetIndex(INMOST_DATA_ENUM_TYPE k) const {assert(k < data.size()); return (data.begin()+k)->first;}
			/// Retrieve a value corresponding to certain position in the array of pairs of index and value.
			/// @param k Position in the array of pairs of index and value.
			/// @return Value corresponding to the position in the array.
			INMOST_DATA_REAL_TYPE   GetValue(INMOST_DATA_ENUM_TYPE k) const {assert(k < data.size()); return (data.begin()+k)->second;}
			/// Retrive interval of nonzeroes
			void                    GetInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end) const;
			/// Retrive indices
			void                    GetIndices(std::vector<Sparse::bit_type>& bitset, std::vector<INMOST_DATA_ENUM_TYPE>& inds) const;
			void                    GetIndices(std::set<INMOST_DATA_ENUM_TYPE>& indset) const;
			/// Merge row indices and values with indices and values in array.
			void                    GetPairs(INMOST_DATA_REAL_TYPE coef, Sparse::Row& inds, Sparse::Row& temp) const;
			/// Merge row indices with indices in array.
			void                    GetIndices(std::vector<INMOST_DATA_ENUM_TYPE>& inds, std::vector<INMOST_DATA_ENUM_TYPE>& temp) const;
			/// Replace indices in array.
			void                    GetIndices(std::vector<INMOST_DATA_ENUM_TYPE>& inds) const;
			/// Retrive values
			void                    GetValues(INMOST_DATA_REAL_TYPE coef, const std::vector<INMOST_DATA_ENUM_TYPE>& inds, std::vector<INMOST_DATA_REAL_TYPE>& vals) const;
			void                    GetValues(INMOST_DATA_REAL_TYPE coef, const std::set<INMOST_DATA_ENUM_TYPE>& indset, std::vector<INMOST_DATA_REAL_TYPE>& vals) const;
			/// An iterator pointing to the first position in the array of pairs of index and value.
			iterator                Begin() {return data.begin();}
			/// An iterator pointing behind the last position in the array of pairs of index and value.
			iterator                End() {return data.end();}
			/// An iterator pointing to the first position in the array of constant pairs of index and value.
			const_iterator          Begin() const {return data.begin();}
			/// An iterator pointing behind the last position in the array of constant pairs of index and value.
			const_iterator          End() const {return data.end();}
			/// An iterator pointing to the last position in the array of pairs of index and value.
			reverse_iterator        rBegin() { return data.rbegin(); }
			/// An iterator pointing before the first position in the array of pairs of index and value.
			reverse_iterator        rEnd() { return data.rend(); }
			/// An iterator pointing to the last position in the array of constant pairs of index and value.
			const_reverse_iterator  rBegin() const { return data.rbegin(); }
			/// An iterator pointing before the first position in the array of constant pairs of index and value.
			const_reverse_iterator  rEnd() const { return data.rend(); }
			/// Last element
			entry&                  Back() {return data.back(); }
			/// Last element
			const entry&            Back() const { return data.back(); }
#if defined(USE_SOLVER)
			/// Return the scalar product of the current sparse row by a dense Vector.
			INMOST_DATA_REAL_TYPE   RowVec(Vector & x) const; // returns A(row) * x
#endif
			/// An optimized assignment of the row, when the content of the source row may not be preserved.
			/// @param source Source raw where to get the contents.
			void                    MoveRow(Row & source) {data = source.data;} //here move constructor and std::move may be used in future
			/// Set the vector entries by zeroes.
			void                    Zero() {for(iterator it = Begin(); it != End(); ++it) it->second = 0;}
			__INLINE void           Insert(iterator it, INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val) { data.insert(it, make_entry(ind, val)); }
			/// Push specified element into sparse row.
			/// This function should be used only if the index is not repeated in the row.
			__INLINE void           Push(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val) {data.push_back(make_entry(ind,val));}
			/// Remove last element.
			__INLINE void           Pop() { data.pop_back(); }
			/// Resize row to specified size.
			/// It is intended to be used together with non-const Row::GetIndex and Row::GetValue
			/// that allow for the modification of individual entries.
			/// @param size New size of the row.
			__INLINE void           Resize(INMOST_DATA_ENUM_TYPE size) {data.resize(size);}
			/// Output all entries of the row.
			void                    Print(double eps = -1, std::ostream & sout = std::cout) const
			{
				int k = 0;
				for(const_iterator it = Begin(); it != End(); ++it) if( fabs(it->second) > eps ) {sout << "(" << it->first << "," << it->second << ") "; k++; }
				if( k ) sout << std::endl;
			}
			/// Sort row
			__INLINE void           Sort() { std::sort(data.begin(), data.end()); }
			void                    Unique();
			/// Check whether the row is sorted.
			bool                    isSorted() const;
			/// Add up two rows. Performs operation output=alpha*left+beta*right.
			/// @param alpha Coefficient to multiply the left row.
			/// @param left The left row.
			/// @param beta Coefficient to multiply the right row.
			/// @param right The right row.
			/// @param output Record result in this vector.
			static void             MergeSortedRows(INMOST_DATA_REAL_TYPE alpha, const Row& left, INMOST_DATA_REAL_TYPE beta, const Row& right, Row& output);
			static void             MergeSortedRows(INMOST_DATA_REAL_TYPE alpha, const Row& left, INMOST_DATA_ENUM_TYPE & ileft, INMOST_DATA_REAL_TYPE beta, const Row& right, INMOST_DATA_ENUM_TYPE& iright, Row& output, INMOST_DATA_ENUM_TYPE & ind);
		};
		
#endif //defined(USE_SOLVER) || defined(USE_AUTODIFF)
		
#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		/// Class to store the compressed symmetric matrix of a hessian row.
		class HessianRow
		{
		public:
			/// Entry of the sparse matrix row.
			typedef struct hessian_index_s
			{
				INMOST_DATA_ENUM_TYPE first;
				INMOST_DATA_ENUM_TYPE second;
				bool operator < (const hessian_index_s & other) const {return first < other.first || (first == other.first && second < other.second); }
				bool operator ==(const hessian_index_s & other) const {return first == other.first && second == other.second;}
			} index;
			__INLINE static index make_index(INMOST_DATA_ENUM_TYPE _first, INMOST_DATA_ENUM_TYPE _second)
			{
				index ret;
				ret.first = _first >_second ? _second : _first;
				ret.second = _first < _second ? _second : _first;
				return ret;
			}
			typedef struct hessian_entry_s
			{
				index first;  ///< the column number of the row element.
				INMOST_DATA_REAL_TYPE second; ///< the real value of the row element.
				bool operator < (const hessian_entry_s & other) const { return first < other.first || (first == other.first && second < other.second); }
			} entry;
			__INLINE static entry make_entry(index ind, INMOST_DATA_REAL_TYPE val)
			{
				entry ret;
				ret.first = ind;
				ret.second = val;
				return ret;
			}
		private:
			typedef std::vector<entry> Entries; //replace later with more memory-efficient chunk_array, with first chunk in stack
		public:
			typedef Entries::iterator iterator;
			typedef Entries::const_iterator const_iterator;
			typedef Entries::reverse_iterator reverse_iterator;
			typedef Entries::const_reverse_iterator const_reverse_iterator;
		private:
			Entries data;
		public:
			HessianRow() : data() {}
			HessianRow(const HessianRow & other) : data(other.data) {}
			HessianRow(entry * pbegin, entry * pend) : data(pbegin,pend) {}
			~HessianRow() {}
			HessianRow & operator = (HessianRow const & other) {data = other.data; return *this;}
			/// The operator [] used to fill the sparse matrix row, but not to access individual elements of the row.
			INMOST_DATA_REAL_TYPE & operator [](index i) // use to fill matrix, not to access individual elements
			{
				for(Entries::size_type it = 0; it < data.size(); ++it)
					if( data[it].first == i ) return data[it].second;
				data.push_back(make_entry(i,0));
				return data.back().second;
			}
			/// The operator [] used to access individual elements of the row.
			INMOST_DATA_REAL_TYPE operator[](index i) const
			{
				for (Entries::size_type it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;
				//you should not come here
				assert(false);
				return 1.0e20;
			}
			INMOST_DATA_REAL_TYPE get_safe(index i) const
			{
				for (Entries::size_type it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;
				return 0.0;
			}
			//void           Reserve(INMOST_DATA_ENUM_TYPE num) { data.reserve(num);}
			/// Clear all data of the current row.
			void                    Clear() { data.clear(); }
			void                    Swap(HessianRow & other) { data.swap(other.data); }
			/// The size of the sparse row, i.e. the total number of nonzero elements.
			INMOST_DATA_ENUM_TYPE   Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			bool                    Empty() const { return data.empty(); }
			index                 & GetIndex(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->first;}
			INMOST_DATA_REAL_TYPE & GetValue(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->second;}
			index                   GetIndex(INMOST_DATA_ENUM_TYPE k) const {assert(k < data.size()); return (data.begin()+k)->first;}
			INMOST_DATA_REAL_TYPE   GetValue(INMOST_DATA_ENUM_TYPE k) const {assert(k < data.size()); return (data.begin()+k)->second;}
			iterator                Begin() {return data.begin();}
			iterator                End() {return data.end();}
			const_iterator          Begin() const {return data.begin();}
			const_iterator          End() const {return data.end();}
			reverse_iterator        rBegin() { return data.rbegin(); }
			reverse_iterator        rEnd() { return data.rend(); }
			const_reverse_iterator  rBegin() const { return data.rbegin(); }
			const_reverse_iterator  rEnd() const { return data.rend(); }
			/// Return the scalar product of the current sparse row by a dense Vector.
			void                    RowVec(INMOST_DATA_REAL_TYPE alpha, const Row & rU, INMOST_DATA_REAL_TYPE beta, Row & rJ) const; // returns A(row) * x
			void                    MoveRow(HessianRow & new_pos) {data = new_pos.data;} //here move constructor and std::move may be used in future
			/// Set the vector entries by zeroes.
			void                    Zero() {for(iterator it = Begin(); it != End(); ++it) it->second = 0;}
			/// Push specified element into sparse row.
			/// This function should be used only if the index is not repeated in the row.
			void                    Push(index ind, INMOST_DATA_REAL_TYPE val) {data.push_back(make_entry(ind,val));}
			/// Resize row to specified size.
			/// It is intended to be used together with non-const Row::GetIndex and Row::GetValue
			/// that allow for the modification of individual entries.
			/// @param size New size of the row.
			void                    Resize(INMOST_DATA_ENUM_TYPE size) {data.resize(size);}
			void                    Print() const
			{
				for(const_iterator it = Begin(); it != End(); ++it) std::cout << "(" << it->first.first << "," << it->first.second << "," << it->second << ") ";
				std::cout << std::endl;
			}
			bool                    isSorted() const;
			/// output = alpha * left + beta *right
			static void             MergeSortedRows(INMOST_DATA_REAL_TYPE alpha, const HessianRow & left, INMOST_DATA_REAL_TYPE beta, const HessianRow & right, HessianRow & output);
			static void             MergeJacobianHessian(INMOST_DATA_REAL_TYPE a, const Row & JL, const Row & JR, INMOST_DATA_REAL_TYPE b, const HessianRow & HL, INMOST_DATA_REAL_TYPE c, const HessianRow & HR, HessianRow & output);
			static void             MergeJacobianHessian(INMOST_DATA_REAL_TYPE a, const Row & JL, const Row & JR, INMOST_DATA_REAL_TYPE b, const HessianRow & H, HessianRow & output);
		};
		
#endif //defined(USE_SOLVER) || defined(USE_AUTODIFF)
		
#if defined(USE_SOLVER)
		
		/// This class can be used for shared access to matrix with OpenMP.
		class LockService
		{
			interval<INMOST_DATA_ENUM_TYPE,INMOST_OMP_LOCK_T> locks;
			void DestroyLocks();
		public:
			LockService(INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0) { if( end != start ) SetInterval(start,end); }
			LockService(const LockService & other) { SetInterval(other.GetFirstIndex(),other.GetLastIndex()); }
			LockService & operator = (LockService const & other) { SetInterval(other.GetFirstIndex(),other.GetLastIndex()); return *this;}
			~LockService() {DestroyLocks();}
			/// Get the first row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const;
			/// Get the last row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetLastIndex() const;
			bool Empty() const;
			void SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end);
			void GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const;
			bool HaveLocks() const;
			bool Lock(INMOST_DATA_ENUM_TYPE row);
			bool TestLock(INMOST_DATA_ENUM_TYPE row);
			bool UnLock(INMOST_DATA_ENUM_TYPE row);
		};
		
		/// Class to store the distributed sparse matrix by compressed rows.
		/// The format used to store sparse matrix is analogous to Compressed Row Storage format (CRS).
		/// @see http://netlib.org/linalg/html_templates/node91.html
		class Matrix
		{
			typedef interval<INMOST_DATA_ENUM_TYPE,Row> Rows;
		public:
			typedef Rows::iterator iterator;
			typedef Rows::const_iterator const_iterator;
		private:
			INMOST_MPI_Comm comm;
			Rows data;
			std::string name;
			bool is_parallel;
		public:
			/// Main constructor of the Matrix class.
			/// @param _name Name of the matrix, empty string by default.
			/// @param start Start of the local data interval.
			/// @param end End of the local data interval.
			/// @param _comm Communicator for parallel data exchanges, MPI_COMM_WORLD by default.
			Matrix(std::string _name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
			Matrix(const Matrix & other);
			Matrix & operator =(Matrix const & other);
			~Matrix();
			/// Return reference to i-th Row of the matrix.
			Row & operator [](INMOST_DATA_ENUM_TYPE i) {return data[i];}
			/// Return reference to i-th Row of the matrix.
			const Row & operator [](INMOST_DATA_ENUM_TYPE i) const {return data[i];}
			/// Return the total number of rows in the matrix.
			INMOST_DATA_ENUM_TYPE  Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			bool                   Empty() const {return data.empty();}
			iterator               Begin() {return data.begin();}
			iterator               End()   {return data.end();}
			const_iterator         Begin() const {return data.begin();}
			const_iterator         End() const   {return data.end();}
			/// Set the start and the end row numbers of the distributed matrix interval.
			void                   SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end) {data.set_interval_beg(start); data.set_interval_end(end);}
			/// Get the start and the end row numbers of the distributed matrix interval.
			void                   GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = data.get_interval_beg(); end = data.get_interval_end();}
			void                   ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
			/// Get the first row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
			/// Get the last row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetLastIndex() const {return data.get_interval_end();}
			/// Get the communicator which the matrix is associated with.
			INMOST_MPI_Comm        GetCommunicator() const {return comm;}
			void                   MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size); //for parallel
			void                   Swap(Matrix & other);
			/// Matrix-vector product of the form: y = alpha*A*x + beta * y.
			/// @param y Input/output vector.
			void MatVec(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & y) const;
			/// Matrix-vector product with transposed matrix of the form: y = alpha*A^T*x + beta * y.
			/// @param y Input/output vector.
			void MatVecTranspose(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & y) const;
			/// Clear all data of the matrix.
			void Clear() {for(Matrix::iterator it = Begin(); it != End(); ++it) it->Clear(); data.clear();}
			/// Count number of nonzeros in matrix.
			INMOST_DATA_ENUM_TYPE Nonzeros() {INMOST_DATA_ENUM_TYPE nnz = 0; for(Matrix::iterator it = Begin(); it != End(); ++it) nnz += it->Size(); return nnz;}
			/// Load the matrix from a single data file in MTX format using the specified interval.
			/// If interval is not specified, then it will be automatically constructed,
			/// with the about equal block size (the last block may has larger dimension).
			void				 Load(std::string file, INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, INMOST_DATA_ENUM_TYPE end = ENUMUNDEF, std::string file_ord = "");
			/// Save the distributed matrix to a single data file in MTX format using parallel MPI I/O.
			/// @see http://math.nist.gov/MatrixMarket/formats.html
			void                 Save(std::string file, const AnnotationService * annotation = NULL);
			/// Check that matrix is in parallel state.
			bool &               isParallel() { return is_parallel; }
			const bool &         isParallel() const { return is_parallel; }
			/// Get the matrix name specified in the main constructor.
			std::string          GetName() const {return name;}
			/// Sort rows
			void                 Sort() { for (iterator it = Begin(); it != End(); ++it) it->Sort(); }
			/// Calculate the real residual.
			///
			/// @param RHS The right-hand side Vector b.
			/// @param SOL The initial guess to the solution on input and the solution Vector x on return.
			/// @return ||A*x-b||
			///
			/// It is assumed that the coefficient matrix A have been set
			/// and the preconditioner have been already constructed.
			INMOST_DATA_REAL_TYPE Residual(Sparse::Vector &RHS, Sparse::Vector &SOL);
			/// ZAXPBY operation for sparse matrices
			/// Z = alpha * X + beta * Y
			static void ZAXPBY(INMOST_DATA_REAL_TYPE alpha, const Sparse::Matrix& X, INMOST_DATA_REAL_TYPE beta, const Sparse::Matrix& Y, Sparse::Matrix& Z);
		};
		
#endif //defined(USE_SOLVER)
		
#if defined(USE_SOLVER)
		
		/// Class to store the distributed sparse hessian hyper matrix by compressed symmetric matrices.
		class HessianMatrix
		{
		public:
			typedef interval<INMOST_DATA_ENUM_TYPE,HessianRow> HessianRows;
			typedef HessianRows::iterator iterator;
			typedef HessianRows::const_iterator const_iterator;
		private:
			INMOST_MPI_Comm comm;
			HessianRows data;
			std::string name;
			bool is_parallel;
		public:
			/// Main constructor of the Matrix class.
			/// @param _name Name of the matrix, empty string by default.
			/// @param start Start of the local data interval.
			/// @param end End of the local data interval.
			/// @param _comm Communicator for parallel data exchanges, MPI_COMM_WORLD by default.
			HessianMatrix(std::string _name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
			HessianMatrix(const HessianMatrix & other);
			HessianMatrix & operator =(HessianMatrix const & other);
			~HessianMatrix();
			/// Return reference to i-th Row of the matrix.
			HessianRow & operator [](INMOST_DATA_ENUM_TYPE i) {return data[i];}
			/// Return reference to i-th Row of the matrix.
			const HessianRow & operator [](INMOST_DATA_ENUM_TYPE i) const {return data[i];}
			/// Return the total number of rows in the matrix.
			INMOST_DATA_ENUM_TYPE  Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			bool                   Empty() const {return data.empty();}
			iterator               Begin() {return data.begin();}
			iterator               End()   {return data.end();}
			const_iterator         Begin() const {return data.begin();}
			const_iterator         End() const   {return data.end();}
			/// Set the start and the end row numbers of the distributed matrix interval.
			void                   SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end) {data.set_interval_beg(start); data.set_interval_end(end);}
			/// Get the start and the end row numbers of the distributed matrix interval.
			void                   GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = data.get_interval_beg(); end = data.get_interval_end();}
			void                   ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
			/// Get the first row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
			/// Get the last row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetLastIndex() const {return data.get_interval_end();}
			/// Get the communicator which the matrix is associated with.
			INMOST_MPI_Comm        GetCommunicator() const {return comm;}
			void                   MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size); //for parallel
			void                   Swap(HessianMatrix & other);
			/// HyperMatrix-Matrix product of the form: J = alpha*H*U + beta * J.
			/// J - jacobian sparse matrix
			/// H - hessian symmetric sparse hypermatrix
			/// U - update sparse matrix
			/// @param J Input/output Matrix.
			void MatVec(INMOST_DATA_REAL_TYPE alpha, const Matrix & U, INMOST_DATA_REAL_TYPE beta, Matrix & J) const;
			/// Clear all data of the matrix.
			void Clear() {for(HessianMatrix::iterator it = Begin(); it != End(); ++it) it->Clear(); data.clear();}
			/// Load the matrix from a single data file in MTX format using the specified interval.
			/// If interval is not specified, then it will be automatically constructed,
			/// with the about equal block size (the last block may has larger dimension).
			void				 Load(std::string file, INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, INMOST_DATA_ENUM_TYPE end = ENUMUNDEF);
			/// Save the distributed matrix to a single data file in MTX format using parallel MPI I/O.
			/// @see http://math.nist.gov/MatrixMarket/formats.html
			void                 Save(std::string file, const AnnotationService * annotation = NULL);
			/// Check that matrix is in parallel state.
			bool &               isParallel() { return is_parallel; }
			const bool &         isParallel() const { return is_parallel; }
			/// Get the matrix name specified in the main constructor.
			std::string          GetName() const {return name;}
		};
		
#endif //defined(USE_SOLVER)
		
#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		struct RowMerger2
		{
			//std::set<INMOST_DATA_ENUM_TYPE> indset;
			std::vector<Sparse::bit_type> bitset;
			std::vector<INMOST_DATA_REAL_TYPE> vals;
			std::vector<INMOST_DATA_ENUM_TYPE> inds, temp;
			RowMerger2() { bitset.resize(1048576, 0); }
			inline void radix_sort();
			inline void naive_sort();
			void clear();
			void set_bitset(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end);
			void set_vals();
			void get_row(Sparse::Row& r);
		};
		struct RowMerger3
		{
			std::vector<INMOST_DATA_REAL_TYPE> vals;
			std::vector<INMOST_DATA_ENUM_TYPE> inds;
			std::vector<INMOST_DATA_ENUM_TYPE> temp;
			void set_vals();
			void clear();
			void get_row(Sparse::Row& r);
		};
		struct RowMerger4
		{
			Sparse::Row inds;
			Sparse::Row temp;
			void clear();
			void get_row(Sparse::Row& r);
		};
		struct RowMerger5
		{
			//BinaryHeapCustom<INMOST_DATA_ENUM_TYPE, std::less<INMOST_DATA_ENUM_TYPE> > heap;
			//typedef std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> queue_t;
			//std::priority_queue< queue_t, std::vector<queue_t>, std::greater<queue_t> > heap;
			std::vector< unsigned short > list;
			std::vector<Sparse::Row> merge;
			std::vector< std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > heap;
			Sparse::Row leafs, store;
			std::vector<INMOST_DATA_ENUM_TYPE> pos;
			std::vector<INMOST_DATA_REAL_TYPE> coefs;
			std::vector<const Sparse::Row*> links;
			RowMerger5();
			void add_row(const Sparse::Row* r, INMOST_DATA_REAL_TYPE coef);
			void add_value(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val);
			void clear();
			void get_row(Sparse::Row& r);
		};
		/// This class may be used to sum multiple sparse rows.
		/// \warning
		/// In parallel column indices of the matrix may span wider then
		/// local row indices, to prevent any problem you can safely
		/// set total size of the matrix as interval of the RowMerger.
		class RowMerger
		{
#if defined(TEST_HASHTABLE)
			typedef HashTable map_container;
#else
			//typedef std::unordered_map<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> map_container;
			typedef robin_hood::unordered_map<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> map_container;
			//typedef judyLArray<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> map_container;
			//typedef tsl::hopscotch_map<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> map_container;
			//typedef tsl::bhopscotch_map<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> map_container;
#endif
		public:
			static const INMOST_DATA_ENUM_TYPE EOL = ENUMUNDEF-1; ///< End of linked list.
			//const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF; ///< Value not defined in linked list.
			class iterator
			{
			private:
				typedef typename RowMerger::map_container::iterator it_type;
				it_type it;
				RowMerger * merger; ///< Link to associated storage for linked list.
				iterator(RowMerger* pmerger) : it(pmerger->pos.begin()), merger(pmerger) {}
				iterator(INMOST_DATA_ENUM_TYPE pos, RowMerger* pmerger) : it(pmerger->pos.find(pos)), merger(pmerger) {}
				iterator(it_type pit, RowMerger* pmerger) : it(pit), merger(pmerger) {}
			public:
				iterator(const iterator & other) : it(other.it), merger(other.merger) {}
				~iterator() {}
				INMOST_DATA_REAL_TYPE& operator *() { return merger->vals[it->second]; }
				INMOST_DATA_REAL_TYPE operator *() const { return merger->vals[it->second]; }
				INMOST_DATA_REAL_TYPE* operator ->() { return &merger->vals[it->second]; }
				const INMOST_DATA_REAL_TYPE* operator ->() const { return &merger->vals[it->second]; }
				iterator & operator ++(){ it++; return *this;}
				iterator operator ++(int) { return iterator(it++, merger); }
				iterator & operator = (const iterator & other) {merger = other.merger; it = other.it; return *this;}
				bool operator ==(const iterator & other) const {return merger == other.merger && it == other.it;}
				bool operator !=(const iterator & other) const {return merger != other.merger || it != other.it;}
				friend class RowMerger;
			};
			//typedef container::iterator iterator;
		private:
			INMOST_DATA_ENUM_TYPE Nonzeros; ///< Number of nonzero in linked list.
			//INMOST_DATA_ENUM_TYPE First; ///< First position.
			//std::unordered_map<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> pos; //Position in vals and next array (huge array)
			map_container pos;  //Position in vals and next array (huge array)
			//mutable void* judy_array;
			std::vector<INMOST_DATA_REAL_TYPE> vals; //Values at the position (small array)
#if 0
			std::vector<INMOST_DATA_ENUM_TYPE> next; //Next nonzero position (small array)
			//interval< INMOST_DATA_ENUM_TYPE, Row::entry > LinkedList; ///< Storage for linked list.
#endif
			//container data;
			INMOST_DATA_ENUM_TYPE get_pos(INMOST_DATA_ENUM_TYPE pos) const;
			INMOST_DATA_ENUM_TYPE get_pos(INMOST_DATA_ENUM_TYPE pos);
			void ins_pos(INMOST_DATA_ENUM_TYPE pos, INMOST_DATA_ENUM_TYPE val);
		public:
			/// Default constructor without size specified.
			RowMerger();
			/// Constructor with size specified.
			/// @param interval_begin First index in linked list.
			/// @param interval_end Last index in linked list.
			/// @param Sorted Result should be sorted or not.
			//RowMerger(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end);
			/// Destructor.
			~RowMerger();
			/// Resize linked list for new interval.
			/// \warning
			/// All contents of linked list will be lost after resize.
			/// @param interval_begin First index in linked list.
			/// @param interval_end Last index in linked list.
			/// @param Sorted Result should be sorted or not.
			//void Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end);
			void Resize(INMOST_DATA_ENUM_TYPE size);
#if defined(USE_SOLVER)
			/// Constructor that gets sizes from the matrix, including non-local mapping.
			/// @param A Matrix to get sizes from.
			/// @param Sorted Result should be sorted.
			RowMerger(const Matrix & A);
			/// Resize linked list for new matrix, including non-local mapping.
			/// \warning
			/// All contents of linked list will be lost after resize.
			/// @param A Matrix to get sizes from.
			/// @param Sorted Result should be sorted or not.
			void Resize(const Matrix & A);
#endif //USE_SOLVER
			/// Clear linked list.
			void Clear();// { data.clear(); }
			/// Add a row with a coefficient into empty linked list.
			/// This routine should be a bit faster then RowMerger::AddRow
			/// for empty linked list. It may result in an unexpected behavior
			/// for non-empty linked list, asserts will fire in debug mode.
			/// @param coef Coefficient to multiply row values.
			/// @param r A row to be added.
			void PushRow(INMOST_DATA_REAL_TYPE coef, const Row & r);
			/// Add a row with a coefficient into non-empty linked list.
			/// Use RowMerger::PushRow for empty linked list.
			/// @param coef Coefficient to multiply row values.
			/// @param r A row to be added.
			void AddRow(INMOST_DATA_REAL_TYPE coef, const Row & r);
			/// Multiply all entries of linked list by a coefficient.
			/// @param coef A coefficient for multiplication.
			void Multiply(INMOST_DATA_REAL_TYPE coef);
			/// Place entries from linked list into row.
			/// \warning
			/// All contents of the row will be overwritten.
			/// If you want contents of the row to be added
			/// use AddRow with this row in advance.
			/// @param r A row to be filled.
			void RetrieveRow(Row & r) const;
			//INMOST_DATA_REAL_TYPE ScalarProd(RowMerger & other);
			/// Get current number of nonzeros from linked list.
			INMOST_DATA_ENUM_TYPE Size() const {return Nonzeros;}// { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			/// Check if linked list is empty.
			bool Empty() const { return Nonzeros == 0; }
			//bool Empty() const { return data.empty(); }
			/// Retrive/add an entry from/to linked list.
			/// @param pos Position in the list.
			INMOST_DATA_REAL_TYPE& operator [] (INMOST_DATA_ENUM_TYPE ipos) { return vals[get_pos(ipos)]; } 
			/// Retrive an entry from linked list.
			/// \warning
			/// Will fire an exception if there is no entry.
			/// @param pos Position in the list.
			INMOST_DATA_REAL_TYPE operator [] (INMOST_DATA_ENUM_TYPE ipos) const { return vals[get_pos(ipos)]; }
			/// Operation of the form c = alpha a + beta b
			/// \warning
			/// Linked list must be clear before operation.
			/// @param c Row c. This will be overwritten.
			/// @param alpha Multiplier for row a.
			/// @param a Row a.
			/// @param beta Multiplier for row b.
			/// @param b Row b.
			void Merge(Row & c, INMOST_DATA_REAL_TYPE alpha, const Row & a, INMOST_DATA_REAL_TYPE beta, const Row & b)
			{
				PushRow(alpha,a);
				AddRow(beta,b);
				RetrieveRow(c);
				Clear();
			}
			///Retrive iterator for the first element.
			iterator Begin() {return iterator(this);}
			//iterator Begin() { return data.begin(); }
			///Retrive iterator for the position beyond the last element.
			iterator End() {return iterator(pos.end(), this);}
			//iterator End() { return data.end(); }
		};
#endif //defined(USE_SOLVER) || defined(USE_AUTODIFF)
	} //namespace Sparse
} //namespace INMOST

#endif //INMOST_SPARSE_INCLUDED
