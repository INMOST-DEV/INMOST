#ifndef INMOST_SPARSE_INCLUDED
#define INMOST_SPARSE_INCLUDED


#include "inmost_common.h"


#if defined(USE_SOLVER)

#define MTX_ALLOW_ANNOTATION
#define MTX_ANNOTATION_SIZE 1024

namespace INMOST
{
    namespace Sparse
    {
      INMOST_MPI_Type GetRowEntryType();
      void CreateRowEntryType();
      void DestroyRowEntryType();
      bool HaveRowEntryType();

      /// Distributed vector class.
		  /// This class can be used to store both local and distributed dense data of real type.
		  /// For example, to form the right-hand side or initial guess to the solution.
		  /// @see Solve
		  class Vector
		  {
		  public:
			  typedef interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> Entries;
			  typedef Entries::iterator iterator;
			  typedef Entries::const_iterator const_iterator;
		  private:
			  INMOST_MPI_Comm comm;
			  Entries data;
			  std::string name;
			  bool is_parallel;
		  public:
			  /// Main constructor of the Vector class.
			  /// @param _name Name of the vector, empty string by default.
			  /// @param start Start of the local data interval.
			  /// @param end End of the local data interval.
			  /// @param _comm Communicator for parallel data exchanges, MPI_COMM_WORLD by default.
			  Vector(std::string _name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
			  Vector(const Vector & other);
			  Vector & operator =(Vector const & other);
			  ~Vector();
			  /// Return reference to i-th element of the vector.
			  INMOST_DATA_REAL_TYPE & operator [](INMOST_DATA_ENUM_TYPE i) {return data[i];}
			  /// Return i-th element of the vector.
			  INMOST_DATA_REAL_TYPE   operator [](INMOST_DATA_ENUM_TYPE i) const {return data[i];}
			  /// Return the global size of the vector.
			  INMOST_DATA_ENUM_TYPE  Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			  iterator             Begin() {return data.begin();}
			  const_iterator       Begin() const {return data.begin();}
			  iterator             End() {return data.end();}
			  const_iterator       End() const {return data.end();}
			  bool                 Empty() const {return data.empty();}
			  /// Set the start and the end of the distributed vector interval.
			  void                 SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end)       {assert(start<=end); data.set_interval_beg(start); data.set_interval_end(end);}
			  /// Get the start and the end of the distributed vector interval.
			  void                 GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = data.get_interval_beg(); end = data.get_interval_end();}
			  void                 ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
			  /// Get the first index of the distributed vector interval.
			  INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
			  /// Get the communicator which the vector is associated with.
			  INMOST_MPI_Comm        GetCommunicator() const {return comm;}

			  void Swap(Vector & other) {data.swap(other.data); name.swap(other.name); std::swap(is_parallel,other.is_parallel); std::swap(comm,other.comm);}


			  /// Save the distributed vector to a single data file using parallel MPI I/O.
			  void                 Save(std::string file);
			  /// Load the vector from a single data file using the specified interval.
			  /// If interval is not specified, then it will be automatically constructed,
			  /// with the about equal block size (the last block may has larger dimension).
			  void                 Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg = ENUMUNDEF, INMOST_DATA_ENUM_TYPE mend = ENUMUNDEF);
			
			  bool                 & isParallel() {return is_parallel;}
			  /// Get the vector name specified in the main constructor.
			  std::string          GetName() {return name;}

			  /// Clear all data of the current vector.
			  void Clear() {data.clear();}
			
		  };



		  /// Class to store the sparse matrix row.
		  class Row 
		  {
		  public:
  #if defined(MTX_ALLOW_ANNOTATION)
			  char annotation[MTX_ANNOTATION_SIZE];
  #endif
			  /// Entry of the sparse matrix row.
			  typedef struct entry_s 
			  {
				  INMOST_DATA_ENUM_TYPE first;  ///< the column number of the row element.
				  INMOST_DATA_REAL_TYPE second; ///< the real value of the row element.
				  //entry_s() :first(0), second(0.0) {}
				  //entry_s(const entry_s & other) :first(other.first), second(other.second) {}//{std::cout << __FUNCTION__ << " " << first << " " << second << std::endl;}
				  //entry_s(INMOST_DATA_ENUM_TYPE first, INMOST_DATA_REAL_TYPE second):first(first),second(second){}
				  //entry_s & operator =(entry_s const & other) {first = other.first, second = other.second; return *this;}
				  bool operator < (const entry_s & other) const { return first < other.first || (first == other.first && second < other.second); }
			  } entry;
			  __INLINE static entry make_entry(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val)
			  {
				  entry ret;
				  ret.first = ind;
				  ret.second = val;
				  return ret;
			  }

		  private:
			  typedef dynarray<entry,16> Entries; //replace later with more memory-efficient chunk_array, with first chunk in stack
			  //typedef array<entry> Entries;
			  //typedef std::vector<entry> Entries;
			  //typedef sparse_data<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> Entries;
			  //typedef Entries::pair entry; //for sparse_data
			
		  public:
			  typedef Entries::iterator iterator;
			  typedef Entries::const_iterator const_iterator;
			  typedef Entries::reverse_iterator reverse_iterator;
			  typedef Entries::const_reverse_iterator const_reverse_iterator;

			  bool modified_pattern; //remove this in future
		  private:
  #if defined(USE_OMP)
			  omp_lock_t lock;
  #endif
			  bool marker;
			  Entries data;
		  public:

        std::string GetAnnotation();
        void SetAnnotation(std::string input);
			  void Report() {data.report_addr();}
			  void SetMarker() { marker = true; }
			  void RemMarker() { marker = false; }
			  bool GetMarker() { return marker; }
			  Row();
			  Row(const Row & other);
			  Row(entry * pbegin, entry * pend);
        bool HaveLock()
        {
  #if defined(USE_OMP)
				  return true;
  #else
          return false;
  #endif
        }
			  void Lock()
			  {
  #if defined(USE_OMP)
				  omp_set_lock(&lock); 
  #endif
			  }
			  void Unlock() 
			  { 
  #if defined(USE_OMP)
				  omp_unset_lock(&lock); 
  #endif
			  }
			  ~Row();
			  Row & operator = (Row const & other);
			  /// The operator [] used to fill the sparse matrix row, but not to access individual elements of the row.
			  INMOST_DATA_REAL_TYPE & operator [](INMOST_DATA_ENUM_TYPE i) // use to fill matrix, not to access individual elements
			  {
				  //for sparse_data type
				  //return data[i];
				  //for dynarray or array
				
				  for(Entries::size_type it = 0; it < data.size(); ++it)
					  if( data[it].first == i ) return data[it].second;
				  entry new_entry;
				  new_entry.first = i;
				  new_entry.second = 0;
				  data.push_back(new_entry);
				  modified_pattern = true;
				  return data.back().second;
				
			  }
			  /// The operator [] used to access individual elements of the row.
			  INMOST_DATA_REAL_TYPE operator[](INMOST_DATA_ENUM_TYPE i) const
			  {
				  //for sparse data type
				  //return data[i];

				  for (Entries::size_type it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;

				  //you should not come here
				  assert(false);
				  return 1.0e20;
			  }
			  //void           Reserve(INMOST_DATA_ENUM_TYPE num) { data.reserve(num);}
			  /// Clear all data of the current row.
			  void                    Clear() { data.clear(); }
			  void                    Swap(Row & other);
			  /// The size of the sparse row, i.e. the total number of nonzero elements.
			  INMOST_DATA_ENUM_TYPE   Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
        bool                    Empty() const { return data.empty(); }
			  INMOST_DATA_ENUM_TYPE & GetIndex(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->first;}
			  INMOST_DATA_REAL_TYPE & GetValue(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->second;}
			  INMOST_DATA_ENUM_TYPE   GetIndex(INMOST_DATA_ENUM_TYPE k) const {assert(k < data.size()); return (data.begin()+k)->first;}
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
			  INMOST_DATA_REAL_TYPE   RowVec(Vector & x) const; // returns A(row) * x
			  void                    MoveRow(Row & new_pos) {data = new_pos.data;} //here move constructor and std::move may be used in future
			  /// Set the vector entries by zeroes.
			  void                    Zero() {for(iterator it = Begin(); it != End(); ++it) it->second = 0;}
			  /// Push specified element into sparse row.
			  /// This function should be used only if the index is not repeated in the row.
			  void                    Push(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val) {data.push_back(make_entry(ind,val));}
			  /// Resize row to specified size. 
			  /// It is intended to be used together with non-const Row::GetIndex and Row::GetValue
			  /// that allow for the modification of individual entries.
			  /// @param size New size of the row.
			  void                    Resize(INMOST_DATA_ENUM_TYPE size) {data.resize(size);}

			  void                    Print() 
			  {
				  for(iterator it = Begin(); it != End(); ++it) std::cout << "(" << it->first << "," << it->second << ") ";
				  std::cout << std::endl;
			  }
		  };


		
		
		  /// Class to store the distributed sparse matrix by compressed rows.
		  /// The format used to store sparse matrix is analogous to Compressed Row Storage format (CRS).
		  /// @see http://netlib.org/linalg/html_templates/node91.html
		  class Matrix
		  {
		  public:
			  typedef interval<INMOST_DATA_ENUM_TYPE,Row> Rows;
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
			  bool                 Empty() const {return data.empty();}
			  iterator             Begin() {return data.begin();}
			  iterator             End()   {return data.end();}
			  const_iterator       Begin() const {return data.begin();}
			  const_iterator       End() const   {return data.end();}
			  /// Set the start and the end row numbers of the distributed matrix interval.
			  void                 SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end)       {data.set_interval_beg(start); data.set_interval_end(end);}
			  /// Get the start and the end row numbers of the distributed matrix interval.
			  void                 GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = data.get_interval_beg(); end = data.get_interval_end();}
			  void                 ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
			  /// Get the first row index of the distributed matrix interval.
			  INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
			  /// Get the communicator which the matrix is associated with.
			  INMOST_MPI_Comm        GetCommunicator() const {return comm;}
			  void                 MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size); //for parallel
			  void                 Swap(Matrix & other) 
			  {
				  data.swap(other.data);
				  name.swap(other.name);
				  INMOST_MPI_Comm ctmp = comm;
				  comm = other.comm;
				  other.comm = ctmp;
				  bool ptmp = is_parallel;
				  is_parallel = other.is_parallel;
				  other.is_parallel = ptmp;
			  }
			  /// Matrix-vector product of the form: y = alpha*A*x + beta * y.
			  /// @param y Input/output vector.
			  void MatVec(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & y) const;
        /// Matrix-vector product with transposed matrix of the form: y = alpha*A^T*x + beta * y.
			  /// @param y Input/output vector.
			  void MatVecTranspose(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & y) const;
        /// Clear all data of the matrix.
			  void Clear() {for(Matrix::iterator it = Begin(); it != End(); ++it) it->Clear(); data.clear();}
        /// Load the matrix from a single data file in MTX format using the specified interval.
			  /// If interval is not specified, then it will be automatically constructed,
			  /// with the about equal block size (the last block may has larger dimension).
			  void				 Load(std::string file, INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, INMOST_DATA_ENUM_TYPE end = ENUMUNDEF);
			  /// Save the distributed matrix to a single data file in MTX format using parallel MPI I/O.
			  /// @see http://math.nist.gov/MatrixMarket/formats.html
			  void                 Save(std::string file);
        /// Check that matrix is in parallel state.
			  bool &               isParallel() { return is_parallel; }
			  /// Get the matrix name specified in the main constructor.
			  std::string          GetName() {return name;}
		  };

		  /// This class may be used to sum multiple sparse rows.
		  /// \warning
		  /// In parallel column indices of the matrix may span wider then 
		  /// local row indices, to prevent any problem you are currently
		  /// advised to set total size of the matrix as interval of the
		  /// RowMerger. In future this may change, see todo 2 below.
		  /// \todo
		  /// 1. (testing!) Add iterators over entries.
		  /// 2. Implement multiple intervals for distributed computation,
		  ///    then in parallel the user may specify additional range of indexes
		  ///    for elements that lay on the borders between each pair of processors.
      ///    Or even better implement mapping that will remap nonlocal entries to the 
      ///    end of current linked list when added and put them back in correct places
      ///    when retrived. May use algorithm from class OrderInfo.
		  class RowMerger
		  {
		  public:
			  static const INMOST_DATA_ENUM_TYPE EOL = ENUMUNDEF-1; ///< End of linked list.
			  static const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF; ///< Value not defined in linked list.
        class iterator
        {
        private:
          INMOST_DATA_ENUM_TYPE pos;
          interval< INMOST_DATA_ENUM_TYPE, Row::entry > * LinkedList; ///< Link to associated storage for linked list.
          iterator(interval< INMOST_DATA_ENUM_TYPE, Row::entry > * pLinkedList) : pos(pLinkedList->get_interval_beg()), LinkedList(pLinkedList) {}
        public:
          iterator(const iterator & other) : pos(other.pos), LinkedList(other.LinkedList) {}
          ~iterator() {}
          INMOST_DATA_REAL_TYPE & operator *() {return (*LinkedList)[pos].second;}
          INMOST_DATA_REAL_TYPE operator *() const {return (*LinkedList)[pos].second;}
          INMOST_DATA_REAL_TYPE * operator ->() {return &(*LinkedList)[pos].second;}
          const INMOST_DATA_REAL_TYPE * operator ->() const {return &(*LinkedList)[pos].second;}
          iterator & operator ++(){ pos = (*LinkedList)[pos].first; return *this;}
          iterator operator ++(int){ iterator ret(LinkedList); ret.pos = (*LinkedList)[pos].first; return ret; }
          iterator & operator = (const iterator & other) {LinkedList = other.LinkedList; pos = other.pos;}
          bool operator ==(const iterator & other) const {return LinkedList == other.LinkedList && pos == other.pos;}
          bool operator !=(const iterator & other) const {return LinkedList != other.LinkedList || pos != other.pos;}
          bool operator < (const iterator & other) const {return LinkedList == other.LinkedList && pos < other.pos;}
          bool operator <=(const iterator & other) const {return LinkedList == other.LinkedList && pos <= other.pos;}
          bool operator > (const iterator & other) const {return LinkedList == other.LinkedList && pos > other.pos;}
          bool operator >=(const iterator & other) const {return LinkedList == other.LinkedList && pos >= other.pos;}
          friend class RowMerger;
        };
		  private:
			  bool Sorted; ///< Contents of linked list should be sorted.
			  INMOST_DATA_ENUM_TYPE Nonzeros; ///< Number of nonzero in linked list.
			  interval< INMOST_DATA_ENUM_TYPE, Row::entry > LinkedList; ///< Storage for linked list.
		  public:
			  /// Default constructor without size specified.
			  RowMerger();
			  /// Constructor with size specified.
			  /// @param interval_begin First index in linked list.
			  /// @param interval_end Last index in linked list.
			  /// @param Sorted Result should be sorted or not.
			  RowMerger(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool Sorted = true);
			  /// Constructor that gets sizes from the matrix.
			  /// @param A Matrix to get sizes from.
			  /// @param Sorted Result should be sorted.
			  RowMerger(Matrix & A, bool Sorted = true);
			  /// Destructor.
			  ~RowMerger();
			  /// Resize linked list for new interval.
			  /// \warning
			  /// All contents of linked list will be lost after resize.
			  /// This behavior may be changed in future.
			  /// @param interval_begin First index in linked list.
			  /// @param interval_end Last index in linked list.
			  /// @param Sorted Result should be sorted or not.
			  void Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool Sorted = true);
			  /// Resize linked list for new matrix.
			  /// \warning
			  /// All contents of linked list will be lost after resize.
			  /// This behavior may be changed in future.
			  /// @param A Matrix to get sizes from.
			  /// @param Sorted Result should be sorted or not.
			  void Resize(Matrix & A, bool Sorted = true);
			  /// Clear linked list.
			  void Clear();
			  /// Add a row with a coefficient into empty linked list.
			  /// This routine should be a bit faster then RowMerger::AddRow
			  /// for empty linked list. It may result in an unexpected behavior
			  /// for non-empty linked list, asserts will fire in debug mode.
			  /// @param coef Coefficient to multiply row values.
			  /// @param r A row to be added.
			  /// @param PreSortRow Sort values of the row before adding. Will be activated only for sorted linked lists.
			  void PushRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow = false);
			  /// Add a row with a coefficient into non-empty linked list.
			  /// Use RowMerger::PushRow for empty linked list.
			  /// @param coef Coefficient to multiply row values.
			  /// @param r A row to be added.
			  /// @param PreSortRow Sort values of the row before adding. Will be activated only for sorted linked lists.
			  void AddRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow = false);
			  /// Multiply all entries of linked list by a coefficient.
			  /// @param coef A coefficient for multiplication.
			  void Multiply(INMOST_DATA_REAL_TYPE coef);
			  /// Place entries from linked list into row.
			  /// \warning
			  /// All contents of the row will be overwritten.
			  /// If you want contents of the row to be added
			  /// use AddRow with this row in advance.
			  /// @param r A row to be filled.
			  void RetriveRow(Row & r);
			  //INMOST_DATA_REAL_TYPE ScalarProd(RowMerger & other);
			  /// Get current number of nonzeros from linked list.
			  INMOST_DATA_ENUM_TYPE Size() {return Nonzeros;}
        /// Retrive/add an entry from/to linked list.
        /// @param pos Position in the list.
        INMOST_DATA_REAL_TYPE & operator [] (INMOST_DATA_ENUM_TYPE pos);
        /// Retrive an entry from linked list.
        /// \warning
        /// Will fire an exception if there is no entry.
        /// @param pos Position in the list.
        INMOST_DATA_REAL_TYPE operator [] (INMOST_DATA_ENUM_TYPE pos) const;
        /// Operation of the form c = alpha a + beta b
        /// \warning
        /// Linked list must be clear before operation.
        /// @param c Row c. This will be overwritten.
        /// @param alpha Multiplier for row a.
        /// @param a Row a.
        /// @param beta Multiplier for row b.
        /// @param b Row b.
        void Merge(Row & c, INMOST_DATA_REAL_TYPE alpha, Row & a, INMOST_DATA_REAL_TYPE beta, Row & b)
        {
          PushRow(alpha,a);
          AddRow(beta,b);
          RetriveRow(c);
          Clear();
        }
        iterator Begin() {return iterator(&LinkedList);}
        iterator End() {iterator ret(&LinkedList); ret.pos = EOL; return ret;}
		  };
   }
}
#endif

#endif