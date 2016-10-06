#ifndef INMOST_SPARSE_INCLUDED
#define INMOST_SPARSE_INCLUDED


#include "inmost_common.h"





namespace INMOST
{
    namespace Sparse
    {
#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		/// Retrive MPI type for row entry type
		INMOST_MPI_Type GetRowEntryType();
		/// Create MPI type for row entry type
		void CreateRowEntryType();
		/// Release MPI type
		void DestroyRowEntryType();
		/// Check whether MPI type was created
		bool HaveRowEntryType();
		/// This class can be used to annotate the matrix
		class AnnotationService
		{
			interval<INMOST_DATA_ENUM_TYPE,std::string> text;
		public:
			AnnotationService(INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0) { if( end != start ) SetInterval(start,end); }
			AnnotationService(const AnnotationService & other) : text(other.text) { }
			AnnotationService & operator = (AnnotationService const & other) {text = other.text; return *this;}
			~AnnotationService() {}
			/// Get the first row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return text.get_interval_beg();}
			/// Get the last row index of the distributed matrix interval.
			INMOST_DATA_ENUM_TYPE  GetLastIndex() const {return text.get_interval_end();}
			bool                   Empty() const {return text.empty();}
			void                   SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end) { text.set_interval_beg(beg); text.set_interval_end(end); }
			void                   GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = text.get_interval_beg(); end = text.get_interval_end();}
			std::string &          GetAnnotation(INMOST_DATA_ENUM_TYPE row) {assert(!text.empty()); return text[row];}
			const std::string &    GetAnnotation(INMOST_DATA_ENUM_TYPE row) const {assert(!text.empty()); return text[row];}
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
            /// Move starting position of local indexes
            void                 ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
            /// Get the first index of the distributed vector interval.
            INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
            /// Get the last index of the distributed vector interval.
            INMOST_DATA_ENUM_TYPE  GetLastIndex() const {return data.get_interval_end();}
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
        
#endif //defined(USE_SOLVER)
		
#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		
        /// Class to store the sparse matrix row.
        class Row
        {
        public:
            /// Entry of the sparse matrix row.
            typedef struct entry_s
            {
                INMOST_DATA_ENUM_TYPE first;  ///< the column number of the row element.
                INMOST_DATA_REAL_TYPE second; ///< the real value of the row element.
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
            //typedef dynarray<entry,16> Entries; //replace later with more memory-efficient chunk_array, with first chunk in stack
            typedef array<entry> Entries;
        public:
            typedef Entries::iterator iterator;
            typedef Entries::const_iterator const_iterator;
            typedef Entries::reverse_iterator reverse_iterator;
            typedef Entries::const_reverse_iterator const_reverse_iterator;
        private:
            Entries data;
        public:
            Row() : data() {}
            Row(const Row & other) : data(other.data) {}
            Row(entry * pbegin, entry * pend)  :data(pbegin, pend) {}
            ~Row() {}
            Row & operator = (Row const & other) { data = other.data; return *this; }
            /// The operator [] used to fill the sparse matrix row, but not to access individual elements of the row.
            INMOST_DATA_REAL_TYPE & operator [](INMOST_DATA_ENUM_TYPE i) // use to fill matrix, not to access individual elements
            {
                for(Entries::size_type it = 0; it < data.size(); ++it)
                    if( data[it].first == i ) return data[it].second;
                data.push_back(make_entry(i,0));
                return data.back().second;
            }
            /// The operator [] used to access individual elements of the row.
            INMOST_DATA_REAL_TYPE operator[](INMOST_DATA_ENUM_TYPE i) const
            {
                for (Entries::size_type it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;
                //you should not come here
                assert(false);
                return 1.0e20;
            }
            /// Returns zero if no entry was found.
            INMOST_DATA_REAL_TYPE get_safe(INMOST_DATA_ENUM_TYPE i) const
            {
                for (Entries::size_type it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;
                return 0.0;
            }
            //void           Reserve(INMOST_DATA_ENUM_TYPE num) { data.reserve(num);}
            /// Clear all data of the current row.
            void                    Clear() { data.clear(); }
            void                    Swap(Row & other) { data.swap(other.data); }
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
#if defined(USE_SOLVER)
            /// Return the scalar product of the current sparse row by a dense Vector.
            INMOST_DATA_REAL_TYPE   RowVec(Vector & x) const; // returns A(row) * x
#endif
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
            void                    Print() const
            {
                for(const_iterator it = Begin(); it != End(); ++it) std::cout << "(" << it->first << "," << it->second << ") ";
                std::cout << std::endl;
            }
            bool                    isSorted() const;
            /// output = alpha * left + beta *right
            static void             MergeSortedRows(INMOST_DATA_REAL_TYPE alpha, const Row & left, INMOST_DATA_REAL_TYPE beta, const Row & right, Row & output);
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
                ret.first = std::min(_first,_second);
                ret.second = std::max(_first,_second);
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
            typedef array<entry> Entries; //replace later with more memory-efficient chunk_array, with first chunk in stack
            //typedef dynarray<entry,8> Entries;
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
            void                    Print()
            {
                for(iterator it = Begin(); it != End(); ++it) std::cout << "(" << it->first.first << "," << it->first.second << "," << it->second << ") ";
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
                iterator & operator = (const iterator & other) {LinkedList = other.LinkedList; pos = other.pos; return *this;}
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
			INMOST_DATA_ENUM_TYPE IntervalBeg; ///< Begin of global interval of owned index interval
			INMOST_DATA_ENUM_TYPE IntervalEnd;  ///< End of global interval of owned index interval
            interval< INMOST_DATA_ENUM_TYPE, Row::entry > LinkedList; ///< Storage for linked list.
			std::vector< INMOST_DATA_ENUM_TYPE > NonlocalPre; ///< List of global indices, that are to the left of owned index interval
			std::vector< INMOST_DATA_ENUM_TYPE > NonlocalPost; ///< List of global indices, that are to the right of owned index interval
        public:
			/// This function converts global index into local index.
			/// @param pos Global index.
			/// @return Local index.
			INMOST_DATA_ENUM_TYPE MapIndex(INMOST_DATA_ENUM_TYPE pos) const;
			/// This function converts local index into global index.
			/// @param pos Local index.
			/// @return Global index.
			INMOST_DATA_ENUM_TYPE UnmapIndex(INMOST_DATA_ENUM_TYPE pos) const;
			/// This function provides information about additional non-local indices.
			/// \warning
			/// All contents of linked list will be lost.
			/// @param Pre Non-local indices that go before IntervalBegin.
			/// @param Post Non-local indices that follow IntervalEnd.
			void SetNonlocal(const std::vector<INMOST_DATA_ENUM_TYPE> & Pre, const std::vector<INMOST_DATA_ENUM_TYPE> & Post);
            /// Default constructor without size specified.
            RowMerger();
            /// Constructor with size specified.
            /// @param interval_begin First index in linked list.
            /// @param interval_end Last index in linked list.
            /// @param Sorted Result should be sorted or not.
            RowMerger(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool Sorted = true);
			/// Constructor with size and non-local mapping specified.
			/// @param interval_begin First index in linked list.
			/// @param interval_end Last index in linked list.
			/// @param Pre Nonlocal indices before First index in linked list.
			/// @param Post Nonlocal indices after Last index in linked list.
			/// @param Sorted Result should be sorted or not.
			RowMerger(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, const std::vector<INMOST_DATA_ENUM_TYPE> & Pre, const std::vector<INMOST_DATA_ENUM_TYPE> & Post, bool Sorted = true);
			/// Destructor.
            ~RowMerger();
            /// Resize linked list for new interval.
            /// \warning
            /// All contents of linked list will be lost after resize.
            /// @param interval_begin First index in linked list.
            /// @param interval_end Last index in linked list.
            /// @param Sorted Result should be sorted or not.
            void Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool Sorted = true);
			/// Resize linked list for new interval with non-local mapping.
			/// \warning
			/// All contents of linked list will be lost after resize.
			/// @param interval_begin First index in linked list.
			/// @param interval_end Last index in linked list.
			/// @param Pre Nonlocal indices before First index in linked list.
			/// @param Post Nonlocal indices after Last index in linked list.
			/// @param Sorted Result should be sorted or not.
			void Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, const std::vector<INMOST_DATA_ENUM_TYPE> & Pre, const std::vector<INMOST_DATA_ENUM_TYPE> & Post, bool Sorted = true);
#if defined(USE_SOLVER)
			/// Constructor that gets sizes from the matrix, including non-local mapping.
			/// @param A Matrix to get sizes from.
			/// @param Sorted Result should be sorted.
			RowMerger(const Matrix & A, bool Sorted = true);
			/// Resize linked list for new matrix, including non-local mapping.
            /// \warning
            /// All contents of linked list will be lost after resize.
            /// @param A Matrix to get sizes from.
            /// @param Sorted Result should be sorted or not.
            void Resize(const Matrix & A, bool Sorted = true);
#endif //USE_SOLVER
            /// Clear linked list.
            void Clear();
            /// Add a row with a coefficient into empty linked list.
            /// This routine should be a bit faster then RowMerger::AddRow
            /// for empty linked list. It may result in an unexpected behavior
            /// for non-empty linked list, asserts will fire in debug mode.
            /// @param coef Coefficient to multiply row values.
            /// @param r A row to be added.
            /// @param PreSortRow Sort values of the row before adding. Will be activated only for sorted linked lists.
            void PushRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow);
            /// Add a row with a coefficient into non-empty linked list.
            /// Use RowMerger::PushRow for empty linked list.
            /// @param coef Coefficient to multiply row values.
            /// @param r A row to be added.
            /// @param PreSortRow Sort values of the row before adding. Will be activated only for sorted linked lists.
            void AddRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow);
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
            void Merge(Row & c, INMOST_DATA_REAL_TYPE alpha, const Row & a, INMOST_DATA_REAL_TYPE beta, const Row & b)
            {
                PushRow(alpha,a);
                AddRow(beta,b);
                RetriveRow(c);
                Clear();
            }
			///Retrive iterator for the first element.
            iterator Begin() {return iterator(&LinkedList);}
			///Retrive iterator for the position beyond the last element.
            iterator End() {iterator ret(&LinkedList); ret.pos = EOL; return ret;}
        };
#endif //defined(USE_SOLVER) || defined(USE_AUTODIFF)
    } //namespace Sparse
} //namespace INMOST

#endif //INMOST_SPARSE_INCLUDED