#pragma once
#ifndef INMOST_SOLVER_INCLUDED
#define INMOST_SOLVER_INCLUDED


#include "inmost_common.h"
//#include "solver_prototypes.hpp"

#if defined(USE_SOLVER)
namespace INMOST
{
	class Solver
	{
	private:
		static INMOST_MPI_Type RowEntryType; //prepared in Initialize
	public:
		enum Type
		{
			INNER_ILU2,
			INNER_MLILUC,
			PETSC,
			ANI //sequential, fortran
		};

		static INMOST_MPI_Type & GetRowEntryType() {return RowEntryType;}
		
		//solver.cpp::::::::::::::::::::::::::::::::::::::::::::::::::::
	public:
		class Matrix;
		class Vector;

		class OrderInfo
		{
		private:
			typedef std::vector<INMOST_DATA_ENUM_TYPE> storage_type;
			storage_type global_to_proc; //stores ends of all non-overlapping intervals of elements, owned by this processor
			storage_type global_overlap; //stores pairs: [begin,end) of overlapping intervals of rows
			std::vector<INMOST_DATA_ENUM_TYPE> vector_exchange_recv, vector_exchange_send;
			std::vector<INMOST_DATA_REAL_TYPE> send_storage, recv_storage;
			std::vector<INMOST_MPI_Request> send_requests, recv_requests;
			std::vector<INMOST_DATA_ENUM_TYPE> extended_indexes;

			//remote indexes
			INMOST_DATA_ENUM_TYPE local_vector_begin, local_vector_end;
			INMOST_DATA_ENUM_TYPE initial_matrix_begin, initial_matrix_end; //local interval of matrix
			INMOST_DATA_ENUM_TYPE local_matrix_begin, local_matrix_end; //local interval of matrix

			bool have_matrix;
			INMOST_MPI_Comm comm;
			int rank,size;
		public:
			bool & HaveMatrix() { return have_matrix; }
			OrderInfo();
			OrderInfo(const OrderInfo & other);
			OrderInfo & operator =(OrderInfo const & other);
			~OrderInfo();
			void PrepareMatrix(Matrix & m, INMOST_DATA_ENUM_TYPE overlap);
			void RestoreMatrix(Matrix & m);
			void PrepareVector(Vector & v);
			void RestoreVector(Vector & v); //Restore initial nonparallel state of the vector
			INMOST_DATA_ENUM_TYPE GetProcessor(INMOST_DATA_ENUM_TYPE gind) const; //retrive processor by binary search in global_to_proc
			void      GetOverlapRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const;
			void        GetLocalRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const;
			void       GetVectorRegion(INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const {mbeg = local_vector_begin; mend = local_vector_end;}
			INMOST_DATA_ENUM_TYPE GetRank() const {return rank;}
			INMOST_DATA_ENUM_TYPE GetSize() const {return size;}
			void                  Update    (Vector & x); // update parallel vector
			void                  Accumulate(Vector & x); // sum shared values in parallel vector
			void                  Integrate(INMOST_DATA_REAL_TYPE * inout, INMOST_DATA_ENUM_TYPE num);
			INMOST_MPI_Comm         GetComm() {return comm;}
			// Access to arrays below allows to organize manual exchange
			INMOST_MPI_Request    * GetSendRequests() {assert(!send_requests.empty()); return &send_requests[0];}
			INMOST_MPI_Request    * GetRecvRequests() {assert(!recv_requests.empty()); return &recv_requests[0];}
			INMOST_DATA_ENUM_TYPE   GetSendRequestsSize() {return static_cast<INMOST_DATA_ENUM_TYPE>(send_requests.size());}
			INMOST_DATA_ENUM_TYPE   GetRecvRequestsSize() {return static_cast<INMOST_DATA_ENUM_TYPE>(recv_requests.size());}
			INMOST_DATA_ENUM_TYPE * GetSendExchangeArray() {assert(!vector_exchange_send.empty()); return &vector_exchange_send[0];}
			INMOST_DATA_ENUM_TYPE   GetSendExchangeSize() {return static_cast<INMOST_DATA_ENUM_TYPE>(send_storage.size());}
			INMOST_DATA_ENUM_TYPE * GetRecvExchangeArray() {assert(!vector_exchange_recv.empty()); return &vector_exchange_recv[0];}
			INMOST_DATA_ENUM_TYPE   GetRecvExchangeSize() {return static_cast<INMOST_DATA_ENUM_TYPE>(recv_storage.size());}
			//for debug
			//~ void                 BeginSequentialCode() {for(int i = 0; i < rank; i++) MPI_Barrier(comm);}
			//~ void                   EndSequentialCode() {for(int i = rank; i < size; i++) MPI_Barrier(comm);}
		};

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
			Vector(std::string _name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
			Vector(const Vector & other);
			Vector & operator =(Vector const & other);
			~Vector();
			INMOST_DATA_REAL_TYPE & operator [](INMOST_DATA_ENUM_TYPE i) {return data[i];}
			INMOST_DATA_REAL_TYPE   operator [](INMOST_DATA_ENUM_TYPE i) const {return data[i];}
			INMOST_DATA_ENUM_TYPE  Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			iterator             Begin() {return data.begin();}
			const_iterator       Begin() const {return data.begin();}
			iterator             End() {return data.end();}
			const_iterator       End() const {return data.end();}
			bool                 Empty() const {return data.empty();}
			void                 SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end)       {data.set_interval_beg(start); data.set_interval_end(end);}
			void                 GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = data.get_interval_beg(); end = data.get_interval_end();}
			void                 ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
			INMOST_MPI_Comm        GetCommunicator() const {return comm;}

			INMOST_DATA_REAL_TYPE ScalarProd(Vector const & other, INMOST_DATA_ENUM_TYPE index_begin, INMOST_DATA_ENUM_TYPE index_end) const;


			void                 Save(std::string file);
			void                 Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg = ENUMUNDEF, INMOST_DATA_ENUM_TYPE mend = ENUMUNDEF);
			
			bool                 & isParallel() {return is_parallel;}
			std::string          GetName() {return name;}
			//~ friend class Solver;
		};



		class Row 
		{
		public:
			
			typedef struct entry_s 
			{
				INMOST_DATA_ENUM_TYPE first; 
				INMOST_DATA_REAL_TYPE second;
				entry_s() :first(0), second(0.0) {}
				entry_s(const entry_s & other) :first(other.first), second(other.second) {}
				entry_s(INMOST_DATA_ENUM_TYPE first, INMOST_DATA_REAL_TYPE second):first(first),second(second){}
				entry_s & operator =(entry_s const & other) {first = other.first, second = other.second; return *this;}
				bool operator < (const entry_s & other) const { return first < other.first || (first == other.first && second < other.second); }
			} entry;
			
			typedef dynarray<entry,16> Entries; //replace later with more memory-efficient chunk_array, with first chunk in stack
			//typedef array<entry> Entries;
			//typedef std::vector<entry> Entries;
			//typedef sparse_data<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> Entries;
			//typedef Entries::pair entry; //for sparse_data
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
			void SetMarker() { marker = true; }
			void RemMarker() { marker = false; }
			bool GetMarker() { return marker; }
			Row() :data() 
			{
#if defined(USE_OMP)
				omp_init_lock(&lock);
#endif
				modified_pattern = marker = false;
			}
			Row(const Row & other) :marker(other.marker),data(other.data) 
			{ 
#if defined(USE_OMP)
				omp_init_lock(&lock);
#endif
				modified_pattern = other.modified_pattern; 
			}
			Row(entry * pbegin, entry * pend) :data(pbegin, pend) 
			{ 
#if defined(USE_OMP)
				omp_init_lock(&lock);
#endif
				modified_pattern = true; marker = false; 
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
			~Row() {}
			Row & operator = (Row const & other) { data = other.data; marker = other.marker; return *this; }
			INMOST_DATA_REAL_TYPE & operator [](INMOST_DATA_ENUM_TYPE i) // use to fill matrix, not to access individual elements
			{
				//for sparse_data type
				//return data[i];
				//for dynarray or array
				
				for(unsigned it = 0; it < data.size(); ++it)
					if( data[it].first == i ) return data[it].second;
				entry new_entry;
				new_entry.first = i;
				new_entry.second = 0;
				data.push_back(new_entry);
				modified_pattern = true;
				return data.back().second;
				
			}
			INMOST_DATA_REAL_TYPE operator[](INMOST_DATA_ENUM_TYPE i) const
			{
				//for sparse data type
				//return data[i];

				for (unsigned it = 0; it < data.size(); ++it) if (data[it].first == i) return data[it].second;

				//you should not come here
				assert(false);
				return 1.0e20;
			}
			//void           Reserve(INMOST_DATA_ENUM_TYPE num) { data.reserve(num);}
			void                  Clear() { data.clear(); }
			void                  Swap(Solver::Row & other) { data.swap(other.data); bool tmp = marker; marker = other.marker; other.marker = tmp; }
			INMOST_DATA_ENUM_TYPE   Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			INMOST_DATA_ENUM_TYPE & GetIndex(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->first;}
			INMOST_DATA_REAL_TYPE & GetValue(INMOST_DATA_ENUM_TYPE k) {assert(k < data.size()); return (data.begin()+k)->second;}
			INMOST_DATA_ENUM_TYPE   GetIndex(INMOST_DATA_ENUM_TYPE k) const {assert(k < data.size()); return (data.begin()+k)->first;}
			INMOST_DATA_REAL_TYPE   GetValue(INMOST_DATA_ENUM_TYPE k) const {assert(k < data.size()); return (data.begin()+k)->second;}
			
			iterator              Begin() {return data.begin();}
			iterator              End() {return data.end();}
			const_iterator        Begin() const {return data.begin();}
			const_iterator        End() const {return data.end();}
			reverse_iterator      rBegin() { return data.rbegin(); }
			reverse_iterator      rEnd() { return data.rend(); }
			const_reverse_iterator rBegin() const { return data.rbegin(); }
			const_reverse_iterator rEnd() const { return data.rend(); }
			INMOST_DATA_REAL_TYPE   RowVec(Vector & x) const; // returns A(row) * x
			void                  MoveRow(Row & new_pos) {data = new_pos.data;} //here move constructor and std::move may be used in future
			void                  Zero() {for(iterator it = Begin(); it != End(); ++it) it->second = 0;}
			//! Use this function if you know that this index is not repeated in the row (and indexes are inserted in increasing order (for sparse_data))
			void                  Push(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val) {data.push_back(entry(ind,val));}
		};
		
		class Matrix
		{
		public:
			typedef interval<INMOST_DATA_ENUM_TYPE,Solver::Row> Rows;
			typedef Rows::iterator iterator;
			typedef Rows::const_iterator const_iterator;
		private:
			INMOST_MPI_Comm comm;
			Rows data;
			std::string name;
			bool is_parallel;
		public:
			Matrix(std::string _name = "", INMOST_DATA_ENUM_TYPE start = 0, INMOST_DATA_ENUM_TYPE end = 0, INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
			Matrix(const Matrix & other);
			Matrix & operator =(Matrix const & other);
			~Matrix();
			
			Row & operator [](INMOST_DATA_ENUM_TYPE i) {return data[i];}
			const Row & operator [](INMOST_DATA_ENUM_TYPE i) const {return data[i];}
			INMOST_DATA_ENUM_TYPE  Size() const { return static_cast<INMOST_DATA_ENUM_TYPE>(data.size()); }
			bool                 Empty() const {return data.empty();}
			iterator             Begin() {return data.begin();}
			iterator             End()   {return data.end();}
			const_iterator       Begin() const {return data.begin();}
			const_iterator       End() const   {return data.end();}
			void                 SetInterval(INMOST_DATA_ENUM_TYPE   start, INMOST_DATA_ENUM_TYPE   end)       {data.set_interval_beg(start); data.set_interval_end(end);}
			void                 GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = data.get_interval_beg(); end = data.get_interval_end();}
			void                 ShiftInterval(INMOST_DATA_ENUM_TYPE shift) {data.shift_interval(shift);}
			INMOST_DATA_ENUM_TYPE  GetFirstIndex() const {return data.get_interval_beg();}
			INMOST_MPI_Comm        GetCommunicator() const {return comm;}
			void                 MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size); //for parallel
			void                 Swap(Solver::Matrix & other) 
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
			
			void MatVec(INMOST_DATA_REAL_TYPE alpha, Solver::Vector & x, INMOST_DATA_REAL_TYPE beta, Solver::Vector & out) const; //y = alpha*A*x + beta * y


			void				 Load(std::string file, INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, INMOST_DATA_ENUM_TYPE end = ENUMUNDEF);
			void                 Save(std::string file);
			bool &               isParallel() { return is_parallel; }
			
			std::string          GetName() {return name;}
			//~ friend class Solver;
		};
		
	private:
		INMOST_MPI_Comm comm;
		std::string name;
		INMOST_DATA_ENUM_TYPE local_size, global_size;
		INMOST_DATA_ENUM_TYPE last_it;
		INMOST_DATA_REAL_TYPE last_resid;
		OrderInfo info;
		INMOST_DATA_ENUM_TYPE overlap;
		
		void * solver_data;
		void * matrix_data;
		void * precond_data;
		
		void * rhs_data;
		void * solution_data;
		
		Type _pack;
		Solver(const Solver & other);// prohibit copy
		Solver & operator =(Solver const & other); //prohibit assignment
	public:
		void SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value);
		void SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value);
		std::string          GetName() {return name;}
	
		Type GetPackage() const {return _pack;}
		
		INMOST_DATA_ENUM_TYPE Iterations();
		INMOST_DATA_REAL_TYPE Residual();
		void SetMatrix(Matrix & A, bool OldPreconditioner = false);
		bool Solve(Vector & RHS, Vector & SOL);


		Solver(Type pack, std::string _name = "", INMOST_MPI_Comm comm = INMOST_MPI_COMM_WORLD);
		~Solver();
		static void Initialize(int * argc, char *** argv, const char * database = "");
		static void Finalize();
	};
}

#endif // USE_SOLVER

#endif // INMOST_SOLVER_INCLUDED
