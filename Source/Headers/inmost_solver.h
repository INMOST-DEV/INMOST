
#ifndef INMOST_SOLVER_INCLUDED
#define INMOST_SOLVER_INCLUDED


#include "inmost_common.h"
#include "inmost_sparse.h"
//#include "solver_prototypes.hpp"


#define DEFAULT_ADDITIVE_SCHWARTZ_OVERLAP             1
#define DEFAULT_ABSOLUTE_TOLERANCE                    1.0e-5
#define DEFAULT_RELATIVE_TOLERANCE                    1.0e-12
#define DEFAULT_DIVERGENCE_TOLERANCE                  1.0e+100
#define DEFAULT_MAXIMUM_ITERATIONS                    2500
#define DEFAULT_SOLVER_GMRES_SUBSTEPS                 2
#define DEFAULT_PRECONDITIONER_DROP_TOLERANCE         0.005
#define DEFAULT_PRECONDITIONER_REUSE_TOLERANCE        0.00005
#define DEFAULT_PRECONDITIONER_FILL_LEVEL             3
#define DEFAULT_PRECONDITIONER_DDPQ_TOLERANCE         0.75
#define DEFAULT_PRECONDITIONER_REORDER_NONZEROS       1
#define DEFAULT_PRECONDITIONER_RESCALE_ITERS          6
#define DEFAULT_PRECONDITIONER_CONDITION_ESTIMATION   1
#define DEFAULT_PRECONDITIONER_ADAPT_DDPQ_TOLERANCE   1

#if defined(USE_SOLVER)
namespace INMOST
{
	/// Main class to set and solve linear system.
	/// Solver class is used to set the coefficient Matrix, the right-hand side Vector
	/// and the initial guess Vector, construct the preconditioner and Solve
	/// the linear system.
	///
	/// Formally, Solver class is independent of INMOST::Mesh class.
	/// @see Sparse::Matrix
	/// @see Sparse::Vector
	/// @see Sparse::Solve
	class Solver
	{
	private:
		static INMOST_MPI_Type RowEntryType; //prepared in Initialize
	public:
		/// Type of the Solver can be currently used in this version of INMOST.
		enum Type
		{
			INNER_ILU2,     ///< inner Solver based on BiCGStab(L) solver with second order ILU factorization as preconditioner.
			INNER_DDPQILUC, ///< inner Solver based on BiCGStab(L) solver with second order Crout-ILU with inversed-based condition estimation and unsymmetric reordering for diagonal dominance as preconditioner.
			INNER_MPTILUC,  ///< inner Solver based on BiCGStab(L) solver with second order Crout-ILU with inversed-based condition estimation and maximum product transversal reordering as preconditioner.
			INNER_MPTILU2,  ///< inner Solver based on BiCGStab(L) solver with second order ILU and maximum product transversal reordering as preconditioner.
			Trilinos_Aztec, ///< external Solver AztecOO from Trilinos package.
			Trilinos_Belos, ///< external Solver Belos from Trilinos package, currently without preconditioner.
			Trilinos_ML,    ///< external Solver AztecOO with ML preconditioner.
			Trilinos_Ifpack,///< external Solver AztecOO with Ifpack preconditioner.
			PETSc,          ///< external Solver PETSc, @see http://www.mcs.anl.gov/petsc/
			ANI,            ///< external Solver from ANI3D based on ILU2 (sequential Fortran version), @see http://ani3d.sourceforge.net/
			FCBIILU2,       ///< external FCBIILU2 Solver (BIILU2 parallel F2C version).
			K3BIILU2        ///< inner    K3BIILU2 Solver (BIILU2 parallel version).
		};

    static std::string TypeName(Type t);
		
		
		//solver.cpp::::::::::::::::::::::::::::::::::::::::::::::::::::
	public:
		
		/// Base class for low level operations with objects of Solver class.
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
			void Clear();
			/// Return true if Matrix data have already been specified.
			bool & HaveMatrix() { return have_matrix; }
			OrderInfo();
			OrderInfo(const OrderInfo & other);
			OrderInfo & operator =(OrderInfo const & other);
			~OrderInfo();
			/// Prepare parallel state of the Matrix with specified overlap size.
			/// This state of the matrix can be used, for instance, to construct
			/// the preconditioner for Additive Swartz method.
			/// @param m Matrix to be expanded.
			/// @param overlap Overlap size, viz. the number of overlap layers.
			void                    PrepareMatrix(Sparse::Matrix & m, INMOST_DATA_ENUM_TYPE overlap);
			/// Restore initial nonparallel state of the Matrix with no overlap.
			void                    RestoreMatrix(Sparse::Matrix & m);
			/// Prepare parallel state of the Vector.
			void                    PrepareVector(Sparse::Vector & v) const;
			/// Restore initial nonparallel state of the Vector.
			void                    RestoreVector(Sparse::Vector & v) const;
			/// Retrieve the processor number by binary search for the specified global index.
			INMOST_DATA_ENUM_TYPE   GetProcessor(INMOST_DATA_ENUM_TYPE gind) const; //retrieve processor by binary search in global_to_proc
			void                    GetOverlapRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const;
			/// Get the local index region for the specified process.
			void                    GetLocalRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const;
			/// Get the local index region for the current process.
			void                    GetVectorRegion(INMOST_DATA_ENUM_TYPE & mbeg, INMOST_DATA_ENUM_TYPE & mend) const {mbeg = local_vector_begin; mend = local_vector_end;}
			/// Get the rank of the current communicator, i.e. the current process index.
			INMOST_DATA_ENUM_TYPE   GetRank() const {return rank;}
			/// Get the size of the current communicator, i.e. the total number of processes used.
			INMOST_DATA_ENUM_TYPE   GetSize() const {return size;}
			/// Update the shared data in parallel vector.
			void                    Update    (Sparse::Vector & x); // update parallel vector
			/// Sum shared values in parallel vector.
			void                    Accumulate(Sparse::Vector & x); // sum shared values in parallel vector
			/// Get the sum of num elements of real array on all processes.
			void                    Integrate(INMOST_DATA_REAL_TYPE * inout, INMOST_DATA_ENUM_TYPE num) const;
			/// Get the communicator which the solver is associated with.
			INMOST_MPI_Comm         GetComm() const {return comm;}
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

			// Get the scalar product of the specified interval of the distributed vector.
			// Conflicts with OpenMP, should not be used in future
			//void ScalarProd(Vector const & left, Vector const & right, INMOST_DATA_ENUM_TYPE index_begin, INMOST_DATA_ENUM_TYPE index_end, INMOST_DATA_REAL_TYPE & sum) const;
		};

		
		
	private:
		static bool is_initialized, is_finalized;
		INMOST_MPI_Comm comm;
		std::string name;
		INMOST_DATA_ENUM_TYPE local_size, global_size;
		INMOST_DATA_ENUM_TYPE last_it;
		INMOST_DATA_REAL_TYPE last_resid;
		OrderInfo info;

		INMOST_DATA_ENUM_TYPE additive_schwartz_overlap;

		INMOST_DATA_ENUM_TYPE maximum_iterations;
		INMOST_DATA_REAL_TYPE absolute_tolerance;
		INMOST_DATA_REAL_TYPE relative_tolerance;
		INMOST_DATA_REAL_TYPE divergence_tolerance;

		INMOST_DATA_REAL_TYPE preconditioner_drop_tolerance;
		INMOST_DATA_REAL_TYPE preconditioner_reuse_tolerance;
		INMOST_DATA_REAL_TYPE preconditioner_ddpq_tolerance;
		INMOST_DATA_ENUM_TYPE preconditioner_reorder_nonzero;
		INMOST_DATA_REAL_TYPE preconditioner_fill_level;
		INMOST_DATA_ENUM_TYPE preconditioner_rescale_iterations;
		INMOST_DATA_ENUM_TYPE preconditioner_condition_estimation;
		INMOST_DATA_ENUM_TYPE preconditioner_adapt_ddpq_tolerance;

		INMOST_DATA_ENUM_TYPE solver_gmres_substeps;

		std::string return_reason;
		
		void * solver_data;
		void * matrix_data;
		void * precond_data;
		
		void * rhs_data;
		void * solution_data;
		
		Type _pack;
		Solver(const Solver & other);// prohibit copy
		Solver & operator =(Solver const & other); //prohibit assignment
	public:
    /// Retrive approximate condition number produced by INNER_MPTILUC.
    /// The number is cond(L^-1).
    INMOST_DATA_REAL_TYPE GetConditionNumberL();
    /// Retrive approximate condition number produced by INNER_MPTILUC.
    /// The number is cond(U^-1).
    INMOST_DATA_REAL_TYPE GetConditionNumberU();
		/// Set the solver parameter of the integer type.
		/// You can find defaults for parameters in the top of the file inmost_solver.h.
		///
		/// Parameters:
		/// - "maximum_iterations" - total number of iterations
		/// - "schwartz_overlap"   - number of overlapping levels for additive schwartz method,
		///                          works for: 
		///                          INNER_ILU2, INNER_MLILUC
		///                          Trilinos_Aztec, Trilinos_Belos, Trilinos_ML, Trilinos_Ifpack
		///                          PETSc
		/// - "gmres_substeps"     - number of gmres steps performed after each bicgstab step,
		///                          works for:
		///                          INNER_ILU2, INNER_MLILUC
		/// - "reorder_nonzeros"   - place sparser rows at the beggining of matrix during reordering,
		///                          works for:
		///                          INNER_MLILUC
		/// - "rescale_iterations" - number of iterations for two-side matrix rescaling,
		///                          works for:
		///                          INNER_ILU2, INNER_MLILUC
		/// - "condition_estimation" - exploit condition estimation of inversed factors to adapt
		///                          drop and reuse tolerances,
		///                          works for:
		///                          INNER_MLILUC
		/// - "adapt_ddpq_tolerance" - adapt ddpq tolerance depending from the complexity 
		///                          of calculation of Schur complement,
		///                          works for:
		///                          INNER_MLILUC
		void SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value);
		/// Set the solver parameter of the real type.
		/// You can find defaults for parameters in the top of the file inmost_solver.h.
		///
		/// Parameters:
		/// - "absolute_tolerance" - iterative method will stop on i-th iteration
		///                          if ||A x(i)-b|| < absolute_tolerance
		/// - "relative_tolerance" - iterative method will stop on i-th iteration
		///                          if ||A x(i)-b||/||A x(0) - b||
		/// - "divergence_tolerance" - iterative method will fail if
		///                          ||A x(i) - b|| > divergence_tolerance
		/// - "drop_tolerance"     - tolerance for dropping values during incomplete factorization,
		///                          works for:
		///                          INNER_ILU2, INNER_MLILUC
		///                          Trilinos_Aztec, Trilinos_Ifpack
		///                          PETSc
		/// - "reuse_tolerance"    - tolerance for reusing values during incomplete factorization,
		///                          these values are used only during calculation of L and U factors
		///                          and/or Schur complement and discarded once factorization is done,
		///                          value should be less then "drop_tolerance",
		///                          typical value is drop_tolerance^2,
		///                          works for:
		///                          INNER_ILU2, INNER_MLILUC
		/// - "ddpq_tolerance"     - by this tolerance most diagonnaly-dominant elements will be selected
		///                          to form the next level of factorization, the closer the tolerance
		///                          is to one the smaller will be the level. Actual rule is:
		///                          A(i,j)/(sum(A(i,:))+sum(A(:,j))-A(i,j)) > ddpq_tolerance *
		///                          A(imax,jmax)/(sum(A(imax,:))+sum(A(:,jmax))-A(imax,jmax))
		///                          where on imax, jmax maximum is reached.
		///                          works for:
		///                          INNER_MLILUC
		/// - "fill_level"         - level of fill for ILU-type preconditioners,
		///                          works for:
		///                          INNER_ILU2 (if LFILL is defined in solver_ilu2.hpp)
		///                          Trilinos, Trilinos_Ifpack
		void SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value);
		/// Get the used defined name of the Solver.
		std::string          GetName() {return name;}
		/// Get the package Type.
		Type GetPackage() const {return _pack;}
		/// Set the matrix and construct the preconditioner.
		/// @param A Matrix A in linear problem Ax = b
    /// @param ModifiedPattern Indicates whether the structure of the matrix have 
    /// changed since last call to Solver::SetMatrix.
		/// @param OldPreconditioner If this parameter is set to true,
		/// then the previous preconditioner will be used,
		/// otherwise the new preconditioner will be constructed. 
		///
		/// Preconditioner will be constructed on call to this function
		/// - for INNER_*, PETSc and ANI packages
		/// - for Trilinos preconditioner will be constructed each time Sparse::Solve is called
		///
		/// Any changes to preconditioner parameters should happen before that point.
		/// If you increase gmres_substep after this point, inner methods most likely will fail
		void SetMatrix(Sparse::Matrix & A, bool ModifiedPattern = true, bool OldPreconditioner = false);
		/// Solver the linear system: A*x = b.
		/// Prior to this call you should call SetMatrix
		///
		/// @param RHS The right-hand side Vector b.
		/// @param SOL The initial guess to the solution on input and the solution Vector x on return.
		///
		/// It is assumed that the coefficient matrix A have been set
		/// and the preconditioner have been already constructed.
		///
		/// @see Sparse::SetMatrix
		bool Solve(Sparse::Vector & RHS, Sparse::Vector & SOL);
		/// Return the number of iterations performed by the last solution.
		/// @see Sparse::Solve
		INMOST_DATA_ENUM_TYPE Iterations();
		/// Return the final residual achieved by the last solution.
		/// @see Sparse::Solve
		INMOST_DATA_REAL_TYPE Residual();
		/// Get the reason of convergence or divergence of the last solution.
		/// @see Sparse::Solve
		std::string GetReason();
    /// Computes the smallest and the largest eigenvalue with the power method.
    /// Requires SetMatrix to be called to compute the preconditioner.
    /// Currently works for internal methods only, since it uses internal matrix-vector multiplication.
    /// Largest eigenvalue: vprev = 0; v = rand(); while( |v|-|vprev| > tol ) {vprev = v; v = A*v; v /= |v|;}
    ///                     lambda_max = |v|;
    /// Smallest eigenvalue: vprev = 0; v = rand(); while( |v|-|vprev| > tol ){vprev = v; solve(A*v = v); v /= |v|;}
    ///                     lambda_min = 1.0/|v|;
    /// See answer by Blair Perot in:
    /// https://www.researchgate.net/post/What_is_the_best_way_to_estimate_the_condition_number_of_a_sparse_matrix.
    /// @param tol Tolerance used for power series.
    /// @param maxits Maximum number of iterations allowed.
    /// @return Condition number or 1.0e100 if not converged.
    INMOST_DATA_REAL_TYPE Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxits = 100);
		/// Main constructor of the solver.
		/// Solver name provided here is used to extract options from database file
		/// for PETSc and Trilinos packages.
		/// @param pack The package Type to be used for solution.
		/// @param _name The user specified name of the current solver.
		/// @param comm Communicator for parallel data exchanges, MPI_COMM_WORLD by default.
		/// @see Solver::Initialize
		/// @see Solver::SetMatrix
		/// @see Solver::Solve
		/// @see Solver::Finalize
		Solver(Type pack, std::string _name = "", INMOST_MPI_Comm comm = INMOST_MPI_COMM_WORLD);
		~Solver();
		/// Initialize the stage of parallel solution.
		/// If MPI is not initialized yet, then it will be initialized.
		///
		/// database file is used to pass parameters to PETSc and Trilinos packages.
		/// if database file for is provided any changes through SetParameterEnum,
		/// SetParameterReal would not be effective for PETSc and Trilinos packages.
		/// Currently this database file provides directions for package-specific
		/// files. In future it is supposed to set up parameters for internal solvers.
		/// @param argc The number of arguments transmitted to the function main.
		/// @param argv The pointer to arguments transmitted to the function main.
		/// @param database Usually the name of the file with the Solver parameters.
		///
		/// The shortest call to this function with the default solver parameters is the following: Initialize(NULL,NULL,"");
		/// @see Solver::Finalize
		/// @see Solver::isInitialized
		///
		/// Example of contents of the database file:
		///
		/// 	PETSc: petsc_options.txt
		/// 	Trilinos_Ifpack: trilinos_ifpack_options.xml
		/// 	Trilinos_ML: trilinos_ml_options.xml
		/// 	Trilinos_Aztec: trilinos_aztec_options.xml
		/// 	Trilinos_Belos: trilinos_belos_options.xml
		static void Initialize(int * argc, char *** argv, const char * database = "");
		/// Finalize the stage of parallel solution.
		/// If MPI was initialized in Solver::Initialize, then it will be finalized.
		/// By this reason, do not use any MPI function after call to this function.
		/// @see Solver::Initialize
		/// @see Solver::isFinalized
		static void Finalize();
		static bool isInitialized() {return is_initialized;}
		static bool isFinalized() {return is_finalized;}
		/// Clear all internal data of the current solver including matrix, preconditioner etc.
		void Clear();
	};
}

#endif // USE_SOLVER

#endif // INMOST_SOLVER_INCLUDED
