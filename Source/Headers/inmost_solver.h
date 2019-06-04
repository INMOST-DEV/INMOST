#ifndef INMOST_SOLVER_INCLUDED
#define INMOST_SOLVER_INCLUDED

#include "inmost_common.h"
#include "inmost_sparse.h"

#if defined(USE_SOLVER)
namespace INMOST
{

	class SolverInterface;
	class SolverParameters;

	enum SolverVerbosityLevel {
		Level0 = 0,
		Level1 = 1,
		Level2 = 2
	};

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
	public:
		//Backward-compatibility
		typedef std::string Type;
		static const Type INNER_ILU2;     ///< inner Solver based on BiCGStab(L) solver with second order ILU factorization as preconditioner.
		static const Type INNER_DDPQILUC; ///< inner Solver based on BiCGStab(L) solver with second order Crout-ILU with inversed-based condition estimation and unsymmetric reordering for diagonal dominance as preconditioner.
		static const Type INNER_MPTILUC;  ///< inner Solver based on BiCGStab(L) solver with second order Crout-ILU with inversed-based condition estimation and maximum product transversal reordering as preconditioner.
		static const Type INNER_MLMPTILUC;  ///< inner Solver based on BiCGStab(L) solver with second order Crout-ILU with inversed-based condition estimation and maximum product transversal reordering as preconditioner.
		static const Type INNER_MPTILU2;  ///< inner Solver based on BiCGStab(L) solver with second order ILU and maximum product transversal reordering as preconditioner.
		static const Type Trilinos_Aztec; ///< external Solver AztecOO from Trilinos package.
		static const Type Trilinos_Belos; ///< external Solver Belos from Trilinos package, currently without preconditioner.
		static const Type Trilinos_ML;    ///< external Solver AztecOO with ML preconditioner.
		static const Type Trilinos_Ifpack;///< external Solver AztecOO with Ifpack preconditioner.
		static const Type PETSc;          ///< external Solver PETSc, @see http://www.mcs.anl.gov/petsc/
		static const Type ANI;            ///< external Solver from ANI3D based on ILU2 (sequential Fortran version), @see http://ani3d.sourceforge.net/
		static const Type FCBIILU2;       ///< external FCBIILU2 Solver (BIILU2 parallel F2C version).
		static const Type K3BIILU2;       ///< inner    K3BIILU2 Solver (BIILU2 parallel version).
		static const Type SUPERLU;        ///< external Solver SuperLU @see https://github.com/starseeker/SuperLU
		/// Backward-compatibility access to integer-typed parameters.
		/// Check Solver::SetParameter for availible parameter names.
		/// @param name Name of the parameter.
		/// @param value Value for the parameter.
		/// @see Solver::SetParameter
		void SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value);
		/// Backward-compatibility access to real-typed parameters.
		/// Check Solver::SetParameter for availible parameter names.
		/// @param name Name of the parameter.
		/// @param value Value for the parameter.
		/// @see Solver::SetParameter
		void SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value);
		/// Backward-compatibility acces to return reason.
		std::string GetReason() {return ReturnReason();}
		//Backward-compatibility
	private:
		static std::vector<SolverParameters> parameters; ///< A global array of solver parameters.
		static int *argc; ///< A global number of command line arguments.
		static char ***argv; ///< A global array of strings corresponding to command line arguments.
		static bool is_initialized; ///< Indicator that solvers were initialized for MPI operations.
		static bool is_finalized; ///< Indicator that solvers were finalized for MPI operations.
		SolverInterface *solver; ///< Pointer to abstract linear solver.
		std::string prefix; ///< Prescribed solver name that is used to assign solver parameters.
		SolverVerbosityLevel versbosity;
		double preconditioner_time;
		double iterations_time;
		bool is_solved;
		/// Reads the parameters from database file stored in xml format.
		/// @param xml_database Path to xml file with solver parameters.
		static void parseXMLDatabase(const char* xml_database);
	public:
		/// Base class for low level parallel operations with objects of Solver class.
		/// This class auguments the matrix with additional layers of overlap for
		/// Additive Schwartz method to work. Also allows to expand vectors
		/// according to new matrix size and to organize vector data exchange in
		/// both directions.
		class OrderInfo
		{
		private:
			typedef std::vector<INMOST_DATA_ENUM_TYPE> storage_type; ///< A type for storing arrays.
			storage_type global_to_proc; ///< Stores ends of all non-overlapping intervals of elements, owned by this processor.
			storage_type global_overlap; ///< Stores pairs: [begin,end) of overlapping intervals of rows
			std::vector<INMOST_DATA_ENUM_TYPE> vector_exchange_recv; ///< Store packed indices to recieve data. Format: # of processors {proc #, # indices, indices}
			std::vector<INMOST_DATA_ENUM_TYPE> vector_exchange_send; ///< Store packed indices to send data. Format: # of processors {proc #, # indices, indices}
			std::vector<INMOST_DATA_REAL_TYPE> send_storage; ///< Storage used to send data of the vector.
			std::vector<INMOST_DATA_REAL_TYPE> recv_storage; ///< Storage used to receive data of the vector.
			std::vector<INMOST_MPI_Request> send_requests; ///< Sturctures used to wait complition of send operations.
			std::vector<INMOST_MPI_Request> recv_requests; ///< Sturctures used to detect complition of receive operations.
			std::vector<INMOST_DATA_ENUM_TYPE> extended_indexes; ///< All the indices that appear outside of the local
			INMOST_DATA_ENUM_TYPE local_vector_begin; ///< Remember the first index in expanded vector.
			INMOST_DATA_ENUM_TYPE local_vector_end; ///< Remember the last index in expanded vector.
			INMOST_DATA_ENUM_TYPE initial_matrix_begin; ///< Initial the first index of the matrix before overlapping was computed.
			INMOST_DATA_ENUM_TYPE initial_matrix_end; ///< Initial the last index of the matrix before overlapping was computed.
			INMOST_DATA_ENUM_TYPE local_matrix_begin; ///< The first index of the matrix with overlapping.
			INMOST_DATA_ENUM_TYPE local_matrix_end; ///< The last index of the matrix with overlapping.
			bool have_matrix; ///< Indicates whether the matrix was specified and all structures were initialized.
			INMOST_MPI_Comm comm; ///< Parallel communicator.
			int rank; ///< Index of the processor in communicator.
			int size; ///< Total number of processors in communicator.
		public:
			/// Clear all structures and data.
			void Clear();
			/// Return true if Matrix data have already been specified.
			bool &HaveMatrix() { return have_matrix; }
			/// Default initialize all structures.
			OrderInfo();
			/// Copy all structures.
			/// @param other Right hand side.
			OrderInfo(const OrderInfo &other);
			/// Assign all structures.
			/// @param other Right hand side.
			OrderInfo &operator=(OrderInfo const &other);
			/// Clears all data.
			~OrderInfo();
			/// Prepare parallel state of the Matrix with specified overlap size.
			/// This state of the matrix can be used, for instance, to construct
			/// the preconditioner for Additive Swartz method.
			/// @param m Matrix to be expanded.
			/// @param overlap Overlap size, viz. the number of overlap layers.
			void PrepareMatrix(Sparse::Matrix &m, INMOST_DATA_ENUM_TYPE overlap);
			/// Restore initial nonparallel state of the Matrix with no overlap.
			/// @param m Matrix to be restored.
			void RestoreMatrix(Sparse::Matrix &m);
			/// Prepare parallel state of the Vector.
			/// @param v Vector to be expanded.
			void PrepareVector(Sparse::Vector &v) const;
			/// Restore initial nonparallel state of the Vector.
			/// @param v Vector to be restored.
			void RestoreVector(Sparse::Vector &v) const;
			/// Retrieve the processor number by binary search for the specified global index.
			/// Finds the processor that contains provided index in it's interval of indices.
			/// @param gind Global index in the matrix.
			INMOST_DATA_ENUM_TYPE GetProcessor(INMOST_DATA_ENUM_TYPE gind) const;
			/// Get the overlap index region for the specified processor.
			/// @param proc Processor number.
			/// @param mbeg Record the first index here.
			/// @param mbeg Record the last index here.
			void GetOverlapRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE &mbeg, INMOST_DATA_ENUM_TYPE &mend) const;
			/// Get the local index region for the specified processor.
			/// @param proc Processor number.
			/// @param mbeg Record the first index here.
			/// @param mbeg Record the last index here.
			void GetLocalRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE &mbeg, INMOST_DATA_ENUM_TYPE &mend) const;
			/// Get the local index region for the current processor.
			/// @param mbeg Record the first index here.
			/// @param mbeg Record the last index here.
			void GetVectorRegion(INMOST_DATA_ENUM_TYPE &mbeg, INMOST_DATA_ENUM_TYPE &mend) const { mbeg = local_vector_begin; mend = local_vector_end; }
			/// Get the rank of the current communicator, i.e. the current process index.
			INMOST_DATA_ENUM_TYPE GetRank() const { return rank; }
			/// Get the size of the current communicator, i.e. the total number of processes used.
			INMOST_DATA_ENUM_TYPE GetSize() const { return size; }
			/// Update the shared data in parallel vector.
			/// @param x Vector for which the shared data should be updated.
			void Update(Sparse::Vector &x);
			/// Sum shared values in parallel vector.
			/// @param x Vector for which the shared data should be accumulated.
			void Accumulate(Sparse::Vector &x);
			/// Get the sum of num elements of real array on all processes.
			/// @param inout Data that should be integrated.
			/// @param num Number of entries in inout array.
			void Integrate(INMOST_DATA_REAL_TYPE *inout, INMOST_DATA_ENUM_TYPE num) const;
			/// Get the communicator which the solver is associated with.
			INMOST_MPI_Comm GetComm() const { return comm; }
			/// MPI structures that hold information on sent data.
			/// Access to MPI structures that allow for asynchronous exchange of data.
			INMOST_MPI_Request *GetSendRequests() { assert(!send_requests.empty()); return &send_requests[0]; }
			/// MPI structures that hold information on posted receive requests.
			/// Access to MPI structures that allow for asynchronous exchange of data.
			INMOST_MPI_Request *GetRecvRequests() { assert(!recv_requests.empty()); return &recv_requests[0]; }
			/// Total number of send requests.
			INMOST_DATA_ENUM_TYPE GetSendRequestsSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(send_requests.size()); }
			/// Total number of posted recieve requests.
			INMOST_DATA_ENUM_TYPE GetRecvRequestsSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(recv_requests.size()); }
			/// Request raw access to array used to send data.
			INMOST_DATA_ENUM_TYPE *GetSendExchangeArray() { assert(!vector_exchange_send.empty()); return &vector_exchange_send[0]; }
			/// Size of the array used to send data.
			INMOST_DATA_ENUM_TYPE GetSendExchangeSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(send_storage.size()); }
			/// Request raw acces to array used to recieve data.
			INMOST_DATA_ENUM_TYPE *GetRecvExchangeArray() { assert(!vector_exchange_recv.empty()); return &vector_exchange_recv[0]; }
			/// Size of the array used to receive data.
			INMOST_DATA_ENUM_TYPE GetRecvExchangeSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(recv_storage.size()); }
			/// Place this function before the code that should be executed sequatially by each processor.
			/// Used for debugging purposes in parallel.
			//void BeginSequentialCode() {for(int i = 0; i < rank; i++) MPI_Barrier(comm);}
			/// Place this function after the code that should be executed sequatially by each processor.
			/// Used for debugging purposes in parallel.
			//void EndSequentialCode() {for(int i = rank; i < size; i++) MPI_Barrier(comm);}
		};
		/// Main constructor of the solver.
		/// @param solverName The solver name to be used for solution.
		/// @param prefix The user specified name of the current solver.
		/// @param comm Communicator for parallel data exchanges, MPI_COMM_WORLD by default.
		/// @see Solver::Initialize
		/// @see Solver::SetMatrix
		/// @see Solver::Solve
		/// @see Solver::Finalize
		Solver(std::string solverName, std::string prefix = "", INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
		/// Copy a solver.
		/// \warning Not all solvers support assignment. This operation may be very expensive.
		Solver(const Solver &other);
		/// Assign a solver.
		/// \warning Not all solvers support assignment. This operation may be very expensive.
		Solver &operator=(const Solver &other);
		/// Return the solver name
		/// @see Sparse::Solve
		std::string SolverName() const;
		/// Return the solver user specified name of the current solver
		/// @see Sparse::Solve
		std::string SolverPrefix() const;
		/// Return the solver verbosity specified level of the current solver
		/// @see Sparse::Solve
		SolverVerbosityLevel VerbosityLevel() const;
		///
		void SetVerbosityLevel(SolverVerbosityLevel level);
		/// Initialize the stage of parallel solution.
		/// If MPI is not initialized yet, then it will be initialized.
		///
		/// database file is used to pass parameters to Inner solvers, PETSc and Trilinos packages.
		/// if database file for is provided any changes through SetParameter
		/// would not be effective for PETSc and Trilinos packages.
		/// @param argc The number of arguments transmitted to the function main.
		/// @param argv The pointer to arguments transmitted to the function main.
		/// @param database Usually the name of the file with the Solver parameters.
		///
		/// The shortest call to this function with the default solver parameters is the following: Initialize(NULL,NULL,"");
		/// @see Solver::Finalize
		/// @see Solver::isInitialized
		///
		/// Example of contents of the database file:
		///     Main: database.xml
		///     PETSc: petsc_options.txt
		///     Trilinos_Ifpack: trilinos_ifpack_options.xml
		///     Trilinos_ML: trilinos_ml_options.xml
		///     Trilinos_Aztec: trilinos_aztec_options.xml
		///     Trilinos_Belos: trilinos_belos_options.xml
		static void Initialize(int *argc, char ***argv, const char *database = NULL);
		/// Finalize the stage of parallel solution.
		/// If MPI was initialized in Solver::Initialize, then it will be finalized.
		/// By this reason, do not use any MPI function after call to this function.
		/// @see Solver::Initialize
		/// @see Solver::isFinalized
		static void Finalize();
		/// Checks the stage of parallel solution is initialized
		static bool isInitialized();
		/// Checks the stage of parallel solution is finalized
		static bool isFinalized();
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
		void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern = true, bool OldPreconditioner = false);
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
		bool Solve(Sparse::Vector &RHS, Sparse::Vector &SOL);
		/// Clear all internal data of the current solver including matrix, preconditioner etc.
		bool Clear();
		/// Get the solver output parameter
		/// @param name The name of solver's output parameter
		/// @see Solver::SetParameter
		std::string GetParameter(std::string name) const;
		/// @param name The name of parameter
		/// @param value The value of parameter
		/// Set the solver parameter of the integer type.
		///
		/// Parameters:
		/// - "maximum_iterations" - total number of iterations
		/// - "schwartz_overlap"   - number of overlapping levels for additive schwartz method,
		///                          works for:
		///                          INNER_ILU2, INNER_MPTILU2, INNER_MPTILUC, INNER_MLMPTILUC, INNER_DDPQILUC
		///                          Trilinos_Aztec, Trilinos_Belos, Trilinos_ML, Trilinos_Ifpack
		///                          PETSc
		/// - "gmres_substeps"     - number of gmres steps performed after each bicgstab step,
		///                          works for:
		///                          INNER_ILU2, INNER_MPTILU2, INNER_MPTILUC, INNER_MLMPTILUC, INNER_DDPQILUC
		/// - "reorder_nonzeros"   - place sparser rows at the beggining of matrix during reordering,
		///                          works for:
		///                          INNER_DDPQILUC
		/// - "rescale_iterations" - number of iterations for two-side matrix rescaling,
		///                          works for:
		///                          INNER_ILU2, INNER_MPTILU2, INNER_MPTILUC, INNER_MLMPTILUC, INNER_DDPQILUC
		/// - "condition_estimation" - exploit condition estimation of inversed factors to adapt
		///                          drop and reuse tolerances,
		///                          works for:
		///                          INNER_MPTILUC, INNER_MLMPTILUC, INNER_DDPQILUC
		/// - "adapt_ddpq_tolerance" - adapt ddpq tolerance depending from the complexity
		///                          of calculation of Schur complement,
		///                          works for:
		///                          INNER_DDPQILUC
		/// Set the solver parameter of the real type.
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
		///                          INNER_ILU2, INNER_MPTILU2, INNER_MPTILUC, INNER_MLMPTILUC, INNER_DDPQILUC
		///                          Trilinos_Aztec, Trilinos_Ifpack
		///                          PETSc
		/// - "reuse_tolerance"    - tolerance for reusing values during incomplete factorization,
		///                          these values are used only during calculation of L and U factors
		///                          and/or Schur complement and discarded once factorization is done,
		///                          value should be less then "drop_tolerance",
		///                          typical value is drop_tolerance^2,
		///                          works for:
		///                          INNER_ILU2, INNER_MPTILU2, INNER_MPTILUC, INNER_MLMPTILUC, INNER_DDPQILUC
		/// - "ddpq_tolerance"     - by this tolerance most diagonnaly-dominant elements will be selected
		///                          to form the next level of factorization, the closer the tolerance
		///                          is to one the smaller will be the level. Actual rule is:
		///                          A(i,j)/(sum(A(i,:))+sum(A(:,j))-A(i,j)) > ddpq_tolerance *
		///                          A(imax,jmax)/(sum(A(imax,:))+sum(A(:,jmax))-A(imax,jmax))
		///                          where on imax, jmax maximum is reached.
		///                          works for:
		///                          INNER_DDPQILUC
		/// - "fill_level"         - level of fill for ILU-type preconditioners,
		///                          works for:
		///                          INNER_ILU2 (if LFILL is defined in solver_ilu2.hpp)
		///                          Trilinos, Trilinos_Ifpack
		/// - "pivot_condition"    - delay factoriation of the row of the matrix to the next level, if
		///                          the estimated condition number is above prescribed value.
		///                          works for:
		///                          INNER_MLMPTILUC
		/// - "pivot_diag"         - delay factoriation of the row of the matrix to the next level, if
		///                          the inverted diagonal value is above the prescribed value.
		///                          works for:
		///                          INNER_MLMPTILUC
		/// @see Solver::GetParameter
		void SetParameter(std::string name, std::string value);
		/// Return the number of iterations performed by the last solution.
		/// @see Sparse::Solve
		INMOST_DATA_ENUM_TYPE Iterations() const;
		/// Return the final precondioned residual achieved by the last solution.
		/// @see Sparse::Solve
		INMOST_DATA_REAL_TYPE Residual() const;
		/// Get the reason of convergence or divergence of the last solution.
		/// @see Sparse::Solve
		const std::string ReturnReason() const;
		/// Return the time of preconditioner evaluation by the last solution
		/// @see Sparse::Solve
		double PreconditionerTime() const;
		/// Return the time of iteration stage by the last solution
		/// @see Sparse::Solve
		double IterationsTime() const;
        /// Return the last solution successfullness
        /// @see Sparse::Solve
		bool IsSolved() const;
        /// Get the useful info of the last solution.
        /// @see Sparse::Solve
        const std::string SolutionMetadataLine(const std::string &delimiter) const;
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
		/// Checks if solver available
		/// @param name Solver name
		/// @see Solver::getAvailableSolvers
		static bool isSolverAvailable(std::string name);
		/// Return the list of all available solvers
		/// @see Solver::isSolverAvailable
		static std::vector<std::string> getAvailableSolvers();
		/// Delete solver and associated data.
		~Solver();
	};

	typedef std::vector<std::string>::iterator solvers_names_iterator_t;
	typedef std::vector<SolverParameters>::iterator solver_parameters_iterator_t;

	/// Helper functions for Optimizer module and Solver
	namespace TTSP {

	    void Enable();

	    void Disable();

	    bool isEnabled();

	    bool isDisabled();

	    void SolverOptimize(const std::string &name, const std::string &type, Solver &solver);

	    void SolverOptimizeSaveResult(const std::string &name, const std::string &type, double metrics, bool is_good);

	    void DestroySavedOptimizer(const std::string &name, const std::string &type);

	}

}

#endif

#endif //INMOST_SOLVER_INCLUDED
