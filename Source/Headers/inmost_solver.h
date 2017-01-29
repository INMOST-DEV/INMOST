#ifndef INMOST_SOLVER_INCLUDED
#define INMOST_SOLVER_INCLUDED

#include "inmost_common.h"
#include "inmost_sparse.h"

#if defined(USE_SOLVER)
namespace INMOST
{

    class SolverInterface;
    class SolverParameters;

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
		void SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value);
		void SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value);
		std::string GetReason() {return ReturnReason();}
		//Backward-compatibility
    private:
        static std::vector<SolverParameters> parameters;
        static int *argc;
        static char ***argv;
        static bool is_initialized;
        static bool is_finalized;

        //Actual solver using for solving system
        SolverInterface *solver;
        std::string prefix;
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
            int rank, size;
        public:
            void Clear();

            /// Return true if Matrix data have already been specified.
            bool &HaveMatrix() { return have_matrix; }

            OrderInfo();

            OrderInfo(const OrderInfo &other);

            OrderInfo &operator=(OrderInfo const &other);

            ~OrderInfo();

            /// Prepare parallel state of the Matrix with specified overlap size.
            /// This state of the matrix can be used, for instance, to construct
            /// the preconditioner for Additive Swartz method.
            /// @param m Matrix to be expanded.
            /// @param overlap Overlap size, viz. the number of overlap layers.
            void PrepareMatrix(Sparse::Matrix &m, INMOST_DATA_ENUM_TYPE overlap);

            /// Restore initial nonparallel state of the Matrix with no overlap.
            void RestoreMatrix(Sparse::Matrix &m);

            /// Prepare parallel state of the Vector.
            void PrepareVector(Sparse::Vector &v) const;

            /// Restore initial nonparallel state of the Vector.
            void RestoreVector(Sparse::Vector &v) const;
            /// Retrieve the processor number by binary search for the specified global index.
            INMOST_DATA_ENUM_TYPE GetProcessor(
                    INMOST_DATA_ENUM_TYPE gind) const; //retrieve processor by binary search in global_to_proc
            void GetOverlapRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE &mbeg,
                                  INMOST_DATA_ENUM_TYPE &mend) const;

            /// Get the local index region for the specified process.
            void
            GetLocalRegion(INMOST_DATA_ENUM_TYPE proc, INMOST_DATA_ENUM_TYPE &mbeg, INMOST_DATA_ENUM_TYPE &mend) const;

            /// Get the local index region for the current process.
            void GetVectorRegion(INMOST_DATA_ENUM_TYPE &mbeg, INMOST_DATA_ENUM_TYPE &mend) const
            {
                mbeg = local_vector_begin;
                mend = local_vector_end;
            }
            /// Get the rank of the current communicator, i.e. the current process index.
            INMOST_DATA_ENUM_TYPE GetRank() const { return rank; }
            /// Get the size of the current communicator, i.e. the total number of processes used.
            INMOST_DATA_ENUM_TYPE GetSize() const { return size; }

            /// Update the shared data in parallel vector.
            void Update(Sparse::Vector &x); // update parallel vector
            /// Sum shared values in parallel vector.
            void Accumulate(Sparse::Vector &x); // sum shared values in parallel vector
            /// Get the sum of num elements of real array on all processes.
            void Integrate(INMOST_DATA_REAL_TYPE *inout, INMOST_DATA_ENUM_TYPE num) const;
            /// Get the communicator which the solver is associated with.
            INMOST_MPI_Comm GetComm() const { return comm; }
            // Access to arrays below allows to organize manual exchange
            INMOST_MPI_Request *GetSendRequests()
            {
                assert(!send_requests.empty());
                return &send_requests[0];
            }

            INMOST_MPI_Request *GetRecvRequests()
            {
                assert(!recv_requests.empty());
                return &recv_requests[0];
            }

            INMOST_DATA_ENUM_TYPE GetSendRequestsSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(send_requests.size()); }

            INMOST_DATA_ENUM_TYPE GetRecvRequestsSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(recv_requests.size()); }

            INMOST_DATA_ENUM_TYPE *GetSendExchangeArray()
            {
                assert(!vector_exchange_send.empty());
                return &vector_exchange_send[0];
            }

            INMOST_DATA_ENUM_TYPE GetSendExchangeSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(send_storage.size()); }

            INMOST_DATA_ENUM_TYPE *GetRecvExchangeArray()
            {
                assert(!vector_exchange_recv.empty());
                return &vector_exchange_recv[0];
            }

            INMOST_DATA_ENUM_TYPE GetRecvExchangeSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(recv_storage.size()); }
            //for debug
            //~ void                 BeginSequentialCode() {for(int i = 0; i < rank; i++) MPI_Barrier(comm);}
            //~ void                   EndSequentialCode() {for(int i = rank; i < size; i++) MPI_Barrier(comm);}

            // Get the scalar product of the specified interval of the distributed vector.
            // Conflicts with OpenMP, should not be used in future
            //void ScalarProd(Vector const & left, Vector const & right, INMOST_DATA_ENUM_TYPE index_begin, INMOST_DATA_ENUM_TYPE index_end, INMOST_DATA_REAL_TYPE & sum) const;
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

        Solver(const Solver &other);

        Solver &operator=(const Solver &other);

        /// Return the solver name 
        /// @see Sparse::Solve
        std::string SolverName() const;

        /// Return the solver user specified name of the current solver
        /// @see Sparse::Solve
        std::string SolverPrefix() const;

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
        bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

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
        /// @see Solver::GetParameter
        void SetParameter(std::string name, std::string value);

        /// Return the number of iterations performed by the last solution.
        /// @see Sparse::Solve
        const INMOST_DATA_ENUM_TYPE Iterations() const;

        /// Return the final residual achieved by the last solution.
        /// @see Sparse::Solve
        const INMOST_DATA_REAL_TYPE Residual() const;

        /// Get the reason of convergence or divergence of the last solution.
        /// @see Sparse::Solve
        const std::string ReturnReason() const;

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

        ~Solver();

    private:

        /// Reads the parameters from database file stored in xml format.
        static void parseXMLDatabase(const char* xml_database);
    };

    typedef std::vector<std::string>::iterator solvers_names_iterator_t;
    typedef std::vector<SolverParameters>::iterator solver_parameters_iterator_t;

}

#endif

#endif //INMOST_SOLVER_INCLUDED
