//
// Created by Dmitri Bagaev on 22.09.16.
//

#ifndef INMOST_SOLVER_INCLUDED
#define INMOST_SOLVER_INCLUDED

#include "inmost_solver_interface.h"
#include "inmost_common.h"
#include "inmost_sparse.h"

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
namespace INMOST {

#define GUARD_MPI(x) {ierr = x; if( ierr != MPI_SUCCESS ) {char str[4096]; int len; MPI_Error_string(ierr,str,&len); std::cout << #x << " not successfull: " << str << std::endl; MPI_Abort(comm,-1000);}}
#define HASH_TABLE_SIZE 2048

    class Solver {
    private:
        static int *argc;
        static char ***argv;
        static const char *database;
        static bool is_initialized;
        static bool is_finalized;

        //Actual solver using for solving system
        SolverInterface *solver;
        std::string prefix;
    public:

        /// Base class for low level operations with objects of Solver class.
        class OrderInfo {
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
            void GetVectorRegion(INMOST_DATA_ENUM_TYPE &mbeg, INMOST_DATA_ENUM_TYPE &mend) const {
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
            INMOST_MPI_Request *GetSendRequests() {
                assert(!send_requests.empty());
                return &send_requests[0];
            }

            INMOST_MPI_Request *GetRecvRequests() {
                assert(!recv_requests.empty());
                return &recv_requests[0];
            }

            INMOST_DATA_ENUM_TYPE GetSendRequestsSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(send_requests.size()); }

            INMOST_DATA_ENUM_TYPE GetRecvRequestsSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(recv_requests.size()); }

            INMOST_DATA_ENUM_TYPE *GetSendExchangeArray() {
                assert(!vector_exchange_send.empty());
                return &vector_exchange_send[0];
            }

            INMOST_DATA_ENUM_TYPE GetSendExchangeSize() { return static_cast<INMOST_DATA_ENUM_TYPE>(send_storage.size()); }

            INMOST_DATA_ENUM_TYPE *GetRecvExchangeArray() {
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

        // Main constructor of the solver.
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
        ///     PETSc: petsc_options.txt
        ///     Trilinos_Ifpack: trilinos_ifpack_options.xml
        ///     Trilinos_ML: trilinos_ml_options.xml
        ///     Trilinos_Aztec: trilinos_aztec_options.xml
        ///     Trilinos_Belos: trilinos_belos_options.xml
        static void Initialize(int *argc, char ***argv, const char *database);

        /// Finalize the stage of parallel solution.
        /// If MPI was initialized in Solver::Initialize, then it will be finalized.
        /// By this reason, do not use any MPI function after call to this function.
        /// @see Solver::Initialize
        /// @see Solver::isFinalized
        static void Finalize();

        static bool isInitialized();

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

        void SetDefaultParameters();

        SolverParameter GetParameter(std::string property) const;

        void SetParameter(std::string property, std::string value);

        /// Return the number of iterations performed by the last solution.
        /// @see Sparse::Solve
        const INMOST_DATA_ENUM_TYPE Iterations() const;

        /// Return the final residual achieved by the last solution.
        /// @see Sparse::Solve
        const INMOST_DATA_REAL_TYPE Residual() const;

        /// Get the reason of convergence or divergence of the last solution.
        /// @see Sparse::Solve
        const std::string ReturnReason() const;

        ~Solver();

    private:
        static std::string parseDatabase(std::string solverName);
    };
}

#endif

#endif //INMOST_SOLVER_INCLUDED
