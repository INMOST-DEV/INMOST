//
// Created by Dmitri Bagaev on 22.09.16.
//

#ifndef INMOST_SOLVERCONTAINER_H
#define INMOST_SOLVERCONTAINER_H

#include <string>
#include "SolverInterface.h"

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

namespace INMOST {

    class Solver2 {
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

        Solver2(std::string solverName, std::string prefix = "", INMOST_MPI_Comm _comm = INMOST_MPI_COMM_WORLD);
        Solver2(const Solver2& other);
        Solver2& operator =(const Solver2& other);

        std::string getSolverName() const;
        std::string getSolverPrefix() const;

        static void Initialize(int *argc, char ***argv, const char *database);
        static bool isInitialized();
        static bool isFinalized();
        static void Finalize();

        void SetMatrix(Sparse::Matrix & A, bool ModifiedPattern = true, bool OldPreconditioner = false);
        bool Solve(INMOST::Sparse::Vector & RHS, INMOST::Sparse::Vector & SOL);

        const INMOST_DATA_ENUM_TYPE Iterations() const;
        const INMOST_DATA_REAL_TYPE Residual() const;
        const std::string ReturnReason() const;

        ~Solver2();

        static std::string parseDatabase(std::string solverName);
    };
}


#endif //INMOST_SOLVERCONTAINER_H
