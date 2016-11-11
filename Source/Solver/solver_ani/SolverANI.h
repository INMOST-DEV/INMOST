#ifndef INMOST_SOLVERANI_H
#define INMOST_SOLVERANI_H

#include "inmost_solver_interface.h"
#include "solver_ani.h"

namespace INMOST {

    class SolverANI : public SolverInterface {
    private:
        bcg solver;
        matrix m;

        INMOST_DATA_ENUM_TYPE local_size;
    public:
        SolverANI(SolverParameters &parameters);

        SolverANI(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual const INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual const INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const;

        virtual void Finalize();

        virtual ~SolverANI();

    };

}


#endif //INMOST_SOLVERANI_H
