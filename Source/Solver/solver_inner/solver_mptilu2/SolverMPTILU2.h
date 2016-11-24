#ifndef INMOST_SOLVERMPTILU2_H
#define INMOST_SOLVERMPTILU2_H


#include <inmost.h>
#include "solver_mtilu2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverMPTILU2 : public SolverInner {
        INMOST_DATA_ENUM_TYPE rescale_iterations, schwartz_overlap, gmres_substeps, reorder_nnz, fill_level;
        INMOST_DATA_REAL_TYPE drop_tolerance, reuse_tolerance;
    public:
        SolverMPTILU2();

        SolverMPTILU2(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual void SetParameter(std::string name, std::string value);

        virtual const std::string SolverName() const;

        virtual ~SolverMPTILU2();
    };

}


#endif //INMOST_SOLVERMPTILU2_H
