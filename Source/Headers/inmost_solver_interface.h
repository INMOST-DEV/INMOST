#ifndef INMOST_SOLVERINTERFACE_H
#define INMOST_SOLVERINTERFACE_H

#include <string>
#include "inmost_sparse.h"

namespace INMOST {

#if defined(USE_SOLVER)

    class SolverInterface {
    protected:
        INMOST_MPI_Comm communicator;
    public:
        SolverInterface() {};
        SolverInterface(const SolverInterface* other) {};
        virtual void Assign(const SolverInterface* other) = 0;

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) = 0;
        virtual void SetMatrix(Sparse::Matrix & A, bool ModifiedPattern, bool OldPreconditioner) = 0;
        virtual bool Solve(INMOST::Sparse::Vector & RHS, INMOST::Sparse::Vector & SOL) = 0;
        virtual bool Clear() = 0;

        virtual bool isMatrixSet() = 0;

        virtual INMOST_DATA_REAL_TYPE GetPropertyReal(std::string property) const = 0;
        virtual INMOST_DATA_ENUM_TYPE GetPropertyEnum(std::string property) const = 0;

        virtual void SetPropertyReal(std::string property, INMOST_DATA_REAL_TYPE value) = 0;
        virtual void SetPropertyEnum(std::string property, INMOST_DATA_ENUM_TYPE value) = 0;

        virtual const INMOST_DATA_ENUM_TYPE Iterations() const = 0;
        virtual const INMOST_DATA_REAL_TYPE Residual() const = 0;
        virtual const std::string ReturnReason() const = 0;

        virtual const std::string SolverName() const = 0;

        virtual void Finalize() = 0;
        virtual ~SolverInterface() {};


        void SetCommunicator(INMOST_MPI_Comm _communicator) {
            communicator = _communicator;
        }

        INMOST_MPI_Comm GetCommunicator() {
            return communicator;
        }
    };

}
#endif


#endif //INMOST_SOLVERINTERFACE_H
