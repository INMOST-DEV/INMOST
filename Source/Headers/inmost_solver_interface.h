#ifndef INMOST_SOLVER_INTERFACE_H
#define INMOST_SOLVER_INTERFACE_H

#include <string>
#include <sstream>
#include "inmost_solver_parameters.h"
#include "inmost_sparse.h"

namespace INMOST {

#if defined(USE_SOLVER)

    class SolverInterface {
    protected:
        INMOST_MPI_Comm communicator;
        SolverParameters &parameters;
    public:
        SolverInterface(SolverParameters &parameters): parameters(parameters) {};

        SolverInterface(const SolverInterface *other): parameters(other->parameters) {
            parameters = SolverParameters(other->parameters);
        };

        virtual void Assign(const SolverInterface *other) = 0;

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) = 0;

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) = 0;

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) = 0;

        virtual bool Clear() = 0;

        virtual bool isMatrixSet() = 0;

        virtual SolverParameters &GetParameters() const {
            return parameters;
        };

        virtual void SetParameters(SolverParameters &_parameters) {
            parameters = _parameters;
        }

        virtual const INMOST_DATA_ENUM_TYPE Iterations() const = 0;

        virtual const INMOST_DATA_REAL_TYPE Residual() const = 0;

        virtual const std::string ReturnReason() const = 0;

        virtual const std::string SolverName() const = 0;

        virtual void Finalize() = 0;

        virtual const INMOST_DATA_REAL_TYPE Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter) {
            throw INMOST::SolverUnsupportedOperation;
        };

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


#endif //INMOST_SOLVER_INTERFACE_H
