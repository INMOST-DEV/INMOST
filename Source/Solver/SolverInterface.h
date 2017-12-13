#ifndef INMOST_SOLVER_INTERFACE_H
#define INMOST_SOLVER_INTERFACE_H

#include "inmost_sparse.h"
#include <string>
#include <sstream>


#define SILENCE_SET_PARAMETER
#if defined(USE_SOLVER)
namespace INMOST {

    class SolverParameters {
    public:
        std::string solverName;
        std::string solverPrefix;
        std::string internalFile;

        std::vector<std::pair<std::string, std::string> > parameters;

        SolverParameters(std::string solverName, std::string solverPrefix, std::string internalFile);
        SolverParameters(const SolverParameters &other);

        ~SolverParameters();
    };

    typedef std::vector<std::pair<std::string, std::string> >::iterator parameters_iterator_t;

    class SolverInterface {
    protected:
        INMOST_MPI_Comm communicator;
    public:
        SolverInterface() {};

        virtual SolverInterface *Copy(const SolverInterface *other) = 0;

        virtual void Assign(const SolverInterface *other) = 0;

        virtual void Setup(int *argc, char ***argv, SolverParameters &p) = 0;

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) = 0;

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) = 0;

        virtual bool Clear() = 0;

        virtual bool isMatrixSet() = 0;

        virtual std::string GetParameter(std::string name) const = 0;

        virtual void SetParameter(std::string name, std::string value) = 0;

        virtual INMOST_DATA_ENUM_TYPE Iterations() const = 0;

        virtual INMOST_DATA_REAL_TYPE Residual() const = 0;

        virtual const std::string ReturnReason() const = 0;

        virtual const std::string SolverName() const = 0;

        virtual void Finalize() = 0;

        virtual INMOST_DATA_REAL_TYPE Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter) {
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
#endif //USE_SOLVER


#endif //INMOST_SOLVER_INTERFACE_H
