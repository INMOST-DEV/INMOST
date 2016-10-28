#ifndef INMOST_SOLVERINTERFACE_H
#define INMOST_SOLVERINTERFACE_H

#include <string>
#include <sstream>
#include "inmost_sparse.h"

namespace INMOST {

#if defined(USE_SOLVER)

    class SolverParameter {
    private:
        std::string value;
    public:
        SolverParameter() {};

        SolverParameter(std::string value) : value(value) {};

        template<typename T>
        const T get() const {
            std::istringstream stream(value);
            T parameter;
            stream >> parameter;
            return parameter;
        }
    };

    class SolverParameters {
    private:
        std::map<std::string, SolverParameter> parameters;
    public:
        SolverParameters();

        void SetParameter(std::string name, std::string value);

        const SolverParameter GetParameter(std::string name) const;

        template<typename T>
        const T get(std::string name) const {
            return GetParameter(name).get<T>();
        }
    };

    class SolverInterface {
    protected:
        INMOST_MPI_Comm communicator;
        SolverParameters parameters;
    public:
        SolverInterface() {};

        SolverInterface(const SolverInterface *other) {};

        virtual void Assign(const SolverInterface *other) = 0;

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) = 0;

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) = 0;

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) = 0;

        virtual bool Clear() = 0;

        virtual bool isMatrixSet() = 0;

        virtual void SetDefaultParameters() = 0;

        virtual SolverParameter GetParameter(std::string name) const {
            return parameters.GetParameter(name);
        };

        virtual void SetParameter(std::string name, std::string value) {
            parameters.SetParameter(name, value);
        };

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


#endif //INMOST_SOLVERINTERFACE_H
