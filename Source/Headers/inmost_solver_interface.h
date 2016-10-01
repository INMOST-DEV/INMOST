#ifndef INMOST_SOLVERINTERFACE_H
#define INMOST_SOLVERINTERFACE_H

#include <string>
#include "inmost_sparse.h"

namespace INMOST {

#if defined(USE_SOLVER)

    class SolverParameters {
    private:
        std::map<std::string, INMOST_DATA_ENUM_TYPE> enums;
        std::map<std::string, INMOST_DATA_REAL_TYPE> reals;
    public:
        SolverParameters();
        void SetDefaults();
        void SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value);
        void SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value);

        const INMOST_DATA_REAL_TYPE GetParameterReal(std::string name) const;
        const INMOST_DATA_ENUM_TYPE GetParameterEnum(std::string name) const;
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

        virtual INMOST_DATA_REAL_TYPE GetParameterReal(std::string name) const {
            return parameters.GetParameterReal(name);
        };

        virtual INMOST_DATA_ENUM_TYPE GetParameterEnum(std::string name) const {
            return parameters.GetParameterEnum(name);
        };

        virtual void SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value) {
            parameters.SetParameterReal(name, value);
        };

        virtual void SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value) {
            parameters.SetParameterEnum(name, value);
        };

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
