#include "SolverTrilinosIfpack.h"

namespace INMOST {

    SolverTrilinosIfpack::SolverTrilinosIfpack() {

    }

    SolverTrilinosIfpack::SolverTrilinosIfpack(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    bool SolverTrilinosIfpack::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
        std::string name = Epetra_problem->first;
        Epetra_LinearProblem *Epetra_linear_problem = &Epetra_problem->second;
        Epetra_Vector VectorRHS(View, matrix->Map(), &*RHS.Begin());
        Epetra_Vector VectorSOL(View, matrix->Map(), &*SOL.Begin());
        Epetra_linear_problem->SetRHS(&VectorRHS);
        Epetra_linear_problem->SetLHS(&VectorSOL);

        bool have_params = parameters_file != "";
        const Teuchos::RCP<Teuchos::ParameterList> top_level_params = Teuchos::createParameterList();
        Teuchos::ParameterList local_list;
        if (have_params) {
            Teuchos::updateParametersFromXmlFileAndBroadcast(parameters_file, top_level_params.ptr(),
                                                             Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(
                                                                     RHS.GetCommunicator())));
            if (!top_level_params->isSublist(name))
                have_params = false;
            else {
                local_list = top_level_params->sublist(name);
            }
        }

        AztecOO AztecSolver(*Epetra_linear_problem);

        if (have_params && local_list.isSublist("AztecOO")) {
            Teuchos::ParameterList AztecOOParams = local_list.sublist("AztecOO");
            if (AztecOOParams.isParameter("Max Iterations")) {
                maximum_iterations = AztecOOParams.get<int>("Max Iterations");
            }
            if (AztecOOParams.isParameter("Tolerance")) {
                relative_tolerance = AztecOOParams.get<double>("Tolerance");
            }
            if (AztecOOParams.isSublist("AztecOO Settings")) {
                AztecSolver.SetParameters(AztecOOParams.sublist("AztecOO Settings"));
            }
        } else {
            AztecSolver.SetAztecOption(AZ_diagnostics, AZ_none);
            AztecSolver.SetAztecOption(AZ_output, AZ_none);
            AztecSolver.SetAztecOption(AZ_solver, AZ_bicgstab);
            AztecSolver.SetAztecOption(AZ_overlap, additive_schwartz_overlap);
        }

        Ifpack *Factory = new Ifpack();
        Ifpack_Preconditioner *Prec;
        std::string PrecType = "ILU";
        if (have_params && local_list.isSublist("Ifpack")) {
            Teuchos::ParameterList ifpacklist = local_list.sublist("Ifpack");
            if (ifpacklist.isParameter("Prec Type")) {
                PrecType = ifpacklist.get<std::string>("Prec Type");
            }
            if (ifpacklist.isParameter("Overlap")) {
                additive_schwartz_overlap = ifpacklist.get<int>("Overlap");
            }
        }
        Prec = Factory->Create(PrecType, matrix, additive_schwartz_overlap);
        Teuchos::ParameterList List;
        if (have_params && local_list.isSublist("Ifpack") &&
            local_list.sublist("Ifpack").isSublist("Ifpack Settings")) {
            List = local_list.sublist("Ifpack").sublist("Ifpack Settings");
        } else
            List.set("fact: level-of-fill", static_cast<int>(preconditioner_fill_level));
        Prec->SetParameters(List);
        Prec->Initialize();
        Prec->Compute();
        AztecSolver.SetPrecOperator(Prec);

        AztecSolver.Iterate(maximum_iterations, relative_tolerance);
        const double *stats = AztecSolver.GetAztecStatus();
        bool success = true;
        std::string reason = "";
        checkStatus(static_cast<int>(stats[AZ_why]), success, reason);
        lastIterations = AztecSolver.NumIters();
        lastResidual = AztecSolver.TrueResidual();
        returnReason = reason;
        delete Prec;
        return success;
    }

    const std::string SolverTrilinosIfpack::SolverName() const {
        return "trilinos_ifpack";
    }

    SolverTrilinosIfpack::~SolverTrilinosIfpack() {
        this->Clear();
    }

}