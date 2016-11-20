#include "SolverTrilinosML.h"

namespace INMOST {

    SolverTrilinosML::SolverTrilinosML() {

    }

    SolverTrilinosML::SolverTrilinosML(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    bool SolverTrilinosML::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
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
                iters = AztecOOParams.get<int>("Max Iterations");

            }
            if (AztecOOParams.isParameter("Tolerance")) {
                rtol = AztecOOParams.get<double>("Tolerance");
            }
            if (AztecOOParams.isSublist("AztecOO Settings")) {
                AztecSolver.SetParameters(AztecOOParams.sublist("AztecOO Settings"));
            }
        } else {
            AztecSolver.SetAztecOption(AZ_diagnostics, AZ_none);
            AztecSolver.SetAztecOption(AZ_output, AZ_none);
            AztecSolver.SetAztecOption(AZ_solver, AZ_bicgstab);
            AztecSolver.SetAztecOption(AZ_overlap, overlap);
        }

        Teuchos::ParameterList List;
        if (have_params && local_list.isSublist("ML") && local_list.sublist("ML").isSublist("ML Settings")) {
            List = local_list.sublist("ML").sublist("ML Settings");
        } else {
            ML_Epetra::SetDefaults("SA", List);
            List.set("max levels", 6);
            List.set("increasing or decreasing", "decreasing");
        }

        ML_Epetra::MultiLevelPreconditioner *Prec = new ML_Epetra::MultiLevelPreconditioner(*matrix, List, true);
        AztecSolver.SetPrecOperator(Prec);

        AztecSolver.Iterate(iters, rtol);
        const double *stats = AztecSolver.GetAztecStatus();
        bool success = true;
        std::string reason = "";
        TrilinosCheckStatus(static_cast<int>(stats[AZ_why]), success, reason);
        lastIterations = static_cast<INMOST_DATA_ENUM_TYPE>(AztecSolver.NumIters());
        lastResidual = AztecSolver.TrueResidual();
        returnReason = reason;
        delete Prec;
        return success;
    }

    const std::string SolverTrilinosML::SolverName() const {
        return "trilinos_ml";
    }

    SolverTrilinosML::~SolverTrilinosML() {
        this->Clear();
    }

}