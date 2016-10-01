#include "SolverTrilinosBelos.h"

namespace INMOST {

    SolverTrilinosBelos::SolverTrilinosBelos() {

    }

    SolverTrilinosBelos::SolverTrilinosBelos(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    bool SolverTrilinosBelos::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
        std::string name = Epetra_problem->first;
        Epetra_LinearProblem *Epetra_linear_problem = &Epetra_problem->second;
        Epetra_Vector VectorRHS(View, matrix->Map(), &*RHS.Begin());
        Epetra_Vector VectorSOL(View, matrix->Map(), &*SOL.Begin());
        Epetra_linear_problem->SetRHS(&VectorRHS);
        Epetra_linear_problem->SetLHS(&VectorSOL);

        bool have_params = parameters_file != "";
        const Teuchos::RCP<Teuchos::ParameterList> top_level_params = Teuchos::createParameterList();
        if (have_params)
            Teuchos::updateParametersFromXmlFileAndBroadcast(parameters_file, top_level_params.ptr(),
                                                             Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(
                                                                     RHS.GetCommunicator())));

        Teuchos::RCP<Teuchos::ParameterList> List = Teuchos::rcp(new Teuchos::ParameterList);

        if (have_params && top_level_params->isSublist(name) && top_level_params->sublist(name).isSublist("Belos"))
            *List = top_level_params->sublist(name).sublist("Belos");
        else {
            List->set("Num Blocks", 100);
            List->set("Block Size", 1);
            List->set("Maximum Iterations", static_cast<int>(parameters.GetParameterEnum("maximum_iterations")));
            List->set("Maximum Restarts", 20);
            List->set("Convergence Tolerance", static_cast<double>(parameters.GetParameterReal("relative_tolerance")));
        }

        Teuchos::RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > Problem =
                Teuchos::rcp(new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>);
        Problem->setLHS(
                Teuchos::rcp_implicit_cast<Epetra_MultiVector, Epetra_Vector>(Teuchos::rcp(&VectorSOL, false)));
        Problem->setRHS(
                Teuchos::rcp_implicit_cast<Epetra_MultiVector, Epetra_Vector>(Teuchos::rcp(&VectorRHS, false)));
        Problem->setOperator(
                Teuchos::rcp_implicit_cast<Epetra_Operator, Epetra_CrsMatrix>(Teuchos::rcp(matrix, false)));

        Problem->setProblem();

        Teuchos::RCP<Belos::SolverManager<double, Epetra_MultiVector, Epetra_Operator> > BelosSolver =
                Teuchos::rcp(new Belos::BlockGmresSolMgr<double, Epetra_MultiVector, Epetra_Operator>);

        BelosSolver->setParameters(List);
        BelosSolver->setProblem(Problem);
        Belos::ReturnType ret = BelosSolver->solve();
        lastIterations = BelosSolver->getNumIters();
        lastResidual = BelosSolver->achievedTol();

        if (ret == Belos::Converged) {
            returnReason = "Converged";
            return true;
        } else {
            std::stringstream reason_stream;
            reason_stream << "Diverged ";
            if (BelosSolver->isLOADetected())
                reason_stream << "loss of accuracy detected ";
            if (BelosSolver->getNumIters() >= 1000)
                reason_stream << "maximum iteration number reached ";
            if (BelosSolver->achievedTol() > 1.0e+6)
                reason_stream << "divergence tolerance reached ";
            returnReason = reason_stream.str();
            return false;
        }

    }

    const std::string SolverTrilinosBelos::SolverName() const {
        return "trilinos_belos";
    }

    SolverTrilinosBelos::~SolverTrilinosBelos() {
        this->Clear();
    }
}