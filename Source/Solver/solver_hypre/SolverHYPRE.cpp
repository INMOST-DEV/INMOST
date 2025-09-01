#include "SolverHYPRE.h"

namespace INMOST 
{
	

    SolverHYPRE::SolverHYPRE(HypreType type)
    {

        size = 0;
    }

    SolverInterface *SolverHYPRE::Copy(const SolverInterface *other) 
    {
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverHYPRE::Assign(const SolverInterface *other)
    {
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverHYPRE::Setup(int *argc, char ***argv, SolverParameters &p)
    {
        //read options from file and arguments
    }

    void SolverHYPRE::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner)
    {
        
    }

    bool SolverHYPRE::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL)
    {
		bool ret = false;
        return ret;
    }

    bool SolverHYPRE::Clear()
    {
        return true;
    }

    bool SolverHYPRE::isMatrixSet()
    {
        return size != 0;
    }

    std::string SolverHYPRE::GetParameter(std::string name) const
    {
        return "";
    }

    void SolverHYPRE::SetParameter(std::string name, std::string value)
    {
    }

    INMOST_DATA_ENUM_TYPE SolverHYPRE::Iterations() const
    {
        return 1;
    }

    INMOST_DATA_REAL_TYPE SolverHYPRE::Residual() const
    {
        return 0;
    }

    const std::string SolverHYPRE::ReturnReason() const
    {
        char reason_str[256] = "";
        return std::string(reason_str);
    }

    const std::string SolverHYPRE::SolverName() const
    {
        if (type == BoomerAMG)
            return "hypre_boomeramg";
        else if (type == ParaSails)
            return "hypre_parasails";
        else if (type == AMS)
            return "hypre_ams";
        else if (type == PILUT)
            return "hypre_pilut";
        else if(type == Euclid)
            return "hypre_euclid";
        return "hypre";
    }

    void SolverHYPRE::Finalize()
    {
        //nothing to do here
    }

    SolverHYPRE::~SolverHYPRE()
    {
        this->Clear();
    }


}
