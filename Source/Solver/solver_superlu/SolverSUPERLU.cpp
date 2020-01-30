#include "SolverSUPERLU.h"

namespace INMOST 
{
	

    SolverSUPERLU::SolverSUPERLU() 
    {
#if defined(USE_SOLVER_SUPERLU_DIST)
		set_default_options_dist(&options_);
		options_.PrintStat          = NO;
		options_.IterRefine         = SLU_DOUBLE;
		options_.ReplaceTinyPivot   = YES;
		options_.Equil              = YES;
		options_.Fact               = DOFACT;
		PStatInit(&stat_);
#else // USE_SOLVER_SUPERLU_DIST
        set_default_options(&options);
        StatInit(&stat);
        perm_c = NULL;
        perm_r = NULL;
#endif //USE_SOLVER_SUPERLU_DIST
        g_size = a_size = 0;
    }

    SolverInterface *SolverSUPERLU::Copy(const SolverInterface *other) 
    {
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverSUPERLU::Assign(const SolverInterface *other) 
    {
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverSUPERLU::Setup(int *argc, char ***argv, SolverParameters &p) 
    {
        //read options from file and arguments
    }

    void SolverSUPERLU::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) 
    {
        //check that the run is serial!
        int *ia, *ja, nnz = 0;
        double *a;
        int mbeg = A.GetFirstIndex();
        int mend = A.GetLastIndex();
        int size = 0;
#if defined(USE_SOLVER_SUPERLU_DIST)
		int mpisize = 1;
		MPI_Comm_size(GetCommunicator(),&mpisize);
		nnz = 0;
		size = mend-mbeg;
		for (int k = 0; k < size; ++k) 
			nnz += A[k + mbeg].Size();
		ia = (int *) malloc(sizeof(int) * (size + 1));
		ja = (int *) malloc(sizeof(int) * nnz);
		a = (double *) malloc(sizeof(double) * nnz);
		int q = 0, f = 0;
		ia[0] = 0;
		for (int k = 0; k < size; ++k) 
		{
			Sparse::Row &r = A[k + mbeg];
			for (INMOST_DATA_ENUM_TYPE l = 0; l < r.Size(); ++l) 
			{
				ja[q] = r.GetIndex(l);
				a[q] = r.GetValue(l);
				++q;
			}
			ia[f + 1] = q;
			f++;
		}
		int dims[2] = {0,0};
		MPI_Dims_create(mpisize,2,dims);
		//~ std::cout << "Size: " << mpisize << " dims: " << dims[0] << "," << dims[1] << std::endl;
		superlu_gridinit(GetCommunicator(),dims[0],dims[1],&grid);
		int offset = 0;
		int local_size = size;
		int global_size = size;
		MPI_Exscan(&local_size,&offset,1,MPI_INT,MPI_SUM,GetCommunicator());
		MPI_Allreduce(&local_size,&global_size,1,MPI_INT,MPI_SUM,GetCommunicator());
	  //std::cout << "M = " << m << " N = " << n << " local offset " << offset << std::endl;
		//~ dCreate_CompRowLoc_Matrix_dist(&A_,m, n, pCSRMatrixFormat->nonzeros,local_size[0],offset,pCSRMatrixFormat->pData,pCSRMatrixFormat->pCols,pCSRMatrixFormat->pRows,SLU_NR_loc, SLU_D, SLU_GE);
		dCreate_CompRowLoc_Matrix_dist(&(this->A),global_size, global_size, nnz,local_size,offset,a,ja,ia,SLU_NR_loc, SLU_D, SLU_GE);
		ScalePermstructInit(global_size,global_size,&ScalePermstruct);
		LUstructInit(global_size,&LUstruct);
		g_size = global_size;
		a_size = size;
#else //USE_SOLVER_SUPERLU_DIST
		remap = new int[mend - mbeg];
        for (int k = 0; k < mend - mbeg; ++k) {
            Sparse::Row &r = A[k + mbeg];
            if (r.Size()) {
                double nrm = 0;
                for (INMOST_DATA_ENUM_TYPE l = 0; l < r.Size(); ++l) {
                    nrm += r.GetValue(l) * r.GetValue(l);
                }
                if (nrm) {
                    std::sort(r.Begin(), r.End());
                    nnz += r.Size();
                    remap[k] = size;
                    size++;
                } else {
                    remap[k] = -1;
                }
            } else {
                remap[k] = -1;
            }
        }
        ia = (int *) malloc(sizeof(int) * (size + 1));
        ja = (int *) malloc(sizeof(int) * nnz);
        a = (double *) malloc(sizeof(double) * nnz);
        int q = 0, f = 0;
        ia[0] = 0;
        for (int k = 0; k < mend - mbeg; ++k) {
            if (remap[k] != -1) {
                Sparse::Row &r = A[k + mbeg];
                for (INMOST_DATA_ENUM_TYPE l = 0; l < r.Size(); ++l) {
                    if (remap[r.GetIndex(l) - mbeg] != -1) {
                        ja[q] = remap[r.GetIndex(l) - mbeg];
                        a[q] = r.GetValue(l);
                        ++q;
                    } else { //if( fabs(a[q]) > 1.0e-9 )
                        std::cout << "Matrix has connections to singular rows" << std::endl;
                    }
                }
                ia[f + 1] = q;
                f++;
            }
        }
        a_size = size;
        dCreate_CompCol_Matrix(&(this->A), size, size, nnz, a, ja, ia, SLU_NR, SLU_D, SLU_GE);
        
        perm_c = new int[size];
        perm_r = new int[size];
#endif //USE_SOLVER_SUPERLU_DIST
    }

    bool SolverSUPERLU::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) 
    {
		bool ret = false;
#if defined(USE_SOLVER_SUPERLU_DIST)
			//~ std::cout << __FILE__ << ":" << __LINE__ << " AStype " << (A.Stype == SLU_NR_loc) << " ADtype " << (A.Dtype == SLU_D? 1 : 0) << " AGEtype " << (A.Mtype == SLU_GE? 1 : 0) << std::endl;
		double pBerr_[2];
		int mbeg = RHS.GetFirstIndex(), mend = RHS.GetLastIndex();
		SOL.SetInterval(mbeg,mend);
		memcpy(&SOL[mbeg],&RHS[mbeg], (mend-mbeg)*sizeof(double));
		//set_default_options_dist(&options_);
		pdgssvx(&options_, &A, &ScalePermstruct, &SOL[mbeg], a_size, 1, &grid, &LUstruct, &SOLVEstruct, pBerr_, &stat_, &info);
		//~ PStatPrint(&options_, &stat_, &grid);
		options_.Fact = FACTORED;
#else //USE_SOLVER_SUPERLU_DIST
		//~ std::cout << __FILE__ << ":" << __LINE__ << std::endl;
		double *inout = new double[a_size];
		int mbeg = RHS.GetFirstIndex(), mend = RHS.GetLastIndex();
		for (int k = 0; k < mend - mbeg; ++k) 
			if (remap[k] != -1) 
				inout[remap[k]] = RHS[k + mbeg];
			
		
		SuperMatrix B;
		dCreate_Dense_Matrix(&B, a_size, 1, inout, a_size, SLU_DN, SLU_D, SLU_GE);
		dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
		Destroy_SuperMatrix_Store(&B);
		for (int k = 0; k < mend - mbeg; ++k) 
			if (remap[k] != -1) 
				SOL[k + mbeg] = inout[remap[k]];
	
		delete[] inout;
#endif //USE_SOLVER_SUPERLU_DIST
            
        
        ret = (info == 0);
        return ret;
    }

    bool SolverSUPERLU::Clear() 
    {
#if defined(USE_SOLVER_SUPERLU_DIST)
		//~ std::cout << "call clear!" << std::endl;
		PStatFree(&stat_);
		ScalePermstructFree(&ScalePermstruct);
		Destroy_LU(g_size,&grid,&LUstruct);
		LUstructFree(&LUstruct);
		if( options_.SolveInitialized )
		{
			//~ std::cout << "call dSolveFinalize!" << std::endl;
			dSolveFinalize(&options_,&SOLVEstruct);
		}
		Destroy_CompRowLoc_Matrix_dist(&A);
		superlu_gridexit(&grid);
#else //USE_SOLVER_SUPERLU_DIST
		Destroy_CompCol_Matrix(&A);
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		StatFree(&stat);
		if (perm_c != NULL) delete[] perm_c;
		if (perm_r != NULL) delete[] perm_r;
		if (remap != NULL) delete[] remap; //allocated outside
#endif
		//~ SUPERLU_FREE(A.Store);
        return true;
    }

    bool SolverSUPERLU::isMatrixSet() 
    {
        return a_size != 0;
    }

    std::string SolverSUPERLU::GetParameter(std::string name) const 
    {
#if !defined(SILENCE_SET_PARAMETER)
        std::cout << "SolverSUPERLU::GetParameter unsupported operation" << std::endl;
#endif
        //throw INMOST::SolverUnsupportedOperation;
        return "";
    }

    void SolverSUPERLU::SetParameter(std::string name, std::string value) 
    {
#if !defined(SILENCE_SET_PARAMETER)
        std::cout << "SolverSUPERLU::SetParameter unsupported operation" << std::endl;
#endif
        //throw INMOST::SolverUnsupportedOperation;
    }

    INMOST_DATA_ENUM_TYPE SolverSUPERLU::Iterations() const 
    {
        return 1;
    }

    INMOST_DATA_REAL_TYPE SolverSUPERLU::Residual() const 
    {
        return 0;
    }

    const std::string SolverSUPERLU::ReturnReason() const 
    {
        char reason_str[256];
        if (info <= a_size) 
            sprintf(reason_str, "diagonal entry of U-factor is exactly singular at %d/%d", info, a_size);
        else if (info > a_size) 
            sprintf(reason_str, "memory allocation failed after %d bytes were allocated", info - a_size);
        else if (info == 0) 
            strcpy(reason_str, "factorization was successfull");
        else
            sprintf(reason_str, "unknown exit code %d", info);
        return std::string(reason_str);
    }

    const std::string SolverSUPERLU::SolverName() const 
    {
        return "superlu";
    }

    void SolverSUPERLU::Finalize() 
    {
        //nothing to do here
    }

    SolverSUPERLU::~SolverSUPERLU() 
    {
        this->Clear();
    }


}
