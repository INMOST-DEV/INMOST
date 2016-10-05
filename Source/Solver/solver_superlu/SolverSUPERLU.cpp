#include "SolverSUPERLU.h"

namespace INMOST {

		SolverSUPERLU::SolverSUPERLU() {
			set_default_options(&options);
			StatInit(&stat);
			perm_c = NULL;
			perm_r = NULL;
			a_size = 0;
		}

        SolverSUPERLU::SolverSUPERLU(const SolverInterface *other) {
        	throw INMOST::SolverUnsupportedOperation;
        }

        void SolverSUPERLU::Assign(const SolverInterface *other) {
        	throw INMOST::SolverUnsupportedOperation;
        }

        void SolverSUPERLU::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {
        	//read options from file and arguments
        }

        void SolverSUPERLU::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        	int *ia, *ja, nnz = 0;
			double *a;
			int mbeg = A.GetFirstIndex();
			int mend = A.GetLastIndex();
			int size = 0;
			remap = new int[mend - mbeg];
			for(int k = 0; k < mend - mbeg; ++k) {
				Sparse::Row & r = A[k + mbeg];
				if (r.Size()) {
					double nrm = 0;
					for(INMOST_DATA_ENUM_TYPE l = 0; l < r.Size(); ++l) {
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
			for(int k = 0; k < mend - mbeg; ++k) {
				if (remap[k] != -1) {
				Sparse::Row & r = A[k + mbeg];
				for(INMOST_DATA_ENUM_TYPE l = 0; l < r.Size(); ++l) {
					if (remap[r.GetIndex(l) - mbeg] != -1) {
							ja[q] = remap[r.GetIndex(l) - mbeg];
							a[q] = r.GetValue(l);
							++q;
						} else { //if( fabs(a[q]) > 1.0e-9 ) 
							std::cout << "Matrix has connections to singular rows" << std::endl;
						}
					}
					ia[f+1] = q;
					f++;
				}
			}
			dCreate_CompCol_Matrix(&(this->A), size, size, nnz, a, ja, ia, SLU_NR, SLU_D, SLU_GE);
			a_size = size;
			perm_c = new int [size];
			perm_r = new int [size];
        }

        bool SolverSUPERLU::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
			double *inout = new double[a_size];
			int mbeg = RHS.GetFirstIndex(), mend = RHS.GetLastIndex();
			for(int k = 0; k < mend - mbeg; ++k) {
				if (remap[k] != -1) {
					inout[remap[k]] = RHS[k + mbeg];
				}
			}
			SuperMatrix B;
			dCreate_Dense_Matrix(&B, a_size, 1, inout, a_size, SLU_DN, SLU_D, SLU_GE);
			dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
			Destroy_SuperMatrix_Store(&B);
			bool ret = (info == 0);
			for (int k = 0; k < mend - mbeg; ++k) {
				if (remap[k] != -1) {
					SOL[k + mbeg] = inout[remap[k]];
				}
			}
			delete [] inout;
			return ret;
        }

        bool SolverSUPERLU::Clear() {
        	Destroy_CompCol_Matrix(&A);
			Destroy_SuperNode_Matrix(&L);
			Destroy_CompCol_Matrix(&U);
			StatFree(&stat);
			if( perm_c != NULL ) delete [] perm_c;
			if( perm_r != NULL ) delete [] perm_r;
			if( remap != NULL ) delete [] remap; //allocated outside
			return true;
        }

        bool SolverSUPERLU::isMatrixSet() {
        	return a_size != 0;
        }

        void SolverSUPERLU::SetDefaultParameters() {

        }

        SolverParameter SolverSUPERLU::GetParameter(std::string name) const {
        	std::cout << "SolverSUPERLU::GetParameter unsupported operation" << std::endl;
        	throw INMOST::SolverUnsupportedOperation;
        }

        void SolverSUPERLU::SetParameter(std::string name, std::string value) {
        	//throw INMOST::SolverUnsupportedOperation;
        }

        const INMOST_DATA_ENUM_TYPE SolverSUPERLU::Iterations() const {
        	return 1;
        }

        const INMOST_DATA_REAL_TYPE SolverSUPERLU::Residual() const {
        	return 0;
        }

        const std::string SolverSUPERLU::ReturnReason() const {
			char reason_str[256];
			if (info <= a_size) {
				sprintf(reason_str, "diagonal entry of U-factor is exactly singular at %d/%d", info, a_size);
			} else if (info > a_size) {
				sprintf(reason_str,"memory allocation failed after %d bytes were allocated", info - a_size);
			} else if (info == 0) {
				strcpy(reason_str, "factorization was successfull");
			} else {
				sprintf(reason_str,"unknown exit code %d", info);
			}
			return std::string(reason_str);
        }

        const std::string SolverSUPERLU::SolverName() const {
        	return "superlu";
        }

        void SolverSUPERLU::Finalize() {
        	//nothing to do here
        }

        SolverSUPERLU::~SolverSUPERLU() {
        	this->Clear();
        }


}