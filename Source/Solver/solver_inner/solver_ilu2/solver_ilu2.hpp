
#ifndef __SOLVER_ILU2__
#define __SOLVER_ILU2__

#include <iomanip>

#include "inmost_solver.h"
#include "../solver_prototypes.hpp"
//#define REPORT_ILU
//#define REPORT_ILU_PROGRESS
using namespace INMOST;

#define DEFAULT_TAU 0.005
#define DEFAULT_TAU2 0.00001
//#define LFILL //control, that factorization is not less then fill for ilu2



class ILU2_preconditioner : public Method {
private:

    Sparse::Matrix *Alink;
    Solver::OrderInfo *info;
    //Sparse::Matrix L,U;
    //Sparse::Vector div;
    std::vector<INMOST_DATA_REAL_TYPE> luv;
    std::vector<INMOST_DATA_ENUM_TYPE> lui;
    interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ilu, iu;
    INMOST_DATA_ENUM_TYPE Lfill;
    INMOST_DATA_REAL_TYPE tau, tau2;
    Sparse::Vector DL, DR;
    INMOST_DATA_ENUM_TYPE nnz, sciters;
    bool init;
public:
    INMOST_DATA_REAL_TYPE &RealParameter(std::string name) {
        if (name == "tau") return tau;
        else if (name == "tau2") return tau2;
        throw -1;
    }

    INMOST_DATA_ENUM_TYPE &EnumParameter(std::string name) {
        if (name == "fill") return Lfill;
        else if (name == "scale_iters") return sciters;
        throw -1;
    }

    ILU2_preconditioner(Solver::OrderInfo &info)
            : info(&info), tau(DEFAULT_TAU), tau2(DEFAULT_TAU2) {
        Alink = NULL;
        init = false;
        sciters = 12;
        Lfill = 1;
    }

    bool Initialize() {
        if (isInitialized()) Finalize();
        assert(Alink != NULL);
        nnz = 0;
        for (Sparse::Matrix::iterator it = (*Alink).Begin(); it != (*Alink).End(); ++it) nnz += it->Size();
#if defined(LFILL)
        std::vector<INMOST_DATA_ENUM_TYPE> lfill;
        lfill.reserve(nnz * 4);
#endif
        luv.reserve(nnz * 4);
        lui.reserve(nnz * 4);


        std::vector<INMOST_DATA_REAL_TYPE> rv;
        std::vector<INMOST_DATA_ENUM_TYPE> ri;
        rv.reserve(nnz * 16);
        ri.reserve(nnz * 16);
        INMOST_DATA_ENUM_TYPE mobeg, moend, vlocbeg, vlocend, vbeg, vend, k, r, end, iter, j;
        INMOST_DATA_REAL_TYPE leabs, flin, ldiag, udiag, mva;
        INMOST_DATA_ENUM_TYPE curr, foll;
        INMOST_DATA_ENUM_TYPE ind, jn;
        INMOST_DATA_INTEGER_TYPE prev, ipred;
        const INMOST_DATA_REAL_TYPE tol_modif = 1e-12, eps = 1.0e-54, subst = 1.0;
        const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;
        //Calculate scaling vectors for matrix (from genebs)
        info->GetOverlapRegion(info->GetRank(), mobeg, moend);
        info->GetLocalRegion(info->GetRank(), vlocbeg, vlocend);
        info->GetVectorRegion(vbeg, vend);
        interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ir(mobeg, moend + 1);
        interval<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_REAL_TYPE> RowValues(vbeg, vend);
#if defined(LFILL)
        interval<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_ENUM_TYPE> RowFill(vbeg, vend);
        //std::fill(RowFill.begin(),RowFill.end(),ENUMUNDEF);
#endif
        interval<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_ENUM_TYPE> RowIndeces(vbeg - 1, vend);

        ilu.set_interval_beg(mobeg);
        ilu.set_interval_end(moend + 1);
        iu.set_interval_beg(mobeg);
        iu.set_interval_end(moend);
        ilu[mobeg] = 0;
        ir[mobeg] = 0;
#if defined(REPORT_ILU)
        std::cout << "Matrix overlap    " << mobeg << ".." << moend << std::endl;
        std::cout << "Local vector part " << vlocbeg << ".." << vlocend << std::endl;
        std::cout << "Entire vector     " << vbeg << ".." << vend << std::endl;
#endif
        //Rescale Matrix
        DL.SetInterval(mobeg, moend);
        info->PrepareVector(DR);
        for (k = mobeg; k < moend; k++) {
            for (Sparse::Row::iterator rit = (*Alink)[k].Begin(); rit != (*Alink)[k].End(); ++rit)
                DL[k] += rit->second * rit->second;
            if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
        }
        for (iter = 0; iter < sciters; iter++) {
            for (Sparse::Vector::iterator rit = DR.Begin(); rit != DR.End(); ++rit) *rit = 0.0;
            for (k = vlocbeg; k < vlocend; k++)
                for (Sparse::Row::iterator rit = (*Alink)[k].Begin(); rit != (*Alink)[k].End(); ++rit)
                    DR[rit->first] += DL[k] * rit->second * rit->second;
            info->Accumulate(DR);
            info->Update(DR);
            for (k = vlocbeg; k < vlocend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
            for (Sparse::Vector::iterator rit = DL.Begin(); rit != DL.End(); ++rit) *rit = 0.0;
            for (k = mobeg; k < moend; k++)
                for (Sparse::Row::iterator rit = (*Alink)[k].Begin(); rit != (*Alink)[k].End(); ++rit)
                    DL[k] += DR[rit->first] * rit->second * rit->second;
            for (k = mobeg; k < moend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
        }
        for (k = mobeg; k < moend; k++) DL[k] = sqrt(DL[k]);
        for (k = vbeg; k < vend; k++) DR[k] = sqrt(DR[k]);
        for (k = mobeg; k < moend; k++)
            for (Sparse::Row::iterator rit = (*Alink)[k].Begin(); rit != (*Alink)[k].End(); ++rit)
                rit->second = DL[k] * rit->second * DR[rit->first];

        //timer = Timer();
        for (interval<INMOST_DATA_INTEGER_TYPE, INMOST_DATA_ENUM_TYPE>::iterator it = RowIndeces.begin();
             it != RowIndeces.end(); ++it)
            *it = UNDEF;
        std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
        //INMOST_DATA_ENUM_TYPE nza = 0, nzl = 0, nzu = 0, nzu2 = 0;
        //for(k = mobeg; k != moend; k++) nza += A[k].Size();
        for (k = mobeg; k != moend; k++) {
#if defined(REPORT_ILU_PROGRESS)
            if (k % 1000 == 0)
            {
                //std::cout << "precond: " << (double)(k-mobeg)/(double)(moend-mobeg)*100 << "\r";
                //printf("%6.2f nza %12d nzl %12d nzu %12d nzu2 %12d\r", (double)(k-mobeg)/(double)(moend-mobeg)*100,nza,nzl,nzu,nzu2);
                printf("precond: %6.2f\r", (double)(k - mobeg) / (double)(moend - mobeg) * 100);
                fflush(stdout);
            }
#endif
            //Uncompress row
            //row_uncompr
            Sparse::Row &Ak = (*Alink)[k];
            end = Ak.Size();
            sort_indeces.clear();
            for (r = 0; r < end; r++)
                if (fabs(Ak.GetValue(r)) > eps) {
                    RowValues[Ak.GetIndex(r)] = Ak.GetValue(r);
#if defined(LFILL)
                    RowFill[Ak.GetIndex(r)] = 0;
#endif
                    ind = Ak.GetIndex(r);
                    sort_indeces.push_back(ind);
                }
            std::sort(sort_indeces.begin(), sort_indeces.end());
            prev = static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg) - 1;
            ipred = static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg) - 1;
            for (r = 0; r < sort_indeces.size(); r++) {
                ind = sort_indeces[r];
                RowIndeces[prev] = ind;
                prev = static_cast<INMOST_DATA_INTEGER_TYPE>(ind);
                if (ind <= k) ipred = ind;
            }
            RowIndeces[prev] = EOL;

            if (ipred != static_cast<INMOST_DATA_INTEGER_TYPE>(k)) {
                RowValues[k] = 0.0;
#if defined(LFILL)
                RowFill[k] = 0;
#endif
                ind = RowIndeces[ipred];
                RowIndeces[ipred] = k;
                RowIndeces[k] = ind;
            }
#if defined(DIAGONAL_PERTURBATION)
            RowValues[k] = RowValues[k]*(1.0+DIAGONAL_PERTURBATION_REL) + (RowValues[k] < 0.0? -1.0 : 1.0)*DIAGONAL_PERTURBATION_REL;
#endif
            //Eliminate lower part
            //elim_lpart
            j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg) - 1];
            while (j < k) //until diagonal entry
            {
                assert(lui[iu[j]] == j);
                RowValues[j] *= luv[iu[j]]; //scale by diagonal
                leabs = fabs(RowValues[j]);
#if defined(LFILL)
                if (leabs > tau2*tau2)// introduce a non-zero, if threshold permits
#else
                if (leabs > tau2)// introduce a non-zero, if threshold permits
#endif
                {
                    curr = j;
                    for (r = iu[j] + 1; r < ilu[j + 1]; r++) {
                        ind = lui[r];
                        if (RowIndeces[ind] != UNDEF) //update without pondering on thresholds
                        {
                            RowValues[ind] -= RowValues[j] * luv[r];
#if defined(LFILL)
                            RowFill[ind] = std::min(lfill[r]+1,RowFill[ind]);
#endif
                        } else {
                            flin = -RowValues[j] * luv[r];
                            //insert new value
                            foll = curr;
                            while (foll < ind) {
                                curr = foll;
                                foll = RowIndeces[curr];
                            }
                            assert(curr < ind);
                            assert(ind < foll);
                            RowIndeces[curr] = ind;
                            RowIndeces[ind] = foll;
                            RowValues[ind] = flin;
#if defined(LFILL)
                            RowFill[ind] = lfill[r] + 1;
#endif
                        }
                        curr = ind;
                    }

                    if (leabs > tau) {
                        curr = j;
                        for (r = ir[j]; r < ir[j + 1]; r++) {
                            //ind = U2j.GetIndex(r);
                            ind = ri[r];
                            if (RowIndeces[ind] != UNDEF) //update without pondering on thresholds
                                RowValues[ind] -= RowValues[j] * rv[r];
                            else // introduce a non-zero if threshold permits
                            {
                                flin = -RowValues[j] * rv[r];
                                //insert new value
                                foll = curr;
                                while (foll < ind) {
                                    curr = foll;
                                    foll = RowIndeces[curr];
                                }
                                assert(curr < ind);
                                assert(ind < foll);
                                RowIndeces[curr] = ind;
                                RowIndeces[ind] = foll;
                                RowValues[ind] = flin;
#if defined(LFILL)
                                RowFill[ind] = ENUMUNDEF;
#endif
                            }
                            curr = ind;
                        }
                    }
                }
                j = RowIndeces[j];
            }
            // Compress row
            //row_compr
            j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg) - 1];
            //find minimum value in row
            ldiag = 0;
            while (j != EOL) {
                INMOST_DATA_REAL_TYPE temp = fabs(RowValues[j]);
                ldiag = std::max(ldiag, temp);
                j = RowIndeces[j];
            }
            if (ldiag < tau2) {
                ldiag = 1.0 / tau2;
                //std::cout << "ldiag too small " << ldiag << std::endl;
            } else
                ldiag = 1.0 / ldiag;

            //if (ldiag > 1000) std::cout << "ldiag is big " << k << " " << ldiag << std::endl;
            //divide all entries on right from the diagonal
            j = k;
            while (j != EOL) {
                RowValues[j] *= ldiag;
                j = RowIndeces[j];
            }
            j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg) - 1];
            while (j < k) {
                mva = fabs(RowValues[j]);
                if (mva > tau2 * tau2) {
                    if (mva > tau
#if defined(LFILL)
                        || RowFill[j] <= Lfill
#endif
                            ) {
                        //L[k][j] = RowValues[j];
                        lui.push_back(j); //lui indicates column index of L matrix
                        luv.push_back(RowValues[j]); //luv indicates corresponding value
#if defined(LFILL)
                        lfill.push_back(RowFill[j]);
#endif
                        //nzl++;
                    }
                }
                jn = RowIndeces[j];
                RowIndeces[j] = UNDEF;
                j = jn;
            }
            //add last diagonal entry to L matrix
            lui.push_back(k);
            luv.push_back(ldiag);
#if defined(LFILL)
            lfill.push_back(0);
#endif
            //nzl++;

            iu[k] = static_cast<INMOST_DATA_ENUM_TYPE>(luv.size()); //iu points to the first entry of current line of U matrix
            // END of L-part
            if (fabs(RowValues[j]) > tol_modif)
                udiag = 1.0 / RowValues[j];
            else {
                //std::cout << "udiag too small " << RowValues[j] << std::endl;
                udiag = (RowValues[j] < 0.0 ? -1.0 : 1.0) / tol_modif;
            }

            //if (fabs(udiag) > 1000) std::cout << "udiag is big " << k << " " << udiag << std::endl;

            jn = RowIndeces[j];
            RowIndeces[j] = UNDEF;
            j = jn;
            //start of U matrix entries
            //add diagonal value for U matrix
            lui.push_back(k);
            luv.push_back(udiag);
#if defined(LFILL)
            lfill.push_back(RowFill[k]);
#endif
            //nzu++;
            while (j != EOL) {
                mva = fabs(RowValues[j]);
                if (mva > tau2 * tau2) {
                    if (mva > tau
#if defined(LFILL)
                        || RowFill[j] <= Lfill
#endif
                            ) {
                        //add values to U matrix
                        lui.push_back(j);
                        luv.push_back(RowValues[j]);
#if defined(LFILL)
                        lfill.push_back(RowFill[j]);
#endif
                        //nzu++;
                    } else if (mva > tau2) {
                        //add values to U2 matrix
                        ri.push_back(j);
                        rv.push_back(RowValues[j]);
                        //nzu2++;
                    }
                }
                jn = RowIndeces[j];
                RowIndeces[j] = UNDEF;
                j = jn;
            }
            ilu[k + 1] = static_cast<INMOST_DATA_ENUM_TYPE>(luv.size()); //next first entry for L

            ir[k + 1] = static_cast<INMOST_DATA_ENUM_TYPE>(rv.size()); //next first entry for U2
            //END U-part
        }
        //printf("\n");
        //std::cout << "iluoo_solve: " << Timer() - timer << std::endl;
        //timer = Timer();
        //Rescale LU
        //xxlusc
        for (k = mobeg; k < moend; k++) {
            for (r = iu[k] - 1; r > ilu[k]; r--) {
                luv[r - 1] /= DL[k];
                //LFNORM += luv[r-1]*luv[r-1];
            }
            luv[iu[k] - 1] *= DL[k]; // L diagonal entry
            //LFNORM += luv[iu[k]-1]*luv[iu[k]-1];
        }
        for (k = mobeg; k < moend; k++) {
            for (r = iu[k] + 1; r < ilu[k + 1]; r++) {
                luv[r] /= DR[lui[r]];
                //UFNORM += luv[r]*luv[r];
            }
            luv[iu[k]] *= DR[k]; // U diagonal entry
            //UFNORM += luv[iu[k]]*luv[iu[k]];
        }
        //std::cout << "xxlusc: " << Timer() - timer << " LFNORM " << sqrt(LFNORM) << " UFNORM " << sqrt(UFNORM) << std::endl;
        //timer = Timer();


        //Rescale matrix back
        //matisc
        for (k = mobeg; k < moend; k++)
            for (Sparse::Row::iterator rit = (*Alink)[k].Begin(); rit != (*Alink)[k].End(); ++rit)
                rit->second = rit->second / DL[k] / DR[rit->first];
        //std::cout << "matisc: " << Timer() - timer << std::endl;

#if defined(REPORT_ILU)
        INMOST_DATA_ENUM_TYPE nzu,nzl, nza;
        nzu = 0;
        nzl = 0;
        nza = 0;
        for(INMOST_DATA_ENUM_TYPE k = mobeg; k < moend; k++)
        {
            nzl += iu[k] - ilu[k];
            nzu += ilu[k+1] - iu[k] - 1;
            nza += (*Alink)[k].Size();
        }
        std::cout << "      nonzeros in A = " << nza << std::endl;
        std::cout << "      nonzeros in L = " << nzl - (moend-mobeg) << std::endl;
        std::cout << "      nonzeros in U = " << nzu << std::endl;
        std::cout << "     nonzeros in LU = " << ilu[moend] - 1 << std::endl;
        std::cout << "     nonzeros in U2 = " << ir[moend] - 1 << std::endl;
        //std::cout << __FUNCTION__ << " done" << std::endl;
#endif

        /*
        info.PrepareVector(div);
        std::fill(div.Begin(),div.End(),0);
        for(k = mobeg; k < moend; k++) div[k] = 1.0;
        info.Accumulate(div);
        for(k = mobeg; k < moend; k++) div[k] = 1.0/div[k];
        */
        init = true;
        return true;
    }

    bool isInitialized() { return init; }

    bool Finalize() {
        if (!isFinalized()) {
            luv.clear();
            lui.clear();
            init = false;
        }
        return true;
    }

    bool isFinalized() { return !init; }

    ~ILU2_preconditioner() {
        if (!isFinalized()) Finalize();
    }

    void Copy(const Method *other) {
        const ILU2_preconditioner *b = dynamic_cast<const ILU2_preconditioner *>(other);
        assert(b != NULL);
        info = b->info;
        Alink = b->Alink;
        DL = b->DL;
        DR = b->DR;
        nnz = b->nnz;
        luv = b->luv;
        lui = b->lui;
        iu = b->iu;
        ilu = b->ilu;
    }

    ILU2_preconditioner(const ILU2_preconditioner &other)
            : Method(other) {
        Copy(&other);
    }

    ILU2_preconditioner &operator=(ILU2_preconditioner const &other) {
        Copy(&other);
        return *this;
    }

    bool Solve(Sparse::Vector &input, Sparse::Vector &output) {
        assert(isInitialized());
#if defined(USE_OMP)
#pragma omp single
#endif
        {
            INMOST_DATA_ENUM_TYPE mobeg, moend, r, k, vbeg, vend; //, end;
            info->GetOverlapRegion(info->GetRank(), mobeg, moend);
            info->GetVectorRegion(vbeg, vend);
            for (k = vbeg; k < mobeg; k++) output[k] = 0; //Restrict additive schwartz (maybe do it outside?)
            for (k = mobeg; k < moend; k++) output[k] = input[k];
            for (k = moend; k < vend; k++) output[k] = 0; //Restrict additive schwartz (maybe do it outside?)
            for (k = mobeg; k < moend; k++) //iterate over L part
            {
                for (r = iu[k] - 1; r > ilu[k]; r--)
                    output[k] -= luv[r - 1] * output[lui[r - 1]];
                output[k] *= luv[iu[k] - 1];
            }
            for (k = moend; k > mobeg; k--) //iterate over U part
            {
                for (r = iu[k - 1] + 1; r < ilu[k]; r++)
                    output[k - 1] -= luv[r] * output[lui[r]];
                output[k - 1] *= luv[iu[k - 1]];
            }
        }
        //May assemble partition of unity instead of restriction before accumulation
        //assembly should be done instead of initialization
        info->Accumulate(output);
        return true;
    }

    bool ReplaceMAT(Sparse::Matrix &A) {
        if (isInitialized()) Finalize();
        Alink = &A;
        return true;
    };

    bool ReplaceSOL(Sparse::Vector &x) {
        (void) x;
        return true;
    }

    bool ReplaceRHS(Sparse::Vector &b) {
        (void) b;
        return true;
    }

    Method *Duplicate() { return new ILU2_preconditioner(*this); }
};


#endif //__SOLVER_ILU2__
