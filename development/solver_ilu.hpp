#pragma once
#ifndef __SOLVER_ILU2__
#define __SOLVER_ILU2__



//TODO:
// how to calculate diagonal perturbation?
// how to calculate schur complement faster
// 
// done! implement crout-ILU from 
//   Documents/Read/solver/crout/crout.pdf to ILUC_preconditioner
// done! implement condition estimation for |L^{-1}| and |U^{-1}| and adaptive tau for ILUC_preconditioner from 
//   Documents\Read\solver\Read_now\bollhofer_siam.pdf
// try to make ILUC2_preconditioner - second order ILU from ILUC_preconditioner
// implement diagonal pivoting for ILUC - maybe need to update diagonal at every step
//   goto references [7],[10]-(data structure!) from
//   Documents\Read\solver\crout\crout.pdf
// return dropped values to diagonal if control vector is provided from
//   Documents\Read\solver\stabilization\milut.pdf
// try to apply dropping while forming linked list,should correctly calculate condition estimates 

// Calculate schur complement faster:
//   Documents\Read\solver\sparse_matmul\sparse.pdf

// in ILUC_preconditioner, replace matrix structures by CSR, estimate number of nonzeros in rows/cols
// before filling, if necessery
#include "inmost_solver.h"
#include "solver_prototypes.hpp"
#define REPORT_ILU
#define REPORT_ILU_PROGRESS
using namespace INMOST;

#define DEFAULT_TAU 0.005
#define DEFAULT_TAU2 0.00001

#define REORDER_DDPQ
#define DDPQ_TAU_INV 3
#define ESTIMATOR
#define RESCALE_B
#define PIVOT_THRESHOLD
//#define DIAGONAL_PERTURBATION
#define DIAGONAL_PERTURBATION_REL 0.005
#define DIAGONAL_PERTURBATION_ABS 1.0e-12
#define SCHUR_DROPPING
#define SCHUR_TAU tau
#define ILUC2
#define ILUC2_TAU pow(tau,2.5)



#define LFILL //control, that factorization is not less then fill for ilu2

#define ADDR(row,col) ((row)*size+(col))

class ILU2_preconditioner : public Method
{
private:
	
	Solver::Matrix * Alink;
	Solver::OrderInfo * info;
	//Solver::Matrix L,U;
	//Solver::Vector div;
	std::vector<INMOST_DATA_REAL_TYPE> luv;
	std::vector<INMOST_DATA_ENUM_TYPE> lui;
	interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> ilu,iu;
	INMOST_DATA_ENUM_TYPE nnzA, Lfill;
	INMOST_DATA_REAL_TYPE tau, tau2;
	Solver::Vector DL, DR;
	INMOST_DATA_ENUM_TYPE nnz, sciters;
	bool init;
public:
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name)
	{
		if( name == "tau" ) return tau;
		else if( name == "tau2" ) return tau2;
		throw -1;
	}
	INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name)
	{
		if (name == "fill") return Lfill;
		else if (name == "scale_iters") return sciters;
		throw -1;
	}
	ILU2_preconditioner(Solver::OrderInfo & info)
		:tau(DEFAULT_TAU), tau2(DEFAULT_TAU2), info(&info)
	{
		Alink = NULL;
		init = false;
		sciters = 8;
		Lfill = 1;
	}
	bool Initialize()
	{
		if (isInitialized()) Finalize();
		assert(Alink != NULL);
		nnz = 0;
		for (Solver::Matrix::iterator it = (*Alink).Begin(); it != (*Alink).End(); ++it) nnz += it->Size();
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
		INMOST_DATA_ENUM_TYPE  ind, jn;
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
		for (k = mobeg; k < moend; k++)
		{
			for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r) DL[k] += r->second*r->second;
			if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
		}
		for (iter = 0; iter < sciters; iter++)
		{
			std::fill(DR.Begin(), DR.End(), 0.0);
			for (k = vlocbeg; k < vlocend; k++)
				for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r) DR[r->first] += DL[k] * r->second*r->second;
			info->Accumulate(DR);
			info->Update(DR);
			for (k = vlocbeg; k < vlocend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
			std::fill(DL.Begin(), DL.End(), 0.0);
			for (k = mobeg; k < moend; k++)
				for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r) DL[k] += DR[r->first] * r->second*r->second;
			for (k = mobeg; k < moend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
		}
		for (k = mobeg; k < moend; k++) DL[k] = sqrt(DL[k]);
		for (k =  vbeg; k <  vend; k++) DR[k] = sqrt(DR[k]);
		for (k = mobeg; k < moend; k++)
			for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				r->second = DL[k] * r->second * DR[r->first];

		//timer = Timer();
		std::fill(RowIndeces.begin(), RowIndeces.end(), UNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		//INMOST_DATA_ENUM_TYPE nza = 0, nzl = 0, nzu = 0, nzu2 = 0;
		//for(k = mobeg; k != moend; k++) nza += A[k].Size();
		for (k = mobeg; k != moend; k++)
		{
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
			Solver::Row & Ak = (*Alink)[k];
			end = Ak.Size();
			sort_indeces.clear();
			for (r = 0; r < end; r++) if (fabs(Ak.GetValue(r)) > eps)
			{
				RowValues[Ak.GetIndex(r)] = Ak.GetValue(r);
#if defined(LFILL)
				RowFill[Ak.GetIndex(r)] = 0;
#endif
				ind = Ak.GetIndex(r);
				sort_indeces.push_back(ind);
			}
			std::sort(sort_indeces.begin(), sort_indeces.end());
			prev = static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1;
			ipred = static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1;
			for (r = 0; r < sort_indeces.size(); r++)
			{
				ind = sort_indeces[r];
				RowIndeces[prev] = ind;
				prev = static_cast<INMOST_DATA_INTEGER_TYPE>(ind);
				if (ind <= k) ipred = ind;
			}
			RowIndeces[prev] = EOL;

			if (ipred != static_cast<INMOST_DATA_INTEGER_TYPE>(k))
			{
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
			j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1];
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
					for (r = iu[j] + 1; r < ilu[j + 1]; r++)
					{
						ind = lui[r];
						if (RowIndeces[ind] != UNDEF) //update without pondering on thresholds
						{
							RowValues[ind] -= RowValues[j] * luv[r];
#if defined(LFILL)
							RowFill[ind] = std::min(lfill[r]+1,RowFill[ind]);
#endif
						}
						else 
						{
							flin = -RowValues[j] * luv[r];
							//insert new value
							foll = curr;
							while (foll < ind)
							{
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

					if (leabs > tau)
					{
						curr = j;
						for (r = ir[j]; r < ir[j + 1]; r++)
						{
							//ind = U2j.GetIndex(r);
							ind = ri[r];
							if (RowIndeces[ind] != UNDEF) //update without pondering on thresholds
								RowValues[ind] -= RowValues[j] * rv[r];
							else // introduce a non-zero if threshold permits
							{
								flin = -RowValues[j] * rv[r];
								//insert new value
								foll = curr;
								while (foll < ind)
								{
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
			j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1];
			//find minimum value in row
			ldiag = 0;
			while (j != EOL)
			{
				ldiag = std::max(ldiag, fabs(RowValues[j]));
				j = RowIndeces[j];
			}
			if (ldiag < tau2)
			{
				ldiag = 1.0 / tau2;
				//std::cout << "ldiag too small " << ldiag << std::endl;
			}
			else
				ldiag = 1.0 / ldiag;

			//if (ldiag > 1000) std::cout << "ldiag is big " << k << " " << ldiag << std::endl;
			//divide all entries on right from the diagonal
			j = k;
			while (j != EOL)
			{
				RowValues[j] *= ldiag;
				j = RowIndeces[j];
			}
			j = RowIndeces[static_cast<INMOST_DATA_INTEGER_TYPE>(vbeg)-1];
			while (j < k)
			{
				mva = fabs(RowValues[j]);
				if (mva > tau2*tau2 )
				{
					if (mva > tau
#if defined(LFILL)
						|| RowFill[j] <= Lfill
#endif
						)
					{
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

			iu[k] = luv.size(); //iu points to the first entry of current line of U matrix
			// END of L-part
			if (fabs(RowValues[j]) > tol_modif)
				udiag = 1.0 / RowValues[j];
			else
			{
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
			while (j != EOL)
			{
				mva = fabs(RowValues[j]);
				if (mva > tau2*tau2)
				{
					if (mva > tau
#if defined(LFILL)
						|| RowFill[j] <= Lfill
#endif
						)
					{
						//add values to U matrix
						lui.push_back(j);
						luv.push_back(RowValues[j]);
#if defined(LFILL)
						lfill.push_back(RowFill[j]);
#endif
						//nzu++;
					}
					else if (mva > tau2)
					{
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
			ilu[k + 1] = luv.size(); //next first entry for L

			ir[k + 1] = rv.size(); //next first entry for U2
			//END U-part
		}
		//printf("\n");
		//std::cout << "iluoo_solve: " << Timer() - timer << std::endl;
		//timer = Timer();
		//Rescale LU
		//xxlusc
		for (k = mobeg; k < moend; k++)
		{
			for (r = iu[k] - 1; r > ilu[k]; r--)
			{
				luv[r - 1] /= DL[k];
				//LFNORM += luv[r-1]*luv[r-1];
			}
			luv[iu[k] - 1] *= DL[k]; // L diagonal entry
			//LFNORM += luv[iu[k]-1]*luv[iu[k]-1];
		}
		for (k = mobeg; k < moend; k++)
		{
			for (r = iu[k] + 1; r < ilu[k + 1]; r++)
			{
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
			for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				r->second = r->second / DL[k] / DR[r->first];
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
	bool isInitialized(){ return init; }
	bool Finalize()
	{
		if (!isFinalized())
		{
			info->RestoreMatrix(*Alink);
			luv.clear();
			lui.clear();
			init = false;
		}
		return true;
	}
	bool isFinalized() { return !init; }
	~ILU2_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}
	void Copy(const Method * other)
	{
		const ILU2_preconditioner * b = dynamic_cast<const ILU2_preconditioner *>(other);
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
	ILU2_preconditioner(const ILU2_preconditioner & other)
	{
		Copy(&other);
	}
	ILU2_preconditioner & operator =(ILU2_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	bool Solve(Solver::Vector & input, Solver::Vector & output)
	{
		assert(isInitialized());
		INMOST_DATA_ENUM_TYPE mobeg, moend, r, k, vbeg,vend; //, end;
		info->GetOverlapRegion(info->GetRank(),mobeg,moend);
		info->GetVectorRegion(vbeg,vend);
		for(k = vbeg; k < mobeg; k++) output[k] = 0;
		for(k = mobeg; k < moend; k++) output[k] = input[k];
		for(k = moend; k < vend; k++) output[k] = 0;
		for(k = mobeg; k < moend; k++) //iterate over L part
		{
			for(r = iu[k]-1; r > ilu[k]; r--) 
				output[k] -= luv[r-1]*output[lui[r-1]];
			output[k] *= luv[iu[k]-1];
		}
		for(k = moend; k > mobeg; k--) //iterate over U part
		{
			for(r = iu[k-1]+1; r < ilu[k]; r++)
				output[k-1] -= luv[r]*output[lui[r]];
			output[k-1] *= luv[iu[k-1]];
		}
		info->Accumulate(output);
		return true;
	}
	bool ReplaceMAT(Solver::Matrix & A) { if (isInitialized()) Finalize();  Alink = &A; return true; };
	bool ReplaceSOL(Solver::Vector & x) {return true;}
	bool ReplaceRHS(Solver::Vector & b) {return true;}
	Method * Duplicate() { return new ILU2_preconditioner(*this); }
};

class ILUC_preconditioner : public Method
{
	typedef struct Interval_t
	{
		INMOST_DATA_ENUM_TYPE first, last;
		INMOST_DATA_ENUM_TYPE Size() { return last - first; }
	} Interval;
	typedef struct wgt_t
	{
		INMOST_DATA_REAL_TYPE max, sum;
		INMOST_DATA_ENUM_TYPE ind, nnz;
	} wgt;
	typedef struct wgt_pair_t
	{
		wgt_t row, col;
	} wgt_pair;
	typedef std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> coord;
	typedef std::pair<INMOST_DATA_REAL_TYPE, coord > wgt_coord;
	typedef std::vector< wgt_coord > wgt_coords;
	typedef struct row_col_t
	{
		Solver::Row row, col;
		INMOST_DATA_REAL_TYPE diag;
	} row_col;
	typedef dynarray<INMOST_DATA_ENUM_TYPE,256> sort_inds;
	typedef dynarray<INMOST_DATA_ENUM_TYPE,256> levels_t;
	//result of multilevel preconditioner
	//        |LU   F |
	//  A  =  |       |
	//        | E   C |
	// LU decomposition is stored in skyline format with reversed direction
	std::vector<Solver::Row::entry> LU_Entries;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LU_Diag;
	interval<INMOST_DATA_ENUM_TYPE, Interval> U_Address, L_Address;

	
	interval<INMOST_DATA_ENUM_TYPE, row_col> EF; //remember supplementary reordered E,F blocks
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Efirst, Ffirst; //remember start positions for additional blocks
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> temp; // temporal place for solve phase
	levels_t level_size; //remember size of each level
	//reordering information
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > ddP,ddQ;
	

	INMOST_DATA_REAL_TYPE tau, eps;
	INMOST_DATA_ENUM_TYPE sciters;
	Solver::Matrix * Alink;
	Solver::OrderInfo * info;
	bool init;
	void ReorderEF(INMOST_DATA_ENUM_TYPE mobeg, 
					INMOST_DATA_ENUM_TYPE cbeg,
					INMOST_DATA_ENUM_TYPE wend, 
					interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ, 
					interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
					interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
	{
		INMOST_DATA_ENUM_TYPE i, k;
		//Reorder E rows block by P
		std::fill(donePQ.begin() + cbeg - mobeg, donePQ.begin() + wend - mobeg, false);
		i = cbeg;
		while (i < wend)
		{
			if (donePQ[i]) i++;
			else
			{
				if (localP[i] != i)
				{
					k = i;
					do
					{
						EF[i].row.Swap(EF[localP[k]].row);
						donePQ[localP[k]] = true;
						k = localP[k];
					} while (k != i);
				}
				donePQ[i] = true;
			}
		}
		//Reorder F columns by  Q
		std::fill(donePQ.begin() + cbeg - mobeg, donePQ.begin() + wend - mobeg, false);
		i = cbeg;
		while (i < wend)
		{
			if (donePQ[i]) i++;
			else
			{
				if (localQ[i] != i)
				{
					k = i;
					do
					{
						EF[i].col.Swap(EF[localQ[k]].col);
						donePQ[localQ[k]] = true;
						k = localQ[k];
					} while (k != i);
				}
				donePQ[i] = true;
			}
		}
	}
	void inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
				   INMOST_DATA_ENUM_TYPE wend, 
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		//inverse reordering
		// in invPQ numbers indicate where to get current column
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k != wend; ++k) invP[localP[k]] = k;
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k != wend; ++k) invQ[localQ[k]] = k;
	}
	void applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
				 INMOST_DATA_ENUM_TYPE wend,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		INMOST_DATA_ENUM_TYPE k;
		// compute reordering in global P,Q, we need it to compute reordering in vector during solve phase
		for (k = wbeg; k != wend; ++k)
		{
			localP[k] = ddP[invP[k]];
			localQ[k] = ddQ[invQ[k]];
		}
		// update reordering in global P,Q
		for (k = wbeg; k != wend; ++k)
		{
			ddP[k] = localP[k];
			ddQ[k] = localQ[k];
		}
	}
public:
	INMOST_DATA_ENUM_TYPE & EnumParameter(std::string name)
	{
		if (name == "scale_iters") return sciters;
		throw - 1;
	}
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name)
	{
		if (name == "tau") return tau;
		throw - 1;
	}
	void Copy(const Method * other)
	{
		const ILUC_preconditioner * b = dynamic_cast<const ILUC_preconditioner *>(other);
		assert(b != NULL);
		tau = b->tau;
		Alink = b->Alink;
		info = b->info;
		sciters = b->sciters;
		eps = b->eps;
	}
	ILUC_preconditioner(const ILUC_preconditioner & other)
	{
		Copy(&other);
	}
	ILUC_preconditioner & operator =(ILUC_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	ILUC_preconditioner(Solver::OrderInfo & info)
		:tau(DEFAULT_TAU),info(&info)
	{
		Alink = NULL;
		init = false;
		sciters = 8;
		eps = 1e-54;
	}
	bool isInitialized() { return init; }
	bool isFinalized() { return !init; }
	bool Initialize()
	{
		const INMOST_DATA_REAL_TYPE tol_modif = 1e-12;
		const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;
		INMOST_DATA_ENUM_TYPE wbeg, wend; //working interval
		INMOST_DATA_ENUM_TYPE mobeg, moend; // total interval
		INMOST_DATA_ENUM_TYPE k, i, j, Li, Ui, curr, next, mi, mj;
		INMOST_DATA_REAL_TYPE l,u,udiag, max_diag, min_diag;
#if defined(SCHUR_DROPPING)
		INMOST_DATA_REAL_TYPE Emeannorm, Fmeannorm, Enorm, Fnorm, Dnorm, drop, Cmax;
#endif
		INMOST_DATA_ENUM_TYPE nzA, nzS, nzA2, nzLU = 0, nzEF = 0, nzL, nzU;
		Solver::Vector DL, DR;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		//prepare temporal array
		temp.set_interval_beg(mobeg);
		temp.set_interval_end(moend);
		std::fill(temp.begin(), temp.end(), 0.0);
		//prepare reordering vectors
		ddP.set_interval_beg(mobeg);
		ddP.set_interval_end(moend);
		ddQ.set_interval_beg(mobeg);
		ddQ.set_interval_end(moend);
		//prepare rescaling vectors
		DL.SetInterval(mobeg, moend);
		DR.SetInterval(mobeg, moend);
		std::fill(DL.Begin(), DL.End(), 0.0);
		for (k = mobeg; k != moend; k++) ddP[k] = ddQ[k] = k;
		//preinitialize data
		interval<INMOST_DATA_ENUM_TYPE, row_col> LU; //remember blocks of ILU decompositions, actual sizes are in level_size array
		LU.set_interval_beg(mobeg);
		LU.set_interval_end(moend);
		level_size.clear();
		// supplementary data for ddPQ reordering
		interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(mobeg, moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(mobeg, moend), invQ(mobeg,moend), localP(mobeg, moend), localQ(mobeg,moend);
#if defined(REORDER_DDPQ)
		interval<INMOST_DATA_ENUM_TYPE, wgt_pair> wgts(mobeg, moend);
		wgt_coords sort_wgts(2 * (moend - mobeg));
#endif
		//supplimentary data structures for ILUC
		INMOST_DATA_ENUM_TYPE LU_Beg;
		U_Address.set_interval_beg(mobeg);
		U_Address.set_interval_end(moend);
		L_Address.set_interval_beg(mobeg);
		L_Address.set_interval_end(moend);
		LU_Diag.set_interval_beg(mobeg);
		LU_Diag.set_interval_end(moend);
		std::vector<Solver::Row::entry> A_Entries, B_Entries, C_Entries, S_Entries;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend), B_Address(mobeg, moend), C_Address(mobeg,moend), S_Address(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ulist(mobeg, moend,UNDEF), Llist(mobeg, moend,UNDEF), Blist(mobeg,moend,UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(mobeg, moend,EOL), Lbeg(mobeg, moend,EOL), Bbeg(mobeg,moend,EOL);

		

		//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
		INMOST_DATA_REAL_TYPE mup, mum, NuU, NuL, NuD, vp, vm, v;
		INMOST_DATA_ENUM_TYPE np, nm;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL(mobeg, moend,0.0), EstU(mobeg, moend,0.0);
		//supplimentary data structures for returning values of dropped elements
		//INMOST_DATA_REAL_TYPE DropLk, DropUk;
		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
		//data structure for linked list
		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValues(mobeg, moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndeces(mobeg, moend+1,UNDEF);
		double tfactor = 0.0, trescale = 0.0, tschur = 0.0, treorder = 0.0, treassamble = 0.0, ttotal, tt;
		ttotal = Timer();
		
		//calculate number of nonzeros
		nzA = 0;
		for (k = mobeg; k < moend; ++k)
		{
			for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0) nzA++;
		}

		
		A_Entries.resize(nzA);
		B_Entries.reserve(nzA);
		C_Entries.reserve(nzA);
		S_Entries.reserve(nzA);
		LU_Entries.reserve(nzA);

		j = 0;
		for (k = mobeg; k < moend; ++k)
		{
			A_Address[k].first = j;
			for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
			{
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
					A_Entries[j++] = Solver::Row::entry(r->first, r->second);
			}
			A_Address[k].last = j;
			assert(A_Address[k].Size() != 0); //singular matrix
		}

#if defined(ILUC2)
		INMOST_DATA_ENUM_TYPE nzLU2 = 0;
		INMOST_DATA_REAL_TYPE tau2 = ILUC2_TAU;
		std::vector<Solver::Row::entry> LU2_Entries;
		interval<INMOST_DATA_ENUM_TYPE, Interval> L2_Address(mobeg, moend+1), U2_Address(mobeg, moend+1);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2list(mobeg, moend, UNDEF), L2list(mobeg, moend, UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2beg(mobeg, moend, EOL), L2beg(mobeg, moend, EOL);
		LU2_Entries.reserve(nzA * 4);
#endif

#if defined(REPORT_ILU)
		std::cout << "nonzeros in A " << nzA << std::endl;
#endif
#if defined(REORDER_DDPQ)
		wgt init_wgt;
		wgt_pair init_wgts;
		init_wgt.ind = ENUMUNDEF;
		init_wgt.max = 0.0;
		init_wgt.nnz = 0;
		init_wgt.sum = 0;
		init_wgts.row = init_wgt;
		init_wgts.col = init_wgt;
#endif
		//set up working interval
		wbeg = mobeg;
		wend = moend;
		while (wbeg < wend) //iterate into levels until all matrix is factored
		{
			//ddPQ reordering on current Schur complement
			INMOST_DATA_ENUM_TYPE cbeg = wbeg, cend = wend; //next size of factored B block
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// REORDERING                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
			std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
			std::fill(localQ.begin() + (wbeg - mobeg), localQ.begin() + (wend - mobeg), ENUMUNDEF);

#if defined(REORDER_DDPQ)

#if defined(REPORT_ILU)
			std::cout << level_size.size() << " calculate weights" << std::endl;
#endif
			
			std::fill(wgts.begin()+wbeg-mobeg, wgts.begin()+wend-mobeg, init_wgts);
			for (k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
				{
					i = A_Entries[it].first;
					u = fabs(A_Entries[it].second);
					wgts[k].row.nnz++;
					wgts[k].row.sum += u;
					if (wgts[k].row.max < u)
					{
						wgts[k].row.max = u;
						wgts[k].row.ind = i;
					}
					wgts[i].col.nnz++;
					wgts[i].col.sum += u;
					if (wgts[i].col.max < u)
					{
						wgts[i].col.max = u;
						wgts[i].col.ind = k;
					}
				}
			}
			//first select those elements that are largest over row and column simultaneously
			sort_wgts.clear();
			for (k = wbeg; k < wend; ++k)
			{
				if( wgts[k].row.ind != ENUMUNDEF )
					sort_wgts.push_back(wgt_coord(wgts[k].row.sum  * (wgts[k].row.nnz + wgts[wgts[k].row.ind].col.nnz) / wgts[k].row.max / (1.0 + (wgts[wgts[k].row.ind].col.ind == k ? 1.0 : 0.0)) , coord(k, wgts[k].row.ind)));
			}
			std::sort(sort_wgts.begin(), sort_wgts.end());
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " fill reordering" << std::endl;
#endif
			i = wbeg;
			for (wgt_coords::iterator it = sort_wgts.begin(); it != sort_wgts.end() && it->first < sort_wgts.front().first * DDPQ_TAU_INV; ++it)
			{
				if (localP[it->second.first] == ENUMUNDEF && localQ[it->second.second] == ENUMUNDEF)
				{
					localP[it->second.first] = i;
					localQ[it->second.second] = i;
					i++;
				}
			}
			cend = i;
#else
			cend = wend;
			i = wbeg;
#endif
			if (cbeg == cend)
			{
				std::cout << __FILE__ << ":" << __LINE__ << " singular result" << std::endl;
				break;
			}
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " new level " << cbeg << ".." << cend << " out of " << wbeg << ".." << wend << std::endl;
#endif
			//finish reordering
			if (i < wend)
			{
				j = i;
				for (k = wbeg; k < wend; ++k)
				{
					if (localP[k] == ENUMUNDEF) localP[k] = i++;
					if (localQ[k] == ENUMUNDEF) localQ[k] = j++;
				}
			}
			treorder += Timer() - tt;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// REASSAMBLE             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
			double tt1, tt2, tt3,ttt;
			//in localPQ numbers indicate where to put current row/column
			//reorder E,F blocks by swaps
			tt1 = Timer();
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " applying permutations to EF blocks " << std::endl;
#endif
			//allocate space for E,F blocks
			if (level_size.empty())
			{
				EF.set_interval_beg(cend);
				EF.set_interval_end(wend);
				Efirst.set_interval_beg(cend);
				Efirst.set_interval_end(wend);
				Ffirst.set_interval_beg(cend);
				Ffirst.set_interval_end(wend);
			}
			else ReorderEF(mobeg, cbeg, wend, donePQ, localP,localQ);
			//update E and F indices, so that next time we access only needed values
			for (k = cend; k < wend; ++k)
			{
				Efirst[k] = EF[k].row.Size();
				Ffirst[k] = EF[k].col.Size();
			}
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			tt1 = Timer() - tt1;
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " rebuilding matrix, B,F blocks " << std::endl;
#endif
			tt2 = Timer();
			//Create B block from current Schur complement for ILUC
			//      | B   F |
			//  A = |       |
			//      | E   C |
			//compute B,F
			for (k = cbeg; k < cend; ++k)
			{
				bool meet_diag = false;
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt != A_Address[invP[k]].last; ++jt)
				{
					i = localQ[A_Entries[jt].first];
					u = A_Entries[jt].second;
					if (i == k) //put diagonal element to diagonal of B
					{
						meet_diag = true;
						LU[k].diag = u;
					}
					else if (i > k)
					{
						assert(i < wend);
						if( i < cend ) LU[k].row.Push(i, u); //put to row of B
						else EF[i].col.Push(k, u); //put to column of F
					}
					else LU[i].col.Push(k, u); //put to column of B
				}
				std::sort(LU[k].row.Begin(), LU[k].row.End());
				if (!meet_diag) LU[k].diag = 0.0;
				//assert(meet_diag);
#if defined(DIAGONAL_PERTURBATION)
				LU[k].diag = LU[k].diag * (1.0 + DIAGONAL_PERTURBATION_REL) + ((LU[k].diag < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS);
#endif
			}
			for (k = cend; k < wend; ++k) nzEF += EF[k].col.Size() - Ffirst[k];

			tt2 = Timer() - tt2;
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " rebuilding matrix, E,C blocks " << std::endl;
#endif

			tt3 = Timer();
			//compute E,C
			{
				C_Entries.clear();
				for (k = cend; k < wend; ++k)
				{
					C_Address[k].first = C_Entries.size();
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt != A_Address[invP[k]].last; ++jt)
					{
						i = localQ[A_Entries[jt].first];
						u = A_Entries[jt].second;
						if (i < cend) EF[k].row.Push(i, u); //put to row of E
						else C_Entries.push_back(Solver::Row::entry(i, u)); //form new Schur complement
					}
					C_Address[k].last = C_Entries.size();
					// sort entries added to E, because it will be very usefull later
					std::sort(EF[k].row.Begin() + Efirst[k], EF[k].row.End());
					nzEF += EF[k].row.Size() - Efirst[k];
				}
				Cmax = 0.0;
				for (k = 0; k < C_Entries.size(); ++k) if (Cmax < fabs(C_Entries[k].second)) Cmax = fabs(C_Entries[k].second);
			}
			tt3 = Timer() - tt3;
			ttt = Timer();
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);
			tt1 += Timer() - ttt;

#if defined(REPORT_ILU)
			std::cout << level_size.size() << " finished permutation, E,F " << tt1 << " B,F " << tt2 << " E,C " << tt3 << std::endl;
#endif
			treassamble += Timer() - tt;			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING                 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Rescale current B block by Row-Column alternating scaling
			tt = Timer();
#if defined(RESCALE_B)
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " rescaling block B " << std::endl;
#endif
			const INMOST_DATA_REAL_TYPE subst = 1.0;
			for (k = cbeg; k < cend; k++)
			{
				DL[k] += LU[k].diag*LU[k].diag;
				for (Solver::Row::iterator r = LU[k].row.Begin(); r != LU[k].row.End(); ++r) DL[k] += r->second*r->second;
				for (Solver::Row::iterator r = LU[k].col.Begin(); r != LU[k].col.End(); ++r) DL[r->first] += r->second*r->second;
			}
			for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					DR[k] += DL[k] * LU[k].diag*LU[k].diag;
					for (Solver::Row::iterator r = LU[k].row.Begin(); r != LU[k].row.End(); ++r) DR[r->first] += DL[k] * r->second*r->second;
					for (Solver::Row::iterator r = LU[k].col.Begin(); r != LU[k].col.End(); ++r) DR[k] += DL[r->first] * r->second*r->second;
				}
				for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					DL[k] += DR[k] * LU[k].diag*LU[k].diag;
					for (Solver::Row::iterator r = LU[k].row.Begin(); r != LU[k].row.End(); ++r) DL[k] += DR[r->first] * r->second*r->second;
					for (Solver::Row::iterator r = LU[k].col.Begin(); r != LU[k].col.End(); ++r) DL[r->first] += DR[k] * r->second*r->second;
				}
				for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}
			for (k = cbeg; k < cend; k++) DL[k] = sqrt(DL[k]);
			for (k = cbeg; k < cend; k++) DR[k] = sqrt(DR[k]);
			for (k = cbeg; k < cend; k++)
			{
				LU[k].diag *= DL[k] * DR[k];
				for (Solver::Row::iterator r = LU[k].row.Begin(); r != LU[k].row.End(); ++r) r->second = r->second * DL[k] * DR[r->first];
				for (Solver::Row::iterator r = LU[k].col.Begin(); r != LU[k].col.End(); ++r) r->second = r->second * DL[r->first] * DR[k];
			}
#endif
			//End rescale B block
			trescale += Timer() - tt;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// FACTORIZATION BEGIN ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
#if defined(ILUC2)
				nzLU2 = 0;
				LU2_Entries.clear();
				U2_Address[cbeg].first = U2_Address[cbeg].last = 0;
				L2_Address[cbeg].first = L2_Address[cbeg].last = 0;
#endif
				//initialize first row and column
				k = cbeg;
				LU_Diag[k] = LU[k].diag;
#if defined(PIVOT_THRESHOLD)
				if (fabs(LU_Diag[k]) < tol_modif)
				{
					LU_Diag[k] = LU_Diag[k] < 0.0 ? -tol_modif : tol_modif;
				}
#endif
				LU_Beg = LU_Entries.size();
				U_Address[cbeg].first = LU_Beg;
				for (Solver::Row::iterator r = LU[k].row.Begin(); r != LU[k].row.End(); ++r) LU_Entries.push_back(Solver::Row::entry(r->first,r->second / LU_Diag[cbeg]));
				L_Address[cbeg].first = LU_Entries.size();
				U_Address[cbeg].last = LU_Entries.size();
				for (Solver::Row::iterator r = LU[k].col.Begin(); r != LU[k].col.End(); ++r) LU_Entries.push_back(Solver::Row::entry(r->first, r->second / LU_Diag[cbeg]));
				L_Address[cbeg].last = LU_Entries.size();
				//form partial list from this row
				//NEW
				if (U_Address[cbeg].Size() > 0)
				{
					i = LU_Entries[U_Address[cbeg].first].first;
					Ubeg[i] = cbeg;
					Ulist[cbeg] = EOL;
				}
				if (L_Address[cbeg].Size() > 0)
				{
					i = LU_Entries[L_Address[cbeg].first].first;
					Lbeg[i] = cbeg;
					Llist[cbeg] = EOL;
				}
				//Initialize data for condition estimator
				tt = Timer();
				NuU = NuL = 1.0;
#if defined(ESTIMATOR)
				EstU[cbeg] = EstL[cbeg] = 0.0;
				for (INMOST_DATA_ENUM_TYPE it = U_Address[cbeg].first; it != U_Address[cbeg].last; ++it) EstU[LU_Entries[it].first] = LU_Entries[it].second;
				for (INMOST_DATA_ENUM_TYPE it = L_Address[cbeg].first; it != L_Address[cbeg].last; ++it) EstL[LU_Entries[it].first] = LU_Entries[it].second;
#endif
				//Factorize the rest
				nzLU += LU[cbeg].col.Size() + LU[cbeg].row.Size() + 1;
				max_diag = min_diag = fabs(LU[cbeg].diag);
				NuD = 1.0;
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " starting factorization " << std::endl;
#endif
				for (k = cbeg + 1; k < cend; ++k)
				{
					//DropLk = DropUk = 0.0;
					//uncompress k-th row
					// add diagonal value first, there shouldn't be values on left from diagonal
					LineIndeces[cbeg] = k;
					LineValues[k] = LU[k].diag;
					Ui = k;
					for (Solver::Row::iterator it = LU[k].row.Begin(); it != LU[k].row.End(); ++it)
					{
						LineValues[it->first] = it->second;
						Ui = LineIndeces[Ui] = it->first;
					}
					LineIndeces[Ui] = EOL;
					// perform the elimination from row
					i = Lbeg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU_Entries[L_Address[i].first].first == k);
						l = LU_Entries[L_Address[i].first].second*LU_Diag[i];
						curr = cbeg;
						for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
						{
							j = LU_Entries[it].first;
							u = LU_Entries[it].second;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= l*u;
							else //add new entry
							{
								next = curr;
								while (next < j)
								{
									curr = next;
									next = LineIndeces[curr];
								}
								assert(curr < j);
								assert(j < next);
								LineIndeces[curr] = j;
								LineIndeces[j] = next;
								LineValues[j] = -l*u;
							}
						}
#if defined(ILUC2)
						if (fabs(l) > tau * NuU)
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								u = LU2_Entries[it].second;
								//eliminate values
								if (LineIndeces[j] != UNDEF) //there is an entry in the list
									LineValues[j] -= l*u;
								else //add new entry
								{
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndeces[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndeces[curr] = j;
									LineIndeces[j] = next;
									LineValues[j] = -l*u;
								}
							}
						}
#endif
						i = Llist[i];
					}
#if defined(ILUC2)
					i = L2beg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU2_Entries[L2_Address[i].first].first == k);
						l = LU2_Entries[L2_Address[i].first].second*LU_Diag[i];
						curr = cbeg;
						for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it) 
						{
							u = LU_Entries[it].second;
							j = LU_Entries[it].first;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= l*u;
							else if (fabs(u) > tau * NuU)//add new entry
							{
								next = curr;
								while (next < j)
								{
									curr = next;
									next = LineIndeces[curr];
								}
								assert(curr < j);
								assert(j < next);
								LineIndeces[curr] = j;
								LineIndeces[j] = next;
								LineValues[j] = -l*u;
							}
						}
						i = L2list[i];
					}
					U2_Address[k].first = LU2_Entries.size();
#endif
					//std::cout << std::endl;
					//define the diagonal value
					udiag = LineValues[k];
#if defined(PIVOT_THRESHOLD)
					if (fabs(udiag) < tol_modif)
					{
						//std::cout << "udiag too small " << 1.0 /udiag << std::endl;
						udiag = udiag < 0.0 ? -tol_modif : tol_modif;
					}
#endif
					//if (1.0 / fabs(udiag) > 1000) std::cout << "udiag is big " << 1.0 / udiag << " " << k << std::endl;
					if (fabs(udiag) > max_diag) max_diag = fabs(udiag);
					if (fabs(udiag) < min_diag) min_diag = fabs(udiag);
					NuD = max_diag / min_diag;
					//rescale elements next after diagonal position
					Ui = LineIndeces[k];
					nzU = 0;
					while (Ui != EOL)
					{
						LineValues[Ui] /= udiag;
						Ui = LineIndeces[Ui];
						nzU++;
					}
					//estimate condition for U^{-1}
#if defined(ESTIMATOR)
					mup = 1.0 - EstU[k];
					mum = -1.0 - EstU[k];
					np = nm = 0;
					//start from the element next after diagonal position
					Ui = LineIndeces[k];
					while (Ui != EOL)
					{
						v = EstU[Ui];
						vp = fabs(v + LineValues[Ui] * mup);
						vm = fabs(v + LineValues[Ui] * mum);
						v = fabs(v);
						if (vp > std::max(2 * v, 0.5)) np++;
						if (std::max(2 * vp, 0.5) < v) np--;
						if (vm > std::max(2 * v, 0.5)) nm++;
						if (std::max(2 * vm, 0.5) < v) nm--;
						Ui = LineIndeces[Ui];
					}
					if (np > nm) NuU = mup; else NuU = mum;
					Ui = LineIndeces[k];
					while (Ui != EOL)
					{
						EstU[Ui] += LineValues[Ui] * NuU;
						Ui = LineIndeces[Ui];
					}
					NuU = std::max(fabs(mup), fabs(mum));
#endif
					//insert line to U part
					U_Address[k].first = LU_Entries.size();
					Ui = LineIndeces[k];
					while (Ui != EOL)
					{
						u = fabs(LineValues[Ui])*NuU;
						if (u > tau) // apply dropping rule
							LU_Entries.push_back(Solver::Row::entry(Ui, LineValues[Ui]));
#if defined(ILUC2)
						else if (u > tau2)
							LU2_Entries.push_back(Solver::Row::entry(Ui,LineValues[Ui]));
#endif
						Ui = LineIndeces[Ui];
					}
					U_Address[k].last = LU_Entries.size();

					//Clean after use
					Ui = cbeg;
					while (Ui != EOL)
					{
						i = LineIndeces[Ui];
						LineValues[Ui] = 0.0; // clean values after use
						LineIndeces[Ui] = UNDEF; //clean indeces after use
						Ui = i;
					}
					if (U_Address[k].Size() != 0)
					{
						i = LU_Entries[U_Address[k].first].first;
						assert(i > k);
						Ulist[k] = Ubeg[i];
						Ubeg[i] = k;
					}
#if defined(ILUC2)
					U2_Address[k].last = LU2_Entries.size();
					if (U2_Address[k].Size() != 0)
					{
						i = LU2_Entries[U2_Address[k].first].first;
						assert(i > k);
						U2list[k] = U2beg[i];
						U2beg[i] = k;
					}
#endif
					//uncompress column
					//insert diagonal value first
					LineIndeces[cbeg] = k;
					LineValues[k] = LU[k].diag;
					//start from diagonal
					Ui = k;
					//column values must be sorted due to construction algorithm
					//so we don't need to sort indeces
					for (Solver::Row::iterator it = LU[k].col.Begin(); it != LU[k].col.End(); ++it)
					{
						LineValues[it->first] = it->second;
						Ui = LineIndeces[Ui] = it->first;
					}
					LineIndeces[Ui] = EOL;
					// perform the elimination from column
					i = Ubeg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU_Entries[U_Address[i].first].first == k);
						u = LU_Entries[U_Address[i].first].second*LU_Diag[i];
						curr = cbeg;
						for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it < L_Address[i].last; ++it)
						{
							j = LU_Entries[it].first;
							l = LU_Entries[it].second;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= l*u;
							else //add new entry
							{
								next = curr;
								while (next < j)
								{
									curr = next;
									next = LineIndeces[curr];
								}
								assert(curr < j);
								assert(j < next);
								LineIndeces[curr] = j;
								LineIndeces[j] = next;
								LineValues[j] = -l*u;
							}
						}
#if defined(ILUC2)
						if (fabs(u)*NuL > tau)
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								l = LU2_Entries[it].second;
								if (LineIndeces[j] != UNDEF) //there is an entry in the list
									LineValues[j] -= l * u;
								else //add new entry
								{
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndeces[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndeces[curr] = j;
									LineIndeces[j] = next;
									LineValues[j] = -l*u;
								}
							}
						}
#endif
						i = Ulist[i];
					}
#if defined(ILUC2)
					i = U2beg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU2_Entries[U2_Address[i].first].first == k);
						u = LU2_Entries[U2_Address[i].first].second*LU_Diag[i];
						curr = cbeg;
						for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it != L_Address[i].last; ++it)
						{
							j = LU_Entries[it].first;
							l = LU_Entries[it].second;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= l*u;
							else if (fabs(l)*NuL > tau)//add new entry
							{
								next = curr;
								while (next < j)
								{
									curr = next;
									next = LineIndeces[curr];
								}
								assert(curr < j);
								assert(j < next);
								LineIndeces[curr] = j;
								LineIndeces[j] = next;
								LineValues[j] = -u*l;
							}
							
						}
						i = U2list[i];
					}
					L2_Address[k].first = LU2_Entries.size();
#endif
					//check that diagonal value is the same(must be!)
					//assert(fabs(LineValues[k] / udiag - 1.0) < 1.0e-10);
					//rescale line by diagonal
					Li = LineIndeces[k];
					nzL = 0;
					while (Li != EOL)
					{
						nzL++;
						LineValues[Li] /= udiag;
						Li = LineIndeces[Li];
					}
					//estimate condition for L^{-1}
#if defined(ESTIMATOR)
					mup = 1.0 - EstL[k];
					mum = -1.0 - EstL[k];
					np = nm = 0;
					//start from the element next after diagonal position
					Li = LineIndeces[k];
					while (Li != EOL)
					{
						v = EstL[Li];
						vp = fabs(v + LineValues[Li] * mup);
						vm = fabs(v + LineValues[Li] * mum);
						v = fabs(v);
						if (vp > std::max(2 * v, 0.5)) np++;
						if (std::max(2 * vp, 0.5) < v) np--;
						if (vm > std::max(2 * v, 0.5)) nm++;
						if (std::max(2 * vm, 0.5) < v) nm--;
						Li = LineIndeces[Li];
					}
					if (np > nm) NuL = mup; else NuL = mum;
					Li = LineIndeces[k];
					while (Li != EOL)
					{
						EstL[Li] += LineValues[Li] * NuL;
						Li = LineIndeces[Li];
					}
					NuL = std::max(fabs(mup), fabs(mum));
#endif
					//insert column to L part
					Li = LineIndeces[k];
					L_Address[k].first = LU_Entries.size();
					while (Li != EOL)
					{
						u = fabs(LineValues[Li]) * NuL;
						if (u > tau) //apply dropping 
							LU_Entries.push_back(Solver::Row::entry(Li, LineValues[Li]));
#if defined(ILUC2)
						else if (u > tau2)
							LU2_Entries.push_back(Solver::Row::entry(Li,LineValues[Li]));
#endif
						Li = LineIndeces[Li];
					}
					L_Address[k].last = LU_Entries.size();
					Li = cbeg;
					while (Li != EOL)
					{
						i = LineIndeces[Li];
						LineValues[Li] = 0.0; //clean values after use
						LineIndeces[Li] = UNDEF; //clean indeces after use
						Li = i;
					}
					if (L_Address[k].Size() > 0)
					{
						i = LU_Entries[L_Address[k].first].first;
						assert(i > k);
						Llist[k] = Lbeg[i];
						Lbeg[i] = k;
					}
#if defined(ILUC2)
					L2_Address[k].last = LU2_Entries.size();
					if (L2_Address[k].Size() != 0)
					{
						i = LU2_Entries[L2_Address[k].first].first;
						assert(i > k);
						L2list[k] = L2beg[i];
						L2beg[i] = k;
					}
#endif
					LU_Diag[k] = udiag;
					//Update indeces in linked lists
					i = Ubeg[k];
					while (i != EOL)
					{
						U_Address[i].first++;
						Li = Ulist[i];
						if (U_Address[i].Size() > 0)
						{
							Ui = LU_Entries[U_Address[i].first].first;
							Ulist[i] = Ubeg[Ui];
							Ubeg[Ui] = i;
						}
						i = Li;
					}
					i = Lbeg[k];
					while (i != EOL)
					{
						L_Address[i].first++;
						Li = Llist[i];
						if (L_Address[i].Size() > 0)
						{
							Ui = LU_Entries[L_Address[i].first].first;
							Llist[i] = Lbeg[Ui];
							Lbeg[Ui] = i;
						}
						i = Li;
					}
#if defined(ILUC2)
					i = U2beg[k];
					while (i != EOL)
					{
						U2_Address[i].first++;
						Li = U2list[i];
						if (U2_Address[i].Size() > 0)
						{
							Ui = LU2_Entries[U2_Address[i].first].first;
							U2list[i] = U2beg[Ui];
							U2beg[Ui] = i;
						}
						i = Li;
					}
					i = L2beg[k];
					while (i != EOL)
					{
						L2_Address[i].first++;
						Li = L2list[i];
						if (L2_Address[i].Size())
						{
							Ui = LU2_Entries[L2_Address[i].first].first;
							L2list[i] = L2beg[Ui];
							L2beg[Ui] = i;
						}
						i = Li;
					}
#endif
					//number of nonzeros
					nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
					nzLU2 += U2_Address[k].Size() + L2_Address[k].Size();
#if defined(REPORT_ILU)
					if (k % 1000 == 0)
					{

						printf("%6.2f%% nnz LU %8d LU2 %8d condition L %10f D %10f U %10f\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, nzLU2, NuL, NuD, NuU);
						fflush(stdout);

					}
#elif defined(REPORT_ILU_PROGRESS)
					if (k % 1000 == 0)
					{
						printf("level %d factor %6.2f%%\t\t\r", level_size.size(), 100.0f*(k - cbeg) / (float)(cend - cbeg));
						fflush(stdout);
					}
#endif
				}
#if defined(REPORT_ILU)
				std::cout << "new level " << cend - cbeg << " total nonzeros in LU " << nzLU;
#if defined(ILUC2)
				std::cout << " in LU2 " << nzLU2;
#endif
				std::cout << " in EF " << nzEF << " conditions L " << NuL << " D " << NuD << " U " << NuU << std::endl;
#endif
			}
			tfactor += Timer() - tt;
			//restore indexes
			U_Address[cbeg].first = LU_Beg;
			L_Address[cbeg].first = U_Address[cbeg].last;
			for (k = cbeg+1; k < cend; ++k)
			{
				U_Address[k].first = L_Address[k - 1].last;
				L_Address[k].first = U_Address[k].last;
			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// FACTORIZATION COMPLETE ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//After we have factored the rescaled matrix we must rescale obtained factors
			tt = Timer();

#if defined(RESCALE_B)
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " rescaling block B back " << std::endl;
#endif
			for (k = cbeg; k < cend; ++k)
			{
				LU_Diag[k] /= DL[k] * DR[k];
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k].first; r < U_Address[k].last; ++r) LU_Entries[r].second *= DR[k] / DR[LU_Entries[r].first];
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r) LU_Entries[r].second *= DL[k] / DL[LU_Entries[r].first];
			}
#endif
			trescale += Timer() - tt;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING DONE   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			level_size.push_back(cend - cbeg);
			//Now we have to calculate new Schur complement
			// S = A - E U^{-1} L^{-1} F
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " computing Schur complement " << std::endl;
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// SCHUR COMPUTATION BEGIN ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
			double t0, t1, t2;
			t0 = Timer();
#if defined(SCHUR_DROPPING)
			Dnorm = 0.0;
			for (k = cbeg; k < cend; ++k)
			{
				Dnorm += LU_Diag[k]*LU_Diag[k];
			}
			Dnorm = sqrt(Dnorm);
			//estimate maximum of schur complement
			Emeannorm = Fmeannorm = 0.0;
			for (k = cend; k < wend; ++k)
			{
				Enorm = Fnorm = 0.0;
				for (j = Efirst[k]; j < EF[k].row.Size(); j++)
					Enorm += EF[k].row.GetValue(j)*EF[k].row.GetValue(j);
				for (j = Ffirst[k]; j < EF[k].col.Size(); j++)
					Fnorm += EF[k].col.GetValue(j)*EF[k].col.GetValue(j);
				Emeannorm += Enorm;
				Fmeannorm += Fnorm;
			}
			Emeannorm = sqrt(Emeannorm);
			Fmeannorm = sqrt(Fmeannorm);
			drop = Cmax / (Emeannorm*NuU*Dnorm*NuL*Emeannorm)*SCHUR_TAU;
			std::cout << "Emeannorm " << Emeannorm << " Fmeannorm " << Fmeannorm << std::endl;
#endif
			
			{
				Solver::Matrix LF;
				LF.SetInterval(cbeg, cend);
				//Compute entire LF matrix with row-major order
				t1 = Timer();
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// LF block    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//calculate drop for values
				
				for (i = cend; i < wend; ++i) //iterate over columns of F
				{
					//unpack row to linked list
					Li = cbeg;
					for (j = Ffirst[i]; j < EF[i].col.Size(); ++j) //iterate over values of i-th column of F
					{
						LineValues[EF[i].col.GetIndex(j)] = EF[i].col.GetValue(j);
						Li = LineIndeces[Li] = EF[i].col.GetIndex(j) + 1;
					}
					LineIndeces[Li] = EOL;

					//compute ~F_i = L^{-1} F_i
					Li = LineIndeces[cbeg];
					while (Li != EOL)
					{
						curr = Li;
						
						for (INMOST_DATA_ENUM_TYPE ru = L_Address[Li - 1].first; ru < L_Address[Li - 1].last; ++ru)
						{
							u = LineValues[Li - 1] * LU_Entries[ru].second;
							j = LU_Entries[ru].first;
							if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
								LineValues[j] -= u;
							else 
#if defined(SCHUR_DROPPING)
							if (fabs(u) > drop)
#endif
							{
								next = curr;
								while (next < j + 1)
								{
									curr = next;
									next = LineIndeces[curr];
								}
								assert(curr < j + 1);
								assert(j + 1 < next);
								LineIndeces[curr] = j + 1;
								LineIndeces[j + 1] = next;
								LineValues[j] = -u;
							}
						}
						Li = LineIndeces[Li];
					}
					
					//Assemble column into matrix
					Li = LineIndeces[cbeg];
					while (Li != EOL)
					{
#if defined(SCHUR_DROPPING)
						if (fabs(LineValues[Li - 1]) > drop)
#endif
							LF[Li - 1].Push(i, LineValues[Li - 1]);
						Li = LineIndeces[Li];
					}
					//Clean after use
					Li = LineIndeces[cbeg];
					while (Li != EOL)
					{
						j = LineIndeces[Li];
						LineIndeces[Li] = UNDEF;
						LineValues[Li - 1] = 0.0;
						Li = j;
					}
#if defined(REPORT_ILU_PROGRESS)
					if (i % 10 == 0)
					{
						printf("%d schur1 %6.2f%%\t\t\r", level_size.size(), 100.f*(i - cend) / (1.f*(wend - cend)));
						fflush(stdout);
					}
#endif
				}
				t1 = Timer() - t1;
				
				
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// EU and Schur blocks  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				t2 = Timer();
				std::fill(Ulist.begin() + cend - mobeg, Ulist.begin() + wend - mobeg, UNDEF);
				S_Entries.clear();
				for (i = cend; i < wend; ++i) //iterate over rows of E
				{
					double tt0, tt1, tt2;
					tt0 = Timer();
					// uncompress sprase i-th row of E to linked list
					Li = cbeg;
					for (j = Efirst[i]; j < EF[i].row.Size(); ++j) //iterate over values of i-th row of E
					{
						LineValues[EF[i].row.GetIndex(j)] = EF[i].row.GetValue(j);
						Li = LineIndeces[Li] = EF[i].row.GetIndex(j) + 1;
					}
					LineIndeces[Li] = EOL;
					// compute ~E_i = E_i U^{-1}
					Li = LineIndeces[cbeg];
					while (Li != EOL)
					{
						curr = Li;
						for (INMOST_DATA_ENUM_TYPE ru = U_Address[Li - 1].first; ru != U_Address[Li - 1].last; ++ru)
						{
							u = LineValues[Li - 1] * LU_Entries[ru].second;
							j = LU_Entries[ru].first;
							if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
								LineValues[j] -= u;
							else 
#if defined(SCHUR_DROPPING)
							if (fabs(u) > drop)
#endif
							{
								next = curr;
								while (next < j + 1)
								{
									curr = next;
									next = LineIndeces[curr];
								}
								assert(curr < j + 1);
								assert(j + 1 < next);
								LineIndeces[curr] = j + 1;
								LineIndeces[j + 1] = next;
								LineValues[j] = -u;
							}
						}
						Li = LineIndeces[Li];
					}
					
#if defined(SCHUR_DROPPING)
					//drop values
					Li = LineIndeces[cbeg];
					Ui = cbeg;
					while (Li != EOL)
					{
						j = LineIndeces[Li];
						if (fabs(LineValues[Li - 1]) < drop)
						{
							LineIndeces[Ui] = j;
							LineIndeces[Li] = UNDEF;
							LineValues[Li - 1] = 0.0;
						}
						else Ui = Li;
						Li = j;
					}
#endif
					//rescale vector by diagonal *E = ~E_i D^{-1} = E_i U^{-1} D^{-1}
					Li = LineIndeces[cbeg];
					while (Li != EOL)
					{
						LineValues[Li - 1] /= LU_Diag[Li - 1];
						Li = LineIndeces[Li];
					}
					//unpack row of C
					tt0 = Timer() - tt0;

					tt1 = Timer();
					Ubeg[cbeg] = EOL;
					for (INMOST_DATA_ENUM_TYPE r = C_Address[i].first; r != C_Address[i].last; ++r)
					{
						temp[C_Entries[r].first] = C_Entries[r].second;
						Ulist[C_Entries[r].first] = Ubeg[cbeg];
						Ubeg[cbeg] = C_Entries[r].first;
					}
					//multiply current row of EU by LF and add to list
					Li = LineIndeces[cbeg];
					while (Li != EOL)
					{
						//iterate over corresponding row of LF, add multiplication to list
						for (Solver::Row::iterator r = LF[Li - 1].Begin(); r != LF[Li - 1].End(); ++r)
						{
							u = LineValues[Li - 1] * r->second;
							if (Ulist[r->first] != UNDEF)
								temp[r->first] -= u;
							else 
#if defined(SCHUR_DROPPING)
							if (fabs(u) > drop)
#endif
							{
								temp[r->first] = -u;
								Ulist[r->first] = Ubeg[cbeg];
								Ubeg[cbeg] = r->first;
							}
						}
						Li = LineIndeces[Li];
					}
					tt1 = Timer() - tt1;

					tt2 = Timer();
					//insert obtained row
					Ui = Ubeg[cbeg];
					S_Address[i].first = S_Entries.size();
					while (Ui != EOL)
					{
						//drop values
						S_Entries.push_back(Solver::Row::entry(Ui, temp[Ui]));
						Ui = Ulist[Ui];
					}
					S_Address[i].last = S_Entries.size();
					//clean after use
					Ui = Ubeg[cbeg];
					while (Ui != EOL)
					{
						j = Ulist[Ui];
						temp[Ui] = 0.0;
						Ulist[Ui] = UNDEF;
						Ui = j;
					}
					Li = LineIndeces[cbeg];
					while (Li != EOL)
					{
						j = LineIndeces[Li];
						LineIndeces[Li] = UNDEF;
						LineValues[Li - 1] = 0.0;
						Li = j;
					}
					tt2 = Timer() - tt2;
#if defined(REPORT_ILU)
					if (i % 10 == 0)
					{
						printf("%6.2f%%  solve %f multiply %f assemble %f\r", 100.f*(i - cend) / (1.f*(wend - cend)), tt0,tt1,tt2);
						fflush(stdout);
					}
#elif defined(REPORT_ILU_PROGRESS)
					if (i % 100 == 0)
					{
						printf("%d schur2 %6.2f%%\t\t\r", level_size.size(),100.f*(i - cend) / (1.f*(wend - cend)));
						fflush(stdout);
					}
#endif
				}
				t2 = Timer() - t2;
			}
			A_Entries.swap(S_Entries);
			A_Address.swap(S_Address);
			t0 = Timer() - t0;
#if defined(REPORT_ILU)
			printf("\ntotal %f phase1 %f (%6.2f%%) phase2 %f (%6.2f%%)\n", t0, t1, 100.0*t1 / t0, t2, 100.0*t2 / t0);
#endif
			tschur += Timer() - tt;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// SCHUR COMPLETE ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//renew working interval
			wbeg = cend;
			//compute number of nonzeros in Schur complement
			nzA = 0;
			for (k = wbeg; k < wend; ++k) nzA += A_Address[k].Size();
			//check that Schur complement is too small or too dense
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// DENSE FACTORIZATION ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
			if ( (wend-wbeg > 0 && wend - wbeg < 32) || (1.0 * nzA / ((wend - wbeg)*(wend - wbeg))) > 0.75)
			//if ( wend-wbeg > 0 && wend - wbeg < 2 )
			{
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " starting dense factorization "  << wend - wbeg << std::endl;
#endif
				for (k = wbeg; k < wend; ++k)
				{
					localP[k] = k;
					localQ[k] = k;
				}
				//Perform dense LU with complete pivoting based on one-rank updates
				// (read siam_bollhofer)
				//Allocate matrix
				INMOST_DATA_ENUM_TYPE size = wend - wbeg;
				INMOST_DATA_REAL_TYPE * entries = new INMOST_DATA_REAL_TYPE[size*size];
				std::fill(entries, entries + size*size, 0.0);
				//Fill the matrix
				for (k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r != A_Address[k].last; ++r)
						entries[ADDR(k - wbeg,A_Entries[r].first - wbeg)] = A_Entries[r].second;
				}
				for (i = 0; i < size; ++i) //iterate over diagonal elements
				{
					//first step - pivoting, select largest element over entire schur complement
					mi = mj = i; //predefine maximum element
					
					for (Li = i; Li < size; Li++) //iterate rows
					for (Ui = i; Ui < size; Ui++) //iterate columns
					{
						//test for maximum absolute value
						if (fabs(entries[ADDR(Li,Ui)]) > fabs(entries[ADDR(mi,mj)])) 
						{
							//choose maximum absolute value
							mi = Li;
							mj = Ui;
						}
					}
					
					
					//remember the choice to update global ordering
					Li = wbeg + i;
					Ui = wbeg + mi;
					while (localP[Li] != wbeg+i) Li = localP[Li];
					while (localP[Ui] != wbeg+mi) Ui = localP[Ui];
					localP[Li] = wbeg + mi;
					localP[Ui] = wbeg + i;

					Li = wbeg + i;
					Ui = wbeg + mj;
					while (localQ[Li] != wbeg+i) Li = localQ[Li];
					while (localQ[Ui] != wbeg+mj) Ui = localQ[Ui];
					localQ[Li] = wbeg + mj;
					localQ[Ui] = wbeg + i;
					//swap current row in entries
					if( i != mi ) for (Li = 0; Li < size; Li++)
					{
						u = entries[ADDR(i,Li)]; //remember current row value
						entries[ADDR(i,Li)] = entries[ADDR(mi,Li)]; //replace current row value
						entries[ADDR(mi,Li)] = u;
					}
					//swap current column in entries
					if(i != mj ) for (Li = 0; Li < size; Li++)
					{
						l = entries[ADDR(Li, i)]; //remember current column value
						entries[ADDR(Li, i)] = entries[ADDR(Li, mj)]; //replace current column value
						entries[ADDR(Li, mj)] = l;
					}
					//TODO
					//check that diagonal element is not bad
					if (fabs(entries[ADDR(i, i)]) < 1e-15)
					{
						//std::cout << "warning: tiny diagonal!" << std::endl;
						entries[ADDR(i, i)] = entries[ADDR(i, i)] >= 0.0 ? 1.0e-15 : -1.0e-15;
					}
					//assert(fabs(entries[ADDR(i, i)])>1e-20);
					//scale off-diagonal elements by diagonal
					for (Li = i+1; Li < size; Li++) 
					{
						entries[ADDR(i, Li)] /= entries[ADDR(i, i)];
						entries[ADDR(Li, i)] /= entries[ADDR(i, i)];
					}
					//compute Schur complement
					// S = E - l d u
					for (Li = i + 1; Li < size; Li++)
					for (Ui = i + 1; Ui < size; Ui++)
					{
						entries[ADDR(Li, Ui)] -= entries[ADDR(Li, i)] * entries[ADDR(i, i)] * entries[ADDR(i, Ui)];
					}
					//run condition estimator
					//TODO
					//can it be done better on dense factors?
					//estimate condition for U^{-1}
#if defined(ESTIMATOR)
					mup = 1.0 - EstU[wbeg + i];
					mum = -1.0 - EstU[wbeg + i];
					np = nm = 0;
					//start from the element next after diagonal position
					for (Ui = i + 1; Ui < size; Ui++)
					{
						v = EstU[wbeg + Ui];
						vp = fabs(v + entries[ADDR(i,Ui)] * mup);
						vm = fabs(v + entries[ADDR(i,Ui)] * mum);
						v = fabs(v);
						if (vp > std::max(2 * v, 0.5)) np++;
						if (std::max(2 * vp, 0.5) < v) np--;
						if (vm > std::max(2 * v, 0.5)) nm++;
						if (std::max(2 * vm, 0.5) < v) nm--;
					}
					if (np > nm) NuU = mup; else NuU = mum;
					for (Ui = i + 1; Ui < size; Ui++) EstU[wbeg + Ui] += entries[ADDR(i,Ui)] * NuU;
					NuU = std::max(fabs(mup), fabs(mum));
					//estimate condition for L^{-1}
					mup = 1.0 - EstL[wbeg+i];
					mum = -1.0 - EstL[wbeg+i];
					np = nm = 0;
					//start from the element next after diagonal position
					for (Li = i + 1; Li < size; Li++)
					{
						v = EstL[wbeg+Li];
						vp = fabs(v + entries[ADDR(Li,i)] * mup);
						vm = fabs(v + entries[ADDR(Li,i)] * mum);
						v = fabs(v);
						if (vp > std::max(2 * v, 0.5)) np++;
						if (std::max(2 * vp, 0.5) < v) np--;
						if (vm > std::max(2 * v, 0.5)) nm++;
						if (std::max(2 * vm, 0.5) < v) nm--;
					}
					if (np > nm) NuL = mup; else NuL = mum;
					for (Li = i + 1; Li < size; Li++) EstL[wbeg+Li] += entries[ADDR(Li,i)] * NuL;
					NuL = std::max(fabs(mup), fabs(mum));
#endif
				}
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " finishing reorder " << std::endl;
#endif
				//apply PQ reordering to other part of matrices
				ReorderEF(mobeg,wbeg,wend,donePQ,localP,localQ);
				inversePQ(wbeg, wend, localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);
				max_diag = min_diag = fabs(entries[ADDR(0, 0)]);
				// copy computed decomposition to stored ILU, apply dropping rule
				for (i = 0; i < size; i++)
				{
					if (fabs(entries[ADDR(i, i)]) > max_diag) max_diag = fabs(entries[ADDR(i, i)]);
					if (fabs(entries[ADDR(i, i)]) < min_diag) min_diag = fabs(entries[ADDR(i, i)]);
				}
				NuD = max_diag / min_diag;
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " filling factors " << std::endl;
#endif
				for (i = 0; i < size; i++)
				{
					LU_Diag[wbeg + i] = entries[ADDR(i, i)];
					U_Address[wbeg + i].first = LU_Entries.size();
					for (j = i + 1; j < size; j++)
					{
						if (fabs(entries[ADDR(i, j)]) > tau / std::max(1.0, NuU))
						{
							LU_Entries.push_back(Solver::Row::entry(j + wbeg, entries[ADDR(i, j)])); //Add to U
						}
					}
					U_Address[wbeg + i].last = LU_Entries.size();
					L_Address[wbeg + i].first = LU_Entries.size();
					for (j = i + 1; j < size; j++)
					{
						if (fabs(entries[ADDR(j, i)]) > tau / std::max(1.0, NuL))
						{
							LU_Entries.push_back(Solver::Row::entry(j + wbeg, entries[ADDR(j, i)])); //Add to L
						}
					}
					L_Address[wbeg + i].last = LU_Entries.size();
					nzLU += 1 + L_Address[wbeg + i].Size() + U_Address[wbeg + i].Size();
				}
				delete[] entries;
				//Remember the size of the block
				level_size.push_back(wend - wbeg);
#if defined(REPORT_ILU)
				std::cout << "new level " << wend - wbeg << " total nonzeros in LU " << nzLU << " in EF " << nzEF << " conditions L " << NuL << " D " << NuD << " U " << NuU << std::endl;
#endif
				//FINISH!!!
				break;
			}
			tfactor += Timer() - tt;
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// DENSE FACTORIZATION COMPLETE //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		ttotal = Timer() - ttotal;
#if defined(REPORT_ILU)
		printf("total      %f\n",ttotal);
		printf("reorder    %f (%6.2f%%)\n", treorder, 100.0*treorder/ttotal);
		printf("reassamble %f (%6.2f%%)\n", treassamble, 100.0*treassamble / ttotal);
		printf("rescale    %f (%6.2f%%)\n", trescale, 100.0*trescale / ttotal);
		printf("factor     %f (%6.2f%%)\n", tfactor, 100.0*tfactor / ttotal);
		printf("schur      %f (%6.2f%%)\n", tschur, 100.0*tschur / ttotal);
#endif
		init = true;
		return true;
	}
	bool Finalize()
	{
		init = false;
		L_Address.clear();
		U_Address.clear();
		LU_Entries.clear();
		EF.clear();
		return true;
	}
	bool Solve(Solver::Vector & input, Solver::Vector & output)
	{
		assert(&input != &output);
		//
		INMOST_DATA_ENUM_TYPE k, mobeg, moend, wbeg, wend, cbeg, cend, vbeg, vend;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetVectorRegion(vbeg, vend);
		for (k = vbeg; k < mobeg; k++) temp[k] = 0;
		for (k = mobeg; k < moend; k++) temp[k] = input[k];
		for (k = moend; k < vend; k++) temp[k] = 0;
		//the original matrix A was separated by multilevel algorithm to the following form
		//     | B  F |
		// A = |      |
		//     | E  C |
		// In order to apply solve phase on this matrix, we consider the matrix in following sense:
		//     | I         0 | | B   0 | | I    B^{-1} F |
		// A = |             | |       | |               | = L D U
		//     | E B^{-1}  I | | 0   S | | 0           I |
		// where S = C - E B^{-1} F
		// consider the system 
		//  | B  F | | u |   | f |
		//  |      | |   | = |   |
		//  | E  C | | y |   | g |
		// then solution is obtained by steps:
		// 1) ~f = B^{-1} f
		// 2) ~g = g - E ~f
		// 3) y = S^{-1} ~g
		// 4) u = ~f - B^{-1} F y
		//first reorder input as prescribed by P
		// ~y = Pt y
		for (k = mobeg; k < moend; ++k) output[k] = temp[ddP[k]];
		for (k = level_size[0]; k < moend; ++k) Efirst[k] = 0;
		wbeg = mobeg;
		wend = moend;
		//perform recursively first two steps of solve phase
		for (levels_t::iterator it = level_size.begin(); it != level_size.end(); ++it)
		{
			//current B block interval
			cbeg = wbeg;
			cend = wbeg + *it;
			// step 1) obtain ~f
			//apply B^{-1} on f
			//Solve with L first
			for (k = cbeg; k < cend; ++k) //iterator over columns of L
			{
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r != L_Address[k].last; ++r)
					output[LU_Entries[r].first] -= output[k] * LU_Entries[r].second;
			}
			//Solve with diagonal
			for (k = cbeg; k < cend; ++k) output[k] /= LU_Diag[k];
			//Solve with U
			for (k = cend; k > cbeg; --k) //iterator over rows of U
			{
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k-1].first; r != U_Address[k-1].last; ++r)
					output[k-1] -= output[LU_Entries[r].first] * LU_Entries[r].second;
			}
			//now can calculate ~g
			//multiply ~f by row of E
			for (k = cend; k < wend; ++k)
			{
				for (Efirst[k]; Efirst[k] < EF[k].row.Size() && EF[k].row.GetIndex(Efirst[k]) < cend; ++Efirst[k])
					output[k] -= output[EF[k].row.GetIndex(Efirst[k])] * EF[k].row.GetValue(Efirst[k]);
			}
			//goto the next level
			wbeg = cend;
		}
		//on the last recursion level third step was accomplished
		//finish last step in unwinding recursion
		for (k = level_size[0]; k < moend; ++k) Ffirst[k] = EF[k].col.Size();
		for (levels_t::reverse_iterator it = level_size.rbegin(); it != level_size.rend(); ++it)
		{
			wbeg -= *it;
			cbeg = wbeg;
			//calculate Fy, iterate over ~g, write to unused input vector
			std::fill(temp.begin() + cbeg - mobeg, temp.begin() + cend - mobeg,0.0);
			for (k = cend; k < wend; ++k)
			{
				for (Ffirst[k]; Ffirst[k] > 0 && EF[k].col.GetIndex(Ffirst[k] - 1) >= cbeg; --Ffirst[k])
					temp[EF[k].col.GetIndex(Ffirst[k] - 1)] += output[k] * EF[k].col.GetValue(Ffirst[k] - 1);
			}
			//perform solve over calculated vector
			//Solve with L first
			for (k = cbeg; k < cend; ++k) //iterator over columns of L
			{
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r != L_Address[k].last; ++r)
					temp[LU_Entries[r].first] -= temp[k] * LU_Entries[r].second; //r->first alwayse > k
			}
			//Solve with diagonal
			for (k = cbeg; k < cend; ++k) temp[k] /= LU_Diag[k];
			//Solve with U
			for (k = cend; k > cbeg; --k) //iterator over rows of U
			{
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k-1].first; r != U_Address[k-1].last; ++r)
					temp[k-1] -= temp[LU_Entries[r].first] * LU_Entries[r].second; // r->first always > k
			}
			//substract obtained B^{-1} F y from ~f to obtain final solution
			for (k = cbeg; k < cend; ++k) output[k] -= temp[k];
			cend = wbeg;
		}
		//System was solved with reordered matrix
		// PAQ x = y
		// AQ x = Pt y
		// Q x = A^{-1} Pt y <- we have now
		// x = Qt A^{-1} Pt y
		//reorder output by Q
		for (k = mobeg; k < moend; ++k) temp[ddQ[k]] = output[k];
		for (k = mobeg; k < moend; ++k) output[k] = temp[k];
		info->Accumulate(output);
		return true;
	}
	bool ReplaceMAT(Solver::Matrix & A){ if (isInitialized()) Finalize(); Alink = &A; return true; }
	bool ReplaceSOL(Solver::Vector & x) { return true; }
	bool ReplaceRHS(Solver::Vector & b) { return true; }
	Method * Duplicate() { return new ILUC_preconditioner(*this); }
};

#endif //__SOLVER_ILU2__
