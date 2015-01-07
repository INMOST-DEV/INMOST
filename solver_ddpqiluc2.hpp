
#ifndef __SOLVER_DDPQILUC2__
#define __SOLVER_DDPQILUC2__
#include <iomanip>

#include "inmost_solver.h"
#include "solver_prototypes.hpp"
//#define REPORT_ILU
//#undef REPORT_ILU
//#define REPORT_ILU_PROGRESS
//#undef REPORT_ILU_PROGRESS
using namespace INMOST;

#ifndef DEFAULT_TAU
#define DEFAULT_TAU 0.01
#endif

#define REORDER_DDPQ
#define REORDER_NNZ
#define DDPQ_TAU 0.8 //ratio is time consumed for factorization by time for schur complement calculation
//#define DDPQ_TAU std::max(0.1,1.0/(1.0+exp(-ratio*5)))
#define DDPQ_ADAPT_TAU
#define ESTIMATOR
#define RESCALE_B
#define PIVOT_THRESHOLD
#define DIAGONAL_PERTURBATION
#define DIAGONAL_PERTURBATION_REL 0.001
#define DIAGONAL_PERTURBATION_ABS 1.0e-15
//#define DIAGONAL_PIVOT //probably there is some bug
#define DIAGONAL_PIVOT_TAU 0.01
#define ILUC2
#define ILUC2_TAU pow(tau,1.75)



#define ADDR(row,col) ((row)*size+(col))

class ILUC_preconditioner : public Method
{
	typedef struct Interval_t
	{
		INMOST_DATA_ENUM_TYPE first, last;
		INMOST_DATA_ENUM_TYPE Size() { return last - first; }
	} Interval;
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
	std::vector<Solver::Row::entry> LU_Entries, B_Entries;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LU_Diag;
	interval<INMOST_DATA_ENUM_TYPE, Interval> U_Address, L_Address, B_Address;
	std::vector<Solver::Row::entry> E_Entries, F_Entries;
	std::vector<interval<INMOST_DATA_ENUM_TYPE, Interval> *> E_Address;
	interval<INMOST_DATA_ENUM_TYPE, Interval> F_Address;
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> temp; // temporal place for solve phase
	levels_t level_size; //remember size of each level
	std::vector<Interval> level_interval;
	//reordering information
	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > ddP,ddQ;
	INMOST_DATA_ENUM_TYPE reorder_nnz, ddpq_tau_adapt, estimator;
	INMOST_DATA_REAL_TYPE ddpq_tau, iluc2_tau;
	INMOST_DATA_REAL_TYPE tau, eps;
	INMOST_DATA_ENUM_TYPE sciters;
	Solver::Matrix * Alink;
	Solver::OrderInfo * info;
	bool init;
	void SwapEntries(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 std::vector<Solver::Row::entry> & Entries, 
					 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend, 
					 INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE i;
		INMOST_DATA_REAL_TYPE u;
		for (i = rbeg; i < rend; i++)
		{
			INMOST_DATA_ENUM_TYPE pk = Address[i].last, pj = Address[i].last, r;
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				if (Entries[r].first >= k) { pk = r; break; }
			}
			for (r = pk; r < Address[i].last; r++)
			{
				if (Entries[r].first >= j) { pj = r; break; }
			}
			if (pk < Address[i].last)
			{
				if (pj < Address[i].last && Entries[pj].first == j)
				{
					if (Entries[pk].first == k)
					{
						u = Entries[pj].second;
						Entries[pj].second = Entries[pk].second;
						Entries[pk].second = u;
					}
					else
					{
						u = Entries[pj].second;
						for (r = pj; r > pk; r--) Entries[r] = Entries[r - 1];
						Entries[pk].first = k;
						Entries[pk].second = u;
					}
				}
				else
				{
					if (Entries[pk].first == k)
					{
						u = Entries[pk].second;
						for (r = pk; r < pj - 1; r++) Entries[r] = Entries[r + 1];
						Entries[pj - 1].first = j;
						Entries[pj - 1].second = u;
					}
				}
			}
		}
	}
	void SwapLine(interval<INMOST_DATA_ENUM_TYPE, Interval> & Line, INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE temp;
		temp = Line[i].first;
		Line[i].first = Line[j].first;
		Line[j].first = temp;
		temp = Line[i].last;
		Line[i].last = Line[j].last;
		Line[j].last = temp;
	}
	void SwapE(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE k;
		for (k = 0; k < E_Address.size(); ++k)
			SwapLine(*E_Address[k], i, j);
	}
	
	void ReorderEF(INMOST_DATA_ENUM_TYPE mobeg, 
					INMOST_DATA_ENUM_TYPE cbeg,
					INMOST_DATA_ENUM_TYPE cend,
					INMOST_DATA_ENUM_TYPE wend, 
					interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ, 
					interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
					interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
	{
		(void)cend;
		INMOST_DATA_ENUM_TYPE i, k;
		//Reorder E rows block by P
		//for(interval<INMOST_DATA_ENUM_TYPE, bool>::iterator it = donePQ.begin() + cbeg - mobeg;
		//	it != donePQ.begin() + wend - mobeg; ++it) *it = false;
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
						SwapE(i, localP[k]);
						donePQ[localP[k]] = true;
						k = localP[k];
					} while (k != i);
				}
				donePQ[i] = true;
			}
		}
		//Reorder F columns by  Q
		/*
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
						SwapF(i, localQ[k]);
						donePQ[localQ[k]] = true;
						k = localQ[k];
					} while (k != i);
				}
				donePQ[i] = true;
			}
		}
		*/
		for (i = mobeg; i < cbeg; i++)
		{
			for (k = F_Address[i].first; k < F_Address[i].last; ++k)
			{
				if (F_Entries[k].first >= cbeg)
					F_Entries[k].first = localQ[F_Entries[k].first];
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
		else if( name == "reorder_nnz") return reorder_nnz;
		else if( name == "ddpq_tau_adapt" ) return ddpq_tau_adapt;
		else if( name == "estimator" ) return estimator;
		throw - 1;
	}
	INMOST_DATA_REAL_TYPE & RealParameter(std::string name)
	{
		if (name == "tau") return tau;
		else if( name == "ddpq_tau" ) return ddpq_tau;
		else if( name == "tau2" ) return iluc2_tau;
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
	ILUC_preconditioner(const ILUC_preconditioner & other) :Method(other)
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
#if defined(DDPQ_ADAPT_TAU)
		ddpq_tau_adapt = 1;
#else
		ddpq_tau_adapt = 0;
#endif
#if defined(ESTIMATOR)
		estimator = 1;
#else
		estimator = 0;
#endif
		ddpq_tau = DDPQ_TAU;
		reorder_nnz = 1;
		iluc2_tau = ILUC2_TAU;
	}
	bool isInitialized() { return init; }
	bool isFinalized() { return !init; }
	bool Initialize()
	{
		const INMOST_DATA_REAL_TYPE subst = 1.0;
		const INMOST_DATA_REAL_TYPE tol_modif = 1e-12;
		const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;

		
		INMOST_DATA_ENUM_TYPE wbeg, wend; //working interval
		INMOST_DATA_ENUM_TYPE mobeg, moend; // total interval
		
		INMOST_DATA_ENUM_TYPE k, i, j, Li, Ui, curr, next, mi, mj;
		INMOST_DATA_REAL_TYPE l,u,udiag, max_diag, min_diag, mean_diag;
		INMOST_DATA_ENUM_TYPE nzA, nzLU = 0, nzEF = 0, nzL;
		Solver::Vector DL, DR;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		
		
		//prepare temporal array
		temp.set_interval_beg(mobeg);
		temp.set_interval_end(moend);
		//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE>::iterator rit = temp.begin();
		//	rit != temp.end(); ++rit) *rit = 0.0;
		std::fill(temp.begin(), temp.end(), 0.0);

		
		//prepare reordering vectors
		ddP.set_interval_beg(mobeg);
		ddP.set_interval_end(moend);
		ddQ.set_interval_beg(mobeg);
		ddQ.set_interval_end(moend);

		

		//prepare rescaling vectors
		DL.SetInterval(mobeg, moend);
		DR.SetInterval(mobeg, moend);
		for(Solver::Vector::iterator ri = DL.Begin(); ri != DL.End(); ++ri) *ri = 0.0;

		

		for (k = mobeg; k != moend; k++) ddP[k] = ddQ[k] = k;


		
		//preinitialize data
		level_size.clear();

		

		// supplementary data for ddPQ reordering
		interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(mobeg, moend);

		

		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(mobeg, moend), invQ(mobeg,moend), localP(mobeg, moend), localQ(mobeg,moend);


		
#if defined(REORDER_DDPQ)
		wgt_coords sort_wgts(2 * (moend - mobeg));
#endif
		

		//supplimentary data structures for ILUC
		INMOST_DATA_ENUM_TYPE LU_Beg;
		U_Address.set_interval_beg(mobeg);
		U_Address.set_interval_end(moend);
		L_Address.set_interval_beg(mobeg);
		L_Address.set_interval_end(moend);
		F_Address.set_interval_beg(mobeg);
		F_Address.set_interval_end(moend);
		B_Address.set_interval_beg(mobeg);
		B_Address.set_interval_end(moend);
		LU_Diag.set_interval_beg(mobeg);
		LU_Diag.set_interval_end(moend);
		std::vector<Solver::Row::entry> A_Entries, C_Entries, S_Entries, LF_Entries, LFt_Entries;
		//std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend), C_Address(mobeg,moend), S_Address(mobeg,moend), LF_Address(mobeg,moend), LFt_Address(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ulist(mobeg, moend), Llist(mobeg, moend), Blist(mobeg,moend), Flist(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(mobeg, moend,EOL), Lbeg(mobeg, moend,EOL), Bbeg(mobeg,moend,EOL), Fbeg(mobeg,moend,EOL);

		

		//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
		INMOST_DATA_REAL_TYPE mup, mum, NuU, NuL, NuD, NuU_old, NuL_old, vp, vm, v;
		INMOST_DATA_ENUM_TYPE np, nm;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL(mobeg, moend,0.0), EstU(mobeg, moend,0.0);
		//supplimentary data structures for returning values of dropped elements
		//INMOST_DATA_REAL_TYPE DropLk, DropUk;
		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
		//data structure for linked list

		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValues(mobeg, moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndeces(mobeg, moend+1,UNDEF);
		double tfactor = 0.0, trescale = 0.0, tschur = 0.0, treorder = 0.0, treassamble = 0.0, ttotal, tt;
		double tlfactor, tlrescale, tlschur, tlreorder, tlreassamble, ratio = 1.0, ratio2 = 1.0;
		ttotal = Timer();

		
		//calculate number of nonzeros
		nzA = 0;
		for (k = mobeg; k < moend; ++k)
		{
			for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0) nzA++;
		}

		//sort_indeces.reserve(256);
		A_Entries.resize(nzA);
		B_Entries.reserve(nzA);
		C_Entries.reserve(nzA);
		S_Entries.reserve(nzA);
		E_Entries.reserve(nzA*2);
		F_Entries.reserve(nzA*2);
		LU_Entries.reserve(nzA);
		LF_Entries.reserve(nzA);
		LFt_Entries.reserve(nzA);

		j = 0;
		for (k = mobeg; k < moend; ++k)
		{
			A_Address[k].first = j;
			for (Solver::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
			{
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
					A_Entries[j++] = Solver::Row::make_entry(r->first, r->second);
			}
			A_Address[k].last = j;
			assert(A_Address[k].Size() != 0); //singular matrix
		}

#if defined(ILUC2)
		INMOST_DATA_ENUM_TYPE nzLU2 = 0;
		INMOST_DATA_REAL_TYPE tau2 = iluc2_tau;
		std::vector<Solver::Row::entry> LU2_Entries;
		interval<INMOST_DATA_ENUM_TYPE, Interval> L2_Address(mobeg, moend+1), U2_Address(mobeg, moend+1);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2list(mobeg, moend, UNDEF), L2list(mobeg, moend, UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2beg(mobeg, moend, EOL), L2beg(mobeg, moend, EOL);
		LU2_Entries.reserve(nzA * 4);
#else
		INMOST_DATA_REAL_TYPE tau2 = tau;
#endif

#if defined(REPORT_ILU)
		std::cout << "nonzeros in A " << nzA << std::endl;
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
			//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = localP.begin()+(wbeg-mobeg);
			//	it != localP.begin() + (wend-mobeg); ++it) *it = ENUMUNDEF;
			//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = localQ.begin()+(wbeg-mobeg);
			//	it != localQ.begin() + (wend-mobeg); ++it) *it = ENUMUNDEF;
			std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
			std::fill(localQ.begin() + (wbeg - mobeg), localQ.begin() + (wend - mobeg), ENUMUNDEF);

#if defined(REORDER_DDPQ)

#if defined(REPORT_ILU)
			std::cout << level_size.size() << " calculate weights" << std::endl;
#endif
			//for(Solver::Vector::iterator rit = DL.Begin()+wbeg-mobeg;
			//	rit != DL.Begin()+wend-mobeg;++rit) *rit = 0.0;
			std::fill(DL.Begin() + wbeg - mobeg, DL.Begin() + wend - mobeg, 0.0);
			//for(Solver::Vector::iterator rit = DR.Begin()+wbeg-mobeg;
			//	rit != DR.Begin()+wend-mobeg;++rit) *rit = 0.0;
			std::fill(DR.Begin() + wbeg - mobeg, DR.Begin() + wend - mobeg, 0.0);
			//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = Ulist.begin()+wbeg-mobeg;
			//	it != Ulist.begin()+wend-mobeg; ++it) *it = 0;
			std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, 0);
			//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = Llist.begin()+wbeg-mobeg;
			//	it != Llist.begin()+wend-mobeg; ++it) *it = 0;
			std::fill(Llist.begin() + wbeg - mobeg, Llist.begin() + wend - mobeg, 0);
			for (k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
				{
					i = A_Entries[it].first;
					u = fabs(A_Entries[it].second);
					//row-col norms
					DL[k] += u; //row norm
					DR[i] += u; //col norm
					//nonzeros
					Ulist[k]++; //row nnz
					Llist[i]++; //col nnz
				}
			}

			INMOST_DATA_REAL_TYPE ddmaxall = 0, ddmaxcur, ddcur, ddmaxmean = 0.0, ddmaxmin = 1.0e+20;
			for (k = wbeg; k < wend; ++k)
			{
				ddmaxcur = 0.0;
				Blist[k] = ENUMUNDEF;
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
				{
					i = A_Entries[it].first;
					u = fabs(A_Entries[it].second);
					ddcur = u / (DL[k] + DR[i] - u);
					if (ddmaxcur < ddcur)
					{
						ddmaxcur = ddcur;
						Blist[k] = i;
					}
				}
				temp[k] = ddmaxcur;
				if (ddmaxall < ddmaxcur) ddmaxall = ddmaxcur;
				if (ddmaxmin > ddmaxcur) ddmaxmin = ddmaxcur;
				ddmaxmean += ddmaxcur*ddmaxcur;
			}

			ddmaxmean = sqrt(ddmaxmean / (1.0*(wend - wbeg)));
			
			//first select those elements that are largest over row and column simultaneously
			sort_wgts.clear();
			for (k = wbeg; k < wend; ++k)
			{
				if (Blist[k] != ENUMUNDEF)
				{
					if (temp[k] > ddmaxall * ddpq_tau)
						sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
				}
			}
			if(reorder_nnz) std::sort(sort_wgts.begin(), sort_wgts.end());
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " fill reordering" << std::endl;
#endif
			i = wbeg;
			for (wgt_coords::iterator it = sort_wgts.begin(); it != sort_wgts.end(); ++it)
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
			if (cbeg == cend && cbeg != wend)
			{
				//std::cout << __FILE__ << ":" << __LINE__ << " singular matrix, factored " << mobeg << ".." << cend << " out of " << mobeg << ".." << moend << std::endl;
				for (k = cbeg; k < cend; ++k)
					LU_Diag[k] = tol_modif;
				break;
			}
#if defined(REORDER_DDPQ)
#if defined(REPORT_ILU)
			//std::cout << level_size.size() << " new level " << cbeg << ".." << cend << " out of " << wbeg << ".." << wend << std::endl;
			std::cout << std::setw(8) << level_size.size() << " total " << std::setw(8) << wend - wbeg << " after tau filtering " << std::setw(8) << sort_wgts.size() << " selected " << std::setw(8) << cend - wbeg;
			std::cout << std::scientific << " max " << std::setw(8) << ddmaxall << " mean " << std::setw(8) << ddmaxmean << " min " << std::setw(8) << ddmaxmin << " ratio " << ddmaxall/ddmaxmean;
			std::cout << std::endl;
#endif
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
			tlreorder = Timer() - tt;
			treorder += tlreorder;
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
			//reorder blocks that we have
			ReorderEF(mobeg, cbeg, cend, wend, donePQ, localP, localQ);
			//allocate space for E,F blocks for current level
			E_Address.push_back(new interval<INMOST_DATA_ENUM_TYPE, Interval>(cend, wend));
			//inverse ordering
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			tt1 = Timer() - tt1;

			//std::cout << "reorder: " << tt1 << std::endl;
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
				LU_Diag[k] = 0;
				F_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(F_Entries.size());
				B_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt != A_Address[invP[k]].last; ++jt)
				{
					i = localQ[A_Entries[jt].first];
					u = A_Entries[jt].second;
					if (i < cend) B_Entries.push_back(Solver::Row::make_entry(i, u));
					else F_Entries.push_back(Solver::Row::make_entry(i,u)); //count nonzero in column of F
				}
				F_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(F_Entries.size());
				B_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				std::sort(B_Entries.begin() + B_Address[k].first, B_Entries.end());
				std::sort(F_Entries.begin() + F_Address[k].first, F_Entries.end());
				nzEF += F_Address[k].Size();
			}
			//Setup column addressing for B,F, in descending order to keep list ordered
			for (k = cend; k > cbeg; --k)
			{
				if (B_Address[k - 1].Size() > 0)
				{
					i = B_Entries[B_Address[k - 1].first].first;
					if (i < k-1)
					{
						Blist[k-1] = Bbeg[i];
						Bbeg[i] = k-1;
					}
				}
			}

			tt2 = Timer() - tt2;
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " rebuilding matrix, E,C blocks " << std::endl;
#endif

			tt3 = Timer();
			//compute E,C
			C_Entries.clear();
			for (k = cend; k < wend; ++k)
			{
				E_Address.back()->at(k).first = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries.size());
				C_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries.size());
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt != A_Address[invP[k]].last; ++jt)
				{
					i = localQ[A_Entries[jt].first];
					u = A_Entries[jt].second;
					if (i < cend) E_Entries.push_back(Solver::Row::make_entry(i, u)); //put to row of E
					else C_Entries.push_back(Solver::Row::make_entry(i, u)); //form new Schur complement
				}
				E_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries.size());
				C_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries.size());
				// sort entries added to E, because it will be very usefull later
				std::sort(E_Entries.begin() + E_Address.back()->at(k).first, E_Entries.end());
				nzEF += E_Address.back()->at(k).Size();
			}
			tt3 = Timer() - tt3;
			ttt = Timer();
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);
			tt1 += Timer() - ttt;

#if defined(REPORT_ILU)
			std::cout << level_size.size() << " finished permutation, E,F " << tt1 << " B,F " << tt2 << " E,C " << tt3 << std::endl;
#endif
			tlreassamble = Timer() - tt;
			treassamble += tlreassamble;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING                 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Rescale current B block by Row-Column alternating scaling
			tt = Timer();
#if defined(RESCALE_B)
#if defined(REPORT_ILU)
			std::cout << level_size.size() << " rescaling block B " << std::endl;
#endif
			//for(Solver::Vector::iterator it = DL.Begin() + cbeg - mobeg;
			//	it != DL.Begin() + cend - mobeg; ++it) *it = 0.0;
			std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
					DL[k] += B_Entries[r].second*B_Entries[r].second;
			}
			for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				//for(Solver::Vector::iterator it = DR.Begin()+cbeg-mobeg; 
				//	it != DR.Begin()+cend-mobeg; ++it) *it = 0.0;
				std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DR[B_Entries[r].first] += DL[k] * B_Entries[r].second*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				//for(Solver::Vector::iterator it = DL.Begin()+cbeg-mobeg;
				//	it != DL.Begin()+cend-mobeg; ++it) *it = 0.0;
				std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DL[k] += B_Entries[r].second*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			}
			for (k = cbeg; k < cend; k++) DL[k] = sqrt(DL[k]);
			for (k = cbeg; k < cend; k++) DR[k] = sqrt(DR[k]);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];
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
				mean_diag = 0.0;
				//initialize diagonal values
				for (k = cbeg; k < cend; k++)
				{
					localP[k] = F_Address[k].first; //remember F block starts permutation
					localQ[k] = B_Address[k].first; //remember B block starts permutation for pivoting
					for (i = B_Address[k].first; i < B_Address[k].last; i++)
					{
						if (B_Entries[i].first == k)
						{
							LU_Diag[k] = B_Entries[i].second;
							mean_diag += fabs(LU_Diag[k]);
							break;
						}
					}
				}
				//initialize first row and column
				mean_diag -= fabs(LU_Diag[cbeg]);
				assert(B_Entries[B_Address[cbeg].first].first == cbeg);
				if (B_Entries[B_Address[cbeg].first].first == cbeg)
					LU_Diag[cbeg] = B_Entries[B_Address[cbeg].first].second;
				else
				{
					std::cout << "No diagonal value! " << cbeg << " " << B_Entries[B_Address[cbeg].first].first << std::endl;
					LU_Diag[cbeg] = 0.0;
				}
#if defined(PIVOT_THRESHOLD)
				if (fabs(LU_Diag[cbeg]) < tol_modif)
				{
					LU_Diag[cbeg] = LU_Diag[cbeg] < 0.0 ? -tol_modif : tol_modif;
				}
				else
#endif
#if defined(DIAGONAL_PERTURBATION)
				{
					LU_Diag[cbeg] = LU_Diag[cbeg] * (1.0 + DIAGONAL_PERTURBATION_REL) + (LU_Diag[cbeg] < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS;
				}
#endif
				LU_Beg = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
				U_Address[cbeg].first = LU_Beg;
				for (INMOST_DATA_ENUM_TYPE r = B_Address[cbeg].first + (B_Entries[B_Address[cbeg].first].first == cbeg ? 1 : 0); r != B_Address[cbeg].last; ++r) 
					LU_Entries.push_back(Solver::Row::make_entry(B_Entries[r].first, B_Entries[r].second / LU_Diag[cbeg]));
				L_Address[cbeg].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
				U_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
				Li = Bbeg[cbeg];
				while (Li != EOL)
				{
					LU_Entries.push_back(Solver::Row::make_entry(Li, B_Entries[B_Address[Li].first].second / LU_Diag[cbeg]));
					Li = Blist[Li];
				}
				L_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
				//update diagonal by factors
				i = U_Address[cbeg].first;
				j = L_Address[cbeg].first;
				while (i != U_Address[cbeg].last && j != L_Address[cbeg].last)
				{
					Ui = LU_Entries[i].first;
					Li = LU_Entries[j].first;
					if (Ui > Li) j++;
					else if (Ui < Li) i++;
					else
					{
						LU_Diag[Ui] -= LU_Entries[i].second * LU_Entries[j].second * LU_Diag[cbeg];
						i++;
						j++;
					}
				}
				//update column indexes
				Li = Bbeg[cbeg];
				while (Li != EOL)
				{
					B_Address[Li].first++;
					Ui = Blist[Li];
					if (B_Address[Li].Size() > 0)
					{
						i = B_Entries[B_Address[Li].first].first;
						if (i < Li)
						{
							if (Bbeg[i] > Li)
							{
								Blist[Li] = Bbeg[i];
								Bbeg[i] = Li;
							}
							else
							{
								curr = next = Bbeg[i];
								while (next < Li)
								{
									curr = next;
									next = Blist[curr];
								}
								assert(curr < Li);
								assert(Li < next);
								Blist[curr] = Li;
								Blist[Li] = next;
							}
						}
					}
					Li = Ui;
				}
				//form partial list from this row
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
				if( estimator )
				{
					EstU[cbeg] = EstL[cbeg] = 0.0;
					for (INMOST_DATA_ENUM_TYPE it = U_Address[cbeg].first; it != U_Address[cbeg].last; ++it) EstU[LU_Entries[it].first] = LU_Entries[it].second;
					for (INMOST_DATA_ENUM_TYPE it = L_Address[cbeg].first; it != L_Address[cbeg].last; ++it) EstL[LU_Entries[it].first] = LU_Entries[it].second;
				}				
				nzLU += L_Address[cbeg].Size() + U_Address[cbeg].Size() + 1;
				max_diag = min_diag = fabs(LU_Diag[cbeg]);
				NuD = 1.0;
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " starting factorization " << std::endl;
#endif
				for (k = cbeg + 1; k < cend; ++k)
				{
					//observe that there is a better diagonal pivot
#if defined(DIAGONAL_PIVOT)
					/*
					double mean_diag2 = 0.0;
					for (i = k; i < cend; i++)
					{
						mean_diag2 += fabs(LU_Diag[i]);
					}
					std::cout << mean_diag << " " << mean_diag2 << std::endl;
					*/

					if (fabs(LU_Diag[k])*(cend-k) < DIAGONAL_PIVOT_TAU*mean_diag)
					{
						j = k;
						for (i = k + 1; i < cend; i++) if (fabs(LU_Diag[j]) < fabs(LU_Diag[i])) j = i;

						//if (fabs(LU_Diag[k]) / fabs(LU_Diag[j]) < tau)
						{
							//TODO!!!
							std::cout << "Detected, that there is a much better pivot, i'm " << k << " " << LU_Diag[k] << " other " << j << " " << LU_Diag[j] << std::endl;
							scanf("%*c");
							//std::cout << "But repivoting algorithm not implemented" << std::endl;
							//This algorithm may be quite costly, but the effect is miraculous
							//First correct E
							SwapEntries(*E_Address.back(), E_Entries, cend, wend, k, j);
							//Then correct F
							SwapLine(F_Address,k, j);
							//Swap lines of B
							SwapLine(B_Address, k, j);
							//Correct B the same way E is corrected
							SwapEntries(B_Address, B_Entries, cbeg, cend, k, j);
							//now reintialize entrie column linked list
							for (i = k; i < cend; ++i) Bbeg[i] = EOL;
							for (i = cend; i > k; --i)
							{
								Li = B_Entries[B_Address[i - 1].first].first;
								if (Li < i - 1)
								{
									Blist[i - 1] = Bbeg[Li];
									Bbeg[Li] = i - 1;
								}
							}
							//now update L
							SwapEntries(L_Address, LU_Entries, cbeg, cend, k, j);
							//reinitialize linked list
							//now update U
							SwapEntries(U_Address, LU_Entries, cbeg, cend, k, j);
							//reinitialize linked list
#if defined(ILUC2)
							//now update L2
							SwapEntries(U2_Address, LU2_Entries, cbeg, cend, k, j);
							//reinitialize linked list
							//now update U2
							SwapEntries(L2_Address, LU2_Entries, cbeg, cend, k, j);
							//reinitialize linked list
#endif
							for (i = k; i < cend; i++)
							{
								Ubeg[i] = EOL;
								Lbeg[i] = EOL;
#if defined(ILUC2)
								U2beg[i] = EOL;
								L2beg[i] = EOL;
#endif
							}
							for (i = cbeg; i < k; i++)
							{
								if (U_Address[i].Size() > 0)
								{
									Li = LU_Entries[U_Address[i].first].first;
									Ulist[i] = Ubeg[Li];
									Ubeg[Li] = i;
								}
								if (L_Address[i].Size() > 0)
								{
									Li = LU_Entries[L_Address[i].first].first;
									Llist[i] = Lbeg[Li];
									Lbeg[Li] = i;
								}
#if defined(ILUC2)
								if (U2_Address[i].Size() > 0)
								{
									Li = LU2_Entries[U2_Address[i].first].first;
									U2list[i] = U2beg[Li];
									U2beg[Li] = i;
								}
								if (L2_Address[i].Size() > 0)
								{
									Li = LU2_Entries[L2_Address[i].first].first;
									L2list[i] = L2beg[Li];
									L2beg[Li] = i;
								}
#endif
							}
							//swap diagonal values
							u = LU_Diag[k];
							LU_Diag[k] = LU_Diag[j];
							LU_Diag[j] = u;
							//swap values in estimator
							u = EstL[k];
							EstL[k] = EstL[j];
							EstL[j] = u;
							u = EstU[k];
							EstU[k] = EstU[j];
							EstU[j] = u;

							//swap rescaling vectors
							u = DL[k];
							DL[k] = DL[j];
							DL[j] = u;

							u = DR[k];
							DR[k] = DR[j];
							DR[j] = u;
							//update indeces in ddPQ
							Ui = ddP[k];
							ddP[k] = ddP[j];
							ddP[j] = Ui;
							Ui = ddQ[k];
							ddQ[k] = ddQ[j];
							ddQ[j] = Ui;

							//update interval starts for B block to correctly restore intervals
							Ui = localQ[k];
							localQ[k] = localQ[j];
							localQ[j] = Ui;

							Ui = localP[k];
							localP[k] = localP[j];
							localP[j] = Ui;
						}
						/*
						else
						{
							std::cout << "Detected that there are no good pivots anymore! interval " << cbeg << "..." << cend << " current element " << k << std::endl;
							std::cout << "But level rejection algorithm is not implemented" << std::endl;
							//probably should reject all unfactored parts for the next level
						}
						*/
					}
					mean_diag -= fabs(LU_Diag[k]); //remove current value
#endif
					
					//DropLk = DropUk = 0.0;
					//uncompress k-th row
					// add diagonal value first, there shouldn't be values on left from diagonal
					assert(B_Entries[B_Address[k].first].first == k);
					LineIndeces[cbeg] = k;
					if (B_Entries[B_Address[k].first].first == k)
						LineValues[k] = B_Entries[B_Address[k].first].second;
					else
					{
						std::cout << "No diagonal value! " << k << " " << B_Entries[B_Address[k].first].first << std::endl;
						LineValues[k] = 0.0;
					}
#if defined(PIVOT_THRESHOLD)
					if (fabs(LineValues[k]) < tol_modif)
					{
						LineValues[k] = LineValues[k] < 0.0 ? -tol_modif : tol_modif;
					}
#endif
#if defined(DIAGONAL_PERTURBATION)
					LineValues[k] = LineValues[k] * (1.0 + DIAGONAL_PERTURBATION_REL) + (LineValues[k] < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS;
#endif
					Ui = k;
					for (INMOST_DATA_ENUM_TYPE it = B_Address[k].first + (B_Entries[B_Address[k].first].first == k ? 1 : 0); it != B_Address[k].last; ++it)
					{
						LineValues[B_Entries[it].first] = B_Entries[it].second;
						Ui = LineIndeces[Ui] = B_Entries[it].first;
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
						
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								u = l*LU2_Entries[it].second;
								//eliminate values
								if (LineIndeces[j] != UNDEF) //there is an entry in the list
									LineValues[j] -= u;
								else if (fabs(u)*NuU > tau)//add new entry
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
									LineValues[j] = -u;
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
							j = LU_Entries[it].first;
							u = l*LU_Entries[it].second;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= u;
							else if (fabs(u)*NuU > tau)//add new entry
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
								LineValues[j] = -u;
							}
						}

						for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
						{
							u = LU2_Entries[it].second;
							j = LU2_Entries[it].first;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= l*u;
						}

						i = L2list[i];
					}
					U2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					//std::cout << std::endl;
					//define the diagonal value
#if defined(DIAGONAL_PIVOT)
					if (fabs(LU_Diag[k] - LineValues[k]) > 1e-9) std::cout << __LINE__ << " Diagonal value went wrong, good: " << LineValues[k] << " have " << LU_Diag[k] << " line " << k << std::endl;
#endif
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
					while (Ui != EOL)
					{
						LineValues[Ui] /= udiag;
						Ui = LineIndeces[Ui];
					}
					//estimate condition for U^{-1}
					if( estimator )
					{
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
						NuU_old = NuU;
						if (np > nm) NuU = mup; else NuU = mum;
						Ui = LineIndeces[k];
						while (Ui != EOL)
						{
							EstU[Ui] += LineValues[Ui] * NuU;
							Ui = LineIndeces[Ui];
						}
						NuU = std::max(fabs(mup), fabs(mum));
						NuU = std::max(NuU_old,NuU);
					}
					//insert line to U part
					U_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					Ui = LineIndeces[k];
					while (Ui != EOL)
					{
						u = fabs(LineValues[Ui]);
						if (u > tau/NuU) // apply dropping rule
							LU_Entries.push_back(Solver::Row::make_entry(Ui, LineValues[Ui]));
#if defined(ILUC2)
						else if (u > tau2)
							LU2_Entries.push_back(Solver::Row::make_entry(Ui, LineValues[Ui]));
#endif
						Ui = LineIndeces[Ui];
					}
					U_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());

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
					U2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
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
					if (B_Entries[B_Address[k].first].first == k)
						LineValues[k] = B_Entries[B_Address[k].first].second;
					else
					{
						std::cout << "No diagonal value! " << k << std::endl;
						LineValues[k] = 0.0;
					}
					//start from diagonal
					//sort_indeces.clear();
					Ui = k;
					Li = Bbeg[k];
					while (Li != EOL)
					{
						assert(B_Entries[B_Address[Li].first].first == k);
						//sort_indeces.push_back(Li);
						LineValues[Li] = B_Entries[B_Address[Li].first].second;
						Ui = LineIndeces[Ui] = Li;
						Li = Blist[Li];
					}
					/*
					std::sort(sort_indeces.begin(), sort_indeces.end());
					for (Li = 0; Li < sort_indeces.size(); Li++)
					Ui = LineIndeces[Ui] = sort_indeces[Li];
					*/
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
						
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								l = u*LU2_Entries[it].second;
								if (LineIndeces[j] != UNDEF) //there is an entry in the list
									LineValues[j] -= l;
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
									LineValues[j] = -l;
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
							l = u*LU_Entries[it].second;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= l;
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
								LineValues[j] = -l;
							}

						}
						for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it != L2_Address[i].last; ++it)
						{
							j = LU2_Entries[it].first;
							l = LU2_Entries[it].second;
							//eliminate values
							if (LineIndeces[j] != UNDEF) //there is an entry in the list
								LineValues[j] -= l*u;
						}
						i = U2list[i];
					}
					L2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					//check that diagonal value is the same(must be!)
					//assert(fabs(LineValues[k] / udiag - 1.0) < 1.0e-10);
#if defined(DIAGONAL_PIVOT)
					if (fabs(LineValues[k] - udiag) > 1e-9) std::cout << __LINE__ << " Diagonal value went wrong, good: " << udiag << " have " << LineValues[k] << " line " << k << std::endl;
#endif
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
					if( estimator )
					{
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
						NuL_old = NuL;
						if (np > nm) NuL = mup; else NuL = mum;
						Li = LineIndeces[k];
						while (Li != EOL)
						{
							EstL[Li] += LineValues[Li] * NuL;
							Li = LineIndeces[Li];
						}
						NuL = std::max(fabs(mup), fabs(mum));
						NuL = std::max(NuL,NuL_old);
					}
					//insert column to L part
					Li = LineIndeces[k];
					L_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					while (Li != EOL)
					{
						u = fabs(LineValues[Li]);
						if (u > tau/NuL) //apply dropping 
							LU_Entries.push_back(Solver::Row::make_entry(Li, LineValues[Li]));
#if defined(ILUC2)
						else if (u > tau2)
							LU2_Entries.push_back(Solver::Row::make_entry(Li, LineValues[Li]));
#endif
						Li = LineIndeces[Li];
					}
					L_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
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
					L2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
					if (L2_Address[k].Size() != 0)
					{
						i = LU2_Entries[L2_Address[k].first].first;
						assert(i > k);
						L2list[k] = L2beg[i];
						L2beg[i] = k;
					}
#endif
					//assert(fabs(LU_Diag[k] - udiag) < 1e-5);

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
					nzLU2 += U2_Address[k].Size() + L2_Address[k].Size();
#endif
					//update column
					Li = Bbeg[k];
					while (Li != EOL)
					{
						B_Address[Li].first++;
						Ui = Blist[Li];
						if (B_Address[Li].Size() > 0)
						{
							i = B_Entries[B_Address[Li].first].first;
							if (i < Li)
							{
								if (Bbeg[i] > Li)
								{
									Blist[Li] = Bbeg[i];
									Bbeg[i] = Li;
								}
								else
								{
									curr = next = Bbeg[i];
									while (next < Li)
									{
										curr = next;
										next = Blist[curr];
									}
									assert(curr < Li);
									assert(Li < next);
									Blist[curr] = Li;
									Blist[Li] = next;
								}
							}
						}
						Li = Ui;
					}
#if defined(DIAGONAL_PIVOT)
					//update diagonal by optained factors
					i = U_Address[k].first;
					j = L_Address[k].first;
					while (i != U_Address[k].last && j != L_Address[k].last)
					{
						Ui = LU_Entries[i].first;
						Li = LU_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							mean_diag -= fabs(LU_Diag[Ui]);
							LU_Diag[Ui] -= LU_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
							mean_diag += fabs(LU_Diag[Ui]);
							i++;
							j++;
						}
					}
#if defined(ILUC2)
					i = U_Address[k].first;
					j = L2_Address[k].first;
					while (i != U_Address[k].last && j != L2_Address[k].last)
					{
						Ui = LU_Entries[i].first;
						Li = LU2_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							mean_diag -= fabs(LU_Diag[Ui]);
							LU_Diag[Ui] -= LU_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
							mean_diag += fabs(LU_Diag[Ui]);
							i++;
							j++;
						}
					}
					i = U2_Address[k].first;
					j = L_Address[k].first;
					while (i != U2_Address[k].last && j != L_Address[k].last)
					{
						Ui = LU2_Entries[i].first;
						Li = LU_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							mean_diag -= fabs(LU_Diag[Ui]);
							LU_Diag[Ui] -= LU2_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
							mean_diag += fabs(LU_Diag[Ui]);
							i++;
							j++;
						}
					}
					i = U2_Address[k].first;
					j = L2_Address[k].first;
					while (i != U2_Address[k].last && j != L2_Address[k].last)
					{
						Ui = LU2_Entries[i].first;
						Li = LU2_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							mean_diag -= fabs(LU_Diag[Ui]);
							LU_Diag[Ui] -= LU2_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
							mean_diag += fabs(LU_Diag[Ui]);
							i++;
							j++;
						}
					}
#endif
#endif

					//number of nonzeros
					nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
					
#if defined(REPORT_ILU)
					if (k % 100 == 0)
					{
						
						std::cout << std::fixed << std::setprecision(2) << std::setw(6) << 100.0f*(k - cbeg) / (float)(cend - cbeg) << "%";
						std::cout << " nnz U " << std::setprecision(10) << std::setw(10) << nzLU;
#if defined(ILUC2)
						std::cout << " U2 " << std::setw(10) << nzLU2;
#endif
						std::cout << std::scientific << std::setprecision(5) << " condition L " << NuL << " D " << NuD << " U " << NuU;
						std::cout << "\r";
						std::cout.flush();
						std::cout << std::setprecision(0);
						//std::cout.setf(0,std::ios::floatfield);
						
						//printf("%6.2f%% nnz LU %8d condition L %10f D %10f U %10f\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, NuL, NuD, NuU);
						//fflush(stdout);
					}
#elif defined(REPORT_ILU_PROGRESS)
					if (k % 100 == 0)
					{
						printf("%lu %d/%d factor %6.2f%%\t\t\r",static_cast<unsigned long>(level_size.size()), cend,moend, 100.0f*(k - cbeg) / (float)(cend - cbeg));
						//printf("%6.2f%% nnz LU %8d condition L %10f D %10f U %10f\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, NuL, NuD, NuU);
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
			tlfactor = Timer() - tt;
			tfactor += tlfactor;
			//restore indexes
			B_Address[cbeg].first = localQ[cbeg];
			U_Address[cbeg].first = LU_Beg;
			L_Address[cbeg].first = U_Address[cbeg].last;
			for (k = cbeg+1; k < cend; ++k)
			{
				B_Address[k].first = localQ[k];
				U_Address[k].first = L_Address[k - 1].last;
				L_Address[k].first = U_Address[k].last;
#if defined(ILUC2)
				U2_Address[k].first = L2_Address[k - 1].last;
				L2_Address[k].first = U2_Address[k].last;
#endif
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
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) B_Entries[r].second /= DL[k] * DR[B_Entries[r].first];
#if defined(ILUC2)
				for (INMOST_DATA_ENUM_TYPE r = U2_Address[k].first; r < U2_Address[k].last; ++r) LU2_Entries[r].second *= DR[k] / DR[LU2_Entries[r].first];
				for (INMOST_DATA_ENUM_TYPE r = L2_Address[k].first; r < L2_Address[k].last; ++r) LU2_Entries[r].second *= DL[k] / DL[LU2_Entries[r].first];
#endif
			}
#endif
			tlrescale = Timer() - tt;
			trescale += tlrescale;
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
			//Setup column addressing for F block
			//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = Fbeg.begin() + cend - mobeg;
			//	it != Fbeg.begin() + wend - mobeg; ++it) *it = EOL;
			std::fill(Fbeg.begin() + cend - mobeg, Fbeg.begin() + wend - mobeg, EOL);
			for (k = cend; k > cbeg; --k)
			{
				if (F_Address[k - 1].Size() > 0)
				{
					i = F_Entries[F_Address[k - 1].first].first;
					assert(i >= cend);
					assert(i < wend);
					Flist[k - 1] = Fbeg[i];
					Fbeg[i] = k - 1;
				}
			}

			{
				LF_Entries.clear();
				//Compute entire LF matrix with row-major order
				t1 = Timer();
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// LF block    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (i = cend; i < wend; ++i) 
				{
					if (Fbeg[i] == EOL)//iterate over columns of F
					{
						LF_Address[i].first = LF_Address[i].last = 0;
						continue;
					}
					//unpack row to linked list
					Li = cbeg;
					Ui = Fbeg[i];
					while (Ui != EOL)  //iterate over values of i-th column of F
					{
						assert(F_Entries[F_Address[Ui].first].first == i);
						LineValues[Ui] = F_Entries[F_Address[Ui].first].second;
						Li = LineIndeces[Li] = Ui + 1;
						Ui = Flist[Ui];
					}
					
					LineIndeces[Li] = EOL;
					
					//Update column indexing
					Li = Fbeg[i];
					while (Li != EOL)
					{
						Ui = Flist[Li];
						F_Address[Li].first++;
						if (F_Address[Li].Size() > 0)
						{
							k = F_Entries[F_Address[Li].first].first;
							assert(k >= cend);
							assert(k < wend);
							if (Fbeg[k] > Li)
							{
								Flist[Li] = Fbeg[k];
								Fbeg[k] = Li;
							}
							else
							{
								curr = next = Fbeg[k];
								while (next < Li)
								{
									curr = next;
									next = Flist[curr];
								}
								assert(curr < Li);
								assert(Li < next);
								Flist[curr] = Li;
								Flist[Li] = next;
							}
						}
						Li = Ui;
					}

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
							else if( fabs(u) > std::min(tau/NuU,tau2) )
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
#if defined(ILUC2)
						if (fabs(LineValues[Li - 1]) > tau)
						{
							curr = Li;
							for (INMOST_DATA_ENUM_TYPE ru = L2_Address[Li - 1].first; ru < L2_Address[Li - 1].last; ++ru)
							{
								u = LineValues[Li - 1] * LU2_Entries[ru].second;
								j = LU2_Entries[ru].first;
								if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
									LineValues[j] -= u;
								else if( fabs(u) > std::min(tau/NuU,tau2) )
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
						}
#endif
						Li = LineIndeces[Li];
					}
					
					//Assemble column into matrix
					Li = LineIndeces[cbeg];
					LF_Address[i].first = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
					while (Li != EOL)
					{
						if (fabs(LineValues[Li - 1]) > std::min(tau/NuU,tau2) )
							LF_Entries.push_back(Solver::Row::make_entry(Li - 1, LineValues[Li - 1]));
						Li = LineIndeces[Li];
					}
					LF_Address[i].last = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
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
					if (i % 100 == 0)
					{
						printf("%lu %d/%d schur1 %6.2f%%\t\t\r", static_cast<unsigned long>(level_size.size()), cend,moend, 100.f*(i - cend) / (1.f*(wend - cend)));
						fflush(stdout);
					}
#endif
				}
				//Restore F indexing
				for (i = cbeg; i < cend; ++i)
					F_Address[i].first = localP[i];
				//setup row indeces for LF block
				//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = Fbeg.begin() + cbeg - mobeg;
				//	it != Fbeg.begin() + cend - mobeg; ++it ) *it = EOL;
				std::fill(Fbeg.begin() + cbeg - mobeg, Fbeg.begin() + cend - mobeg, EOL);
				for (k = wend; k > cend; --k)
				{
					if (LF_Address[k-1].Size() > 0)
					{
						i = LF_Entries[LF_Address[k-1].first].first;
						Flist[k-1] = Fbeg[i];
						Fbeg[i] = k-1;
					}
				}
				//transpone LF matrix
				LFt_Entries.resize(LF_Entries.size());
				j = 0;
				for (k = cbeg; k < cend; k++)
				{
					LFt_Address[k].first = j;
					Li = Fbeg[k];
					while (Li != EOL)
					{
						LFt_Entries[j].first = Li;
						LFt_Entries[j].second = LF_Entries[LF_Address[Li].first].second;
						j++;
						Li = Flist[Li];
					}
					LFt_Address[k].last = j;
					//Update row indexing
					Li = Fbeg[k];
					while (Li != EOL)
					{
						Ui = Flist[Li];
						LF_Address[Li].first++;
						if (LF_Address[Li].Size() > 0)
						{
							i = LF_Entries[LF_Address[Li].first].first;
							Flist[Li] = Fbeg[i];
							Fbeg[i] = Li;
						}
						Li = Ui;
					}
				}
				t1 = Timer() - t1;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// EU and Schur blocks  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				t2 = Timer();
				//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = Ulist.begin() + cend - mobeg;
				//	it != Ulist.begin() + wend - mobeg; ++it) *it = UNDEF;
				std::fill(Ulist.begin() + cend - mobeg, Ulist.begin() + wend - mobeg, UNDEF);
				S_Entries.clear();
				
				for (i = cend; i < wend; ++i) 
				{
					if (E_Address.back()->at(i).Size() == 0) //iterate over rows of E
					{
						S_Address[i].first = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
						for (k = C_Address[i].first; k < C_Address[i].last; k++)
							S_Entries.push_back(C_Entries[k]);
						S_Address[i].last = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
						continue;
					}
					double tt0, tt1, tt2;
					tt0 = Timer();
					// uncompress sprase i-th row of E to linked list
					Li = cbeg;
					for (j = E_Address.back()->at(i).first; j < E_Address.back()->at(i).last; ++j) //iterate over values of i-th row of E
					{
						LineValues[E_Entries[j].first] = E_Entries[j].second;
						Li = LineIndeces[Li] = E_Entries[j].first + 1;
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
							else if (fabs(u)  > std::min(tau/NuL,tau2))
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
#if defined(ILUC2)
						if (fabs(LineValues[Li - 1]) > tau)
						{
							curr = Li;
							for (INMOST_DATA_ENUM_TYPE ru = U2_Address[Li - 1].first; ru != U2_Address[Li - 1].last; ++ru)
							{
								u = LineValues[Li - 1] * LU2_Entries[ru].second;
								j = LU2_Entries[ru].first;
								if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
									LineValues[j] -= u;
								else if (fabs(u)  > std::min(tau/NuL,tau2))
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
						}
#endif
						Li = LineIndeces[Li];
					}
					
					//drop values
					Li = LineIndeces[cbeg];
					Ui = cbeg;
					while (Li != EOL)
					{
						j = LineIndeces[Li];
						if (fabs(LineValues[Li - 1]) < std::min(tau/NuL,tau2) )//tau2)
						{
							LineIndeces[Ui] = j;
							LineIndeces[Li] = UNDEF;
							LineValues[Li - 1] = 0.0;
						}
						else Ui = Li;
						Li = j;
					}
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
						for (INMOST_DATA_ENUM_TYPE r = LFt_Address[Li - 1].first; r < LFt_Address[Li - 1].last; ++r)
						{
							j = LFt_Entries[r].first;
							u = LFt_Entries[r].second;
							l = LineValues[Li - 1];
							if (Ulist[j] != UNDEF)
								temp[j] -= u*l;
							else if( fabs(u) > std::max(tau/NuU,tau2) || fabs(l) > std::max(tau/NuL,tau2) )
							{
								temp[j] = -u*l;
								Ulist[j] = Ubeg[cbeg];
								Ubeg[cbeg] = j;
							}
						}
						Li = LineIndeces[Li];
					}
					tt1 = Timer() - tt1;

					tt2 = Timer();
					INMOST_DATA_REAL_TYPE Smax = 1.0e-54;
					INMOST_DATA_REAL_TYPE Snorm = 0.0, Snum = 0.0;
					Ui = Ubeg[cbeg];
					while (Ui != EOL)
					{
						u = fabs(temp[Ui]);
						if( u > Smax ) Smax = u;
						Snorm += u*u;
						Snum++;
						Ui = Ulist[Ui];
					}
					Snorm = sqrt(Snorm/Snum);
					//insert obtained row
					Ui = Ubeg[cbeg];
					S_Address[i].first = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
					while (Ui != EOL)
					{
						//drop values
						if( fabs(temp[Ui])/Smax > std::max(tau2,std::min(tau/NuU,tau/NuL)) * Smax / Snorm)
						//if( fabs(temp[Ui]) > tau * Smax )// Snorm)
							S_Entries.push_back(Solver::Row::make_entry(Ui, temp[Ui]));
						Ui = Ulist[Ui];
					}
					S_Address[i].last = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
					//clean after use
					Ui = Ubeg[cbeg];
					while (Ui != EOL)
					{
						j = Ulist[Ui];
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
						printf("%lu %d/%d schur2 %6.2f%%\t\t\r", static_cast<unsigned long>(level_size.size()),cend,moend,100.f*(i - cend) / (1.f*(wend - cend)));
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
			INMOST_DATA_REAL_TYPE SnormF = 0.0;
			for(k = 0; k < A_Entries.size(); ++k) SnormF += A_Entries[k].second*A_Entries[k].second;
			SnormF = sqrt(SnormF);
			std::cout << std::endl << " ||S||f: " << SnormF << std::endl;


			printf("\ntotal %f phase1 %f (%6.2f%%) phase2 %f (%6.2f%%) nnz %d\n", t0, t1, 100.0*t1 / t0, t2, 100.0*t2 / t0, A_Entries.size());
#endif
			tlschur = Timer() - tt;
			tschur += tlschur;
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
			
			if ((1.0 * nzA / ((wend - wbeg)*(wend - wbeg))) > 0.75 && (1.0*(wend - wbeg)) / (1.0 *(moend - mobeg)) > 0.1)
			{
				std::cout << "Try to sparsify schur complement!!!" << std::endl;
				std::cout << "Sparsity: " << 1.0 * nzA / ((wend - wbeg)*(wend - wbeg)) << std::endl;
				std::cout << "Schur size: " << wend - wbeg << " total size: " << moend - mobeg << std::endl;
			}
			
			if ( (wend-wbeg > 0 && wend - wbeg < 64) || (1.0 * nzA / ((wend - wbeg)*(wend - wbeg))) > 0.75)
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
				
				if( entries == NULL )
				{
					std::cout << __FILE__ << ":" << __LINE__ << " array of size " << size*size << " was not allocated " << std::endl;
					throw -1;
				}
				//memset(entries,0,sizeof(INMOST_DATA_REAL_TYPE)*size*size);
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
					if( estimator )
					{
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
						NuU_old = NuU;
						if (np > nm) NuU = mup; else NuU = mum;
						for (Ui = i + 1; Ui < size; Ui++) EstU[wbeg + Ui] += entries[ADDR(i,Ui)] * NuU;
						NuU = std::max(fabs(mup), fabs(mum));
						NuU = std::max(NuU,NuU_old);
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
						NuL_old = NuL;
						if (np > nm) NuL = mup; else NuL = mum;
						for (Li = i + 1; Li < size; Li++) EstL[wbeg+Li] += entries[ADDR(Li,i)] * NuL;
						NuL = std::max(fabs(mup), fabs(mum));
						NuL = std::max(NuL,NuL_old);
					}
				}
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " finishing reorder " << std::endl;
#endif
				//apply PQ reordering to other part of matrices
				ReorderEF(mobeg,wbeg,wend,wend,donePQ,localP,localQ);
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
					U_Address[wbeg + i].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					for (j = i + 1; j < size; j++)
					{
						if (fabs(entries[ADDR(i, j)]) > tau / std::max(1.0, NuU))
						{
							LU_Entries.push_back(Solver::Row::make_entry(j + wbeg, entries[ADDR(i, j)])); //Add to U
						}
					}
					U_Address[wbeg + i].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					L_Address[wbeg + i].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					for (j = i + 1; j < size; j++)
					{
						if (fabs(entries[ADDR(j, i)]) > tau / std::max(1.0, NuL))
						{
							LU_Entries.push_back(Solver::Row::make_entry(j + wbeg, entries[ADDR(j, i)])); //Add to L
						}
					}
					L_Address[wbeg + i].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					nzLU += 1 + L_Address[wbeg + i].Size() + U_Address[wbeg + i].Size();
				}
				delete[] entries;
				//Remember the size of the block
				level_size.push_back(wend - wbeg);
#if defined(REPORT_ILU)
				std::cout << "new level " << wend - wbeg << " total nonzeros in LU " << nzLU << " in EF " << nzEF << " conditions L " << NuL << " D " << NuD << " U " << NuU << std::endl;
#endif
				//FINISH!!!
				tfactor += Timer() - tt;
				break;
			}
			

			ratio = (tlrescale + tlfactor)/(tlreorder + tlreassamble + tlschur);
			ratio2 = (trescale + tfactor)/(treorder + treassamble + tschur);
			if( ddpq_tau_adapt )
			{
				if( ratio < 1.0 && ratio2 < 1.0 )
					ddpq_tau /= 1.015;
				else
					ddpq_tau *= 1.015;
				ddpq_tau = std::max(0.01,ddpq_tau);
				ddpq_tau = std::min(0.99,ddpq_tau);
			}
			/*
//#if defined(REPORT_ILU)
			std::cout << "============================================================================" << std::endl;
			//std::cout << "factor: " << tlrescale + tlfactor << " trescale: " << tlrescale << " tfactor " << tlfactor << std::endl;
			//std::cout << "schur:  " << tlreorder + tlreassamble + tlschur << " treorder: " << tlreorder << " treassamble " << tlreassamble << " tschur: " << tlschur << std::endl;
			std::cout << " time ratio " << ratio << std::endl;
			std::cout << " time ratio2 " << ratio2 << std::endl;
			std::cout << " ddpq_tau " << ddpq_tau << std::endl;
			std::cout <<  cbeg << "..." << cend << " of " << mobeg << "..." << moend << " level " << level_size.size() << std::endl;
			std::cout << "============================================================================" << std::endl;
//#endif
			*/
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// DENSE FACTORIZATION COMPLETE //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		ttotal = Timer() - ttotal;
#if defined(REPORT_ILU)
		printf("total      %f   levels %d\n",ttotal,level_size.size());
		printf("reorder    %f (%6.2f%%)\n", treorder, 100.0*treorder/ttotal);
		printf("reassamble %f (%6.2f%%)\n", treassamble, 100.0*treassamble / ttotal);
		printf("rescale    %f (%6.2f%%)\n", trescale, 100.0*trescale / ttotal);
		printf("factor     %f (%6.2f%%)\n", tfactor, 100.0*tfactor / ttotal);
		printf("schur      %f (%6.2f%%)\n", tschur, 100.0*tschur / ttotal);
#endif
		init = true;
		level_interval.resize(level_size.size());
		level_interval[0].first = mobeg;
		level_interval[0].last = mobeg + level_size[0];
		for (k = 1; k < level_size.size(); k++)
		{
			level_interval[k].first = level_interval[k - 1].last;
			level_interval[k].last = level_interval[k].first + level_size[k];
		}
		return true;
	}
	bool Finalize()
	{
		init = false;
		L_Address.clear();
		U_Address.clear();
		LU_Entries.clear();
		for (INMOST_DATA_ENUM_TYPE k = 0; k < E_Address.size(); k++)
			delete E_Address[k];
		E_Address.clear();
		E_Entries.clear();
		F_Entries.clear();
		F_Address.clear();
		B_Entries.clear();
		B_Address.clear();
		temp.clear();
		ddP.clear();
		ddQ.clear();
		level_size.clear();
		LU_Diag.clear();
		return true;
	}
	void Multiply(int level, Solver::Vector & input, Solver::Vector & output)
	{
		INMOST_DATA_ENUM_TYPE k, m, cbeg, cend;
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		for (k = cbeg; k < cend; k++)
		{
			output[k] = 0.0;
			for (m = B_Address[k].first; m < B_Address[k].last; ++m)
				output[k] += input[B_Entries[m].first] * B_Entries[m].second;
		}
	}
	int Descend(int level, Solver::Vector & input, Solver::Vector & output)
	{
		INMOST_DATA_ENUM_TYPE k, m, cbeg, cend;
		//current B block interval
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		// step 1) obtain ~f
		//apply B^{-1} on f
		for (k = cbeg; k < cend; ++k) output[k] = input[k];
		for (k = cbeg; k < cend; ++k) temp[k] = input[k];
		//Solve with L first
		for (k = cbeg; k < cend; ++k) //iterator over columns of L
		{
			for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r != L_Address[k].last; ++r)
				temp[LU_Entries[r].first] -= temp[k] * LU_Entries[r].second;
		}
		//Solve with diagonal
		for (k = cbeg; k < cend; ++k) temp[k] /= LU_Diag[k];
		//Solve with U
		for (k = cend; k > cbeg; --k) //iterator over rows of U
		{
			for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r != U_Address[k - 1].last; ++r)
				temp[k - 1] -= temp[LU_Entries[r].first] * LU_Entries[r].second;
		}
		//now can calculate ~g
		//multiply ~f by row of E
		for (k = cend; k < level_interval.back().last; ++k)
		{
			for (m = E_Address[level]->at(k).first; m < E_Address[level]->at(k).last; ++m)
				output[k] -= temp[E_Entries[m].first] * E_Entries[m].second;
		}
		return level + 1;
	}
	int Ascend(int level, Solver::Vector & input, Solver::Vector & output)
	{
		level = level - 1;
		INMOST_DATA_ENUM_TYPE cbeg,cend, k, m;
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		for (k = cbeg; k < cend; ++k) output[k] = input[k];
		//calculate Fy, iterate over ~g, write to unused input vector
		for (k = cbeg; k < cend; ++k)
		{
			for (m = F_Address[k].first; m < F_Address[k].last; m++)
				output[k] -= output[F_Entries[m].first] * F_Entries[m].second;
		}
		//perform solve over calculated vector
		//Solve with L first
		for (k = cbeg; k < cend; ++k) //iterator over columns of L
		{
			for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r != L_Address[k].last; ++r)
				output[LU_Entries[r].first] -= output[k] * LU_Entries[r].second; //r->first alwayse > k
		}
		//Solve with diagonal
		for (k = cbeg; k < cend; ++k) output[k] /= LU_Diag[k];
		//Solve with U
		for (k = cend; k > cbeg; --k) //iterator over rows of U
		{
			for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r != U_Address[k - 1].last; ++r)
				output[k - 1] -= output[LU_Entries[r].first] * LU_Entries[r].second; // r->first always > k
		}
		return level;
	}
	bool Solve(Solver::Vector & input, Solver::Vector & output)
	{
		assert(&input != &output);
		//
		INMOST_DATA_ENUM_TYPE k, level, mobeg, moend, vbeg, vend;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetVectorRegion(vbeg, vend);
		for (k = mobeg; k < moend; k++) temp[k] = input[k];
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
		for (k = vbeg; k < mobeg; k++) output[k] = 0;
		for (k = mobeg; k < moend; ++k) output[k] = temp[ddP[k]];
		for (k = moend; k < vend; k++) output[k] = 0;

		level = 0;
		//perform recursively first two steps of solve phase
		while(level < level_size.size()) level = Descend(level, output, output);
		

		//level = Ascend(level, output);
		//level = Ascend(level, output);
		//level = Descend(level, output);
		//level = Descend(level, output);

		//on the last recursion level third step was accomplished
		//finish last step in unwinding recursion
		
		while (level > 0) level = Ascend(level, output, output);

		
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
	bool ReplaceSOL(Solver::Vector & x) {(void)x; return true; }
	bool ReplaceRHS(Solver::Vector & b) {(void)b; return true; }
	Method * Duplicate() { return new ILUC_preconditioner(*this); }
	~ILUC_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}
};

#endif //__SOLVER_DDPQILUC2__
