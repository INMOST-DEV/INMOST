#define _CRT_SECURE_NO_WARNINGS
#include "inmost_solver.h"
#if defined(USE_SOLVER)
#include "solver_mlmtiluc2.hpp"
#include <sstream>
#include <deque>
#include <iomanip>

#include "../../../Misc/utils.h"
//#define USE_OMP

using namespace INMOST;

#ifndef DEFAULT_TAU
#define DEFAULT_TAU 0.01
#endif


#define REORDER_RCM
//#define REORDER_NNZ
#if defined(USE_SOLVER_METIS)
#define REORDER_METIS_ND
#endif
#if defined(USE_SOLVER_MONDRIAAN)
//#define REORDER_MONDRIAAN
#endif
//#define REORDER_ZOLTAN_HUND


static bool run_mpt = true;
static bool rescale_b = true;
static bool allow_pivot = true;

#define ESTIMATOR
#define ESTIMATOR_REFINE

//~ #define PREMATURE_DROPPING

//#define EQUALIZE_1NORM
//~ #define EQUALIZE_2NORM
#define EQUALIZE_IDOMINANCE

#define PIVOT_THRESHOLD
#define PIVOT_THRESHOLD_VALUE 1.0e-12
//#define DIAGONAL_PERTURBATION
#define DIAGONAL_PERTURBATION_REL 1.0e-7
#define DIAGONAL_PERTURBATION_ABS 1.0e-10
#define ILUC2
#define ILUC2_SCHUR
#define TRACK_DIAGONAL
//#define PIVOT_COND_DEFAULT 0.1/tau
#define PIVOT_COND_DEFAULT 1.0e+2
#define PIVOT_DIAG_DEFAULT 1.0e+5
//~ #define SCHUR_DROPPING_LF
//~ #define SCHUR_DROPPING_EU
//~ #define SCHUR_DROPPING_S
#define DIAGONAL_PIVOT
#define CONDITION_PIVOT

#if defined(REORDER_METIS_ND)
#define METIS_EXPORT
#include "metis.h"
#endif
#if defined(REORDER_ZOLTAN_HUND)
#include <zoltan.h>
#endif
#if defined(REORDER_MONDRIAAN)
#include <Mondriaan.h>
#endif




	void MLMTILUC_preconditioner::ReorderEF(INMOST_DATA_ENUM_TYPE wbeg,
											INMOST_DATA_ENUM_TYPE wend,
											interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
	{
		INMOST_DATA_ENUM_TYPE i, k, l;
		if( !E_Address.empty() )
		{
			for(k = wbeg; k < wend; ++k) donePQ[k] = false;
			i = wbeg;
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
							for(l = 0; l < (int)E_Address.size(); ++l)
							{
								Interval t = E_Address[l]->at(i);
								E_Address[l]->at(i) = E_Address[l]->at(localP[k]);
								E_Address[l]->at(localP[k]) = t;
							}
							donePQ[localP[k]] = true;
							k = localP[k];
						} while (k != i);
					}
					donePQ[i] = true;
				}
			}
		}
		if( !F_Address.empty() )
		{
			for(k = wbeg; k < wend; ++k) donePQ[k] = false;
			i = wbeg;
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
							for(l = 0; l < (int)F_Address.size(); ++l)
							{
								Interval t = F_Address[l]->at(i);
								F_Address[l]->at(i) = F_Address[l]->at(localQ[k]);
								F_Address[l]->at(localQ[k]) = t;
							}
							donePQ[localQ[k]] = true;
							k = localQ[k];
						} while (k != i);
					}
					donePQ[i] = true;
				}
			}
		}
	}


	void MLMTILUC_preconditioner::DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
											 std::vector<Sparse::Row::entry> & Entries,
											 INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
											 std::string file_name)
	{
		INMOST_DATA_REAL_TYPE norm1 = 0, norm2 = 0, max = -1.0e54, min = 1.0e54, minabs = 1.0e54;
		INMOST_DATA_REAL_TYPE vrowmax, diag, mindiag = 1.0e54, maxdiag = -1.0e54, maxabsdiag = -1.0e54, minabsdiag = 1.0e54;
		INMOST_DATA_ENUM_TYPE nnz = 0, dominant_rows = 0, dominant_cols = 0, irowmax = 0, nodiag = 0;
		INMOST_DATA_ENUM_TYPE addrbeg = Address.get_interval_beg();
		INMOST_DATA_ENUM_TYPE addrend = Address.get_interval_end();
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> vcolmax(wmbeg,wmend,0);
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> icolmax(wmbeg,wmend,ENUMUNDEF);
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k) 
		{
			if( k < addrbeg || k >= addrend) continue;
			nnz += Address[k].Size();
			vrowmax = 0;

			bool diag_found = false;
			diag = 0;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it)
			{
				norm1 += fabs(Entries[it].second);
				norm2 += Entries[it].second*Entries[it].second;
				if( fabs(Entries[it].second) > vrowmax ) 
				{
					vrowmax = fabs(Entries[it].second);
					irowmax = Entries[it].first;
				}

				if( fabs(Entries[it].second) > vcolmax[Entries[it].first] ) 
				{
					vcolmax[Entries[it].first] = fabs(Entries[it].second);
					icolmax[Entries[it].first] = k;
				}
				if( Entries[it].second > max ) max = Entries[it].second;
				if( Entries[it].second < min ) min = Entries[it].second;
				if( fabs(Entries[it].second) < fabs(minabs) ) minabs = Entries[it].second;

				if( Entries[it].first == k )
				{
					diag_found = true;
					diag = Entries[it].second;
				}
			}

			if( diag_found )
			{
				if( mindiag > diag ) mindiag = diag;
				if( maxdiag < diag ) maxdiag = diag;
				if( minabsdiag > fabs(diag) ) minabsdiag = fabs(diag);
				if( maxabsdiag < fabs(diag) ) maxabsdiag = fabs(diag);
			}
			else nodiag++;

			if( irowmax == k ) ++dominant_rows;
		}

		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k) if( icolmax[k] == k ) ++dominant_cols;

		std::cout << "Writing matrix to " << file_name.c_str() << std::endl;
		std::fstream fout(file_name.c_str(),std::ios::out);
		fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
		fout << "% maximum " << max << std::endl;
		fout << "% minimum " << min << std::endl;
		fout << "% absolute minimum " << minabs << std::endl;
		fout << "% A 1-norm  " << norm1 << std::endl;
		fout << "% A 2-norm  " << sqrt(norm2) << std::endl;
		fout << "% mean 1-norm  " << norm1/(wmend-wmbeg) << std::endl;
		fout << "% mean 2-norm  " << sqrt(norm2/(wmend-wmbeg)) << std::endl;
		fout << "% dominant rows  " << dominant_rows << std::endl;
		fout << "% dominant cols  " << dominant_cols << std::endl;
		fout << "% maximal diagonal value " << maxdiag << std::endl;
		fout << "% minimal diagonal value " << mindiag << std::endl;
		fout << "% absolute maximal diagonal value " << maxabsdiag << std::endl;
		fout << "% absolute minimal diagonal value " << minabsdiag << std::endl;
		fout << "% no diagonal value  " << nodiag << std::endl;
		fout << "% true matrix indices interval " << wmbeg << ":" << wmend << std::endl;
		fout << std::scientific;
		
		//fout.close(); return;
		
		fout << wmend-wmbeg << " " << wmend-wmbeg << " " << nnz << std::endl;;
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k)
		{
			if( k < addrbeg || k >= addrend) continue;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it)
				fout << (k-wmbeg+1) << " " << (Entries[it].first-wmbeg+1) << " " << Entries[it].second << std::endl;
		}
		fout.close();
	}
	void MLMTILUC_preconditioner::CheckOrder(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
											 std::vector<Sparse::Row::entry> & Entries,
											 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend)
	{
		INMOST_DATA_ENUM_TYPE i,r;
		for (i = rbeg; i < rend; i++) if( Address[i].first < Address[i].last )
		{
			//check ordered on entry
			bool ordered = true;
			for (r = Address[i].first; r < Address[i].last-1; r++)
			{
				if( !(Entries[r].first < Entries[r+1].first) )
				{
					ordered = false;
					break;
				}
			}
			if( !ordered )
			{
				std::cout << "Row " << i << " not ordered: " << std::endl;
				std::cout << "Interval: " << Address[i].first << ":" << Address[i].last << std::endl;
				for (r = Address[i].first; r < Address[i].last; r++)
				{
					std::cout << "(" << Entries[r].first << "," << Entries[r].second << ") ";
				}
				std::cout << std::endl;
				throw -1;
			}
			bool nan = false;
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				if( Entries[r].second != Entries[r].second )
				{
					nan = true;
					break;
				}
			}
			if( nan )
			{
				std::cout << "Row " << i << " contains nan: " << std::endl;
				std::cout << "Interval: " << Address[i].first << ":" << Address[i].last << std::endl;
				for (r = Address[i].first; r < Address[i].last; r++)
				{
					std::cout << "(" << Entries[r].first << "," << Entries[r].second << ") ";
				}
				std::cout << std::endl;
				throw -1;
			}
		}
	}
	
	void MLMTILUC_preconditioner::inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
											INMOST_DATA_ENUM_TYPE wend,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		//inverse reordering
		// in invPQ numbers indicate where to get current column
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invP[localP[k]] = k;
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invQ[localQ[k]] = k;
	}
	void MLMTILUC_preconditioner::applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
										  INMOST_DATA_ENUM_TYPE wend,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		INMOST_DATA_ENUM_TYPE k;
		// compute reordering in global P,Q, we need it to compute reordering in vector during solve phase
		for (k = wbeg; k < wend; ++k)
		{
			localP[k] = ddP[invP[k]];
			localQ[k] = ddQ[invQ[k]];
		}
		// update reordering in global P,Q
		for (k = wbeg; k < wend; ++k)
		{
			ddP[k] = localP[k];
			ddQ[k] = localQ[k];
		}
	}
	INMOST_DATA_ENUM_TYPE & MLMTILUC_preconditioner::EnumParameter(std::string name)
	{
		if (name == "scale_iters") return sciters;
		else if( name == "estimator" ) return estimator;
		else if( name == "verbosity" ) return verbosity;
		throw - 1;
	}
	INMOST_DATA_REAL_TYPE & MLMTILUC_preconditioner::RealParameter(std::string name)
	{
		if (name == "tau") return tau;
		else if( name == "tau2" ) return iluc2_tau;
		else if( name == "pivot_cond" ) return pivot_cond;
		else if( name == "pivot_diag" ) return pivot_diag;
		else if( name == "condition_number_L" ) return condestL;
		else if( name == "condition_number_U" ) return condestU;
		throw - 1;
	}
	void MLMTILUC_preconditioner::Copy(const Method * other)
	{
		const MLMTILUC_preconditioner * b = dynamic_cast<const MLMTILUC_preconditioner *>(other);
		assert(b != NULL);
		tau = b->tau;
		iluc2_tau = b->iluc2_tau;
		pivot_cond = b->pivot_cond;
		pivot_diag = b->pivot_diag;
		Alink = b->Alink;
		info = b->info;
		sciters = b->sciters;
		eps = b->eps;
	}
	MLMTILUC_preconditioner::MLMTILUC_preconditioner(const MLMTILUC_preconditioner & other) :Method(other)
	{
		Copy(&other);
	}
	MLMTILUC_preconditioner & MLMTILUC_preconditioner::operator =(MLMTILUC_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	MLMTILUC_preconditioner::MLMTILUC_preconditioner(Solver::OrderInfo & info)
		:tau(DEFAULT_TAU),info(&info)
	{
		Alink = NULL;
		init = false;
		sciters = 8;
		eps = 1e-54;
#if defined(ESTIMATOR)
		estimator = 1;
#else
		estimator = 0;
#endif
		tau = 1.0e-3;
		iluc2_tau = tau*tau;
		pivot_cond = PIVOT_COND_DEFAULT;
		pivot_diag = PIVOT_DIAG_DEFAULT;
	}
	bool MLMTILUC_preconditioner::isInitialized() { return init; }
	bool MLMTILUC_preconditioner::isFinalized() { return !init; }
	bool MLMTILUC_preconditioner::Initialize()
	{
		
        const INMOST_DATA_REAL_TYPE subst = 1.0; (void)subst;
		const INMOST_DATA_REAL_TYPE tol_modif = PIVOT_THRESHOLD_VALUE;
		const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;
		bool block_pivot = false;
		
		INMOST_DATA_ENUM_TYPE wbeg, wend; //working interval
		INMOST_DATA_ENUM_TYPE mobeg, moend; // total interval
		INMOST_DATA_ENUM_TYPE vbeg, vend; // vector interval
		
		INMOST_DATA_ENUM_TYPE k, i, j, Li, Ui, curr, next;
		INMOST_DATA_REAL_TYPE l,u,udiag, abs_udiag, max_diag = 0, min_diag = 0;
		INMOST_DATA_REAL_TYPE max_diag_old, min_diag_old;
		INMOST_DATA_ENUM_TYPE nzA, nzLU = 0, nzA0;
		Sparse::Vector DL, DR;
		Sparse::Vector DL0,DR0;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetVectorRegion(vbeg,vend);
		
		//prepare temporal array
		temp.set_interval_beg(vbeg);
		temp.set_interval_end(vend);
		//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE>::iterator rit = temp.begin();
		//	rit != temp.end(); ++rit) *rit = 0.0;
		for(k = vbeg; k < vend; ++k) temp[k] = 0.0;

		
		//prepare reordering vectors
		ddP.set_interval_beg(mobeg);
		ddP.set_interval_end(moend);
		ddQ.set_interval_beg(mobeg);
		ddQ.set_interval_end(moend);

		

		//prepare rescaling vectors
		DL.SetInterval(mobeg, moend);
		DR.SetInterval(mobeg, moend);
		DL0.SetInterval(mobeg, moend);
		DR0.SetInterval(mobeg, moend);
		for(k = mobeg; k < moend; ++k) DL[k] = DR[k] = DL0[k] = DR0[k] = 1.0;
		//for(Sparse::Vector::iterator ri = DL.Begin(); ri != DL.End(); ++ri) *ri = 1.0;
		//for(Sparse::Vector::iterator ri = DR.Begin(); ri != DR.End(); ++ri) *ri = 1.0;
		//for(Sparse::Vector::iterator ri = DL0.Begin(); ri != DL0.End(); ++ri) *ri = 1.0;
		//for(Sparse::Vector::iterator ri = DR0.Begin(); ri != DR0.End(); ++ri) *ri = 1.0;

		

		for (k = mobeg; k < moend; k++) ddP[k] = ddQ[k] = k;
		// supplementary data for ddPQ reordering
		interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(mobeg, moend, false);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(mobeg, moend), invQ(mobeg,moend), localP(mobeg, moend,ENUMUNDEF), localQ(mobeg,moend,ENUMUNDEF);		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Bstart(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, bool> Pivot(mobeg,moend,false);

		//supplimentary data structures for ILUC
		INMOST_DATA_ENUM_TYPE LU_Beg, Sbeg;
		U_Address.set_interval_beg(mobeg);
		U_Address.set_interval_end(moend);
		L_Address.set_interval_beg(mobeg);
		L_Address.set_interval_end(moend);
		B_Address.set_interval_beg(mobeg);
		B_Address.set_interval_end(moend);
		LU_Diag.set_interval_beg(mobeg);
		LU_Diag.set_interval_end(moend);
		std::vector<Sparse::Row::entry> A_Entries, S_Entries, LF_Entries, LFt_Entries;
		std::vector<Sparse::Row::entry> EU_Entries;
		//std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend), S_Address(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, Interval> LF_Address(mobeg,moend), LFt_Address(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, Interval> EU_Address(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ulist(mobeg, moend), Llist(mobeg, moend), Blist(mobeg,moend), Flist(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(mobeg, moend,EOL), Lbeg(mobeg, moend,EOL), Bbeg(mobeg,moend,EOL), Fbeg(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Scolnorm(mobeg,moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Scolnum(mobeg,moend,0.0);

		INMOST_DATA_REAL_TYPE NuU = 1, NuL = 1, NuD = 1, NuU_max = 1.0, NuL_max = 1.0;
		INMOST_DATA_REAL_TYPE NuU_acc = 1, NuL_acc = 1, NuD_acc = 1;
#if defined(ESTIMATOR)
		//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
		INMOST_DATA_REAL_TYPE mup, mum, smup, smum, NuL1 = 1, NuL2 = 1, NuU1 = 1, NuU2 = 1;
		INMOST_DATA_REAL_TYPE NuU1_old = 1, NuL1_old = 1, NuU2_old = 1, NuL2_old = 1, NuD_old = 1;
		INMOST_DATA_REAL_TYPE NuU1_new = 1, NuL1_new = 1, vp, vm, v;
#if defined(ESTIMATOR_REFINE)
		INMOST_DATA_ENUM_TYPE np, nm;
		INMOST_DATA_REAL_TYPE NuU2_new, NuL2_new;
#endif
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL1(mobeg, moend,0.0), EstU1(mobeg, moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL2(mobeg, moend,0.0), EstU2(mobeg, moend,0.0);
#endif
#if defined(ESTIMATOR)
		INMOST_DATA_REAL_TYPE NuU_tmp, NuL_tmp;
#endif
		//supplimentary data structures for returning values of dropped elements
		//INMOST_DATA_REAL_TYPE DropLk, DropUk;
		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
		//data structure for linked list

		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Dist(mobeg,moend+1,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValuesU(mobeg, moend,0.0), LineValuesL(mobeg,moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndecesU(mobeg, moend+1,UNDEF), LineIndecesL(mobeg,moend+1,UNDEF);
		
		std::vector<INMOST_DATA_ENUM_TYPE> indicesU, indicesL, indicesS;
        double tfactor = 0.0, trescale = 0.0, treorder = 0.0, ttransversal = 0.0, treassamble = 0.0, ttotal, tt, testimator = 0.0, tschur = 0.0, tlocal;
#if defined(REORDER_METIS_ND)
		double tmetisgraph = 0, tmetisnd = 0;
#endif
#if defined(REORDER_RCM)
		double trcmgraph = 0, trcmorder = 0;
#endif
        double tlfactor, tlrescale, tlreorder, tlreassamble, tlschur;
		ttotal = Timer();

		//(*Alink).Save("M.mtx");
		
		//calculate number of nonzeros
		nzA = 0;
		for (k = mobeg; k < moend; ++k)
		{
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0) nzA++;
		}
		nzA0 = nzA;

		//sort_indeces.reserve(256);
		A_Entries.resize(nzA);
		//B_Entries.reserve(nzA);
		//LU_Entries.reserve(nzA*4);

#if defined(REORDER_NNZ)
		wgt_coords sort_wgts(2*(moend - mobeg));
#endif

#if defined(REORDER_ZOLTAN_HUND)
#endif

#if defined(REORDER_MONDRIAAN)
		{
			sparsematrix sp;
			int numparts = 8;
			int * cols = new int[nzA];
			int * rows = new int[Alink->Size()+1];
			double * values = new double[nzA];
			int cnt = 0;
			rows[0] = 0;
			for (k = mobeg; k < moend; ++k)
			{
				for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				{
					if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
					{
						cols[cnt] = r->first;
						values[cnt] = r->second;
						++cnt;
					}
				}
				rows[k-mobeg+1] = cnt;
			}
			CRSSparseMatrixInit(&sp,Alink->Size(),Alink->Size(),nzA,cols,rows,values,0);
			
			PstartInit(&sp,numparts);
			opts options;
			
			SetDefaultOptions(&options);
			//SetOption(&options,"SplitMethod","KLFM");
			SetOption(&options,"Permute","BBD");
			SetOption(&options,"SplitStrategy","onedimrow");
			DistributeMatrixMondriaan(&sp,numparts,0.1,&options,NULL);
			sp.

			j = 0;
			for (k = mobeg; k < moend; ++k)
			{
				A_Address[k].first = j;
				for (Sparse::Row::iterator r = (*Alink)[sp.row_perm[k-mobeg]].Begin(); r != (*Alink)[sp.row_perm[k-mobeg]].End(); ++r)
				{
					if (r->first >= mobeg && r->first < moend && r->second != 0.0)
						A_Entries[j++] = Sparse::Row::make_entry(sp.col_perm_inv[r->first-mobeg], r->second);
				}
				A_Address[k].last = j;
				assert(A_Address[k].Size() != 0); //singular matrix
			}
			//DumpMatrix(A_Address,A_Entries,mobeg,moend,"mondriaan.mtx");
			/*
			sp.MMTypeCode[0] = 'M';
			FILE * f = fopen("mondriaan.mtx","w");
			MMWriteSparseMatrix(&sp,f,"mondriaan",&options);
			fclose(f);
			*/
			MMDeleteSparseMatrix(&sp);
			delete [] cols;
			delete [] rows;
			delete [] values;
		}
#else
		
		j = 0;
		for (k = mobeg; k < moend; ++k)
		{
			A_Address[k].first = j;
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
			{
				if (r->first >= mobeg && r->first < moend && r->second != 0.0)
					A_Entries[j++] = Sparse::Row::make_entry(r->first, r->second);
			}
			A_Address[k].last = j;
			//assert(A_Address[k].Size() != 0); //singular matrix
		}
#endif
		//DumpMatrix(A_Address, A_Entries, mobeg, moend, "A.mtx");

		std::vector<INMOST_DATA_REAL_TYPE> C_Entries(A_Entries.size());
		INMOST_DATA_REAL_TYPE Unorm, Unum, Lnorm, Lnum;
		INMOST_DATA_REAL_TYPE LFnorm, LFnum, LFmax, LFmin, LFtau, LFdrop;
		INMOST_DATA_REAL_TYPE EUnorm, EUnum, EUmax, EUmin, EUtau, EUdrop;
		INMOST_DATA_REAL_TYPE Snorm, Snum, Smax, Smin, Stau, Sdrop;

		INMOST_DATA_ENUM_TYPE nzLU2 = 0, nzLU2tot = 0, ndrops = 0;
#if defined(ILUC2)
		
		INMOST_DATA_REAL_TYPE tau2 = iluc2_tau;
		std::vector<Sparse::Row::entry> LU2_Entries;
		interval<INMOST_DATA_ENUM_TYPE, Interval> L2_Address(mobeg, moend+1), U2_Address(mobeg, moend+1);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2list(mobeg, moend, UNDEF), L2list(mobeg, moend, UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2beg(mobeg, moend, EOL), L2beg(mobeg, moend, EOL);
		//LU2_Entries.reserve(nzA*4);
#else
		INMOST_DATA_REAL_TYPE tau2 = tau;
#endif

		for(k = mobeg; k < moend; ++k)
		{
			U_Address[k].first = U_Address[k].last = ENUMUNDEF;
			L_Address[k].first = L_Address[k].last = ENUMUNDEF;
			U2_Address[k].first = U2_Address[k].last = ENUMUNDEF;
			L2_Address[k].first = L2_Address[k].last = ENUMUNDEF;
		}

		if( verbosity > 1 ) 
			std::cout << "nonzeros in A " << nzA << " pivot cond " << pivot_cond << " diag " << pivot_diag << " tau " << tau << " tau2 " << tau2 << std::endl;

		
		int swaps = 0, totswaps = 0;
		//set up working interval
		wbeg = mobeg;
		wend = moend;
		while (wbeg < wend) //iterate into levels until all matrix is factored
		{
			//ddPQ reordering on current Schur complement
			INMOST_DATA_ENUM_TYPE cbeg = wbeg, cend = wend; //next size of factored B block
			
		//DumpMatrix(A_Address,A_Entries,cbeg,cend,"mat"+to_string(level_size.size())+".mtx");
///////////////////////////////////////////////////////////////////////////////////
///   MAXIMUM TRANSVERSE REORDERING                                             ///
///////////////////////////////////////////////////////////////////////////////////
			if( run_mpt )
			{
				if( verbosity > 1 )
					printf("Reordering with MPT\n");

				ttransversal = Timer();
				INMOST_DATA_ENUM_TYPE ColumnBegin;
				INMOST_DATA_ENUM_TYPE pop_heap_pos;
				INMOST_DATA_REAL_TYPE ShortestPath, AugmentPath;
				INMOST_DATA_ENUM_TYPE PathEnd, Trace, IPermPrev;
				//Sparse::Vector & U = DL;
				//Sparse::Vector & V = DR;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Cmax(wbeg,wend,0.0);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Perm = localQ;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & IPerm = localP;
				//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & UpdateStack = Ulist;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & ColumnList = Llist;//(wbeg,wend,ENUMUNDEF);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Parent = Blist;
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> AugmentPosition(wbeg,wend,ENUMUNDEF);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ColumnPosition(wbeg,wend,ENUMUNDEF);
				
				
				BinaryHeap Heap(&Dist[wbeg],wend-wbeg);
///////////////////////////////////////////////////////////////////////////////////
///  Arrays initialization                                                      ///
///////////////////////////////////////////////////////////////////////////////////
				// arrays U,V,Dist are cleared at the end of schur complement calculation
				//std::fill(U.begin() + wbeg - mobeg, U.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(V.begin() + wbeg - mobeg, V.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(Cmax.begin() + wbeg - mobeg, Cmax.begin() + wend - mobeg, 0.0);
				//std::fill(Dist.begin() + wbeg - mobeg, Dist.begin() + wend + 1 - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				for(k = wbeg; k < wend; ++k)
				{
					Perm[k] = ENUMUNDEF;
					IPerm[k] = ENUMUNDEF;
					Parent[k] = ENUMUNDEF;
					ColumnList[k] = ENUMUNDEF;
				}
				C_Entries.resize(A_Entries.size());
///////////////////////////////////////////////////////////////////////////////////
///       Initial LOG transformation to dual problem and initial extreme match  ///
///////////////////////////////////////////////////////////////////////////////////
				//double T = Timer();
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						u = C_Entries[it] = fabs(A_Entries[it].second);
						if( u > Cmax[i] ) Cmax[i] = u;
						//C_Entries.push_back(Sparse::Row::make_entry(i,u));
					}
				}

				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						if( Cmax[i] == 0 || C_Entries[it] == 0 )
							C_Entries[it] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
						else
						{
							C_Entries[it] = log(Cmax[i])-log(C_Entries[it]);
							//C_Entries[it] = log10(Cmax[i])-log10(C_Entries[it]);
							//C_Entries[it] = fabs(log10(Cmax[i]/C_Entries[it]));
							//C_Entries[it] = fabs(log(Cmax[i]/C_Entries[it]));
							if( C_Entries[it] < U[i] ) U[i] = C_Entries[it];
						}
					}
				}
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						u = C_Entries[it] - U[A_Entries[it].first];
						if( u < V[k] ) V[k] = u;
					}
				}

///////////////////////////////////////////////////////////////////////////////////
///                  Update cost and match                                      ///
///////////////////////////////////////////////////////////////////////////////////
				
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
					{
						u = fabs(C_Entries[it] - V[k] - U[A_Entries[it].first]);
						if( u < 1.0e-30 && Perm[A_Entries[it].first] == ENUMUNDEF && IPerm[k] == ENUMUNDEF )
						{
							 Perm[A_Entries[it].first] = k;
							 IPerm[k] = A_Entries[it].first;
							 ColumnPosition[k] = it;
						}
					}
				}
				
///////////////////////////////////////////////////////////////////////////////////
/// 1-step augmentation                                                         ///
///////////////////////////////////////////////////////////////////////////////////

				for(k = wbeg; k < wend; ++k)
				{
					if( IPerm[k] == ENUMUNDEF ) //unmatched row
					{
						for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last && IPerm[k] == ENUMUNDEF; ++it)
						{
							u = fabs(C_Entries[it] - V[k] - U[A_Entries[it].first]);
							if( u <= 1.0e-30 )
							{
								Li = Perm[A_Entries[it].first];
								assert(Li != ENUMUNDEF);
								// Search other row in C for 0
								for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit < A_Address[Li].last; ++Lit)
								{
									u = fabs(C_Entries[Lit]- V[Li] - U[A_Entries[Lit].first]);
									if( u <= 1.0e-30 && Perm[A_Entries[Lit].first] == ENUMUNDEF )
									{
										Perm[A_Entries[it].first] = k;
										IPerm[k] = A_Entries[it].first;
										ColumnPosition[k] = it;
										Perm[A_Entries[Lit].first] = Li;
										IPerm[Li] = A_Entries[Lit].first;
										ColumnPosition[Li] = Lit;
										break;
									}
								}
							}
						}
					}
				}

///////////////////////////////////////////////////////////////////////////////////
///             Weighted bipartite matching                                     ///
///////////////////////////////////////////////////////////////////////////////////
				for(k = wbeg; k < wend; ++k)
				{
					if( IPerm[k] != ENUMUNDEF )
						continue;
					Li = k;
					ColumnBegin = EOL;
					Parent[Li] = ENUMUNDEF;
					PathEnd = ENUMUNDEF;
					Trace = k;
					ShortestPath = 0;
					AugmentPath = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					while(true)
					{
						for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit < A_Address[Li].last; ++Lit)
						{
							Ui = A_Entries[Lit].first;
							//if( ColumnList[Ui] == k ) continue;
							if( ColumnList[Ui] != ENUMUNDEF ) continue;
							l = fabs(ShortestPath + C_Entries[Lit] - V[Li] - U[Ui]);
							//if( l < 0.0 ) printf("row %d col %d negative l %g Augment %lf Shortest %lf C %lf V %lf U %lf\n",k,Ui,l,AugmentPath,ShortestPath,C_Entries[Lit],V[Li],U[Ui]);
							if( l < 0.0 && l > -1.0e-8 ) l = 0;
							if( l < 0.0 ) continue;
							if( l < AugmentPath )
							{
								if( Perm[Ui] == ENUMUNDEF )
								{
									PathEnd = Ui;
									Trace = Li;
									AugmentPath = l;
									AugmentPosition[Ui] = Lit;
								}
								else if( l < Dist[Ui] )
								{
									Parent[Perm[Ui]] = Li;
									AugmentPosition[Ui] = Lit;
									if( Heap.Contains(Ui-wbeg) )
										Heap.DecreaseKey(Ui-wbeg,l);
									else
										Heap.PushHeap(Ui-wbeg,l);
								}
							}
						}

						pop_heap_pos = Heap.PopHeap();
						if( pop_heap_pos == ENUMUNDEF ) break;
					
						Ui = pop_heap_pos+wbeg;
						ShortestPath = Dist[Ui];

						if( AugmentPath <= ShortestPath ) 
						{
							Dist[Ui] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
							//Heap.increaseKey(Ui,Dist[Ui]);
							break;
						}


						ColumnList[Ui] = ColumnBegin;
						ColumnBegin = Ui;

						Li = Perm[Ui];
						
					}
					if( PathEnd != ENUMUNDEF )
					{
						Ui = ColumnBegin;
						while(Ui != EOL)
						{
							U[Ui] += Dist[Ui] - AugmentPath;
							if( Perm[Ui] != ENUMUNDEF ) V[Perm[Ui]] = C_Entries[ColumnPosition[Perm[Ui]]] - U[Ui];
							Dist[Ui] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
							Li = ColumnList[Ui];
							ColumnList[Ui] = ENUMUNDEF;
							Ui = Li;
						}

						Ui = PathEnd;
						while(Trace != ENUMUNDEF)
						{
							IPermPrev = IPerm[Trace];
							Perm[Ui] = Trace;
							IPerm[Trace] = Ui;

							ColumnPosition[Trace] = AugmentPosition[Ui];
							V[Trace] = C_Entries[ColumnPosition[Trace]] - U[Ui];

							Ui = IPermPrev;
							Trace = Parent[Trace];

						}
						Heap.Clear();
					}
				}
				//printf("Maximum product transversal %lf\n",Timer()-T);
				//fclose(rec);

				//T = Timer();
				for (k = cbeg; k < cend; ++k)
				{
					if( V[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() ) l = 1;
					else l = exp(V[k]);
					if( U[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || Cmax[k] == 0 ) u = 1;
					else u = exp(U[k])/Cmax[k];
					//if( l != l || fabs(l) < 1.0e-12 || isnan(l) || isinf(l) ) std::cout << "k " << k << " l is " << l << " V " << V[k] << std::endl;
					//if( u != u || fabs(u) < 1.0e-12 || isnan(u) || isinf(u) ) std::cout << "k " << k << " u is " << u << " U " << U[k] << " Cmax " << Cmax[k] << std::endl;
					DL[k] = l;
					DR[k] = u;
					

					bool flip_sign = false;
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
					{
						i = A_Entries[jt].first;
						j = Perm[A_Entries[jt].first];
						
						if( U[i] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || Cmax[i] == 0 ) u = 1;
						else u = exp(U[i])/Cmax[i];
						
						if( fabs(l*A_Entries[jt].second*u) > 1 + 1.0e-7 )
							std::cout << "element on row " << k << " col " << A_Entries[jt].first << " value " << A_Entries[jt].second << " u " << u << " l " << l << " U " << U[i] << " V " << V[k] << " Cmax " << Cmax[k] << " scaled " << l*A_Entries[jt].second*u << std::endl;
						if( j == k )
						{
							
							if( l*A_Entries[jt].second*u < 0.0 ) flip_sign = true;
						}
					}

					if( flip_sign ) DL[k] *= -1;
				}
				//printf("Reorder matrix %lf\n",Timer()-T);
				/*
				{
					std::stringstream s;
					s << "DL" << level_size.size() << ".mtx";
					DL.Save(s.str());
				}
				{
					std::stringstream s;
					s << "DR" << level_size.size() << ".mtx";
					DR.Save(s.str());
				}
				 */
				/*
				{
					std::stringstream s;
					s << "MPT" << level_size.size() << ".txt";
					std::fstream fout(s.str().c_str(),std::ios::out);
					for(k = cbeg; k < cend; ++k)
					{
						fout << k << " " << DL[k] << " " << DR[k] << " " << localQ[k] << std::endl;
					}
					fout.close();
				}
				*/
				//std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
				for(k = wbeg; k < wend; ++k) localP[k] = ENUMUNDEF;
                
                { //check that there are no gaps in Perm
                    for(k = cbeg; k < cend; ++k)
                    {
                        if( Perm[k] != ENUMUNDEF )
						{
							assert(localP[Perm[k]] == ENUMUNDEF);
                            localP[Perm[k]] = 0;
						}
                    }
					std::vector<INMOST_DATA_ENUM_TYPE> gaps;
                    for(k = cbeg; k < cend; ++k)
                        if( localP[k] == ENUMUNDEF )
                            gaps.push_back(k);
					
					//std::cout << "@ gaps: " << gaps.size() << std::endl;
                    
                    for(k = cbeg; k < cend; ++k)
                        if( Perm[k] == ENUMUNDEF )
                        {
                            Perm[k] = gaps.back();
                            gaps.pop_back();
                        }
                    //std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
                    for(k = wbeg; k < wend; ++k) localP[k] = ENUMUNDEF;
                }
				
				
				
				//exit(-1);

				ttransversal = Timer() - ttransversal;

				treorder += ttransversal;
				if( verbosity > 1 )
					printf("Time %g\n",ttransversal);

			}
			else
			{
				for(k = wbeg; k < wend; ++k)
				{
					localQ[k] = k;
					DL[k] = DR[k] = 1;
				}
			}
///////////////////////////////////////////////////////////////////////////////////
///                  END MAXIMUM TRANSVERSE REORDERING                          ///
///////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
#if defined(REORDER_METIS_ND)
			idx_t nvtxs = wend-wbeg;
			std::vector<idx_t> xadj(nvtxs+1), adjncy, perm(nvtxs),iperm(nvtxs);
			//adjncy.reserve(nzA*2);
			if( verbosity > 1 )
				printf("Reordering with METIS\n");

			
			for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) printf("%s:%d No column permutation for row %d. Matrix is structurally singular\n",__FILE__,__LINE__,k);

			for (k = cbeg; k < cend; ++k) localP[k] = k;
			for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
			for (k = cbeg; k < cend; ++k) DR[k] = V[k];
			
			for (k = cbeg; k < cend; ++k)
			{
				if( !(localQ[k] >= cbeg && localQ[k] < cend) )
					std::cout << "Bad permutation: " << localQ[k] << " interval [" << cbeg << ":" << cend << "]" << std::endl;
				V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
			}
			for (k = cbeg; k < cend; ++k) DR0[k] = V[k];

			
			
			for(k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
				{
					A_Entries[jt].first = localQ[A_Entries[jt].first];
					//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					
				}
				std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
			}
			ReorderEF(wbeg, wend, donePQ, localP, localQ);
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);

			//DumpMatrix(A_Address,A_Entries,wbeg,wend,"a.mtx");

///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of B block                       //
///////////////////////////////////////////////////////////////////////////////////
			
			tmetisgraph = Timer();

			for (k = wend; k > wbeg; --k)
			{
				//vwgt[k-1] = 1;
				if (A_Address[k-1].Size() > 0)
				{
					i = A_Entries[A_Address[k-1].first].first;
					Blist[k-1] = Bbeg[i];
					Bbeg[i] = k-1;
				}
				Ulist[k-1] = A_Address[k-1].first;
			}
			

			xadj[0] = 0;
			for(i = wbeg; i < wend; ++i)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
				{
					if( A_Entries[jt].first != i )
						adjncy.push_back(static_cast<idx_t>(A_Entries[jt].first-wbeg));
				}
				
				Li = Bbeg[i];
				while (Li != EOL)
				{
					if( Li != i )
						adjncy.push_back(static_cast<idx_t>(Li-wbeg));
					Li = Blist[Li];
				}
				
				Li = Bbeg[i];
				while (Li != EOL)
				{
					Ui = Blist[Li];
					Ulist[Li]++;
					if (A_Address[Li].last - Ulist[Li] > 0)
					{
						k = A_Entries[Ulist[Li]].first;
						Blist[Li] = Bbeg[k];
						Bbeg[k] = Li;
					}

					Li = Ui;
				}
				
				std::sort(adjncy.begin()+xadj[i-wbeg],adjncy.end());
				adjncy.resize(std::unique(adjncy.begin()+xadj[i-wbeg],adjncy.end())-adjncy.begin());
				
				xadj[i-wbeg+1] = static_cast<idx_t>(adjncy.size());
			}

			//std::fill(Bbeg.begin(),Bbeg.end(),EOL);
			for(i = mobeg; i < moend; ++i) Bbeg[i] = EOL;

			tmetisgraph = Timer() - tmetisgraph;
			/*
			{
				FILE * f = fopen("metis.mtx","w");
				fprintf(f,"%%MatrixMarket matrix coordinate real general\n");
				fprintf(f,"%d %d %d\n",nvtxs,nvtxs,adjncy.size());
				for(i = wbeg; i < wend; ++i)
				{
					for(j = xadj[i]; j < xadj[i+1]; ++j)
						fprintf(f,"%d %d 1\n",i+1,adjncy[j]+1);
				}
				fclose(f);
			}
			*/
			//A_Address[0].first = 0;
			//for(k = wbeg+1; k < wend; ++k)
			//	A_Address[k].first = A_Address[k-1].last;
			if( !adjncy.empty() )
			{
				idx_t options[METIS_NOPTIONS];
				METIS_SetDefaultOptions(options);
				options[METIS_OPTION_NUMBERING] = 0;
				
				//options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_TIME | METIS_DBG_COARSEN | METIS_DBG_CONNINFO | METIS_DBG_CONTIGINFO | METIS_DBG_IPART | METIS_DBG_MEMORY | METIS_DBG_MOVEINFO | METIS_DBG_REFINE | METIS_DBG_SEPINFO;
				//options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
				//options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED;
				//options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
				//options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
				//options[METIS_OPTION_NO2HOP] = 0;
				//options[METIS_OPTION_MINCONN] = 0;
				//options[METIS_OPTION_NITER] = 4;
				//printf("before metis\n");
				//METIS_NodeNDP(nvtxs,&xadj[0],&adjncy[0],NULL,4,options,&perm[0],&iperm[0],&sizes[0]);
				tmetisnd = Timer();
				METIS_NodeND(&nvtxs,&xadj[0],&adjncy[0],NULL,options,&perm[0],&iperm[0]);
				tmetisnd = Timer()-tmetisnd;
				//printf("after metis\n");
				for(k = wbeg; k < wend; ++k)
				{
					//localP[k] = iperm[k-wbeg]+wbeg;
					//localQ[k] = iperm[k-wbeg]+wbeg;
					localP[perm[k-wbeg]+wbeg] = k;
					localQ[perm[k-wbeg]+wbeg] = k;
				}
			}
			else
			{
				for(k = wbeg; k < wend; ++k)
				{
					localP[k] = k;
					localQ[k] = k;
				}
			}
			cend = wend;
			i = wend;
#elif defined(REORDER_RCM)
			{
				if( verbosity > 1 ) 
					printf("Reordering with RCM\n");
				//create a symmetric graph of the matrix A + A^T
				std::vector<INMOST_DATA_ENUM_TYPE> xadj(wend-wbeg+1), adjncy;
				//adjncy.reserve(nzA*2);

				for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) printf("%s:%d No column permutation for row %d. Matrix is structurally singular\n",__FILE__,__LINE__,k);

				for (k = cbeg; k < cend; ++k) localP[k] = k;
				for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < cend; ++k) DR[k] = V[k];
				
				for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < cend; ++k) DR0[k] = V[k];

				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
					{
						A_Entries[jt].first = localQ[A_Entries[jt].first];
						//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					}
					std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
				}
				
				//DumpMatrix(A_Address,A_Entries,cbeg,cend,"A.mtx");
				
				ReorderEF(wbeg, wend, donePQ, localP, localQ);
				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);

				trcmgraph = Timer();
				
				for(k = wbeg; k < wend; ++k)
				{
					Ulist[k] = EOL;
					Bbeg[k] = EOL;
					Blist[k] = EOL;
				}
				//std::fill(Ulist.begin()+wbeg-mobeg, Ulist.begin()+wend-mobeg,EOL);
				//std::fill(Bbeg.begin()+wbeg-mobeg, Bbeg.begin()+wend-mobeg,EOL);
				//std::fill(Blist.begin()+wbeg-mobeg, Blist.begin()+wend-mobeg,EOL);
				for (k = wend; k > wbeg; --k)
				{
					//vwgt[k-1] = 1;
					if (A_Address[k-1].Size() > 0)
					{
						i = A_Entries[A_Address[k-1].first].first;
						Blist[k-1] = Bbeg[i];
						Bbeg[i] = k-1;
					}
					Ulist[k-1] = A_Address[k-1].first;
				}
				xadj[0] = 0;
				for(i = wbeg; i < wend; ++i)
				{
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
					{
						if( A_Entries[jt].first != i )
							adjncy.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(A_Entries[jt].first));
					}

					Li = Bbeg[i];
					while (Li != EOL)
					{
						if( Li != i )
							adjncy.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(Li));
						Li = Blist[Li];
					}

					Li = Bbeg[i];
					while (Li != EOL)
					{
						Ui = Blist[Li];
						Ulist[Li]++;
						if (A_Address[Li].Size() && A_Address[Li].last - Ulist[Li] > 0)
						{
							k = A_Entries[Ulist[Li]].first;
							Blist[Li] = Bbeg[k];
							Bbeg[k] = Li;
						}

						Li = Ui;
					}

					std::sort(adjncy.begin()+xadj[i-wbeg],adjncy.end());
					adjncy.resize(std::unique(adjncy.begin()+xadj[i-wbeg],adjncy.end())-adjncy.begin());

					xadj[i-wbeg+1] = static_cast<INMOST_DATA_ENUM_TYPE>(adjncy.size());
				}

				//std::fill(Bbeg.begin()+wbeg - mobeg, Bbeg.begin()+wend - mobeg,EOL);
				for(k = wbeg; k < wend; ++k) Bbeg[k] = EOL;

				trcmgraph = Timer()-trcmgraph;

				trcmorder = Timer();
				//std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, ENUMUNDEF);
				for(k = wbeg; k < wend; ++k) Ulist[k] = ENUMUNDEF;
				//find node with the lowest order
                INMOST_DATA_ENUM_TYPE index = wbeg;
				INMOST_DATA_ENUM_TYPE cur = ENUMUNDEF;
				std::deque<INMOST_DATA_ENUM_TYPE> q;
				std::vector<INMOST_DATA_ENUM_TYPE> conns;
                do
				{
					cur = ENUMUNDEF;
					for(k = wbeg; k < wend && cur == ENUMUNDEF; ++k)
					{
						if( Ulist[k] == ENUMUNDEF )
							cur = k;
					}
					assert(cur != ENUMUNDEF);
					for(k = cur+1; k < wend; ++k) if( Ulist[k] == ENUMUNDEF )
					{
						if( RCM_Comparator(wbeg,xadj)(k,cur) )
							cur = k;
					}
					q.push_back(cur);
					Ulist[cur] = index++;
					while(!q.empty())
					{
						cur = q.front();
						q.pop_front();
						for (INMOST_DATA_ENUM_TYPE it = xadj[cur-wbeg]; it < xadj[cur-wbeg+1]; ++it)
							if( Ulist[adjncy[it]] == ENUMUNDEF )
								conns.push_back(adjncy[it]);
						std::sort(conns.begin(),conns.end(),RCM_Comparator(wbeg,xadj));
						for (k = 0; k < static_cast<INMOST_DATA_ENUM_TYPE>(conns.size()); ++k)
						{
							Ulist[conns[k]] = index++;
							q.push_back(conns[k]);
						}
						conns.clear();
					}
				}
				while( index < wend );

				for(k = wbeg; k < wend; ++k)
					Ulist[k] = wend-(Ulist[k]-wbeg)-1;

				for(k = wbeg; k < wend; ++k)
				{
					localP[k] = Ulist[k];//wend - Ulist[k] + 1;
					localQ[k] = Ulist[k];//wend - Ulist[k] + 1;
					//localP[Ulist[k]] = k;
					//localQ[Ulist[k]] = k;
				}
				cend = wend;
				i = wend;


				trcmorder = Timer() - trcmorder;
			}
#else
			cend = wend;
			i = wbeg;
#endif
			tt = Timer() - tt;
			treorder += tt;
			if( verbosity > 1 )
			{
				printf("Time %g\n",tt);
				printf("Reorder\n");
			}
			tt = Timer();
			if (cbeg == cend && cbeg != wend)
			{
				//std::cout << __FILE__ << ":" << __LINE__ << " singular matrix, factored " << mobeg << ".." << cend << " out of " << mobeg << ".." << moend << std::endl;
				for (k = cbeg; k < wend; ++k)
					LU_Diag[k] = tol_modif;
				break;
			}
			
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

			for (k = cbeg; k < cend; ++k) 
			{
				U[localP[k]] = DL[k];
				V[localQ[k]] = DR[k];
			}
			for (k = cbeg; k < cend; ++k) 
			{
				DL[k] = U[k];
				DR[k] = V[k];
			}
			
			for (k = cbeg; k < cend; ++k) 
			{
				if( !(localQ[k] >= cbeg && localQ[k] < cend) ||
				    !(localP[k] >= cbeg && localP[k] < cend))
					std::cout << "Bad permutations P: " << localP[k] << " Q: " << localQ[k] << " interval [" << cbeg << ":" << cend << "]" << std::endl;
				U[localP[k]] = DL0[k];
				V[localQ[k]] = DR0[k];
			}
			for (k = cbeg; k < cend; ++k) 
			{
				DL0[k] = U[k];
				DR0[k] = V[k];
			}
			
			tlreorder = Timer() - tt;
			treorder += tlreorder;
			if( verbosity > 1 )
				printf("Time %g\n",tlreorder);
///////////////////////////////////////////////////////////////////////////////////
///                  REASSAMBLE                                                 ///
///////////////////////////////////////////////////////////////////////////////////
			if( verbosity > 1 ) 
				printf("Reassemble\n");
			//tt = Timer();
			double tt1, tt2,ttt;
			//in localPQ numbers indicate where to put current row/column
			//reorder E,F blocks by swaps
			tt1 = Timer();
			//inverse ordering
			ReorderEF(wbeg,wend, donePQ, localP, localQ);
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			tt1 = Timer() - tt1;

			//std::cout << "reorder: " << tt1 << std::endl;
			tt2 = Timer();
			nzA = 0;
			B_Entries.clear();
			for (k = cbeg; k < cend; ++k)
			{
				B_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt < A_Address[invP[k]].last; ++jt)
				{
					i = localQ[A_Entries[jt].first];
					u = A_Entries[jt].second;
					B_Entries.push_back(Sparse::Row::make_entry(i, u));
				}
				B_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				std::sort(B_Entries.begin() + B_Address[k].first, B_Entries.end());
				nzA += A_Address[k].Size();
			}
			B_Entries.push_back(Sparse::Row::make_entry(-1,0.0));
			
			//std::stringstream s;
			//s << "B" << level_size.size() << ".mtx";
			//DumpMatrix(B_Address,B_Entries,cbeg,cend,s.str());
			/*
			{
				int rank = 0;
#if defined(USE_MPI)
				MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
				std::stringstream str;
				str << "mat_reorder" << rank << ".mtx";
				DumpMatrix(B_Address,B_Entries,cbeg,cend,str.str());
			}
#if defined(USE_MPI)
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			exit(0);
			*/
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

			ttt = Timer();
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);
			tt1 += Timer() - ttt;
			
			//reset localPQ
			
			//for(int k = mobeg; k < moend; ++k) localP[k] = localQ[k] = k;
			/*
			{
				std::stringstream s;
				s << "C" << level_size.size() << ".mtx";
				DumpMatrix(B_Address,B_Entries,cbeg,cend,s.str());
			}
			*/
			tlreassamble = Timer() - tt;
			treassamble += tlreassamble;
			if( verbosity > 1 )
				printf("Time %g\n",tlreassamble);
///////////////////////////////////////////////////////////////////////////////////
///                  RESCALING                                                  ///
///////////////////////////////////////////////////////////////////////////////////
			//Rescale current B block by Row-Column alternating scaling
			tt = Timer();
			if( rescale_b )
			{
				if( verbosity > 1 )
					std::cout << " rescaling block B " << std::endl;

				
#if defined(EQUALIZE_1NORM) || defined(EQUALIZE_2NORM)
				for (k = cbeg; k < cend; k++) DL[k] = DR[k] = 1.0;
#elif defined(EQUALIZE_IDOMINANCE)
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						B_Entries[r].second *= (DL[k] * DR[B_Entries[r].first]);
					
					U[k] = DL[k];
					V[k] = DR[k];
				}
#else
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						B_Entries[r].second *= (DL[k] * DR[B_Entries[r].first]);
				}
#endif
				/*
				{
					std::stringstream s;
					s << "CC" << level_size.size() << ".mtx";
					DumpMatrix(B_Address,B_Entries,cbeg,cend,s.str());
				}
				*/
				//DumpMatrix(B_Address,B_Entries,cbeg,cend,"MC64_MPTILUC2.mtx");
				/*
				{
					std::stringstream name;
					name << "mc64_" << info->GetRank() << ".mtx";
					DumpMatrix(B_Address,B_Entries,cbeg,cend,name.str());
				}
				*/
				
///////////////////////////////////////////////////////////////////////////////////
///                COMPUTE GIRSCHGORIN RADIUS                                   ///
///////////////////////////////////////////////////////////////////////////////////
				/*
				INMOST_DATA_REAL_TYPE radii = 0;
				for (k = cbeg; k < cend; k++)
				{
					INMOST_DATA_REAL_TYPE local_radii = 0;
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						if( B_Entries[r].first != k ) local_radii += fabs(B_Entries[r].second);
					if( radii < local_radii ) radii = local_radii;
				}

				printf("Gershgorin's radius after mc64: %e\n",radii);
				*/


				//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_mc64.mtx");
#if defined(EQUALIZE_1NORM)
///////////////////////////////////////////////////////////////////////////////////
///         ROW-COLUMN ALTERNATING SCALING FOR 1-NORM BALANCING                 ///
///////////////////////////////////////////////////////////////////////////////////
				//std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
				for(k = cbeg; k < cend; ++k) DL[k] = 0.0;
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						DL[k] += fabs(B_Entries[r].second);//*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
				for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
				{
					//std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
					for(k = cbeg; k < cend; ++k) DR[k] = 0.0;
					for (k = cbeg; k < cend; k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
							DR[B_Entries[r].first] += DL[k] * fabs(B_Entries[r].second);//*B_Entries[r].second;
					}
					for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
					//std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
					for(k = cbeg; k < cend; ++k) DL[k] = 0.0;
					for (k = cbeg; k < cend; k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
							DL[k] += DR[B_Entries[r].first] * fabs(B_Entries[r].second);//*B_Entries[r].second;
					}
					for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
				}
				//for (k = cbeg; k < cend; k++) DL[k] = sqrt(DL[k]);
				//for (k = cbeg; k < cend; k++) DR[k] = sqrt(DR[k]);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];
				}
#elif defined(EQUALIZE_2NORM)
				//std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
				for(k = cbeg; k < cend; ++k) DL[k] = 0.0;
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						DL[k] += B_Entries[r].second*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
				for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
				{
					//std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
					for(k = cbeg; k < cend; ++k) DR[k] = 0.0;
					for (k = cbeg; k < cend; k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
							DR[B_Entries[r].first] += DL[k] * B_Entries[r].second*B_Entries[r].second;
					}
					for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
					//std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
					for(k = cbeg; k < cend; ++k) DL[k] = 0.0;
					for (k = cbeg; k < cend; k++)
					{
						for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
							DL[k] += DR[B_Entries[r].first] * B_Entries[r].second*B_Entries[r].second;
					}
					for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
				}
				//for (k = cbeg; k < cend; k++) DL[k] = sqrt(DL[k]);
				//for (k = cbeg; k < cend; k++) DR[k] = sqrt(DR[k]);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];
				}
#elif defined(EQUALIZE_IDOMINANCE)
///////////////////////////////////////////////////////////////////////////////////
///        THIS VERSION OF RESCALING INCREASES DIAGONAL DOMINANCE               ///
///////////////////////////////////////////////////////////////////////////////////
				if( C_Entries.size() < B_Entries.size() )
					C_Entries.resize(B_Entries.size());
				
				for(k = cbeg; k < cend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					{
						u = fabs(B_Entries[r].second);
						if( u > 1.0e-12 )
							C_Entries[r] = -log(u);
						else
							C_Entries[r] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					}
				}

				
				//std::fill(temp.begin() + cbeg - mobeg, temp.begin() + cend - mobeg, 0.0);
				for(k = cbeg; k < cend; ++k) temp[k] = 0.0;

				for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
				{
					//std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
					//std::fill(DR.Begin() + cbeg - mobeg, DR.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
					for(k = cbeg; k < cend; ++k) DL[k] = DR[k] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					for (k = cbeg; k < cend; k++) //row number
					{
						for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						{
							i = B_Entries[r].first; //column number
							if( i != k ) //out of diagonal
							{
								u = C_Entries[r] + temp[k] - temp[i];
								//if( isnan(u) || u != u ) std::cout << __FILE__ << ":" << __LINE__ << " u is " << u << std::endl;
								DL[k] = std::min(DL[k],u);// update Y1
								DR[i] = std::min(DR[i],u);// update Y2
							}
						}
					}

					for (k = cbeg; k < cend; k++)
					{
						if( DR[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() &&
							DL[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() )
							temp[k] += (DR[k]-DL[k])*0.5;
					}

				}

				for (k = cbeg; k < cend; k++)
				{
					DL[k] = exp(-temp[k]);
					DR[k] = exp(temp[k]);
					//if( isnan(DL[k]) || DL[k] != DL[k] || fabs(DL[k]) < 1.0e-12 ) std::cout << __FILE__ << ":" << __LINE__ << " DL[" << k << "] is " << DL[k] << std::endl;
					//if( isnan(DR[k]) || DR[k] != DR[k] || fabs(DR[k]) < 1.0e-12 ) std::cout << __FILE__ << ":" << __LINE__ << " DR[" << k << "] is " << DR[k] << std::endl;
				}

				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						B_Entries[r].second *= (DL[k] * DR[B_Entries[r].first]);
				}
#endif
///////////////////////////////////////////////////////////////////////////////////
///              COMPUTE GIRSCHGORIN RADIUS                                     ///
///////////////////////////////////////////////////////////////////////////////////
				/*
				radii = 0;
				for (k = cbeg; k < cend; k++)
				{
					INMOST_DATA_REAL_TYPE local_radii = 0;
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
						if( B_Entries[r].first != k ) local_radii += fabs(B_Entries[r].second);
					if( radii < local_radii ) radii = local_radii;
				}

				printf("Gershgorin's radius after equilibration: %e\n",radii);
				*/
#if defined(EQUALIZE_IDOMINANCE)
				for (k = cbeg; k < cend; k++)
				{
					DL[k] *= U[k];
					DR[k] *= V[k];
				}
#endif
				//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_equilibration.mtx");
///////////////////////////////////////////////////////////////////////////////////
///               RESCALING DONE                                                ///
///////////////////////////////////////////////////////////////////////////////////
				
				//stack scaling
				
				for (k = cbeg; k < cend; k++)
				{
					DL0[k] *= DL[k];
					DR0[k] *= DR[k];
				}
				/*
				DL.Save("DL_"+to_string(level_size.size())+".mtx");
				DR.Save("DR_"+to_string(level_size.size())+".mtx");
				DL0.Save("DL0_"+to_string(level_size.size())+".mtx");
				DR0.Save("DR0_"+to_string(level_size.size())+".mtx");
				DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_equil"+to_string(level_size.size())+".mtx");
				 */
			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"rescale"+to_string(level_size.size())+".mtx");
				//exit(-1);
				//rescale EF
				if( !level_size.empty() )
				{
					INMOST_DATA_ENUM_TYPE first = mobeg, last;
					INMOST_DATA_ENUM_TYPE kbeg, kend, r;
					for(size_t level = 0; level < level_size.size(); ++level)
					{
						last = first + level_size[level];
						kbeg = std::max(cbeg,last);
						kend = std::min(cend,moend);
						//std::cout << "Rescale EF level " << level << " ";
						//std::cout << "interval [" << first << "," << last << "] ";
						//std::cout << "matrix [" << mobeg << "," << moend << "] ";
						//std::cout << "rescaling [" << cbeg << "," << cend << "] ";
						//std::cout << "selected [" << kbeg << ":" << kend << "]" << std::endl;
						for (k = kbeg; k < kend; ++k) 
						{
							if( F_Address[level]->at(k).Size() )
							{
								for (r = F_Address[level]->at(k).first; r < F_Address[level]->at(k).last; ++r)
									F_Entries[r].second *= DR[k];
							}
							if( E_Address[level]->at(k).Size() )
							{
								for (r = E_Address[level]->at(k).first; r < E_Address[level]->at(k).last; ++r)
									E_Entries[r].second *= DL[k];
							}
						}
						first = last;
					}
				}
				//End rescale B block
				tlrescale = Timer() - tt;
				trescale += tlrescale;
				if( verbosity > 1 )
					printf("Time %g\n",tlrescale);

			}
			/*
			{
				std::stringstream s;
				s << "B" << level_size.size() << ".mtx";
				DumpMatrix(B_Address,B_Entries,cbeg,cend,s.str());
			}
			*/
///////////////////////////////////////////////////////////////////////////////////
///           FACTORIZATION BEGIN                                               ///
///////////////////////////////////////////////////////////////////////////////////
			{
				
#if defined(ILUC2)
				nzLU2tot += nzLU2;
				nzLU2 = 0;
				LU_Beg = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
				LU2_Entries.clear();
#endif
///////////////////////////////////////////////////////////////////////////////////
//       setup diagonal values and memorize indices                              //
///////////////////////////////////////////////////////////////////////////////////
				//initialize diagonal values
				for (k = cbeg; k < cend; k++)
				{
					LU_Diag[k] = 0.0;
					Bstart[k] = B_Address[k].first; //remember B block starts permutation for pivoting
					for (i = B_Address[k].first; i < B_Address[k].last; i++)
					{
						if (B_Entries[i].first == k)
						{
#if defined(TRACK_DIAGONAL)
							LU_Diag[k] = B_Entries[i].second;
#endif
							break;
						}
					}
				}
				//Initialize data for condition estimator
				tt = Timer();
///////////////////////////////////////////////////////////////////////////////////
//              initialize condition estimator                                   //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
				if( estimator )
				{
					NuU_acc *= NuU;
					NuL_acc *= NuL;
					NuD_acc *= NuD;
					NuU = NuL = 1.0;
					NuU1 = NuL1 = 1.0;
					NuU2 = NuL2 = 1.0;
					NuU_max = NuL_max = 1.0;
				}
#endif
				max_diag = 0;
				min_diag = 1.0e300;
				NuD = 1.0;
				if( verbosity > 1 )
					std::cout << " starting factorization " << std::endl;

				for (k = cbeg; k < cend; ++k)
				{
///////////////////////////////////////////////////////////////////////////////////
//              starting factorization step                                      //
///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
//            Prepare diagonal value                                             //
///////////////////////////////////////////////////////////////////////////////////
#if defined(TRACK_DIAGONAL)
					udiag = LU_Diag[k]; // LU_Diag is calculated with less dropping
#if defined(DIAGONAL_PIVOT)
					if( allow_pivot && !block_pivot && fabs(udiag)*pivot_diag < 1 )
					{
						//std::cout << "k " << k << " diag " << udiag << " cond " << pivot_diag << std::endl;
						Pivot[k] = true;
						swaps++;
						Li = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
						Ui = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
						U_Address[k].first = U_Address[k].last = Li;
						L_Address[k].first = L_Address[k].last = Li;
						U2_Address[k].first = U2_Address[k].last = Ui;
						L2_Address[k].first = L2_Address[k].last = Ui;
						LineIndecesU[cbeg] = LineIndecesL[cbeg] = EOL;
						goto exit_point;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
// Prevent small pivot, hopefully it would work with correct right hand side     //
///////////////////////////////////////////////////////////////////////////////////
#if defined(PIVOT_THRESHOLD)
					if (fabs(udiag) < tol_modif)
					{
						udiag = udiag < 0.0 ? -tol_modif : tol_modif;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
// Perturb diagonal value, this will shift eigenvalues                           //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PERTURBATION)
					udiag = udiag * (1.0 + DIAGONAL_PERTURBATION_REL) + (udiag < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS;
#endif
					//abs_udiag = mean_diag/static_cast<INMOST_DATA_REAL_TYPE>(cend-k);//fabs(LU_Diag[k]);
					abs_udiag = fabs(LU_Diag[k]);
#else
					abs_udiag = 1.0;
#endif
///////////////////////////////////////////////////////////////////////////////////
//            Prepare linked list for U-part                                     //
///////////////////////////////////////////////////////////////////////////////////

					//DropLk = DropUk = 0.0;
					//uncompress k-th row
					// add diagonal value first, there shouldn't be values on left from diagonal
					//assert(B_Entries[B_Address[k].first].first == k);
					LineIndecesU[cbeg] = k;
					LineIndecesU[k] = EOL;
					LineValuesU[k] = 0.0;
					if( B_Address[k].Size() )
					{
						if (B_Entries[B_Address[k].first].first == k)
							LineValuesU[k] = B_Entries[B_Address[k].first].second;
						else
							LineValuesU[k] = 0.0;
						Ui = k;
						for (INMOST_DATA_ENUM_TYPE it = B_Address[k].first + (B_Entries[B_Address[k].first].first == k ? 1 : 0); it < B_Address[k].last; ++it)
						{
							LineValuesU[B_Entries[it].first] = B_Entries[it].second;
							Ui = LineIndecesU[Ui] = B_Entries[it].first;
						}
						LineIndecesU[Ui] = EOL;
					}
					Sbeg = k;
///////////////////////////////////////////////////////////////////////////////////
//                    U-part elimination with L                                  //
///////////////////////////////////////////////////////////////////////////////////
					// perform the elimination from row
					i = Lbeg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU_Entries[L_Address[i].first].first == k);
						l = LU_Entries[L_Address[i].first].second*LU_Diag[i];
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > 1.0e-9 )
#endif
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
							{
								assert(!Pivot[i]); // pivoted rows should be empty
								j = LU_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								u = l*LU_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else if( u )
								{
									/*
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesU[curr] = j;
									LineIndecesU[j] = next;
									LineValuesU[j] = -u;
									*/
									
									LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									 
								}
							}
						}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with L                        //
///////////////////////////////////////////////////////////////////////////////////
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > 1.0e-9 )
#endif
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								assert(!Pivot[i]); // pivoted rows should be empty
								j = LU2_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								u = l*LU2_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else if( u )
								{
									/*
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesU[curr] = j;
									LineIndecesU[j] = next;
									LineValuesU[j] = -u;
									*/
									
									LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									 
								}
							}
						}
#endif
						i = Llist[i];
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 U-part elimination with second-order L                        //
///////////////////////////////////////////////////////////////////////////////////
					i = L2beg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU2_Entries[L2_Address[i].first].first == k);
						l = LU2_Entries[L2_Address[i].first].second*LU_Diag[i];
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > 1.0e-9 )
#endif
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
							{
								assert(!Pivot[i]); // pivoted rows should be empty
								j = LU_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								u = l*LU_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else if( u )
								{
									/*
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesU[curr] = j;
									LineIndecesU[j] = next;
									LineValuesU[j] = -u;
									*/
									
									LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									 
								}
							}
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with second-order L           //
///////////////////////////////////////////////////////////////////////////////////
#if 0
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								assert(!Pivot[i]); // pivoted rows should be empty
								j = LU2_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								u = l*LU2_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else if( u )
								{
									
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									/*
									assert(curr < j);
									assert(j < next);
									LineIndecesU[curr] = j;
									LineIndecesU[j] = next;
									LineValuesU[j] = -u;
									*/
									
									LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									 
								}
							}
#endif
						}
						i = L2list[i];
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//                  Order contents of linked list                                //
///////////////////////////////////////////////////////////////////////////////////
					
					indicesU.clear();
					Ui = Sbeg;
					while(Ui != EOL)
					{
						indicesU.push_back(Ui);
						Ui = LineIndecesU[Ui];
					}
					if( !indicesU.empty() )
					{
						std::sort(indicesU.begin(),indicesU.end());
						//Sbeg = indices[0];
						for(size_t qt = 1; qt < indicesU.size(); ++qt)
							LineIndecesU[indicesU[qt-1]] = indicesU[qt];
						LineIndecesU[indicesU.back()] = EOL;
					}
					
///////////////////////////////////////////////////////////////////////////////////
//                  Retrieve diagonal value                                      //
///////////////////////////////////////////////////////////////////////////////////
#if !defined(TRACK_DIAGONAL)
					udiag = LineValuesU[k];
#if defined(DIAGONAL_PIVOT)
					if( allow_pivot && !block_pivot && fabs(udiag)*pivot_diag < 1 )
					{
						Pivot[k] = true;
						swaps++;
						Li = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
						Ui = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
						U_Address[k].first = U_Address[k].last = Li;
						L_Address[k].first = L_Address[k].last = Li;
						U2_Address[k].first = U2_Address[k].last = Ui;
						L2_Address[k].first = L2_Address[k].last = Ui;
						LineIndecesL[cbeg] = EOL;
						goto exit_point;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
// Prevent small pivot, hopefully it would work with correct right hand side     //
///////////////////////////////////////////////////////////////////////////////////
#if defined(PIVOT_THRESHOLD)
					if (fabs(udiag) < tol_modif)
					{
						//std::cout << "udiag too small " << 1.0 /udiag << std::endl;
						udiag = udiag < 0.0 ? -tol_modif : tol_modif;
					}
#endif
#if defined(DIAGONAL_PERTURBATION)
					udiag = udiag * (1.0 + DIAGONAL_PERTURBATION_REL) + (udiag < 0.0 ? -1.0 : 1.0)*DIAGONAL_PERTURBATION_ABS;
#endif
#endif
///////////////////////////////////////////////////////////////////////////////////
//                    Condition estimator for diagonal D                         //
///////////////////////////////////////////////////////////////////////////////////
					max_diag_old = max_diag;
					min_diag_old = min_diag;
					if (fabs(udiag) > max_diag) max_diag = fabs(udiag);
					if (fabs(udiag) < min_diag) min_diag = fabs(udiag);
					NuD_old = NuD;
					NuD = max_diag / min_diag;
///////////////////////////////////////////////////////////////////////////////////
//                    Rescale with diagonal                                     //
///////////////////////////////////////////////////////////////////////////////////
					//rescale elements next after diagonal position
					Ui = LineIndecesU[k];
					while (Ui != EOL)
					{
						LineValuesU[Ui] /= udiag;
						Ui = LineIndecesU[Ui];
					}
///////////////////////////////////////////////////////////////////////////////////
//                 prepearing linked list for L-part                             //
///////////////////////////////////////////////////////////////////////////////////
					//uncompress column
					//insert diagonal value first
					LineIndecesL[cbeg] = k;
					LineIndecesL[k] = EOL;
					LineValuesL[k] = 0.0;
					if( B_Address[k].Size() )
					{
						if (B_Entries[B_Address[k].first].first == k)
							LineValuesL[k] = B_Entries[B_Address[k].first].second;
						else
							LineValuesL[k] = 0.0;
						//start from diagonal
						Ui = k;
						Li = Bbeg[k];
						while (Li != EOL)
						{
							assert(B_Entries[B_Address[Li].first].first == k);
							LineValuesL[Li] = B_Entries[B_Address[Li].first].second;
							Ui = LineIndecesL[Ui] = Li;
							Li = Blist[Li];
						}
						LineIndecesL[Ui] = EOL;
					}
					Sbeg = k;
///////////////////////////////////////////////////////////////////////////////////
//                 L-part elimination with U                                     //
///////////////////////////////////////////////////////////////////////////////////
					// perform the elimination from column
					i = Ubeg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU_Entries[U_Address[i].first].first == k);
						u = LU_Entries[U_Address[i].first].second*LU_Diag[i];
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > 1.0e-9 )
#endif
						{
							curr = cbeg;
							assert(!Pivot[i]);
							for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it < L_Address[i].last; ++it)
							{
								j = LU_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								if( Pivot[j] ) continue;
								l = u*LU_Entries[it].second;
								//eliminate values
								if (LineIndecesL[j] != UNDEF) //there is an entry in the list
									LineValuesL[j] -= l;
								else if( l )
								{
									/*
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesL[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesL[curr] = j;
									LineIndecesL[j] = next;
									LineValuesL[j] = -l;
									*/
									
									LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									
								}
							}
						}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second-order L-part elimination with U                        //
///////////////////////////////////////////////////////////////////////////////////
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > 1.0e-9 )
#endif
						{
							curr = cbeg;
							assert(!Pivot[i]);
							for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								if( Pivot[j] ) continue;
								l = u*LU2_Entries[it].second;
								if (LineIndecesL[j] != UNDEF) //there is an entry in the list
									LineValuesL[j] -= l;
								else if( l )
								{
									/*
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesL[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesL[curr] = j;
									LineIndecesL[j] = next;
									LineValuesL[j] = -l;
									*/
									
									LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									
								}
							}
						}
#endif
						i = Ulist[i];
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 L-part elimination with second-order U                        //
///////////////////////////////////////////////////////////////////////////////////
					i = U2beg[k];
					while (i != EOL)
					{
						assert(i != UNDEF);
						assert(LU2_Entries[U2_Address[i].first].first == k);
						u = LU2_Entries[U2_Address[i].first].second*LU_Diag[i];
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > 1.0e-9 )
#endif
						{
							curr = cbeg;
							assert(!Pivot[i]);
							for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it < L_Address[i].last; ++it)
							{
								j = LU_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								if( Pivot[j] ) continue;
								l = u*LU_Entries[it].second;
								//eliminate values
								if (LineIndecesL[j] != UNDEF) //there is an entry in the list
									LineValuesL[j] -= l;
								else if( l )
								{
									/*
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesL[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesL[curr] = j;
									LineIndecesL[j] = next;
									LineValuesL[j] = -l;
									*/
									
									LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									
								}
							}
///////////////////////////////////////////////////////////////////////////////////
//           second-order L-part elimination with second-order U                 //
///////////////////////////////////////////////////////////////////////////////////
#if 0
							curr = cbeg;
							assert(!Pivot[i]);
							for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								assert(j >= k); // all indices are supposed to be to the right
								if( Pivot[j] ) continue;
								l = u*LU2_Entries[it].second;
								//eliminate values
								if (LineIndecesL[j] != UNDEF) //there is an entry in the list
									LineValuesL[j] -= l;
								else if( l )
								{
									/*
									next = curr;
									while (next < j)
									{
										curr = next;
										next = LineIndecesL[curr];
									}
									assert(curr < j);
									assert(j < next);
									LineIndecesL[curr] = j;
									LineIndecesL[j] = next;
									LineValuesL[j] = -l;
									*/
									
									LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									
								}
							}
#endif
						}
						i = U2list[i];
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//                  Order contents of linked list                                //
///////////////////////////////////////////////////////////////////////////////////
					
					indicesL.clear();
					Ui = Sbeg;
					while(Ui != EOL)
					{
						indicesL.push_back(Ui);
						Ui = LineIndecesL[Ui];
					}
					if( !indicesL.empty() )
					{
						std::sort(indicesL.begin(),indicesL.end());
						//Sbeg = indices[0];
						for(size_t qt = 1; qt < indicesL.size(); ++qt)
							LineIndecesL[indicesL[qt-1]] = indicesL[qt];
						LineIndecesL[indicesL.back()] = EOL;
					}
					
					//check that diagonal value is the same(must be!)
					//assert(fabs(LineValues[k] / udiag - 1.0) < 1.0e-10);
///////////////////////////////////////////////////////////////////////////////////
//                    Rescale with diagonal                                      //
///////////////////////////////////////////////////////////////////////////////////
					//rescale line by diagonal
					Li = LineIndecesL[k];
					while (Li != EOL)
					{
						LineValuesL[Li] /= udiag;
						Li = LineIndecesL[Li];
					}
///////////////////////////////////////////////////////////////////////////////////
//                    Condition estimator for U part                             //
///////////////////////////////////////////////////////////////////////////////////
					//estimate condition for U^{-1}
#if defined(ESTIMATOR)
					if( estimator )
					{
						tlocal = Timer();
						NuU1_old = NuU1;
						mup = 1.0 - EstU1[k];
						mum = -1.0 - EstU1[k];
						smup = smum = 0;
						Ui = LineIndecesU[k];
						while (Ui != EOL)
						{
							v = EstU1[Ui];
							vp = fabs(v + LineValuesU[Ui] * mup);
							vm = fabs(v + LineValuesU[Ui] * mum);
							smup += vp;
							smum += vm;
							Ui = LineIndecesU[Ui];
						}
						if (smup > smum) NuU1 = mup; else NuU1 = mum;
						NuU1_new = std::max(fabs(mup),fabs(mum));
#if defined(ESTIMATOR_REFINE)
						NuU2_old = NuU2;
						mup = 1.0 - EstU2[k];
						mum = -1.0 - EstU2[k];
						np = nm = 0;
						//start from the element next after diagonal position
						Ui = LineIndecesU[k];
						while (Ui != EOL)
						{
							v = EstU2[Ui];
							vp = fabs(v + LineValuesU[Ui] * mup);
							vm = fabs(v + LineValuesU[Ui] * mum);
							v = fabs(v);
							if (vp > std::max(2 * v, 0.5)) np++;
							if (std::max(2 * vp, 0.5) < v) np--;
							if (vm > std::max(2 * v, 0.5)) nm++;
							if (std::max(2 * vm, 0.5) < v) nm--;
							Ui = LineIndecesU[Ui];
						}
						if (np > nm) NuU2 = mup; else NuU2 = mum;
						NuU2_new = std::max(fabs(mup), fabs(mum));
#endif
						testimator += Timer()-tlocal;
					} //if( estimator )
#endif //ESTIMATOR
///////////////////////////////////////////////////////////////////////////////////
//                 condition estimation for L part                               //
///////////////////////////////////////////////////////////////////////////////////
					//estimate condition for L^{-1}
#if defined(ESTIMATOR)
					if( estimator )
					{
						tlocal = Timer();
						NuL1_old = NuL1;
						mup = 1.0 - EstL1[k];
						mum = -1.0 - EstL1[k];
						smup = smum = 0;
						Li = LineIndecesL[k];
						while (Li != EOL)
						{
							v = EstL1[Li];
							vp = fabs(v + LineValuesL[Li] * mup);
							vm = fabs(v + LineValuesL[Li] * mum);
							smup += vp;
							smum += vm;
							Li = LineIndecesL[Li];
						}
						if (smup > smum) NuL1 = mup; else NuL1 = mum;
						NuL1_new = std::max(fabs(mup),fabs(mum));
#if defined(ESTIMATOR_REFINE)
						NuL2_old = NuL2;
						mup = 1.0 - EstL2[k];
						mum = -1.0 - EstL2[k];
						np = nm = 0;
						//start from the element next after diagonal position
						Li = LineIndecesL[k];
						while (Li != EOL)
						{
							v = EstL2[Li];
							vp = fabs(v + LineValuesL[Li] * mup);
							vm = fabs(v + LineValuesL[Li] * mum);
							v = fabs(v);
							if (vp > std::max(2 * v, 0.5)) np++;
							if (std::max(2 * vp, 0.5) < v) np--;
							if (vm > std::max(2 * v, 0.5)) nm++;
							if (std::max(2 * vm, 0.5) < v) nm--;
							Li = LineIndecesL[Li];
						}
						if (np > nm) NuL2 = mup; else NuL2 = mum;
						NuL2_new = std::max(fabs(mup), fabs(mum));
#endif
						testimator += Timer()-tlocal;
					}
#endif //ESTIMATOR
///////////////////////////////////////////////////////////////////////////////////
//         discarding current iteration based on condition numbers               //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
					//discard line if condition number is too large
#if defined(ESTIMATOR_REFINE)
					NuU_tmp = std::max(NuU1_new,NuU2_new);
					NuL_tmp = std::max(NuL1_new,NuL2_new);
#else
					NuU_tmp = NuU1_new;
					NuL_tmp = NuL1_new;
#endif
#if defined(CONDITION_PIVOT)
					if(allow_pivot && !block_pivot && (NuU_tmp > pivot_cond || NuL_tmp > pivot_cond || NuD > pivot_diag) )
					{
						//restore condition number
						NuL1 = NuU1_old;
						NuU1 = NuU1_old;
						EstU1[k] = EstL1[k] = 0;
#if defined(ESTIMATOR_REFINE)
						NuL2 = NuU2_old;
						NuU2 = NuU2_old;
						EstU2[k] = EstL2[k] = 0;
#endif
						NuD = NuD_old;
						max_diag = max_diag_old;
						min_diag = min_diag_old;
						//Cancel this row
						Pivot[k] = true;
						swaps++;
						Li = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
						Ui = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
						U_Address[k].first = U_Address[k].last = Li;
						L_Address[k].first = L_Address[k].last = Li;
						U2_Address[k].first = U2_Address[k].last = Ui;
						L2_Address[k].first = L2_Address[k].last = Ui;
						//Proceed to the next row
						goto exit_point;
					}
					else
#endif
					{
#if defined(ESTIMATOR_REFINE)
						NuL = std::max(NuL1_new,NuL2_new);
						NuU = std::max(NuU1_new,NuU2_new);
#else
						NuL = NuL1_new;
						NuU = NuU1_new;
#endif
						NuL_max = std::max(NuL,NuL_max);
						NuU_max = std::max(NuU,NuU_max);
					}
#endif //ESTIMATOR
					
///////////////////////////////////////////////////////////////////////////////////
//                 reconstruct U-part from linked list                           //
///////////////////////////////////////////////////////////////////////////////////
					//insert line to U part
					U_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
					U2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					
					Unum = Unorm = 0;
					
					Ui = LineIndecesU[k];
					while (Ui != EOL)
					{
						u = LineValuesU[Ui];
						Unorm += u*u;
						Unum++;
						Ui = LineIndecesU[Ui];
					}
					if( Unum ) Unorm = sqrt(Unorm/Unum);
					
					//Unorm = std::min(1.0,Unorm);
					
					Ui = LineIndecesU[k];
					while (Ui != EOL)
					{
						u = fabs(LineValuesU[Ui]);
						if (u*NuU > tau )//*Unorm) // apply dropping rule
						//if (u*NuU*NuU_acc*NuD*NuD_acc > tau) // apply dropping rule
						//if (u*NuU_acc*NuD_acc > tau) // apply dropping rule
						//if( u > tau*Unorm )
							LU_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#if defined(ILUC2)
						else if (u*NuU > tau2 )//*Unorm)
						//else if (u*NuU*NuU_acc*NuD*NuD_acc > tau2)
						//else if (u*NuU_acc*NuD_acc > tau2)
						//else if( u > tau2*Unorm )
							LU2_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#endif
						else ndrops++;
						Ui = LineIndecesU[Ui];
					}
					U_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
					U2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 reconstruct L-part from linked list                           //
///////////////////////////////////////////////////////////////////////////////////
					//insert column to L part
					L_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
					L2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					
					Lnum = Lnorm = 0;
					
					Li = LineIndecesL[k];
					while (Li != EOL)
					{
						u = LineValuesL[Li];
						Lnorm += u*u;
						Lnum++;
						Li = LineIndecesL[Li];
					}
					if( Lnum ) Lnorm = sqrt(Lnorm/Lnum);
					
					//Lnorm = std::min(1.0,Lnorm);
					
					Li = LineIndecesL[k];
					while (Li != EOL)
					{
						u = fabs(LineValuesL[Li]);
						if (u*NuL > tau )//*Lnorm) //apply dropping
						//if (u*NuL*NuL_acc*NuD*NuD_acc > tau) //apply dropping
						//if (u*NuL_acc*NuD_acc > tau) //apply dropping
						//if( u > tau*Lnorm )
							LU_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#if defined(ILUC2)
						else if (u*NuL > tau2 )//*Lnorm)
						//else if (u*NuL*NuL_acc*NuD*NuD_acc > tau2)
						//else if (u*NuL_acc*NuD_acc > tau2)
						//else if( u > tau2*Lnorm )
							LU2_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#endif
						else ndrops++;
						Li = LineIndecesL[Li];
					}
					L_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
					L2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 Update estimator vectors                                      //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
					if( estimator )
					{
						tlocal = Timer();
						{ //U-estimator
							Ui = LineIndecesU[k];
							while (Ui != EOL)
							{
								EstU1[Ui] += LineValuesU[Ui] * NuU1;
								Ui = LineIndecesU[Ui];
							}
							//choose new condition number
							NuU1 = NuU1_new;//std::max(NuU_new,NuU_old);
							//update estimator vector
#if defined(ESTIMATOR_REFINE)
							Ui = LineIndecesU[k];
							while (Ui != EOL)
							{
								EstU2[Ui] += LineValuesU[Ui] * NuU2;
								Ui = LineIndecesU[Ui];
							}
							//choose new condition number
							NuU2 = NuU2_new;//std::max(NuU_new,NuU_old);
#endif
						}
						{ //L-estimator
							Li = LineIndecesL[k];
							while (Li != EOL)
							{
								EstL1[Li] += LineValuesL[Li] * NuL1;
								Li = LineIndecesL[Li];
							}
							NuL1 = NuL1_new;
							//update estimator vector
#if defined(ESTIMATOR_REFINE)
							Li = LineIndecesL[k];
							while (Li != EOL)
							{
								EstL2[Li] += LineValuesL[Li] * NuL2;
								Li = LineIndecesL[Li];
							}
							//choose new condition number
							NuL2 = NuL2_new;//std::max(NuL_new,NuL_old);
#endif
						}
						testimator += Timer() - tlocal;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//     Update diagonal                                                           //
///////////////////////////////////////////////////////////////////////////////////
					LU_Diag[k] = udiag;
///////////////////////////////////////////////////////////////////////////////////
//      Iteration exit point with cleanups and updates                           //
///////////////////////////////////////////////////////////////////////////////////
				exit_point:
///////////////////////////////////////////////////////////////////////////////////
//                 Cleaning data structures                                      //
///////////////////////////////////////////////////////////////////////////////////
					Li = cbeg;
					while (Li != EOL)
					{
						i = LineIndecesL[Li];
						LineValuesL[Li] = 0.0; //clean values after use
						LineIndecesL[Li] = UNDEF; //clean indeces after use
						Li = i;
					}
					Ui = cbeg;
					while (Ui != EOL)
					{
						i = LineIndecesU[Ui];
						LineValuesU[Ui] = 0.0; // clean values after use
						LineIndecesU[Ui] = UNDEF; //clean indeces after use
						Ui = i;
					}
///////////////////////////////////////////////////////////////////////////////////
//                 Insert indexes for transposed U-part traversal                //
///////////////////////////////////////////////////////////////////////////////////
					//Update linked list for factors
					if (U_Address[k].Size() != 0)
					{
						i = LU_Entries[U_Address[k].first].first;
						assert(i > k);
						Ulist[k] = Ubeg[i];
						Ubeg[i] = k;
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         Insert indexes for transposed second-order U-part traversal           //
///////////////////////////////////////////////////////////////////////////////////
					if (U2_Address[k].Size() != 0)
					{
						i = LU2_Entries[U2_Address[k].first].first;
						assert(i > k);
						U2list[k] = U2beg[i];
						U2beg[i] = k;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 Insert indexes for transposed L-part traversal                //
///////////////////////////////////////////////////////////////////////////////////
					//update linked list for factors
					if (L_Address[k].Size() > 0)
					{
						i = LU_Entries[L_Address[k].first].first;
						assert(i > k);
						Llist[k] = Lbeg[i];
						Lbeg[i] = k;
					}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         Insert indexes for transposed second-order L-part traversal           //
///////////////////////////////////////////////////////////////////////////////////
					if (L2_Address[k].Size() != 0)
					{
						i = LU2_Entries[L2_Address[k].first].first;
						assert(i > k);
						L2list[k] = L2beg[i];
						L2beg[i] = k;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//     Shift indexes for transposed U-part traversal for next iteration          //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//     Shift indexes for transposed L-part traversal for next iteration          //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//         Shift indexes for transposed second-order U-part traversal            //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//         Shift indexes for transposed second-order L-part traversal            //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//         Shift indexes for transposed B matrix traversal                       //
///////////////////////////////////////////////////////////////////////////////////
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
#if defined(TRACK_DIAGONAL)
					if( !Pivot[k] )
					{
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with L and U                      //
///////////////////////////////////////////////////////////////////////////////////
						//update diagonal by optained factors
						i = U_Address[k].first;
						j = L_Address[k].first;
						while (i < U_Address[k].last && j < L_Address[k].last)
						{
							Ui = LU_Entries[i].first;
							Li = LU_Entries[j].first;
							if (Ui > Li) j++;
							else if (Ui < Li) i++;
							else
							{
								assert(Ui > k);
								LU_Diag[Ui] -= LU_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
								i++;
								j++;
							}
						}
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with L and second-order U         //
///////////////////////////////////////////////////////////////////////////////////
						i = U_Address[k].first;
						j = L2_Address[k].first;
						while (i < U_Address[k].last && j < L2_Address[k].last)
						{
							Ui = LU_Entries[i].first;
							Li = LU2_Entries[j].first;
							if (Ui > Li) j++;
							else if (Ui < Li) i++;
							else
							{
								assert(Ui > k);
								LU_Diag[Ui] -= LU_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
								i++;
								j++;
							}
						}
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with second-order L and U         //
///////////////////////////////////////////////////////////////////////////////////
						i = U2_Address[k].first;
						j = L_Address[k].first;
						while (i < U2_Address[k].last && j < L_Address[k].last)
						{
							Ui = LU2_Entries[i].first;
							Li = LU_Entries[j].first;
							if (Ui > Li) j++;
							else if (Ui < Li) i++;
							else
							{
								assert(Ui > k);
								LU_Diag[Ui] -= LU2_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
								i++;
								j++;
							}
						}
///////////////////////////////////////////////////////////////////////////////////
// Update values on the whole diagonal with second-order L and second-order U    //
///////////////////////////////////////////////////////////////////////////////////
#if 0
						i = U2_Address[k].first;
						j = L2_Address[k].first;
						while (i < U2_Address[k].last && j < L2_Address[k].last)
						{
							Ui = LU2_Entries[i].first;
							Li = LU2_Entries[j].first;
							if (Ui > Li) j++;
							else if (Ui < Li) i++;
							else
							{
								assert(Ui > k);
								LU_Diag[Ui] -= LU2_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
								i++;
								j++;
							}
						}
#endif
#endif //ILUC2
					}
#endif //TRACK_DIAGONAL
					//number of nonzeros
					nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
					if( verbosity )
					{
						if( k % 100 == 0 )
						{
							if (verbosity == 1)
								printf("%d/%d factor %6.2f%%\t\t\r",cend,moend, 100.0f*(k - cbeg) / (float)(cend - cbeg));
							else 
								printf("%6.2f%% nnz LU %8d LU2 %8d condition L %10f D %10f U %10f swaps %4d drops %4d\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, nzLU2, NuL, NuD, NuU, swaps,ndrops);
							fflush(stdout);
						}
					}
///////////////////////////////////////////////////////////////////////////////////
//                       iteration done                                          //
///////////////////////////////////////////////////////////////////////////////////
				}
				if( verbosity > 1 )
				{
					std::cout << "size " << cend-cbeg;
					std::cout << " total nonzeros in A " << nzA << " (sparsity " << nzA/(double)(cend-cbeg)/(double)(cend-cbeg)*100.0 << "%) in LU " << nzLU << " (fill " << nzLU/(double)nzA0*100.0 << "%)";
					std::cout << " in LU2 " << nzLU2;
					std::cout << " conditions L " << NuL_max << " D " << NuD << " U " << NuU_max << " pivot swaps " << swaps << " drops " << ndrops << "            " << std::endl;
				}
			}

			//restore indexes
			B_Address[cbeg].first = Bstart[cbeg];
			U_Address[cbeg].first = LU_Beg;
			L_Address[cbeg].first = U_Address[cbeg].last;
#if defined(ILUC2)
			U2_Address[cbeg].first = 0;
			L2_Address[cbeg].first = U2_Address[cbeg].last;
#endif
			for (k = cbeg+1; k < cend; ++k)
			{
				B_Address[k].first = Bstart[k];
				U_Address[k].first = L_Address[k - 1].last;
				L_Address[k].first = U_Address[k].last;
#if defined(ILUC2)
				U2_Address[k].first = L2_Address[k - 1].last;
				L2_Address[k].first = U2_Address[k].last;
#endif
			}
///////////////////////////////////////////////////////////////////////////////////
///             FACTORIZATION COMPLETE                                          ///
///////////////////////////////////////////////////////////////////////////////////
			tlfactor = Timer() - tt;
			tfactor += tlfactor;
			if( verbosity > 1 )
				printf("Time %g\n",tlfactor);

			//After we have factored the rescaled matrix we must rescale obtained factors
			/*
			tt = Timer();
			if( rescale_b )
			{
#if defined(REPORT_ILU)
				std::cout << std::endl << " rescaling block B back " << std::endl;
#endif
				for (k = cbeg; k < cend; ++k)
				{
					LU_Diag[k] /= (DL[k] * DR[k]);
					for (INMOST_DATA_ENUM_TYPE r = U_Address[k].first; r < U_Address[k].last; ++r) LU_Entries[r].second *= (DR[k] / DR[LU_Entries[r].first]);
					for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r) LU_Entries[r].second *= (DL[k] / DL[LU_Entries[r].first]);
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) B_Entries[r].second /= (DL[k] * DR[B_Entries[r].first]);
				}
#if defined(ILUC2)
				if( swaps )
				{
					for (k = cbeg; k < cend; ++k)
					{
						for (INMOST_DATA_ENUM_TYPE r = U2_Address[k].first; r < U2_Address[k].last; ++r) LU2_Entries[r].second *= (DR[k] / DR[LU2_Entries[r].first]);
						for (INMOST_DATA_ENUM_TYPE r = L2_Address[k].first; r < L2_Address[k].last; ++r) LU2_Entries[r].second *= (DL[k] / DL[LU2_Entries[r].first]);
					}
				}
#endif
			}
			tlrescale = Timer() - tt;
			trescale += tlrescale;
			 */
///////////////////////////////////////////////////////////////////////////////////
///                RESCALING DONE                                               ///
///////////////////////////////////////////////////////////////////////////////////
			
			//reset the scaling vectors
			for (k = cbeg; k < cend; ++k) DL[k] = DR[k] = 1.0;
///////////////////////////////////////////////////////////////////////////////////
//  Check that we have to enter the next level due to pivoting                   //
///////////////////////////////////////////////////////////////////////////////////
			if( swaps )
			{
				tlocal = Timer();
				cend = wend-swaps;
				tt = Timer();
				if( verbosity > 1 )
				{
					std::cout << "Total swaps: " << swaps << " interval: " << cend << " " << wend << std::endl;
					printf("Reassemble pivots\n");
				}

				level_size.push_back(cend - wbeg);
				i = wbeg;
				//enumerate entries that we keep first
				for(k = wbeg; k < wend; ++k) if( !Pivot[k] )
				{
					localP[k] = i;
					localQ[k] = i;
					++i;
				}
				//enumerate entries that we dropped at the end
				for(k = wbeg; k < wend; ++k) if( Pivot[k] )
				{
					localP[k] = i;
					localQ[k] = i;
					++i;
				}
				
				for (k = wbeg; k < wend; ++k)
				{
					U[localP[k]] = DL0[k];
					V[localQ[k]] = DR0[k];
				}
				for (k = wbeg; k < wend; ++k)
				{
					DL0[k] = U[k];
					DR0[k] = V[k];
				}
				
				//inverse the ordering
				ReorderEF(wbeg, wend, donePQ, localP, localQ);
				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
///////////////////////////////////////////////////////////////////////////////////
//          Reassamble matrix                                                    //
///////////////////////////////////////////////////////////////////////////////////
				//here F is stored by columns and E by rows
				A_Entries.clear();
				E_Address.push_back(new interval<INMOST_DATA_ENUM_TYPE, Interval>(cend, wend));
				F_Address.push_back(new interval<INMOST_DATA_ENUM_TYPE, Interval>(cend, wend));
				//find out the space required for each column of F
				//and setup the first and last pointers
				for(k = cend; k < wend; ++k)
					F_Address.back()->at(k).first = F_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(F_Entries.size());
				for(k = wbeg; k < wend; ++k)
				{
					j = invP[k];
					if( k < cend )
					{
						for (INMOST_DATA_ENUM_TYPE jt = B_Address[j].first; jt < B_Address[j].last; ++jt)
						{
							i = localQ[B_Entries[jt].first];
							if( i >= cend ) //put into F
								F_Address.back()->at(i).last++;
						}
					}
				}
				//Establish the addresses
				INMOST_DATA_ENUM_TYPE offset = 0;
				for(k = cend; k < wend; ++k)
				{
					F_Address.back()->at(k).first += offset;
					F_Address.back()->at(k).last  += offset;
					offset += F_Address.back()->at(k).Size();
				}
				F_Entries.resize(F_Address.back()->at(wend-1).last);
				//Reset addresses to advance current fill position
				for(k = cend; k < wend; ++k)
					F_Address.back()->at(k).last = F_Address.back()->at(k).first;
				//reassamble lower-right corner of B into A
				//put lower-left into E and upper-right into F
				for(k = wbeg; k < wend; ++k)
				{
					j = invP[k];
					if( k >= cend )
					{
						E_Address.back()->at(k).first = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries.size());
						A_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(A_Entries.size());
						for (INMOST_DATA_ENUM_TYPE jt = B_Address[j].first; jt < B_Address[j].last; ++jt)
						{
							i = localQ[B_Entries[jt].first];
							u = B_Entries[jt].second;
							if( i >= cend ) //put into A
								A_Entries.push_back(Sparse::Row::make_entry(i, u));
							else //put into E
								E_Entries.push_back(Sparse::Row::make_entry(i, u));
						}
						A_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(A_Entries.size());
						E_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries.size());
						std::sort(A_Entries.begin() + A_Address[k].first, A_Entries.end());
						std::sort(E_Entries.begin() + E_Address.back()->at(k).first, E_Entries.end());
					}
					else
					{
						for (INMOST_DATA_ENUM_TYPE jt = B_Address[j].first; jt < B_Address[j].last; ++jt)
						{
							i = localQ[B_Entries[jt].first];
							u = B_Entries[jt].second;
							if( i >= cend ) //put into F
								F_Entries[F_Address.back()->at(i).last++] = Sparse::Row::make_entry(k, u);
						}
					}
				}
				//restore order in F
				/*
				for(k = cend; k < wend; ++k)
				{
					std::sort(F_Entries.begin()+F_Address.back()->at(k).first, F_Entries.begin()+F_Address.back()->at(k).last);
				}
				 */
				A_Entries.push_back(Sparse::Row::make_entry(-1,0.0));
				
				//DumpMatrix(A_Address,A_Entries,cend,wend,"sA.mtx");
				//DumpMatrix(*E_Address.back(),E_Entries,wbeg,wend,"E.mtx");
				//DumpMatrix(*F_Address.back(),F_Entries,wbeg,wend,"F.mtx");
				//reorder LU according to pivoting
///////////////////////////////////////////////////////////////////////////////////
//  Reassamble L and U                                                           //
///////////////////////////////////////////////////////////////////////////////////
				//DumpMatrix(L_Address,LU_Entries,wbeg,wend,"preL.mtx");
				//DumpMatrix(U_Address,LU_Entries,wbeg,wend,"preU.mtx");
				//DumpMatrix(L2_Address,LU2_Entries,wbeg,wend,"preL2.mtx");
				//DumpMatrix(U2_Address,LU2_Entries,wbeg,wend,"preU2.mtx");
				/*
				{
					std::fstream fout("preDiag.mtx",std::ios::out);
					fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
					fout << std::scientific;
					fout << wend-wbeg << " " << wend-wbeg << " " << wend-wbeg << std::endl;
					for(k = wbeg; k < wend; ++k)
						fout << k+1 << " " << k+1 << " " << LU_Diag[k] << std::endl;
					fout.close();
				}
				*/
				//first change positions
				for(k = wbeg; k < wend; ++k)
				{
					j = invP[k];
					//if( j < cend )
					{
						INMOST_DATA_ENUM_TYPE cnt;
						// L-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = L_Address[j].first; jt < L_Address[j].last; ++jt)
						{
							i = localQ[LU_Entries[jt].first];
							if( i < cend )
								LU_Entries[jt].first = i;
							else
							{
								LU_Entries[jt].first = ENUMUNDEF;
								LU_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(LU_Entries.begin()+L_Address[j].first,LU_Entries.begin()+L_Address[j].last);
						//this will remove references to undefined
						L_Address[j].last -= cnt;
						// U-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = U_Address[j].first; jt < U_Address[j].last; ++jt)
						{
							Ui = LU_Entries[jt].first;
							i = localQ[Ui];
							if( i < cend )
								LU_Entries[jt].first = i;
							else
							{
								LU_Entries[jt].first = ENUMUNDEF;
								LU_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(LU_Entries.begin()+U_Address[j].first,LU_Entries.begin()+U_Address[j].last);
						//this will remove references to undefined
						U_Address[j].last -= cnt;
					}
				}
				//now the same for second-order LU
				for(k = wbeg; k < wend; ++k)
				{
					j = invP[k];
					//if( j < cend )
					{
						INMOST_DATA_ENUM_TYPE cnt;
						//L-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = L2_Address[j].first; jt < L2_Address[j].last; ++jt)
						{
							i = localQ[LU2_Entries[jt].first];
							if( i < cend )
								LU2_Entries[jt].first = i;
							else
							{
								LU2_Entries[jt].first = ENUMUNDEF;
								LU2_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(LU2_Entries.begin()+L2_Address[j].first,LU2_Entries.begin()+L2_Address[j].last);
						//this will remove references to undefined
						L2_Address[j].last -= cnt;
						//U-part
						cnt = 0;
						for (INMOST_DATA_ENUM_TYPE jt = U2_Address[j].first; jt < U2_Address[j].last; ++jt)
						{
							i = localQ[LU2_Entries[jt].first];
							if( i < cend )
								LU2_Entries[jt].first = i;
							else
							{
								LU2_Entries[jt].first = ENUMUNDEF;
								LU2_Entries[jt].second = 0.0;
								cnt++;
							}
						}
						//this will move all undefined to the end
						std::sort(LU2_Entries.begin()+U2_Address[j].first,LU2_Entries.begin()+U2_Address[j].last);
						//this will remove references to undefined
						U2_Address[j].last -= cnt;
					}
				}
				//first establish intervals in arrays
				{
					interval<INMOST_DATA_ENUM_TYPE, Interval> tmpL(wbeg,cend), tmpU(wbeg,cend), tmpL2(wbeg,cend), tmpU2(wbeg,cend);
					interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> tmpD(wbeg,cend);
					for(k = wbeg; k < cend; ++k)
					{
						j = invP[k];
						tmpL[k] = L_Address[j];
						tmpU[k] = U_Address[j];
						tmpL2[k] = L2_Address[j];
						tmpU2[k] = U2_Address[j];
						tmpD[k] = LU_Diag[j];
					}
					for(k = wbeg; k < cend; ++k)
					{
						L_Address[k] = tmpL[k];
						U_Address[k] = tmpU[k];
						L2_Address[k] = tmpL2[k];
						U2_Address[k] = tmpU2[k];
						LU_Diag[k] = tmpD[k];
					}
				}
				//DumpMatrix(L_Address,LU_Entries,wbeg,cend,"L.mtx");
				//DumpMatrix(U_Address,LU_Entries,wbeg,cend,"U.mtx");
				//DumpMatrix(L2_Address,LU2_Entries,wbeg,cend,"L2.mtx");
				//DumpMatrix(U2_Address,LU2_Entries,wbeg,cend,"U2.mtx");
				/*
				{
					std::fstream fout("Diag.mtx",std::ios::out);
					fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
					fout << std::scientific;
					fout << cend-wbeg << " " << cend-wbeg << " " << cend-wbeg << std::endl;
					for(k = wbeg; k < cend; ++k)
						fout << k+1 << " " << k+1 << " " << LU_Diag[k] << std::endl;
					fout.close();
				}
				 */
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);
				tt = Timer()-tt;
				if( verbosity > 1 )
					printf("Time %g\n",tt);

				int ndrops_lf = 0, ndrops_eu = 0, ndrops_s = 0;
///////////////////////////////////////////////////////////////////////////////////
//  Compute Schur                                                                //
///////////////////////////////////////////////////////////////////////////////////
				//result of multilevel preconditioner
				//        |LDU  F |
				//  A  =  |       |
				//        | E   A |
				//
				// For the next step of multilevel factorization we have to compute
				// S = C - (E U^{-1}) D^{-1} (L^{-1} F)
				//
				//first precompute LF block
				tt = Timer();
				if( verbosity > 1 )
					printf("Assemble LF\n");

				
				for (k = cend; k < wend; ++k)
				{
					Li = cbeg;
					for(INMOST_DATA_ENUM_TYPE jt = F_Address.back()->at(k).first; jt < F_Address.back()->at(k).last; ++jt)
					{
						u = F_Entries[jt].second;
						j = F_Entries[jt].first;
						LineValuesU[j] = u;
						Li = LineIndecesU[Li] = j+1;
					}
					LineIndecesU[Li] = EOL;
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with L part                                             //
///////////////////////////////////////////////////////////////////////////////////
					//compute ~F_i = L^{-1} F_i
					Li = LineIndecesU[cbeg];
					Sbeg = cbeg;
					while (Li != EOL)
					{
						curr = Li;
						if( 1 + LineValuesU[Li - 1] != 1 )
						{
							for (INMOST_DATA_ENUM_TYPE ru = L_Address[Li - 1].first; ru < L_Address[Li - 1].last; ++ru)
							{
								u = LineValuesU[Li - 1] * LU_Entries[ru].second;
								j = LU_Entries[ru].first;
								if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
									LineValuesU[j] -= u;
								else if( 1+u != 1 )
								{
									
									next = curr;
									while (next < j + 1)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j + 1);
									assert(j + 1 < next);
									LineIndecesU[curr] = j + 1;
									LineIndecesU[j + 1] = next;
									LineValuesU[j] = -u;
									
									//LineValuesU[j] = -u;
									//LineIndecesU[j + 1] = Sbeg;
									//Sbeg = j+1;
								}
							}
#if defined(ILUC2) && defined(ILUC2_SCHUR)
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with second-order L part                                //
///////////////////////////////////////////////////////////////////////////////////
							curr = Li;
							for (INMOST_DATA_ENUM_TYPE ru = L2_Address[Li - 1].first; ru < L2_Address[Li - 1].last; ++ru)
							{
								u = LineValuesU[Li - 1] * LU2_Entries[ru].second;
								j = LU2_Entries[ru].first;
								if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
									LineValuesU[j] -= u;
								else if( 1+u != 1 )
								{
									
									next = curr;
									while (next < j + 1)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j + 1);
									assert(j + 1 < next);
									LineIndecesU[curr] = j + 1;
									LineIndecesU[j + 1] = next;
									LineValuesU[j] = -u;
									
									//LineValuesU[j] = -u;
									//LineIndecesU[j + 1] = Sbeg;
									//Sbeg = j+1;
								}
							}
#endif
						}
///////////////////////////////////////////////////////////////////////////////////
//             solve iteration done                                              //
///////////////////////////////////////////////////////////////////////////////////						
						Li = LineIndecesU[Li];
					}
///////////////////////////////////////////////////////////////////////////////////
//                 Rescale by diagonal                                           //
///////////////////////////////////////////////////////////////////////////////////
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						LineValuesU[Li - 1] /= LU_Diag[Li-1]; //Li - 1 or k?
						Li = LineIndecesU[Li];
					}
///////////////////////////////////////////////////////////////////////////////////
//                 Compute norm for dropping                                     //
///////////////////////////////////////////////////////////////////////////////////
#if defined(SCHUR_DROPPING_LF)
					LFnorm = LFnum = 0;
					Li = LineIndecesU[cbeg];
					LFmax = 0;
					LFmin = 1.0e+54;
					while (Li != EOL)
					{
						u = fabs(LineValuesU[Li - 1]);
						LFnorm += u*u;
						if( u > LFmax ) LFmax = u;
						if( u < LFmin ) LFmin = u;
						LFnum++;
						Li = LineIndecesU[Li];
					}
					if( LFnum ) LFnorm = sqrt(LFnorm/LFnum);
					//LFnorm = std::min(1.0,LFnorm);
					//LFtau = LFmin + (LFmax - LFmin)*std::min(tau*tau,tau2) / NuL_max;
					LFtau = LFmin + std::min(tau*tau,tau2)*LFnorm / NuL_max;
					//LFtau = LFmin + (LFmax - LFmin)*pow(LFnorm/LFmax,4);
					//LFtau = LFmin + (LFmax - LFmin)*pow(LFnorm/LFmax,4);
#endif
///////////////////////////////////////////////////////////////////////////////////
//                 Assemble column into matrix                                   //
///////////////////////////////////////////////////////////////////////////////////
					LFdrop = 0;
					Li = LineIndecesU[cbeg];
					LF_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
					while (Li != EOL)
					{
						u = LineValuesU[Li - 1];
						assert(Li - 1 >= wbeg && Li - 1 < cend);
#if defined(SCHUR_DROPPING_LF)
						if( fabs(u) >= LFtau )
#else
						if( 1+u != 1 )
#endif
							LF_Entries.push_back(Sparse::Row::make_entry(Li - 1, u));
						else 
						{
							LFdrop += u;
							ndrops_lf++;
						}
						Li = LineIndecesU[Li];
					}
					LF_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
#if defined(SCHUR_DROPPING_LF)					
					if( LFdrop )
					{
						LFdrop /= (double) LF_Address[k].Size();
						for(j = LF_Address[k].first; j < LF_Address[k].last; ++j)
							LF_Entries[j].second += LFdrop;
					}
#endif
					//assert(std::is_sorted(LF_Entries.begin()+LF_Address[k].first,LF_Entries.end()));
///////////////////////////////////////////////////////////////////////////////////
//              clean linked list                                                //
///////////////////////////////////////////////////////////////////////////////////
					//Clean after use
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						j = LineIndecesU[Li];
						LineIndecesU[Li] = UNDEF;
						LineValuesU[Li - 1] = 0.0;
						Li = j;
					}
					LineIndecesU[cbeg] = UNDEF;
					if( verbosity )
					{
						if (k % 100 == 0)
						{
							printf("LF %6.2f%% nnz %lu drops %d\t\t\r", 100.f*(k - cend+1) / (1.f*(wend - cend)),LF_Entries.size(),ndrops_lf);
							fflush(stdout);
						}
					}

///////////////////////////////////////////////////////////////////////////////////
//             iteration done!                                                   //
///////////////////////////////////////////////////////////////////////////////////
				}
				tt = Timer()-tt;
				if( verbosity > 1 )
				{
					printf("LF nnz %lu drops %d         \t\t\n", LF_Entries.size(),ndrops_lf);
					printf("Time %g\n",tt);
				}

				//DumpMatrix(LF_Address,LF_Entries,wbeg,wend,"LF.mtx");
///////////////////////////////////////////////////////////////////////////////////
//         prepearing LF block for transposed traversal                          //
///////////////////////////////////////////////////////////////////////////////////
				tt = Timer();
				if( verbosity > 1 )
					printf("Transpose LF\n");
				//std::fill(Fbeg.begin() + wbeg - mobeg, Fbeg.begin() + cend - mobeg, EOL);
				for(k = wbeg; k < cend; ++k) Fbeg[k] = EOL;
				for(k = wend; k > cend; --k)
				{
					if (LF_Address[k-1].Size() > 0)
					{
						i = LF_Entries[LF_Address[k-1].first].first;
						Flist[k-1] = Fbeg[i];
						Fbeg[i] = k-1;
					}
					else Flist[k-1] = EOL;
				}
				//transpone LF matrix
///////////////////////////////////////////////////////////////////////////////////
//          transpose LF block                                                   //
///////////////////////////////////////////////////////////////////////////////////
				LFt_Entries.resize(LF_Entries.size());
				j = 0;
				for (k = wbeg; k < cend; k++)
				{
///////////////////////////////////////////////////////////////////////////////////
//         assemble line of transposed LF                                        //
///////////////////////////////////////////////////////////////////////////////////
					LFt_Address[k].first = j;
					Li = Fbeg[k];
					//Fnorms[k] = 0;
					while (Li != EOL)
					{
						LFt_Entries[j].first = Li;
						LFt_Entries[j].second = LF_Entries[LF_Address[Li].first].second;
						j++;
						Li = Flist[Li];
					}
					LFt_Address[k].last = j;
///////////////////////////////////////////////////////////////////////////////////
//          shift indices for transposed traversal                               //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//              iteration done                                                   //
///////////////////////////////////////////////////////////////////////////////////
				}
				LF_Entries.clear();
				tt = Timer()-tt;
				if( verbosity > 1 )
					printf("Time %g\n",tt);

				//DumpMatrix(LFt_Address,LFt_Entries,wbeg,wend,"LFt.mtx");
///////////////////////////////////////////////////////////////////////////////////
//             EU and Schur                                                      //
///////////////////////////////////////////////////////////////////////////////////
				tt = Timer();
				if( verbosity > 1 )
					printf("Assemble EU\n");

				for(k = cend; k < wend; ++k)
				{
///////////////////////////////////////////////////////////////////////////////////
//          no values at E block - unpack C block into Schur                     //
///////////////////////////////////////////////////////////////////////////////////
					Li = cbeg;
					for (j = E_Address.back()->at(k).first; j < E_Address.back()->at(k).last; ++j) //iterate over values of k-th row of E
					{
						u = E_Entries[j].second;
						i = E_Entries[j].first;
						LineValuesU[i] = E_Entries[j].second;
						Li = LineIndecesU[Li] = i + 1;
					}
					LineIndecesU[Li] = EOL;
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with U part                                             //
///////////////////////////////////////////////////////////////////////////////////
					// compute ~E_i = E_i U^{-1}
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						curr = Li;
						if( 1 + LineValuesU[Li - 1] != 1 )
						{
							for (INMOST_DATA_ENUM_TYPE ru = U_Address[Li - 1].first; ru < U_Address[Li - 1].last; ++ru)
							{
								u = LineValuesU[Li - 1] * LU_Entries[ru].second;
								j = LU_Entries[ru].first;
								if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
									LineValuesU[j] -= u;
								else if ( 1+u != 1 )
								{
									next = curr;
									while (next < j + 1)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j + 1);
									assert(j + 1 < next);
									LineIndecesU[curr] = j + 1;
									LineIndecesU[j + 1] = next;
									LineValuesU[j] = -u;
								}
							}
#if defined(ILUC2) && defined(ILUC2_SCHUR)
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with second-order U part                                //
///////////////////////////////////////////////////////////////////////////////////
							curr = Li;
							for (INMOST_DATA_ENUM_TYPE ru = U2_Address[Li - 1].first; ru < U2_Address[Li - 1].last; ++ru)
							{
								u = LineValuesU[Li - 1] * LU2_Entries[ru].second;
								j = LU2_Entries[ru].first;
								if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
									LineValuesU[j] -= u;
								else if ( 1+u != 1 )
								{
									next = curr;
									while (next < j + 1)
									{
										curr = next;
										next = LineIndecesU[curr];
									}
									assert(curr < j + 1);
									assert(j + 1 < next);
									LineIndecesU[curr] = j + 1;
									LineIndecesU[j + 1] = next;
									LineValuesU[j] = -u;
								}
							}
#endif
						}
						Li = LineIndecesU[Li];
///////////////////////////////////////////////////////////////////////////////////
//         solve iteration done                                                  //
///////////////////////////////////////////////////////////////////////////////////
					}
///////////////////////////////////////////////////////////////////////////////////
//         rescale by diagonal                                                   //
///////////////////////////////////////////////////////////////////////////////////
					//rescale vector by diagonal *E = ~E_i D^{-1} = E_i U^{-1} D^{-1}
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						LineValuesU[Li - 1] /= LU_Diag[Li - 1];
						Li = LineIndecesU[Li];
					}
///////////////////////////////////////////////////////////////////////////////////
//       compute norm                                                            //
///////////////////////////////////////////////////////////////////////////////////
#if defined(SCHUR_DROPPING_EU)
					EUnorm = EUnum = 0;
					EUmax = 0;
					EUmin = 1.0e+54;
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						u = fabs(LineValuesU[Li - 1]);
						EUnorm += u*u;
						if( u > EUmax ) EUmax = u;
						if( u < EUmin ) EUmin = u;
						EUnum++;
						Li = LineIndecesU[Li];
					}
					if( EUnum ) EUnorm = sqrt(EUnorm/EUnum);
					EUtau = EUmin + std::min(tau*tau,tau2)*EUnorm / NuU_max;
					//EUnorm = std::min(1.0,EUnorm);
					//EUtau = EUmin + (EUmax - EUmin)*pow(EUnorm/EUmax,4);
					//EUtau = EUmin + (EUmax - EUmin)*pow(EUnorm/EUmax,4);
#endif
///////////////////////////////////////////////////////////////////////////////////
//drop values that do not satisfy tolerances from linked list of line of EU block//
///////////////////////////////////////////////////////////////////////////////////
					//drop values
					EUdrop = 0;
					//Assemble column into matrix
					Li = LineIndecesU[cbeg];
					EU_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(EU_Entries.size());
					while (Li != EOL)
					{
						u = LineValuesU[Li - 1];
						assert(Li - 1 >= wbeg && Li - 1 < cend);
#if defined(SCHUR_DROPPING_EU)
						if( fabs(u) >= EUtau )
#else
						if( 1 + u != u )
#endif
							EU_Entries.push_back(Sparse::Row::make_entry(Li - 1, u));
						else
						{
							EUdrop += u;
							ndrops_eu++;
						}
						Li = LineIndecesU[Li];
					}
					EU_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(EU_Entries.size());
#if defined(SCHUR_DROPPING_EU)					
					if( EUdrop )
					{
						EUdrop /= (double)EU_Address[k].Size();
						for(INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
							EU_Entries[r].second += EUdrop;
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//         clean up linked list                                                  //
///////////////////////////////////////////////////////////////////////////////////
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						j = LineIndecesU[Li];
						LineIndecesU[Li] = UNDEF;
						LineValuesU[Li - 1] = 0.0;
						Li = j;
					}
					LineIndecesU[cbeg] = UNDEF;
					if( verbosity )
					{
						if (k % 100 == 0)
						{
							printf("EU %6.2f%% nnz %lu drops %d\t\t\r", 100.f*(k - cend+1) / (1.f*(wend - cend)),EU_Entries.size(),ndrops_eu);
							fflush(stdout);
						}
					}
				}
				tt = Timer()-tt;
				if(verbosity > 1 )
				{
					printf("EU nnz %lu drops %d          \t\t\n",EU_Entries.size(),ndrops_eu);
					printf("Time %g\n",tt);
				}

#if defined(SCHUR_DROPPING_S)
///////////////////////////////////////////////////////////////////////////////////
//         Compute column-norms of schur                                         //
///////////////////////////////////////////////////////////////////////////////////
				tt = Timer();
				if( verbosity > 1 )
					printf("Compute Schur column norm\n");

				for(k = cend; k < wend; ++k) Scolnorm[k] = Scolnum[k] = 0.0;
				for(k = cend; k < wend; ++k)
				{
					Sbeg = EOL;
					for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
					{
						Scolnorm[A_Entries[r].first] += A_Entries[r].second*A_Entries[r].second;
						LineIndecesU[A_Entries[r].first] = Sbeg;
						Sbeg = A_Entries[r].first;
					}
					for(INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
					{
						Li = EU_Entries[r].first;
						l = EU_Entries[r].second;
						for (INMOST_DATA_ENUM_TYPE r = LFt_Address[Li].first; r < LFt_Address[Li].last; ++r)
						{
							j = LFt_Entries[r].first;
							u = LFt_Entries[r].second;
							v = -l*u*LU_Diag[Li];
							Scolnorm[j] += v*v;
							if( LineIndecesU[j] == UNDEF )
							{
								LineIndecesU[j] = Sbeg;
								Sbeg = j;
							}
						}
					}
					Ui = Sbeg;
					while(Ui != EOL)
					{
						Scolnum[Ui]++;
						Li = Ui;
						Ui = LineIndecesU[Ui];
						LineIndecesU[Li] = UNDEF;
					}
					if( verbosity )
					{
						if (k % 100 == 0)
						{
							printf("Schur column norm %6.2f%%\t\t\r", 100.f*(k - cend+1) / (1.f*(wend - cend)));
							fflush(stdout);
						}
					}
				}
				for(k = cend; k < wend; ++k) Scolnorm[k] = Scolnum[k] ? sqrt(Scolnorm[k]/(double)Scolnum[k]) : 1.0;
				tt = Timer()-tt;
				if( verbosity > 1 )
					printf("\nTime %g\n",tt);
#endif //SCHUR_DROPPING_S
///////////////////////////////////////////////////////////////////////////////////
//         Construction of Schur complement                                      //
///////////////////////////////////////////////////////////////////////////////////
				tt = Timer();
				if( verbosity > 1 )
					printf("Construct Schur\n");

				for(k = cend; k < wend; ++k)
				{
///////////////////////////////////////////////////////////////////////////////////
//         unpack line of A matrix to temporary linked list                      //
///////////////////////////////////////////////////////////////////////////////////
					Sbeg = EOL;
					for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
					{
						LineValuesL[A_Entries[r].first] = A_Entries[r].second;
						LineIndecesL[A_Entries[r].first] = Sbeg;
						Sbeg = A_Entries[r].first;
					}
					//multiply current row of EU by LF and add to list
///////////////////////////////////////////////////////////////////////////////////
//        multiply EU and LF blocks and add to temporary linked list             //
///////////////////////////////////////////////////////////////////////////////////
					for(INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
					{
						//iterate over corresponding row of LF, add multiplication to list
						Li = EU_Entries[r].first;
						l = EU_Entries[r].second;
						for (INMOST_DATA_ENUM_TYPE r = LFt_Address[Li].first; r < LFt_Address[Li].last; ++r)
						{
							j = LFt_Entries[r].first;
							u = LFt_Entries[r].second;
							v = -l*u*LU_Diag[Li];
							if (LineIndecesL[j] != UNDEF)
								LineValuesL[j] += v;
							else if( 1+v != 1 )
							{
								LineValuesL[j] = v;
								LineIndecesL[j] = Sbeg;
								Sbeg = j;
							}
						}
						//Li = LineIndecesU[Li];
					}
///////////////////////////////////////////////////////////////////////////////////
//         sort contents of linked list                                          //
///////////////////////////////////////////////////////////////////////////////////
					indicesS.clear();
					Ui = Sbeg;
					while(Ui != EOL)
					{
						indicesS.push_back(Ui);
						Ui = LineIndecesL[Ui];
					}
					if( !indicesS.empty() )
					{
						std::sort(indicesS.begin(),indicesS.end());
						
						Sbeg = indicesS[0];
						for(size_t qt = 1; qt < indicesS.size(); ++qt)
							LineIndecesL[indicesS[qt-1]] = indicesS[qt];
						LineIndecesL[indicesS.back()] = EOL;
					}
///////////////////////////////////////////////////////////////////////////////////
//         prepare row norm for dropping                                         //
///////////////////////////////////////////////////////////////////////////////////
#if defined(SCHUR_DROPPING_S)
					Snorm = Snum = 0;
					Smax = 0;
					Smin = 1.0e+54;
					
					Ui = Sbeg;
					while (Ui != EOL)
					//for(size_t qt = 0; qt < indicesS.size(); ++qt)
					{
						//Ui = indicesS[qt];
						u = fabs(LineValuesL[Ui]);
						if( u > Smax ) Smax = u;
						if( u < Smin ) Smin = u;
						Snorm += u*u;
						Snum++;
						Ui = LineIndecesL[Ui];
					}
					if( Snum ) Snorm = sqrt(Snorm/Snum);
					//Snorm = std::min(1.0,Snorm);
					
					//tauS = std::min(tau*tau,tau2)*Smax;
					Stau = std::min(tau*tau,tau2)/std::max(NuU_max,NuL_max);
					//Stau = std::min(tau*tau,tau2);
					//Stau = Smin + (Smax-Smin)*std::min(tau*tau,tau2);
					//tauS = Smin + (Smax-Smin)*std::min(std::min(tau2,tau*tau),pow(Snorm/Smax,4));
					//exp(log((fabs(u)-Smin)/(Smax-Smin)*127+1)) > exp(0.5*log(128))
					//fabs(u) > exp(0.5)/127*(Smax-Smin) + Smax
					//Stau = Smin + (Smax-Smin)*pow(Snorm/Smax,4);
#endif
///////////////////////////////////////////////////////////////////////////////////
//         put calculated row to Schur complement                                //
///////////////////////////////////////////////////////////////////////////////////
					Sdrop = 0;
					Ui = Sbeg;
					S_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
					while (Ui != EOL)
					//for(size_t qt = 0; qt < indicesS.size(); ++qt)
					{
						//Ui = indicesS[qt];
						u = LineValuesL[Ui];
#if defined(SCHUR_DROPPING_S)
						//if( fabs(Smax-Smin) < 1.0e-5 || log((fabs(u)-Smin)/(Smax-Smin)*127+1)/log(128) > tau )
						if( fabs(u) >= Stau*std::min(Snorm,Scolnorm[Ui]) )
#else
						if( 1+u != 1 )
#endif
							S_Entries.push_back(Sparse::Row::make_entry(Ui,u));
						else
						{
							Sdrop += u;
							ndrops_s++;
						}
						Ui = LineIndecesL[Ui];
					}
					S_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
///////////////////////////////////////////////////////////////////////////////////
//         Substract dropped elements from maximal element to keep sum of row    //
//    (TODO: may collect dropping per row/column and preserve row-column sum)    //
///////////////////////////////////////////////////////////////////////////////////
#if defined(SCHUR_DROPPING_S)
					if( Sdrop )
					{
						INMOST_DATA_ENUM_TYPE rmax = S_Address[k].first, r;
						Sdrop /= (double)S_Address[k].Size();
						for(r = S_Address[k].first; r < S_Address[k].last; ++r)
							S_Entries[r].second += Sdrop;
					}
#endif
					//assert(std::is_sorted(S_Entries.begin()+S_Addres[k].first,S_Entries.end()));
///////////////////////////////////////////////////////////////////////////////////
//         clean up temporary linked list                                        //
///////////////////////////////////////////////////////////////////////////////////
					//clean after use
					Ui = Sbeg;
					while (Ui != EOL)
					{
						j = LineIndecesL[Ui];
						LineIndecesL[Ui] = UNDEF;
						LineValuesL[Ui] = 0.0;
						Ui = j;
					}
///////////////////////////////////////////////////////////////////////////////////
//         clean up linked list                                                 //
///////////////////////////////////////////////////////////////////////////////////
					/*
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						j = LineIndecesU[Li];
						LineIndecesU[Li] = UNDEF;
						LineValuesU[Li - 1] = 0.0;
						Li = j;
					}
					LineIndecesU[cbeg] = UNDEF;
					*/
					if( verbosity )
					{
						if (k % 100 == 0)
						{
							printf("Schur %6.2f%% nnz %lu drop S %d\t\t\r",  100.f*(k - cend+1) / (1.f*(wend - cend)),S_Entries.size(),ndrops_s);
							fflush(stdout);
						}
					}
///////////////////////////////////////////////////////////////////////////////////
//         Schur complement row done                                             //
///////////////////////////////////////////////////////////////////////////////////
				}
				tt = Timer()-tt;
				if( verbosity > 1 )
				{
					printf("Schur nnz %lu drop S %d          \t\t\n", S_Entries.size(),ndrops_s);
					printf("Time %g\n",tt);
				}

				//DumpMatrix(EU_Address,EU_Entries,wbeg,wend,"EU.mtx");
				LFt_Entries.clear();
				EU_Entries.clear();
				
				A_Entries.swap(S_Entries);
				A_Address.swap(S_Address);
				
				//DumpMatrix(S_Address,S_Entries,cend,wend,"S.mtx");
				
				//std::cout << "Press enter" << std::endl;
				//scanf("%*c");
				
///////////////////////////////////////////////////////////////////////////////////
//         Cleanup arrays for the next iteration                                 //
///////////////////////////////////////////////////////////////////////////////////
				tt = Timer();
				if( verbosity > 1 )
					printf("Cleanup\n");

				for(k = wbeg; k < wend; ++k)
				{
					localP[k] = localQ[k] = ENUMUNDEF;
					Lbeg[k] = Ubeg[k] = L2beg[k] = U2beg[k] = Bbeg[k] = EOL;
					EstL1[k] = EstU1[k] = EstL2[k] = EstU2[k] = 0.0;
					Dist[k] = U[k] = V[k] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					Pivot[k] = false;
				}
				//std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
				//std::fill(localQ.begin() + (wbeg - mobeg), localQ.begin() + (wend - mobeg), ENUMUNDEF);
				//std::fill(Lbeg.begin() + (wbeg-mobeg), Lbeg.begin() + (wend-mobeg),EOL);
				//std::fill(Ubeg.begin() + (wbeg-mobeg), Ubeg.begin() + (wend-mobeg),EOL);
				//std::fill(L2beg.begin() + (wbeg-mobeg), L2beg.begin() + (wend-mobeg),EOL);
				//std::fill(U2beg.begin() + (wbeg-mobeg), U2beg.begin() + (wend-mobeg),EOL);
				//std::fill(Bbeg.begin() + (wbeg-mobeg), Bbeg.begin() + (wend-mobeg),EOL);
				//std::fill(EstL1.begin() + (wbeg-mobeg), EstL1.begin() + (wend-mobeg),0.0);
				//std::fill(EstU1.begin() + (wbeg-mobeg), EstU1.begin() + (wend-mobeg),0.0);
				//std::fill(EstL2.begin() + (wbeg-mobeg), EstL2.begin() + (wend-mobeg),0.0);
				//std::fill(EstU2.begin() + (wbeg-mobeg), EstU2.begin() + (wend-mobeg),0.0);
				//std::fill(Dist.begin() + (wbeg-mobeg), Dist.begin() + (wend-mobeg),std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(U.begin() + (wbeg-mobeg), U.begin() + (wend-mobeg),std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(V.begin() + (wbeg-mobeg), V.begin() + (wend-mobeg),std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(Pivot.begin() + (wbeg-mobeg), Pivot.begin() + (wend-mobeg),false);
				for(k = cend; k < wend; ++k)
				{
					U_Address[k].first = U_Address[k].last = ENUMUNDEF;
					L_Address[k].first = L_Address[k].last = ENUMUNDEF;
					U2_Address[k].first = U2_Address[k].last = ENUMUNDEF;
					L2_Address[k].first = L2_Address[k].last = ENUMUNDEF;
				}
				LU2_Entries.clear();
				//no longer needed
				S_Entries.clear();
				LF_Entries.clear();
				LFt_Entries.clear();
				for(k = wbeg; k < wend; ++k)
				{
					S_Address[k].first = S_Address[k].last = ENUMUNDEF;
					LF_Address[k].first = LF_Address[k].last = ENUMUNDEF;
					LFt_Address[k].first = LFt_Address[k].last = ENUMUNDEF;
					//EU_Address[k].first = EU_Address[k].last = ENUMUNDEF;
				}
				
				if( /*wend-cend < 16 ||*/ A_Entries.empty() )
					block_pivot = true;
				//rescale_b = false;
				//run_mpt = false;
				tt = Timer()-tt;
				tlschur = Timer()-tlocal;
				tschur += tlschur;
				if( verbosity > 1 )
				{
					printf("Time %g\n",tt);
					printf("Total Schur %g\n",tlschur);
				}

				wbeg = cend; //there is more work to do
			}
			else
			{
				level_size.push_back(wend - wbeg);
				wbeg = wend;
			}
			totswaps += swaps;
			swaps = 0;
		}
		
		tt = Timer();
		
		if( rescale_b )
		{
			if( verbosity > 1 )
				std::cout << std::endl << " rescaling block B back " << std::endl;

			for (k = mobeg; k < moend; ++k)
			{
				LU_Diag[k] /= (DL0[k] * DR0[k]);
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k].first; r < U_Address[k].last; ++r) LU_Entries[r].second *= (DR0[k] / DR0[LU_Entries[r].first]);
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r) LU_Entries[r].second *= (DL0[k] / DL0[LU_Entries[r].first]);
				//for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) B_Entries[r].second /= (DL0[k] * DR0[B_Entries[r].first]);
			}
			
			//rescale EF
			if( !level_size.empty() )
			{
				INMOST_DATA_ENUM_TYPE first = mobeg, last;
				INMOST_DATA_ENUM_TYPE kbeg, kend, r;
				for(size_t level = 0; level < level_size.size(); ++level)
				{
					last = first + level_size[level];
					kbeg = last;//std::max(cbeg,last);
					kend = moend;//std::min(cend,moend);
					//std::cout << "Rescale EF back level " << level << " ";
					//std::cout << "size " << level_size[level] << " ";
					//std::cout << "interval [" << first << "," << last << "] ";
					//std::cout << "selected [" << kbeg << ":" << kend << "]" << std::endl;
					for (k = kbeg; k < kend; ++k) 
					{
						if( F_Address[level]->at(k).Size() )
						{
							for (r = F_Address[level]->at(k).first; r < F_Address[level]->at(k).last; ++r)
								F_Entries[r].second /= DR0[k]*DL0[F_Entries[r].first];
						}
						if( E_Address[level]->at(k).Size() )
						{
							for (r = E_Address[level]->at(k).first; r < E_Address[level]->at(k).last; ++r)
								E_Entries[r].second /= DL0[k]*DR0[E_Entries[r].first];
						}
					}
					first = last;
				}
			}
		}
		
///////////////////////////////////////////////////////////////////////////////////
///                RESCALING DONE                                               ///
///////////////////////////////////////////////////////////////////////////////////
		tlrescale = Timer() - tt;
		trescale += tlrescale;
		
		level_interval.resize(level_size.size());
		if( !level_size.empty() )
		{
			level_interval[0].first = mobeg;
			level_interval[0].last = mobeg + level_size[0];
			for (k = 1; k < level_size.size(); k++)
			{
				level_interval[k].first = level_interval[k - 1].last;
				level_interval[k].last = level_interval[k].first + level_size[k];
			}
		}
///////////////////////////////////////////////////////////////////////////////////
///  FACTORIZATION COMPLETE                                                     ///
///////////////////////////////////////////////////////////////////////////////////
		ttotal = Timer() - ttotal;
		if( verbosity > 1 )
		{
			printf("total      %f\n",ttotal);
			printf("reorder    %f (%6.2f%%)\n", treorder, 100.0*treorder/ttotal);
			printf("   mpt     %f (%6.2f%%)\n",ttransversal, 100.0*ttransversal/ttotal);
#if defined(REORDER_METIS_ND)
			printf("   graph   %f (%6.2f%%)\n",tmetisgraph, 100.0*tmetisgraph/ttotal);
			printf("   nd      %f (%6.2f%%)\n", tmetisnd, 100.0*tmetisnd/ttotal);
#endif
#if defined(REORDER_RCM)
			printf("   graph   %f (%6.2f%%)\n",trcmgraph, 100.0*trcmgraph/ttotal);
			printf("   rcm     %f (%6.2f%%)\n", trcmorder, 100.0*trcmorder/ttotal);
#endif
			printf("reassamble %f (%6.2f%%)\n", treassamble, 100.0*treassamble / ttotal);
			printf("rescale    %f (%6.2f%%)\n", trescale, 100.0*trescale / ttotal);
			printf("factor     %f (%6.2f%%)\n", tfactor, 100.0*tfactor / ttotal);
			printf("   cond    %f (%6.2f%%)\n", testimator, 100.0*testimator / ttotal);
			printf("  schur    %f (%6.2f%%)\n", tschur, 100.0*tschur / ttotal);
			printf("nnz A %d LU %d LU2 %d swaps %d levels %d\n",nzA,nzLU,nzLU2tot, totswaps, (int)level_size.size());
		}
		init = true;
		/*
		{
			std::fstream fout("Diag.mtx",std::ios::out);
			fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
			fout << std::scientific;
			fout << moend-mobeg << " " << moend-mobeg << " " << moend-mobeg << std::endl;
			for(k = mobeg; k < moend; ++k)
				fout << k+1 << " " << k+1 << " " << LU_Diag[k] << std::endl;
			fout.close();
		}
		*/
		condestL = NuL;
		condestU = NuU;
		return true;
	}
	bool MLMTILUC_preconditioner::Finalize()
	{
		init = false;
		L_Address.clear();
		U_Address.clear();
		LU_Entries.clear();
		B_Entries.clear();
		B_Address.clear();
		temp.clear();
		ddP.clear();
		ddQ.clear();
		LU_Diag.clear();
		for(int k = 0; k < (int)E_Address.size(); ++k)
			delete E_Address[k];
		for(int k = 0; k < (int)F_Address.size(); ++k)
			delete F_Address[k];
		E_Address.clear();
		F_Address.clear();
		level_size.clear();
		level_interval.clear();
		E_Entries.clear();
		F_Entries.clear();
		return true;
	}
	void MLMTILUC_preconditioner::ApplyB(double alpha, Sparse::Vector & x, double beta, Sparse::Vector & y) // y = alpha A x + beta y
	{
		INMOST_DATA_ENUM_TYPE k, m, cbeg, cend;
		INMOST_DATA_REAL_TYPE temp;
		info->GetOverlapRegion(info->GetRank(), cbeg, cend);
		for(k = cbeg; k < cend; ++k)
		{
			temp = 0.0;
			for(m = B_Address[k].first; m < B_Address[k].last; ++m)
			{
				temp += B_Entries[m].second*x[B_Entries[m].first];
			}
			y[k] = beta*y[k] + alpha * temp;
		}
	}
	int MLMTILUC_preconditioner::Descend(int level, Sparse::Vector & inout)
	{
		INMOST_DATA_ENUM_TYPE k, m, cbeg, cend;
		//current B block interval
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		//check that it makes sense to make restriction step
		if( cend < level_interval.back().last )
		{
			// step 1) obtain ~f
			//apply B^{-1} on f
			for (k = cbeg; k < cend; ++k) temp[k] = inout[k];
			//Solve with L first
			for (k = cbeg; k < cend; ++k) if( L_Address[k].Size() )//iterator over columns of L
			{
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r)
					temp[LU_Entries[r].first] -= temp[k] * LU_Entries[r].second;
			}
			//Solve with diagonal
			for (k = cbeg; k < cend; ++k) temp[k] /= LU_Diag[k];
			//Solve with U
			for (k = cend; k > cbeg; --k) if( U_Address[k - 1].Size() ) //iterator over rows of U
			{
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r < U_Address[k - 1].last; ++r)
					temp[k - 1] -= temp[LU_Entries[r].first] * LU_Entries[r].second;
			}
			//now can calculate ~g
			//multiply ~f by row of E
			for (k = cend; k < level_interval.back().last; ++k) if( E_Address[level]->at(k).Size() )
			{
				for (m = E_Address[level]->at(k).first; m < E_Address[level]->at(k).last; ++m)
				{
					assert(E_Entries[m].first < cend);
					inout[k] -= temp[E_Entries[m].first] * E_Entries[m].second;
				}
			}
		}
		return level + 1;
	}
	int MLMTILUC_preconditioner::Ascend(int level, Sparse::Vector & inout)
	{
		level = level - 1;
		INMOST_DATA_ENUM_TYPE cbeg,cend, k, m;
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		//calculate Fy, iterate over ~g, write to unused input vector
		//for (k = cbeg; k < cend; ++k)
		//{
		//	for (m = F_Address[k].first; m < F_Address[k].last; m++)
		//		output[k] -= output[F_Entries[m].first] * F_Entries[m].second;
		//}
		//check that it makes sense to make prolongation step
		for (k = cend; k < level_interval.back().last; ++k) if( F_Address[level]->at(k).Size() )
		{
				for (m = F_Address[level]->at(k).first; m < F_Address[level]->at(k).last; ++m)
				{
					assert(F_Entries[m].first < cend);
					inout[F_Entries[m].first] -= inout[k] * F_Entries[m].second;
				}
		}
		//perform solve over calculated vector
		//Solve with L first
		for (k = cbeg; k < cend; ++k) if( L_Address[k].Size() )//iterator over columns of L
		{
			for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r)
				inout[LU_Entries[r].first] -= inout[k] * LU_Entries[r].second; //r->first alwayse > k
		}
		//Solve with diagonal
		for (k = cbeg; k < cend; ++k) inout[k] /= LU_Diag[k];
		//Solve with U
		for (k = cend; k > cbeg; --k) if( U_Address[k - 1].Size() )//iterator over rows of U
		{
			for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r < U_Address[k - 1].last; ++r)
				inout[k - 1] -= inout[LU_Entries[r].first] * LU_Entries[r].second; // r->first always > k
		}
		return level;
	}
	bool MLMTILUC_preconditioner::Solve(Sparse::Vector & input, Sparse::Vector & output)
	{
		assert(&input != &output);
		//
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			INMOST_DATA_ENUM_TYPE k, mobeg, moend, vbeg, vend;
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
			
			for (k = vbeg; k < mobeg; k++) output[k] = 0;
			for (k = mobeg; k < moend; ++k) output[k] = temp[ddP[k]];//*DL0[k];//*DL0[k];
			for (k = moend; k < vend; k++) output[k] = 0;
			
			int level = 0;
			while (level < level_size.size()) level = Descend(level, output);
			while (level > 0) level = Ascend(level, output);
			
			for (k = mobeg; k < moend; ++k) temp[ddQ[k]] = output[k];//*DR0[k];//*DR0[k];
			for (k = mobeg; k < moend; ++k) output[k] = temp[k];
			
			

			//Restrict additive schwartz (maybe do it outside?)
			//May assamble partition of unity instead of restriction before accumulation
			//assembly should be done instead of initialization
			for (k = vbeg; k < mobeg; ++k) output[k] = 0;
			for (k = moend; k < vend; ++k) output[k] = 0;
		}
		info->Accumulate(output);
		return true;
	}
	bool MLMTILUC_preconditioner::ReplaceMAT(Sparse::Matrix & A){ if (isInitialized()) Finalize(); Alink = &A; return true; }
	bool MLMTILUC_preconditioner::ReplaceSOL(Sparse::Vector & x) {(void)x; return true; }
	bool MLMTILUC_preconditioner::ReplaceRHS(Sparse::Vector & b) {(void)b; return true; }
	Method * MLMTILUC_preconditioner::Duplicate() { return new MLMTILUC_preconditioner(*this); }
	MLMTILUC_preconditioner::~MLMTILUC_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}
#endif

