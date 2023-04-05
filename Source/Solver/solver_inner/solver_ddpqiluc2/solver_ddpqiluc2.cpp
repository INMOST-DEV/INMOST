#define _CRT_SECURE_NO_WARNINGS

#include "inmost_solver.h"
#if defined(USE_SOLVER)
#include "solver_ddpqiluc2.hpp"
#include <iomanip>

/// \todo
/// 1. Correct estimator (as in MPTILUC)
/// 2. Remove premature dropping.

//#define PRINT_DDPQ
//#define PRINT_HISTOGRAM
//#define REPORT_ILU
//#undef REPORT_ILU
//#define REPORT_ILU_PROGRESS
//#define REPORT_ILU_END
//#define REPORT_ILU_SUMMARY
//#undef REPORT_ILU_PROGRESS
//#define REPORT_SCHUR
using namespace INMOST;

#ifndef DEFAULT_TAU
#define DEFAULT_TAU 0.01
#endif

//#define REORDER_METIS_ND
#define REORDER_DDPQ
#define DDPQ_TAU 0.8 //ratio is time consumed for factorization by time for schur complement calculation
//#define DDPQ_TAU std::max(0.1,1.0/(1.0+exp(-ratio*5)))
#define DDPQ_ADAPT_TAU
#define ESTIMATOR
#define RESCALE_B
#define PIVOT_THRESHOLD
#define DIAGONAL_PERTURBATION
#define DIAGONAL_PERTURBATION_REL 1.0e-6
#define DIAGONAL_PERTURBATION_ABS 1.0e-15
//#define DIAGONAL_PIVOT //probably there is some bug
#define DIAGONAL_PIVOT_TAU 0.01
//#define DIAGONAL_PIVOT_COND 5.0
#define ILUC2

#if defined(REORDER_METIS_ND)
#include "metis.h"
#endif


#define ADDR(row,col) ((row)*size+(col))


	void ILUC_preconditioner::DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
									std::vector<Sparse::Row::entry> & Entries,
									INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
									std::string file_name)
	{
		INMOST_DATA_REAL_TYPE norm1 = 0, norm2 = 0, max = -1.0e54, min = 1.0e54, minabs = 1.0e54;
		INMOST_DATA_ENUM_TYPE nnz = 0;
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k) 
		{
			nnz += Address[k].Size();

			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it != Address[k].last; ++it)
			{
				norm1 += fabs(Entries[it].second);
				norm2 += Entries[it].second*Entries[it].second;
				if( Entries[it].second > max ) max = Entries[it].second;
				if( Entries[it].second < min ) min = Entries[it].second;
				if( fabs(Entries[it].second) < minabs ) minabs = Entries[it].second;
			}
		}

		std::cout << "Writing matrix to " << file_name.c_str() << std::endl;
		std::fstream fout(file_name.c_str(),std::ios::out);
		fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
		fout << "% maximum " << max << std::endl;
		fout << "% minimum " << min << std::endl;
		fout << "% absolute minimum " << minabs << std::endl;
		fout << "% 1-norm  " << norm1 << std::endl;
		fout << "% 2-norm  " << sqrt(norm2) << std::endl;
		fout << "% mean 1-norm  " << norm1/nnz << std::endl;
		fout << "% mean 2-norm  " << sqrt(norm2/nnz) << std::endl;
		fout << std::scientific;
		
		fout << wmend-wmbeg << " " << wmend-wmbeg << " " << nnz << std::endl;;
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k)
		{
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it != Address[k].last; ++it)
				fout << (k-wmbeg+1) << " " << (Entries[it].first-wmbeg+1) << " " << Entries[it].second << std::endl;
		}
		fout.close();
	}
	void ILUC_preconditioner::SwapEntries(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 std::vector<Sparse::Row::entry> & Entries, 
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
	void ILUC_preconditioner::SwapLine(interval<INMOST_DATA_ENUM_TYPE, Interval> & Line, INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE temp;
		temp = Line[i].first;
		Line[i].first = Line[j].first;
		Line[j].first = temp;
		temp = Line[i].last;
		Line[i].last = Line[j].last;
		Line[j].last = temp;
	}
	void ILUC_preconditioner::SwapE(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE k;
		for (k = 0; k < E_Address.size(); ++k)
			SwapLine(*E_Address[k], i, j);
	}
	
	void ILUC_preconditioner::ReorderEF(INMOST_DATA_ENUM_TYPE mobeg, 
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
	void ILUC_preconditioner::inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
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
	void ILUC_preconditioner::applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
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
	INMOST_DATA_ENUM_TYPE & ILUC_preconditioner::EnumParameter(std::string name)
	{
		if (name == "scale_iters") return sciters;
		else if( name == "reorder_nnz") return reorder_nnz;
		else if( name == "ddpq_tau_adapt" ) return ddpq_tau_adapt;
		else if( name == "estimator" ) return estimator;
		throw - 1;
	}
	INMOST_DATA_REAL_TYPE & ILUC_preconditioner::RealParameter(std::string name)
	{
		if (name == "tau") return tau;
		else if( name == "ddpq_tau" ) return ddpq_tau;
		else if( name == "tau2" ) return iluc2_tau;
		throw - 1;
	}
	void ILUC_preconditioner::Copy(const Method * other)
	{
		const ILUC_preconditioner * b = dynamic_cast<const ILUC_preconditioner *>(other);
		assert(b != NULL);
		tau = b->tau;
		Alink = b->Alink;
		info = b->info;
		sciters = b->sciters;
		eps = b->eps;
	}
	ILUC_preconditioner::ILUC_preconditioner(const ILUC_preconditioner & other) :Method(other)
	{
		Copy(&other);
	}
	ILUC_preconditioner & ILUC_preconditioner::operator =(ILUC_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	ILUC_preconditioner::ILUC_preconditioner(Solver::OrderInfo & info)
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
		tau = 1.0e-3;
		iluc2_tau = tau*tau;
	}
	bool ILUC_preconditioner::isInitialized() { return init; }
	bool ILUC_preconditioner::isFinalized() { return !init; }
	bool ILUC_preconditioner::Initialize()
	{

#if defined(PRINT_HISTOGRAM)
		FILE * hist = fopen("histogram.csv","w");
#endif
#if defined(PRINT_DDPQ)
		FILE * fddpq = fopen("ddpq.csv","w");
#endif
#if defined(REPORT_SCHUR)
		FILE * fschur = fopen("schur.csv","w");
#endif
		const INMOST_DATA_REAL_TYPE subst = 1.0;
		const INMOST_DATA_REAL_TYPE tol_modif = 1e-12;
		const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;

		
		INMOST_DATA_ENUM_TYPE wbeg, wend; //working interval
		INMOST_DATA_ENUM_TYPE mobeg, moend; // total interval
		
		INMOST_DATA_ENUM_TYPE k, i, j, Li, Ui, curr, next, mi, mj, Sbeg;
		INMOST_DATA_REAL_TYPE l,u,udiag, abs_udiag, max_diag, min_diag, mean_diag, tol_schur;
		INMOST_DATA_ENUM_TYPE nzA, nzLU = 0, nzEF = 0, nzS;
		Sparse::Vector DL, DR;
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
		for(Sparse::Vector::iterator ri = DL.Begin(); ri != DL.End(); ++ri) *ri = 0.0;

		

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
		std::vector<Sparse::Row::entry> A_Entries, C_Entries, S_Entries, LF_Entries, LFt_Entries;
		//std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend), C_Address(mobeg,moend), S_Address(mobeg,moend), LF_Address(mobeg,moend), LFt_Address(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ulist(mobeg, moend), Llist(mobeg, moend), Blist(mobeg,moend),Alist(mobeg,moend),Afirst(mobeg,moend), Flist(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(mobeg, moend,EOL), Lbeg(mobeg, moend,EOL), Bbeg(mobeg,moend,EOL),Abeg(mobeg,moend,EOL), Fbeg(mobeg,moend,EOL);


		//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
#if defined(ESTIMATOR)
		INMOST_DATA_REAL_TYPE mup, mum, smup, smum, NuU = 1, NuL = 1, NuD, NuL1 = 1, NuL2 = 1, NuU1 = 1, NuU2 = 1;
		INMOST_DATA_REAL_TYPE NuU1_old = 1, NuL1_old = 1, NuU2_old = 1, NuL2_old = 1;
		INMOST_DATA_REAL_TYPE NuU1_new = 1, NuU2_new = 1, NuL1_new = 1, NuL2_new = 1, vp, vm, v;
		INMOST_DATA_ENUM_TYPE np, nm;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL1(mobeg, moend,0.0), EstU1(mobeg, moend,0.0), CondU(mobeg,moend,1.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL2(mobeg, moend,0.0), EstU2(mobeg, moend,0.0), CondL(mobeg,moend,1.0);
#endif
		//supplimentary data structures for returning values of dropped elements
		//INMOST_DATA_REAL_TYPE DropLk, DropUk;
		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
		//data structure for linked list

		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValuesL(mobeg, moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndecesL(mobeg, moend+1,UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValuesU(mobeg, moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndecesU(mobeg, moend+1,UNDEF);
		double tfactor = 0.0, trescale = 0.0, tschur = 0.0, treorder = 0.0, treassamble = 0.0, ttotal, tt;
		double tlfactor, tlrescale, tlschur, tlreorder, tlreassamble, ratio = 1.0, ratio2 = 1.0;
		ttotal = Timer();

		
		//calculate number of nonzeros
		nzA = 0;
		for (k = mobeg; k < moend; ++k)
		{
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
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
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
			{
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
					A_Entries[j++] = Sparse::Row::make_entry(r->first, r->second);
			}
			A_Address[k].last = j;
			//assert(A_Address[k].Size() != 0); //singular matrix
		}
		INMOST_DATA_REAL_TYPE  Anorm, Enorm, Fnorm;
		INMOST_DATA_ENUM_TYPE nzE,nzF;
#if defined(ILUC2)
		INMOST_DATA_ENUM_TYPE nzLU2 = 0;
		INMOST_DATA_REAL_TYPE tau2 = iluc2_tau;
		std::vector<Sparse::Row::entry> LU2_Entries;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Fnorms(mobeg,moend);
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
			std::fill(DL.Begin() + wbeg - mobeg, DL.Begin() + wend - mobeg, 0.0);
			std::fill(DR.Begin() + wbeg - mobeg, DR.Begin() + wend - mobeg, 0.0);
			std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, 0);
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

#if defined(PRINT_DDPQ)
			int histogram[32];
			memset(histogram,0,sizeof(int)*32);
			for (k = wbeg; k < wend; ++k)
			{
				//histogram[std::max(std::min(31-(int)floor((temp[k]-ddmaxmin)/(ddmaxall-ddmaxmin)*31),31),0)]++;
				int pos = std::min((int)ceil(-log((fabs(temp[k])-ddmaxmin)/(ddmaxall-ddmaxmin))),31);
				if( pos < 0 ) pos = 31;
				++histogram[ pos ];
			}
			fprintf(fddpq,"%g; %g; %g; %d",ddmaxall,ddmaxmin,ddmaxmean,histogram[0]);
			for(int hi = 1; hi < 32; ++hi)
				fprintf(fddpq,"; %d",histogram[hi]);
			fprintf(fddpq,"\n");
#endif
			
			//first select those elements that are largest over row and column simultaneously
			sort_wgts.clear();
			for (k = wbeg; k < wend; ++k)
			{
				if (Blist[k] != ENUMUNDEF)
				{
					//if ((temp[k]-ddmaxmin)/(ddmaxall-ddmaxmin) > ddpq_tau)
					//	sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					sort_wgts.push_back(wgt_coord(temp[k], coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(pow(temp[k],0.25)/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(1.0/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(temp[k]/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
				}
			}
			//if(reorder_nnz) 
				std::sort(sort_wgts.rbegin(), sort_wgts.rend());
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

				if( it->first / sort_wgts.begin()->first < ddpq_tau ) break;
				//if( it-sort_wgts.begin() > 0.01*(moend-mobeg) ) break; //5% of all values
			}
			cend = i;
#elif defined(REORDER_METIS_ND)
			idx_t nvtxs = wend-wbeg;
			std::vector<idx_t> xadj(nvtxs+1), adjncy, perm(nvtxs),iperm(nvtxs), sizes(nvtxs), vwgt(nvtxs);
			adjncy.reserve(nvtxs*5);
			xadj[0] = 0;

///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of B block                       //
///////////////////////////////////////////////////////////////////////////////////
			
			for (k = wend; k > wbeg; --k)
			{
				vwgt[k-1] = 1;
				if (A_Address[k-1].Size() > 0)
				{
					i = A_Entries[A_Address[k-1].first].first;
					Alist[k-1] = Abeg[i];
					Abeg[i] = k-1;
				}
				Afirst[k-1] = A_Address[k-1].first;
			}
			
			for(i = wbeg; i < wend; ++i)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt != A_Address[i].last; ++jt)
				{
					adjncy.push_back(static_cast<idx_t>(A_Entries[jt].first-wbeg));
				}
				
				Li = Abeg[i];
				while (Li != EOL)
				{
					adjncy.push_back(static_cast<idx_t>(Li-wbeg));
					Li = Alist[Li];
				}
				Li = Abeg[i];
				while (Li != EOL)
				{
					Ui = Alist[Li];
					Afirst[Li]++;
					if (A_Address[Li].last - Afirst[Li] > 0)
					{
						k = A_Entries[Afirst[Li]].first;
						Alist[Li] = Abeg[k];
						Abeg[k] = Li;
					}

					Li = Ui;
				}
				std::sort(adjncy.begin()+xadj[i],adjncy.end());
				adjncy.resize(std::unique(adjncy.begin()+xadj[i],adjncy.end())-adjncy.begin());
				
				xadj[i+1] = adjncy.size();
			}
			//A_Address[0].first = 0;
			//for(k = wbeg+1; k < wend; ++k)
			//	A_Address[k].first = A_Address[k-1].last;
			idx_t options[METIS_NOPTIONS];
			METIS_SetDefaultOptions(options);
			options[METIS_OPTION_NUMBERING] = 0;
			options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
			options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_TIME | METIS_DBG_COARSEN | METIS_DBG_CONNINFO | METIS_DBG_CONTIGINFO | METIS_DBG_IPART | METIS_DBG_MEMORY | METIS_DBG_MOVEINFO | METIS_DBG_REFINE | METIS_DBG_SEPINFO;
			printf("before metis\n");
			METIS_NodeNDP(nvtxs,&xadj[0],&adjncy[0],&vwgt[0],4,options,&perm[0],&iperm[0],&sizes[0]);
			//METIS_NodeND(&nvtxs,&xadj[0],&adjncy[0],&vwgt[0],options,&perm[0],&iperm[0]);
			printf("after metis\n");
			for(k = wbeg; k < wend; ++k)
			{
				localP[k] = perm[k-wbeg]+wbeg;
				localQ[k] = perm[k-wbeg]+wbeg;
			}
			cend = wend;
			i = wend;
#else
			cend = wend;
			i = wbeg;
#endif
			if (cbeg == cend && cbeg != wend)
			{
				//std::cout << __FILE__ << ":" << __LINE__ << " singular matrix, factored " << mobeg << ".." << cend << " out of " << mobeg << ".." << moend << std::endl;
				for (k = cbeg; k < wend; ++k)
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
			Anorm = 0;
			nzA = 0;
			for (k = cbeg; k < cend; ++k)
			{
				F_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(F_Entries.size());
				B_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt != A_Address[invP[k]].last; ++jt)
				{
					i = localQ[A_Entries[jt].first];
					u = A_Entries[jt].second;
					if (i < cend) B_Entries.push_back(Sparse::Row::make_entry(i, u));
					else F_Entries.push_back(Sparse::Row::make_entry(i,u)); //count nonzero in column of F
					Anorm += u*u;
				}
				F_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(F_Entries.size());
				B_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				std::sort(B_Entries.begin() + B_Address[k].first, B_Entries.end());
				std::sort(F_Entries.begin() + F_Address[k].first, F_Entries.end());
				nzEF += F_Address[k].Size();
				nzA += A_Address[k].Size();
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
					if (i < cend) E_Entries.push_back(Sparse::Row::make_entry(i, u)); //put to row of E
					else C_Entries.push_back(Sparse::Row::make_entry(i, u)); //form new Schur complement
					Anorm += u*u;
				}
				E_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries.size());
				C_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries.size());
				// sort entries added to E, because it will be very useful later
				std::sort(E_Entries.begin() + E_Address.back()->at(k).first, E_Entries.end());
				nzEF += E_Address.back()->at(k).Size();
				nzA += A_Address[k].Size();
			}
			Anorm = sqrt(Anorm/static_cast<INMOST_DATA_REAL_TYPE>(nzA));
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
			//for(Sparse::Vector::iterator it = DL.Begin() + cbeg - mobeg;
			//	it != DL.Begin() + cend - mobeg; ++it) *it = 0.0;
			std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
					DL[k] += fabs(B_Entries[r].second);
			}
			for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				//for(Sparse::Vector::iterator it = DR.Begin()+cbeg-mobeg; 
				//	it != DR.Begin()+cend-mobeg; ++it) *it = 0.0;
				std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DR[B_Entries[r].first] += DL[k] * fabs(B_Entries[r].second);
				}
				for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				//for(Sparse::Vector::iterator it = DL.Begin()+cbeg-mobeg;
				//	it != DL.Begin()+cend-mobeg; ++it) *it = 0.0;
				std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DL[k] += DR[B_Entries[r].first] * fabs(B_Entries[r].second);
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
#endif
			//End rescale B block
			trescale += Timer() - tt;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// FACTORIZATION BEGIN ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
#if defined(REPORT_ILU)
        int swaps = 0;
#endif
#if defined(ILUC2)
				nzLU2 = 0;
				LU2_Entries.clear();
#endif
				mean_diag = 0.0;
				(void)mean_diag;
///////////////////////////////////////////////////////////////////////////////////
//       setup diagonal values and memorize indices                              //
///////////////////////////////////////////////////////////////////////////////////
				//initialize diagonal values
				for (k = cbeg; k < cend; k++)
				{
          LU_Diag[k] = 0.0;
					localP[k] = F_Address[k].first; //remember F block starts permutation
					localQ[k] = B_Address[k].first; //remember B block starts permutation for pivoting
					for (i = B_Address[k].first; i < B_Address[k].last; i++)
					{
						if (B_Entries[i].first == k)
						{
							LU_Diag[k] = B_Entries[i].second;
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[k]);
#endif
							break;
						}
					}
				}
				//initialize first row and column
///////////////////////////////////////////////////////////////////////////////////
//       get diagonal value                                                      //
///////////////////////////////////////////////////////////////////////////////////
				assert(B_Entries[B_Address[cbeg].first].first == cbeg);
				if (B_Entries[B_Address[cbeg].first].first == cbeg)
					LU_Diag[cbeg] = B_Entries[B_Address[cbeg].first].second;
				else
				{
					std::cout << __LINE__ << " No diagonal value! " << cbeg << " " << B_Entries[B_Address[cbeg].first].first << std::endl;
					LU_Diag[cbeg] = 0.0;
				}
///////////////////////////////////////////////////////////////////////////////////
//       update for diagonal pivoting                                            //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PIVOT)
				mean_diag -= fabs(LU_Diag[cbeg]);
#endif
///////////////////////////////////////////////////////////////////////////////////
//       modify diagonal value according to rules                                //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//      form first line of U and second-order U                                  //
///////////////////////////////////////////////////////////////////////////////////
				LU_Beg = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
				U_Address[cbeg].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
        U2_Address[cbeg].first = 0;
#endif
				for (INMOST_DATA_ENUM_TYPE r = B_Address[cbeg].first + (B_Entries[B_Address[cbeg].first].first == cbeg ? 1 : 0); r != B_Address[cbeg].last; ++r) 
        {
          u = B_Entries[r].second / LU_Diag[cbeg];
          if( fabs(u) > tau )
					  LU_Entries.push_back(Sparse::Row::make_entry(B_Entries[r].first, u));
#if defined(ILUC2)
          else if( fabs(u) > tau2 )
            LU2_Entries.push_back(Sparse::Row::make_entry(B_Entries[r].first, u));
#endif
        }
        U_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());				
#if defined(ILUC2)
        U2_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());				
#endif
///////////////////////////////////////////////////////////////////////////////////
//       form first line of L and second-order L                                 //
///////////////////////////////////////////////////////////////////////////////////
        L_Address[cbeg].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
        L2_Address[cbeg].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());				
#endif
				Li = Bbeg[cbeg];
				while (Li != EOL)
				{
          l = B_Entries[B_Address[Li].first].second / LU_Diag[cbeg];
          if( fabs(l) > tau )
					  LU_Entries.push_back(Sparse::Row::make_entry(Li, l));
#if defined(ILUC2)
          else if( fabs(l) > tau2 )
            LU2_Entries.push_back(Sparse::Row::make_entry(Li, l));
#endif
					Li = Blist[Li];
				}
				L_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
        L2_Address[cbeg].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());				
#endif
				//update diagonal by factors
///////////////////////////////////////////////////////////////////////////////////
//       update whole diagonal with calculated L and U lines                     //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of B block                       //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of U                             //
///////////////////////////////////////////////////////////////////////////////////
				//form partial list from this row
				if (U_Address[cbeg].Size() > 0)
				{
					i = LU_Entries[U_Address[cbeg].first].first;
					Ubeg[i] = cbeg;
					Ulist[cbeg] = EOL;
				}
///////////////////////////////////////////////////////////////////////////////////
//       setup indices for transposed traversal of L                             //
///////////////////////////////////////////////////////////////////////////////////
				if (L_Address[cbeg].Size() > 0)
				{
					i = LU_Entries[L_Address[cbeg].first].first;
					Lbeg[i] = cbeg;
					Llist[cbeg] = EOL;
				}
				//Initialize data for condition estimator
				tt = Timer();
///////////////////////////////////////////////////////////////////////////////////
//              initialize condition estimator                                   //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
        if( estimator )
        {
				  NuU = NuL = 1.0;
				  if( estimator )
				  {
					  EstU1[cbeg] = EstL1[cbeg] = 0.0;
						EstU2[cbeg] = EstL2[cbeg] = 0.0;
					  for (INMOST_DATA_ENUM_TYPE it = U_Address[cbeg].first; it != U_Address[cbeg].last; ++it) 
						{
							EstU1[LU_Entries[it].first] = LU_Entries[it].second;
							EstU2[LU_Entries[it].first] = LU_Entries[it].second;
						}
					  for (INMOST_DATA_ENUM_TYPE it = L_Address[cbeg].first; it != L_Address[cbeg].last; ++it) 
						{
							EstL1[LU_Entries[it].first] = LU_Entries[it].second;
							EstL2[LU_Entries[it].first] = LU_Entries[it].second;
						}
				  }
        }
#endif
				nzLU += L_Address[cbeg].Size() + U_Address[cbeg].Size() + 1;
				max_diag = min_diag = fabs(LU_Diag[cbeg]);
				NuD = 1.0;
				(void)NuD;
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " starting factorization " << std::endl;
#endif
				for (k = cbeg + 1; k < cend; ++k)
				{
///////////////////////////////////////////////////////////////////////////////////
//              starting factorization step                                      //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//              diagonal pivoting algorithm                                      //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PIVOT_COND)
          bool no_swap_algorithm = false;
#endif
					if ( k < cend-1 &&  fabs(LU_Diag[k])*(cend-k) < DIAGONAL_PIVOT_TAU*mean_diag)
					{
#if defined(DIAGONAL_PIVOT_COND)
            if( false )
            {
swap_algorithm:
              no_swap_algorithm = true;
            }
#endif
						j = k+1;
						for (i = k + 1; i < cend; i++) if (fabs(LU_Diag[j]) < fabs(LU_Diag[i])) j = i;

            if( j > cend-1 )
            {
              std::cout << "Oops! asking to exchange beyond matrix! k = " << k << " j = " << j << " cend = " << cend << std::endl;
            }

            if( fabs(LU_Diag[j]) < fabs(LU_Diag[k]) )
            {
              std::cout << "Oops! I'm the best pivot!" << std::endl;
              
            }
            else
						{
							//TODO!!!
              ++swaps;
#if defined(REPORT_ILU)
							//std::cout << "Detected, that there is a much better pivot, i'm " << k << " " << LU_Diag[k] << " other " << j << " " << LU_Diag[j] << std::endl;
							//std::cout << "Condition numbers: L " << NuL << " D " << NuD << " U " << NuU << std::endl;
#endif
							//scanf("%*c");
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
#if defined(ESTIMATOR)
							if( estimator )
							{
								u = EstL1[k];
								EstL1[k] = EstL1[j];
								EstL1[j] = u;
								u = EstL2[k];
								EstL2[k] = EstL2[j];
								EstL2[j] = u;
								u = EstU1[k];
								EstU1[k] = EstU1[j];
								EstU1[j] = u;
								u = EstU2[k];
								EstU2[k] = EstU2[j];
								EstU2[j] = u;
							}
#endif
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
					
#endif
///////////////////////////////////////////////////////////////////////////////////
//            Prepare diagonal value                                             //
///////////////////////////////////////////////////////////////////////////////////
					udiag = LU_Diag[k]; // LU_Diag is calculated with less dropping
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
					abs_udiag = fabs(LU_Diag[k]);
					(void)abs_udiag;
///////////////////////////////////////////////////////////////////////////////////
//            Prepare linked list for U-part                                     //
///////////////////////////////////////////////////////////////////////////////////
					//DropLk = DropUk = 0.0;
					//uncompress k-th row
					// add diagonal value first, there shouldn't be values on left from diagonal
					assert(B_Entries[B_Address[k].first].first == k);
					LineIndecesU[cbeg] = k;
					if (B_Entries[B_Address[k].first].first == k)
						LineValuesU[k] = B_Entries[B_Address[k].first].second;
					else
					{
						std::cout << __LINE__ << " No diagonal value! " << k << " " << B_Entries[B_Address[k].first].first << std::endl;
						LineValuesU[k] = 0.0;
					}

					Ui = k;
					for (INMOST_DATA_ENUM_TYPE it = B_Address[k].first + (B_Entries[B_Address[k].first].first == k ? 1 : 0); it != B_Address[k].last; ++it)
					{
						LineValuesU[B_Entries[it].first] = B_Entries[it].second;
						Ui = LineIndecesU[Ui] = B_Entries[it].first;
					}
					LineIndecesU[Ui] = EOL;
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
            //if( fabs(l)*NuU*NuL_old > tau2*abs_udiag )//global check
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  u = LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= l*u;
							  else //if (fabs(u)*NuU*NuL_old > tau2*abs_udiag )//add new entry
							  {
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
								  LineValuesU[j] = -l*u;
							  }
						  }
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with L                        //
///////////////////////////////////////////////////////////////////////////////////
						//if( fabs(l)*NuU*NuL_old > tau*abs_udiag )//global check
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								u = l*LU2_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else //if (fabs(u)*NuU*NuL_old > tau2*abs_udiag)//add new entry
								{
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
            //if( fabs(l)*NuU*NuL_old > tau2*abs_udiag )//global check
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  u = l*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= u;
							  else //if (fabs(u)*NuU*NuL_old > tau2*abs_udiag)//add new entry
							  {
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
							  }
						  }
            }
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with second-order L           //
///////////////////////////////////////////////////////////////////////////////////
            //if( fabs(l)*NuU*NuL_old > tau*abs_udiag ) //global check
            {
              curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
						  {
                j = LU2_Entries[it].first;
							  u = l*LU2_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= u;
                else //if (fabs(u)*NuU*NuL_old > tau2*abs_udiag)//add new entry
							  {
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
							  }
						  }
            }

						i = L2list[i];
					}
#endif
					//std::cout << std::endl;
					//define the diagonal value
//this check will fail due to global check of tolerances, if global check will be removed then this will never fail
//#if defined(DIAGONAL_PIVOT) && !defined(DIAGONAL_PERTURBATION)
//					if (fabs(LU_Diag[k] - LineValues[k]) > 1e-6) std::cout << __LINE__ << " Diagonal value went wrong, good: " << LineValues[k] << " have " << LU_Diag[k] << " line " << k << std::endl;
//#endif
///////////////////////////////////////////////////////////////////////////////////
//                  Retrieve diagonal value                                      //
///////////////////////////////////////////////////////////////////////////////////
/*          
//#if defined(DIAGONAL_PIVOT) 
          udiag = LU_Diag[k]; // LU_Diag is calculated with less dropping
//#else
//					udiag = LineValues[k];
//#endif
          
          //udiag = LineValues[k];

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
					*/
///////////////////////////////////////////////////////////////////////////////////
//                    Condition estimator for diagonal D                         //
///////////////////////////////////////////////////////////////////////////////////
          if (fabs(udiag) > max_diag) max_diag = fabs(udiag);
					if (fabs(udiag) < min_diag) min_diag = fabs(udiag);
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
					if (B_Entries[B_Address[k].first].first == k)
						LineValuesL[k] = B_Entries[B_Address[k].first].second;
					else
					{
						std::cout << __LINE__ << " No diagonal value! " << k << std::endl;
						LineValuesL[k] = 0.0;
					}
					//start from diagonal
					//sort_indeces.clear();
					Ui = k;
					Li = Bbeg[k];
					while (Li != EOL)
					{
						assert(B_Entries[B_Address[Li].first].first == k);
						//sort_indeces.push_back(Li);
						LineValuesL[Li] = B_Entries[B_Address[Li].first].second;
						Ui = LineIndecesL[Ui] = Li;
						Li = Blist[Li];
					}
					/*
					std::sort(sort_indeces.begin(), sort_indeces.end());
					for (Li = 0; Li < sort_indeces.size(); Li++)
					Ui = LineIndeces[Ui] = sort_indeces[Li];
					*/
					LineIndecesL[Ui] = EOL;
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
            //if( fabs(u)*NuL*NuU_old > tau2*abs_udiag )//global check
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it < L_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  l = LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l*u;
							  else //if (fabs(l)*NuL*NuU_old > tau2*abs_udiag)//add new entry
							  {
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
								  LineValuesL[j] = -l*u;
							  }
						  }
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second-order L-part elimination with U                        //
///////////////////////////////////////////////////////////////////////////////////
						//if( fabs(u)*NuL*NuU_old > tau*abs_udiag )//global check
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								l = u*LU2_Entries[it].second;
								if (LineIndecesL[j] != UNDEF) //there is an entry in the list
									LineValuesL[j] -= l;
								else //if (fabs(l)*NuL*NuU_old > tau2*abs_udiag)//add new entry
								{
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
            //if( fabs(u)*NuL*NuU_old > tau2*abs_udiag )//global check
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it != L_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  l = u*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
							  else //if (fabs(l)*NuL*NuU_old > tau2*abs_udiag)//add new entry
							  {
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
                }
						  }
            }
///////////////////////////////////////////////////////////////////////////////////
//           second-order L-part elimination with second-order U                 //
///////////////////////////////////////////////////////////////////////////////////
            //if( fabs(u)*NuL*NuU_old > tau*abs_udiag )//global check
            {
              curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it != L2_Address[i].last; ++it)
						  {
							  j = LU2_Entries[it].first;
							  l = u*LU2_Entries[it].second;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
                else //if (fabs(l)*NuL*NuU_old > tau2*abs_udiag)//add new entry
							  {
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
                }
						  }
            }
						i = U2list[i];
					}
#endif
					//check that diagonal value is the same(must be!)
					//assert(fabs(LineValues[k] / udiag - 1.0) < 1.0e-10);
//this check will fail due to global check of tolerances, if global check will be removed then this will never fail
//#if defined(DIAGONAL_PIVOT) && !defined(DIAGONAL_PERTURBATION)
//					if (fabs(LineValues[k] - udiag) > 1e-9) std::cout << __LINE__ << " Diagonal value went wrong, good: " << udiag << " have " << LineValues[k] << " line " << k << std::endl;
//#endif
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
						NuU1_old = NuU1;
						(void)NuU1_old;
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
							if (vp > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) np++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vp, 0.5) < v) np--;
							if (vm > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) nm++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vm, 0.5) < v) nm--;
							Ui = LineIndecesU[Ui];
						}
						NuU2_old = NuU2;
						(void)NuU2_old;
						if (np > nm) NuU2 = mup; else NuU2 = mum;
						NuU2_new = std::max(fabs(mup), fabs(mum));
					} //if( estimator )
#endif //ESTIMATOR
///////////////////////////////////////////////////////////////////////////////////
//                 condition estimation for L part                               //
///////////////////////////////////////////////////////////////////////////////////
					//estimate condition for L^{-1}
#if defined(ESTIMATOR)
					if( estimator )
					{
						NuL1_old = NuL1;
						(void)NuL1_old;
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

						NuL2_old = NuL2;
						(void)NuL2_old;
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
							if (vp > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) np++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vp, 0.5) < v) np--;
							if (vm > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) nm++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vm, 0.5) < v) nm--;
							Li = LineIndecesL[Li];
						}
						if (np > nm) NuL2 = mup; else NuL2 = mum;
						NuL2_new = std::max(fabs(mup), fabs(mum));
					}
#endif //ESTIMATOR

///////////////////////////////////////////////////////////////////////////////////
//         discarding current iteration based on condition numbers               //
///////////////////////////////////////////////////////////////////////////////////
          //discard line if condition number is too large
#if defined(DIAGONAL_PIVOT) && defined(DIAGONAL_PIVOT_COND)
          if( k != cend-1 && no_swap_algorithm && 
						(std::max(NuU1_new,NuU2_new) > DIAGONAL_PIVOT_COND || std::max(NuL1_new,NuL2_new) > DIAGONAL_PIVOT_COND) )
          {
#if defined(REPORT_ILU)
            std::cout << "Requested pivoting based on condition estimator (L" << NuL1_new << " " << NuL2_new << " U " << NuU1_new << " " << NuU2_new << ")! row " << k << "/" << cend << std::endl;
#endif //REPORT_ILU
            //restore condition number
						NuL1 = NuU1_old;
            NuL2 = NuU2_old;
						NuU1 = NuU1_old;
            NuU2 = NuU2_old;
            //clean up values
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
            goto swap_algorithm;
          }
          else
#endif //DIAGONAL_PIVOT / DIAGONAL_PIVOT_COND
          {
						NuL = std::max(NuL1_new,NuL2_new);
						NuU = std::max(NuU1_new,NuU2_new);
          }
///////////////////////////////////////////////////////////////////////////////////
//                 reconstruct U-part from linked list                           //
///////////////////////////////////////////////////////////////////////////////////
          //insert line to U part
					U_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
          U2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					Ui = LineIndecesU[k];
					while (Ui != EOL)
					{
						u = fabs(LineValuesU[Ui]);
						if (u*NuU > tau) // apply dropping rule
							LU_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#if defined(ILUC2)
						else if (u*NuU > tau2)
							LU2_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#endif
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
					Li = LineIndecesL[k];
					L_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
#if defined(ILUC2)
          L2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU2_Entries.size());
#endif
					while (Li != EOL)
					{
						u = fabs(LineValuesL[Li]);
						if (u*NuL > tau) //apply dropping 
							LU_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#if defined(ILUC2)
						else if (u*NuL > tau2)
							LU2_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#endif
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
							Ui = LineIndecesU[k];
							while (Ui != EOL)
							{
								EstU2[Ui] += LineValuesU[Ui] * NuU2;
								Ui = LineIndecesU[Ui];
							}
							//choose new condition number
							NuU2 = NuU2_new;//std::max(NuU_new,NuU_old);
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
              Li = LineIndecesL[k];
						  while (Li != EOL)
						  {
							  EstL2[Li] += LineValuesL[Li] * NuL2;
							  Li = LineIndecesL[Li];
						  }
              //choose new condition number
              NuL2 = NuL2_new;//std::max(NuL_new,NuL_old);
						}
						CondU[k] = NuU;
						CondL[k] = NuL;
					}
#endif
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
//     Update diagonal                                                           //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PIVOT)
          mean_diag -= fabs(LU_Diag[k]); //remove current value
#endif
					//LU_Diag[k] = udiag;
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
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with L and U                      //
///////////////////////////////////////////////////////////////////////////////////
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
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
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
					while (i != U_Address[k].last && j != L2_Address[k].last)
					{
						Ui = LU_Entries[i].first;
						Li = LU2_Entries[j].first;
						if (Ui > Li) j++;
						else if (Ui < Li) i++;
						else
						{
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
							i++;
							j++;
						}
					}
///////////////////////////////////////////////////////////////////////////////////
//         Update values on the whole diagonal with second-order L and U         //
///////////////////////////////////////////////////////////////////////////////////
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
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU2_Entries[i].second * LU_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
							i++;
							j++;
						}
					}
///////////////////////////////////////////////////////////////////////////////////
// Update values on the whole diagonal with second-order L and second-order U    //
///////////////////////////////////////////////////////////////////////////////////
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
							assert(Ui > k);
#if defined(DIAGONAL_PIVOT)
							mean_diag -= fabs(LU_Diag[Ui]);
#endif
							LU_Diag[Ui] -= LU2_Entries[i].second * LU2_Entries[j].second * LU_Diag[k];
#if defined(DIAGONAL_PIVOT)
							mean_diag += fabs(LU_Diag[Ui]);
#endif
							i++;
							j++;
						}
					}
#endif //ILUC2
//#endif //DIAGONAL_PIVOT

					//number of nonzeros
					nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
					
#if defined(REPORT_ILU)
					if (k % 100 == 0)
					{
						
						std::cout << std::fixed << std::setprecision(2) << std::setw(6) << 100.0f*(k - cbeg) / (float)(cend - cbeg) << "%";
						std::cout << " nnz LU " << std::setprecision(10) << std::setw(10) << nzLU;
#if defined(ILUC2)
						std::cout << " LU2 " << std::setw(10) << nzLU2;
#endif
						std::cout << std::scientific << std::setprecision(5) << " condition L " << NuL << " D " << NuD << " U " << NuU;
            std::cout << " pivot swaps " << swaps;
						std::cout << "\r"; //carrige return
						std::cout.flush();
            //std::cout << std::endl;
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
///////////////////////////////////////////////////////////////////////////////////
//                       iteration done                                          //
///////////////////////////////////////////////////////////////////////////////////
				}
#if defined(REPORT_ILU_END)
#if defined(ILUC2)
				printf("%05.2f%% new level %6d nnz LU %10d LU2 %10d EF %10d cond L %8e D %8e U %8e pivot %6d\n",(cend-mobeg)/(1.0*(moend-mobeg))*100.0,cend-cbeg,nzLU,nzLU2,nzEF,NuL,NuD,NuU,swaps);
#else
				printf("%05.2f%% new level %6d nnz LU %10d EF %10d cond L %8e D %8e U %8e pivot %6d\n",(wend-mobeg)/(1.0*(cend-mobeg))*100.0,wend-wbeg,nzLU,nzEF,NuL,NuD,NuU,swaps);
#endif
        //std::cout << std::endl;
				//std::cout << "new level " << cend - cbeg << " total nonzeros in LU " << nzLU;
//#if defined(ILUC2)
//				std::cout << " in LU2 " << nzLU2;
//#endif
//				std::cout << " in EF " << nzEF << " conditions L " << NuL << " D " << NuD << " U " << NuU << " pivot swaps " << swaps << std::endl;
#endif
			}
			tlfactor = Timer() - tt;
			tfactor += tlfactor;
			//restore indexes
			B_Address[cbeg].first = localQ[cbeg];
			U_Address[cbeg].first = LU_Beg;
			L_Address[cbeg].first = U_Address[cbeg].last;
#if defined(ILUC2)
      U2_Address[cbeg].first = 0;
      L2_Address[cbeg].first = U2_Address[cbeg].last;
#endif
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
///////////////////////////////////////////////////////////////////////////////////
//         prepearing F block for transposed traversal                           //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//         put F block to linked list                                            //
///////////////////////////////////////////////////////////////////////////////////
					Li = cbeg;
					Ui = Fbeg[i];
					Fnorm = 0;
					nzF = 0;
					while (Ui != EOL)  //iterate over values of i-th column of F
					{
						assert(F_Entries[F_Address[Ui].first].first == i);
						LineValuesU[Ui] = F_Entries[F_Address[Ui].first].second;
						Fnorm += fabs(F_Entries[F_Address[Ui].first].second);
						//if( Fnorm < fabs(F_Entries[F_Address[Ui].first].second) ) Fnorm = fabs(F_Entries[F_Address[Ui].first].second);
						++nzF;
						Li = LineIndecesU[Li] = Ui + 1;
						Ui = Flist[Ui];
					}
					LineIndecesU[Li] = EOL;
///////////////////////////////////////////////////////////////////////////////////
//         shift indices for transposed traversal                                //
///////////////////////////////////////////////////////////////////////////////////					
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
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with L part                                            //
///////////////////////////////////////////////////////////////////////////////////
					//compute ~F_i = L^{-1} F_i
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						//if(fabs(LineValues[Li - 1])*NuU*NuL > tau2*Fnorm )
            {
              curr = Li;
						  for (INMOST_DATA_ENUM_TYPE ru = L_Address[Li - 1].first; ru < L_Address[Li - 1].last; ++ru)
						  {
							  u = LineValuesU[Li - 1] * LU_Entries[ru].second;
							  j = LU_Entries[ru].first;
							  if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
								  LineValuesU[j] -= u;
							  else //if( fabs(u)*NuU*NuL > tau2*Fnorm )
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
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with second-order L part                                //
///////////////////////////////////////////////////////////////////////////////////
						//if (fabs(LineValues[Li - 1])*NuU*NuL > tau*Fnorm)
						{
							curr = Li;
							for (INMOST_DATA_ENUM_TYPE ru = L2_Address[Li - 1].first; ru < L2_Address[Li - 1].last; ++ru)
							{
								u = LineValuesU[Li - 1] * LU2_Entries[ru].second;
								j = LU2_Entries[ru].first;
								if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
									LineValuesU[j] -= u;
								else //if( fabs(u)*NuU*NuL > tau2*Fnorm)
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
						}
#endif
						Li = LineIndecesU[Li];
///////////////////////////////////////////////////////////////////////////////////
//             solve iteration done                                              //
///////////////////////////////////////////////////////////////////////////////////
					}
///////////////////////////////////////////////////////////////////////////////////
//        assemble line of LF from linked list                                  //
///////////////////////////////////////////////////////////////////////////////////
					//Assemble column into matrix
					Li = LineIndecesU[cbeg];
					LF_Address[i].first = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
					while (Li != EOL)
					{
						if (fabs(LineValuesU[Li - 1]) > tau2*Fnorm)
							LF_Entries.push_back(Sparse::Row::make_entry(Li - 1, LineValuesU[Li - 1]));
						Li = LineIndecesU[Li];
					}
					LF_Address[i].last = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries.size());
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
#if defined(REPORT_ILU_PROGRESS)
					if (i % 100 == 0)
					{
						printf("%lu %d/%d schur1 %6.2f%%\t\t\r", static_cast<unsigned long>(level_size.size()), cend,moend, 100.f*(i - cend) / (1.f*(wend - cend)));
						fflush(stdout);
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//             iteration done!                                                   //
///////////////////////////////////////////////////////////////////////////////////
				}
				//Restore F indexing
				for (i = cbeg; i < cend; ++i)
					F_Address[i].first = localP[i];
				//setup row indeces for LF block
				//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = Fbeg.begin() + cbeg - mobeg;
				//	it != Fbeg.begin() + cend - mobeg; ++it ) *it = EOL;
///////////////////////////////////////////////////////////////////////////////////
//         prepearing LF block for transposed traversal                          //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//          transpose LF block                                                   //
///////////////////////////////////////////////////////////////////////////////////
				LFt_Entries.resize(LF_Entries.size());
				j = 0;
				for (k = cbeg; k < cend; k++)
				{
///////////////////////////////////////////////////////////////////////////////////
//         assemble line of transposed LF                                        //
///////////////////////////////////////////////////////////////////////////////////
					LFt_Address[k].first = j;
					Li = Fbeg[k];
					//Fnorms[k] = 0;
					nzF = 0;
					while (Li != EOL)
					{
						LFt_Entries[j].first = Li;
						LFt_Entries[j].second = LF_Entries[LF_Address[Li].first].second;
						//Fnorms[k] += fabs(LF_Entries[LF_Address[Li].first].second);
						++nzF;
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
				t1 = Timer() - t1;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// EU and Schur blocks  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				t2 = Timer();
        nzS = 0;
				//for(interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE>::iterator it = Ulist.begin() + cend - mobeg;
				//	it != Ulist.begin() + wend - mobeg; ++it) *it = UNDEF;
				
				S_Entries.clear();
				///
#if defined(REPORT_SCHUR)
				fprintf(fschur,"level %d\n",level_size.size());
				fprintf(fschur,"Sdiag; Smin; Smax; Snorm; Smean; Snum; Accepted; NuU; NuL; NuD; tol; dominance\n");
#endif
#if defined(PRINT_HISTOGRAM)
				//fprintf(hist,"level %d\n",level_size.size());
				int histogram[32];
				memset(histogram,0,sizeof(int)*32);
#endif
				///
				for (i = cend; i < wend; ++i) 
				{
///////////////////////////////////////////////////////////////////////////////////
//          no values at E block - unpack C block into Schur                     //
///////////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////////
//         put line of E block to linked list                                    //
///////////////////////////////////////////////////////////////////////////////////
					Enorm = 0;
					nzE = 0;
					Li = cbeg;
					for (j = E_Address.back()->at(i).first; j < E_Address.back()->at(i).last; ++j) //iterate over values of i-th row of E
					{
						LineValuesU[E_Entries[j].first] = E_Entries[j].second;
						Enorm += fabs(E_Entries[j].second);
						//if( Enorm < fabs(E_Entries[j].second) ) Enorm = fabs(E_Entries[j].second);
						++nzE;
						Li = LineIndecesU[Li] = E_Entries[j].first + 1;
					}
					LineIndecesU[Li] = EOL;
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with U part               //
///////////////////////////////////////////////////////////////////////////////////
					// compute ~E_i = E_i U^{-1}
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
            //if (fabs(LineValues[Li - 1])*NuL*NuU > tau2*Enorm*fabs(LU_Diag[Li-1]))
						{
						  curr = Li;
						  for (INMOST_DATA_ENUM_TYPE ru = U_Address[Li - 1].first; ru != U_Address[Li - 1].last; ++ru)
						  {
							  u = LineValuesU[Li - 1] * LU_Entries[ru].second;
							  j = LU_Entries[ru].first;
							  if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
								  LineValuesU[j] -= u;
							  else //if (fabs(u)*NuL*NuU > tau2*Enorm*fabs(LU_Diag[Li-1]))
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
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//         perform solve with second-order U part                                //
///////////////////////////////////////////////////////////////////////////////////
						//if (fabs(LineValues[Li - 1])*NuL*NuU > tau*Enorm*fabs(LU_Diag[Li-1]))
						{
							curr = Li;
							for (INMOST_DATA_ENUM_TYPE ru = U2_Address[Li - 1].first; ru != U2_Address[Li - 1].last; ++ru)
							{
								u = LineValuesU[Li - 1] * LU2_Entries[ru].second;
								j = LU2_Entries[ru].first;
								if (LineIndecesU[j + 1] != UNDEF) // There is an entry in the list
									LineValuesU[j] -= u;
								else //if (fabs(u)*NuL*NuU  > tau2*Enorm*fabs(LU_Diag[Li-1]) )
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
						}
#endif
						Li = LineIndecesU[Li];
///////////////////////////////////////////////////////////////////////////////////
//         solve iteration done                                                  //
///////////////////////////////////////////////////////////////////////////////////
					}
					
///////////////////////////////////////////////////////////////////////////////////
//drop values that do not satisfy tolerances from linked list of line of EU block//
///////////////////////////////////////////////////////////////////////////////////
					//drop values
					Li = LineIndecesU[cbeg];
					Ui = cbeg;
					while (Li != EOL)
					{
						j = LineIndecesU[Li];
						//if (fabs(LineValuesU[Li - 1])*NuL*NuU < tau2*Enorm*fabs(LU_Diag[Li-1]) )//tau2)
						if (fabs(LineValuesU[Li - 1]) < tau2*Enorm*fabs(LU_Diag[Li-1]) )
						{
							LineIndecesU[Ui] = j;
							LineIndecesU[Li] = UNDEF;
							LineValuesU[Li - 1] = 0.0;
						}
						else Ui = Li;
						Li = j;
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
					//unpack row of C
					tt0 = Timer() - tt0;

					tt1 = Timer();
///////////////////////////////////////////////////////////////////////////////////
//         unpack line of C matrix to temporary linked list                      //
///////////////////////////////////////////////////////////////////////////////////
					Sbeg = EOL;
					for (INMOST_DATA_ENUM_TYPE r = C_Address[i].first; r != C_Address[i].last; ++r)
					{
						LineValuesL[C_Entries[r].first] = C_Entries[r].second;
						LineIndecesL[C_Entries[r].first] = Sbeg;
						Sbeg = C_Entries[r].first;
					}
					//multiply current row of EU by LF and add to list
///////////////////////////////////////////////////////////////////////////////////
//        multiply EU and LF blocks and add to temporary linked list             //
///////////////////////////////////////////////////////////////////////////////////
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						//iterate over corresponding row of LF, add multiplication to list
						l = LineValuesU[Li - 1];
						for (INMOST_DATA_ENUM_TYPE r = LFt_Address[Li - 1].first; r < LFt_Address[Li - 1].last; ++r)
						{
							j = LFt_Entries[r].first;
							u = LFt_Entries[r].second;
							v = -l*u;
							if (LineIndecesL[j] != UNDEF)
								LineValuesL[j] += v;
							//else if( fabs(u)*NuU > tau || fabs(l)*NuL > tau )
							else //if( fabs(v) > 1.0e-54 )//if( fabs(v)*std::max(NuU,NuL) > tau2*std::max(Enorm,Fnorms[Li-1]) )//*Enorm*Fnorms[Li-1] )
							{
								LineValuesL[j] = v;
								LineIndecesL[j] = Sbeg;
								Sbeg = j;
							}
						}
						Li = LineIndecesU[Li];
					}
					tt1 = Timer() - tt1;

					tt2 = Timer();
///////////////////////////////////////////////////////////////////////////////////
//         calculate norms for dropping                                          //
///////////////////////////////////////////////////////////////////////////////////
					INMOST_DATA_REAL_TYPE Smax = 1.0e-54;
					INMOST_DATA_REAL_TYPE Snorm = 0.0, Snum = 0.0, Smin = 1.0e+54, Smean;
					Ui = Sbeg;
					while (Ui != EOL)
					{
						u = fabs(LineValuesL[Ui]);
						if( u > Smax ) Smax = u;
						if( u < Smin ) Smin = u;
						Snorm += u;
						Snum++;
						Ui = LineIndecesL[Ui];
					}
					Smean = Snorm/Snum;
					(void)Smean;
					//insert obtained row
///////////////////////////////////////////////////////////////////////////////////
//         put calculated row to Schur complement                                //
///////////////////////////////////////////////////////////////////////////////////
					//tol_schur = std::max(tau2,tau / (NuU*NuL*Enorm*Fnorm)) * Snorm ;
					//tol_schur = tau2 * Snorm * Snorm / Smax * NuU * NuL;
					
          //tol_schur = exp(log(tau2)*0.6+log(tau)*0.4) * Snorm * Snorm / Smax;// * sqrt(Smax*Snum/Snorm);//tau2 * Smax / std::max(NuU,NuL);
					//tol_schur = (Smax-Smin)*(std::min(tau2*sqrt(Smax/Smean),tau)) + Smin*0.999 ;
					//tol_schur = tau2;//*Snorm* sqrt(Smax/Smean);
					tol_schur = (Smax-Smin)*tau2 + Smin*0.99;
					//tol_schur = tau2*Smax*Smax/Snorm;
					//tol_schur = (Smax-Smin)*exp(log(tau2)*0.6+log(tau)*0.4) + Smin;
					//tol_schur = exp(log(tau2)*0.8+log(tau)*0.2)*Smax*Smax/Smean;
					//tol_schur = (Smax-Smin)*exp((log(tau2)*Smax+log(tau)*Snorm)/(Snorm+Smax)) + Smin;
					Ui = Sbeg;
					S_Address[i].first = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
					while (Ui != EOL)
					{
#if defined(PRINT_HISTOGRAM)
						//printf("val %g logval %g intlogval %d\n",(fabs(LineValuesL[Ui])-Smin)/(Smax-Smin), log((fabs(LineValuesL[Ui])-Smin)/(Smax-Smin)), (int)ceil(-log((fabs(LineValuesL[Ui])-Smin)/(Smax-Smin))));
						int pos = std::min((int)ceil(-log((fabs(LineValuesL[Ui])-Smin)/(Smax-Smin))),31);
						if( pos < 0 ) pos = 31;
						++histogram[ pos ];
						//++histogram[31- (int)floor((fabs(LineValuesL[Ui])-Smin)/(Smax-Smin)*31) ];
#endif
						//drop values
						if( fabs(LineValuesL[Ui]) >= tol_schur )
						{
						//if( fabs(LineValuesL[Ui]) > tau * Smax )// Snorm)
							S_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesL[Ui]));
							//fprintf(fschur, "%4d %20g accepted\n",Ui,LineValuesL[Ui]);
						}
						//else fprintf(fschur, "%4d %20g dropped\n",Ui,LineValuesL[Ui]);
						Ui = LineIndecesL[Ui];
					}
					S_Address[i].last = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries.size());
          nzS += S_Address[i].Size();

#if defined(REPORT_SCHUR)
					fprintf(fschur,"%e; %e; %e; %e; %e; %e; %d; %e; %e; %e; %e; %e\n",
						LineValuesL[i],Smin,Smax,Snorm,Smean,Snum,S_Address[i].Size(), NuU,NuL,NuD, tol_schur,Smax / (Snorm/sqrt(Snum)));
#endif
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
					Li = LineIndecesU[cbeg];
					while (Li != EOL)
					{
						j = LineIndecesU[Li];
						LineIndecesU[Li] = UNDEF;
						LineValuesU[Li - 1] = 0.0;
						Li = j;
					}
					tt2 = Timer() - tt2;
#if defined(REPORT_ILU)
					if (i % 10 == 0)
					{
            //carrige return
						printf("%6.2f%%  nzS %8d Smax %20g Snorm %20g tol %20g\r",
              100.f*(i - cend) / (1.f*(wend - cend)), nzS,Smax,Snorm,tol_schur);
            fflush(stdout);
					}
#elif defined(REPORT_ILU_PROGRESS)
					if (i % 100 == 0)
					{
						printf("%lu %d/%d schur2 %6.2f%%\t\t\r", static_cast<unsigned long>(level_size.size()),cend,moend,100.f*(i - cend) / (1.f*(wend - cend)));
						fflush(stdout);
					}
#endif
///////////////////////////////////////////////////////////////////////////////////
//         Schur complement row done                                             //
///////////////////////////////////////////////////////////////////////////////////
				}

#if defined(PRINT_HISTOGRAM)
				fprintf(hist,"%d",histogram[0]);
				for(int hi = 1; hi < 32; ++hi )
					fprintf(hist,"; %d",histogram[hi]);
				fprintf(hist,"\n");
#endif

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
			
			if ((1.0 * nzA / std::max<INMOST_DATA_ENUM_TYPE>(1,(wend - wbeg)*(wend - wbeg))) > 0.75 && (1.0*(wend - wbeg)) / (1.0 * std::max<INMOST_DATA_ENUM_TYPE>(1,moend - mobeg)) > 0.1)
			{
				std::cout << "Try to sparsify schur complement!!!" << std::endl;
				std::cout << "Sparsity: " << 1.0 * nzA / ((wend - wbeg)*(wend - wbeg)) << std::endl;
				std::cout << "Schur size: " << wend - wbeg << " total size: " << moend - mobeg << std::endl;

				//DumpMatrix(A_Address,A_Entries,wbeg,wend,"dense_schur.mtx");
				
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
        int swaps = 0;
				for (i = 0; i < size; ++i) //iterate over diagonal elements
				{
					//first step - pivoting, select largest element over entire schur complement
					mi = mj = i; //predefine maximum element
///////////////////////////////////////////////////////////////////////////////////
//                    Choose pivot                                               //
///////////////////////////////////////////////////////////////////////////////////
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
					
///////////////////////////////////////////////////////////////////////////////////
//                    Reorder for pivot                                          //
///////////////////////////////////////////////////////////////////////////////////
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
          if( i != mi || i != mj ) ++swaps;
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
///////////////////////////////////////////////////////////////////////////////////
//                    Apply threshold on pivot                                   //
///////////////////////////////////////////////////////////////////////////////////
					//TODO
					//check that diagonal element is not bad
					if (fabs(entries[ADDR(i, i)]) < 1e-15)
					{
						//std::cout << "warning: tiny diagonal!" << std::endl;
						entries[ADDR(i, i)] = entries[ADDR(i, i)] >= 0.0 ? 1.0e-15 : -1.0e-15;
					}
					//assert(fabs(entries[ADDR(i, i)])>1e-20);
					//scale off-diagonal elements by diagonal
///////////////////////////////////////////////////////////////////////////////////
//                   Rescale by pivot                                            //
///////////////////////////////////////////////////////////////////////////////////
					for (Li = i+1; Li < size; Li++) 
					{
						entries[ADDR(i, Li)] /= entries[ADDR(i, i)];
						entries[ADDR(Li, i)] /= entries[ADDR(i, i)];
					}
///////////////////////////////////////////////////////////////////////////////////
//                    Perform elimination from whole matrix                      //
///////////////////////////////////////////////////////////////////////////////////
					//compute Schur complement
					// S = E - l d u
					for (Li = i + 1; Li < size; Li++)
					for (Ui = i + 1; Ui < size; Ui++)
					{
						entries[ADDR(Li, Ui)] -= entries[ADDR(Li, i)] * entries[ADDR(i, i)] * entries[ADDR(i, Ui)];
					}
///////////////////////////////////////////////////////////////////////////////////
//                    Estimate condition numbers                                 //
///////////////////////////////////////////////////////////////////////////////////
					//run condition estimator
					//TODO
					//can it be done better on dense factors?
					//estimate condition for U^{-1}
#if defined(ESTIMATOR)
					if( estimator )
					{
						//estimate condition for U^{-1}
						NuU1_old = NuU1;
						mup = 1.0 - EstU1[wbeg+i];
						mum = -1.0 - EstU1[wbeg+i];
						smup = smum = 0;
						for (Ui = i + 1; Ui < size; Ui++)
						{
							v = EstU1[Ui];
							vp = fabs(v + entries[ADDR(i,Ui)] * mup);
							vm = fabs(v + entries[ADDR(i,Ui)] * mum);
							smup += vp;
							smum += vm;
						}
						if (smup > smum) NuU1 = mup; else NuU1 = mum;
						for (Ui = i + 1; Ui < size; Ui++) EstU1[wbeg + Ui] += entries[ADDR(i,Ui)] * NuU1;
						NuU1_new = std::max(fabs(mup),fabs(mum));

						NuU2_old = NuU2;
						mup = 1.0 - EstU2[wbeg + i];
						mum = -1.0 - EstU2[wbeg + i];
						np = nm = 0;
						//start from the element next after diagonal position
						for (Ui = i + 1; Ui < size; Ui++)
						{
							v = EstU2[wbeg + Ui];
							vp = fabs(v + entries[ADDR(i,Ui)] * mup);
							vm = fabs(v + entries[ADDR(i,Ui)] * mum);
							v = fabs(v);
							if (vp > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) np++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vp, 0.5) < v) np--;
							if (vm > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) nm++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vm, 0.5) < v) nm--;
						}
						if (np > nm) NuU2 = mup; else NuU2 = mum;
						NuU2_new = std::max(fabs(mup), fabs(mum));
						for (Ui = i + 1; Ui < size; Ui++) EstU2[wbeg + Ui] += entries[ADDR(i,Ui)] * NuU2;

						NuU = std::max(NuU1_new,NuU2_new);

						//estimate condition for L^{-1}
						NuL1_old = NuL1;
						mup = 1.0 - EstL1[wbeg+i];
						mum = -1.0 - EstL1[wbeg+i];
						smup = smum = 0;
						for (Li = i + 1; Li < size; Li++)
						{
							v = EstL1[Li];
							vp = fabs(v + entries[ADDR(Li,i)] * mup);
							vm = fabs(v + entries[ADDR(Li,i)] * mum);
							smup += vp;
							smum += vm;
						}
						if (smup > smum) NuL1 = mup; else NuL1 = mum;
						for (Li = i + 1; Li < size; Li++) EstL1[wbeg+Li] += entries[ADDR(Li,i)] * NuL1;
						NuL1_new = std::max(fabs(mup),fabs(mum));

						NuL2_old = NuL2;
						mup = 1.0 - EstL2[wbeg+i];
						mum = -1.0 - EstL2[wbeg+i];
						np = nm = 0;
						//start from the element next after diagonal position
						for (Li = i + 1; Li < size; Li++)
						{
							v = EstL2[wbeg+Li];
							vp = fabs(v + entries[ADDR(Li,i)] * mup);
							vm = fabs(v + entries[ADDR(Li,i)] * mum);
							v = fabs(v);
							if (vp > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) np++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vp, 0.5) < v) np--;
							if (vm > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) nm++;
							if (std::max<INMOST_DATA_REAL_TYPE>(2 * vm, 0.5) < v) nm--;
						}
						
						if (np > nm) NuL2 = mup; else NuL2 = mum;
						for (Li = i + 1; Li < size; Li++) EstL2[wbeg+Li] += entries[ADDR(Li,i)] * NuL2;
						NuL2 = std::max(fabs(mup), fabs(mum));

						NuL = std::max(NuL1_new,NuL2_new);

						CondU[i+wbeg] = NuU;
						CondL[i+wbeg] = NuL;
					}
#endif //ESTIMATOR
#if defined(REPORT_ILU)
					//if (k % 100 == 0)
					{
						
						std::cout << std::fixed << std::setprecision(2) << std::setw(6) << 100.0f*i / (float)size << "%";
						std::cout << std::scientific << std::setprecision(5) << " condition L " << NuL << " D " << NuD << " U " << NuU;
            std::cout << " pivot swaps " << swaps;
						std::cout << "\r"; //carrige return
						std::cout.flush();
            //std::cout << std::endl;

						std::cout << std::setprecision(0);
						//std::cout.setf(0,std::ios::floatfield);
						
						//printf("%6.2f%% nnz LU %8d condition L %10f D %10f U %10f\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, NuL, NuD, NuU);
						//fflush(stdout);
					}
#elif defined(REPORT_ILU_PROGRESS)
					if (k % 100 == 0)
					{
						printf("%lu %d/%d factor %6.2f%%\t\t\r",static_cast<unsigned long>(level_size.size()), wend,moend, 100.0f*i / (float)size);
						//printf("%6.2f%% nnz LU %8d condition L %10f D %10f U %10f\r", 100.0f*(k - cbeg) / (float)(cend - cbeg), nzLU, NuL, NuD, NuU);
						fflush(stdout);
					}
#endif
				}
#if defined(REPORT_ILU)
        std::cout << std::endl;
				std::cout << level_size.size() << " finishing reorder " << std::endl;
#endif
///////////////////////////////////////////////////////////////////////////////////
//                   finish reordering                                           //
///////////////////////////////////////////////////////////////////////////////////
				//apply PQ reordering to other part of matrices
				ReorderEF(mobeg,wbeg,wend,wend,donePQ,localP,localQ);
				inversePQ(wbeg, wend, localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);
				max_diag = min_diag = fabs(entries[ADDR(0, 0)]);
///////////////////////////////////////////////////////////////////////////////////
//                    Estimate condition number for diagonal                     //
///////////////////////////////////////////////////////////////////////////////////
				// copy computed decomposition to stored ILU, apply dropping rule
				for (i = 0; i < size; i++)
				{
					if (fabs(entries[ADDR(i, i)]) > max_diag) max_diag = fabs(entries[ADDR(i, i)]);
					if (fabs(entries[ADDR(i, i)]) < min_diag) min_diag = fabs(entries[ADDR(i, i)]);
				}
				NuD = max_diag / min_diag;
///////////////////////////////////////////////////////////////////////////////////
//                    Fill factors with account for tolerances                   //
///////////////////////////////////////////////////////////////////////////////////
#if defined(REPORT_ILU)
				std::cout << level_size.size() << " filling factors " << std::endl;
#endif
				for (i = 0; i < size; i++)
				{
					LU_Diag[wbeg + i] = entries[ADDR(i, i)];
					U_Address[wbeg + i].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					for (j = i + 1; j < size; j++)
					{
						if (fabs(entries[ADDR(i, j)])*CondU[i] > tau2)
						{
							LU_Entries.push_back(Sparse::Row::make_entry(j + wbeg, entries[ADDR(i, j)])); //Add to U
						}
					}
					U_Address[wbeg + i].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					L_Address[wbeg + i].first = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					for (j = i + 1; j < size; j++)
					{
						if (fabs(entries[ADDR(j, i)])*CondL[i] > tau2)
						{
							LU_Entries.push_back(Sparse::Row::make_entry(j + wbeg, entries[ADDR(j, i)])); //Add to L
						}
					}
					L_Address[wbeg + i].last = static_cast<INMOST_DATA_ENUM_TYPE>(LU_Entries.size());
					nzLU += 1 + L_Address[wbeg + i].Size() + U_Address[wbeg + i].Size();
				}
				delete[] entries;
				//Remember the size of the block
				level_size.push_back(wend - wbeg);
#if defined(REPORT_ILU_END)
				printf("%05.2f%% new level %6d nnz LU %10d EF %10d cond L %8e D %8e U %8e pivot %6d\n",100,wend-wbeg,nzLU,nzEF,NuL,NuD,NuU,swaps);
				//std::cout << "new level " << wend - wbeg << " total nonzeros in LU " << nzLU << " in EF " << nzEF << " conditions L " << NuL << " D " << NuD << " U " << NuU << std::endl;
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
				ddpq_tau = std::max<INMOST_DATA_REAL_TYPE>(0.01,ddpq_tau);
				ddpq_tau = std::min<INMOST_DATA_REAL_TYPE>(0.99,ddpq_tau);
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
#if defined(REPORT_ILU_SUMMARY)
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
#if defined(PRINT_HISTOGRAM)
		fclose(hist);
#endif
#if defined(PRINT_DDPQ)
		fclose(fddpq);
#endif
#if defined(REPORT_SCHUR)
		fclose(fschur);
#endif
		return true;
	}
	bool ILUC_preconditioner::Finalize()
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
	void ILUC_preconditioner::Multiply(int level, Sparse::Vector & input, Sparse::Vector & output)
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
	void ILUC_preconditioner::ApplyLU(int level, Sparse::Vector & input,Sparse::Vector & output)
	{
		INMOST_DATA_ENUM_TYPE k, cbeg, cend;
		//current B block interval
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;

		for (k = cbeg; k < cend; ++k) output[k] = input[k];
		//Solve with L first
		for (k = cbeg; k < cend; ++k) //iterator over columns of L
		{
			for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r != L_Address[k].last; ++r)
				output[LU_Entries[r].first] -= output[k] * LU_Entries[r].second;
		}
		//Solve with diagonal
		for (k = cbeg; k < cend; ++k) temp[k] /= LU_Diag[k];
		//Solve with U
		for (k = cend; k > cbeg; --k) //iterator over rows of U
		{
			for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r != U_Address[k - 1].last; ++r)
				output[k - 1] -= output[LU_Entries[r].first] * LU_Entries[r].second;
		}
	}
	void ILUC_preconditioner::ApplyB(int level, double alpha, Sparse::Vector & x, double beta, Sparse::Vector & y) // y = alpha A x + beta y
	{
		INMOST_DATA_ENUM_TYPE k, m, cbeg, cend;
		INMOST_DATA_REAL_TYPE temp;
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		for(k = cbeg; k < cend; ++k)
		{
			temp = 0.0;
			for(m = B_Address[k].first; m != B_Address[k].last; ++m)
			{
				temp += B_Entries[m].second*x[B_Entries[m].first];
			}
			y[k] = beta*y[k] + alpha * temp;
		}
	}
	int ILUC_preconditioner::Descend(int level, Sparse::Vector & input, Sparse::Vector & output)
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
	int ILUC_preconditioner::Ascend(int level, Sparse::Vector & input, Sparse::Vector & output)
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
	bool ILUC_preconditioner::Solve(Sparse::Vector & input, Sparse::Vector & output)
	{
		assert(&input != &output);
		//
#if defined(USE_OMP)
#pragma omp single
#endif
		{
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
			for (k = vbeg; k < mobeg; k++) output[k] = 0;//Restrict additive schwartz (maybe do it outside?)
			for (k = mobeg; k < moend; ++k) output[k] = temp[ddP[k]];
			for (k = moend; k < vend; k++) output[k] = 0;//Restrict additive schwartz (maybe do it outside?)

			level = 0;
			//perform recursively first two steps of solve phase
			// outer V-cycle

			while(level < level_size.size()) level = Descend(level, output, output);
		
			// W-cycle, should do ApplyB between moves
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
		}
		//May assemble partition of unity instead of restriction before accumulation
		//assembly should be done instead of initialization
		info->Accumulate(output);
		return true;
	}
	bool ILUC_preconditioner::ReplaceMAT(Sparse::Matrix & A){ if (isInitialized()) Finalize(); Alink = &A; return true; }
	bool ILUC_preconditioner::ReplaceSOL(Sparse::Vector & x) {(void)x; return true; }
	bool ILUC_preconditioner::ReplaceRHS(Sparse::Vector & b) {(void)b; return true; }
	Method * ILUC_preconditioner::Duplicate() { return new ILUC_preconditioner(*this); }
	ILUC_preconditioner::~ILUC_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}
#endif

