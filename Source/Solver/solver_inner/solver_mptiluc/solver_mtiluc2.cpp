#define _CRT_SECURE_NO_WARNINGS

#include "inmost_solver.h"
#if defined(USE_SOLVER)
#include "solver_mtiluc2.hpp"
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
//#define REORDER_DDPQ
//#define REORDER_ZOLTAN_HUND

#define ESTIMATOR
#define ESTIMATOR_REFINE
#define RESCALE_B

//#define PREMATURE_DROPPING

//#define EQUALIZE_1NORM
#define EQUALIZE_2NORM
//#define EQUALIZE_IDOMINANCE

#define PIVOT_THRESHOLD
#define PIVOT_THRESHOLD_VALUE 1.0e-9
//#define DIAGONAL_PERTURBATION
#define DIAGONAL_PERTURBATION_REL 1.0e-7
#define DIAGONAL_PERTURBATION_ABS 1.0e-10
//#define DIAGONAL_PIVOT //probably there is some bug
#define DIAGONAL_PIVOT_TAU 0.01
//#define DIAGONAL_PIVOT_COND 100
#define ILUC2
//#define TRACK_DIAGONAL

#if defined(DIAGONAL_PIVOT) && defined(DIAGONAL_PIVOT_TAU) && !defined(TRACK_DIAGONAL)
#define TRACK_DIAGONAL
#endif


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



	template<typename T>
	std::ostream & fmt(std::ostream & in, const T & value, int width, int precision)
	{
		std::ios save(NULL);
		save.copyfmt(in);
		in << std::fixed << std::setw(width) << std::setprecision(precision) << value;
		std::cout.copyfmt(save);
		return in;
	}

	void MTILUC_preconditioner::DumpMatrix(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
									std::vector<Sparse::Row::entry> & Entries,
									INMOST_DATA_ENUM_TYPE wmbeg, INMOST_DATA_ENUM_TYPE wmend,
									std::string file_name)
	{
		INMOST_DATA_REAL_TYPE norm1 = 0, norm2 = 0, max = -1.0e54, min = 1.0e54, minabs = 1.0e54, vrowmax, diag, mindiag = 1.0e54, maxdiag = -1.0e54, maxabsdiag = -1.0e54, minabsdiag = 1.0e54;
		INMOST_DATA_ENUM_TYPE nnz = 0, dominant_rows = 0, dominant_cols = 0, irowmax = 0, nodiag = 0;
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> vcolmax(wmbeg,wmend,0);
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> icolmax(wmbeg,wmend,ENUMUNDEF);
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k) 
		{
			nnz += Address[k].Size();
			vrowmax = 0;

			bool diag_found = false;
			diag = 0;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it != Address[k].last; ++it)
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
		
		fout << wmend-wmbeg << " " << wmend-wmbeg << " " << nnz << std::endl;;
		for (INMOST_DATA_ENUM_TYPE k = wmbeg; k < wmend; ++k)
		{
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it != Address[k].last; ++it)
				fout << (k-wmbeg+1) << " " << (Entries[it].first-wmbeg+1) << " " << Entries[it].second << std::endl;
		}
		fout.close();
	}
	void MTILUC_preconditioner::CheckOrder(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 std::vector<Sparse::Row::entry> & Entries, 
					 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend)
	{
		INMOST_DATA_ENUM_TYPE i,r;
		for (i = rbeg; i < rend; i++) if( Address[i].first != Address[i].last )
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
	void MTILUC_preconditioner::SwapEntries(interval<INMOST_DATA_ENUM_TYPE, Interval> & Address, 
					 std::vector<Sparse::Row::entry> & Entries, 
					 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend, 
					 INMOST_DATA_ENUM_TYPE k, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE i;
		INMOST_DATA_REAL_TYPE u;
		assert(k < j);
		for (i = rbeg; i < rend; i++) if( Address[i].last != Address[i].last )
		{
			INMOST_DATA_ENUM_TYPE pk = Address[i].last, pj = Address[i].last, r;
			//check ordered on entry
			for (r = Address[i].first; r < Address[i].last-1; r++)
			{
				assert(Entries[r].first < Entries[r+1].first);
			}
			
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				if (Entries[r].first >= k) { pk = r; break; }
			}
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				if (Entries[r].first >= j) { pj = r; break; }
			}
			bool found_k = pk < Address[i].last && Entries[pk].first == k;
			bool found_j = pj < Address[i].last && Entries[pj].first == j;
			
			if( found_k && found_j )
			{
				//just swap values, indices are intact
				u = Entries[pj].second;
				Entries[pj].second = Entries[pk].second;
				Entries[pk].second = u;
			}
			else if( found_k )
			{
				//we have to place value of k into position of j and give it an index of j
				//don't move at pj, since it's index is greater then j
				u = Entries[pk].second;
				for (r = pk; r < pj - 1; r++) Entries[r] = Entries[r + 1];
				Entries[pj - 1].first = j;
				Entries[pj - 1].second = u;
				
				//check ordered
				for (r = Address[i].first; r < Address[i].last-1; r++)
				{
					assert(Entries[r].first < Entries[r+1].first);
				}
				
			}
			else if( found_j )
			{
				//we have to place value of j into position of k and give it an index of k
				// move at pk, since it's index is greater then k
				u = Entries[pj].second;
				for (r = pj; r > pk; r--) Entries[r] = Entries[r - 1];
				Entries[pk].first = k;
				Entries[pk].second = u;
				
				//check ordered
				for (r = Address[i].first; r < Address[i].last-1; r++)
				{
					assert(Entries[r].first < Entries[r+1].first);
				}
			}
			//check nan
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				assert(Entries[r].second == Entries[r].second);
			}
			
		}
	}
	void MTILUC_preconditioner::SwapLine(interval<INMOST_DATA_ENUM_TYPE, Interval> & Line, INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
	{
		INMOST_DATA_ENUM_TYPE temp;
		temp = Line[i].first;
		Line[i].first = Line[j].first;
		Line[j].first = temp;
		temp = Line[i].last;
		Line[i].last = Line[j].last;
		Line[j].last = temp;
	}
	
	void MTILUC_preconditioner::inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
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
	void MTILUC_preconditioner::applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
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
	INMOST_DATA_ENUM_TYPE & MTILUC_preconditioner::EnumParameter(std::string name)
	{
		if (name == "scale_iters") return sciters;
		else if( name == "estimator" ) return estimator;
		else if( name == "verbosity" ) return verbosity;
		throw - 1;
	}
	INMOST_DATA_REAL_TYPE & MTILUC_preconditioner::RealParameter(std::string name)
	{
		if (name == "tau") return tau;
		else if( name == "tau2" ) return iluc2_tau;
    else if( name == "condition_number_L" ) return condestL;
    else if( name == "condition_number_U" ) return condestU;
		throw - 1;
	}
	void MTILUC_preconditioner::Copy(const Method * other)
	{
		const MTILUC_preconditioner * b = dynamic_cast<const MTILUC_preconditioner *>(other);
		assert(b != NULL);
		tau = b->tau;
		Alink = b->Alink;
		info = b->info;
		sciters = b->sciters;
		eps = b->eps;
		verbosity = b->verbosity;
	}
	MTILUC_preconditioner::MTILUC_preconditioner(const MTILUC_preconditioner & other) :Method(other)
	{
		Copy(&other);
	}
	MTILUC_preconditioner & MTILUC_preconditioner::operator =(MTILUC_preconditioner const & other)
	{
		Copy(&other);
		return *this;
	}
	MTILUC_preconditioner::MTILUC_preconditioner(Solver::OrderInfo & info)
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
		verbosity = 0;
		tau = 1.0e-3;
		iluc2_tau = tau*tau;
	}
	bool MTILUC_preconditioner::isInitialized() { return init; }
	bool MTILUC_preconditioner::isFinalized() { return !init; }
	bool MTILUC_preconditioner::Initialize()
	{

		const INMOST_DATA_REAL_TYPE subst = 1.0;
		const INMOST_DATA_REAL_TYPE tol_modif = PIVOT_THRESHOLD_VALUE;
		const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;

		
		INMOST_DATA_ENUM_TYPE wbeg = 0, wend = 0; //working interval
		INMOST_DATA_ENUM_TYPE mobeg = 0, moend = 0; // total interval
		INMOST_DATA_ENUM_TYPE vbeg = 0, vend = 0; // vector interval
		
		INMOST_DATA_ENUM_TYPE k = 0, i = 0, j = 0, Li = 0, Ui = 0, curr = 0, next = 0;
		INMOST_DATA_REAL_TYPE l = 0, u = 0,udiag = 0, abs_udiag = 0, max_diag = -1.0e+20, min_diag = 1.0e+20, mean_diag = 0;
		INMOST_DATA_ENUM_TYPE nzA = 0, nzLU = 0;
		Sparse::Vector DL, DR;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetVectorRegion(vbeg,vend);
		
		//prepare temporal array
		temp.set_interval_beg(vbeg);
		temp.set_interval_end(vend);
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
		for(Sparse::Vector::iterator ri = DL.Begin(); ri != DL.End(); ++ri) *ri = 1.0;
		for(Sparse::Vector::iterator ri = DR.Begin(); ri != DR.End(); ++ri) *ri = 1.0;

		

		for (k = mobeg; k != moend; k++) ddP[k] = ddQ[k] = k;
		// supplementary data for ddPQ reordering
		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(mobeg, moend), invQ(mobeg,moend), localP(mobeg, moend,ENUMUNDEF), localQ(mobeg,moend,ENUMUNDEF);		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Bstart(mobeg,moend);

		//supplimentary data structures for ILUC
		INMOST_DATA_ENUM_TYPE LU_Beg = 0;
		U_Address.set_interval_beg(mobeg);
		U_Address.set_interval_end(moend);
		L_Address.set_interval_beg(mobeg);
		L_Address.set_interval_end(moend);
		B_Address.set_interval_beg(mobeg);
		B_Address.set_interval_end(moend);
		LU_Diag.set_interval_beg(mobeg);
		LU_Diag.set_interval_end(moend);
		std::vector<Sparse::Row::entry> A_Entries;
		//std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ulist(mobeg, moend), Llist(mobeg, moend), Blist(mobeg,moend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(mobeg, moend,EOL), Lbeg(mobeg, moend,EOL), Bbeg(mobeg,moend,EOL);

		INMOST_DATA_REAL_TYPE NuU = 1, NuL = 1, NuD, NuU_max = 1.0, NuL_max = 1.0;
#if defined(ESTIMATOR)
		//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
		INMOST_DATA_REAL_TYPE mup = 0, mum = 0, smup = 0, smum = 0, NuL1 = 1, NuL2 = 1, NuU1 = 1, NuU2 = 1;
		INMOST_DATA_REAL_TYPE NuU1_old = 1, NuL1_old = 1, NuU2_old = 1, NuL2_old = 1;
		INMOST_DATA_REAL_TYPE NuU1_new = 1, NuL1_new = 1, vp = 0, vm = 0, v = 0;
#if defined(ESTIMATOR_REFINE)
		INMOST_DATA_ENUM_TYPE np = 0, nm = 0;
		INMOST_DATA_REAL_TYPE NuU2_new = 0, NuL2_new = 0;
#endif
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL1(mobeg, moend,0.0), EstU1(mobeg, moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL2(mobeg, moend,0.0), EstU2(mobeg, moend,0.0);
#endif
#if defined(ESTIMATOR) && defined(DIAGONAL_PIVOT_COND)
		INMOST_DATA_REAL_TYPE NuU_tmp = 0, NuL_tmp = 0;
#endif
		//supplimentary data structures for returning values of dropped elements
		//INMOST_DATA_REAL_TYPE DropLk, DropUk;
		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
		//data structure for linked list

		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Dist(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());		
		
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValuesU(mobeg, moend,0.0), LineValuesL(mobeg,moend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndecesU(mobeg, moend+1,UNDEF), LineIndecesL(mobeg,moend+1,UNDEF);
		double tfactor = 0.0, tswap = 0.0, trescale = 0.0, treorder = 0.0, ttransversal = 0.0, treassamble = 0.0, ttotal, tt = 0.0, testimator = 0.0, tlocal = 0.0;
        (void) tswap; //can be unused for some defines
#if defined(REORDER_METIS_ND)
		double tmetisgraph = 0, tmetisnd = 0;
#endif
#if defined(REORDER_RCM)
		double trcmgraph = 0, trcmorder = 0;
#endif
        double tlfactor, tlrescale, tlreorder, tlreassamble;
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
					if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
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
				if (r->first >= mobeg && r->first < moend && fabs(r->second) > 0.0)
					A_Entries[j++] = Sparse::Row::make_entry(r->first, r->second);
			}
			A_Address[k].last = j;
			//assert(A_Address[k].Size() != 0); //singular matrix
		}
#endif
		//DumpMatrix(A_Address, A_Entries, mobeg, moend, "A.mtx");

		std::vector<INMOST_DATA_REAL_TYPE> C_Entries(A_Entries.size());

		//~ INMOST_DATA_ENUM_TYPE Sbeg;
		std::vector<INMOST_DATA_ENUM_TYPE> indices;
#if defined(ILUC2)
		INMOST_DATA_ENUM_TYPE nzLU2 = 0;
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
#if defined(ILUC2)
			U2_Address[k].first = U2_Address[k].last = ENUMUNDEF;
			L2_Address[k].first = L2_Address[k].last = ENUMUNDEF;
#endif
		}


		if( verbosity > 1 ) 
			std::cout << "nonzeros in A " << nzA << std::endl;
		int swaps = 0;
		//set up working interval
		wbeg = mobeg;
		wend = moend;
		while (wbeg < wend) //iterate into levels until all matrix is factored
		{
			//std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
			//std::fill(localQ.begin() + (wbeg - mobeg), localQ.begin() + (wend - mobeg), ENUMUNDEF);
			//ddPQ reordering on current Schur complement
			INMOST_DATA_ENUM_TYPE cbeg = wbeg, cend = wend; //next size of factored B block
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// MAXIMUM TRANSVERSE REORDERING /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Arrays initialization   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
				//std::fill(U.Begin() + wbeg - mobeg, U.Begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(V.Begin() + wbeg - mobeg, V.Begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				//std::fill(Cmax.begin() + wbeg - mobeg, Cmax.begin() + wend - mobeg, 0.0);
				//std::fill(Dist.begin() + wbeg - mobeg, Dist.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				std::fill(Perm.begin() + wbeg - mobeg, Perm.begin() + wend - mobeg, ENUMUNDEF);
				std::fill(IPerm.begin() + wbeg - mobeg, IPerm.begin() + wend - mobeg, ENUMUNDEF);
				std::fill(Parent.begin() + wbeg - mobeg, Parent.begin() + wend - mobeg, ENUMUNDEF);
				std::fill(ColumnList.begin() + wbeg - mobeg, ColumnList.begin() + wend - mobeg, ENUMUNDEF);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Initial LOG transformation to dual problem and initial extreme match /////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//double T = Timer();
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						u = C_Entries[it] = fabs(A_Entries[it].second);
						if( u > Cmax[i] ) Cmax[i] = u;
						//C_Entries.push_back(Sparse::Row::make_entry(i,u));
					}
				}

				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
					{
						i = A_Entries[it].first;
						if( Cmax[i] == 0.0 || C_Entries[it] == 0.0 )
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
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
					{
						u = C_Entries[it] - U[A_Entries[it].first];
						if( u < V[k] ) V[k] = u;
					}
				}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Update cost and match ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// 1-step augmentation   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				for(k = wbeg; k < wend; ++k)
				{
					if( IPerm[k] == ENUMUNDEF ) //unmatched row
					{
						for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last && IPerm[k] == ENUMUNDEF; ++it)
						{
							u = fabs(C_Entries[it] - V[k] - U[A_Entries[it].first]);
							if( u <= 1.0e-30 )
							{
								Li = Perm[A_Entries[it].first];
								assert(Li != ENUMUNDEF);
								// Search other row in C for 0
								for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit != A_Address[Li].last; ++Lit)
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

				
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Weighted bipartite matching ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
						for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit != A_Address[Li].last; ++Lit)
						{
							Ui = A_Entries[Lit].first;
							//if( ColumnList[Ui] == k ) continue;
							if( ColumnList[Ui] != ENUMUNDEF ) continue;
							l = ShortestPath + C_Entries[Lit] - V[Li] - U[Ui];
							//if( l < 0.0 ) printf("row %d col %d negative l %g Shortest %lf C %lf V %lf U %lf\n",k,Ui,l,C_Entries[Lit],V[Li],U[Ui]);
							if( l < 0.0 && fabs(l) < 1.0e-12 ) l = 0;
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
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt != A_Address[k].last; ++jt)
					{
						i = A_Entries[jt].first;
						j = Perm[A_Entries[jt].first];
						if( U[i] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || Cmax[i] == 0 ) u = 1;
						else u = exp(U[i])/Cmax[i];
						if( j == k && (l*A_Entries[jt].second*u) < 0.0 ) flip_sign = true;
					}
					
					if( flip_sign ) DL[k] *= -1;
				}
				
				//DumpMatrix(A_Address,A_Entries,wbeg,wend,"A.mtx");
				//DumpMatrix(B_Address,B_Entries,wbeg,wend,"MC64.mtx");

				std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
                
                { //check that there are no gaps in Perm
                    for(k = cbeg; k < cend; ++k)
                    {
                        if( Perm[k] != ENUMUNDEF )
                            localP[Perm[k]] = 0;
                    }
                    std::vector<INMOST_DATA_ENUM_TYPE> gaps;
                    for(k = cbeg; k < cend; ++k)
                        if( localP[k] == ENUMUNDEF )
                            gaps.push_back(k);
                    
                    for(k = cbeg; k < cend; ++k)
                        if( Perm[k] == ENUMUNDEF )
                        {
                            Perm[k] = gaps.back();
                            gaps.pop_back();
                        }
                    std::fill(localP.begin() + (wbeg - mobeg), localP.begin() + (wend - mobeg), ENUMUNDEF);
                }
				
				//exit(-1);

				ttransversal = Timer() - ttransversal;

				treorder += ttransversal;
			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// END MAXIMUM TRANSVERSE REORDERING /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
#if defined(REORDER_DDPQ)

			if( verbosity > 1 )
				std::cout << level_size.size() << " calculate weights" << std::endl;
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

			
			//first select those elements that are largest over row and column simultaneously
			sort_wgts.clear();
			for (k = wbeg; k < wend; ++k)
			{
				if (Blist[k] != ENUMUNDEF)
				{
					//if ((temp[k]-ddmaxmin)/(ddmaxall-ddmaxmin) > ddpq_tau)
					//	sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(temp[k], coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(pow(temp[k],0.25)/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					//sort_wgts.push_back(wgt_coord(1.0/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
					sort_wgts.push_back(wgt_coord(temp[k]/(Ulist[k]+Llist[Blist[k]]-1), coord(k, Blist[k])));
				}
			}
			//if(reorder_nnz) 
				std::sort(sort_wgts.rbegin(), sort_wgts.rend());
			if( verbosity > 1 )
				std::cout << level_size.size() << " fill reordering" << std::endl;

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
				if( it-sort_wgts.begin() > 0.01*(moend-mobeg) ) break; //5% of all values
			}
			cend = i;
#elif defined(REORDER_METIS_ND)
			tt = Timer();
			idx_t nvtxs = wend-wbeg;
			std::vector<idx_t> xadj(nvtxs+1), adjncy, perm(nvtxs),iperm(nvtxs);
			adjncy.reserve(nzA*2);
			
			for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";

			for (k = cbeg; k < cend; ++k) localP[k] = k;
			for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
			for (k = cbeg; k < cend; ++k) DR[k] = V[k];

			
			
			for(k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt != A_Address[k].last; ++jt)
				{
					A_Entries[jt].first = localQ[A_Entries[jt].first];
					//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					
				}
				std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
			}

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
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt != A_Address[i].last; ++jt)
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

			std::fill(Bbeg.begin(),Bbeg.end(),EOL);

			tmetisgraph = Timer() - tmetisgraph;
			
			//A_Address[0].first = 0;
			//for(k = wbeg+1; k < wend; ++k)
			//	A_Address[k].first = A_Address[k-1].last;
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
			//METIS_NodeNDP(nvtxs,&xadj[0],&adjncy[0],NULL,4,options,&perm[0],&iperm[0],&sizes[0]);
			tmetisnd = Timer();
			METIS_NodeND(&nvtxs,&xadj[0],&adjncy[0],NULL,options,&perm[0],&iperm[0]);
			tmetisnd = Timer()-tmetisnd;
			for(k = wbeg; k < wend; ++k)
			{
				//localP[k] = iperm[k-wbeg]+wbeg;
				//localQ[k] = iperm[k-wbeg]+wbeg;
				localP[perm[k-wbeg]+wbeg] = k;
				localQ[perm[k-wbeg]+wbeg] = k;
			}
			cend = wend;
			i = wend;

#elif defined(REORDER_NNZ)
			//for (k = cbeg; k < cend; ++k) V[k] = DR[localQ[k]];
			//for (k = cbeg; k < cend; ++k) DR[k] = V[k];

			//DumpMatrix(A_Address,A_Entries,cbeg,cend,"mat_original.mtx");
			std::fill(U.begin() + wbeg - mobeg, U.begin() + wend - mobeg, 0.0);
			std::fill(V.begin() + wbeg - mobeg, V.begin() + wend - mobeg, 0.0);
			

			std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, 0);
			std::fill(Llist.begin() + wbeg - mobeg, Llist.begin() + wend - mobeg, 0);

			for (k = wbeg; k < wend; ++k) Blist[localQ[k]] = k;

			for (k = wbeg; k < wend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it != A_Address[k].last; ++it)
				{
					i = A_Entries[it].first;
					u = fabs(A_Entries[it].second);
					//row-col norms
					U[k] += u; //row norm
					V[i] += u; //col norm
					//nonzeros
					Ulist[k]++; //row nnz
					Llist[i]++; //col nnz

					if( Blist[i] == k ) temp[k] = u;
				}
			}

			sort_wgts.clear();
			for (k = wbeg; k < wend; ++k)
                sort_wgts.push_back(wgt_coord(Ulist[k]+Llist[Blist[k]]-1, coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1)*(U[k]+V[Blist[k]]), coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord((Ulist[k]+Llist[Blist[k]]-1)*(U[k]+V[Blist[k]]-temp[k])/temp[k], coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord(1.0/temp[k], coord(k, Blist[k])));
				//sort_wgts.push_back(wgt_coord((U[k]+V[Blist[k]])*(U[k]+V[Blist[k]])/(temp[k]), coord(k, Blist[k])));
			std::sort(sort_wgts.begin(), sort_wgts.end());

			i = wbeg;
			for (wgt_coords::iterator it = sort_wgts.begin(); it != sort_wgts.end(); ++it)
			{
				localP[it->second.first] = i;
				localQ[it->second.second] = i;
				++i;
			}
			cend = i;
#elif defined(REORDER_RCM)
			{
				//create a symmetric graph of the matrix A + A^T
				std::vector<INMOST_DATA_ENUM_TYPE> xadj(wend-wbeg+1), adjncy;
				adjncy.reserve(nzA*2);

				for (k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF ) std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";

				for (k = cbeg; k < cend; ++k) localP[k] = k;
				for (k = cbeg; k < cend; ++k) V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
				for (k = cbeg; k < cend; ++k) DR[k] = V[k];

				for(k = wbeg; k < wend; ++k)
				{
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt != A_Address[k].last; ++jt)
					{
						A_Entries[jt].first = localQ[A_Entries[jt].first];
						//A_Entries[jt].second *= DL[k]*DR[A_Entries[jt].first];
					}
					std::sort(A_Entries.begin()+A_Address[k].first,A_Entries.begin()+A_Address[k].last);
				}

				inversePQ(wbeg,wend,localP,localQ, invP,invQ);
				applyPQ(wbeg, wend, localP, localQ, invP, invQ);

				trcmgraph = Timer();

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
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt != A_Address[i].last; ++jt)
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

					xadj[i-wbeg+1] = static_cast<INMOST_DATA_ENUM_TYPE>(adjncy.size());
				}

				std::fill(Bbeg.begin(),Bbeg.end(),EOL);

				trcmgraph = Timer()-trcmgraph;

				trcmorder = Timer();
				std::fill(Ulist.begin() + wbeg - mobeg, Ulist.begin() + wend - mobeg, ENUMUNDEF);
				//find node with the lowest order
				//INMOST_DATA_ENUM_TYPE start = wbeg;
				INMOST_DATA_ENUM_TYPE index = wbeg;
				INMOST_DATA_ENUM_TYPE cur = ENUMUNDEF;
				std::deque<INMOST_DATA_ENUM_TYPE> q;
				std::vector<INMOST_DATA_ENUM_TYPE> conns;
				//for(k = wbeg+1; k < wend; ++k)
				//    if( RCM_Comparator(wbeg,xadj)(k,start) )
				//        start = k;
				//Ulist[start] = index++;
				//enumerate unknowns with no connections
				for (k = wbeg; k < wend; ++k) if (Ulist[k] == ENUMUNDEF)
				{
					if (xadj[k + 1 - wbeg] - xadj[k - wbeg] == 0)
						Ulist[k] = index++;
				}
				while (index < wend)
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
			if (cbeg == cend && cbeg != wend)
			{
				//std::cout << __FILE__ << ":" << __LINE__ << " singular matrix, factored " << mobeg << ".." << cend << " out of " << mobeg << ".." << moend << std::endl;
				for (k = cbeg; k < wend; ++k)
					LU_Diag[k] = tol_modif;
				break;
			}
#if defined(REORDER_DDPQ)
			if( verbosity > 1 )
			{
				//std::cout << level_size.size() << " new level " << cbeg << ".." << cend << " out of " << wbeg << ".." << wend << std::endl;
				std::cout << std::setw(8) << level_size.size() << " total " << std::setw(8) << wend - wbeg << " after tau filtering " << std::setw(8) << sort_wgts.size() << " selected " << std::setw(8) << cend - wbeg;
				std::cout << std::scientific << " max " << std::setw(8) << ddmaxall << " mean " << std::setw(8) << ddmaxmean << " min " << std::setw(8) << ddmaxmin << " ratio " << ddmaxall/ddmaxmean;
				std::cout << std::endl;
			}
#endif
			tt = Timer();
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

				//U[k] = DL[localP[k]];
				//V[k] = DR[localQ[k]];
			}
			for (k = cbeg; k < cend; ++k) 
			{
				DL[k] = U[k];
				DR[k] = V[k];
			}

			tlreorder = Timer() - tt;
			treorder += tlreorder;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// REASSAMBLE             ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			tt = Timer();
			double tt1, tt2,ttt;
			//in localPQ numbers indicate where to put current row/column
			//reorder E,F blocks by swaps
			tt1 = Timer();
			//inverse ordering
			inversePQ(wbeg,wend,localP,localQ, invP,invQ);
			tt1 = Timer() - tt1;

			//std::cout << "reorder: " << tt1 << std::endl;
			tt2 = Timer();
			nzA = 0;
			for (k = cbeg; k < cend; ++k)
			{
				B_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(B_Entries.size());
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt != A_Address[invP[k]].last; ++jt)
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
			

			tlreassamble = Timer() - tt;
			treassamble += tlreassamble;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING                 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Rescale current B block by Row-Column alternating scaling
			tt = Timer();
#if defined(RESCALE_B)
			if( verbosity > 1 )
				std::cout << " rescaling block B, iters " << sciters << std::endl;

			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];

				U[k] = DL[k];
				V[k] = DR[k];
			}
			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"MC64_MPTILUC2.mtx");
			/*
			{
				std::stringstream name;
				name << "mc64_" << info->GetRank() << ".mtx";
				DumpMatrix(B_Address,B_Entries,cbeg,cend,name.str());
			}
			*/
			/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// COMPUTE GIRSCHGORIN RADIUS ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			INMOST_DATA_REAL_TYPE radii = 0;
			for (k = cbeg; k < cend; k++)
			{
				INMOST_DATA_REAL_TYPE local_radii = 0;
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					if( B_Entries[r].first != k ) local_radii += fabs(B_Entries[r].second);
				if( radii < local_radii ) radii = local_radii;
			}

			*/


			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_mc64.mtx");
#if defined(EQUALIZE_1NORM)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// ROW-COLUMN ALTERNATING SCALING FOR 1-NORM BALANCING////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
					DL[k] += fabs(B_Entries[r].second);//*B_Entries[r].second;
			}
			for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DR[B_Entries[r].first] += DL[k] * fabs(B_Entries[r].second);//*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
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
			std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, 0.0);
			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
					DL[k] += B_Entries[r].second*B_Entries[r].second;
			}
			for (k = cbeg; k < cend; k++) if (DL[k] < eps) DL[k] = 1.0 / subst; else DL[k] = 1.0 / DL[k];
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DR.Begin()+cbeg-mobeg, DR.Begin()+cend-mobeg, 0.0);
				for (k = cbeg; k < cend; k++)
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
						DR[B_Entries[r].first] += DL[k] * B_Entries[r].second*B_Entries[r].second;
				}
				for (k = cbeg; k < cend; k++) if (DR[k] < eps) DR[k] = 1.0 / subst; else DR[k] = 1.0 / DR[k];
				std::fill(DL.Begin()+cbeg-mobeg, DL.Begin()+cend-mobeg, 0.0);
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////THIS VERSION OF RESCALING INCREASES DIAGONAL DOMINANCE///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(k = cbeg; k < cend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r) 
				{
					u = fabs(B_Entries[r].second);
					if( u > 0.0 )
						C_Entries[r] = -log(u);
					else 
						C_Entries[r] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
				}
			}

			
			std::fill(temp.begin() + cbeg - mobeg, temp.begin() + cend - mobeg, 0.0);

			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
				std::fill(DL.Begin() + cbeg - mobeg, DL.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
				std::fill(DR.Begin() + cbeg - mobeg, DR.Begin() + cend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());	
				for (k = cbeg; k < cend; k++) //row number
				{
					for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					{
						i = B_Entries[r].first; //column number
						if( i != k ) //out of diagonal
						{
							u = C_Entries[r] + temp[k] - temp[i];
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
			}

			for (k = cbeg; k < cend; k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					B_Entries[r].second *= DL[k] * DR[B_Entries[r].first];
			}
#endif
			
			/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// COMPUTE GIRSCHGORIN RADIUS ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			radii = 0;
			for (k = cbeg; k < cend; k++)
			{
				INMOST_DATA_REAL_TYPE local_radii = 0;
				for (INMOST_DATA_ENUM_TYPE r = B_Address[k].first; r < B_Address[k].last; ++r)
					if( B_Entries[r].first != k ) local_radii += fabs(B_Entries[r].second);
				if( radii < local_radii ) radii = local_radii;
			}

			*/

#if defined(EQUALIZE_1NORM) || defined(EQUALIZE_2NORM) || defined(EQUALIZE_IDOMINANCE)
			for (k = cbeg; k < cend; k++)
			{
				DL[k] *= U[k];
				DR[k] *= V[k];
			}
#else
			for (k = cbeg; k < cend; k++)
			{
				DL[k] = U[k];
				DR[k] = V[k];
			}
#endif

			//DumpMatrix(B_Address,B_Entries,cbeg,cend,"mat_equilibration.mtx");
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// RESCALING DONE ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//End rescale B block
			trescale += Timer() - tt;
#endif
			


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// FACTORIZATION BEGIN ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			{
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
					Bstart[k] = B_Address[k].first; //remember B block starts permutation for pivoting
					for (i = B_Address[k].first; i < B_Address[k].last; i++)
					{
						if (B_Entries[i].first == k)
						{
#if defined(TRACK_DIAGONAL)
							LU_Diag[k] = B_Entries[i].second;
#endif
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
				//assert(B_Entries[B_Address[cbeg].first].first == cbeg);
				/*
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
					//if( !(l == l) ) std::cout << __FILE__ << ":" << __LINE__ << " nan " << std::endl;
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
#if defined(TRACK_DIAGONAL)
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
						assert(Ui > cbeg);
#if defined(DIAGONAL_PIVOT)
						mean_diag -= fabs(LU_Diag[Ui]);
#endif
						LU_Diag[Ui] -= LU_Entries[i].second * LU_Entries[j].second * LU_Diag[cbeg];
#if defined(DIAGONAL_PIVOT)
						mean_diag += fabs(LU_Diag[Ui]);
#endif
						i++;
						j++;
					}
				}
#endif
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
					tlocal = Timer();
					NuU = NuL = 1.0;
					if( estimator )
					{
						EstU1[cbeg] = EstL1[cbeg] = 0.0;
#if defined(ESTIMATOR_REFINE)
						EstU2[cbeg] = EstL2[cbeg] = 0.0;
#endif
						for (INMOST_DATA_ENUM_TYPE it = U_Address[cbeg].first; it != U_Address[cbeg].last; ++it) 
						{
							EstU1[LU_Entries[it].first] = LU_Entries[it].second;
#if defined(ESTIMATOR_REFINE)
							EstU2[LU_Entries[it].first] = LU_Entries[it].second;
#endif
						}
						for (INMOST_DATA_ENUM_TYPE it = L_Address[cbeg].first; it < L_Address[cbeg].last; ++it) 
						{
							EstL1[LU_Entries[it].first] = LU_Entries[it].second;
#if defined(ESTIMATOR_REFINE)
							EstL2[LU_Entries[it].first] = LU_Entries[it].second;
#endif
						}
					}
					testimator += Timer() - tlocal;
				}
#endif
				nzLU += L_Address[cbeg].Size() + U_Address[cbeg].Size() + 1;
				max_diag = min_diag = fabs(LU_Diag[cbeg]);
				*/
				NuD = 1.0;
				if( verbosity > 1 )
					std::cout << " starting factorization " << std::endl;
				for (k = cbeg; k < cend; ++k)
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
					std::cout << k << "/" << cend << " " << mean_diag << " " << mean_diag2 << std::endl;
					*/
///////////////////////////////////////////////////////////////////////////////////
//              diagonal pivoting algorithm                                      //
///////////////////////////////////////////////////////////////////////////////////
#if defined(DIAGONAL_PIVOT_COND)
					int no_swap_algorithm = 1;
#endif
#if defined(DIAGONAL_PIVOT_TAU)
					if ( k < cend-1 &&  fabs(LU_Diag[k])*(cend-k) < DIAGONAL_PIVOT_TAU*mean_diag)
#else
					if( false )
#endif
					{
#if defined(DIAGONAL_PIVOT_COND)
						if( false )
						{
swap_algorithm:
							no_swap_algorithm--;
						}
#endif
						tlocal = Timer();
						j = k+1;
						for (i = k + 1; i < cend; i++) if (fabs(LU_Diag[j]) < fabs(LU_Diag[i])) j = i;

						if( j > cend-1 )
						{
							std::cout << "Oops! asking to exchange beyond matrix! k = " << k << " j = " << j << " cend = " << cend << std::endl;
						}
						/*
            if( fabs(LU_Diag[j]) < fabs(LU_Diag[k]) )
            {
              std::cout << "Oops! I'm the best pivot! [" << j << "] " << LU_Diag[j] << " [" << k << "] " << LU_Diag[k] << " mean " << mean_diag/static_cast<INMOST_DATA_REAL_TYPE>(cend-k) << "\t\t\t" << std::endl;
              
            }
            else*/
						{
							//TODO!!!
							
							++swaps;
							if( verbosity > 2 )
							{
								std::cout << "Detected, that there is a much better pivot, i'm " << k << " " << LU_Diag[k] << " other " << j << " " << LU_Diag[j] << std::endl;
								std::cout << "Condition numbers: L " << NuL << " D " << NuD << " U " << NuU << std::endl;
							}
							//scanf("%*c");
							//std::cout << "But repivoting algorithm not implemented" << std::endl;
							//This algorithm may be quite costly, but the effect is miraculous
							//First correct E
							//SwapEntries(*E_Address.back(), E_Entries, cend, wend, k, j);
							//Then correct F

							//CheckOrder(B_Address,B_Entries,cbeg,cend);
							//CheckOrder(L_Address,LU_Entries,cbeg,cend);
							//CheckOrder(U_Address,LU_Entries,cbeg,cend);
							//CheckOrder(L2_Address,LU2_Entries,cbeg,cend);
							//CheckOrder(U2_Address,LU2_Entries,cbeg,cend);
							
							//SwapLine(F_Address,k, j);
							//Swap lines of B
							SwapLine(B_Address, k, j);
							//Correct B the same way E is corrected
							SwapEntries(B_Address, B_Entries, cbeg, cend, k, j);
							//now reintialize entrie column linked list
							for (i = cbeg; i < cend; ++i) Bbeg[i] = EOL;
							for (i = cend; i > cbeg; --i)
							{
								if( B_Address[i - 1].Size() > 0 )
								{
									Li = B_Entries[B_Address[i - 1].first].first;
									if (Li < i - 1)
									{
										Blist[i - 1] = Bbeg[Li];
										Bbeg[Li] = i - 1;
									}
								}
							}
							//now update L
							//SwapLine(L_Address, k, j);
							SwapEntries(L_Address, LU_Entries, cbeg, k, k, j);
							//reinitialize linked list
							//now update U
							//SwapLine(U_Address, k, j);
							SwapEntries(U_Address, LU_Entries, cbeg, k, k, j);
							//reinitialize linked list
#if defined(ILUC2)
							//now update L2
							//SwapLine(U2_Address, k, j);
							SwapEntries(U2_Address, LU2_Entries, cbeg, k, k, j);
							//reinitialize linked list
							//now update U2
							//SwapLine(L2_Address, k, j);
							SwapEntries(L2_Address, LU2_Entries, cbeg, k, k, j);
							//reinitialize linked list
#endif
							for (i = cbeg; i < cend; i++)
							{
								Ubeg[i] = EOL;
								Lbeg[i] = EOL;
#if defined(ILUC2)
								U2beg[i] = EOL;
								L2beg[i] = EOL;
#endif
							}
							for (i = cbeg; i < cend; i++)
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
								u = EstU1[k];
								EstU1[k] = EstU1[j];
								EstU1[j] = u;
#if defined(ESTIMATOR_REFINE)
								u = EstU2[k];
								EstU2[k] = EstU2[j];
								EstU2[j] = u;
								u = EstL2[k];
								EstL2[k] = EstL2[j];
								EstL2[j] = u;
#endif
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
							Ui = Bstart[k];
							Bstart[k] = Bstart[j];
							Bstart[j] = Ui;

							//Ui = localQ[k];
							//localQ[k] = localQ[j];
							//localQ[j] = Ui;

							//Ui = localP[k];
							//localP[k] = localP[j];
							//localP[j] = Ui;
							
							
							
							//CheckOrder(B_Address,B_Entries,cbeg,cend);
							//CheckOrder(L_Address,LU_Entries,cbeg,cend);
							//CheckOrder(U_Address,LU_Entries,cbeg,cend);
							//CheckOrder(L2_Address,LU2_Entries,cbeg,cend);
							//CheckOrder(U2_Address,LU2_Entries,cbeg,cend);
						}
						/*
						else
						{
							std::cout << "Detected that there are no good pivots anymore! interval " << cbeg << "..." << cend << " current element " << k << std::endl;
							std::cout << "But level rejection algorithm is not implemented" << std::endl;
							//probably should reject all unfactored parts for the next level
						}
						*/
						tlocal = Timer() - tlocal;
						tswap += tlocal;
					}
					
#endif
///////////////////////////////////////////////////////////////////////////////////
//            Prepare diagonal value                                             //
///////////////////////////////////////////////////////////////////////////////////
#if defined(TRACK_DIAGONAL)
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
					//abs_udiag = mean_diag/static_cast<INMOST_DATA_REAL_TYPE>(cend-k);//fabs(LU_Diag[k]);
					abs_udiag = fabs(LU_Diag[k]);
#else
					abs_udiag = 1.0;
#endif
					(void)abs_udiag;

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
					if (B_Entries[B_Address[k].first].first == k)
						LineValuesU[k] = B_Entries[B_Address[k].first].second;

					Ui = k;
					for (INMOST_DATA_ENUM_TYPE it = B_Address[k].first + (B_Entries[B_Address[k].first].first == k ? 1 : 0); it < B_Address[k].last; ++it)
					{
						LineValuesU[B_Entries[it].first] = B_Entries[it].second;
						Ui = LineIndecesU[Ui] = B_Entries[it].first;
					}
					LineIndecesU[Ui] = EOL;
					//~ Sbeg = k;
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
            //if( fabs(l)*NuU /*NuL_old*/ > tau2 *abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > tau2*tau2 * abs_udiag )
#endif
						{
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  u = l*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= u;
							  else //if (fabs(u)*NuU /*NuL_old*/ > tau2 *abs_udiag )//add new entry
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
								  /*
									LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									*/
							  }
						  }
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with L                        //
///////////////////////////////////////////////////////////////////////////////////
						//if( fabs(l)*NuU /*NuL_old*/ > tau*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(l) > tau2*abs_udiag )
#endif
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								u = l*LU2_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else //if (fabs(u)*NuU /*NuL_old*/ > tau2*abs_udiag)//add new entry
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
									/*
									LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									*/
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
						if( fabs(l) > tau2*abs_udiag )
#endif
						{
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = U_Address[i].first; it < U_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  u = l*LU_Entries[it].second;
							  //eliminate values
							  if (LineIndecesU[j] != UNDEF) //there is an entry in the list
								  LineValuesU[j] -= u;
							  else //if (fabs(u)*NuU /*NuL_old*/ > tau2*abs_udiag)//add new entry
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
								  /*
								  LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									*/
							  }
						  }
///////////////////////////////////////////////////////////////////////////////////
//                 second order U-part elimination with second-order L           //
///////////////////////////////////////////////////////////////////////////////////
#if 0
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = U2_Address[i].first; it < U2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								u = l*LU2_Entries[it].second;
								//eliminate values
								if (LineIndecesU[j] != UNDEF) //there is an entry in the list
									LineValuesU[j] -= u;
								else //if (fabs(u)*NuU /*NuL_old*/ > tau2*abs_udiag)//add new entry
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
									/*
									LineValuesU[j] = -u;
									LineIndecesU[j] = Sbeg;
									Sbeg = j;
									*/
									
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
					/*
					indices.clear();
					Ui = Sbeg;
					while(Ui != EOL)
					{
						indices.push_back(Ui);
						Ui = LineIndecesU[Ui];
					}
					if( !indices.empty() )
					{
						std::sort(indices.begin(),indices.end());
						//Sbeg = indices[0];
						for(size_t qt = 1; qt < indices.size(); ++qt)
							LineIndecesU[indices[qt-1]] = indices[qt];
						LineIndecesU[indices.back()] = EOL;
					}
					*/
					//std::cout << std::endl;
					//define the diagonal value
//this check will fail due to global check of tolerances, if global check will be removed then this will never fail
//#if defined(DIAGONAL_PIVOT) && !defined(DIAGONAL_PERTURBATION)
//					if (fabs(LU_Diag[k] - LineValues[k]) > 1e-6) std::cout << __LINE__ << " Diagonal value went wrong, good: " << LineValues[k] << " have " << LU_Diag[k] << " line " << k << std::endl;
//#endif
///////////////////////////////////////////////////////////////////////////////////
//                  Retrieve diagonal value                                      //
///////////////////////////////////////////////////////////////////////////////////
#if !defined(TRACK_DIAGONAL)
					udiag = LineValuesU[k];
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
					LineIndecesL[k] = EOL;
					LineValuesL[k] = 0.0;
					if( B_Address[k].Size() )
					{
						if (B_Entries[B_Address[k].first].first == k)
							LineValuesL[k] = B_Entries[B_Address[k].first].second;
						Ui = k;
						Li = Bbeg[k];
						while (Li != EOL)
						{
							if( !(B_Entries[B_Address[Li].first].first == k) ) std::cout << "B_Entries[B_Address[Li].first].first: " << B_Entries[B_Address[Li].first].first << " k: " << k << std::endl;
							assert(B_Entries[B_Address[Li].first].first == k);
							LineValuesL[Li] = B_Entries[B_Address[Li].first].second;
							Ui = LineIndecesL[Ui] = Li;
							Li = Blist[Li];
						}
						LineIndecesL[Ui] = EOL;
					}
					//~ Sbeg = k;
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
            //if( fabs(u)*NuL /*NuU_old*/ > tau2*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > tau2*tau2*abs_udiag )
#endif
						{
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it < L_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  l = u*LU_Entries[it].second;
							  //if( !(l == l) ) std::cout << __FILE__ << ":" << __LINE__ << " nan " << std::endl;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
							  else //if (fabs(l)*NuL /*NuU_old*/ > tau2*abs_udiag)//add new entry
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
								  /*
								  LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									*/
							  }
						  }
            }
#if defined(ILUC2)
///////////////////////////////////////////////////////////////////////////////////
//                 second-order L-part elimination with U                        //
///////////////////////////////////////////////////////////////////////////////////
						//if( fabs(u)*NuL /*NuU_old*/ > tau*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > tau2*abs_udiag )
#endif
						{
							curr = cbeg;
							for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
							{
								j = LU2_Entries[it].first;
								l = u*LU2_Entries[it].second;
								//if( !(l == l) ) std::cout << __FILE__ << ":" << __LINE__ << " nan " << std::endl;
								if (LineIndecesL[j] != UNDEF) //there is an entry in the list
									LineValuesL[j] -= l;
								else //if (fabs(l)*NuL /*NuU_old*/ > tau2*abs_udiag)//add new entry
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
									
									/*
									LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									*/
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
            //if( fabs(u)*NuL /*NuU_old*/ > tau2*abs_udiag )//global check
#if defined(PREMATURE_DROPPING)
						if( fabs(u) > tau2*abs_udiag )
#endif
            {
						  curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L_Address[i].first; it < L_Address[i].last; ++it)
						  {
							  j = LU_Entries[it].first;
							  l = u*LU_Entries[it].second;
							  //if( !(l == l) ) std::cout << __FILE__ << ":" << __LINE__ << " nan " << std::endl;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
							  else //if (fabs(l)*NuL /*NuU_old*/ > tau2*abs_udiag)//add new entry
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
								  
								  /*
								  LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									*/
                }
						  }
///////////////////////////////////////////////////////////////////////////////////
//           second-order L-part elimination with second-order U                 //
///////////////////////////////////////////////////////////////////////////////////
#if 0
              curr = cbeg;
						  for (INMOST_DATA_ENUM_TYPE it = L2_Address[i].first; it < L2_Address[i].last; ++it)
						  {
							  j = LU2_Entries[it].first;
							  l = u*LU2_Entries[it].second;
							  //if( !(l == l) ) std::cout << __FILE__ << ":" << __LINE__ << " nan " << std::endl;
							  //eliminate values
							  if (LineIndecesL[j] != UNDEF) //there is an entry in the list
								  LineValuesL[j] -= l;
                else //if (fabs(l)*NuL /*NuU_old*/ > tau2 *abs_udiag)//add new entry
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
								  
								  /*
								  LineValuesL[j] = -l;
									LineIndecesL[j] = Sbeg;
									Sbeg = j;
									*/
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
					/*
					indices.clear();
					Ui = Sbeg;
					while(Ui != EOL)
					{
						indices.push_back(Ui);
						Ui = LineIndecesL[Ui];
					}
					if( !indices.empty() )
					{
						std::sort(indices.begin(),indices.end());
						//Sbeg = indices[0];
						for(size_t qt = 1; qt < indices.size(); ++qt)
							LineIndecesL[indices[qt-1]] = indices[qt];
						LineIndecesL[indices.back()] = EOL;
					}
					*/
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
						//if( !(LineValuesL[Li] == LineValuesL[Li]) ) std::cout << __FILE__ << ":" << __LINE__ << " nan " << std::endl;
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
#if defined(ESTIMATOR_REFINE)
						NuU2_old = NuU2;
						(void)NuU2_old;
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
#if defined(ESTIMATOR_REFINE)
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
#endif
						testimator += Timer()-tlocal;
					}
#endif //ESTIMATOR
///////////////////////////////////////////////////////////////////////////////////
//         discarding current iteration based on condition numbers               //
///////////////////////////////////////////////////////////////////////////////////
#if defined(ESTIMATOR)
          //discard line if condition number is too large
#if defined(DIAGONAL_PIVOT) && defined(DIAGONAL_PIVOT_COND)
#if defined(ESTIMATOR_REFINE)
					NuU_tmp = std::max(NuU1_new,NuU2_new);
					NuL_tmp = std::max(NuL1_new,NuL2_new);
#else
					NuU_tmp = NuU1_new;
					NuL_tmp = NuL1_new;
#endif
					if( k != cend-1 && no_swap_algorithm && (NuU_tmp > DIAGONAL_PIVOT_COND || NuL_tmp > DIAGONAL_PIVOT_COND) )
					{
						if( verbosity > 1 )
							std::cout << "Requested pivoting based on condition estimator (L" << NuL1_new << " " << NuL2_new << " U " << NuU1_new << " " << NuU2_new << ")! row " << k << "/" << cend << std::endl;
						//restore condition number
						NuL1 = NuU1_old;
						NuU1 = NuU1_old;
#if defined(ESTIMATOR_REFINE)
						NuL2 = NuU2_old;
						NuU2 = NuU2_old;
#endif
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
					Ui = LineIndecesU[k];
					while (Ui != EOL)
					{
						u = fabs(LineValuesU[Ui]);
						if (u*NuU /*NuL_old*/ > tau) // apply dropping rule
							LU_Entries.push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#if defined(ILUC2)
						else if (u*NuU /*NuL_old*/ > tau2)
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
						//if( !(u == u) ) std::cout << __FILE__ << ":" << __LINE__ << " nan " << std::endl;
						if (u*NuL /*NuU_old*/ > tau) //apply dropping
							LU_Entries.push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#if defined(ILUC2)
						else if (u*NuL /*NuU_old*/ > tau2)
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
					LU_Diag[k] = udiag;
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
#if defined(TRACK_DIAGONAL)
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
#if 0
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
#endif
#endif //ILUC2
#endif //TRACK_DIAGONAL
//#endif //DIAGONAL_PIVOT

					//CheckOrder(L_Address,LU_Entries,cbeg,k+1);
					//CheckOrder(U_Address,LU_Entries,cbeg,k+1);
					//CheckOrder(L2_Address,LU2_Entries,cbeg,k+1);
					//CheckOrder(U2_Address,LU2_Entries,cbeg,k+1);

					//number of nonzeros
					nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
					

					if( verbosity )
					{
						if( k % 100 == 0 )
						{
							std::ios save(NULL);
							save.copyfmt(std::cout);
							if (verbosity == 1)
								std::cout << cend << "/" << moend << " factor " << std::setw(6) << std::fixed << std::setprecision(2) << 100.0f*(k - cbeg) / (float)(cend - cbeg) << "\r" << std::flush;
							else 
							{
								std::cout << std::setw(6) << std::fixed << std::setprecision(2) << 100.0f*(k - cbeg) / (float)(cend - cbeg);
								std::cout << " nnz LU " << std::setw(8) << nzLU;
#if defined(ILUC2)
								std::cout << " LU2 " << std::setw(8) << nzLU2;
#endif
								std::cout << " cond L " << std::setprecision(6) << std::setw(10) << NuL << " D " << NuD << " U " << NuD << " pivot " << std::setw(10) << swaps;
								std::cout << "\r" << std::flush;
							}
							std::cout.copyfmt(save);
						}
					}
///////////////////////////////////////////////////////////////////////////////////
//                       iteration done                                          //
///////////////////////////////////////////////////////////////////////////////////
				}
				if( verbosity > 1 )
				{
					std::cout << "size " << moend-mobeg;
					std::cout << " total nonzeros in A " << nzA << " in LU " << nzLU;
#if defined(ILUC2)
					std::cout << " in LU2 " << nzLU2;
#endif
					std::cout << " conditions L " << NuL_max << " D " << NuD << " U " << NuU_max/* << " pivot swaps " << swaps*/ << std::endl;
				}
			}
			tlfactor = Timer() - tt;
			tfactor += tlfactor;
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
			
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// FACTORIZATION COMPLETE ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//After we have factored the rescaled matrix we must rescale obtained factors
			tt = Timer();
#if defined(RESCALE_B)
			if( verbosity > 1 )
				std::cout << " rescaling block B back " << std::endl;

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
			
			wend = wbeg;
		}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////  UPDATE ORDERING DUE TO PIVOTING //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//inversePQ(wbeg,wend,localP,localQ,invP,invQ);
		//applyPQ(wbeg,wend,localP,localQ,invP,invQ);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////  FACTORIZATION COMPLETE ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		ttotal = Timer() - ttotal;
		if( verbosity > 2 )
		{
			std::cout << "total      " << ttotal       << "\n";
			std::cout << "reorder    " << treorder     << " ("; fmt(std::cout,100.0*treorder/ttotal    ,6,2) << "%)\n";
			std::cout << "   mpt     " << ttransversal << " ("; fmt(std::cout,100.0*ttransversal/ttotal,6,2) << "%)\n";
#if defined(REORDER_METIS_ND)
			std::cout << "   graph   " << tmetisgraph  << " ("; fmt(std::cout,100.0*tmetisgraph/ttotal ,6,2) << "%)\n";
			std::cout << "   nd      " << tmetisnd     << " ("; fmt(std::cout,100.0*tmetisnd/ttotal    ,6,2) << "%)\n";
#endif
#if defined(REORDER_RCM)
			std::cout << "   graph   " << trcmgraph    << " ("; fmt(std::cout,100.0*trcmgraph/ttotal   ,6,2) << "%)\n";
			std::cout << "   rcm     " << trcmorder    << " ("; fmt(std::cout,100.0*trcmorder/ttotal   ,6,2) << "%)\n";
#endif
			std::cout << "reassamble " << treassamble  << " ("; fmt(std::cout,100.0*treassamble/ttotal ,6,2) << "%)\n";
			std::cout << "rescale    " << trescale     << " ("; fmt(std::cout,100.0*trescale/ttotal    ,6,2) << "%)\n";
			std::cout << "factor     " << tfactor      << " ("; fmt(std::cout,100.0*tfactor/ttotal     ,6,2) << "%)\n";
			std::cout << "   swap    " << tswap        << " ("; fmt(std::cout,100.0*tswap/ttotal       ,6,2) << "%)\n";
			std::cout << "   cond    " << testimator   << " ("; fmt(std::cout,100.0*testimator/ttotal  ,6,2) << "%)\n";
#if defined(ILUC2)
			std::cout << "nnz A " << nzA << " LU " << nzLU << " LU2 " << nzLU2 << " swaps " << swaps << "\n";
#else
			std::cout << "nnz A " << nzA << " LU " << nzLU << " swaps " << swaps << "\n";
#endif
		}
		init = true;
		/*
		level_interval.resize(level_size.size());
		level_interval[0].first = mobeg;
		level_interval[0].last = mobeg + level_size[0];
		for (k = 1; k < level_size.size(); k++)
		{
			level_interval[k].first = level_interval[k - 1].last;
			level_interval[k].last = level_interval[k].first + level_size[k];
		}
		*/
		condestL = NuL;
		condestU = NuU;
		//for partition of unity
		/*
		info->PrepareVector(div);
		std::fill(div.Begin(),div.End(),0);
		for(k = mobeg; k < moend; k++) div[k] = 1.0;
		info->Accumulate(div);
		for(k = mobeg; k < moend; k++) div[k] = 1.0/div[k];
		 */
		return true;
	}
	bool MTILUC_preconditioner::Finalize()
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
		return true;
	}
	void MTILUC_preconditioner::ApplyB(double alpha, Sparse::Vector & x, double beta, Sparse::Vector & y) // y = alpha A x + beta y
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
	bool MTILUC_preconditioner::Solve(Sparse::Vector & input, Sparse::Vector & output)
	{
		assert(&input != &output);
		//
		INMOST_DATA_ENUM_TYPE mobeg, moend, k;
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			INMOST_DATA_ENUM_TYPE vbeg, vend;
			info->GetOverlapRegion(info->GetRank(), mobeg, moend);
			info->GetVectorRegion(vbeg, vend);
		
			for (k = vbeg; k < mobeg; k++) temp[k] = 0;
			for (k = mobeg; k < moend; ++k) temp[k] = input[ddP[k]];
			for (k = moend; k < vend; k++) temp[k] = 0;

			//for (k = mobeg; k < moend; ++k) if( temp[k] != temp[k] ) {std::cout << "Nan in input" << " rank " << info->GetRank() << std::endl; break;}

			//Solve with L first
			for (k = mobeg; k < moend; ++k) if( L_Address[k].first != L_Address[k].last ) //iterator over columns of L
			{
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r)
				{
					temp[LU_Entries[r].first] -= temp[k] * LU_Entries[r].second;
					//if( LU_Entries[r].second != LU_Entries[r].second ) std::cout << "L-row " << k << " column " << LU_Entries[r].first << " is nan " << " rank " << info->GetRank() << std::endl;
				}
			}
			
			//for (k = mobeg; k < moend; ++k) if( temp[k] != temp[k] ) {std::cout << "Nan after L-solve" << " rank " << info->GetRank() << std::endl; break;}
			//Solve with diagonal
			for (k = mobeg; k < moend; ++k) 
			{
				temp[k] /= LU_Diag[k];
				//if( LU_Diag[k] != LU_Diag[k] ) std::cout << "Diag-row " << k << " is nan " << std::endl;
			}
			
			//for (k = mobeg; k < moend; ++k) if( temp[k] != temp[k] ) {std::cout << "Nan after diagonal" << " rank " << info->GetRank() << std::endl; break;}
			//Solve with U
			for (k = moend; k > mobeg; --k) if( U_Address[k - 1].first != U_Address[k - 1].last ) //iterator over rows of U
			{
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first; r < U_Address[k - 1].last; ++r)
				{
					temp[k - 1] -= temp[LU_Entries[r].first] * LU_Entries[r].second;
					//if( LU_Entries[r].second != LU_Entries[r].second ) std::cout << "U-row " << k << " at " << LU_Entries[r].first << " is nan " << " rank " << info->GetRank() << std::endl;
				}
			}
			
			//for (k = mobeg; k < moend; ++k) if( temp[k] != temp[k] ) {std::cout << "Nan after U-solve" << " rank " << info->GetRank() << std::endl; break;}
			
			//for (k = mobeg; k < moend; ++k) if( temp[k] != temp[k] ) {std::cout << "Nan in output" << " rank " << info->GetRank() << std::endl; break;}

		
			for (k = mobeg; k < moend; ++k) output[ddQ[k]] = temp[k];

			//Restrict additive schwartz (maybe do it outside?)
			//May assamble partition of unity instead of restriction before accumulation
			//assembly should be done instead of initialization
			for (k = vbeg; k < mobeg; ++k) output[k] = 0;
			for (k = moend; k < vend; ++k) output[k] = 0;
		}
		info->Accumulate(output);
		/*
#if defined(USE_OMP)
#pragma omp single
#endif
		{
			for (k = mobeg; k < moend; ++k) output[k] *= div[k];
		}
		 */
		return true;
	}
	bool MTILUC_preconditioner::ReplaceMAT(Sparse::Matrix & A){ if (isInitialized()) Finalize(); Alink = &A; return true; }
	bool MTILUC_preconditioner::ReplaceSOL(Sparse::Vector & x) {(void)x; return true; }
	bool MTILUC_preconditioner::ReplaceRHS(Sparse::Vector & b) {(void)b; return true; }
	Method * MTILUC_preconditioner::Duplicate() { return new MTILUC_preconditioner(*this); }
	MTILUC_preconditioner::~MTILUC_preconditioner()
	{
		if (!isFinalized()) Finalize();
	}
#endif

