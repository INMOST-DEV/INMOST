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

#if defined(USE_OMP)
#define USE_OMP_FACT
#endif //USE_OMP

//~ #undef USE_OMP


//#define REORDER_RCM
#define REORDER_WRCM
//~ #define REORDER_BRCM
#if defined(USE_SOLVER_METIS)
//#define REORDER_METIS_ND
#endif
#if defined(USE_SOLVER_MONDRIAAN)
//#define REORDER_MONDRIAAN
#endif

static int run_nd = 0; //1 - unsymmetric dissection before mpt, 2 - symmetric dissection after mpt
static bool run_mpt = true;
static bool reorder_b = true;
static bool rescale_b = true;
static bool allow_pivot = true;
static bool print_mem = false;
static bool show_summary = false;
const INMOST_DATA_ENUM_TYPE UNDEF = ENUMUNDEF, EOL = ENUMUNDEF - 1;


//~ #define EQUALIZE_1NORM
//#define EQUALIZE_2NORM
#define EQUALIZE_IDOMINANCE

//#define PREMATURE_DROPPING
//#define PIVOT_THRESHOLD
//#define DIAGONAL_PERTURBATION
const double rpert = 1.0e-6;
const double apert = 1.0e-8;
#define ILUC2
#define ILUC2_SCHUR

//#define PIVOT_COND_DEFAULT 0.1/tau
#define PIVOT_COND_DEFAULT 1.0e+2
#define PIVOT_DIAG_DEFAULT 1.0e+5
//#define SCHUR_DROPPING_LF
//#define SCHUR_DROPPING_EU
#define SCHUR_DROPPING_S
#define SCHUR_DROPPING_LF_PREMATURE
#define SCHUR_DROPPING_EU_PREMATURE
// #define SCHUR_DROPPING_S_PREMATURE
#define CONDITION_PIVOT

#if defined(USE_SOLVER_METIS)
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


	void MLMTILUC_preconditioner::ReorderEF(INMOST_DATA_ENUM_TYPE wbeg,
											INMOST_DATA_ENUM_TYPE wend,
											interval<INMOST_DATA_ENUM_TYPE, bool> & donePQ,
											const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
											const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
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
							for(l = 0; l < E_Address.size(); ++l)
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
							for(l = 0; l < F_Address.size(); ++l)
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

	


	void MLMTILUC_preconditioner::DumpMatrix(const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
											 const std::vector< std::vector<Sparse::Row::entry> >& Entries,
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
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it)
			{
				norm1 += fabs(Entries[Athr][it].second);
				norm2 += Entries[Athr][it].second*Entries[Athr][it].second;
				if( fabs(Entries[Athr][it].second) > vrowmax )
				{
					vrowmax = fabs(Entries[Athr][it].second);
					irowmax = Entries[Athr][it].first;
				}

				if( fabs(Entries[Athr][it].second) > vcolmax[Entries[Athr][it].first] )
				{
					vcolmax[Entries[Athr][it].first] = fabs(Entries[Athr][it].second);
					icolmax[Entries[Athr][it].first] = k;
				}
				if( Entries[Athr][it].second > max ) max = Entries[Athr][it].second;
				if( Entries[Athr][it].second < min ) min = Entries[Athr][it].second;
				if( fabs(Entries[Athr][it].second) < fabs(minabs) ) minabs = Entries[Athr][it].second;

				if( Entries[Athr][it].first == k )
				{
					diag_found = true;
					diag = Entries[Athr][it].second;
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
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it)
				fout << (k-wmbeg+1) << " " << (Entries[Athr][it].first-wmbeg+1) << " " << Entries[Athr][it].second << std::endl;
		}
		fout.close();
	}
	INMOST_DATA_ENUM_TYPE MLMTILUC_preconditioner::ComputeNonzeroes(const Block& b, const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address, const std::vector<std::vector<Sparse::Row::entry> >& Entries)
	{
		INMOST_DATA_ENUM_TYPE rbeg = b.row_start, rend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		INMOST_DATA_ENUM_TYPE ret = 0;
		for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
		{
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) if (cbeg <= Entries[Athr][it].first && Entries[Athr][it].first < cend)
				ret++;
		}
		return ret;
	}
	INMOST_DATA_ENUM_TYPE MLMTILUC_preconditioner::ComputeNonzeroes(const Block& b, const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address, const std::vector< std::vector<Sparse::Row::entry> >& Entries, const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ)
	{
		INMOST_DATA_ENUM_TYPE rbeg = b.row_start, rend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		INMOST_DATA_ENUM_TYPE ret = 0;
		for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
		{
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) if (cbeg <= localQ[Entries[Athr][it].first] && localQ[Entries[Athr][it].first] < cend)
				ret++;
		}
		return ret;
	}
	void MLMTILUC_preconditioner::DumpMatrixBlock(const Block & b, 
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		std::string file_name)
	{
		INMOST_DATA_ENUM_TYPE rbeg = b.row_start, rend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		INMOST_DATA_REAL_TYPE norm1 = 0, norm2 = 0, max = -1.0e54, min = 1.0e54, minabs = 1.0e54;
		INMOST_DATA_REAL_TYPE vrowmax, diag, mindiag = 1.0e54, maxdiag = -1.0e54, maxabsdiag = -1.0e54, minabsdiag = 1.0e54;
		INMOST_DATA_ENUM_TYPE nnz = 0, dominant_rows = 0, dominant_cols = 0, irowmax = 0, nodiag = 0;
		INMOST_DATA_ENUM_TYPE addrbeg = Address.get_interval_beg();
		INMOST_DATA_ENUM_TYPE addrend = Address.get_interval_end();
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> vcolmax(cbeg, cend, 0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> icolmax(cbeg, cend, ENUMUNDEF);
		for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
		{
			if (k < addrbeg || k >= addrend) continue;
			vrowmax = 0;

			bool diag_found = false;
			diag = 0;
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) if( cbeg <= Entries[Athr][it].first && Entries[Athr][it].first < cend )
			{
				nnz++;
				norm1 += fabs(Entries[Athr][it].second);
				norm2 += Entries[Athr][it].second * Entries[Athr][it].second;
				if (fabs(Entries[Athr][it].second) > vrowmax)
				{
					vrowmax = fabs(Entries[Athr][it].second);
					irowmax = Entries[Athr][it].first;
				}

				if (fabs(Entries[Athr][it].second) > vcolmax[Entries[Athr][it].first])
				{
					vcolmax[Entries[Athr][it].first] = fabs(Entries[Athr][it].second);
					icolmax[Entries[Athr][it].first] = k;
				}
				if (Entries[Athr][it].second > max) max = Entries[Athr][it].second;
				if (Entries[Athr][it].second < min) min = Entries[Athr][it].second;
				if (fabs(Entries[Athr][it].second) < fabs(minabs)) minabs = Entries[Athr][it].second;

				if (Entries[Athr][it].first == k)
				{
					diag_found = true;
					diag = Entries[Athr][it].second;
				}
			}

			if (diag_found)
			{
				if (mindiag > diag) mindiag = diag;
				if (maxdiag < diag) maxdiag = diag;
				if (minabsdiag > fabs(diag)) minabsdiag = fabs(diag);
				if (maxabsdiag < fabs(diag)) maxabsdiag = fabs(diag);
			}
			else nodiag++;

			if (irowmax == k) ++dominant_rows;
		}

		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) if (icolmax[k] == k) ++dominant_cols;

		std::cout << "Writing matrix to " << file_name.c_str() << std::endl;
		std::fstream fout(file_name.c_str(), std::ios::out);
		fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
		fout << "% maximum " << max << std::endl;
		fout << "% minimum " << min << std::endl;
		fout << "% absolute minimum " << minabs << std::endl;
		fout << "% A 1-norm  " << norm1 << std::endl;
		fout << "% A 2-norm  " << sqrt(norm2) << std::endl;
		fout << "% mean 1-norm  " << norm1 / (rend - rbeg) << std::endl;
		fout << "% mean 2-norm  " << sqrt(norm2 / (rend - rbeg)) << std::endl;
		fout << "% dominant rows  " << dominant_rows << std::endl;
		fout << "% dominant cols  " << dominant_cols << std::endl;
		fout << "% maximal diagonal value " << maxdiag << std::endl;
		fout << "% minimal diagonal value " << mindiag << std::endl;
		fout << "% absolute maximal diagonal value " << maxabsdiag << std::endl;
		fout << "% absolute minimal diagonal value " << minabsdiag << std::endl;
		fout << "% no diagonal value  " << nodiag << std::endl;
		fout << "% true matrix row indices interval " << rbeg << ":" << rend << std::endl;
		fout << "% true matrix col indices interval " << cbeg << ":" << cend << std::endl;
		fout << std::scientific;

		//fout.close(); return;

		fout << rend - rbeg << " " << cend - cbeg << " " << nnz << std::endl;;
		for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
		{
			if (k < addrbeg || k >= addrend) continue;
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) if (cbeg <= Entries[Athr][it].first && Entries[Athr][it].first < cend)
				fout << (k - rbeg + 1) << " " << (Entries[Athr][it].first - cbeg + 1) << " " << Entries[Athr][it].second << std::endl;
		}
		fout.close();
	}
	void MLMTILUC_preconditioner::DumpMatrixBlock(const Block& b, 
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector<std::vector<Sparse::Row::entry> >& Entries,
		const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
		std::string file_name)
	{
		INMOST_DATA_ENUM_TYPE rbeg = b.row_start, rend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		INMOST_DATA_REAL_TYPE norm1 = 0, norm2 = 0, max = -1.0e54, min = 1.0e54, minabs = 1.0e54;
		INMOST_DATA_REAL_TYPE vrowmax, diag, mindiag = 1.0e54, maxdiag = -1.0e54, maxabsdiag = -1.0e54, minabsdiag = 1.0e54;
		INMOST_DATA_ENUM_TYPE nnz = 0, dominant_rows = 0, dominant_cols = 0, irowmax = 0, nodiag = 0;
		INMOST_DATA_ENUM_TYPE addrbeg = Address.get_interval_beg();
		INMOST_DATA_ENUM_TYPE addrend = Address.get_interval_end();
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> vcolmax(cbeg, cend, 0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> icolmax(cbeg, cend, ENUMUNDEF);
		for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
		{
			if (k < addrbeg || k >= addrend) continue;
			vrowmax = 0;

			bool diag_found = false;
			diag = 0;
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) if (cbeg <= localQ[Entries[Athr][it].first] && localQ[Entries[Athr][it].first] < cend)
			{
				nnz++;
				norm1 += fabs(Entries[Athr][it].second);
				norm2 += Entries[Athr][it].second * Entries[Athr][it].second;
				if (fabs(Entries[Athr][it].second) > vrowmax)
				{
					vrowmax = fabs(Entries[Athr][it].second);
					irowmax = localQ[Entries[Athr][it].first];
				}

				if (fabs(Entries[Athr][it].second) > vcolmax[localQ[Entries[Athr][it].first]])
				{
					vcolmax[localQ[Entries[Athr][it].first]] = fabs(Entries[Athr][it].second);
					icolmax[localQ[Entries[Athr][it].first]] = k;
				}
				if (Entries[Athr][it].second > max) max = Entries[Athr][it].second;
				if (Entries[Athr][it].second < min) min = Entries[Athr][it].second;
				if (fabs(Entries[Athr][it].second) < fabs(minabs)) minabs = Entries[Athr][it].second;

				if (localQ[Entries[Athr][it].first] == k)
				{
					diag_found = true;
					diag = Entries[Athr][it].second;
				}
			}

			if (diag_found)
			{
				if (mindiag > diag) mindiag = diag;
				if (maxdiag < diag) maxdiag = diag;
				if (minabsdiag > fabs(diag)) minabsdiag = fabs(diag);
				if (maxabsdiag < fabs(diag)) maxabsdiag = fabs(diag);
			}
			else nodiag++;

			if (irowmax == k) ++dominant_rows;
		}

		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) if (icolmax[k] == k) ++dominant_cols;

		std::cout << "Writing matrix to " << file_name.c_str() << std::endl;
		std::fstream fout(file_name.c_str(), std::ios::out);
		fout << "%%MatrixMarket matrix coordinate real general" << std::endl;
		fout << "% maximum " << max << std::endl;
		fout << "% minimum " << min << std::endl;
		fout << "% absolute minimum " << minabs << std::endl;
		fout << "% A 1-norm  " << norm1 << std::endl;
		fout << "% A 2-norm  " << sqrt(norm2) << std::endl;
		fout << "% mean 1-norm  " << norm1 / (rend - rbeg) << std::endl;
		fout << "% mean 2-norm  " << sqrt(norm2 / (rend - rbeg)) << std::endl;
		fout << "% dominant rows  " << dominant_rows << std::endl;
		fout << "% dominant cols  " << dominant_cols << std::endl;
		fout << "% maximal diagonal value " << maxdiag << std::endl;
		fout << "% minimal diagonal value " << mindiag << std::endl;
		fout << "% absolute maximal diagonal value " << maxabsdiag << std::endl;
		fout << "% absolute minimal diagonal value " << minabsdiag << std::endl;
		fout << "% no diagonal value  " << nodiag << std::endl;
		fout << "% true matrix row indices interval " << rbeg << ":" << rend << std::endl;
		fout << "% true matrix col indices interval " << cbeg << ":" << cend << std::endl;
		fout << std::scientific;

		//fout.close(); return;

		fout << rend - rbeg << " " << cend - cbeg << " " << nnz << std::endl;;
		for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
		{
			if (k < addrbeg || k >= addrend) continue;
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) //if (cbeg <= localQ[Entries[Athr][it].first] && localQ[Entries[Athr][it].first] < cend)
			{
				fout << (k - rbeg + 1);
				fout << " original " << Entries[Athr][it].first;
				fout << " new " << (localQ[Entries[Athr][it].first] - cbeg + 1);
				fout << " " << Entries[Athr][it].second;
				if (localQ[Entries[Athr][it].first] - cbeg == k - rbeg)
					fout << " diag ";
				fout << std::endl;
			}
		}
		fout.close();
	}
	void MLMTILUC_preconditioner::CheckOrder(const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
											 const std::vector< std::vector<Sparse::Row::entry> > & Entries,
											 INMOST_DATA_ENUM_TYPE rbeg, INMOST_DATA_ENUM_TYPE rend)
	{
		INMOST_DATA_ENUM_TYPE i,r;
		for (i = rbeg; i < rend; i++) if( Address[i].first < Address[i].last )
		{
			//check ordered on entry
			bool ordered = true;
			for (r = Address[i].first; r < Address[i].last-1; r++)
			{
				if( !(Entries[Address[i].thr][r].first < Entries[Address[i].thr][r+1].first) )
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
					std::cout << "(" << Entries[Address[i].thr][r].first << "," << Entries[Address[i].thr][r].second << ") ";
				}
				std::cout << std::endl;
				throw -1;
			}
			bool nan = false;
			for (r = Address[i].first; r < Address[i].last; r++)
			{
				if( Entries[Address[i].thr][r].second != Entries[Address[i].thr][r].second )
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
					std::cout << "(" << Entries[Address[i].thr][r].first << "," << Entries[Address[i].thr][r].second << ") ";
				}
				std::cout << std::endl;
				throw -1;
			}
		}
	}
	
	void MLMTILUC_preconditioner::inversePQ(INMOST_DATA_ENUM_TYPE wbeg,
											INMOST_DATA_ENUM_TYPE wend,
											const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
											const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
	{
		//inverse reordering
		// in invPQ numbers indicate where to get current column
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invP[k] = ENUMUNDEF;
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			assert(invP[localP[k]] == ENUMUNDEF);
			invP[localP[k]] = k;
		}
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invQ[k] = ENUMUNDEF;
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			assert(invQ[localQ[k]] == ENUMUNDEF);
			invQ[localQ[k]] = k;
		}
	}
	void MLMTILUC_preconditioner::applyPQ(INMOST_DATA_ENUM_TYPE wbeg,
										  INMOST_DATA_ENUM_TYPE wend,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
										  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
										  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
										  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ)
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
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::AddListOrdered(INMOST_DATA_ENUM_TYPE cbeg, 
																  const Interval & Address,
																  const std::vector<Sparse::Row::entry> & Entries,
																  INMOST_DATA_REAL_TYPE coef,
																  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
																  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
																  INMOST_DATA_REAL_TYPE droptol)
	{
		INMOST_DATA_ENUM_TYPE curr = cbeg, next;
		INMOST_DATA_REAL_TYPE drop_sum = 0;
		for (INMOST_DATA_ENUM_TYPE it = Address.first; it < Address.last; ++it)
		{
			//~ assert(!Pivot[i]); // pivoted rows should be empty
			//~ assert(j >= k); // all indices are supposed to be to the right
			INMOST_DATA_ENUM_TYPE j = Entries[it].first;
			INMOST_DATA_REAL_TYPE u = coef*Entries[it].second;
			//eliminate values
			if (LineIndeces[j] != UNDEF) //there is an entry in the list
				LineValues[j] += u;
			else if( fabs(u) > droptol )
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
				LineValues[j] = u; 
			}
			else drop_sum += u;
		}
		return drop_sum;
	}
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::AddListUnordered(INMOST_DATA_ENUM_TYPE & Sbeg,
																	const Interval & Address,
																	const std::vector<Sparse::Row::entry> & Entries,
																	INMOST_DATA_REAL_TYPE coef,
																	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
																	interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
																	INMOST_DATA_REAL_TYPE droptol)
	{
		INMOST_DATA_REAL_TYPE drop_sum = 0;
		for (INMOST_DATA_ENUM_TYPE it = Address.first; it < Address.last; ++it)
		{
			//~ assert(!Pivot[i]); // pivoted rows should be empty
			//~ assert(j >= k); // all indices are supposed to be to the right
			INMOST_DATA_ENUM_TYPE j = Entries[it].first;
			INMOST_DATA_REAL_TYPE u = coef*Entries[it].second;
			//eliminate values
			if (LineIndeces[j] != UNDEF) //there is an entry in the list
				LineValues[j] += u;
			else if( fabs(u) > droptol )
			{
				LineValues[j] = u;
				LineIndeces[j] = Sbeg;
				Sbeg = j;
			}
			else drop_sum += u;
		}
		return drop_sum;
	}
	void MLMTILUC_preconditioner::OrderList(INMOST_DATA_ENUM_TYPE & Sbeg,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
											std::vector<INMOST_DATA_ENUM_TYPE> & indices)
	{
		indices.clear();
		INMOST_DATA_ENUM_TYPE i = Sbeg;
		while(i != EOL)
		{
			indices.push_back(i);
			i = LineIndeces[i];
		}
		if( !indices.empty() )
		{
			std::sort(indices.begin(),indices.end());
			Sbeg = indices[0];
			for(size_t qt = 1; qt < indices.size(); ++qt)
				LineIndeces[indices[qt-1]] = indices[qt];
			LineIndeces[indices.back()] = EOL;
		}
	}
	
	void MLMTILUC_preconditioner::ScaleList(INMOST_DATA_REAL_TYPE coef,
											INMOST_DATA_ENUM_TYPE beg,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues)
	{										
		INMOST_DATA_ENUM_TYPE i = beg;
		while (i != EOL)
		{
			LineValues[i] *= coef;
			i = LineIndeces[i];
		}
	}
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::Estimator1(INMOST_DATA_ENUM_TYPE k,
															  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
															  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
															  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
															  INMOST_DATA_REAL_TYPE & mu_update)
	{
		INMOST_DATA_REAL_TYPE mup = - Est[k] + 1.0;
		INMOST_DATA_REAL_TYPE mum = - Est[k] - 1.0;
		INMOST_DATA_REAL_TYPE smup = 0, smum = 0;
		INMOST_DATA_ENUM_TYPE i = LineIndeces[k];
		while (i != EOL)
		{
			if (i > k)
			{
				smup += fabs(Est[i] + LineValues[i] * mup);
				smum += fabs(Est[i] + LineValues[i] * mum);
			}
			i = LineIndeces[i];
		}
		if (fabs(mup) + smup > fabs(mum) + smum) 
			mu_update = mup; 
		else 
			mu_update = mum;
		return std::max(fabs(mup),fabs(mum));
	}
	
	
	INMOST_DATA_REAL_TYPE MLMTILUC_preconditioner::Estimator2(INMOST_DATA_ENUM_TYPE k,
															  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
															  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
															  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
															  INMOST_DATA_REAL_TYPE & mu_update)
	{
		INMOST_DATA_REAL_TYPE mup = - Est[k] + 1.0;
		INMOST_DATA_REAL_TYPE mum = - Est[k] - 1.0;
		INMOST_DATA_ENUM_TYPE np = 0, nm = 0;
		//start from the element next after diagonal position
		INMOST_DATA_ENUM_TYPE i = LineIndeces[k];
		while (i != EOL)
		{
			INMOST_DATA_REAL_TYPE v = Est[i];
			INMOST_DATA_REAL_TYPE vp = fabs(v + LineValues[i] * mup);
			INMOST_DATA_REAL_TYPE vm = fabs(v + LineValues[i] * mum);
			v = fabs(v);
			if (vp > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) np++;
			if (vm > std::max<INMOST_DATA_REAL_TYPE>(2 * v, 0.5)) nm++;
			if (std::max<INMOST_DATA_REAL_TYPE>(2 * vp, 0.5) < v) np--;
			if (std::max<INMOST_DATA_REAL_TYPE>(2 * vm, 0.5) < v) nm--;
			i = LineIndeces[i];
		}
		if (fabs(mup) + np > fabs(mum)+nm) 
			mu_update = mup; 
		else 
			mu_update = mum;
		return std::max(fabs(mup), fabs(mum));
	}
	
	void MLMTILUC_preconditioner::EstimatorUpdate(INMOST_DATA_ENUM_TYPE k,
												  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
												  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues,
												  interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Est,
												  INMOST_DATA_REAL_TYPE & mu_update)
	{
		INMOST_DATA_ENUM_TYPE i = LineIndeces[k];
		while (i != EOL)
		{
			Est[i] += LineValues[i] * mu_update;
			i = LineIndeces[i];
		}
	}
	
	void MLMTILUC_preconditioner::DiagonalUpdate(INMOST_DATA_ENUM_TYPE k,
												 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & Diag,
												 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndecesL,
												 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesL,
												 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndecesU,
												 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValuesU)
	{
		INMOST_DATA_ENUM_TYPE Li = LineIndecesL[k];
		INMOST_DATA_ENUM_TYPE Ui = LineIndecesU[k];
		while (Li != EOL && Ui != EOL)
		{
			if( Ui > Li ) Li = LineIndecesL[Li];
			else if( Li > Ui ) Ui = LineIndecesU[Ui];
			else
			{
				assert(Ui > k && Li > k && Ui == Li);
				//~ std::cout << "Diag[" << Ui << "] " << Diag[Ui] << " - " << LineValuesL[Li]*LineValuesU[Ui]*Diag[k] << " = ";
				Diag[Ui] -= LineValuesL[Li]*LineValuesU[Ui]*Diag[k];
				//~ std::cout << Diag[Ui] << " from " << Li << "," << Ui << std::endl;
				Li = LineIndecesL[Li];
				Ui = LineIndecesU[Ui];
			}
		}
	}
	
	void MLMTILUC_preconditioner::ClearList(INMOST_DATA_ENUM_TYPE beg,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & LineIndeces,
											interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & LineValues)
	{
		INMOST_DATA_ENUM_TYPE i = beg, j;
		while (i != EOL)
		{
			j = LineIndeces[i];
			LineValues[i] = 0.0; //clean values after use
			LineIndeces[i] = UNDEF; //clean indeces after use
			i = j;
		}
	}
	
	void MLMTILUC_preconditioner::PrepareTranspose(INMOST_DATA_ENUM_TYPE cbeg,
												   INMOST_DATA_ENUM_TYPE cend,
												   const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
												   const std::vector< std::vector<Sparse::Row::entry> > & Entries,
												   interval<INMOST_DATA_ENUM_TYPE, Interval> G_Address,
												   std::vector< std::vector< std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > > G_Entries)
												   //interval<INMOST_DATA_ENUM_TYPE, std::vector< std::pair<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > > & Indices )
	{
#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
		{
			int thr = Thread();
			INMOST_DATA_ENUM_TYPE tseg = (cend - cbeg) / Threads();
			INMOST_DATA_ENUM_TYPE tbeg = cbeg + tseg * thr;
			INMOST_DATA_ENUM_TYPE tend = tbeg + tseg;
			if (thr == Threads() - 1) tend = cend;
			if (tend > tbeg)
			{
				for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
				{
					G_Address[k].thr = thr;
					G_Address[k].first = G_Address[k].last = 0;
				}
				for (INMOST_DATA_ENUM_TYPE i = cbeg; i < cend; ++i) //row-index
				{
					for (INMOST_DATA_ENUM_TYPE jt = Address[i].first; jt < Address[i].last; ++jt) //entry index
					{
						INMOST_DATA_ENUM_TYPE j = Entries[Address[i].thr][jt].first;
						if (tbeg <= j && j < tend)
							G_Address[j].last++;
					}
				}
				//Establish the addresses
				INMOST_DATA_ENUM_TYPE offset = 0;
				for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
				{
					G_Address[k].first += offset;
					G_Address[k].last += offset;
					offset += G_Address[k].Size();
				}
				//Allocate size for F storage
				G_Entries[thr].resize(G_Address[tend - 1].last);
				//Reset addresses to advance current fill position
				for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
					G_Address[k].last = G_Address[k].first;
				//Fill data for each thread
				for (INMOST_DATA_ENUM_TYPE i = cbeg; i < cend; ++i) //row-index
				{
					for (INMOST_DATA_ENUM_TYPE jt = Address[i].first; jt < Address[i].last; ++jt) //entry index
					{
						INMOST_DATA_ENUM_TYPE j = Entries[Address[i].thr][jt].first;
						if (tbeg <= j && j < tend)
							G_Entries[thr][G_Address[j].last++] = std::make_pair(i, jt);
							//Indices[j].push_back(std::make_pair(i, jt)); // Entries[j].first is column index, record row-index, entry index
					}
				}
			}
		}
	}

	void MLMTILUC_preconditioner::PrepareGraph(INMOST_DATA_ENUM_TYPE wbeg,
					  						   INMOST_DATA_ENUM_TYPE wend,
					  						   const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
					  						   const std::vector< std::vector<Sparse::Row::entry> >& Entries,
											   interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
											   std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& G_Entries)
					  						   //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G)
	{
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
		{
			//G[k].reserve(Address[k].Size());
			G_Address[k].thr = Thread();
			G_Address[k].first = (INMOST_DATA_ENUM_TYPE)G_Entries[G_Address[k].thr].size();
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) //if( !check_zero(Entries[it].second) )
				G_Entries[G_Address[k].thr].push_back(Entries[Address[k].thr][it].first);
				//G[k].push_back(Entries[Address[k].thr][it].first);
			G_Address[k].last = (INMOST_DATA_ENUM_TYPE)G_Entries[G_Address[k].thr].size();
		}
	}

	void MLMTILUC_preconditioner::PrepareGraphTranspose(INMOST_DATA_ENUM_TYPE wbeg,
							   							INMOST_DATA_ENUM_TYPE wend,
														const interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
														const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& G_Entries,
														interval<INMOST_DATA_ENUM_TYPE, Interval>& tG_Address,
														std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& tG_Entries)
							   							//const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
							   							//interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG)
	{
		//preallocate tG???
#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
		{
			int thr = Thread();
			INMOST_DATA_ENUM_TYPE cbeg = tG_Address.get_interval_beg();
			INMOST_DATA_ENUM_TYPE cend = tG_Address.get_interval_end();
			INMOST_DATA_ENUM_TYPE tseg = (cend - cbeg) / Threads();
			INMOST_DATA_ENUM_TYPE tbeg = cbeg + tseg * thr;
			INMOST_DATA_ENUM_TYPE tend = tbeg + tseg;
			if (thr == Threads() - 1) tend = cend;
			if (tend > tbeg)
			{
				for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
				{
					tG_Address[k].thr = thr;
					tG_Address[k].first = tG_Address[k].last = 0;
				}
				for (INMOST_DATA_ENUM_TYPE i = wbeg; i < wend; ++i) //row-index
				{
					for (INMOST_DATA_ENUM_TYPE jt = G_Address[i].first; jt < G_Address[i].last; ++jt) //entry index
					{
						INMOST_DATA_ENUM_TYPE j = G_Entries[G_Address[i].thr][jt];
						if (tbeg <= j && j < tend)
							tG_Address[j].last++;
					}
				}
				//Establish the addresses
				INMOST_DATA_ENUM_TYPE offset = 0;
				for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
				{
					tG_Address[k].first += offset;
					tG_Address[k].last += offset;
					offset += tG_Address[k].Size();
				}
				//Allocate size for F storage
				tG_Entries[thr].resize(tG_Address[tend - 1].last);
				//Reset addresses to advance current fill position
				for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
					tG_Address[k].last = tG_Address[k].first;
				//Fill data for each thread
				for (INMOST_DATA_ENUM_TYPE i = wbeg; i < wend; ++i) //row-index
				{
					for (INMOST_DATA_ENUM_TYPE jt = G_Address[i].first; jt < G_Address[i].last; ++jt) //entry index
					{
						INMOST_DATA_ENUM_TYPE j = G_Entries[G_Address[i].thr][jt];
						if (tbeg <= j && j < tend)
							tG_Entries[thr][tG_Address[j].last++] = i;
					}
				}
			}
		}
	}

	void MLMTILUC_preconditioner::PrepareGraphProduct(INMOST_DATA_ENUM_TYPE wbeg,
													  INMOST_DATA_ENUM_TYPE wend,
													  const interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
													  const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& G_Entries,
													  const interval<INMOST_DATA_ENUM_TYPE, Interval>& tG_Address,
													  const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& tG_Entries,
													  interval<INMOST_DATA_ENUM_TYPE, Interval>& pG_Address,
													  std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& pG_Entries)
							 						  //const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G,
							 						  //const interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG,
							 						  //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG)
	{
#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
		{
			interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > List(wbeg,wend,ENUMUNDEF);
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
			for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
			{
				INMOST_DATA_ENUM_TYPE Beg = EOL, nnz = 0;
				// go over connection of k-th row
				List[k] = Beg;
				Beg = k;
				//for(INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it)
				for (INMOST_DATA_ENUM_TYPE it = G_Address[k].first; it < G_Address[k].last; ++it)
				{
					//INMOST_DATA_ENUM_TYPE kt = G[k][it];
					INMOST_DATA_ENUM_TYPE kt = G_Entries[G_Address[k].thr][it];
					//for(INMOST_DATA_ENUM_TYPE jt = 0; jt < tG[kt].size(); ++jt)
					for (INMOST_DATA_ENUM_TYPE jt = tG_Address[kt].first; jt < tG_Address[kt].last; ++jt)
					{
						//INMOST_DATA_ENUM_TYPE lt = tG[kt][jt];
						INMOST_DATA_ENUM_TYPE lt = tG_Entries[tG_Address[kt].thr][jt];
						//insert to linked list
						if (List[lt] == ENUMUNDEF)
						{
							List[lt] = Beg;
							Beg = lt;
							nnz++;
						}
					}
				}
				// wadj_sep[k-wbeg] = nnz;
				int thr = Thread();
				pG_Address[k].thr = thr;
				pG_Address[k].first = (INMOST_DATA_ENUM_TYPE)pG_Entries[thr].size();
				//pG[k].reserve(nnz);
				INMOST_DATA_ENUM_TYPE it = Beg, jt;
				while( it != static_cast<INMOST_DATA_ENUM_TYPE>(k) )
				{
					//pG[k].push_back(it);
					pG_Entries[thr].push_back(it);
					jt = List[it]; //save next entry
					List[it] = ENUMUNDEF; //clear list
					it = jt; //restore next
				}
				List[k] = ENUMUNDEF;
				pG_Address[k].last = (INMOST_DATA_ENUM_TYPE)pG_Entries[thr].size();
			}
		}
	}

	int MLMTILUC_preconditioner::Thread()
	{
#if defined(USE_OMP)
		return omp_get_thread_num();
#else
		return 0;
#endif
	}

	int MLMTILUC_preconditioner::Threads()
	{
#if defined(USE_OMP)
		return omp_get_max_threads();
#else
		return 1;
#endif
	}

	void MLMTILUC_preconditioner::ColumnInterval(INMOST_DATA_ENUM_TYPE wbeg, 
						 						 INMOST_DATA_ENUM_TYPE wend,
						 						 const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 						 const std::vector< std::vector<Sparse::Row::entry> > & Entries,
						 						 INMOST_DATA_ENUM_TYPE & cbeg,
						 						 INMOST_DATA_ENUM_TYPE & cend)
	{
		cbeg = ENUMUNDEF, cend = 0;
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		{
			int Athr = Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = Address[k].first; it < Address[k].last; ++it) 
			{
				INMOST_DATA_ENUM_TYPE jt = Entries[Athr][it].first;
				if( jt > cend ) cend = jt;
				if( jt < cbeg ) cbeg = jt;
			}
		}
		cend++;
	}

	void MLMTILUC_preconditioner::FilterGraph(const Block & b,
											  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 		  					  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
											  const interval<INMOST_DATA_ENUM_TYPE, Interval>& in_G_Address,
											  const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > in_G_Entries,
											  interval<INMOST_DATA_ENUM_TYPE, Interval>& out_G_Address,
											  std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > out_G_Entries)
											  //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G_in,
											  //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G_out)
	{
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		out_G_Address.set_interval_beg(wbeg);
		out_G_Address.set_interval_end(wend);
		out_G_Entries.clear();
		out_G_Entries.resize(Threads());
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
		{
			int thr = Thread();
			out_G_Address[k].thr = thr;
			out_G_Address[k].first = (INMOST_DATA_ENUM_TYPE)out_G_Entries[thr].size();
			for (INMOST_DATA_ENUM_TYPE jt = in_G_Address[invP[k]].first; jt < in_G_Address[invP[k]].last; ++jt)
			{
				INMOST_DATA_ENUM_TYPE j = localQ[in_G_Entries[in_G_Address[invP[k]].thr][jt]];
				if (cbeg <= j && j < cend)
					out_G_Entries[thr].push_back(j);
			}
			out_G_Address[k].last = (INMOST_DATA_ENUM_TYPE)out_G_Entries[thr].size();
			/*
			G_out[k] = G_in[invP[k]];
			INMOST_DATA_ENUM_TYPE it = 0, jt = 0;
			while(it < G_out[k].size())
			{
				G_out[k][jt] = localQ[G_out[k][it++]]; //localQ is the new column position
				if( cbeg <= G_out[k][jt] && G_out[k][jt] < cend ) jt++; //column is inside block
			}
			G_out[k].resize(jt);
			*/
		}
	}

	void MLMTILUC_preconditioner::FilterGraphTranspose(const Block & b,
											  		  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		  					  		  const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invQ,
													  const interval<INMOST_DATA_ENUM_TYPE, Interval> & in_tG_Address,
													  const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > in_tG_Entries,
													  interval<INMOST_DATA_ENUM_TYPE, Interval> & out_tG_Address,
													  std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > out_tG_Entries)
											  		  //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG_in,
											  		  //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & tG_out)
	{
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		out_tG_Address.set_interval_beg(cbeg);
		out_tG_Address.set_interval_end(cend);
		out_tG_Entries.clear();
		out_tG_Entries.resize(Threads());
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(INMOST_DATA_INTEGER_TYPE k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
		{
			int thr = Thread();
			out_tG_Address[k].thr = thr;
			out_tG_Address[k].first = (INMOST_DATA_ENUM_TYPE)out_tG_Entries[thr].size();
			for (INMOST_DATA_ENUM_TYPE jt = in_tG_Address[invQ[k]].first; jt < in_tG_Address[invQ[k]].last; ++jt)
			{
				INMOST_DATA_ENUM_TYPE j = localP[in_tG_Entries[in_tG_Address[invQ[k]].thr][jt]];
				if (wbeg <= j && j < wend)
					out_tG_Entries[thr].push_back(j);
			}
			out_tG_Address[k].last = (INMOST_DATA_ENUM_TYPE)out_tG_Entries[thr].size();
			/*
			tG_out[k] = tG_in[invQ[k]];
			INMOST_DATA_ENUM_TYPE it = 0, jt = 0;
			while(it < tG_out[k].size())
			{
				tG_out[k][jt] = localP[tG_out[k][it++]]; //localP is the new row position
				if( wbeg <= tG_out[k][jt] && tG_out[k][jt] < wend ) jt++; //row is inside block
			}
			tG_out[k].resize(jt);
			*/
		}
	}

	void MLMTILUC_preconditioner::FilterGraphProduct(const Block & b,
													 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & invP,
				 		  							 const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
													 const interval<INMOST_DATA_ENUM_TYPE, Interval>& in_pG_Address,
													 const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > in_pG_Entries,
													 interval<INMOST_DATA_ENUM_TYPE, Interval>& out_pG_Address,
													 std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > out_pG_Entries)
													 //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG_in,
													 //interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & pG_out)
	{
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end;
		out_pG_Address.set_interval_beg(wbeg);
		out_pG_Address.set_interval_end(wend);
		out_pG_Entries.clear();
		out_pG_Entries.resize(Threads());
#if defined(USE_OMP)
#pragma omp parallel for
#endif	
		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
		{
			int thr = Thread();
			out_pG_Address[k].thr = thr;
			out_pG_Address[k].first = (INMOST_DATA_ENUM_TYPE)out_pG_Entries[thr].size();
			for (INMOST_DATA_ENUM_TYPE jt = in_pG_Address[invP[k]].first; jt < in_pG_Address[invP[k]].last; ++jt)
			{
				INMOST_DATA_ENUM_TYPE j = localP[in_pG_Entries[in_pG_Address[invP[k]].thr][jt]];
				if (wbeg <= j && j < wend)
					out_pG_Entries[thr].push_back(j);
			}
			out_pG_Address[k].last = (INMOST_DATA_ENUM_TYPE)out_pG_Entries[thr].size();
			/*
			pG_out[k] = pG_in[invP[k]];
			INMOST_DATA_ENUM_TYPE it = 0, jt = 0;
			while(it < pG_out[k].size())
			{
				pG_out[k][jt] = localP[pG_out[k][it++]]; //localP is the new row position
				if( wbeg <= pG_out[k][jt] && pG_out[k][jt] < wend ) jt++; //row is inside block
			}
			pG_out[k].resize(jt);
			*/
		}
	}

	void MLMTILUC_preconditioner::DumpGraph(std::string name, 
											const interval<INMOST_DATA_ENUM_TYPE, Interval>& G_Address,
											const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > G_Entries)
											//interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > & G)
	{
		INMOST_DATA_ENUM_TYPE wbeg, wend, cbeg, cend, nnz = 0, side;
		wbeg = G_Address.get_interval_beg();
		wend = G_Address.get_interval_end();
		cbeg = ENUMUNDEF;
		cend = 0;
		std::cout << __FUNCTION__ << " " << name << std::endl;
		std::cout << "wbeg " << wbeg << " wend " << wend << std::endl;
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			for(INMOST_DATA_ENUM_TYPE it = G_Address[k].first; it < G_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE jt = G_Entries[G_Address[k].thr][it];
				if( cbeg > jt ) cbeg = jt;
				if( cend < jt ) cend = jt;
				nnz++;
			}
		}
		std::cout << "cbeg " << cbeg << " cend " << cend << std::endl;
		side = std::max(wend-wbeg,cend-cbeg);
		std::ofstream file(name.c_str());
		file << "%%MatrixMarket matrix coordinate real general\n";
		file << "%graph only, values are units\n";
		file << "%writing as square graph, some rows/cols may be missing\n";
		file << side << " " << side << " " << nnz << std::endl;
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			for(INMOST_DATA_ENUM_TYPE it = G_Address[k].first; it < G_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE jt = G_Entries[G_Address[k].thr][it];
				file << k-wbeg+1 << " " << jt-cbeg+1 << " 1\n";
			}
		}
		file.close();
	}

	void MLMTILUC_preconditioner::NestedDissection(INMOST_DATA_ENUM_TYPE wbeg,
												   INMOST_DATA_ENUM_TYPE wend,
												   const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 						   const std::vector< std::vector<Sparse::Row::entry> > & Entries,
						 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
												   std::vector<Block> & blocks,
												   INMOST_DATA_ENUM_TYPE max_size)
	{
		const int kway_parts = 2;
		double /*timer = Timer(),*/ total_time = Timer();
		INMOST_DATA_ENUM_TYPE cbeg, cend, sep, blks;
		ColumnInterval(wbeg,wend,Address,Entries,cbeg,cend);

		interval< INMOST_DATA_ENUM_TYPE, Interval > G_Address(wbeg,wend), tG_Address(cbeg,cend), pG_Address(wbeg,wend);
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > G_Entries(Threads()), tG_Entries(Threads()), pG_Entries(Threads());
		interval< INMOST_DATA_ENUM_TYPE, Interval > G2_Address, tG2_Address, pG2_Address;
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > G2_Entries(Threads()), tG2_Entries(Threads()), pG2_Entries(Threads());
		interval< INMOST_DATA_ENUM_TYPE, Interval > * rG_Address = &G_Address, * rtG_Address = &tG_Address, * rpG_Address = &pG_Address;
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >* rG_Entries = &G_Entries, * rtG_Entries = &tG_Entries, * rpG_Entries = &pG_Entries;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>  P(wbeg,wend), Q(cbeg,cend), invP(wbeg,wend), invQ(cbeg,cend);
		std::vector<Block> blocks_new;

		// std::cout << __FUNCTION__ << " row " << wbeg << ":" << wend << " col " << cbeg << ":" << cend << std::endl;
		
		//~ timer = Timer();

		PrepareGraph(wbeg,wend,Address,Entries,G_Address,G_Entries);
		PrepareGraphTranspose(wbeg,wend,G_Address,G_Entries,tG_Address,tG_Entries);

		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) localP[k] = invP[k] = k;
		for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) localQ[k] = invQ[k] = k;
		
		// std::cout << "prepare tG time " << Timer() - timer << std::endl;
		// if( wgt_sep )
		{
			//~ timer = Timer();
			PrepareGraphProduct(wbeg,wend,G_Address,G_Entries,tG_Address,tG_Entries,pG_Address,pG_Entries);
			// std::cout << "prepare pG time " << Timer() - timer << std::endl;
		}

		// DumpGraph("G.mtx",G);
		// DumpGraph("tG.mtx",tG);
		// DumpGraph("pG.mtx",pG);
		blocks.push_back(Block(wbeg,wend,cbeg,cend));

		//find largest non-separator block
		INMOST_DATA_ENUM_TYPE cur = 0;//, cnt = 0;
		do
		{
			// std::cout << "separate block " << cur << " rows " << blocks[cur].RowSize() << " cols " << blocks[cur].ColSize() << std::endl;
			GreedyDissection(blocks[cur],*rG_Address,*rG_Entries,*rtG_Address,*rtG_Entries,*rpG_Address,*rpG_Entries,P,Q,blocks_new,kway_parts);//,1,0,0);
			// compute reordering in global P,Q, we need it to compute reordering in vector during solve phase
			for (INMOST_DATA_ENUM_TYPE k = blocks[cur].row_start; k < blocks[cur].row_end; ++k) localP[invP[k]] = P[k];
			for (INMOST_DATA_ENUM_TYPE k = blocks[cur].col_start; k < blocks[cur].col_end; ++k) localQ[invQ[k]] = Q[k];
			//replace blocks
			blocks.insert(blocks.erase(blocks.begin()+cur),blocks_new.begin(),blocks_new.end());
			blocks_new.clear();
			
			for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k) if( !blocks[k].separator )
			{
				if( cur == ENUMUNDEF || blocks[cur].RowSize() < blocks[k].RowSize() )
					cur = k;
			}
			// std::cout << "next selected block " << cur << " rows " << blocks[cur].RowSize() << " cols " << blocks[cur].ColSize() << std::endl;
			// break;
			//prepare reduced graphs
			//TODO: do not use entire initial graph!
			//inverse permutation

			for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) invP[localP[k]] = k;
			for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) invQ[localQ[k]] = k;

			// if( transpose )
			// {
			// 	FilterGraph(blocks[cur],invP,localQ,tG,G2);
			// 	FilterGraphTranspose(blocks[cur],localP,invQ,G,tG2);
			// 	FilterGraphProduct(blocks[cur],invP,localP,pG,pG2);
			// }
			// else
			{
				FilterGraph(blocks[cur],invP,localQ,G_Address,G_Entries,G2_Address,G2_Entries);
				FilterGraphTranspose(blocks[cur],localP,invQ,tG_Address,tG_Entries,tG2_Address,tG2_Entries);
				FilterGraphProduct(blocks[cur],invP,localP,pG_Address,pG_Entries,pG2_Address,pG2_Entries);
			}


			rG_Address = &G2_Address;
			rG_Entries = &G2_Entries;
			rtG_Address = &tG2_Address;
			rtG_Entries = &tG2_Entries;
			rpG_Address = &pG2_Address;
			rpG_Entries = &pG2_Entries;
			// DumpGraph("rG.mtx",G2);
			// DumpGraph("rtG.mtx",tG2);
			// DumpGraph("rpG.mtx",pG2);
			//break;
			// scanf("%*c");
			// cnt++;
		}
		while( blocks[cur].RowSize() > max_size );//&& cnt < 2);

		blks = sep = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
		{
			if( blocks[k].separator ) sep += blocks[k].RowSize();
			else blks++;
			//std::cout << (blocks[k].separator?"separator":"block") << "[" << k << "] rows " << blocks[k].row_start << ":" << blocks[k].row_end << "(" << blocks[k].RowSize() << ") cols " << blocks[k].col_start << ":" << blocks[k].col_end << "(" << blocks[k].ColSize() << ")" << std::endl;
		}
		//std::cout << "total separator " << sep << "/" << wend-wbeg << " blocks " << blks << std::endl;


		//std::cout << __FUNCTION__ << " time " << Timer() - total_time << std::endl;

	}

	void MLMTILUC_preconditioner::KwayDissection(INMOST_DATA_ENUM_TYPE wbeg,
												 INMOST_DATA_ENUM_TYPE wend,
												 const interval<INMOST_DATA_ENUM_TYPE, Interval> & Address,
						 						 const std::vector< std::vector<Sparse::Row::entry> > & Entries,
						 						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		 						 interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
												 std::vector<Block> & blocks, int parts)
	{
		double timer = Timer(), total_time = Timer();
		INMOST_DATA_ENUM_TYPE cbeg, cend, sep, blks;
		ColumnInterval(wbeg,wend,Address,Entries,cbeg,cend);


		interval< INMOST_DATA_ENUM_TYPE, Interval > G_Address(wbeg, wend), tG_Address(cbeg, cend), pG_Address(wbeg, wend);
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > G_Entries(Threads()), tG_Entries(Threads()), pG_Entries(Threads());

		
		// std::cout << __FUNCTION__ << " row " << wbeg << ":" << wend << " col " << cbeg << ":" << cend << std::endl;
		
		//timer = Timer();
		PrepareGraph(wbeg,wend,Address,Entries,G_Address,G_Entries);
		//std::cout << "prepare G time " << Timer() - timer << std::endl;
		//timer = Timer();
		PrepareGraphTranspose(wbeg,wend,G_Address,G_Entries,tG_Address,tG_Entries);
		//std::cout << "prepare tG time " << Timer() - timer << std::endl;
		//timer = Timer();
		PrepareGraphProduct(wbeg,wend,G_Address,G_Entries,tG_Address,tG_Entries,pG_Address,pG_Entries);
		//std::cout << "prepare pG time " << Timer() - timer << std::endl;
		//timer = Timer();
		GreedyDissection(Block(wbeg,wend,cbeg,cend),G_Address,G_Entries,tG_Address,tG_Entries,pG_Address,pG_Entries,localP,localQ,blocks,parts);
		//std::cout << "greedy dissection " << Timer() - timer << std::endl;
		blks = sep = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
		{
			if( blocks[k].separator ) sep += blocks[k].RowSize();
			else blks++;
			//std::cout << (blocks[k].separator?"separator":"block") << "[" << k << "] rows " << blocks[k].row_start << ":" << blocks[k].row_end << "(" << blocks[k].RowSize() << ") cols " << blocks[k].col_start << ":" << blocks[k].col_end << "(" << blocks[k].ColSize() << ")" << std::endl;
			
		}
		//std::cout << "total separator " << sep << " blocks " << blks << std::endl;


		//std::cout << __FUNCTION__ << " time " << Timer() - total_time << std::endl;
	}
	void MLMTILUC_preconditioner::KwaySymmetricDissection(INMOST_DATA_ENUM_TYPE wbeg,
		INMOST_DATA_ENUM_TYPE wend,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& Address,
		const std::vector< std::vector<Sparse::Row::entry> >& Entries,
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localP,
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
		std::vector<Block>& blocks, int parts)
	{
		double timer = Timer(), total_time = Timer();
		INMOST_DATA_ENUM_TYPE cbeg, cend, sep, blks;
		ColumnInterval(wbeg, wend, Address, Entries, cbeg, cend);

		interval< INMOST_DATA_ENUM_TYPE, Interval > G_Address(wbeg, wend), tG_Address(cbeg, cend), pG_Address(wbeg, wend);
		std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > G_Entries(Threads()), tG_Entries(Threads()), pG_Entries(Threads());
		//interval< INMOST_DATA_ENUM_TYPE, std::vector<INMOST_DATA_ENUM_TYPE> > G(wbeg, wend), tG(cbeg, cend), pG(wbeg, wend);

		// std::cout << __FUNCTION__ << " row " << wbeg << ":" << wend << " col " << cbeg << ":" << cend << std::endl;

		//timer = Timer();
		PrepareGraph(wbeg, wend, Address, Entries, G_Address, G_Entries);
		//std::cout << "prepare G time " << Timer() - timer << std::endl;
		//timer = Timer();
		PrepareGraphTranspose(wbeg, wend, G_Address, G_Entries, tG_Address, tG_Entries);
		//std::cout << "prepare tG time " << Timer() - timer << std::endl;
		//timer = Timer();
		PrepareGraphProduct(wbeg, wend, G_Address, G_Entries, tG_Address, tG_Entries, pG_Address, pG_Entries);
		//std::cout << "prepare pG time " << Timer() - timer << std::endl;
		//timer = Timer();
		GreedyDissection(Block(wbeg, wend, cbeg, cend), G_Address, G_Entries, tG_Address, tG_Entries, pG_Address, pG_Entries, localP, localQ, blocks, parts);
		//std::cout << "greedy dissection " << Timer() - timer << std::endl;

		blks = sep = 0;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
		{
			if (blocks[k].separator) sep += blocks[k].RowSize();
			else blks++;
			//std::cout << (blocks[k].separator?"separator":"block") << "[" << k << "] rows " << blocks[k].row_start << ":" << blocks[k].row_end << "(" << blocks[k].RowSize() << ") cols " << blocks[k].col_start << ":" << blocks[k].col_end << "(" << blocks[k].ColSize() << ")" << std::endl;

		}

		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
			localQ[k] = localP[k];
		for (INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k) if( !blocks[k].separator )
		{
			blocks[k].col_start = blocks[k].row_start;
			blocks[k].col_end = blocks[k].row_end;
		}
		//std::cout << "total separator " << sep << " blocks " << blks << std::endl;


		//std::cout << __FUNCTION__ << " time " << Timer() - total_time << std::endl;
	}

	void MLMTILUC_preconditioner::GreedyDissection(const Block & b,
												   const interval< INMOST_DATA_ENUM_TYPE, Interval > & G_Address,
												   const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > & G_Entries,
						  						   const interval< INMOST_DATA_ENUM_TYPE, Interval > & tG_Address,
												   const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& tG_Entries,
						  						   const interval< INMOST_DATA_ENUM_TYPE, Interval > & pG_Address,
												   const std::vector< std::vector<INMOST_DATA_ENUM_TYPE> >& pG_Entries,
						 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
				 		 						   interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
												   std::vector<Block> & blocks, int kway_parts)
	{
		// double timer = Timer(), total_time = Timer();
		
		// const bool kway = false;
		bool kway = (kway_parts > 1);
		// const int kway_parts = 4;
		const int upd_sep = 1, upd_blk = 1;
		const int wgt_sep = 1, wgt_blk = 1;

		
		// std::cout << __FUNCTION__ <<  " wgt sep " << wgt_sep << " blk " << wgt_blk << " kway " << kway << std::endl;
		
		INMOST_DATA_ENUM_TYPE wbeg = b.row_start, wend = b.row_end, cbeg = b.col_start, cend = b.col_end;
		INMOST_DATA_ENUM_TYPE index_col = cbeg, index_row = wbeg, cur = ENUMUNDEF;//, nxt;
		INMOST_DATA_ENUM_TYPE Beg, Beg_upd, sep_size = 0, row_size = 0, col_size = 0, tot_sep_size = 0;//, nnz;
		INMOST_DATA_ENUM_TYPE index_row_beg = index_row, index_col_beg = index_col;
		INMOST_DATA_REAL_TYPE min_ratio = 1.0e+20, ratio, block_mean, balance_ratio, optsep_ratio;
		INMOST_DATA_ENUM_TYPE min_index = ENUMUNDEF;
		int had_col, had_sep, had_wgt_blk, had_wgt_sep;

		std::vector<int> wadj_blk(wend-wbeg,0), wadj_sep(wend-wbeg,0);
		BinaryHeapCustom<int, std::less<int> > q0(wend-wbeg);
		interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > List(wbeg,wend,ENUMUNDEF);
		interval< INMOST_DATA_ENUM_TYPE, std::pair<int,INMOST_DATA_ENUM_TYPE> > List_upd(wbeg,wend,std::make_pair(0,ENUMUNDEF));
		std::vector<INMOST_DATA_ENUM_TYPE> visits;
		std::vector< std::pair<INMOST_DATA_ENUM_TYPE, bool> > upd;

		if( !kway ) visits.reserve(wend-wbeg);

		// std::cout << "alloc time " << Timer() - timer << std::endl;
		// timer = Timer();

		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		{
			localP[k] = ENUMUNDEF;
			if( wgt_blk ) wadj_blk[k-wbeg] = (int)G_Address[k].Size(); //row connection wgt
			if( wgt_sep ) wadj_sep[k-wbeg] = (int)pG_Address[k].Size();
		}

		for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) localQ[k] = ENUMUNDEF;


		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
			q0.PushHeap(k-wbeg,wgt_sep*wadj_sep[k-wbeg]+wgt_blk*wadj_blk[k-wbeg]); 

		// std::cout << "prepare heap time " << Timer() - timer << std::endl;
		// timer = Timer();
		// exit(-1);
		
		
		// cur = 0;
		// for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) 
		// 	if( wadj_sep[k-wbeg] > wadj_sep[cur-wbeg] ) cur = k;
		// srand(time(NULL));
		// cur = rand()%(wend-wbeg)+wbeg;
		// std::cout << "start at " << cur << std::endl;
		// goto start;
		// std::ofstream rec("rec.txt");
		Beg = EOL;
		Beg_upd = EOL;
		{
			while( !q0.Empty() ) 
			{
				//cur has minimal separator growth size
				cur = q0.PopHeap() + wbeg; //row index
			// start:
				if( localP[cur] != ENUMUNDEF ) continue; //skip separator of other blocks
				had_col = index_col;
				had_sep = sep_size;
				had_wgt_blk = wadj_blk[cur-wbeg]; 
				had_wgt_sep = wadj_sep[cur-wbeg];
				// rec << "(" << cur << "," << wadj_blk[cur-wbeg] << "," << wadj_sep[cur-wbeg] << ")";
				// rec << " row " << index_row << " col " << index_col << " sep " << sep_size;
				localP[cur] = index_row++;
				row_size++;
				if( List[cur] != ENUMUNDEF ) //uncount me from separator
					sep_size--;
				else //include me to linked list to avoid count as separator
				{
					List[cur] = Beg;
					Beg = cur;
					//extract me from rows that depend on me
					if( wgt_sep ) upd.push_back(std::make_pair(cur,false));
				}
				// go over connection of enumerated row
				int col_skip = 0;
				int sep_skip = 0;
				for(INMOST_DATA_ENUM_TYPE it = G_Address[cur].first; it < G_Address[cur].last; ++it)
				{
					INMOST_DATA_ENUM_TYPE kt = G_Entries[G_Address[cur].thr][it]; //column index
					if( localQ[kt] == ENUMUNDEF ) 
					{
						localQ[kt] = index_col++;
						col_size++;
						for(INMOST_DATA_ENUM_TYPE jt = tG_Address[kt].first; jt < tG_Address[kt].last; ++jt)
						{
							INMOST_DATA_ENUM_TYPE lt = tG_Entries[tG_Address[kt].thr][jt]; //row index
							wadj_blk[lt-wbeg]-=upd_blk;
							// if( wgt_blk && q0.Contains(lt-wbeg) )
							// 	q0.ChangeKey(lt-wbeg,wgt_sep*wadj_sep[lt-wbeg]+wgt_blk*wadj_blk[lt-wbeg]);
							if( q0.Contains(lt-wbeg) )
							{
								List_upd[lt].first += wgt_blk*upd_blk;
								if (List_upd[lt].second == ENUMUNDEF)
								{
									List_upd[lt].second = Beg_upd;
									Beg_upd = lt;
								}
							}
							// add dependencies to linked-list of separator
							if( List[lt] == ENUMUNDEF ) //new row into separator
							{
								List[lt] = Beg;
								Beg = lt;
								sep_size++;
								//extract me from rows that depend on me
								if( wgt_sep ) upd.push_back(std::make_pair(lt,true));
							}
							else sep_skip++;
						}
					}
					else col_skip++;
				}
				// collect updates
				for(INMOST_DATA_ENUM_TYPE it = 0; it < upd.size(); ++it)
				{
					INMOST_DATA_ENUM_TYPE lt = upd[it].first; //row index
					for(INMOST_DATA_ENUM_TYPE ut = pG_Address[lt].first; ut < pG_Address[lt].last; ++ut)
					{
						INMOST_DATA_ENUM_TYPE qt = pG_Entries[pG_Address[lt].thr][ut]; //row index
						if( localP[qt] == ENUMUNDEF )
						{
							wadj_sep[qt-wbeg]-=upd_sep;
							if( q0.Contains(qt-wbeg) )
							{
								List_upd[qt].first += wgt_sep*upd_sep;
								if (List_upd[qt].second == ENUMUNDEF)
								{
									List_upd[qt].second = Beg_upd;
									Beg_upd = qt;
								}
							}
						}
					}
					if( upd[it].second )
					{
						wadj_sep[lt-wbeg]-=upd_sep; //taking me will reduce separator size
						if( q0.Contains(lt-wbeg) )
						{
							List_upd[lt].first += wgt_sep*upd_sep;
							if (List_upd[lt].second == ENUMUNDEF)
							{
								List_upd[lt].second = Beg_upd;
								Beg_upd = lt;
							}
						}
					}
				}
				// update weights
				{
					INMOST_DATA_ENUM_TYPE it = Beg_upd,jt;
					while(it != EOL)
					{
						q0.ChangeKey(it-wbeg, q0.GetKey(it-wbeg)-List_upd[it].first);
						jt = List_upd[it].second;
						List_upd[it].first = 0;
						List_upd[it].second = ENUMUNDEF;
						it = jt;
					}
					Beg_upd = EOL;
				}
				upd.clear();
				if( wgt_blk && upd_blk == 1 && had_wgt_blk != index_col - had_col )
					std::cout << "Wrong prediction on column growth: " << had_wgt_blk << " real growth "  << index_col - had_col  << std::endl;

				if( wgt_sep && upd_sep == 1 && had_wgt_sep != ((int)sep_size - had_sep) )
					std::cout << "Wrong prediction on separtor growth: " << had_wgt_sep << " real growth "  << ((int)sep_size - had_sep)  << std::endl;

				//disconnect each row that connects to me to know block growth size
				if( !kway )
				{
					visits.push_back(cur);
					// block_mean = sqrt(((INMOST_DATA_REAL_TYPE)row_size)*((INMOST_DATA_REAL_TYPE)wend-wbeg-row_size-sep_size));
					block_mean = ((INMOST_DATA_REAL_TYPE)std::min(row_size,wend-wbeg-row_size-sep_size));
					// balance_ratio = pow(0.5*((INMOST_DATA_REAL_TYPE)wend-wbeg-sep_size) / block_mean,5);
					balance_ratio = 1;
					optsep_ratio = ((INMOST_DATA_REAL_TYPE)sep_size) / block_mean;
					ratio = optsep_ratio * balance_ratio;
					if( ratio < min_ratio )
					{
						min_ratio = ratio;
						min_index = index_row;
					}
				}
				else
				{
					bool take_next = false;
					// if( !q0.Empty() )
					// {
					// 	nxt = q0.PeekHeap()+wbeg;
					// 	take_next = (wgt_sep && upd_sep == 1 && nxt != ENUMUNDEF && wadj_sep[nxt-wbeg] <= 0);
					// }
					if( (row_size > (wend-wbeg-tot_sep_size-sep_size)/kway_parts) && !take_next )
					// if( col_size > (cend-cbeg)/kway_parts && !take_next )
					// if( row_size > 5*sep_size && !take_next )
					{
						
						INMOST_DATA_ENUM_TYPE it = Beg;//, jt;
						while( it != EOL )
						{
							if( localP[it] == ENUMUNDEF ) 
							{
								localP[it] = ENUMUNDEF-1; //enumerate separator
								tot_sep_size++;
							}
							it = List[it];
						}
						// std::cout << "make new block! total separator " << tot_sep_size << " gained " << sep_size << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
						blocks.push_back(Block(index_row_beg,index_row,index_col_beg,index_col));
						index_row_beg = index_row;
						index_col_beg = index_col;
						col_size = row_size = sep_size = 0;
					}
				}
			}
		}
		// rec.close();
		if( !kway )
		{
			// std::cout << "simulation time " << Timer() - timer << std::endl;
			// timer = Timer();

			//clear linked-list
			{
				//interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > & List = Lists[0];
				INMOST_DATA_ENUM_TYPE it = Beg, jt;
				while( it != EOL )
				{
					jt = List[it]; //save next entry
					List[it] = ENUMUNDEF; //clear list
					it = jt; //restore next
				}
				Beg = EOL;
			}
		
			for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) localP[k] = ENUMUNDEF;
			for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) localQ[k] = ENUMUNDEF;

			index_row_beg = index_row = wbeg;
			index_col_beg = index_col = cbeg;
			//fill first block
			for(INMOST_DATA_ENUM_TYPE k = 0; k < min_index-wbeg; ++k)
			{
				//interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > & List = Lists[0];
				cur = visits[k];
				localP[cur] = index_row++;
				for(INMOST_DATA_ENUM_TYPE it = G_Address[cur].first; it < G_Address[cur].last; ++it)
				{
					INMOST_DATA_ENUM_TYPE kt = G_Entries[G_Address[cur].thr][it];
					if( localQ[kt] == ENUMUNDEF ) 
						localQ[kt] = index_col++;
					for(INMOST_DATA_ENUM_TYPE jt = tG_Address[kt].first; jt < tG_Address[kt].last; ++jt)
					{
						INMOST_DATA_ENUM_TYPE lt = tG_Entries[tG_Address[kt].thr][jt];
						if (List[lt] == ENUMUNDEF) //new row into separator
						{
							List[lt] = Beg;
							Beg = lt;
						}
					}
				}
			}
			//mark separator
			tot_sep_size = 0;
			{
				//interval< INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE > & List = Lists[0];
				INMOST_DATA_ENUM_TYPE it = Beg;//, jt;
				while( it != EOL )
				{
					if( localP[it] == ENUMUNDEF ) 
					{
						localP[it] = ENUMUNDEF-1; //enumerate separator
						tot_sep_size++;
					}
					it = List[it];
				}
			}
			if( index_row_beg != index_row )
			{
				// std::cout << "make new block! total separator " << tot_sep_size  << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
				blocks.push_back(Block(index_row_beg,index_row,index_col_beg,index_col));
				index_row_beg = index_row;
				index_col_beg = index_col;
			}
			//fill second block preserving order and skipping separator
			for(INMOST_DATA_ENUM_TYPE k = min_index-wbeg; k < wend-wbeg; ++k)
			{
				cur = visits[k];
				if( localP[cur] != ENUMUNDEF ) continue;
				localP[cur] = index_row++;
				for(INMOST_DATA_ENUM_TYPE it = G_Address[cur].first; it < G_Address[cur].last; ++it)
				{
					INMOST_DATA_ENUM_TYPE kt = G_Entries[G_Address[cur].thr][it];
					if( localQ[kt] == ENUMUNDEF ) 
						localQ[kt] = index_col++;
				}
			}
		}
		//add last block
		if( index_row_beg != index_row )
		{
			// std::cout << "make new block! total separator " << tot_sep_size << " rows " << index_row-index_row_beg << " cols " << index_col-index_col_beg << std::endl;
			blocks.push_back(Block(index_row_beg,index_row,index_col_beg,index_col));
			index_row_beg = index_row;
			index_col_beg = index_col;
		}

		// std::cout << "blocks: " << blocks.size() << " total separator " << tot_sep_size << "/" << wend-wbeg << std::endl;

		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) if( localP[k] == ENUMUNDEF-1 ) localP[k] = index_row++;
		for(INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k) if( localQ[k] == ENUMUNDEF   ) localQ[k] = index_col++;

		if( index_row_beg != index_row )
		{
			// std::cout << "make separator block! total separator " << tot_sep_size << " rows " << index_row-index_row_beg << " extra cols " << index_col-index_col_beg << std::endl;
			blocks.push_back(Block(index_row_beg,index_row,cbeg,index_col,true));
		}

		//check blocks
#if 0
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			INMOST_DATA_ENUM_TYPE r = localP[k]; //new k-th row number
			INMOST_DATA_ENUM_TYPE b = ENUMUNDEF; //block number
			for(INMOST_DATA_ENUM_TYPE it = 0; it < blocks.size(); ++it)
				if( blocks[it].row_start <= r && r < blocks[it].row_end )
				{
					b = it;
					break;
				}
			if( b == ENUMUNDEF )
			{
				std::cout << "cannot find block number for row " << k << " new position " << localP[k] << std::endl;
				continue;
			}
			for(INMOST_DATA_ENUM_TYPE it = 0; it < G[k].size(); ++it)
			{
				INMOST_DATA_ENUM_TYPE c = localQ[G[k][it]]; //new column position
				if( !(blocks[b].col_start <= c && c < blocks[b].col_end) )
				{
					std::cout << "position (" << k << "," << G[k][it] << ")";
					std::cout << " new (" << r << "," << c << ") outside block ";
					std::cout << "(" << blocks[b].row_start << ":" << blocks[b].row_end << ",";
					std::cout << blocks[b].col_start << ":" << blocks[b].col_end << ")" << std::endl;
				}
			}
		}
#endif
		// std::cout << "construction time " << Timer() - timer << std::endl;
		// std::cout << __FUNCTION__ << " time " << Timer() - total_time << std::endl;

		// exit(-1);
		
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
		condestL = condestU = 1;
		verbosity = 0;
		Alink = NULL;
		init = false;
		sciters = 8;
		eps = 1e-54;
		estimator = 1;
		tau = 1.0e-3;
		iluc2_tau = tau*tau;
		pivot_cond = PIVOT_COND_DEFAULT;
		pivot_diag = PIVOT_DIAG_DEFAULT;
	}
	bool MLMTILUC_preconditioner::isInitialized() { return init; }
	bool MLMTILUC_preconditioner::isFinalized() { return !init; }
	
	
	void MLMTILUC_preconditioner::MaximalTransversal(	INMOST_DATA_ENUM_TYPE rbeg, //matrix row begin
														INMOST_DATA_ENUM_TYPE rend, //matrix row end
														INMOST_DATA_ENUM_TYPE cbeg, //matrix row begin
														INMOST_DATA_ENUM_TYPE cend, //matrix row end
														const interval<INMOST_DATA_ENUM_TYPE, Interval> & A_Address,
														const std::vector< std::vector<Sparse::Row::entry> > & A_Entries,
														interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
														interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & DL,
														interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & DR)
	{
		//Sparse::Vector & U = DL;
		//Sparse::Vector & V = DR;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Cmax(cbeg,cend,0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & Perm = localQ;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> IPerm(rbeg,rend,ENUMUNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ColumnList(cbeg,cend,ENUMUNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Parent(rbeg,rend,ENUMUNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> AugmentPosition(cbeg,cend,ENUMUNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> ColumnPosition(rbeg,rend,ENUMUNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(cbeg, cend, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(rbeg, rend, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Dist(cbeg, cend + 1, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		std::vector<INMOST_DATA_REAL_TYPE> C_Entries;// (A_Entries.size());
		std::vector<INMOST_DATA_ENUM_TYPE> C_Row;
		BinaryHeap Heap(&Dist[cbeg],cend-cbeg);


		//std::fill(U.begin() + wbeg - mobeg, U.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		//std::fill(V.begin() + wbeg - mobeg, V.begin() + wend - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		//std::fill(Cmax.begin() + wbeg - mobeg, Cmax.begin() + wend - mobeg, 0.0);
		//std::fill(Dist.begin() + wbeg - mobeg, Dist.begin() + wend + 1 - mobeg, std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());

		for (INMOST_DATA_INTEGER_TYPE k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
			Perm[k] = ENUMUNDEF;

		if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		
		// Initial LOG transformation to dual problem and initial extreme match
		C_Row.push_back(0);
		for(INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
		{
			int Athr = A_Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE i = A_Entries[Athr][it].first;
				INMOST_DATA_REAL_TYPE u = fabs(A_Entries[Athr][it].second);
				C_Entries.push_back(u);
				Cmax[i] = std::max(Cmax[i],u);
			}
			assert(C_Row.back() != C_Entries.size()); //no empty row
			C_Row.push_back(static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries.size()));
		}

		for(INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
		{
			int Athr = A_Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE i = A_Entries[Athr][it].first;
				INMOST_DATA_ENUM_TYPE Ci = C_Row[k - rbeg] + it - A_Address[k].first;
				if( check_zero(Cmax[i]) || check_zero(C_Entries[Ci]) )
					C_Entries[Ci] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
				else
				{
					INMOST_DATA_REAL_TYPE logCmax = log(Cmax[i]);
					INMOST_DATA_REAL_TYPE logCval = log(C_Entries[Ci]);
					INMOST_DATA_REAL_TYPE logCdif = logCmax - logCval;
					C_Entries[Ci] = logCdif;
					U[i] = std::min(U[i], C_Entries[Ci]);
				}
			}
		}
		for(INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
		{
			int Athr = A_Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE Ci = C_Row[k - rbeg] + it - A_Address[k].first;
				INMOST_DATA_REAL_TYPE u = C_Entries[Ci] - U[A_Entries[Athr][it].first];
				V[k] = std::min(V[k], u);
			}
		}
		/// Update cost and match
		for(INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
		{
			int Athr = A_Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE Ci = C_Row[k - rbeg] + it - A_Address[k].first;
				INMOST_DATA_REAL_TYPE u = fabs(C_Entries[Ci] - V[k] - U[A_Entries[Athr][it].first]);
				if( check_zero(u) && Perm[A_Entries[Athr][it].first] == ENUMUNDEF && IPerm[k] == ENUMUNDEF )
				{
					 Perm[A_Entries[Athr][it].first] = k;
					 IPerm[k] = A_Entries[Athr][it].first;
					 ColumnPosition[k] = Ci;
				}
			}
		}
		/// 1-step augmentation
		for(INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
		{
			if( IPerm[k] == ENUMUNDEF ) //unmatched row
			{
				int Athrk = A_Address[k].thr;
				for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last && IPerm[k] == ENUMUNDEF; ++it)
				{
					INMOST_DATA_ENUM_TYPE Ci = C_Row[k - rbeg] + it - A_Address[k].first;
					INMOST_DATA_REAL_TYPE u = fabs(C_Entries[Ci] - V[k] - U[A_Entries[Athrk][it].first]);
					if( check_zero(u) )
					{
						INMOST_DATA_ENUM_TYPE Li = Perm[A_Entries[Athrk][it].first];
						assert(Li != ENUMUNDEF);
						// Search other row in C for 0
						int AthrLi = A_Address[Li].thr;
						for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit < A_Address[Li].last; ++Lit)
						{
							INMOST_DATA_ENUM_TYPE CLi = C_Row[Li - rbeg] + Lit - A_Address[Li].first;
							u = fabs(C_Entries[CLi]- V[Li] - U[A_Entries[AthrLi][Lit].first]);
							if( check_zero(u) && Perm[A_Entries[AthrLi][Lit].first] == ENUMUNDEF )
							{
								Perm[A_Entries[Athrk][it].first] = k;
								IPerm[k] = A_Entries[Athrk][it].first;
								ColumnPosition[k] = Ci;
								Perm[A_Entries[AthrLi][Lit].first] = Li;
								IPerm[Li] = A_Entries[AthrLi][Lit].first;
								ColumnPosition[Li] = CLi;
								break;
							}
						}
					}
				}
			}
		}
		/// Weighted bipartite matching
		for(INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
		{
			if( IPerm[k] != ENUMUNDEF )
				continue;
			INMOST_DATA_ENUM_TYPE Li = k;
			INMOST_DATA_ENUM_TYPE ColumnBegin = EOL;
			Parent[Li] = ENUMUNDEF;
			INMOST_DATA_ENUM_TYPE PathEnd = ENUMUNDEF;
			INMOST_DATA_ENUM_TYPE Trace = k;
			INMOST_DATA_REAL_TYPE ShortestPath = 0;
			INMOST_DATA_REAL_TYPE AugmentPath = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
			while(true)
			{
				int AthrLi = A_Address[Li].thr;
				for (INMOST_DATA_ENUM_TYPE Lit = A_Address[Li].first; Lit < A_Address[Li].last; ++Lit)
				{
					INMOST_DATA_ENUM_TYPE Ui = A_Entries[AthrLi][Lit].first;
					//if( ColumnList[Ui] == k ) continue;
					if( ColumnList[Ui] != ENUMUNDEF ) continue;
					INMOST_DATA_ENUM_TYPE CLi = C_Row[Li - rbeg] + Lit - A_Address[Li].first;
					INMOST_DATA_REAL_TYPE l = fabs(ShortestPath + C_Entries[CLi] - V[Li] - U[Ui]);
					//if( l < 0.0 ) printf("row %d col %d negative l %g Augment %lf Shortest %lf C %lf V %lf U %lf\n",k,Ui,l,AugmentPath,ShortestPath,C_Entries[Lit],V[Li],U[Ui]);
					if( l < 0.0 && l > -1.0e-8 ) l = 0;
					//if( l < 0.0 ) continue;
					if (l < 0.0) std::cout << "row " << k << " col " << Ui << " negative l " << l << " augment " << AugmentPath << " ShortestPath " << ShortestPath << " C " << C_Entries[CLi] << " V " << V[Li] << " U " << U[Ui] << " l " << ShortestPath + C_Entries[CLi] - V[Li] - U[Ui] << std::endl;
					if( l < AugmentPath )
					{
						if( Perm[Ui] == ENUMUNDEF )
						{
							PathEnd = Ui;
							Trace = Li;
							AugmentPath = l;
							AugmentPosition[Ui] = CLi;//Lit;
						}
						else if( l < Dist[Ui] )
						{
							Parent[Perm[Ui]] = Li;
							AugmentPosition[Ui] = CLi;// Lit;
							if( Heap.Contains(Ui-cbeg) )
								Heap.DecreaseKey(Ui-cbeg,l);
							else
								Heap.PushHeap(Ui-cbeg,l);
						}
					}
				}

				INMOST_DATA_ENUM_TYPE pop_heap_pos = Heap.PopHeap();
				if( pop_heap_pos == ENUMUNDEF ) break;
			
				INMOST_DATA_ENUM_TYPE Ui = pop_heap_pos+cbeg;
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
				INMOST_DATA_ENUM_TYPE Ui = ColumnBegin;
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
					INMOST_DATA_ENUM_TYPE IPermPrev = IPerm[Trace];
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

#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
		for (INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
			DL[k] = (V[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() ? 1 : exp(V[k]));
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
		for (INMOST_DATA_INTEGER_TYPE k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k)
			DR[k] = (U[k] == std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() || check_zero(Cmax[k])  ? 1 : exp(U[k]) / Cmax[k]);

		int err = 0;
#if defined(USE_OMP_FACT)
#pragma omp parallel for reduction(+:err)
#endif
		for (INMOST_DATA_INTEGER_TYPE k = rbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(rend); ++k)
		{
			bool flip_sign = false;
			int Athr = A_Address[k].thr;
			for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
			{
				INMOST_DATA_ENUM_TYPE i = A_Entries[Athr][jt].first;
				INMOST_DATA_ENUM_TYPE j = Perm[i];				
				if( fabs(DL[k]*A_Entries[Athr][jt].second*DR[i]) > 1 + 1.0e-7 )
				{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
					std::cout << "element on row " << k << " col " << A_Entries[Athr][jt].first << " value " << A_Entries[Athr][jt].second << " u " << DR[i] << " l " << DL[k] << " U " << U[i] << " V " << V[k] << " Cmax " << Cmax[i] << " scaled " << DL[k]*A_Entries[Athr][jt].second*DR[i] << std::endl;
					err++;
					//exit(-1);
				}
				if( static_cast<INMOST_DATA_INTEGER_TYPE>(j) == k )
				{
					if( DL[k]*A_Entries[Athr][jt].second*DR[i] < 0.0 ) flip_sign = true;
				}
			}

			if( flip_sign ) DL[k] *= -1;
		}


		if (err) std::cout << "Total scaling errors: " << err << std::endl;
		
		
	}

	void MLMTILUC_preconditioner::FillColumnGaps(	INMOST_DATA_ENUM_TYPE cbeg,
													INMOST_DATA_ENUM_TYPE cend,
													interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ)
	{
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> tmp(cbeg, cend, ENUMUNDEF);

	}

	void MLMTILUC_preconditioner::GraphWeights(	INMOST_DATA_ENUM_TYPE wbeg,
												INMOST_DATA_ENUM_TYPE wend,
												const interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
												const std::vector< std::vector<Sparse::Row::entry> >& A_Entries,
												const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR,
												std::vector<INMOST_DATA_REAL_TYPE>& wadj)
	{
		for (INMOST_DATA_ENUM_TYPE i = wbeg; i < wend; ++i)
		{
			for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
			{
				INMOST_DATA_ENUM_TYPE j = A_Entries[A_Address[i].thr][jt].first;
				if (j != i && j >= wbeg && j < wend && !Pivot[j])
				{
					INMOST_DATA_REAL_TYPE v = fabs(DL[i] * A_Entries[A_Address[i].thr][jt].second * DR[j]);
					wadj[i - wbeg] += v;
					wadj[j - wbeg] += v;
				}
			}
		}
	}

	void MLMTILUC_preconditioner::CheckBlock(const Block& b,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
		const std::vector< std::vector<Sparse::Row::entry> >& A_Entries,
		INMOST_DATA_ENUM_TYPE sepbeg,
		INMOST_DATA_ENUM_TYPE sepend,
		std::string file, int line)
	{

#ifndef NDEBUG
		INMOST_DATA_ENUM_TYPE inblock = 0, insep = 0, outside = 0;
		for (INMOST_DATA_ENUM_TYPE k = b.row_start; k < b.row_end; ++k)
		{
			for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE i = A_Entries[A_Address[k].thr][it].first;
				if (i >= b.col_start && i < b.col_end)
					inblock++;
				else if (i >= sepbeg && i < sepend)
					insep++;
				else
				{
					std::cout << "row " << k << " index " << i << " outside block " << b.col_start << ":" << b.col_end << " and separator " << sepbeg << ":" << sepend << std::endl;
					outside++;
				}
			}
		}
		//if (outside) 
			std::cout << file << ":" << line << " block rows " << b.row_start << ":" << b.row_end << " cols " << b.col_start << ":" << b.col_end << " inside " << inblock << " separator " << insep << " outside " << outside << std::endl;
		assert(!outside);
#endif
	}
	void MLMTILUC_preconditioner::CheckBlock(const Block& b,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
		const std::vector< std::vector<Sparse::Row::entry> >& A_Entries,
		const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
		INMOST_DATA_ENUM_TYPE sepbeg,
		INMOST_DATA_ENUM_TYPE sepend,
		std::string file, int line)
	{
#ifndef NDEBUG
		INMOST_DATA_ENUM_TYPE inblock = 0, insep = 0, outside = 0;
		for (INMOST_DATA_ENUM_TYPE k = b.row_start; k < b.row_end; ++k)
		{
			for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
			{
				INMOST_DATA_ENUM_TYPE i = localQ[A_Entries[A_Address[k].thr][it].first];
				if (i >= b.col_start && i < b.col_end)
					inblock++;
				else if (i >= sepbeg && i < sepend)
					insep++;
				else
				{
					std::cout << "row " << k << " index " << i << " outside block " << b.col_start << ":" << b.col_end << " and separator " << sepbeg << ":" << sepend << std::endl;
					outside++;
				}
			}
		}
		//if (outside) 
		std::cout << file << ":" << line << " block rows " << b.row_start << ":" << b.row_end << " cols " << b.col_start << ":" << b.col_end << " inside " << inblock << " separator " << insep << " outside " << outside << std::endl;
		assert(!outside);
#endif
	}

	void MLMTILUC_preconditioner::SymmetricGraph(INMOST_DATA_ENUM_TYPE wbeg,
												INMOST_DATA_ENUM_TYPE wend,
												const interval<INMOST_DATA_ENUM_TYPE, Interval> & A_Address,
												const std::vector< std::vector<Sparse::Row::entry> > & A_Entries,
												const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
												std::vector<INMOST_DATA_ENUM_TYPE> & xadj,
												std::vector< INMOST_DATA_ENUM_TYPE > & adjncy)
	{
		xadj.resize(wend-wbeg+1);
		{//graph assembly
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Apos(wbeg, wend), Anext(wbeg, wend), Afirst(wbeg,wend);
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
				if (Pivot[k]) std::cout << "pivot row " << k << " appears in graph " << std::endl;
			for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
			{
				Apos[k] = EOL;
				Afirst[k] = EOL;
				Anext[k] = EOL;
			}
			for (INMOST_DATA_INTEGER_TYPE k = wend; k > static_cast<INMOST_DATA_INTEGER_TYPE>(wbeg); --k)
			{
				if (A_Address[k-1].Size() > 0)
				{
					INMOST_DATA_ENUM_TYPE beg = A_Address[k - 1].first;
					while (beg < A_Address[k - 1].last && Pivot[A_Entries[A_Address[k - 1].thr][beg].first]) beg++;
					if (beg < A_Address[k - 1].last)
					{
						INMOST_DATA_ENUM_TYPE i = A_Entries[A_Address[k - 1].thr][beg].first;
						Anext[k - 1] = Afirst[i];
						Afirst[i] = k - 1;
						Apos[k - 1] = beg;
					}
				}
			}
			xadj[0] = 0;
			for(INMOST_DATA_ENUM_TYPE i = wbeg; i < wend; ++i)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[i].first; jt < A_Address[i].last; ++jt)
				{
					INMOST_DATA_ENUM_TYPE j = A_Entries[A_Address[i].thr][jt].first;
					if( j != i && !Pivot[j] ) adjncy.push_back(j);
				}

				INMOST_DATA_ENUM_TYPE Li = Afirst[i];
				while (Li != EOL)
				{
					if( Li != i && !Pivot[Li] ) adjncy.push_back(Li);
					Li = Anext[Li];
				}

				Li = Afirst[i];
				while (Li != EOL)
				{
					INMOST_DATA_ENUM_TYPE Ui = Anext[Li];
					do
					{
						Apos[Li]++;
					} while (Apos[Li] < A_Address[Li].last && Pivot[A_Entries[A_Address[Li].thr][Apos[Li]].first]);
					if (A_Address[Li].Size() && A_Address[Li].last - Apos[Li] > 0)
					{
						INMOST_DATA_ENUM_TYPE k = A_Entries[A_Address[Li].thr][Apos[Li]].first;
						Anext[Li] = Afirst[k];
						Afirst[k] = Li;
					}

					Li = Ui;
				}

				std::sort(adjncy.begin()+xadj[i-wbeg],adjncy.end());
				adjncy.resize(std::unique(adjncy.begin()+xadj[i-wbeg],adjncy.end())-adjncy.begin());

				xadj[i-wbeg+1] = static_cast<INMOST_DATA_ENUM_TYPE>(adjncy.size());
			}
			
			//~ if (verbosity > 1&& print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		}
	}
	
	void MLMTILUC_preconditioner::ReorderColumns(	INMOST_DATA_ENUM_TYPE wbeg,
													INMOST_DATA_ENUM_TYPE wend,
													interval<INMOST_DATA_ENUM_TYPE, Interval> & A_Address,
													std::vector< std::vector<Sparse::Row::entry> > & A_Entries,
													interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ,
													interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> & DR,
													Sparse::Vector & DR0)
	{
		//for debug
#if !defined(NDEBUG)
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) 
			if( localQ[k] == ENUMUNDEF ) 
				std::cout << __FILE__ << ":" << __LINE__ << " No column permutation for row " << k << ". Matrix is structurally singular\n";
#endif
		{
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(wbeg, wend);
#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
			{
#if defined(USE_OMP_FACT)
#pragma omp for 
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					V[localQ[k]] = DR[k]; //if you expirience a bug here then the matrix is structurally singular
#if defined(USE_OMP_FACT)
#pragma omp for 
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					DR[k] = V[k];
#if defined(USE_OMP_FACT)
#pragma omp for 
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					V[localQ[k]] = DR0[k]; //if you expirience a bug here then the matrix is structurally singular
#if defined(USE_OMP_FACT)
#pragma omp for 
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
					DR0[k] = V[k];


#if defined(USE_OMP_FACT)
#pragma omp for 
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					int Athr = A_Address[k].thr;
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt < A_Address[k].last; ++jt)
						A_Entries[Athr][jt].first = localQ[A_Entries[Athr][jt].first];
					std::sort(A_Entries[Athr].begin() + A_Address[k].first, A_Entries[Athr].begin() + A_Address[k].last);
				}
			}
		}
		{
			interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(wbeg, wend, false);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(wbeg, wend), invQ(wbeg, wend);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> localP(wbeg,wend);
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				localP[k] = k;
			ReorderEF(wbeg, wend, donePQ, localP, localQ);
			inversePQ(wbeg, wend, localP, localQ, invP, invQ);
			applyPQ(wbeg, wend, localP, localQ, invP, invQ);
		}
	}


	void MLMTILUC_preconditioner::ReorderSystem(INMOST_DATA_ENUM_TYPE wbeg,
												INMOST_DATA_ENUM_TYPE wend,
												interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
												std::vector< std::vector<Sparse::Row::entry> >& A_Entries,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localP,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& localQ,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR,
												Sparse::Vector& DL0,
												Sparse::Vector& DR0)
	{
		{
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(wbeg, wend), V(wbeg, wend);
#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
			{
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					U[localP[k]] = DL[k];
					V[localQ[k]] = DR[k];
				}
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					DL[k] = U[k];
					DR[k] = V[k];
				}
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					U[localP[k]] = DL0[k];
					V[localQ[k]] = DR0[k];
				}
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					DL0[k] = U[k];
					DR0[k] = V[k];
				}
			}
		}
		{
			interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(wbeg, wend, false);
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(wbeg, wend), invQ(wbeg, wend);
			ReorderEF(wbeg, wend, donePQ, localP, localQ);//reorder E,F blocks by swaps
			inversePQ(wbeg, wend, localP, localQ, invP, invQ);//in localPQ numbers indicate where to put current row/column
			{
				interval<INMOST_DATA_ENUM_TYPE, Interval> B_Address(A_Address.get_interval_beg(), A_Address.get_interval_end());
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				{
					int Athr = A_Address[invP[k]].thr;
					B_Address[k].thr = Athr;
					B_Address[k].first = A_Address[invP[k]].first;
					for (INMOST_DATA_ENUM_TYPE jt = A_Address[invP[k]].first; jt < A_Address[invP[k]].last; ++jt)
						A_Entries[Athr][jt].first = localQ[A_Entries[Athr][jt].first];
					B_Address[k].last = A_Address[invP[k]].last;
				}
				//B_Entries.push_back(Sparse::Row::make_entry(-1, 0.0));
				A_Address.swap(B_Address);
				//A_Entries.swap(B_Entries);
			}
#if defined(USE_OMP_FACT)
#pragma omp parallel for 
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				std::sort(A_Entries[A_Address[k].thr].begin() + A_Address[k].first, A_Entries[A_Address[k].thr].begin() + A_Address[k].last);


			applyPQ(wbeg, wend, localP, localQ, invP, invQ); //reset localPQ
		}
	}

	void MLMTILUC_preconditioner::IdominanceScaling(INMOST_DATA_ENUM_TYPE wbeg,
													INMOST_DATA_ENUM_TYPE wend,
													const interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
													std::vector< std::vector<Sparse::Row::entry> >& A_Entries,
													const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
													interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
													interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR)
	{
		{/// THIS VERSION OF RESCALING INCREASES DIAGONAL DOMINANCE
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(wbeg, wend), V(wbeg, wend);
			std::vector< std::vector<INMOST_DATA_REAL_TYPE> > C_Entries(A_Entries.size());
			interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> temp(wbeg, wend, 0);
#if defined(USE_OMP)
			bool parallel = !omp_in_parallel();
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
				if (Pivot[k]) std::cout << "pivot row " << k << " appears in scaling " << std::endl;
			for (size_t q = 0; q < A_Entries.size(); ++q)
				C_Entries[q].resize(A_Entries[q].size());
#if defined(USE_OMP_FACT)
#pragma omp parallel for if(parallel)
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if( !Pivot[A_Entries[A_Address[k].thr][r].first] )
					A_Entries[A_Address[k].thr][r].second *= (DL[k] * DR[A_Entries[A_Address[k].thr][r].first]);

				U[k] = DL[k];
				V[k] = DR[k];
			}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if(parallel)
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
			{
				for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (!Pivot[A_Entries[A_Address[k].thr][r].first])
				{
					INMOST_DATA_REAL_TYPE u = fabs(A_Entries[A_Address[k].thr][r].second);
					if (check_zero(u))
						C_Entries[A_Address[k].thr][r] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
					else
						C_Entries[A_Address[k].thr][r] = -log(u);
				}
			}
			for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
			{
#if defined(USE_OMP_FACT)
#pragma omp parallel for if(parallel)
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) 
					DL[k] = DR[k] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
//needs atomic for std::min
//#if defined(USE_OMP_FACT)
//#pragma omp parallel for if(parallel)
//#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++) //row number
				{
					for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (!Pivot[A_Entries[A_Address[k].thr][r].first])
					{
						INMOST_DATA_ENUM_TYPE i = A_Entries[A_Address[k].thr][r].first; //column number
						if (static_cast<INMOST_DATA_INTEGER_TYPE>(i) != k) //out of diagonal
						{
							INMOST_DATA_REAL_TYPE u = C_Entries[A_Address[k].thr][r] + temp[k] - temp[i];
							DL[k] = std::min(DL[k], u);// update Y1
							DR[i] = std::min(DR[i], u);// update Y2
						}
					}
				}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if(parallel)
#endif
				for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
				{
					if (DR[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max() &&
						DL[k] != std::numeric_limits<INMOST_DATA_REAL_TYPE>::max())
						temp[k] += (DR[k] - DL[k]) * 0.5;
				}
			}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if(parallel)
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			{
				DL[k] = exp(-temp[k]);
				DR[k] = exp(temp[k]);
				if (DL[k] != DL[k]) DL[k] = 1;
				if (DR[k] != DR[k]) DR[k] = 1;
			}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if(parallel)
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (!Pivot[A_Entries[A_Address[k].thr][r].first])
					A_Entries[A_Address[k].thr][r].second *= (DL[k] * DR[A_Entries[A_Address[k].thr][r].first]);
			}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if(parallel)
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			{
				DL[k] *= U[k];
				DR[k] *= V[k];
			}
		}
	}

	void MLMTILUC_preconditioner::SinkhornScaling(INMOST_DATA_ENUM_TYPE wbeg,
		INMOST_DATA_ENUM_TYPE wend,
		interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
		std::vector<Sparse::Row::entry>& A_Entries,
		interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR,
		INMOST_DATA_REAL_TYPE p)
	{
		/// ROW-COLUMN ALTERNATING SCALING FOR p-NORM BALANCING
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(wbeg, wend), V(wbeg, wend);
#if defined(USE_OMP)
		bool parallel = !omp_in_parallel();
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			if (Pivot[k]) std::cout << "pivot row " << k << " appears in scaling " << std::endl;
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
		{
			for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (!Pivot[A_Entries[r].first])
				A_Entries[r].second *= (DL[k] * DR[A_Entries[r].first]);

			U[k] = DL[k];
			V[k] = DR[k];
		}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
			DL[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
		{
			for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (!Pivot[A_Entries[r].first])
				DL[k] += std::pow(std::fabs(A_Entries[r].second), p);
		}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			DL[k] = 1.0 / std::max(eps, DL[k]);
		for (INMOST_DATA_ENUM_TYPE iter = 0; iter < sciters; iter++)
		{
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				DR[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (!Pivot[A_Entries[r].first])
				{
					INMOST_DATA_REAL_TYPE u = DL[k] * std::pow(std::fabs(A_Entries[r].second), p);
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
					DR[A_Entries[r].first] += u;
				}
			}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
				DR[k] = 1.0 / std::max(eps, DR[k]);
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
				DL[k] = 0.0;
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
			{
				for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (!Pivot[A_Entries[r].first])
					DL[k] += DR[A_Entries[r].first] * std::pow(std::fabs(A_Entries[r].second), p);
			}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
			for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
				DL[k] = 1.0 / std::max(eps, DL[k]);
		}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
		{
			for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r)
				A_Entries[r].second *= DL[k] * DR[A_Entries[r].first];
		}
#if defined(USE_OMP_FACT)
#pragma omp parallel for if( parallel )
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
		{
			DL[k] *= U[k];
			DR[k] *= V[k];
		}
		/// RESCALING DONE
	}

	void MLMTILUC_preconditioner::PivotScaling(	INMOST_DATA_ENUM_TYPE wbeg,
												INMOST_DATA_ENUM_TYPE wend,
												const interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
												std::vector< std::vector<Sparse::Row::entry> > & A_Entries,
												const interval<INMOST_DATA_ENUM_TYPE, bool>& Pivot,
												const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
												const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR)
	{
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
		for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++) 
		{
			for (INMOST_DATA_ENUM_TYPE r = A_Address[k].first; r < A_Address[k].last; ++r) if (Pivot[k] || Pivot[A_Entries[A_Address[k].thr][r].first])
				A_Entries[A_Address[k].thr][r].second *= (DL[k] * DR[A_Entries[A_Address[k].thr][r].first]);
		}
	}

	void MLMTILUC_preconditioner::EFScaling(INMOST_DATA_ENUM_TYPE wbeg,
											INMOST_DATA_ENUM_TYPE wend,
											INMOST_DATA_ENUM_TYPE mobeg,
											INMOST_DATA_ENUM_TYPE moend,
											const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DL,
											const interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& DR)
	{
		//rescale EF
		if (!level_size.empty())
		{
			INMOST_DATA_ENUM_TYPE first = mobeg, last;
			INMOST_DATA_ENUM_TYPE kbeg, kend;
			for (size_t level = 0; level < level_size.size(); ++level)
			{
				last = first + level_size[level];
				kbeg = std::max(wbeg, last);
				kend = std::min(wend, moend);
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
				for (INMOST_DATA_INTEGER_TYPE k = kbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(kend); ++k)
				{
					if (F_Address[level]->at(k).Size())
					{
						int Fthr = F_Address[level]->at(k).thr;
						for (INMOST_DATA_ENUM_TYPE r = F_Address[level]->at(k).first; r < F_Address[level]->at(k).last; ++r)
							F_Entries[Fthr][r].second *= DR[k];
					}
					if (E_Address[level]->at(k).Size())
					{
						int Ethr = E_Address[level]->at(k).thr;
						for (INMOST_DATA_ENUM_TYPE r = E_Address[level]->at(k).first; r < E_Address[level]->at(k).last; ++r)
							E_Entries[Ethr][r].second *= DL[k];
					}
				}
				first = last;
			}
		}
		//End rescale B block
	}
	
	void MLMTILUC_preconditioner::WRCMOrdering(	INMOST_DATA_ENUM_TYPE wbeg,
												INMOST_DATA_ENUM_TYPE wend,
												std::vector<INMOST_DATA_ENUM_TYPE> & xadj, 
												std::vector<INMOST_DATA_ENUM_TYPE> & adjncy,
												std::vector<INMOST_DATA_REAL_TYPE> & wadj,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
	{
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> order(wbeg,wend,ENUMUNDEF);
		//find node with the lowest order
		INMOST_DATA_ENUM_TYPE index = wbeg;
		INMOST_DATA_ENUM_TYPE cur = ENUMUNDEF;
		std::deque<INMOST_DATA_ENUM_TYPE> q;
		//~ std::vector< std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> > conns;
		std::vector< INMOST_DATA_ENUM_TYPE > conns;
		//enumerate elements without connections
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) if (order[k] == ENUMUNDEF)
		{
			if (xadj[k + 1 - wbeg] - xadj[k - wbeg] == 0)
				order[k] = index++;
		}
		while (index < wend)
		{
			cur = ENUMUNDEF;
			for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend) && cur == ENUMUNDEF; ++k)
			{
				if( order[k] == ENUMUNDEF )
					cur = k;
			}
			if( cur == ENUMUNDEF ) break;
			assert(cur != ENUMUNDEF);
			for(INMOST_DATA_INTEGER_TYPE k = cur+1; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k) if( order[k] == ENUMUNDEF )
			{
				if( WRCM_Comparator(wbeg,wadj)(k,cur) )
					cur = k;
			}
			q.push_back(cur);
			order[cur] = index++;
			bool finish_block = false;
			while(!q.empty())
			{
				cur = q.front();
				q.pop_front();
				if( finish_block )
				{
					for (INMOST_DATA_ENUM_TYPE it = xadj[cur-wbeg]; it < xadj[cur-wbeg+1]; ++it)
						if( order[adjncy[it]] == ENUMUNDEF )
							order[adjncy[it]] = ENUMUNDEF-1;
				}
				else
				{
					for (INMOST_DATA_ENUM_TYPE it = xadj[cur-wbeg]; it < xadj[cur-wbeg+1]; ++it)
						if( order[adjncy[it]] == ENUMUNDEF )
							conns.push_back(adjncy[it]);
					std::sort(conns.begin(),conns.end(),WRCM_Comparator(wbeg,wadj));
					for (INMOST_DATA_INTEGER_TYPE k = 0; k < static_cast<INMOST_DATA_INTEGER_TYPE>(conns.size()); ++k)
					{
						order[conns[k]] = index++;
						q.push_back(conns[k]);
						//~ if( (index-wbeg)%512==0 ) finish_block = true;
					}
					conns.clear();
				}
			}
			
		}
		
		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
			if( order[k] == ENUMUNDEF-1 ) order[k] = index++;

		//reverse
		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
			order[k] = wend-(order[k]-wbeg)-1;

		for(INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
		{
			localP[k] = order[k];
			localQ[k] = order[k];
		}		
	}
	
	void MLMTILUC_preconditioner::RCMOrdering(INMOST_DATA_ENUM_TYPE wbeg,
												INMOST_DATA_ENUM_TYPE wend,
												std::vector<INMOST_DATA_ENUM_TYPE> & xadj, 
												std::vector<INMOST_DATA_ENUM_TYPE> & adjncy,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
	{
		interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> order(wbeg,wend,ENUMUNDEF);
		//find node with the lowest order
		INMOST_DATA_ENUM_TYPE index = wbeg;
		INMOST_DATA_ENUM_TYPE cur = ENUMUNDEF;
		std::deque<INMOST_DATA_ENUM_TYPE> q;
		std::vector<INMOST_DATA_ENUM_TYPE> conns;
		//enumerate elements without connections
		for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) if (order[k] == ENUMUNDEF)
		{
			if (xadj[k + 1 - wbeg] - xadj[k - wbeg] == 0)
				order[k] = index++;
		}
		while (index < wend)
		{
			cur = ENUMUNDEF;
			for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend && cur == ENUMUNDEF; ++k)
			{
				if( order[k] == ENUMUNDEF )
					cur = k;
			}
			assert(cur != ENUMUNDEF);
			for(INMOST_DATA_ENUM_TYPE k = cur+1; k < wend; ++k) if( order[k] == ENUMUNDEF )
			{
				if( RCM_Comparator(wbeg,xadj)(k,cur) )
					cur = k;
			}
			q.push_back(cur);
			order[cur] = index++;
			while(!q.empty())
			{
				cur = q.front();
				q.pop_front();
				for (INMOST_DATA_ENUM_TYPE it = xadj[cur-wbeg]; it < xadj[cur-wbeg+1]; ++it)
					if( order[adjncy[it]] == ENUMUNDEF )
						conns.push_back(adjncy[it]);
				std::sort(conns.begin(),conns.end(),RCM_Comparator(wbeg,xadj));
				for (INMOST_DATA_ENUM_TYPE k = 0; k < static_cast<INMOST_DATA_ENUM_TYPE>(conns.size()); ++k)
				{
					order[conns[k]] = index++;
					q.push_back(conns[k]);
				}
				conns.clear();
			}
		}
		//reverse
		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
			order[k] = wend-(order[k]-wbeg)-1;

		for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
		{
			localP[k] = order[k];
			localQ[k] = order[k];
		}
	}

	void MLMTILUC_preconditioner::CheckColumnGaps(const Block& b,
		const interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
		const std::vector< std::vector<Sparse::Row::entry> >& A_Entries)
	{
#ifndef NDEBUG
		INMOST_DATA_ENUM_TYPE rbeg = b.row_start, rend = b.row_end;
		INMOST_DATA_ENUM_TYPE cbeg = b.col_start, cend = b.col_end;
		INMOST_DATA_ENUM_TYPE nocol = 0, norow = 0;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> count(cbeg, cend, 0);
		for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
		{
			INMOST_DATA_ENUM_TYPE nrow = 0;
			for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt != A_Address[k].last; ++jt)// if (!check_zero(A_Entries[jt].second))
			{
				nrow++;
				count[A_Entries[A_Address[k].thr][jt].first]++;
			}
			if (nrow == 0)
			{
				//std::cout << "Empty row at " << k << std::endl;
				norow++;
			}
		}
		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
			if (count[k] == 0)
			{
				//std::cout << "No column index at " << k << std::endl;
				nocol++;
			}
		if( norow || nocol ) std::cout << "  ### No indices at " << nocol << " columns and " << norow << " rows" << std::endl;
#endif //NDEBUG
	}

	void MLMTILUC_preconditioner::Factorize(INMOST_DATA_ENUM_TYPE cbeg,
											INMOST_DATA_ENUM_TYPE cend,
											interval<INMOST_DATA_ENUM_TYPE, Interval>& A_Address,
											std::vector< std::vector<Sparse::Row::entry> > & A_Entries,
											interval<INMOST_DATA_ENUM_TYPE, Interval>& L2_Address,
											std::vector< std::vector<Sparse::Row::entry> > & L2_Entries,
											interval<INMOST_DATA_ENUM_TYPE, Interval>& U2_Address,
											std::vector< std::vector<Sparse::Row::entry> > & U2_Entries,
											interval<INMOST_DATA_ENUM_TYPE, bool> & Pivot,
											bool block_pivot,
											INMOST_DATA_REAL_TYPE& NuLout,
											INMOST_DATA_REAL_TYPE& NuUout,
											double& testimator)
	{
		double tlocal;
		INMOST_DATA_ENUM_TYPE swaps = 0, ndrops = 0, nzLU = 0, nzLU2 = 0;
#if defined(ILUC2)
		INMOST_DATA_REAL_TYPE tau2 = iluc2_tau;
#else
		INMOST_DATA_REAL_TYPE tau2 = tau;
#endif
#if defined(ILUC2)
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2list(cbeg, cend, UNDEF), L2list(cbeg, cend, UNDEF);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2beg(cbeg, cend, EOL), L2beg(cbeg, cend, EOL);
		//L2_Address.set_interval_beg(cbeg);
		//L2_Address.set_interval_end(cend + 1);
		//U2_Address.set_interval_beg(cbeg);
		//U2_Address.set_interval_end(cend + 1);
#endif
#if defined(ILUC2)
		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
		{
			L2_Address[k].first = L2_Address[k].last = ENUMUNDEF;
			U2_Address[k].first = U2_Address[k].last = ENUMUNDEF;
		}
#endif
		//Setup column addressing for B,F, in descending order to keep list ordered
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Anext(cbeg, cend), Afirst(cbeg, cend, EOL);
		//supplimentary data structures for condition estimates of L^{-1}, U^{-1}
		INMOST_DATA_REAL_TYPE mu_update_L1, mu_update_L2, mu_update_U1, mu_update_U2;
		INMOST_DATA_REAL_TYPE NuU_old = 1, NuL_old = 1, NuD_old = 1, NuU_max = 1, NuL_max = 1;
		INMOST_DATA_REAL_TYPE NuU = 1, NuL = 1, NuD = 1, max_diag = 0, min_diag = 1.0e+300;
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> LineValuesU(cbeg, cend, 0.0), LineValuesL(cbeg, cend, 0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> LineIndecesU(cbeg, cend + 1, UNDEF), LineIndecesL(cbeg, cend + 1, UNDEF);
		std::vector<INMOST_DATA_ENUM_TYPE> indicesU, indicesL;
#if defined(PREMATURE_DROPPING)
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Unorms(cbeg, cend, 0.0), Lnorms(cbeg, cend, 0.0);
#endif
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL1(cbeg, cend, 0.0), EstU1(cbeg, cend, 0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> EstL2(cbeg, cend, 0.0), EstU2(cbeg, cend, 0.0);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Apos(cbeg, cend, ENUMUNDEF), Upos(cbeg, cend, ENUMUNDEF), Lpos(cbeg, cend, ENUMUNDEF);
#if defined(ILUC2)
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2pos(cbeg, cend, ENUMUNDEF), L2pos(cbeg, cend, ENUMUNDEF);
#endif
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Llist(cbeg, cend), Ulist(cbeg, cend);
		interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Ubeg(cbeg, cend, EOL), Lbeg(cbeg, cend, EOL);
		int thr = Thread();
		//if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		// //~ Clear TransposeB
		//~ PrepareTranspose(cbeg,cend,A_Address,A_Entries,TransposeA);
		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
			Apos[k] = A_Address[k].first;
		for (INMOST_DATA_ENUM_TYPE k = cend; k > cbeg; --k)
		{
			if (A_Address[k - 1].Size() > 0)
			{
				INMOST_DATA_ENUM_TYPE beg = A_Address[k - 1].first;
				while (beg < A_Address[k - 1].last && Pivot[A_Entries[A_Address[k - 1].thr][beg].first]) beg++;
				if (beg < A_Address[k - 1].last)
				{
					INMOST_DATA_ENUM_TYPE i = A_Entries[A_Address[k - 1].thr][beg].first;
					if (i < k - 1)
					{
						Anext[k - 1] = Afirst[i];
						Afirst[i] = k - 1;
					}
				}
			}
		}
		//setup diagonal values and memorize indices
		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; k++)
		{
			LU_Diag[k] = 0.0;
			for (INMOST_DATA_ENUM_TYPE i = A_Address[k].first; i < A_Address[k].last; i++)
			{
				if (A_Entries[A_Address[k].thr][i].first == k)
				{
					LU_Diag[k] = A_Entries[A_Address[k].thr][i].second;
					break;
				}
			}
		}
		if (verbosity > 1)
		{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
			std::cout << " starting factorization " << cbeg << "<->" << cend << " thread " << thr << std::endl;
		}

		int report_pace = std::max<int>((cend - cbeg) / 100,1);

		if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
		{
			if (!Pivot[k])
			{
//#if defined(USE_OMP_FACT)
//#pragma omp parallel sections num_threads(2)
//#endif
				{
					// Compute U-part
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
					{
						// Prepare linked list for U-part
						INMOST_DATA_ENUM_TYPE Sbeg, i;
						INMOST_DATA_REAL_TYPE coef;
						//uncompress k-th row
						// add diagonal value first, there shouldn't be values on left from diagonal
						//assert(B_Entries[A_Address[k].first].first == k);
						Sbeg = k;
						LineIndecesU[k] = EOL;
						LineValuesU[k] = 0.0;
						if (A_Address[k].Size())
						{
							i = k;
							for (INMOST_DATA_ENUM_TYPE it = A_Address[k].first; it < A_Address[k].last; ++it)
							{
								assert(A_Entries[A_Address[k].thr][it].first >= k);
								if (!Pivot[A_Entries[A_Address[k].thr][it].first])
								{
									LineValuesU[A_Entries[A_Address[k].thr][it].first] = A_Entries[A_Address[k].thr][it].second;
									i = LineIndecesU[i] = A_Entries[A_Address[k].thr][it].first;
								}
							}
							LineIndecesU[i] = EOL;
						}
						//diagonal perturbation
#if defined(DIAGONAL_PERTURBATION)
						LineValuesU[k] = LineValuesU[k] * (1.0 + rpert) + (LineValuesU[k] < 0.0 ? -1.0 : 1.0) * apert;
#endif
						// U-part elimination with L
						i = Lbeg[k];
						while (i != EOL)
						{
							assert(i != UNDEF);
							assert(static_cast<INMOST_DATA_INTEGER_TYPE>(L_Entries[L_Address[i].thr][L_Address[i].first].first) == k);
							coef = -L_Entries[L_Address[i].thr][L_Address[i].first].second * LU_Diag[i];
#if defined(PREMATURE_DROPPING)
							if (std::fabs(coef) * Unorms[i] > tau2)
#endif
							{
								AddListUnordered(Sbeg, U_Address[i], U_Entries[U_Address[i].thr], coef, LineIndecesU, LineValuesU, 0.0);
#if defined(ILUC2)
								AddListUnordered(Sbeg, U2_Address[i], U2_Entries[U2_Address[i].thr], coef, LineIndecesU, LineValuesU, 0.0);
#endif
							}
							i = Llist[i];
						}
#if defined(ILUC2)
						// U-part elimination with second-order L 
						i = L2beg[k];
						while (i != EOL)
						{
							assert(i != UNDEF);
							assert(static_cast<INMOST_DATA_INTEGER_TYPE>(L2_Entries[L2_Address[i].thr][L2_Address[i].first].first) == k);
							coef = -L2_Entries[L2_Address[i].thr][L2_Address[i].first].second * LU_Diag[i];
#if defined(PREMATURE_DROPPING)
							if (std::fabs(coef) * Unorms[i] > tau)
#endif
							{
								AddListUnordered(Sbeg, U_Address[i], U_Entries[U_Address[i].thr], coef, LineIndecesU, LineValuesU, 0.0);
#if 0
								AddListUnordered(Sbeg, U2_Address[i], U2_Entries[U2_Address[i].thr], coef, LineIndecesU, LineValuesU, 0.0);
#endif
							}
							i = L2list[i];
						}
#endif
						// Order contents of linked list
						OrderList(Sbeg, LineIndecesU, indicesU);
					}
					// Compute L-part
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
					{
						// prepearing linked list for L-part
						INMOST_DATA_ENUM_TYPE Sbeg, i, j;
						INMOST_DATA_REAL_TYPE coef;
						//uncompress column
						//insert diagonal value first
						Sbeg = k;
						LineIndecesL[k] = EOL;
						LineValuesL[k] = 0.0;
						if (A_Address[k].Size())
						{
							if (static_cast<INMOST_DATA_INTEGER_TYPE>(A_Entries[A_Address[k].thr][A_Address[k].first].first) == k)
								LineValuesL[k] = A_Entries[A_Address[k].thr][A_Address[k].first].second;
							else
								LineValuesL[k] = 0.0;
							//start from diagonal
							j = k;
							i = Afirst[k];
							while (i != EOL)
							{
								assert(static_cast<INMOST_DATA_INTEGER_TYPE>(A_Entries[A_Address[i].thr][A_Address[i].first].first) == k);
								if (!Pivot[i])
								{
									LineValuesL[i] = A_Entries[A_Address[i].thr][A_Address[i].first].second;
									j = LineIndecesL[j] = i;
								}
								i = Anext[i];
							}
							LineIndecesL[j] = EOL;
						}
						{// Shift indexes for transposed B matrix traversal
							INMOST_DATA_ENUM_TYPE Li = Afirst[k];
							while (Li != EOL)
							{
								INMOST_DATA_ENUM_TYPE Ui = Anext[Li];
								do
								{
									A_Address[Li].first++;
								} while (A_Address[Li].first < A_Address[Li].last && Pivot[A_Entries[A_Address[Li].thr][A_Address[Li].first].first]);
								if (A_Address[Li].Size() > 0)
								{
									INMOST_DATA_ENUM_TYPE i = A_Entries[A_Address[Li].thr][A_Address[Li].first].first, curr, next;
									if (i < Li)
									{
										if (Afirst[i] > Li)
										{
											Anext[Li] = Afirst[i];
											Afirst[i] = Li;
										}
										else
										{
											curr = next = Afirst[i];
											while (next < Li)
											{
												curr = next;
												next = Anext[curr];
											}
											assert(curr < Li);
											assert(Li < next);
											Anext[curr] = Li;
											Anext[Li] = next;
										}
									}
								}
								Li = Ui;
							}
						}
						//diagonal perturbation
#if defined(DIAGONAL_PERTURBATION)
						LineValuesL[k] = LineValuesL[k] * (1.0 + rpert) + (LineValuesU[k] < 0.0 ? -1.0 : 1.0) * apert;
#endif
						// L-part elimination with U
						i = Ubeg[k];
						while (i != EOL)
						{
							assert(i != UNDEF);
							assert(static_cast<INMOST_DATA_INTEGER_TYPE>(U_Entries[U_Address[i].thr][U_Address[i].first].first) == k);
							coef = -U_Entries[U_Address[i].thr][U_Address[i].first].second * LU_Diag[i];
#if defined(PREMATURE_DROPPING)
							if (std::fabs(coef) * Lnorms[i] > tau2)
#endif
							{
								AddListUnordered(Sbeg, L_Address[i], L_Entries[L_Address[i].thr], coef, LineIndecesL, LineValuesL, 0.0);
#if defined(ILUC2)
								AddListUnordered(Sbeg, L2_Address[i], L2_Entries[L2_Address[i].thr], coef, LineIndecesL, LineValuesL, 0.0);
#endif
							}
							i = Ulist[i];
						}
#if defined(ILUC2)
						// L-part elimination with second-order U
						i = U2beg[k];
						while (i != EOL)
						{
							assert(i != UNDEF);
							assert(static_cast<INMOST_DATA_INTEGER_TYPE>(U2_Entries[U2_Address[i].thr][U2_Address[i].first].first) == k);
							coef = -U2_Entries[U2_Address[i].thr][U2_Address[i].first].second * LU_Diag[i];
#if defined(PREMATURE_DROPPING)
							if (std::fabs(coef) * Lnorms[i] > tau)
#endif
							{
								AddListUnordered(Sbeg, L_Address[i], L_Entries[L_Address[i].thr], coef, LineIndecesL, LineValuesL, 0.0);
#if 0
								AddListUnordered(Sbeg, L2_Address[i], L2_Entries[L2_Address[i].thr], coef, LineIndecesL, LineValuesL, 0.0);
#endif
							}
							i = U2list[i];
						}
#endif
						// Order contents of linked list
						OrderList(Sbeg, LineIndecesL, indicesL);
					}
				}
				
				if( false ) if (fabs(LineValuesU[k]) < tau2 || fabs(LineValuesL[k]) < tau2)
				{
//#if defined(USE_OMP_FACT)
//#pragma omp critical
//#endif
					{
						std::cout << "Tiny diagonal " << k << ": U " << LineValuesU[k] << " L " << LineValuesL[k] << std::endl;
						/*
						INMOST_DATA_ENUM_TYPE Li = k;
						while (Li != EOL)
						{
							if( !check_zero(LineValuesU[Li]) )
								std::cout << "(" << Li << "," << LineValuesU[Li] << ") ";
							Li = LineIndecesL[Li];
						}
						std::cout << std::endl;
						Li = k;
						while (Li != EOL)
						{
							if (!check_zero(LineValuesL[Li]))
								std::cout << "(" << Li << "," << LineValuesL[Li] << ") ";
							Li = LineIndecesU[Li];
						}
						std::cout << std::endl;
						*/
					}
				}
				
				
#if defined(PIVOT_THRESHOLD)
				LineValuesU[k] = std::max(fabs(LineValuesU[k]), tau) * (LineValuesU[k] < 0 ? -1.0 : 1.0);
				LineValuesL[k] = std::max(fabs(LineValuesL[k]), tau) * (LineValuesL[k] < 0 ? -1.0 : 1.0);
#endif			
				if (block_pivot)
				{
					LineValuesU[k] = std::max(std::abs(LineValuesU[k]), tau) * (LineValuesU[k] < 0 ? -1.0 : 1.0);
					LineValuesL[k] = std::max(std::abs(LineValuesL[k]), tau) * (LineValuesL[k] < 0 ? -1.0 : 1.0);
				}
				LU_Diag[k] = (LineValuesU[k] + LineValuesL[k]) / 2.0;
				//if (std::fabs(LineValuesU[k]) < std::fabs(LineValuesL[k]))
					//LU_Diag[k] = LineValuesU[k];
				//else
					//LU_Diag[k] = LineValuesL[k];
				//LU_Diag[k] = std::max(fabs(LU_Diag[k]), tau) * (LU_Diag[k] < 0 ? -1.0 : 1.0);

				// Condition estimator for diagonal D
				NuD_old = NuD;
				NuD = std::max<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]), max_diag) / std::min<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]), min_diag);
//#if defined(USE_OMP_FACT)
//#pragma omp parallel sections num_threads(2)
//#endif
				{
					// Compute U-part
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
					{
						// Rescale with diagonal
						ScaleList(1.0 / LU_Diag[k], k, LineIndecesU, LineValuesU);
						// Condition estimator for U part
						if (estimator)
						{
							tlocal = Timer();
							NuU_old = NuU;
							NuU = Estimator1(k, LineIndecesU, LineValuesU, EstU1, mu_update_U1);
							NuU = std::max(NuU, Estimator2(k, LineIndecesU, LineValuesU, EstU2, mu_update_U2));
							testimator += Timer() - tlocal;
						} //if( estimator )
					}
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
					{
						// Rescale with diagonal
						ScaleList(1.0 / LU_Diag[k], k, LineIndecesL, LineValuesL);
						// Estimate condition for L^{-1}
						if (estimator)
						{
							tlocal = Timer();
							NuL_old = NuL;
							NuL = Estimator1(k, LineIndecesL, LineValuesL, EstL1, mu_update_L1);
							NuL = std::max(NuL, Estimator2(k, LineIndecesL, LineValuesL, EstL2, mu_update_L2));
							testimator += Timer() - tlocal;
						}
					}
				}
					
					//std::cout << "k " << k << " NuU " << NuU << " NuL " << NuL << " D " << LU_Diag[k] << "    " << std::endl;
					//discard line if condition number is too large
#if defined(CONDITION_PIVOT)
				if (allow_pivot && !block_pivot && (NuU > pivot_cond || NuL > pivot_cond || NuD > pivot_diag))
				{
					//std::cout << "swap " << k << std::endl;
					//restore condition number
					//EstU1[k] = EstL1[k] = 0;
					//EstU2[k] = EstL2[k] = 0;
					//Cancel this row
					Pivot[k] = true;
					swaps++;
				}
#endif
			}
			if (Pivot[k])
			{
				NuL = NuL_old;
				NuU = NuU_old;
				NuD = NuD_old;
				U_Address[k].thr = 0;
				Upos[k] = U_Address[k].first = U_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U_Entries[U_Address[k].thr].size());
				L_Address[k].thr = 0;
				Lpos[k] = L_Address[k].first = L_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L_Entries[L_Address[k].thr].size());
#if defined(ILUC2)
				U2_Address[k].thr = 0;
				U2pos[k] = U2_Address[k].first = U2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U2_Entries[U2_Address[k].thr].size());
				L2_Address[k].thr = 0;
				L2pos[k] = L2_Address[k].first = L2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L2_Entries[L2_Address[k].thr].size());
#endif
			}
			else
			{
				NuL_max = std::max(NuL, NuL_max);
				NuU_max = std::max(NuU, NuU_max);
				max_diag = std::max<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]), max_diag);
				min_diag = std::min<INMOST_DATA_REAL_TYPE>(fabs(LU_Diag[k]), min_diag);

				if (fabs(LineValuesU[k]) < tau2 || fabs(LineValuesL[k]) < tau2)
				{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
					{
						std::cout << "Tiny diagonal " << k << ": U " << LineValuesU[k] << " L " << LineValuesL[k] << " NuL " << NuL << " NuU " << NuU << " NuD " << NuD << " max_diag " << max_diag << " min_diag " << min_diag << std::endl;
						
						INMOST_DATA_ENUM_TYPE Li = k;
						while (Li != EOL)
						{
							if( !check_zero(LineValuesU[Li]) )
								std::cout << "(" << Li << "," << LineValuesU[Li] << ") ";
							Li = LineIndecesL[Li];
						}
						std::cout << std::endl;
						Li = k;
						while (Li != EOL)
						{
							if (!check_zero(LineValuesL[Li]))
								std::cout << "(" << Li << "," << LineValuesL[Li] << ") ";
							Li = LineIndecesU[Li];
						}
						std::cout << std::endl;
						
					}
				}
				// Update values on the whole diagonal with L and U
				//DiagonalUpdate(k, LU_Diag, LineIndecesL, LineValuesL, LineIndecesU, LineValuesU);
//#if defined(USE_OMP_FACT)
//#pragma omp parallel sections num_threads(2)
//#endif
				{
					// reconstruct U-part from linked list
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
					{
						INMOST_DATA_REAL_TYPE Unorm = 0;
						INMOST_DATA_ENUM_TYPE Ui = LineIndecesU[k];
						while (Ui != EOL)
						{
							Unorm += LineValuesU[Ui] * LineValuesU[Ui];
							Ui = LineIndecesU[Ui];
						}
						Unorm = sqrt(Unorm);
#if defined(PREMATURE_DROPPING)
						Unorms[k] = Unorm;
#endif
						U_Address[k].thr = thr;
						Upos[k] = U_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(U_Entries[U_Address[k].thr].size());
#if defined(ILUC2)
						U2_Address[k].thr = thr;
						U2pos[k] = U2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(U2_Entries[U2_Address[k].thr].size());
#endif
						assert(!(LineValuesU[k] != LineValuesU[k])); //check nan
						U_Entries[U_Address[k].thr].push_back(Sparse::Row::make_entry(k, LineValuesU[k]));
						U_Address[k].first++; //shift from diagonal
						Ui = LineIndecesU[k];
						while (Ui != EOL)
						{
							assert(!(LineValuesU[Ui] != LineValuesU[Ui])); //check nan
							INMOST_DATA_REAL_TYPE u = fabs(LineValuesU[Ui]);
							if (u * NuU > tau * Unorm) // apply dropping rule
								U_Entries[U_Address[k].thr].push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#if defined(ILUC2)
							else if (u * NuU > tau2 * Unorm)
								U2_Entries[U2_Address[k].thr].push_back(Sparse::Row::make_entry(Ui, LineValuesU[Ui]));
#endif
							else ndrops++;
							Ui = LineIndecesU[Ui];
						}
						U_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U_Entries[U_Address[k].thr].size());
#if defined(ILUC2)
						U2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(U2_Entries[U2_Address[k].thr].size());
#endif
						// Update estimator vectors
						if (estimator)
						{
							tlocal = Timer();
							//U-estimator
							EstimatorUpdate(k, LineIndecesU, LineValuesU, EstU1, mu_update_U1);
							EstimatorUpdate(k, LineIndecesU, LineValuesU, EstU2, mu_update_U2);
							testimator += Timer() - tlocal;
						}
					}
					// reconstruct L-part from linked list
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
					{
						INMOST_DATA_REAL_TYPE Lnorm = 0;
						INMOST_DATA_ENUM_TYPE Li = LineIndecesL[k];
						while (Li != EOL)
						{
							Lnorm += LineValuesL[Li] * LineValuesL[Li];
							Li = LineIndecesL[Li];
						}
						Lnorm = sqrt(Lnorm);
#if defined(PREMATURE_DROPPING)
						Lnorms[k] = Lnorm;
#endif
						//insert column to L part
						L_Address[k].thr = thr;
						Lpos[k] = L_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(L_Entries[L_Address[k].thr].size());
#if defined(ILUC2)
						L2_Address[k].thr = thr;
						L2pos[k] = L2_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(L2_Entries[L2_Address[k].thr].size());
#endif
						assert(!(LineValuesL[k] != LineValuesL[k])); //check nan
						L_Entries[L_Address[k].thr].push_back(Sparse::Row::make_entry(k, LineValuesL[k]));
						L_Address[k].first++;
						Li = LineIndecesL[k];
						while (Li != EOL)
						{
							assert(!(LineValuesL[Li] != LineValuesL[Li])); //check nan
							INMOST_DATA_REAL_TYPE u = fabs(LineValuesL[Li]);
							if (u * NuL > tau * Lnorm)
								L_Entries[L_Address[k].thr].push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#if defined(ILUC2)
							else if (u * NuL > tau2 * Lnorm)
								L2_Entries[L2_Address[k].thr].push_back(Sparse::Row::make_entry(Li, LineValuesL[Li]));
#endif
							else ndrops++;
							Li = LineIndecesL[Li];
						}
						L_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L_Entries[L_Address[k].thr].size());
#if defined(ILUC2)
						L2_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(L2_Entries[L2_Address[k].thr].size());
#endif
						// Update estimator vectors
						if (estimator)
						{
							tlocal = Timer();
							//L-estimator
							EstimatorUpdate(k, LineIndecesL, LineValuesL, EstL1, mu_update_L1);
							EstimatorUpdate(k, LineIndecesL, LineValuesL, EstL2, mu_update_L2);
							testimator += Timer() - tlocal;
						}
					}
				}
			}
//#if defined(USE_OMP_FACT)
//#pragma omp parallel sections num_threads(2)
//#endif
			{//clear and advance
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
				{//U-part
					// Cleanup structures
					ClearList(k, LineIndecesU, LineValuesU);
					// Insert indexes for transposed U-part traversal
					//~ PrepareTranspose(k,k+1,U_Address,U_Entries,TransposeU);
					if (U_Address[k].Size() != 0)
					{
						INMOST_DATA_ENUM_TYPE  i = U_Entries[U_Address[k].thr][U_Address[k].first].first;
						assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
						Ulist[k] = Ubeg[i];
						Ubeg[i] = k;
					}
					// Shift indexes for transposed U-part traversal for next iteration
					{
						INMOST_DATA_ENUM_TYPE i = Ubeg[k];
						while (i != EOL)
						{
							U_Address[i].first++;
							INMOST_DATA_ENUM_TYPE Li = Ulist[i], Ui;
							if (U_Address[i].Size() > 0)
							{
								Ui = U_Entries[U_Address[i].thr][U_Address[i].first].first;
								Ulist[i] = Ubeg[Ui];
								Ubeg[Ui] = i;
							}
							i = Li;
						}
					}
#if defined(ILUC2)
					// Insert indexes for transposed second-order U-part traversal
					//~ PrepareTranspose(k,k+1,U2_Address,U2_Entries,TransposeU2);
					if (U2_Address[k].Size() != 0)
					{
						INMOST_DATA_ENUM_TYPE i = U2_Entries[U2_Address[k].thr][U2_Address[k].first].first;
						assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
						U2list[k] = U2beg[i];
						U2beg[i] = k;
					}
					// Shift indexes for transposed second-order U-part traversal
					{
						INMOST_DATA_ENUM_TYPE i = U2beg[k];
						while (i != EOL)
						{
							U2_Address[i].first++;
							INMOST_DATA_ENUM_TYPE Li = U2list[i], Ui;
							if (U2_Address[i].Size() > 0)
							{
								Ui = U2_Entries[U2_Address[i].thr][U2_Address[i].first].first;
								U2list[i] = U2beg[Ui];
								U2beg[Ui] = i;
							}
							i = Li;
						}
					}
#endif
				}
//#if defined(USE_OMP_FACT)
//#pragma omp section
//#endif
				{//L-part
// Cleanup structures
					ClearList(k, LineIndecesL, LineValuesL);
					// Insert indexes for transposed L-part traversal
					//~ PrepareTranspose(k,k+1,L_Address,L_Entries,TransposeU);
					if (L_Address[k].Size() > 0)
					{
						INMOST_DATA_ENUM_TYPE i = L_Entries[L_Address[k].thr][L_Address[k].first].first;
						assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
						Llist[k] = Lbeg[i];
						Lbeg[i] = k;
					}
					// Shift indexes for transposed L-part traversal for next iteration
					{
						INMOST_DATA_ENUM_TYPE  i = Lbeg[k];
						while (i != EOL)
						{
							L_Address[i].first++;
							INMOST_DATA_ENUM_TYPE Li = Llist[i], Ui;
							if (L_Address[i].Size() > 0)
							{
								Ui = L_Entries[L_Address[i].thr][L_Address[i].first].first;
								Llist[i] = Lbeg[Ui];
								Lbeg[Ui] = i;
							}
							i = Li;
						}
					}
#if defined(ILUC2)
					// Insert indexes for transposed second-order L-part traversal
					//~ PrepareTranspose(k,k+1,L2_Address,L2_Entries,TransposeL2);
					if (L2_Address[k].Size() != 0)
					{
						INMOST_DATA_ENUM_TYPE i = L2_Entries[L2_Address[k].thr][L2_Address[k].first].first;
						assert(static_cast<INMOST_DATA_INTEGER_TYPE>(i) > k);
						L2list[k] = L2beg[i];
						L2beg[i] = k;
					}
					// Shift indexes for transposed second-order L-part traversal
					{
						INMOST_DATA_ENUM_TYPE i = L2beg[k];
						while (i != EOL)
						{
							L2_Address[i].first++;
							INMOST_DATA_ENUM_TYPE Li = L2list[i], Ui;
							if (L2_Address[i].Size())
							{
								Ui = L2_Entries[L2_Address[i].thr][L2_Address[i].first].first;
								L2list[i] = L2beg[Ui];
								L2beg[Ui] = i;
							}
							i = Li;
						}
					}
#endif
				}
			}
			if (verbosity)
			{
				//number of nonzeros
				nzLU += U_Address[k].Size() + L_Address[k].Size() + 1;
#if defined(ILUC2)
				nzLU2 += U2_Address[k].Size() + L2_Address[k].Size();
#endif
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
				progress_cur++;
				if (k % report_pace == 0)
				{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
					{
						std::ios save(NULL);
						save.copyfmt(std::cout);
						if (verbosity == 1)
							//std::cout << cbeg << "<->" << cend << " factor " << std::setw(6) << std::fixed << std::setprecision(2) << 100.0f * (k - cbeg) / (float)(cend - cbeg) << "\r" << std::flush;
							std::cout << cbeg << "<->" << cend << " factor " << std::setw(6) << std::fixed << std::setprecision(2) << 100.0f * progress_cur / (float)progress_all << "\r" << std::flush;
						else
						{
							std::cout << std::setw(6) << std::fixed << std::setprecision(2) << 100.0f * progress_cur / (float)progress_all;
							std::cout << " nnz LU " << std::setw(8) << nzLU;
							std::cout << " LU2 " << std::setw(8) << nzLU2;
							std::cout << " cond L " << std::setprecision(6) << std::setw(10) << NuL;
							std::cout << " D " << std::setw(10) << NuD;
							std::cout << " U " << std::setw(10) << NuU;
							std::cout << " swaps " << std::setw(10) << swaps;
							std::cout << " drops " << std::setw(10) << ndrops;
							std::cout << "\r" << std::flush;
						}
						std::cout.copyfmt(save);
					}
				}
			}
			// iteration done
		}
		//restore indexes
		for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
		{
			A_Address[k].first = Apos[k];
			U_Address[k].first = Upos[k];
			L_Address[k].first = Lpos[k];
#if defined(ILUC2)
			U2_Address[k].first = U2pos[k];
			L2_Address[k].first = L2pos[k];
#endif
		}
		if (verbosity > 1)
		{
			INMOST_DATA_ENUM_TYPE nzA = 0;
			for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
			{
				for (INMOST_DATA_ENUM_TYPE jt = A_Address[k].first; jt != A_Address[k].last; ++jt)
				{
					INMOST_DATA_ENUM_TYPE j = A_Entries[A_Address[k].thr][jt].first;
					if (j >= cbeg && j < cend) nzA++;
				}
			}
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
			{
				std::cout << "size " << cend - cbeg;
				std::cout << " total nonzeros in A " << nzA << " (sparsity " << nzA / (double)(cend - cbeg) / (double)(cend - cbeg) * 100.0 << "%) in LU " << nzLU << " (fill " << nzLU / (double)nzA * 100.0 << "%)";
				std::cout << " in LU2 " << nzLU2;
				std::cout << " conditions L " << NuL_max << " D " << NuD << " U " << NuU_max << " pivot swaps " << swaps << " drops " << ndrops << "            " << std::endl;
			}
		}
		NuLout = std::max(NuLout,NuL);
		NuUout = std::max(NuUout,NuU);
	}


	
	void MLMTILUC_preconditioner::MetisOrdering(INMOST_DATA_ENUM_TYPE wbeg,
												INMOST_DATA_ENUM_TYPE wend,
												std::vector<INMOST_DATA_ENUM_TYPE> & xadj, 
												std::vector<INMOST_DATA_ENUM_TYPE> & adjncy,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localP,
												interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> & localQ)
	{
#if defined(USE_SOLVER_METIS)
		if( !adjncy.empty() )
		{
			idx_t nvtxs = wend-wbeg;
			idx_t options[METIS_NOPTIONS];
			METIS_SetDefaultOptions(options);
			options[METIS_OPTION_NUMBERING] = 0;
			std::vector<idx_t> ixadj(xadj.size()), iadjncy(adjncy.size()), perm(nvtxs), iperm(nvtxs);
			std::copy(xadj.begin(),xadj.end(),ixadj.begin());
			std::copy(adjncy.begin(),adjncy.end(),iadjncy.begin());
			
			//options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_TIME | METIS_DBG_COARSEN | METIS_DBG_CONNINFO | METIS_DBG_CONTIGINFO | METIS_DBG_IPART | METIS_DBG_MEMORY | METIS_DBG_MOVEINFO | METIS_DBG_REFINE | METIS_DBG_SEPINFO;
			//options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
			//options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED;
			//options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
			//options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
			//options[METIS_OPTION_NO2HOP] = 0;
			//options[METIS_OPTION_MINCONN] = 0;
			//options[METIS_OPTION_NITER] = 4;
			//METIS_NodeNDP(nvtxs,&xadj[0],&adjncy[0],NULL,4,options,&perm[0],&iperm[0],&sizes[0]);
			METIS_NodeND(&nvtxs,&ixadj[0],&iadjncy[0],NULL,options,&perm[0],&iperm[0]);
			
			
			for(INMOST_DATA_ENUM_TYPE  k = wbeg; k < wend; ++k)
			{
				//localP[k] = iperm[k-wbeg]+wbeg;
				//localQ[k] = iperm[k-wbeg]+wbeg;
				localP[perm[k-wbeg]+wbeg] = k;
				localQ[perm[k-wbeg]+wbeg] = k;
			}
		}
		else
#endif
		{
			for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
			{
				localP[k] = k;
				localQ[k] = k;
			}
		}
	}

	struct sort_pred 
	{
		bool operator()(const Sparse::Row::entry& left, const Sparse::Row::entry& right) 
		{
			return left.second > right.second;
		}
	};
	
	bool MLMTILUC_preconditioner::Initialize()
	{
		if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		bool block_pivot = false;
		
		INMOST_DATA_ENUM_TYPE wbeg, wend; //working interval
		INMOST_DATA_ENUM_TYPE mobeg, moend; // total interval
		INMOST_DATA_ENUM_TYPE vbeg, vend; // vector interval
		
		//Sparse::Vector DL, DR;
		Sparse::Vector DL0,DR0;
		info->GetOverlapRegion(info->GetRank(), mobeg, moend);
		info->GetVectorRegion(vbeg,vend);
		
		

		
		//prepare reordering vectors
		ddP.set_interval_beg(mobeg);
		ddP.set_interval_end(moend);
		ddQ.set_interval_beg(mobeg);
		ddQ.set_interval_end(moend);

		

		//prepare rescaling vectors
		//DL.SetInterval(mobeg, moend);
		//DR.SetInterval(mobeg, moend);
		DL0.SetInterval(mobeg, moend);
		DR0.SetInterval(mobeg, moend);
		for(INMOST_DATA_ENUM_TYPE k = mobeg; k < moend; ++k) 
			DL0[k] = DR0[k] = 1.0;
		
		

		for (INMOST_DATA_ENUM_TYPE k = mobeg; k < moend; k++)
			ddP[k] = ddQ[k] = k;
		
		if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

		//supplimentary data structures for ILUC
		U_Address.set_interval_beg(mobeg);
		U_Address.set_interval_end(moend);
		U_Entries.resize(Threads());
		L_Address.set_interval_beg(mobeg);
		L_Address.set_interval_end(moend);
		L_Entries.resize(Threads());
		LU_Diag.set_interval_beg(mobeg);
		LU_Diag.set_interval_end(moend);
		std::vector< std::vector<Sparse::Row::entry> > A_Entries(Threads());
		E_Entries.resize(Threads());
		F_Entries.resize(Threads());
		
		//std::vector<INMOST_DATA_ENUM_TYPE> sort_indeces;
		interval<INMOST_DATA_ENUM_TYPE, Interval> A_Address(mobeg, moend);
		
		
		INMOST_DATA_ENUM_TYPE nnzE = 0, nnzF = 0, nnzA = 0;
		
		//supplimentary data structures for returning values of dropped elements
		//~ INMOST_DATA_REAL_TYPE DropLk, DropUk, SumLk, SumUk;
		//~ interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU(mobeg,moend,0.0), DropL(mobeg,moend,0.0);
//~ #if defined(ILUC2)
		//~ interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DropU2(mobeg,moend,0.0), DropL2(mobeg,moend,0.0);
//~ #endif
		//data structure for linked list

		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(mobeg,moend,std::numeric_limits<INMOST_DATA_REAL_TYPE>::max());
		//interval<INMOST_DATA_ENUM_TYPE, std::vector< std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_ENUM_TYPE> > > TransposeU(mobeg,moend), TransposeL(mobeg,moend);
		
		
		
        double tfactor = 0.0, trescale = 0.0, treorder = 0.0, ttransversal = 0.0, treassamble = 0.0, ttotal, testimator = 0.0, tschur = 0.0, tnd = 0, torder = 0, tlocal;
	    double tlschur;
		INMOST_DATA_ENUM_TYPE totswaps = 0;
		ttotal = Timer();

		//(*Alink).Save("M.mtx");
		
		//calculate number of nonzeros
		INMOST_DATA_ENUM_TYPE nzA = 0;
		for (INMOST_DATA_ENUM_TYPE k = mobeg; k < moend; ++k)
		{
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
				if (r->first >= mobeg && r->first < moend && !check_zero(fabs(r->second)) ) nzA++;
		}
		INMOST_DATA_ENUM_TYPE nzA0 = nzA;

		//sort_indeces.reserve(256);
		//A_Entries.resize(nzA);
		//LU_Entries.reserve(nzA*4);

#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
		for (INMOST_DATA_INTEGER_TYPE k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k)
		{
			int thr = Thread();
			A_Address[k].thr = thr;
			A_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(A_Entries[thr].size());
			for (Sparse::Row::iterator r = (*Alink)[k].Begin(); r != (*Alink)[k].End(); ++r)
			{
				if (r->first >= mobeg && r->first < moend && (!check_zero(r->second)) )
					A_Entries[thr].push_back(Sparse::Row::make_entry(r->first, r->second));
			}
			A_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(A_Entries[thr].size());
			assert(A_Address[k].Size() != 0); //singular matrix
		}
		

		if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
		
		
		INMOST_DATA_ENUM_TYPE nzLU = 0, nzLU2 = 0, nzLU2tot = 0;
#if defined(ILUC2)
		INMOST_DATA_REAL_TYPE tau2 = iluc2_tau;
#else
		INMOST_DATA_REAL_TYPE tau2 = tau;
#endif


		for(INMOST_DATA_ENUM_TYPE k = mobeg; k < moend; ++k)
		{
			U_Address[k].first = U_Address[k].last = ENUMUNDEF;
			L_Address[k].first = L_Address[k].last = ENUMUNDEF;
		}

		if( verbosity > 1 ) 
			std::cout << "nonzeros in A " << nzA << " pivot cond " << pivot_cond << " diag " << pivot_diag << " tau " << tau << " tau2 " << tau2 << std::endl;

		

		//set up working interval
		wbeg = mobeg;
		wend = moend;
		INMOST_DATA_REAL_TYPE NuL = 1, NuU = 1;
		while (wbeg < wend) //iterate into levels until all matrix is factored
		{
			interval<INMOST_DATA_ENUM_TYPE, bool> Pivot(wbeg, wend, false);
			//reordering on current Schur complement
			std::vector<Block> blocks; //blocks from nested dissection
			INMOST_DATA_ENUM_TYPE cbeg = wbeg, cend = wend; //next size of factored B block
			

			if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
			//this scope is for reordering and rescaling
			{
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> localP(wbeg, wend, ENUMUNDEF), localQ(wbeg, wend, ENUMUNDEF);
				interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> DL(wbeg,wend,1.0), DR(wbeg,wend,1.0);
				//this scope is for reorderings
				{
					if (verbosity > 1) std::cout << "Reorder\n";
					double tlreorder = Timer();
					if (run_nd == 1 && !block_pivot && Threads() > 1 && wend-wbeg > 128)
					{
						if (verbosity > 1) std::cout << "Reordering with k-way dissection, threads = " << Threads() << std::endl;
						tlocal = Timer();
						//NestedDissection(wbeg, wend, A_Address, A_Entries, localP, localQ, blocks, (wend - wbeg) / Threads());
						//KwayDissection(wbeg,wend,A_Address,A_Entries,localP,localQ,blocks,Threads());
						//KwayDissection(wbeg, wend, A_Address, A_Entries, localP, localQ, blocks, Threads());
						KwayDissection(wbeg, wend, A_Address, A_Entries, localP, localQ, blocks, Threads());
						//complete matrix ordering
						tlocal = Timer() - tlocal;
						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";
						tnd += tlocal;
						treorder += tlocal;


						size_t sep = 0;
						for (size_t q = 0; q < blocks.size(); ++q)
							if (blocks[q].separator) sep += blocks[q].row_end - blocks[q].row_start;
						
						if (2*sep < (wend - wbeg) ) //separator is not too big
						{
							if (verbosity > 1) std::cout << "Reassemble\n";
							tlocal = Timer();
							ReorderSystem(wbeg, wend, A_Address, A_Entries, localP, localQ, DL, DR, DL0, DR0);
							tlocal = Timer() - tlocal;
							treassamble += tlocal;
							if (verbosity > 1) std::cout << "Time " << tlocal << "\n";

							if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

							if (false)
							{
								DumpMatrix(A_Address, A_Entries, wbeg, wend, "A_nd" + to_string(level_size.size()) + ".mtx");
								std::ofstream file(("blocks_nd" + to_string(level_size.size()) + ".txt").c_str());
								for (INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
									file << blocks[k].row_start - wbeg << " " << blocks[k].row_end - wbeg << " " << blocks[k].col_start - wbeg << " " << blocks[k].col_end - wbeg << std::endl;
								file.close();
							}
						}
						else
						{
							if (verbosity > 1) std::cout << "Separator is too big " << sep << " matrix size " << wend-wbeg << "\n";
							blocks.clear();
							blocks.push_back(Block(wbeg, wend, wbeg, wend, false));
							for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
							{
								localP[k] = k;
								localQ[k] = k;
							}
						}
						//exit(-1);
					}
					else blocks.push_back(Block(wbeg, wend, wbeg, wend,false));

					//DumpMatrix(A_Address,A_Entries,cbeg,cend,"mat"+to_string(level_size.size())+".mtx");
					/// MAXIMUM TRANSVERSE REORDERING
					if( run_mpt )
					{
						if( verbosity > 1 ) std::cout << "Reordering with MPT\n";
						tlocal = Timer();
						for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
							localQ[k] = ENUMUNDEF;

						for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
							CheckColumnGaps(blocks[q], A_Address, A_Entries);

						//for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
						//	CheckBlock(blocks[q], A_Address, A_Entries, wend, wend, __FILE__, __LINE__); //no separator

#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
						{
							INMOST_DATA_ENUM_TYPE tseg = (wend - wbeg) / Threads();
							INMOST_DATA_ENUM_TYPE tbeg = wbeg + tseg * Thread();
							INMOST_DATA_ENUM_TYPE tend = tbeg + tseg;

							for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
								std::sort(A_Entries[A_Address[k].thr].begin() + A_Address[k].first, A_Entries[A_Address[k].thr].begin() + A_Address[k].last, sort_pred());
						}

#if defined(USE_OMP_FACT)
#pragma omp parallel for schedule(static,1)
#endif
						for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
						{
							if (verbosity > 1)
							{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
								{
									std::cout << "Block " << q;
									std::cout << " rows " << blocks[q].row_start << ":" << blocks[q].row_end;
									std::cout << " cols " << blocks[q].col_start << ":" << blocks[q].col_end;
									std::cout << " nnz " << ComputeNonzeroes(blocks[q], A_Address, A_Entries);
									std::cout << std::endl;
								}
							}

							//if (blocks.size() > 1)
							{
								//sort so that largest values are encountered first
								//maximum transversal algorithm may choose less beneficial path if not sorted
								//for (INMOST_DATA_ENUM_TYPE k = blocks[q].row_start; k < blocks[q].row_end; ++k)
								//	std::sort(A_Entries[A_Address[k].thr].begin() + A_Address[k].first, A_Entries[A_Address[k].thr].begin() + A_Address[k].last, sort_pred());
							}
							//DumpMatrixBlock(blocks[q], A_Address, A_Entries, "block_" + to_string(level_size.size()) + "_" + to_string(q) + ".mtx");

							
							MaximalTransversal(blocks[q].row_start, blocks[q].row_end, blocks[q].col_start, blocks[q].col_end, A_Address, A_Entries, localQ, DL, DR);

							if (verbosity > 1)
							{
								INMOST_DATA_ENUM_TYPE cmin = moend, cmax = mobeg;
								for (INMOST_DATA_ENUM_TYPE k = blocks[q].col_start; k < blocks[q].col_end; ++k) if (localQ[k] != ENUMUNDEF)
								{
									cmin = std::min(cmin, localQ[k]);
									cmax = std::max(cmax, localQ[k]);
								}
								Block b(blocks[q].row_start, blocks[q].row_end, blocks[q].row_start, blocks[q].row_end, false);
								INMOST_DATA_ENUM_TYPE nnz = ComputeNonzeroes(b, A_Address, A_Entries, localQ);
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
								std::cout << "Block " << q << " column indices in [" << cmin << ":" << cmax << "] rows " << blocks[q].row_start << ":" << blocks[q].row_end << " columns " << blocks[q].col_start << ":" << blocks[q].col_end << " nnz " << nnz << std::endl;
							}
							
						}
						tlocal = Timer() - tlocal;
						if( verbosity > 1 ) std::cout << "Time " << tlocal << "\n";
						ttransversal += tlocal;
						treorder += tlocal;

						
						
						{ //collect columns and shift gaps in localQ
							//collect gaps
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
								localP[k] = ENUMUNDEF;
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
							{
								if (localQ[k] != ENUMUNDEF)
								{
									if (localP[localQ[k]] != ENUMUNDEF)
									{
#if defined(USE_OMP)
#pragma omp critical
#endif
										std::cout << "k " << k << " Q[k] " << localQ[k] << " P[Q[k]] " << localP[localQ[k]] << std::endl;
									}
									localP[localQ[k]] = k;
								}
							}
							//debug checks
#ifndef NDEBUG
							for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
							{
								INMOST_DATA_ENUM_TYPE rbeg = blocks[q].row_start;
								INMOST_DATA_ENUM_TYPE rend = blocks[q].row_end;
								INMOST_DATA_ENUM_TYPE cbeg = blocks[q].col_start;
								INMOST_DATA_ENUM_TYPE cend = blocks[q].col_end;
								INMOST_DATA_ENUM_TYPE first = moend;// , firstc = moend;
								INMOST_DATA_ENUM_TYPE last = mobeg;// , lastc = mobeg;
								INMOST_DATA_ENUM_TYPE skip = 0;
								for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
								{
									if (localQ[k] != ENUMUNDEF)
									{
										first = std::min(first, localQ[k]);
										last = std::max(last, localQ[k]);
									}
									else skip++;
								}
								if (first != rbeg) std::cout << " !!!! first column index " << first << " row index " << rbeg << " block " << q << std::endl;
								if (last + 1 != rend) std::cout << " !!!! last column index " << last + 1 << " row index " << rend << " block " << q << std::endl;
								if (blocks[q].row_end - blocks[q].row_start > blocks[q].col_end - blocks[q].col_start - skip)
									std::cout << " @@@@ block " << q << " number of rows: " << blocks[q].row_end - blocks[q].row_start << " assigned columns: " << blocks[q].col_end - blocks[q].col_start - skip << " total columns " << blocks[q].col_end - blocks[q].col_start << " skip " << skip << std::endl;
							}
#endif
							std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > gaps(blocks.size());
#if defined(USE_OMP)
#pragma omp parallel for schedule(static,1)
#endif
							for (int q = 0; q < (int)blocks.size(); ++q)
							{
								for (INMOST_DATA_ENUM_TYPE k = blocks[q].row_end; k > blocks[q].row_start; --k)
									if (localP[k - 1] == ENUMUNDEF)
									{

										if (verbosity > 1 && !blocks[q].separator)
										{
#if defined(USE_OMP)
#pragma omp critical
#endif
											{
												std::cout << "Gap at " << k - 1 << " is from block " << q;
												std::cout << " rows " << blocks[q].row_start << ":" << blocks[q].row_end;
												std::cout << " cols " << blocks[q].col_start << ":" << blocks[q].col_end;
												std::cout << std::endl;
											}
										}
										gaps[q].push_back(k - 1);
									}
								if (verbosity > 1 && !gaps[q].empty() && !blocks[q].separator)
								{
#if defined(USE_OMP)
#pragma omp critical
#endif
									std::cout << "block " << q << " " << (blocks[q].separator?"separator":"internal") << " gaps " << gaps[q].size() << std::endl;
								}
							}
							//fill gaps for internal blocks, some columns fill the gaps
#if defined(USE_OMP)
#pragma omp parallel for schedule(static,1)
#endif
							for (int q = 0; q < (int)blocks.size(); ++q) if( !blocks[q].separator )
							{
								if (!gaps[q].empty())
								{
									if (verbosity > 1)
									{
#if defined(USE_OMP)
#pragma omp critical
#endif
										std::cout << "assigning columns for block " << q << " gaps " << gaps[q].size() << " columns " << blocks[q].col_start << ":" << blocks[q].col_end << std::endl;
									}
									for (INMOST_DATA_ENUM_TYPE k = blocks[q].col_start; k < blocks[q].col_end; ++k)
									{
										if (localQ[k] == ENUMUNDEF)
										{
											if (verbosity > 1)
											{
#if defined(USE_OMP)
#pragma omp critical
#endif
												std::cout << "assign column " << k << " to " << gaps[q].back() << " for block " << q << std::endl;
											}
											localQ[k] = gaps[q].back();
											gaps[q].pop_back();
											if (gaps[q].empty()) 
												break;
										}
									}
								}
								if (verbosity > 1 && !gaps[q].empty())
								{
#if defined(USE_OMP)
#pragma omp critical
#endif
									std::cout << "block " << q << " still have gaps " << gaps[q].size() << std::endl;
								}
							}
							//fill gaps and set pivots for separators
							for (int q = 0; q < (int)blocks.size(); ++q) if( blocks[q].separator )
							{
								INMOST_DATA_ENUM_TYPE rbeg = blocks[q].row_start;
								INMOST_DATA_ENUM_TYPE rend = blocks[q].row_end;
								INMOST_DATA_ENUM_TYPE cbeg = blocks[q].col_start;
								INMOST_DATA_ENUM_TYPE cend = blocks[q].col_end;
								for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k)
									Pivot[k] = true;
								INMOST_DATA_ENUM_TYPE gbeg = wend, gend = wbeg, gtot = 0;
								for (INMOST_DATA_ENUM_TYPE k = cbeg; k < cend; ++k)
									if (localQ[k] == ENUMUNDEF)
									{
										if (!gaps[q].empty())
										{
											localQ[k] = gaps[q].back();
											gaps[q].pop_back();
											gtot++;
											gbeg = std::min(gbeg, localQ[k]);
											gend = std::max(gend, localQ[k]);
										}
										//else localQ[k] = last++;
										else std::cout << __FILE__ << ":" << __LINE__ << " gaps for block " << q << " are empty! k = " << k << " interval " << cbeg << ":" << cend << std::endl;
									}
								if (verbosity > 1)
								{
									{
										std::cout << "Separator " << q;
										std::cout << " rows " << blocks[q].row_start << ":" << blocks[q].row_end;
										std::cout << " cols " << blocks[q].col_start << ":" << blocks[q].col_end;
										std::cout << std::endl;
										std::cout << "gaps at " << gbeg << ":" << gend << " total " << gtot << std::endl;
									}
								}
							}
							if (verbosity > 1)
							{
								for (size_t q = 0; q < gaps.size(); ++q)
									if (!gaps[q].empty()) std::cout << "block " << q << " gaps " << gaps[q].size() << std::endl;
								for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
									if (localQ[k] == ENUMUNDEF)
									{
										std::cout << "Column " << k << " not assigned" << std::endl;
										for (int q = 0; q < (int)blocks.size(); ++q)
											if (blocks[q].col_start <= k && k < blocks[q].col_end)
												std::cout << "in block " << q << " " << (blocks[q].separator ? "separator" : "internal") << std::endl;
									}
							}
							
							for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
							{
								blocks[q].col_start = blocks[q].row_start;
								blocks[q].col_end = blocks[q].row_end;
							}
							for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
								localP[k] = k;

							//check blocks
							{
								INMOST_DATA_ENUM_TYPE sepbeg = wend, sepend = wend;
								for (size_t q = 0; q < blocks.size(); ++q) if (blocks[q].separator)
								{
									sepbeg = blocks[q].row_start;
									sepend = blocks[q].row_end;
								}
								//for (size_t q = 0; q < blocks.size(); ++q) if (!blocks[q].separator)
								//	CheckBlock(blocks[q], A_Address, A_Entries, localQ, sepbeg, sepend, __FILE__, __LINE__); //no separator
							}
						}

						if (verbosity > 1) std::cout << "Reassemble columns\n";
						tlocal = Timer();
						ReorderColumns(wbeg, wend, A_Address, A_Entries, localQ, DR, DR0);
						tlocal = Timer() - tlocal;
						treassamble += tlocal;

						//check blocks
						{
							INMOST_DATA_ENUM_TYPE sepbeg = wend, sepend = wend;
							for (size_t q = 0; q < blocks.size(); ++q) if (blocks[q].separator)
							{
								sepbeg = blocks[q].row_start;
								sepend = blocks[q].row_end;
							}
							//for (size_t q = 0; q < blocks.size(); ++q) if (!blocks[q].separator)
							//	CheckBlock(blocks[q], A_Address, A_Entries, sepbeg, sepend, __FILE__, __LINE__); //no separator
						}

						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";

						if( false )
						{
							DumpMatrix(A_Address, A_Entries, wbeg, wend, "A_mt"+to_string(level_size.size())+".mtx");
							std::ofstream file(("blocks_mt" + to_string(level_size.size()) + ".txt").c_str());
							for (INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
								file << blocks[k].row_start-wbeg << " " << blocks[k].row_end-wbeg << " " << blocks[k].col_start-wbeg << " " << blocks[k].col_end-wbeg << std::endl;
							file.close();
						}
						//exit(-1);

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
					}
					else
					{
						for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
						{
							localQ[k] = k;
							DL[k] = DR[k] = 1;
						}
					}
					///  END MAXIMUM TRANSVERSE REORDERING
					cend = wend;

					if (run_nd == 2 && !block_pivot && Threads() > 1 && wend - wbeg > 128)
					{
						if (verbosity > 1) std::cout << "Reordering with k-way symmetric dissection, threads = " << Threads() << std::endl;
						tlocal = Timer();
						blocks.clear();
						KwaySymmetricDissection(wbeg, wend, A_Address, A_Entries, localP, localQ, blocks, Threads());
						//complete matrix ordering
						tlocal = Timer() - tlocal;
						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";
						tnd += tlocal;
						treorder += tlocal;


						size_t sep = 0;
						for (size_t q = 0; q < blocks.size(); ++q)
							if (blocks[q].separator) sep += blocks[q].row_end - blocks[q].row_start;

						if (2 * sep < (wend - wbeg)) //separator is not too big
						{
							if (verbosity > 1) std::cout << "Reassemble\n";
							tlocal = Timer();
							ReorderSystem(wbeg, wend, A_Address, A_Entries, localP, localQ, DL, DR, DL0, DR0);
							tlocal = Timer() - tlocal;
							treassamble += tlocal;
							if (verbosity > 1) std::cout << "Time " << tlocal << "\n";

							if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

							if (false)
							{
								for (size_t q = 0; q < blocks.size(); ++q) if (blocks[q].separator)
								{
									blocks.push_back(Block(blocks[q].col_start, blocks[q].col_end, blocks[q].row_start, blocks[q].row_end, true));
									break;
								}

								DumpMatrix(A_Address, A_Entries, wbeg, wend, "A_snd" + to_string(level_size.size()) + ".mtx");
								std::ofstream file(("blocks_snd" + to_string(level_size.size()) + ".txt").c_str());
								for (INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
									file << blocks[k].row_start - wbeg << " " << blocks[k].row_end - wbeg << " " << blocks[k].col_start - wbeg << " " << blocks[k].col_end - wbeg << std::endl;
								file.close();

								if( blocks.back().separator ) blocks.pop_back();
							}

							if (verbosity > 1)
							{
								for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
								{
									std::cout << "Block " << q;
									std::cout << " rows " << blocks[q].row_start << ":" << blocks[q].row_end;
									std::cout << " cols " << blocks[q].col_start << ":" << blocks[q].col_end;
									std::cout << " nnz " << ComputeNonzeroes(blocks[q], A_Address, A_Entries);
									std::cout << std::endl;
								}
							}

							for (int q = 0; q < (int)blocks.size(); ++q) if (blocks[q].separator)
								for (INMOST_DATA_ENUM_TYPE k = blocks[q].row_start; k < blocks[q].row_end; ++k)
									Pivot[k] = true;
						}
						else
						{
							if (verbosity > 1) std::cout << "Separator is too big " << sep << " matrix size " << wend - wbeg << "\n";
							blocks.clear();
							blocks.push_back(Block(wbeg, wend, wbeg, wend, false));
							for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
							{
								localP[k] = k;
								localQ[k] = k;
							}
						}
						//exit(-1);
					}
					else if (blocks.empty()) blocks.push_back(Block(wbeg, wend, wbeg, wend, false));


					if (reorder_b)
					{
						tlocal = Timer();

						if (verbosity > 1) std::cout << "Reordering\n";


						for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
						{
							localP[k] = k;
							localQ[k] = k;
						}

						//check blocks
						{
							INMOST_DATA_ENUM_TYPE sepbeg = wend, sepend = wend;
							for (size_t q = 0; q < blocks.size(); ++q) if (blocks[q].separator)
							{
								sepbeg = blocks[q].row_start;
								sepend = blocks[q].row_end;
							}
							//for (size_t q = 0; q < blocks.size(); ++q) if (!blocks[q].separator)
							//	CheckBlock(blocks[q], A_Address, A_Entries, sepbeg, sepend, __FILE__, __LINE__); //no separator
						}
#if defined(USE_OMP_FACT)
#pragma omp parallel for schedule(static,1)
#endif
						for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
						{
							INMOST_DATA_ENUM_TYPE rbeg = blocks[q].row_start;
							INMOST_DATA_ENUM_TYPE rend = blocks[q].row_end;

							std::vector<INMOST_DATA_ENUM_TYPE> xadj(rend - rbeg + 1), adjncy;



							SymmetricGraph(rbeg, rend, A_Address, A_Entries, Pivot, xadj, adjncy);
							
							
#if defined(REORDER_METIS_ND)
							{
								idx_t nvtxs = rend - rbeg;
								MetisOrdering(rbeg, rend, xadj, adjncy, localP, localQ);
							}
#elif defined(REORDER_RCM)
							RCMOrdering(rbeg, rend, xadj, adjncy, localP, localQ);
#elif defined(REORDER_WRCM)
							{
								std::vector<INMOST_DATA_REAL_TYPE> wadj(rend - rbeg);
								GraphWeights(rbeg, rend, A_Address, A_Entries, Pivot, DL, DR, wadj);
								WRCMOrdering(rbeg, rend, xadj, adjncy, wadj, localP, localQ);
							}
#endif
						}
						tlocal = Timer() - tlocal;
						torder += tlocal;
						treorder += tlocal;
						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";
						/// REASSAMBLE
						if (verbosity > 1) std::cout << "Reassemble\n";
						tlocal = Timer();
						ReorderSystem(wbeg, wend, A_Address, A_Entries, localP, localQ, DL, DR, DL0, DR0);
						tlocal = Timer() - tlocal;
						treassamble += tlocal;
						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";

						//check blocks
						{
							INMOST_DATA_ENUM_TYPE sepbeg = wend, sepend = wend;
							for (size_t q = 0; q < blocks.size(); ++q) if (blocks[q].separator)
							{
								sepbeg = blocks[q].row_start;
								sepend = blocks[q].row_end;
							}
							//for (size_t q = 0; q < blocks.size(); ++q) if (!blocks[q].separator)
							//	CheckBlock(blocks[q], A_Address, A_Entries, sepbeg, sepend, __FILE__, __LINE__); //no separator
						}

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;


						if (false)
						{
							for (size_t q = 0; q < blocks.size(); ++q) if (blocks[q].separator)
							{
								blocks.push_back(Block(blocks[q].col_start, blocks[q].col_end, blocks[q].row_start, blocks[q].row_end, true));
								break;
							}
							DumpMatrix(A_Address, A_Entries, wbeg, wend, "A_ord" + to_string(level_size.size()) + ".mtx");
							std::ofstream file(("blocks_ord" + to_string(level_size.size()) + ".txt").c_str());
							for (INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
								file << blocks[k].row_start-wbeg << " " << blocks[k].row_end-wbeg << " " << blocks[k].col_start-wbeg << " " << blocks[k].col_end-wbeg << std::endl;
							file.close();

							if (blocks.back().separator) blocks.pop_back();
						}
						//exit(-1);
					}
					else
					{
						for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
						{
							localP[k] = k;
							localQ[k] = k;
						}
					}
					
					tlreorder = Timer() - tlreorder;
					if( verbosity > 1 ) std::cout << "Reordering time " << tlreorder << "\n";
					
					
					


					{ //report A data
						nzA = 0;
						for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
							nzA += A_Address[k].Size();
						double sparsity = nzA / (double)(wend - wbeg) / (double)(wend - wbeg);
						if (verbosity > 1)
							std::cout << "nzA " << nzA << " size " << wend - wbeg << " sparsity " << sparsity << "\n";

						if ((sparsity > 0.85) || (wend - wbeg < 32))
							block_pivot = true;
					}

					if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
				} //scope for reordering
				/// RESCALING
				{ //scope for rescalings
					//Rescale current B block by Row-Column alternating scaling
					if( rescale_b )
					{
						if( verbosity > 1 ) std::cout << "Rescaling, iters " << sciters << std::endl;
						tlocal = Timer();
#if defined(USE_OMP_FACT)
#pragma omp parallel for schedule(static,1)
#endif
						for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
						{
							INMOST_DATA_ENUM_TYPE rbeg = blocks[q].row_start;
							INMOST_DATA_ENUM_TYPE rend = blocks[q].row_end;

#if defined(EQUALIZE_1NORM)
							SinkhornScaling(rbeg, rend, A_Address, A_Entries, Pivot, DL,  DR0, 1);
#elif defined(EQUALIZE_2NORM)
							SinkhornScaling(rbeg, rend, A_Address, A_Entries, Pivot, DL, DR,  2);
#elif defined(EQUALIZE_IDOMINANCE)
							IdominanceScaling(rbeg, rend, A_Address, A_Entries, Pivot, DL, DR);
#endif

							//DumpMatrixBlock(blocks[q], A_Address, A_Entries, "sc_block_" + to_string(level_size.size()) + "_" + to_string(q) + ".mtx");
						}
						//stack scaling
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
						for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); k++)
						{
							DL0[k] *= DL[k];
							DR0[k] *= DR[k];
						}
						tlocal = Timer() - tlocal;
						trescale += tlocal;
						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";
						if (verbosity > 1) std::cout << "Pivot rescaling" << std::endl;
						tlocal = Timer();
						PivotScaling(wbeg, wend, A_Address, A_Entries, Pivot, DL, DR);
						tlocal = Timer() - tlocal;
						trescale += tlocal;
						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";
						if (verbosity > 1) std::cout << "EF rescaling" << std::endl;
						tlocal = Timer();
						EFScaling(wbeg, wend, mobeg, moend, DL, DR);
						tlocal = Timer() - tlocal;
						trescale += tlocal;


						if (verbosity > 1) std::cout << "Time " << tlocal << "\n";
						
					}
				} //scope for scalings


				if (false)
				{
					DumpMatrix(A_Address, A_Entries, wbeg, wend, "scA.mtx");
					std::ofstream file("blocks_sc.txt");
					for (INMOST_DATA_ENUM_TYPE k = 0; k < blocks.size(); ++k)
						file << blocks[k].row_start << " " << blocks[k].row_end << " " << blocks[k].col_start << " " << blocks[k].col_end << std::endl;
					file.close();
				}
				//exit(-1);

			} //scope for scalings and reorderings

			if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
			/*
			{
				std::stringstream s;
				s << "B" << level_size.size() << ".mtx";
				DumpMatrix(A_Address,A_Entries,cbeg,cend,s.str());
				std::cout << "Press any key" << std::endl;
				scanf("%*c");
			}
			*/
			{ //this scope is for factorization and schur
				
				tlocal = Timer();
				std::vector< std::vector<Sparse::Row::entry> > L2_Entries(Threads()), U2_Entries(Threads());
				interval<INMOST_DATA_ENUM_TYPE, Interval> L2_Address(wbeg,wend), U2_Address(wbeg,wend);
				//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2list(cbeg, cend, UNDEF), L2list(cbeg, cend, UNDEF);
				//interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> U2beg(cbeg, cend, EOL), L2beg(cbeg, cend, EOL);
				progress_all = wend - wbeg;
				progress_cur = 0;
#if defined(USE_OMP_FACT)
				//int nested = omp_get_nested();
				//omp_set_nested(1);
#pragma omp parallel for schedule(static,1)
#endif
				for (int q = 0; q < (int)blocks.size(); ++q) if (!blocks[q].separator)
				{/// FACTORIZATION BEGIN
					INMOST_DATA_ENUM_TYPE rbeg = blocks[q].row_start;
					INMOST_DATA_ENUM_TYPE rend = blocks[q].row_end;
					Factorize(rbeg, rend, A_Address, A_Entries, L2_Address, L2_Entries, U2_Address, U2_Entries, Pivot,block_pivot,NuL,NuU,testimator);
				} //end of factorization scope
#if defined(USE_OMP_FACT)
				//omp_set_nested(nested);
#endif
				nzLU = 0;
				for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
					nzLU += L_Address[k].Size() + U_Address[k].Size();
#if defined(ILUC2)
				nzLU2 = 0;
				for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
					nzLU2 += L2_Address[k].Size() + U2_Address[k].Size();
				nzLU2tot += nzLU2;
				nzLU2 = 0;
#endif
#if defined(ILUC2) && !defined(ILUC2_SCHUR)
				L2_Entries.clear();
				U2_Entries.clear();
				L2_Address.clear();
				U2_Address.clear();
#endif

				/// FACTORIZATION COMPLETE
				tlocal = Timer() - tlocal;
				tfactor += tlocal;
				if( verbosity > 1 ) std::cout << "Time " << tlocal << "\n";

				if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
				
				INMOST_DATA_ENUM_TYPE swaps = 0;
				if( verbosity > 1 ) std::cout << "total blocks: " << blocks.size() << std::endl;
				for (int q = 0; q < (int)blocks.size(); ++q)
				{
					INMOST_DATA_ENUM_TYPE rbeg = blocks[q].row_start;
					INMOST_DATA_ENUM_TYPE rend = blocks[q].row_end;
					INMOST_DATA_ENUM_TYPE swaps0 = swaps;
					for (INMOST_DATA_ENUM_TYPE k = rbeg; k < rend; ++k) 
						if (Pivot[k]) 
							swaps++;
					if (!blocks[q].separator)
					{
						if (verbosity > 1)
						{
							std::cout << "swaps before block " << swaps0 << " after " << swaps << " inside " << swaps - swaps0 << std::endl;
							std::cout << "before,block " << blocks[q].row_start << ":" << blocks[q].row_end << " " << blocks[q].col_start << ":" << blocks[q].col_end << std::endl;
						}
						blocks[q].row_start -= swaps0;
						blocks[q].col_start -= swaps0;
						blocks[q].row_end -= swaps;
						blocks[q].col_end -= swaps;
						if (verbosity > 1)  
							std::cout << "after ,block " << blocks[q].row_start << ":" << blocks[q].row_end << " " << blocks[q].col_start << ":" << blocks[q].col_end << std::endl;
					}
					else
					{
						if (verbosity > 1)
						{
							std::cout << "swaps before separator " << swaps0 << " after " << swaps << " inside " << swaps - swaps0 << std::endl;
							std::cout << "before,separator " << blocks[q].row_start << ":" << blocks[q].row_end << " " << blocks[q].col_start << ":" << blocks[q].col_end << std::endl;
						}
						blocks[q].row_start -= swaps0;
						if (verbosity > 1) 
							std::cout << "after ,separator " << blocks[q].row_start << ":" << blocks[q].row_end << " " << blocks[q].col_start << ":" << blocks[q].col_end << std::endl;
					}
				}
				level_blocks.push_back(blocks);

				//  Check that we have to enter the next level due to pivoting
				if( swaps )
				{
					tlschur = tlocal = Timer();
					cend = wend-swaps;
					if (verbosity > 1)
					{
						std::cout << "S block: " << wbeg << ":" << wend << std::endl;
						std::cout << "B block: " << cbeg << ":" << cend << std::endl;
					}

					if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

					std::vector< std::vector<Sparse::Row::entry> > C_Entries(Threads());
					interval<INMOST_DATA_ENUM_TYPE, Interval> C_Address(cend, wend);
					
					{ //reorderings

						tlocal = Timer();
						interval<INMOST_DATA_ENUM_TYPE, bool> donePQ(wbeg, wend, false);
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> invP(wbeg, wend), invQ(wbeg, wend);
						interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> localP(wbeg, wend, ENUMUNDEF), localQ(wbeg, wend, ENUMUNDEF);


						if( verbosity > 1 )
						{
							std::cout << "Total swaps: " << swaps << " interval: " << cend << " " << wend << std::endl;
							std::cout << "Reassemble pivots\n";
						}

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

						level_size.push_back(cend - wbeg);
						{
							INMOST_DATA_ENUM_TYPE i = wbeg;
							//enumerate entries that we keep first
							for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) if (!Pivot[k])
							{
								localP[k] = i;
								localQ[k] = i;
								++i;
							}
							//enumerate entries that we dropped at the end
							for (INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k) if (Pivot[k])
							{
								localP[k] = i;
								localQ[k] = i;
								++i;
							}
						}
						{
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> U(wbeg,wend);
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> V(wbeg,wend);
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
							for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
							{
								U[localP[k]] = DL0[k];
								V[localQ[k]] = DR0[k];
							}
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
							for (INMOST_DATA_INTEGER_TYPE k = wbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
							{
								DL0[k] = U[k];
								DR0[k] = V[k];
							}
						}
						
						//inverse the ordering
						ReorderEF(wbeg, wend, donePQ, localP, localQ);
						inversePQ(wbeg,wend,localP,localQ, invP,invQ);
						// Reassamble matrix
						//here F is stored by columns and E by rows
						E_Address.push_back(new interval<INMOST_DATA_ENUM_TYPE, Interval>(cend, wend));
						F_Address.push_back(new interval<INMOST_DATA_ENUM_TYPE, Interval>(cend, wend));
						//and setup the first and last pointers
						nnzF = nnzE = 0;
						//reassamble upper-right corner of B into F
#if defined(USE_OMP_FACT)
#pragma omp parallel reduction(+:nnzF)
#endif
						{
							int thr = Thread();
							INMOST_DATA_ENUM_TYPE tseg = (wend - cend) / Threads();
							INMOST_DATA_ENUM_TYPE tbeg = cend + tseg * thr;
							INMOST_DATA_ENUM_TYPE tend = tbeg + tseg;
							if (thr == Threads() - 1) tend = wend;
							if (tend > tbeg)
							{
								for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
								{
									F_Address.back()->at(k).thr = thr;
									F_Address.back()->at(k).first = F_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(F_Entries[thr].size());
								}
								//find out the space required for each column of F
								//compute how many columns fit into each thread
								for (INMOST_DATA_ENUM_TYPE k = wbeg; k < cend; ++k)
								{
									INMOST_DATA_ENUM_TYPE j = invP[k];
									for (INMOST_DATA_ENUM_TYPE jt = A_Address[j].first; jt < A_Address[j].last; ++jt)
									{
										INMOST_DATA_ENUM_TYPE i = localQ[A_Entries[A_Address[j].thr][jt].first];
										if (i >= tbeg && i < tend) //put into F
											F_Address.back()->at(i).last++;
									}
								}
								//Establish the addresses
								INMOST_DATA_ENUM_TYPE offset = 0;
								for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
								{
									F_Address.back()->at(k).first += offset;
									F_Address.back()->at(k).last += offset;
									offset += F_Address.back()->at(k).Size();
								}
								//allocate size for F storage
								F_Entries[thr].resize(F_Address.back()->at(tend - 1).last);
								//Reset addresses to advance current fill position
								for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
									F_Address.back()->at(k).last = F_Address.back()->at(k).first;
								//Fill data for each thread
								for (INMOST_DATA_ENUM_TYPE k = wbeg; k < cend; ++k)
								{
									INMOST_DATA_ENUM_TYPE j = invP[k];
									for (INMOST_DATA_ENUM_TYPE jt = A_Address[j].first; jt < A_Address[j].last; ++jt)
									{
										INMOST_DATA_ENUM_TYPE i = localQ[A_Entries[A_Address[j].thr][jt].first];
										if (i >= tbeg && i < tend) //put into F
										{
											nnzF++;
											F_Entries[thr][F_Address.back()->at(i).last++] = Sparse::Row::make_entry(k, A_Entries[A_Address[j].thr][jt].second);
										}
									}
								}
							}
						}
						//reassamble lower-left corner of B into E
#if defined(USE_OMP_FACT)
#pragma omp parallel for reduction(+:nnzE)
#endif
						for(INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
						{
							int thr = Thread();
							INMOST_DATA_ENUM_TYPE j = invP[k];
							E_Address.back()->at(k).thr = thr;
							E_Address.back()->at(k).first = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries[thr].size());
							for (INMOST_DATA_ENUM_TYPE jt = A_Address[j].first; jt < A_Address[j].last; ++jt)
							{
								INMOST_DATA_ENUM_TYPE i = localQ[A_Entries[A_Address[j].thr][jt].first];
								INMOST_DATA_REAL_TYPE u = A_Entries[A_Address[j].thr][jt].second;
								if( i < cend ) //put into E
								{
									nnzE++;
									E_Entries[thr].push_back(Sparse::Row::make_entry(i, u));
								}
							}
							E_Address.back()->at(k).last = static_cast<INMOST_DATA_ENUM_TYPE>(E_Entries[thr].size());
							std::sort(E_Entries[thr].begin() + E_Address.back()->at(k).first, E_Entries[thr].end());
						}
						//reassamble lower-right corner of B into C
#if defined(USE_OMP_FACT)
#pragma omp parallel for
#endif
						for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
						{
							int thr = Thread();
							INMOST_DATA_ENUM_TYPE j = invP[k];
							C_Address[k].thr = thr;
							C_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries[thr].size());
							for (INMOST_DATA_ENUM_TYPE jt = A_Address[j].first; jt < A_Address[j].last; ++jt)
							{
								INMOST_DATA_ENUM_TYPE i = localQ[A_Entries[A_Address[j].thr][jt].first];
								INMOST_DATA_REAL_TYPE u = A_Entries[A_Address[j].thr][jt].second;
								if (i >= cend) //put into A
									C_Entries[thr].push_back(Sparse::Row::make_entry(i, u));
							}
							C_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(C_Entries[thr].size());
						}

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

						A_Address.clear();
						A_Entries.clear();
						{
							std::vector< std::vector<Sparse::Row::entry> > empty;
							A_Entries.swap(empty); //release memory
						}

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
						//  Reassamble L and U
						//first change positions
						for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
						{
							INMOST_DATA_ENUM_TYPE j;
							//if( j < cend )
							j = invQ[k];
							{// L-part
								INMOST_DATA_ENUM_TYPE cnt = 0;
								for (INMOST_DATA_ENUM_TYPE jt = L_Address[j].first; jt < L_Address[j].last; ++jt)
								{
									INMOST_DATA_ENUM_TYPE i = localP[L_Entries[L_Address[j].thr][jt].first];
									if (i < cend)
										L_Entries[L_Address[j].thr][jt].first = i;
									else
									{
										L_Entries[L_Address[j].thr][jt].first = ENUMUNDEF;
										L_Entries[L_Address[j].thr][jt].second = 0.0;
										cnt++;
									}
								}
								//this will move all undefined to the end
								if (L_Address[j].Size())
									std::sort(L_Entries[L_Address[j].thr].begin() + L_Address[j].first, L_Entries[L_Address[j].thr].begin() + L_Address[j].last);
								//this will remove references to undefined
								L_Address[j].last -= cnt;
							}
							j = invP[k];
							{// U-part
								INMOST_DATA_ENUM_TYPE cnt = 0;
								for (INMOST_DATA_ENUM_TYPE jt = U_Address[j].first; jt < U_Address[j].last; ++jt)
								{
									INMOST_DATA_ENUM_TYPE Ui = U_Entries[U_Address[j].thr][jt].first;
									INMOST_DATA_ENUM_TYPE i = localQ[Ui];
									if( i < cend )
										U_Entries[U_Address[j].thr][jt].first = i;
									else
									{
										U_Entries[U_Address[j].thr][jt].first = ENUMUNDEF;
										U_Entries[U_Address[j].thr][jt].second = 0.0;
										cnt++;
									}
								}
								//this will move all undefined to the end
								if(U_Address[j].Size() )
									std::sort(U_Entries[U_Address[j].thr].begin()+U_Address[j].first,U_Entries[U_Address[j].thr].begin()+U_Address[j].last);
								//this will remove references to undefined
								U_Address[j].last -= cnt;
							}
						}
#if defined(ILUC2) && defined(ILUC2_SCHUR)
						//now the same for second-order LU
						for(INMOST_DATA_ENUM_TYPE k = wbeg; k < wend; ++k)
						{
							INMOST_DATA_ENUM_TYPE j;
							//if( j < cend )
							j = invQ[k];
							{//L-part
								if (L2_Address[j].Size())
								{
									INMOST_DATA_ENUM_TYPE  cnt = 0;
									for (INMOST_DATA_ENUM_TYPE jt = L2_Address[j].first; jt < L2_Address[j].last; ++jt)
									{
										INMOST_DATA_ENUM_TYPE i = localP[L2_Entries[L2_Address[j].thr][jt].first];
										if (i < cend)
											L2_Entries[L2_Address[j].thr][jt].first = i;
										else
										{
											L2_Entries[L2_Address[j].thr][jt].first = ENUMUNDEF;
											L2_Entries[L2_Address[j].thr][jt].second = 0.0;
											cnt++;
										}
									}
									//this will move all undefined to the end
									if (L2_Address[j].Size())
										std::sort(L2_Entries[L2_Address[j].thr].begin() + L2_Address[j].first, L2_Entries[L2_Address[j].thr].begin() + L2_Address[j].last);
									//this will remove references to undefined
									L2_Address[j].last -= cnt;
								}
							}
							j = invP[k];
							{//U-part
								if( U2_Address[j].Size() )
								{
									INMOST_DATA_ENUM_TYPE cnt = 0;
									for (INMOST_DATA_ENUM_TYPE jt = U2_Address[j].first; jt < U2_Address[j].last; ++jt)
									{
										INMOST_DATA_ENUM_TYPE i = localQ[U2_Entries[U2_Address[j].thr][jt].first];
										if( i < cend )
											U2_Entries[U2_Address[j].thr][jt].first = i;
										else
										{
											U2_Entries[U2_Address[j].thr][jt].first = ENUMUNDEF;
											U2_Entries[U2_Address[j].thr][jt].second = 0.0;
											cnt++;
										}
									}
									//this will move all undefined to the end
									if(U2_Address[j].Size() )
										std::sort(U2_Entries[U2_Address[j].thr].begin()+U2_Address[j].first,U2_Entries[U2_Address[j].thr].begin()+U2_Address[j].last);
									//this will remove references to undefined
									U2_Address[j].last -= cnt;
								}
							}
						}
#endif
						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
						//first establish intervals in arrays
						{
							interval<INMOST_DATA_ENUM_TYPE, Interval> tmpL(wbeg,cend), tmpU(wbeg,cend);
#if defined(ILUC2) && defined(ILUC2_SCHUR)
							interval<INMOST_DATA_ENUM_TYPE, Interval> tmpL2(wbeg, cend), tmpU2(wbeg, cend);
#endif
							interval<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> tmpD(wbeg,cend);
							for(INMOST_DATA_ENUM_TYPE k = wbeg; k < cend; ++k)
							{
								INMOST_DATA_ENUM_TYPE j = invP[k];
								tmpL[k] = L_Address[j];
								tmpU[k] = U_Address[j];
#if defined(ILUC2) && defined(ILUC2_SCHUR)
								tmpL2[k] = L2_Address[j];
								tmpU2[k] = U2_Address[j];
#endif
								tmpD[k] = LU_Diag[j];
							}
							for(INMOST_DATA_ENUM_TYPE k = wbeg; k < cend; ++k)
							{
								L_Address[k] = tmpL[k];
								U_Address[k] = tmpU[k];
#if defined(ILUC2) && defined(ILUC2_SCHUR)
								L2_Address[k] = tmpL2[k];
								U2_Address[k] = tmpU2[k];
#endif
								LU_Diag[k] = tmpD[k];
							}
						}
					
						applyPQ(wbeg, wend, localP, localQ, invP, invQ);

						tlocal = Timer() - tlocal;
						if( verbosity > 1 )
							std::cout << "Time " << tlocal << "\n";

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
					} //reorderings
					int ndrops_lf = 0, ndrops_eu = 0, ndrops_s = 0;
					// Compute Schur
					//result of multilevel preconditioner
					//        |LDU  F |
					//  A  =  |       |
					//        | E   C |
					//
					// For the next step of multilevel factorization we have to compute
					// S = C - (E U^{-1}) D^{-1} (L^{-1} F)
					//
					//first precompute LF block
					{
						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

						std::vector< interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> > LineValuesS(Threads(), interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>(mobeg, moend, 0.0));
						std::vector< interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> > LineIndecesS(Threads(), interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>(mobeg, moend + 1, UNDEF));
						

						std::vector< std::vector<Sparse::Row::entry> > LFt_Entries(Threads());
						interval<INMOST_DATA_ENUM_TYPE, Interval> LFt_Address(wbeg, cend);

						{
							tlocal = Timer();
							std::vector< std::vector<Sparse::Row::entry> > LF_Entries(Threads());
							interval<INMOST_DATA_ENUM_TYPE, Interval> LF_Address(cend, wend);

							
							if (verbosity > 1)
								std::cout << "Assemble LF\n";
							if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
							progress_cur = 0;
							progress_all = wend - cend;
#if defined(USE_OMP_FACT)
#pragma omp parallel for reduction(+:ndrops_lf)
#endif
							for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
							{
								int nthr = Thread();
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
								INMOST_DATA_REAL_TYPE Fmax = 0, Fnorm = 0, Fnum = 0;
								INMOST_DATA_ENUM_TYPE Li = cbeg;
								int Fthr = F_Address.back()->at(k).thr;
								for (INMOST_DATA_ENUM_TYPE jt = F_Address.back()->at(k).first; jt < F_Address.back()->at(k).last; ++jt)
								{
									INMOST_DATA_REAL_TYPE u = F_Entries[Fthr][jt].second;
									INMOST_DATA_ENUM_TYPE j = F_Entries[Fthr][jt].first;
									assert(L_Address[j].Size() >= 1);// at least diagonal
									assert(L_Entries[L_Address[j].thr][L_Address[j].first].first == j);// diagonal entry
									INMOST_DATA_REAL_TYPE ldiag = L_Entries[L_Address[j].thr][L_Address[j].first].second;
									LineValues[j] = u/ldiag;
									Li = LineIndeces[Li] = j + 1;
									Fnorm += u*u;
									Fnum  ++;
									Fmax = std::max<INMOST_DATA_REAL_TYPE>(Fmax, fabs(u));
								}
								LineIndeces[Li] = EOL;
								if( Fnum ) Fnorm = sqrt(Fnorm/Fnum);
								//~ LFdrop = 0;
								// perform solve with L part
								//compute ~F_i = L^{-1} F_i
								Li = LineIndeces[cbeg];
								//Sbeg = cbeg;
								while (Li != EOL)
								{
									INMOST_DATA_ENUM_TYPE curr = Li;
									if ( !check_zero(LineValues[Li - 1]) )
									{
										INMOST_DATA_REAL_TYPE dtol = std::min(tau * tau, tau2) * Fnorm / NuL / NuU; //SCHUR_DROPPING_LF_PREMATURE
										assert(L_Address[Li - 1].Size() >= 1); // at least diagonal
										assert(L_Entries[L_Address[Li - 1].thr][L_Address[Li - 1].first].first == Li - 1); //diagonal goes first
										for (INMOST_DATA_ENUM_TYPE ru = L_Address[Li - 1].first+1; ru < L_Address[Li - 1].last; ++ru)
										{
											INMOST_DATA_REAL_TYPE u = LineValues[Li - 1] * L_Entries[L_Address[Li - 1].thr][ru].second;
											INMOST_DATA_ENUM_TYPE j = L_Entries[L_Address[Li - 1].thr][ru].first;
											assert(L_Address[j].Size() >= 1);// at least diagonal
											assert(L_Entries[L_Address[j].thr][L_Address[j].first].first == j);// diagonal entry
											INMOST_DATA_REAL_TYPE ldiag = L_Entries[L_Address[j].thr][L_Address[j].first].second;
											assert(!check_zero(ldiag));
											u /= ldiag;
											assert(!(u != u));//check nan
											if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
												LineValues[j] -= u;
#if defined(SCHUR_DROPPING_LF_PREMATURE)
											else if (fabs(u) > dtol * std::min(LU_Diag[j],INMOST_DATA_REAL_TYPE(1.0)))
#else
											else if (!check_zero(u))
#endif
											{
												INMOST_DATA_ENUM_TYPE next = curr;
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

												//LineValues[j] = -u;
												//LineIndeces[j + 1] = Sbeg;
												//Sbeg = j+1;
											}
											//else LFdrop += u/LU_Diag[j];
										}
#if defined(ILUC2) && defined(ILUC2_SCHUR)
										// perform solve with second-order L part
										curr = Li;
										for (INMOST_DATA_ENUM_TYPE ru = L2_Address[Li - 1].first; ru < L2_Address[Li - 1].last; ++ru)
										{
											INMOST_DATA_REAL_TYPE u = LineValues[Li - 1] * L2_Entries[L2_Address[Li - 1].thr][ru].second;
											INMOST_DATA_ENUM_TYPE j = L2_Entries[L2_Address[Li - 1].thr][ru].first;
											assert(L_Address[j].Size() >= 1);// at least diagonal
											assert(L_Entries[L_Address[j].thr][L_Address[j].first].first == j);// diagonal entry
											INMOST_DATA_REAL_TYPE ldiag = L_Entries[L_Address[j].thr][L_Address[j].first].second;
											assert(!check_zero(ldiag));
											u /= ldiag;
											assert(!(u != u));//check nan
											if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
												LineValues[j] -= u;
#if defined(SCHUR_DROPPING_LF_PREMATURE)
											else if (fabs(u) > dtol * std::min(LU_Diag[j],INMOST_DATA_REAL_TYPE(1.0)))
#else
											else if (!check_zero(u))
#endif
											{

												INMOST_DATA_ENUM_TYPE next = curr;
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

												//LineValues[j] = -u;
												//LineIndeces[j + 1] = Sbeg;
												//Sbeg = j+1;
											}
											//else LFdrop += u/LU_Diag[j];
										}
#endif
									}
									// solve iteration done
									Li = LineIndeces[Li];
								}
								// Rescale by diagonal
								Li = LineIndeces[cbeg];
								while (Li != EOL)
								{
									LineValues[Li - 1] /= LU_Diag[Li - 1]; //Li - 1 or k?
									Li = LineIndeces[Li];
								}
								// Compute norm for dropping
#if defined(SCHUR_DROPPING_LF)
								INMOST_DATA_REAL_TYPE LFnorm, LFnum, LFmax, LFmin, LFtau;
								LFnorm = LFnum = 0;
								Li = LineIndeces[cbeg];
								LFmax = 0;
								LFmin = 1.0e+54;
								while (Li != EOL)
								{
									INMOST_DATA_REAL_TYPE u = fabs(LineValues[Li - 1]);
									LFnorm += u * u;
									if (u > LFmax) LFmax = u;
									if (u < LFmin) LFmin = u;
									LFnum++;
									Li = LineIndeces[Li];
								}
								if (LFnum) LFnorm = sqrt(LFnorm / LFnum);
								//LFnorm = std::min(1.0,LFnorm);
								LFtau = LFmin + LFnorm*std::min(tau*tau,tau2) / NuL;
								//LFtau = LFmin + std::min(tau*tau,tau2)*LFnorm;// / NuL_max;
								//LFtau = LFmin + (LFmax - LFmin) * pow(LFnorm / LFmax, 2);
								//LFtau = std::min(tau * tau, tau2) * LFnorm / NuL / NuU;
								//LFtau = LFmin + (LFmax - LFmin)*pow(LFnorm/LFmax,4);
#endif
								// Assemble column into matrix
								{
									Li = LineIndeces[cbeg];
									LF_Address[k].thr = nthr;
									LF_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries[nthr].size());
									while (Li != EOL)
									{
										INMOST_DATA_REAL_TYPE u = LineValues[Li - 1];
										assert(!(u != u));//check nan
										assert(Li - 1 >= wbeg && Li - 1 < cend);
#if defined(SCHUR_DROPPING_LF)
										if (fabs(u) >= LFtau)
#else
										if (!check_zero(u))
#endif
											LF_Entries[nthr].push_back(Sparse::Row::make_entry(Li - 1, u));
										else
										{
											//LFdrop += u;
											ndrops_lf++;
										}
										Li = LineIndeces[Li];
									}
									LF_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(LF_Entries[nthr].size());
								}
								// if( LFdrop )
								// {
								// 	LFdrop /= (double) LF_Address[k].Size();
								// 	for(j = LF_Address[k].first; j < LF_Address[k].last; ++j)
								// 		LF_Entries[j].second += LFdrop;
								// }
								//assert(std::is_sorted(LF_Entries.begin()+LF_Address[k].first,LF_Entries.end()));
								// clean linked list
								Li = LineIndeces[cbeg];
								while (Li != EOL)
								{
									INMOST_DATA_ENUM_TYPE j = LineIndeces[Li];
									LineIndeces[Li] = UNDEF;
									LineValues[Li - 1] = 0.0;
									Li = j;
								}
								LineIndeces[cbeg] = UNDEF;
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
								progress_cur++;
								if (verbosity)
								{
									if (k % 100 == 0)
									{
										size_t LFnnz = 0;
										for (size_t q = 0; q < LF_Entries.size(); ++q) LFnnz += LF_Entries[q].size();
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
										{
											std::ios save(NULL);
											save.copyfmt(std::cout);
											//std::cout << "LF " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
											std::cout << "LF " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * progress_cur / (1.f * progress_all);
											std::cout << " nnz " << std::setw(10) << LFnnz << " drops " << std::setw(10) << ndrops_lf;
											std::cout << "\r" << std::flush;
											std::cout.copyfmt(save);
										}
									}
								}
								// iteration done!
							}
							tlocal = Timer() - tlocal;
							if (verbosity > 1)
							{
								size_t LFnnz = 0;
								for (size_t q = 0; q < LF_Entries.size(); ++q) LFnnz += LF_Entries[q].size();
								std::cout << "F nnz " << nnzF << " LF nnz " << LFnnz << " drops " << ndrops_lf << "         \t\t\n";
								std::cout << "Time " << tlocal << "\n";
							}

							if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

							// prepearing LF block for transposed traversal
							tlocal = Timer();
							/*
							{
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Flist(cend, wend);
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE> Fbeg(wbeg, cend);
								tlocal = Timer();
								if (verbosity > 1)
									std::cout << "Transpose LF\n";
								//std::fill(Fbeg.begin() + wbeg - mobeg, Fbeg.begin() + cend - mobeg, EOL);
								for (INMOST_DATA_ENUM_TYPE k = wbeg; k < cend; ++k) Fbeg[k] = EOL;
								for (INMOST_DATA_ENUM_TYPE k = wend; k > cend; --k)
								{
									if (LF_Address[k - 1].Size() > 0)
									{
										INMOST_DATA_ENUM_TYPE i = LF_Entries[LF_Address[k - 1].thr][LF_Address[k - 1].first].first;
										Flist[k - 1] = Fbeg[i];
										Fbeg[i] = k - 1;
									}
									else Flist[k - 1] = EOL;
								}
								//transpone LF matrix
								size_t LFnnz = 0;
								for (size_t q = 0; q < LF_Entries.size(); ++q) LFnnz += LF_Entries[q].size();
								LFt_Entries.resize(LFnnz);

								if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

								{
									INMOST_DATA_ENUM_TYPE j = 0;
									for (INMOST_DATA_ENUM_TYPE k = wbeg; k < cend; k++)
									{
										// assemble line of transposed LF
										LFt_Address[k].first = j;
										INMOST_DATA_ENUM_TYPE Li = Fbeg[k];
										//Fnorms[k] = 0;
										while (Li != EOL)
										{
											LFt_Entries[j].first = Li;
											LFt_Entries[j].second = LF_Entries[LF_Address[Li].thr][LF_Address[Li].first].second;
											j++;
											Li = Flist[Li];
										}
										LFt_Address[k].last = j;
										// shift indices for transposed traversal
										Li = Fbeg[k];
										while (Li != EOL)
										{
											INMOST_DATA_ENUM_TYPE Ui = Flist[Li];
											LF_Address[Li].first++;
											if (LF_Address[Li].Size() > 0)
											{
												INMOST_DATA_ENUM_TYPE i = LF_Entries[LF_Address[Li].thr][LF_Address[Li].first].first;
												Flist[Li] = Fbeg[i];
												Fbeg[i] = Li;
											}
											Li = Ui;
										}
										// iteration done
									}
								}
							}
							*/

#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
							{
								int thr = Thread();
								INMOST_DATA_ENUM_TYPE tseg = (cend - wbeg) / Threads();
								INMOST_DATA_ENUM_TYPE tbeg = wbeg + tseg * thr;
								INMOST_DATA_ENUM_TYPE tend = tbeg + tseg;
								if (thr == Threads() - 1) tend = cend;
								if (tend > tbeg)
								{
									for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
									{
										LFt_Address[k].thr = thr;
										LFt_Address[k].first = LFt_Address[k].last = 0;
									}
									//find out the space required for each column of F
									//compute how many columns fit into each thread
									for (INMOST_DATA_ENUM_TYPE k = cend; k < wend; ++k)
									{
										for (INMOST_DATA_ENUM_TYPE jt = LF_Address[k].first; jt < LF_Address[k].last; ++jt)
										{
											INMOST_DATA_ENUM_TYPE i = LF_Entries[LF_Address[k].thr][jt].first;
											if (i >= tbeg && i < tend) //put into F
												LFt_Address[i].last++;
										}
									}
									//Establish the addresses
									INMOST_DATA_ENUM_TYPE offset = 0;
									for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
									{
										LFt_Address[k].first += offset;
										LFt_Address[k].last += offset;
										offset += LFt_Address[k].Size();
									}
									//allocate size for F storage
									LFt_Entries[thr].resize(LFt_Address[tend - 1].last);
									//Reset addresses to advance current fill position
									for (INMOST_DATA_ENUM_TYPE k = tbeg; k < tend; ++k)
										LFt_Address[k].last = LFt_Address[k].first;
									//Fill data for each thread
									for (INMOST_DATA_ENUM_TYPE k = cend; k < wend; ++k)
									{
										for (INMOST_DATA_ENUM_TYPE jt = LF_Address[k].first; jt < LF_Address[k].last; ++jt)
										{
											INMOST_DATA_ENUM_TYPE i = LF_Entries[LF_Address[k].thr][jt].first;
											if (i >= tbeg && i < tend) //put into F
												LFt_Entries[thr][LFt_Address[i].last++] = Sparse::Row::make_entry(k, LF_Entries[LF_Address[k].thr][jt].second);
										}
									}
								}
							}

							tlocal = Timer() - tlocal;
							if (verbosity > 1)
								std::cout << "Time " << tlocal << "\n";
						}

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

						//DumpMatrix(LFt_Address,LFt_Entries,wbeg,wend,"LFt.mtx");
						// EU and Schur
						std::vector< std::vector<Sparse::Row::entry> > EU_Entries(Threads());
						interval<INMOST_DATA_ENUM_TYPE, Interval> EU_Address(cend, wend);
						tlocal = Timer();
						if (verbosity > 1)
							std::cout << "Assemble EU\n";
						progress_cur = 0;
						progress_all = wend - cend;
#if defined(USE_OMP_FACT)
#pragma omp parallel for reduction(+:ndrops_eu)
#endif
						for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
						{
							int nthr = Thread();
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
							INMOST_DATA_REAL_TYPE Emax = 0, Enorm = 0, Enum = 0;
							// no values at E block - unpack C block into Schur 
							INMOST_DATA_ENUM_TYPE Li = cbeg;
							int Ethr = E_Address.back()->at(k).thr;
							for (INMOST_DATA_ENUM_TYPE j = E_Address.back()->at(k).first; j < E_Address.back()->at(k).last; ++j) //iterate over values of k-th row of E
							{
								INMOST_DATA_REAL_TYPE u = E_Entries[Ethr][j].second;
								INMOST_DATA_ENUM_TYPE i = E_Entries[Ethr][j].first;
								assert(U_Address[i].Size() >= 1);// at least diagonal
								assert(U_Entries[U_Address[i].thr][U_Address[i].first].first == i);// diagonal entry
								INMOST_DATA_REAL_TYPE udiag = U_Entries[U_Address[i].thr][U_Address[i].first].second;
								LineValues[i] = u / udiag;
								Li = LineIndeces[Li] = i + 1;
								Enorm += u*u;
								Enum  ++;
								Emax = std::max(Emax, std::fabs(u));
							}
							LineIndeces[Li] = EOL;
							if( Enum ) Enorm = sqrt(Enorm/Enum);
							//~ EUdrop = 0;
							// perform solve with U part
							// compute ~E_i = E_i U^{-1}
							Li = LineIndeces[cbeg];
							while (Li != EOL)
							{
								INMOST_DATA_ENUM_TYPE curr = Li;
								if (!check_zero(LineValues[Li - 1]))
								{
									INMOST_DATA_REAL_TYPE dtol = std::min(tau*tau,tau2) * Enorm / NuU / NuL; //SCHUR_DROPPING_EU_PREMATURE
									assert(U_Address[Li - 1].Size() >= 1); // at least diagonal
									assert(U_Entries[U_Address[Li - 1].thr][U_Address[Li - 1].first].first == Li - 1); //diagonal goes first
									for (INMOST_DATA_ENUM_TYPE ru = U_Address[Li - 1].first+1; ru < U_Address[Li - 1].last; ++ru)
									{
										INMOST_DATA_REAL_TYPE u = LineValues[Li - 1] * U_Entries[U_Address[Li - 1].thr][ru].second;
										INMOST_DATA_ENUM_TYPE j = U_Entries[U_Address[Li - 1].thr][ru].first;
										assert(U_Address[j].Size() >= 1);// at least diagonal
										assert(U_Entries[U_Address[j].thr][U_Address[j].first].first == j);// diagonal entry
										INMOST_DATA_REAL_TYPE udiag = U_Entries[U_Address[j].thr][U_Address[j].first].second;
										assert(!check_zero(udiag));
										u /= udiag;
										assert(!(u != u));//check nan
										if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
											LineValues[j] -= u;
#if defined(SCHUR_DROPPING_EU_PREMATURE)
										else if (fabs(u) > dtol * std::min(LU_Diag[j],INMOST_DATA_REAL_TYPE(1.0)))
#else
										else if (!check_zero(u))
#endif
										{
											INMOST_DATA_ENUM_TYPE next = curr;
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
										//else EUdrop += u/LU_Diag[j];
									}
#if defined(ILUC2) && defined(ILUC2_SCHUR)
									// perform solve with second-order U part
									curr = Li;
									for (INMOST_DATA_ENUM_TYPE ru = U2_Address[Li - 1].first; ru < U2_Address[Li - 1].last; ++ru)
									{
										INMOST_DATA_REAL_TYPE u = LineValues[Li - 1] * U2_Entries[U2_Address[Li - 1].thr][ru].second;
										INMOST_DATA_ENUM_TYPE j = U2_Entries[U2_Address[Li - 1].thr][ru].first;
										assert(U_Address[j].Size() >= 1);// at least diagonal
										assert(U_Entries[U_Address[j].thr][U_Address[j].first].first == j);// diagonal entry
										INMOST_DATA_REAL_TYPE udiag = U_Entries[U_Address[j].thr][U_Address[j].first].second;
										assert(!check_zero(udiag));
										u /= udiag;
										assert(!(u != u));//check nan
										if (LineIndeces[j + 1] != UNDEF) // There is an entry in the list
											LineValues[j] -= u;
#if defined(SCHUR_DROPPING_EU_PREMATURE)
										else if (fabs(u) > dtol * std::min(LU_Diag[j],INMOST_DATA_REAL_TYPE(1.0)))
#else
										else if (!check_zero(u))
#endif
										{
											INMOST_DATA_ENUM_TYPE next = curr;
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
										//else EUdrop += u/LU_Diag[j];
									}
#endif
								}
								Li = LineIndeces[Li];
								// solve iteration done
							}
							// rescale vector by diagonal *E = ~E_i D^{-1} = E_i U^{-1} D^{-1}
							Li = LineIndeces[cbeg];
							while (Li != EOL)
							{
								assert(!check_zero(LU_Diag[Li - 1]));
								LineValues[Li - 1] /= LU_Diag[Li - 1];
								Li = LineIndeces[Li];
							}
							// compute norm
#if defined(SCHUR_DROPPING_EU)
							INMOST_DATA_REAL_TYPE EUnorm, EUnum, EUmax, EUmin, EUtau;
							EUnorm = EUnum = 0;
							EUmax = 0;
							EUmin = 1.0e+54;
							Li = LineIndeces[cbeg];
							while (Li != EOL)
							{
								INMOST_DATA_REAL_TYPE u = fabs(LineValues[Li - 1]);
								EUnorm += u * u;
								if (u > EUmax) EUmax = u;
								if (u < EUmin) EUmin = u;
								EUnum++;
								Li = LineIndeces[Li];
							}
							if (EUnum) EUnorm = sqrt(EUnorm / EUnum);
							//~ EUtau = EUmin + std::min(tau*tau,tau2)*EUnorm;// / NuU_max;
							EUtau = EUmin + EUnorm * std::min(tau * tau, tau2) / NuU;// / NuU_max;
							//EUnorm = std::min(1.0,EUnorm);
							//EUtau = EUmin + (EUmax - EUmin) * pow(EUnorm / EUmax, 2);
							//EUtau = EUmin + (EUmax - EUmin)*pow(EUnorm/EUmax,4);
#endif
							// drop values that do not satisfy tolerances from linked list of line of EU block
							// Assemble column into matrix
							Li = LineIndeces[cbeg];
							{
								EU_Address[k].thr = nthr;
								EU_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(EU_Entries[nthr].size());
								while (Li != EOL)
								{
									INMOST_DATA_REAL_TYPE u = LineValues[Li - 1];
									assert(!(u != u));//check nan
									assert(Li - 1 >= wbeg && Li - 1 < cend);
#if defined(SCHUR_DROPPING_EU)
									if (fabs(u) >= EUtau)
#else
									if (!check_zero(u))
#endif
										EU_Entries[nthr].push_back(Sparse::Row::make_entry(Li - 1, u));
									else
									{
										//EUdrop += u;
										ndrops_eu++;
									}
									Li = LineIndeces[Li];
								}
								EU_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(EU_Entries[nthr].size());
							}
							// if( EUdrop )
							// {
							// 	EUdrop /= (double)EU_Address[k].Size();
							// 	for(INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
							// 		EU_Entries[r].second += EUdrop;
							// }
							//         clean up linked list 
							Li = LineIndeces[cbeg];
							while (Li != EOL)
							{
								INMOST_DATA_ENUM_TYPE j = LineIndeces[Li];
								LineIndeces[Li] = UNDEF;
								LineValues[Li - 1] = 0.0;
								Li = j;
							}
							LineIndeces[cbeg] = UNDEF;
							if (verbosity)
							{
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
								progress_cur++;
								if (k % 100 == 0)
								{
									size_t EUnnz = 0;
									for (size_t q = 0; q < EU_Entries.size(); ++q) EUnnz += EU_Entries[q].size();
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
									{
										std::ios save(NULL);
										save.copyfmt(std::cout);
										//std::cout << "EU " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
										std::cout << "EU " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * progress_cur / (1.f * progress_all);
										std::cout << " nnz " << std::setw(10) << EUnnz << " drops " << std::setw(10) << ndrops_eu;
										std::cout << "\r" << std::flush;
										std::cout.copyfmt(save);
									}
								}
							}
						}
						tlocal = Timer() - tlocal;
						if (verbosity > 1)
						{
							size_t EUnnz = 0;
							for (size_t q = 0; q < EU_Entries.size(); ++q) EUnnz += EU_Entries[q].size();
							std::cout << "E nnz " << nnzE << " EU nnz " << EUnnz << " drops " << ndrops_eu << "          \t\t\n";
							std::cout << "Time " << tlocal << "\n";
						}

						if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

						{
#if defined(SCHUR_DROPPING_S) || defined(SCHUR_DROPPING_S_PREMATURE)
							// Compute column-norms of schur
							tlocal = Timer();
							interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> Scolmax(cend, wend, 0.0), Srowmax(cend, wend, 0.0);
							{
								std::vector< interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE> > Scolmax_local(Threads(), interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>(cend, wend, 0.0));
								if (verbosity > 1)
									std::cout << "Compute Schur column norm\n";
								progress_cur = 0;
								progress_all = wend - cend;

#if defined(USE_OMP_FACT)
#pragma omp parallel
#endif
								{
									int nthr = Thread();
									interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
									interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
									interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& Scolmaxloc = Scolmax_local[nthr];
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
									for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
									{
										INMOST_DATA_REAL_TYPE v, coef;
										INMOST_DATA_ENUM_TYPE Sbeg = EOL, i, j;
										for (INMOST_DATA_ENUM_TYPE r = C_Address[k].first; r < C_Address[k].last; ++r)
										{
											LineValues[C_Entries[C_Address[k].thr][r].first] = C_Entries[C_Address[k].thr][r].second;
											LineIndeces[C_Entries[C_Address[k].thr][r].first] = Sbeg;
											Sbeg = C_Entries[C_Address[k].thr][r].first;
										}
										for (INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
										{
											coef = -EU_Entries[EU_Address[k].thr][r].second * LU_Diag[EU_Entries[EU_Address[k].thr][r].first];
											AddListUnordered(Sbeg, LFt_Address[EU_Entries[EU_Address[k].thr][r].first], LFt_Entries[LFt_Address[EU_Entries[EU_Address[k].thr][r].first].thr], coef, LineIndeces, LineValues, 0);
										}
										i = Sbeg;
										while (i != EOL)
										{
											j = i;
											i = LineIndeces[i];
											v = fabs(LineValues[j]);
											Srowmax[k] = std::max(Srowmax[k], v);
											Scolmaxloc[j] = std::max(Scolmaxloc[j], v);
											LineValues[j] = 0;
											LineIndeces[j] = UNDEF;
										}
										if (verbosity)
										{
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
											progress_cur++;
											if (k % 100 == 0)
											{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
												{
													std::ios save(NULL);
													save.copyfmt(std::cout);
													//std::cout << "Schur column norm " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
													std::cout << "Schur column norm " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * progress_cur / (1.f * progress_all);
													std::cout << "\t\t\r" << std::flush;
													std::cout.copyfmt(save);
												}
											}
										}
									}
#if defined(USE_OMP_FACT)
#pragma omp for
#endif
									for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
									{
										for (unsigned j = 0; j < Scolmax_local.size(); ++j)
											Scolmax[k] = std::max(Scolmax[k], Scolmax_local[j][k]);
									}
								}
							}
#ifndef NDEBUG
							if (verbosity > 1)
							{
								INMOST_DATA_ENUM_TYPE norow = 0, nocol = 0;
								std::cout << "\n";
								for (INMOST_DATA_ENUM_TYPE k = cend; k < wend; ++k)
								{
									if (check_zero(Scolmax[k]))
									{
										std::cout << "Zero column " << k << " norm " << Scolmax[k] << std::endl;
										nocol++;
									}
									if (check_zero(Srowmax[k]))
									{
										std::cout << "Zero row " << k << " norm " << Srowmax[k] << std::endl;
										norow++;
									}
								}
								std::cout << "Zero norm for " << norow << " rows and " << nocol << " columns " << std::endl;
							}
#endif
							{
								INMOST_DATA_ENUM_TYPE norow = 0, nocol = 0;
								for (INMOST_DATA_ENUM_TYPE k = cend; k < wend; ++k)
								{
									if (check_zero(Scolmax[k])) nocol++;
									if (check_zero(Srowmax[k]))	norow++;
								}
								if (nocol || norow) block_pivot = true;
							}

							tlocal = Timer() - tlocal;
							if (verbosity > 1)
								std::cout << "\nTime " << tlocal << "\n";

							if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
#endif //SCHUR_DROPPING_S
							// Construction of Schur complement 
							std::vector< std::vector<Sparse::Row::entry> >  S_Entries(Threads());
							interval<INMOST_DATA_ENUM_TYPE, Interval> S_Address(cend, wend);
							std::vector< std::vector<INMOST_DATA_ENUM_TYPE> > indicesS(Threads());

							tlocal = Timer();
							if (verbosity > 1)
								std::cout << "Construct Schur\n";

							if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

							nnzA = 0;
							for (INMOST_DATA_ENUM_TYPE k = cend; k < wend; ++k) 
								nnzA += C_Address[k].Size();
							progress_cur = 0;
							progress_all = wend - cend;
#if defined(USE_OMP_FACT)
#pragma omp parallel for reduction(+:ndrops_s)
#endif
							for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
							{
								int nthr = Thread();
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE>& LineValues = LineValuesS[nthr];
								interval<INMOST_DATA_ENUM_TYPE, INMOST_DATA_ENUM_TYPE>& LineIndeces = LineIndecesS[nthr];
								std::vector<INMOST_DATA_ENUM_TYPE>& indices = indicesS[nthr];
								// unpack line of A matrix to temporary linked list 
								INMOST_DATA_ENUM_TYPE Sbeg = EOL;
								for (INMOST_DATA_ENUM_TYPE r = C_Address[k].first; r < C_Address[k].last; ++r)
								{
									assert(!(C_Entries[C_Address[k].thr][r].second != C_Entries[C_Address[k].thr][r].second));
									LineValues[C_Entries[C_Address[k].thr][r].first] = C_Entries[C_Address[k].thr][r].second;
									LineIndeces[C_Entries[C_Address[k].thr][r].first] = Sbeg;
									Sbeg = C_Entries[C_Address[k].thr][r].first;
								}
								// multiply EU and LF blocks and add to temporary linked list
								for (INMOST_DATA_ENUM_TYPE r = EU_Address[k].first; r < EU_Address[k].last; ++r)
								{
									//iterate over corresponding row of LF, add multiplication to list
									INMOST_DATA_REAL_TYPE l = -EU_Entries[EU_Address[k].thr][r].second * LU_Diag[EU_Entries[EU_Address[k].thr][r].first];
									assert(!(l != l));
									AddListUnordered(Sbeg, LFt_Address[EU_Entries[EU_Address[k].thr][r].first], LFt_Entries[LFt_Address[EU_Entries[EU_Address[k].thr][r].first].thr], l, LineIndeces, LineValues, 0);
								}
								OrderList(Sbeg, LineIndeces, indices);
								INMOST_DATA_ENUM_TYPE Ui = Sbeg, rownnz = 0, rowannz = 0;

								while (Ui != EOL)
								{
									if (!check_zero(LineValues[Ui])) rowannz++;
									rownnz++;
									Ui = LineIndeces[Ui];
								}

								if (rowannz == 0 || rownnz == 0)
								{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
									std::cout << "empty row " << k << " nnz " << rownnz << " annz " << rowannz << std::endl;
								}

								Ui = Sbeg;

								{
									S_Address[k].thr = nthr;
									S_Address[k].first = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries[nthr].size());
									while (Ui != EOL)
									{
										INMOST_DATA_REAL_TYPE u = LineValues[Ui];
										assert(!(u != u));
#if defined(SCHUR_DROPPING_S)
										if (fabs(u) >= std::min(Srowmax[k], Scolmax[Ui]) * std::min(double(std::pow(tau,2)),double(tau2)) / NuL / NuU )
#else
										if (!check_zero(u))
										//if( fabs(u) > std::pow(std::min(pow(tau,2),tau2),1.5) )
#endif
											S_Entries[nthr].push_back(Sparse::Row::make_entry(Ui, u));
										else ndrops_s++;
										Ui = LineIndeces[Ui];
									}
									S_Address[k].last = static_cast<INMOST_DATA_ENUM_TYPE>(S_Entries[nthr].size());

									if (S_Address[k].Size() == 0)
									{
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
										std::cout << "empty assembled row " << k << " nnz " << rownnz << " annz " << rowannz << std::endl;
									}
								}
								ClearList(Sbeg, LineIndeces, LineValues);
								if (verbosity)
								{
#if defined(USE_OMP_FACT)
#pragma omp atomic
#endif
									progress_cur++;
									if (k % 100 == 0)
									{
										size_t Snnz = 0;
										for (size_t q = 0; q < S_Entries.size(); ++q)
											Snnz += S_Entries[q].size();
#if defined(USE_OMP_FACT)
#pragma omp critical
#endif
										{
											std::ios save(NULL);
											save.copyfmt(std::cout);
											//std::cout << "Schur " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * (k - cend + 1) / (1.f * (wend - cend));
											std::cout << "Schur " << std::setw(6) << std::fixed << std::setprecision(2) << 100.f * progress_cur / (1.f * progress_all);
											std::cout << " nnz " << std::setw(10) << Snnz << " drop S " << ndrops_s;
											std::cout << "\t\t\r" << std::flush;
											std::cout.copyfmt(save);
										}
									}
								}
							}
							// Schur complement row done
							size_t Snnz = 0;
							for (size_t q = 0; q < S_Entries.size(); ++q)
								Snnz += S_Entries[q].size();

							tlocal = Timer() - tlocal;
							if (verbosity > 1)
							{
								std::cout << "Schur nnz " << Snnz << " drop S " << ndrops_s << "          \t\t\n";
								std::cout << "Time " << tlocal << "\n";
							}

							if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;

							A_Entries.swap(S_Entries);
							A_Address.swap(S_Address);
							/*
							A_Entries.clear();
							A_Entries.resize(Snnz);
							A_Address.set_interval_beg(cend);
							A_Address.set_interval_end(wend);
							for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(wend); ++k)
							{
								A_Address[k].first = A_Entries.size();
								A_Entries.insert(A_Entries.end(),S_Entries[S_Address[k].thr].begin() + S_Address[k].first, S_Entries[S_Address[k].thr].begin() + S_Address[k].last);
								A_Address[k].last = A_Entries.size();
							}
							*/
						}

					}
					// Cleanup arrays for the next iteration
					if( verbosity > 1 )
						std::cout << "Cleanup\n";

					if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
					for(INMOST_DATA_ENUM_TYPE k = cend; k < wend; ++k)
					{
						U_Address[k].first = U_Address[k].last = ENUMUNDEF;
						L_Address[k].first = L_Address[k].last = ENUMUNDEF;
					}
					if( /*wend-cend < 16 ||*/ A_Entries.empty() )
						block_pivot = true;
					//rescale_b = false;
					//run_mpt = false;
					tlschur = Timer()-tlschur;
					tschur += tlschur;
					wbeg = cend; //there is more work to do
					if( verbosity > 1 )
					{
						std::cout << "Total Schur " << tlschur << "\n";
						std::cout << "next S block: " << wbeg << ":" << wend << " size " << wend-wbeg << std::endl;
					}

					


				}
				else
				{
					if (verbosity > 1)
						std::cout << "no swaps" << std::endl;

					level_size.push_back(wend - wbeg);
					wbeg = wend;
				}
				totswaps += swaps;
			}	
		} //factorization and schur complement scope
		
		
		
		if( rescale_b )
		{
			if( verbosity > 1 )
				std::cout << std::endl << " rescaling block B back " << std::endl;

			for (INMOST_DATA_ENUM_TYPE k = mobeg; k < moend; ++k)
			{
				LU_Diag[k] /= (DL0[k] * DR0[k]);
				assert(!(LU_Diag[k] != LU_Diag[k]));
				for (INMOST_DATA_ENUM_TYPE r = U_Address[k].first; r < U_Address[k].last; ++r)
				{
					U_Entries[U_Address[k].thr][r].second *= (DR0[k] / DR0[U_Entries[U_Address[k].thr][r].first]);
					assert(!(U_Entries[U_Address[k].thr][r].second != U_Entries[U_Address[k].thr][r].second));
				}
				for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first; r < L_Address[k].last; ++r)
				{
					L_Entries[L_Address[k].thr][r].second *= (DL0[k] / DL0[L_Entries[L_Address[k].thr][r].first]);
					assert(!(L_Entries[L_Address[k].thr][r].second != L_Entries[L_Address[k].thr][r].second));
				}
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
					for (INMOST_DATA_ENUM_TYPE k = kbeg; k < kend; ++k)
					{
						if( F_Address[level]->at(k).Size() )
						{
							int Fthr = F_Address[level]->at(k).thr;
							for (r = F_Address[level]->at(k).first; r < F_Address[level]->at(k).last; ++r)
							{
								F_Entries[Fthr][r].second /= DR0[k] * DL0[F_Entries[Fthr][r].first];
								assert(!(F_Entries[Fthr][r].second != F_Entries[Fthr][r].second));
							}
						}
						if( E_Address[level]->at(k).Size() )
						{
							int Ethr = E_Address[level]->at(k).thr;
							for (r = E_Address[level]->at(k).first; r < E_Address[level]->at(k).last; ++r)
							{
								E_Entries[Ethr][r].second /= DL0[k] * DR0[E_Entries[Ethr][r].first];
								assert(!(E_Entries[Ethr][r].second != E_Entries[Ethr][r].second));
							}
						}
					}
					first = last;
				}
			}
		}
		/// RESCALING DONE
		
		level_interval.resize(level_size.size());
		if( !level_size.empty() )
		{
			level_interval[0].first = mobeg;
			level_interval[0].last = mobeg + level_size[0];
			for (INMOST_DATA_ENUM_TYPE k = 1; k < level_size.size(); k++)
			{
				level_interval[k].first = level_interval[k - 1].last;
				level_interval[k].last = level_interval[k].first + level_size[k];
			}
		}
		/// MULTILEVEL FACTORIZATION COMPLETE
		ttotal = Timer() - ttotal;
		if( verbosity > 1 || show_summary )
		{
			std::cout << "total      " << ttotal       << "\n";
			std::cout << "reorder    " << treorder     << " ("; fmt(std::cout,100.0*treorder/ttotal    ,6,2) << "%)\n";
			std::cout << "   nd      " << tnd          << " ("; fmt(std::cout,100.0*tnd/ttotal,6,2) << "%)\n";
			std::cout << "   mpt     " << ttransversal << " ("; fmt(std::cout,100.0*ttransversal/ttotal,6,2) << "%)\n";
			std::cout << "   order   " << torder       << " ("; fmt(std::cout,100.0*torder/ttotal    ,6,2) << "%)\n";
			std::cout << "reassamble " << treassamble  << " ("; fmt(std::cout,100.0*treassamble/ttotal ,6,2) << "%)\n";
			std::cout << "rescale    " << trescale     << " ("; fmt(std::cout,100.0*trescale/ttotal    ,6,2) << "%)\n";
			std::cout << "factor     " << tfactor      << " ("; fmt(std::cout,100.0*tfactor/ttotal     ,6,2) << "%)\n";
			std::cout << "   schur   " << tschur       << " ("; fmt(std::cout,100.0*tschur/ttotal       ,6,2) << "%)\n";
			std::cout << "   cond    " << testimator   << " ("; fmt(std::cout,100.0*testimator/ttotal  ,6,2) << "%)\n";
			std::cout << "nnz A " << nzA << " LU " << nzLU << " LU2 " << nzLU2tot << " swaps " << totswaps << " levels " << level_size.size();
			std::cout << " blocks: ";
			for(size_t q = 0; q < level_blocks.size(); ++q)
				std::cout << " " << level_blocks[q].size();
			std::cout  << "\n";
		}
		if (verbosity > 1 && print_mem) std::cout << __FILE__ << ":" << __LINE__ << " mem " << getCurrentRSS() << " peak " << getPeakRSS() << std::endl;
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
		// info->PrepareVector(div);
		// std::fill(div.Begin(),div.End(),0);
		// for(k = mobeg; k < moend; k++) div[k] = 1.0;
		// info->Accumulate(div);
		// for(k = mobeg; k < moend; k++) div[k] = 1.0/div[k];

		condestL = NuL;
		condestU = NuU;

		//prepare temporal array
		temp.set_interval_beg(vbeg);
		temp.set_interval_end(vend);
		for (INMOST_DATA_ENUM_TYPE k = vbeg; k < vend; ++k) temp[k] = 0.0;
		return true;
	}
	bool MLMTILUC_preconditioner::Finalize()
	{
		init = false;
		L_Address.clear();
		U_Address.clear();
		L_Entries.clear();
		U_Entries.clear();
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
	
	int MLMTILUC_preconditioner::Descend(INMOST_DATA_ENUM_TYPE level, Sparse::Vector & inout)
	{
		INMOST_DATA_ENUM_TYPE cbeg, cend;
		//current B block interval
		cbeg = level_interval[level].first;
		cend = level_interval[level].last;
		//check that it makes sense to make restriction step
		if( cend < level_interval.back().last )
		{
			// step 1) obtain ~f
			//apply B^{-1} on f
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) temp[k] = inout[k];
			//Solve with L first
#if defined(USE_OMP)
#pragma omp for schedule(static,1)
#endif
			for (int q = 0; q < (int)level_blocks[level].size(); ++q) if( !level_blocks[level][q].separator )
			{
				assert(level_blocks[level][q].row_start == level_blocks[level][q].col_start);
				assert(level_blocks[level][q].row_end == level_blocks[level][q].col_end);
				INMOST_DATA_ENUM_TYPE bbeg = level_blocks[level][q].row_start;
				INMOST_DATA_ENUM_TYPE bend = level_blocks[level][q].row_end;
				assert(bbeg >= cbeg && bbeg < cend);
				assert(bend > cbeg && bend <= cend);
				for (INMOST_DATA_ENUM_TYPE k = bbeg; k < bend; ++k)
				{
					assert(L_Address[k].Size() >= 1);
					assert(L_Entries[L_Address[k].thr][L_Address[k].first].first == k);
					if (L_Address[k].Size())//iterator over columns of L
					{
						temp[k] /= L_Entries[L_Address[k].thr][L_Address[k].first].second; // scale by L diagonal
						for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first + 1; r < L_Address[k].last; ++r)
							temp[L_Entries[L_Address[k].thr][r].first] -= temp[k] * L_Entries[L_Address[k].thr][r].second;
					}
				}
			}
			//Solve with diagonal
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) temp[k] /= LU_Diag[k];
			//Solve with U
#if defined(USE_OMP)
#pragma omp for schedule(static,1)
#endif
			for (int q = 0; q < (int)level_blocks[level].size(); ++q) if (!level_blocks[level][q].separator)
			{
				assert(level_blocks[level][q].row_start == level_blocks[level][q].col_start);
				assert(level_blocks[level][q].row_end == level_blocks[level][q].col_end);
				INMOST_DATA_ENUM_TYPE bbeg = level_blocks[level][q].row_start;
				INMOST_DATA_ENUM_TYPE bend = level_blocks[level][q].row_end;
				assert(bbeg >= cbeg && bbeg < cend);
				assert(bend > cbeg && bend <= cend);
				for (INMOST_DATA_ENUM_TYPE k = bend; k > bbeg; --k)
				{
					assert(U_Address[k - 1].Size() >= 1);
					assert(U_Entries[U_Address[k - 1].thr][U_Address[k - 1].first].first == k - 1);
					if (U_Address[k - 1].Size()) //iterator over rows of U
					{
						for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first + 1; r < U_Address[k - 1].last; ++r)
							temp[k - 1] -= temp[U_Entries[U_Address[k - 1].thr][r].first] * U_Entries[U_Address[k - 1].thr][r].second;
						temp[k - 1] /= U_Entries[U_Address[k - 1].thr][U_Address[k - 1].first].second;
					}
				}
			}
			//now can calculate ~g
			//multiply ~f by row of E
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(level_interval.back().last); ++k) if( E_Address[level]->at(k).Size() )
			{
				int Ethr = E_Address[level]->at(k).thr;
				for (INMOST_DATA_ENUM_TYPE m = E_Address[level]->at(k).first; m < E_Address[level]->at(k).last; ++m)
				{
					assert(E_Entries[Ethr][m].first < cend);
					inout[k] -= temp[E_Entries[Ethr][m].first] * E_Entries[Ethr][m].second;
				}
			}
		}
		return level + 1;
	}
	int MLMTILUC_preconditioner::Ascend(INMOST_DATA_ENUM_TYPE level, Sparse::Vector & inout)
	{
		INMOST_DATA_ENUM_TYPE cbeg,cend;
		cbeg = level_interval[level-1].first;
		cend = level_interval[level-1].last;
		//calculate Fy, iterate over ~g, write to unused input vector
		//for (k = cbeg; k < cend; ++k)
		//{
		//	for (m = F_Address[k].first; m < F_Address[k].last; m++)
		//		output[k] -= output[F_Entries[m].first] * F_Entries[m].second;
		//}
		//check that it makes sense to make prolongation step
		if (cend < static_cast<INMOST_DATA_INTEGER_TYPE>(level_interval.back().last))
		{
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = cend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(level_interval.back().last); ++k) if (F_Address[level-1]->at(k).Size())
			{
				int Fthr = F_Address[level - 1]->at(k).thr;
				for (INMOST_DATA_ENUM_TYPE m = F_Address[level-1]->at(k).first; m < F_Address[level-1]->at(k).last; ++m)
				{
					assert(F_Entries[Fthr][m].first < cend);
					INMOST_DATA_REAL_TYPE v = inout[k] * F_Entries[Fthr][m].second;
					INMOST_DATA_ENUM_TYPE i = F_Entries[Fthr][m].first;
#if defined(USE_OMP)
#pragma omp atomic
#endif
					inout[i] -= v;
				}
			}
		}
		//perform solve over calculated vector
		//Solve with L first
#if defined(USE_OMP)
#pragma omp for schedule(static,1)
#endif
		for (int q = 0; q < (int)level_blocks[level-1].size(); ++q) if (!level_blocks[level-1][q].separator)
		{
			assert(level_blocks[level-1][q].row_start == level_blocks[level-1][q].col_start);
			assert(level_blocks[level-1][q].row_end == level_blocks[level-1][q].col_end);
			INMOST_DATA_ENUM_TYPE bbeg = level_blocks[level-1][q].row_start;
			INMOST_DATA_ENUM_TYPE bend = level_blocks[level-1][q].row_end;
			assert(bbeg >= cbeg && bbeg < cend);
			assert(bend > cbeg && bend <= cend);
			for (INMOST_DATA_ENUM_TYPE k = bbeg; k < bend; ++k)
			{
				assert(L_Address[k].Size() >= 1);
				assert(L_Entries[L_Address[k].thr][L_Address[k].first].first == k);
				inout[k] /= L_Entries[L_Address[k].thr][L_Address[k].first].second; // scale by L diagonal
				if (L_Address[k].Size())//iterator over columns of L
				{
					for (INMOST_DATA_ENUM_TYPE r = L_Address[k].first + 1; r < L_Address[k].last; ++r)
						inout[L_Entries[L_Address[k].thr][r].first] -= inout[k] * L_Entries[L_Address[k].thr][r].second; //r->first alwayse > k
				}
			}
		}
		//Solve with diagonal
#if defined(USE_OMP)
#pragma omp for
#endif
		for (INMOST_DATA_INTEGER_TYPE k = cbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(cend); ++k) inout[k] /= LU_Diag[k];
		//Solve with U
#if defined(USE_OMP)
#pragma omp for schedule(static,1)
#endif
		for (int q = 0; q < (int)level_blocks[level-1].size(); ++q) if (!level_blocks[level-1][q].separator)
		{
			assert(level_blocks[level-1][q].row_start == level_blocks[level-1][q].col_start);
			assert(level_blocks[level-1][q].row_end == level_blocks[level-1][q].col_end);
			INMOST_DATA_ENUM_TYPE bbeg = level_blocks[level-1][q].row_start;
			INMOST_DATA_ENUM_TYPE bend = level_blocks[level-1][q].row_end;
			assert(bbeg >= cbeg && bbeg < cend);
			assert(bend > cbeg && bend <= cend);
			for (INMOST_DATA_ENUM_TYPE k = bend; k > bbeg; --k)
			{
				assert(U_Address[k - 1].Size() >= 1);
				assert(U_Entries[U_Address[k - 1].thr][U_Address[k - 1].first].first == k - 1);
				if (U_Address[k - 1].Size())//iterator over rows of U
				{
					for (INMOST_DATA_ENUM_TYPE r = U_Address[k - 1].first + 1; r < U_Address[k - 1].last; ++r)
						inout[k - 1] -= inout[U_Entries[U_Address[k - 1].thr][r].first] * U_Entries[U_Address[k - 1].thr][r].second; // r->first always > k
					inout[k - 1] /= U_Entries[U_Address[k - 1].thr][U_Address[k - 1].first].second;
				}
			}
		}
		return level-1;
	}
	bool MLMTILUC_preconditioner::Solve(Sparse::Vector & input, Sparse::Vector & output)
	{
		assert(&input != &output);
		//
		{
			INMOST_DATA_ENUM_TYPE mobeg, moend, vbeg, vend;
			info->GetOverlapRegion(info->GetRank(), mobeg, moend);
			info->GetVectorRegion(vbeg, vend);
		
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); k++) temp[k] = input[k];
			
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
			
			for (INMOST_DATA_INTEGER_TYPE k = vbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(mobeg); k++) output[k] = 0;
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k) output[k] = temp[ddP[k]];//*DL0[k];//*DL0[k];
			for (INMOST_DATA_INTEGER_TYPE k = moend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(vend); k++) output[k] = 0;
			
			INMOST_DATA_ENUM_TYPE level = 0;
			while (level < level_size.size()) level = Descend(level, output);
			while (level > 0) level = Ascend(level, output);
			
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k) temp[ddQ[k]] = output[k];//*DR0[k];//*DR0[k];
#if defined(USE_OMP)
#pragma omp for
#endif
			for (INMOST_DATA_INTEGER_TYPE k = mobeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(moend); ++k) output[k] = temp[k];
			
			

			//Restrict additive schwartz (maybe do it outside?)
			//May assamble partition of unity instead of restriction before accumulation
			//assembly should be done instead of initialization
			// for (k = vbeg; k < mobeg; ++k) output[k] = 0;
			// for (k = moend; k < vend; ++k) output[k] = 0;
			for (INMOST_DATA_INTEGER_TYPE k = vbeg; k < static_cast<INMOST_DATA_INTEGER_TYPE>(mobeg); ++k) output[k] = 0;
			for (INMOST_DATA_INTEGER_TYPE k = moend; k < static_cast<INMOST_DATA_INTEGER_TYPE>(vend); ++k) output[k] = 0;
			// for (k = mobeg; k < moend; ++k) output[k] *= div[k];
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

