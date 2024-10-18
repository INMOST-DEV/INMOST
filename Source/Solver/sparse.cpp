#define _CRT_SECURE_NO_WARNINGS
#include "inmost_sparse.h"
#include "inmost_dense.h"
#include <fstream>
#include <sstream>
#include <queue>
#include "../Misc/utils.h"
#include "../IO/io.hpp"

//if both are commented then sparse row sum is computed with merge
//#define USE_UNORDERED_SPA // use sparse accumulator with unordered result
//#define USE_ORDERED_SPA //use sparse accumulator with ordered result

#ifndef USE_UNORDERED_SPA
#define ASSUME_SORTED
#endif //USE_UNRORDERED_SPA

#define LIST_SIZE 65536 //array size for SPA
//#define LIST_SIZE 32 //for testing
#define LIST_UNDEF USHRT_MAX //undefined entry in SPA


namespace INMOST
{
	namespace Sparse
	{

#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		bool _hasRowEntryType = false;

		INMOST_MPI_Type RowEntryType = INMOST_MPI_DATATYPE_NULL;

		INMOST_MPI_Type GetRowEntryType() { return RowEntryType; }

		bool HaveRowEntryType() { return _hasRowEntryType; }

		void CreateRowEntryType()
		{
#if defined(USE_MPI)
			if (!HaveRowEntryType())
			{
				int ierr;
				MPI_Datatype type[2] = { INMOST_MPI_DATA_ENUM_TYPE, INMOST_MPI_DATA_REAL_TYPE };
				int blocklen[2] = { 1, 1 };
				MPI_Aint disp[2];
				disp[0] = offsetof(Sparse::Row::entry, first);
				disp[1] = offsetof(Sparse::Row::entry, second);
				//disp[2] = sizeof(Sparse::Row::entry);
				ierr = MPI_Type_create_struct(2, blocklen, disp, type, &RowEntryType);
				if (ierr != MPI_SUCCESS)
				{
					std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Type_create_struct" << std::endl;
				}
				ierr = MPI_Type_commit(&RowEntryType);
				if (ierr != MPI_SUCCESS)
				{
					std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Type_commit" << std::endl;
				}
			}
#endif
			_hasRowEntryType = true;
		}

		void DestroyRowEntryType()
		{
#if defined(USE_MPI)
			if (HaveRowEntryType())
			{
				MPI_Type_free(&RowEntryType);
				RowEntryType = INMOST_MPI_DATATYPE_NULL;
			}
#endif
			_hasRowEntryType = false;
		}

#define _m0(x)	(x & 0x7FF)
#define _m1(x)	(x >> 11 & 0x7FF)
#define _m2(x)	(x >> 22 )

		inline static  unsigned int flip(unsigned int fp)
		{
			unsigned int mask = -((int)(fp >> 31)) | 0x80000000;
			return fp ^ mask;
		}
		
		
		void RowMerger::AddRow(INMOST_DATA_REAL_TYPE coef, const Row& r)
		{
			if (!r.Empty() && 1.0 + coef != 1.0)
			{
				for(size_t k = 0; k < rows.size(); ++k)
					if (rows[k].first == &r)
					{
						rows[k].second += coef;
						return;
					}
				rows.push_back(std::make_pair(&r,coef));
			}
		}

		RowMerger::RowMerger() 
#if defined(USE_UNORDERED_SPA) || defined(USE_ORDERED_SPA)
			: list(LIST_SIZE, LIST_UNDEF) 
#endif
		{}

		void RowMerger::RetrieveRow(Sparse::Row& r)
		{
			if (rows.empty()) //only leafs are present
			{
				leafs.Sort();
				leafs.Unique();
				INMOST_DATA_ENUM_TYPE s = leafs.Size();
				r.Resize(s);
				for (INMOST_DATA_ENUM_TYPE q = 0; q < s; ++q)
				{
					r.GetIndex(q) = leafs.GetIndex(q);
					r.GetValue(q) = leafs.GetValue(q);
				}
			}
			else
			{
				if (!leafs.Empty())
				{
					leafs.Sort();
					leafs.Unique();
					rows.push_back(std::make_pair(&leafs,1.0));
				}
				if (rows.size() == 1) //should work even if rows[0] points to r
				{
					INMOST_DATA_ENUM_TYPE s = rows[0].first->Size();
					r.Resize(s);
					for (INMOST_DATA_ENUM_TYPE k = 0; k < s; ++k)
					{
						r.GetIndex(k) = rows[0].first->GetIndex(k);
						r.GetValue(k) = rows[0].first->GetValue(k) * rows[0].second;
					}
				}
#if defined(USE_UNORDERED_SPA)
				else if (true) //fixed SPA with fixed window position
				{
					INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end, beg_next;
					INMOST_DATA_ENUM_TYPE size_max = 0, nlists = rows.size();
					INMOST_DATA_ENUM_TYPE ind = 0, ind_last;
					std::vector< INMOST_DATA_ENUM_TYPE> beg_list, end_list;
					pos.resize(rows.size());
					std::fill(pos.begin(), pos.end(), 0);
					for (INMOST_DATA_ENUM_TYPE k = 0; k < nlists; ++k)
					{
						beg = std::min(rows[k].first->GetIndex(0), beg);
						size_max += rows[k].first->Size();
					}
					store.Resize(size_max);
					while (nlists)
					{
						beg = beg - (beg % LIST_SIZE); //shift to fixed position
						end = beg + LIST_SIZE;
						beg_list.push_back(beg);
						end_list.push_back(end);
						beg_next = ENUMUNDEF;
						ind_last = ind;
						//merge within list range
						for (INMOST_DATA_ENUM_TYPE k = 0; k < nlists; ++k) if (pos[k] < rows[k].first->Size())
						{
							INMOST_DATA_ENUM_TYPE q = pos[k];
							while (q < rows[k].first->Size() && rows[k].first->GetIndex(q) < end)
							{
								INMOST_DATA_ENUM_TYPE i = rows[k].first->GetIndex(q) - beg;
								INMOST_DATA_REAL_TYPE v = rows[k].first->GetValue(q) * rows[k].second;
								if (list[i] == LIST_UNDEF)
								{
									if (1.0 + v != 1.0)
									{
										store.GetIndex(ind) = rows[k].first->GetIndex(q);
										store.GetValue(ind) = v;
										list[i] = ind++;
									}
								}
								else store.GetValue(list[i]) += v;
								++q;
							}
							pos[k] = q;
							if (pos[k] < rows[k].first->Size())
								beg_next = std::min(beg_next, rows[k].first->GetIndex(pos[k]));
							else
							{
								if (k != nlists - 1)
								{
									std::swap(pos[k], pos[nlists - 1]);
									std::swap(rows[k], rows[nlists - 1]);
								}
								pos.pop_back();
								rows.pop_back();
								--nlists;
								--k;
							}
						}
						if (beg_next == ENUMUNDEF && nlists)
							std::cout << __FILE__ << ":" << __LINE__ << " oops!" << std::endl;
						//clear list
						for (INMOST_DATA_ENUM_TYPE q = ind_last; q < ind; ++q)
							list[store.GetIndex(q) - beg] = LIST_UNDEF;
						//sort values within range
						//std::sort(store.Begin() + ind_last, store.Begin() + ind);
						beg = beg_next;
					}
					r.Resize(ind);
					for (INMOST_DATA_ENUM_TYPE q = 0; q < ind; ++q)
					{
						r.GetIndex(q) = store.GetIndex(q);
						r.GetValue(q) = store.GetValue(q);
					}
				}
#endif
				else if (rows.size() == 2)
				{
					Row::MergeSortedRows(rows[0].second, *rows[0].first, rows[1].second, *rows[1].first, store);
					INMOST_DATA_ENUM_TYPE s = store.Size();
					r.Resize(s);
					for (INMOST_DATA_ENUM_TYPE q = 0; q < s; ++q)
					{
						r.GetIndex(q) = store.GetIndex(q);
						r.GetValue(q) = store.GetValue(q);
					}
				}
#if defined(USE_ORDERED_SPA)
				else if (true) //fixed SPA with moving window
				{
					INMOST_DATA_ENUM_TYPE beg = ENUMUNDEF, end, beg_next;
					INMOST_DATA_ENUM_TYPE size_max = 0, nlists = rows.size();
					INMOST_DATA_ENUM_TYPE ind = 0, ind_last;
					pos.resize(rows.size());
					std::fill(pos.begin(), pos.end(), 0);
					for (INMOST_DATA_ENUM_TYPE k = 0; k < nlists; ++k)
					{
						beg = std::min(rows[k].first->GetIndex(0), beg);
						size_max += rows[k].first->Size();
					}
					store.Resize(size_max);
					while (nlists)
					{
						end = beg + LIST_SIZE;
						beg_next = ENUMUNDEF;
						ind_last = ind;
						//merge within list range
						for (INMOST_DATA_ENUM_TYPE k = 0; k < nlists; ++k) if (pos[k] < rows[k].first->Size())
						{
							INMOST_DATA_ENUM_TYPE q = pos[k];
							while (q < rows[k].first->Size() && rows[k].first->GetIndex(q) < end)
							{
								INMOST_DATA_ENUM_TYPE i = rows[k].first->GetIndex(q) - beg;
								INMOST_DATA_REAL_TYPE v = rows[k].first->GetValue(q) * rows[k].second;
								if (list[i] == LIST_UNDEF)
								{
									if (1.0 + v != 1.0)
									{
										store.GetIndex(ind) = rows[k].first->GetIndex(q);
										store.GetValue(ind) = v;
										list[i] = ind++;
									}
								}
								else store.GetValue(list[i]) += v;
								++q;
							}
							pos[k] = q;
							if (pos[k] < rows[k].first->Size())
								beg_next = std::min(beg_next, rows[k].first->GetIndex(pos[k]));
							else
							{
								if (k != nlists - 1)
								{
									std::swap(pos[k], pos[nlists - 1]);
									std::swap(rows[k], rows[nlists - 1]);
								}
								pos.pop_back();
								rows.pop_back();
								--nlists;
								--k;
							}
						}
						if (beg_next == ENUMUNDEF && nlists)
							std::cout << __FILE__ << ":" << __LINE__ << " oops!" << std::endl;
						//clear list
						for (INMOST_DATA_ENUM_TYPE q = ind_last; q < ind; ++q)
							list[store.GetIndex(q) - beg] = LIST_UNDEF;
						//sort values within range
						std::sort(store.Begin() + ind_last, store.Begin() + ind);
						beg = beg_next;
					}
					r.Resize(ind);
					for (INMOST_DATA_ENUM_TYPE q = 0; q < ind; ++q)
					{
						r.GetIndex(q) = store.GetIndex(q);
						r.GetValue(q) = store.GetValue(q);
					}
				}
#endif
				else //merge pairs
				{
					INMOST_DATA_ENUM_TYPE k1, k2, nmerge = 0;
					if (merge.size() < rows.size() - 1)
						merge.resize(rows.size() - 1);
					for (INMOST_DATA_ENUM_TYPE k = 0; k < rows.size(); ++k)
						heap.push_back(std::make_pair(rows[k].first->Size(), k));
					std::make_heap(heap.begin(), heap.end(), std::greater<>());
					while (heap.size() > 1)
					{
						std::pop_heap(heap.begin(), heap.end(), std::greater<>());
						k1 = heap.back().second;
						heap.pop_back();
						std::pop_heap(heap.begin(), heap.end(), std::greater<>());
						k2 = heap.back().second;
						heap.pop_back();
						Row::MergeSortedRows(rows[k1].second, *rows[k1].first, rows[k2].second, *rows[k2].first, merge[nmerge]);
						heap.push_back(std::make_pair(merge[nmerge].Size(), static_cast<INMOST_DATA_ENUM_TYPE>(rows.size())));
						std::push_heap(heap.begin(), heap.end(), std::greater<>());
						rows.push_back(std::make_pair(&merge[nmerge], 1.0));
						nmerge++;
					}
					INMOST_DATA_ENUM_TYPE i = heap.back().second;
					INMOST_DATA_ENUM_TYPE s = rows[i].first->Size();
					r.Resize(s);
					for (INMOST_DATA_ENUM_TYPE k = 0; k < s; ++k)
					{
						r.GetIndex(k) = rows[i].first->GetIndex(k);
						r.GetValue(k) = rows[i].first->GetValue(k);
					}
					heap.clear();
				}
			}
		}


		void RowMerger::Clear()
		{
			leafs.Clear();
			rows.clear();
		}

		
		RowMerger::~RowMerger() {}		

		void RowMerger::Multiply(INMOST_DATA_REAL_TYPE coef)
		{
			for (size_t k = 0; k < rows.size(); ++k)
				rows[k].second *= coef;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < leafs.Size(); ++k)
				leafs.GetValue(k) *= coef;
		}

		
////////class HessianRow

		void   HessianRow::RowVec(INMOST_DATA_REAL_TYPE alpha, const Row & rU, INMOST_DATA_REAL_TYPE beta, Row & rJ) const
		{
			INMOST_DATA_ENUM_TYPE end = rU.Size();
			for(INMOST_DATA_ENUM_TYPE i = 0; i < end; i++) rJ.GetValue(i) *= beta;
			for(INMOST_DATA_ENUM_TYPE i = 0; i < end; i++)
			{
				const index & ind = GetIndex(i);
				if( ind.first == ind.second )
					rJ[ind.first] += alpha*rU.get_safe(ind.first)*GetValue(i);
				else
				{
					rJ[ind.first ] += alpha*rU.get_safe(ind.second)*GetValue(i);
					rJ[ind.second] += alpha*rU.get_safe(ind.first )*GetValue(i);
				}
			}
		}

		bool HessianRow::isSorted() const
		{
			for(INMOST_DATA_ENUM_TYPE k = 1; k < Size(); ++k)
				if( GetIndex(k) < GetIndex(k-1) ) return false;
			return true;
		}

		void HessianRow::MergeSortedRows(INMOST_DATA_REAL_TYPE alpha, const HessianRow & left, INMOST_DATA_REAL_TYPE beta, const HessianRow & right, HessianRow & output)
		{
			assert(left.isSorted());
			assert(right.isSorted());
			output.Resize(left.Size()+right.Size());
			INMOST_DATA_ENUM_TYPE i = 0, j = 0, q = 0;
			while( i < left.Size() && j < right.Size() )
			{
				if( left.GetIndex(i) < right.GetIndex(j) )
				{
					output.GetIndex(q) = left.GetIndex(i);
					output.GetValue(q) = alpha*left.GetValue(i);
					++q;
					++i;
				}
				else if( left.GetIndex(i) == right.GetIndex(j) )
				{
					output.GetIndex(q) = left.GetIndex(i);
					output.GetValue(q) = alpha*left.GetValue(i) + beta*right.GetValue(j);
					++q;
					++i;
					++j;
				}
				else //right is smaller
				{
					output.GetIndex(q) = right.GetIndex(j);
					output.GetValue(q) = beta*right.GetValue(j);
					++q;
					++j;
				}
			}
			while( i < left.Size() )
			{
				if( q > 0 && (output.GetIndex(q-1) == left.GetIndex(i)) )
					output.GetValue(q-1) += alpha*left.GetValue(i);
				else
				{
					output.GetIndex(q) = left.GetIndex(i);
					output.GetValue(q) = alpha*left.GetValue(i);
					++q;
				}
				++i;
			}
			while( j < right.Size() )
			{
				if( q > 0 && (output.GetIndex(q-1) == right.GetIndex(j)) )
					output.GetValue(q-1) += beta*right.GetValue(j);
				else
				{
					output.GetIndex(q) = right.GetIndex(j);
					output.GetValue(q) = beta*right.GetValue(j);
					++q;
				}
				++j;
			}
			output.Resize(q);
		}

		void HessianRow::MergeJacobianHessian(INMOST_DATA_REAL_TYPE a, const Row & JL, const Row & JR, INMOST_DATA_REAL_TYPE b, const HessianRow & HL, INMOST_DATA_REAL_TYPE c, const HessianRow & HR, HessianRow & output)
		{
			// merge three sorted arrays at once
			// one of the array is synthesized from JL and JR on the go
			static const entry stub_entry = make_entry(make_index(ENUMUNDEF,ENUMUNDEF),0.0);
			assert(JL.isSorted());
			assert(JR.isSorted());
			assert(HL.isSorted());
			assert(HR.isSorted());
			output.Resize(HL.Size()+HR.Size()+JL.Size()*JR.Size());
			INMOST_DATA_ENUM_TYPE i = 0, j = 0, k1 = 0, l1 = 0, k2 = 0, l2 = 0, q = 0, kk1 = 0, ll2 = 0, r;
			while(kk1 < JL.Size() && JL.GetIndex(kk1) < JR.GetIndex(l1) ) kk1++;
			k1 = kk1;
			while(ll2 < JR.Size() && JL.GetIndex(k2) > JR.GetIndex(ll2) ) ll2++;
			l2 = ll2;
			entry candidate[4] = {stub_entry,stub_entry,stub_entry,stub_entry};
			if( i < HL.Size() )
				candidate[0] = make_entry(HL.GetIndex(i),b*HL.GetValue(i));
			if( j < HR.Size() )
				candidate[1] = make_entry(HR.GetIndex(j),c*HR.GetValue(j));
			if( k1 < JL.Size() && l1 < JR.Size() )
				candidate[2] = make_entry(make_index(JL.GetIndex(k1),JR.GetIndex(l1)),(JL.GetIndex(k1)==JR.GetIndex(l1)?0.5:1)*a*JL.GetValue(k1)*JR.GetValue(l1));
			if( k2 < JL.Size() && l2 < JR.Size() )
				candidate[3] = make_entry(make_index(JL.GetIndex(k2),JR.GetIndex(l2)),(JL.GetIndex(k2)==JR.GetIndex(l2)?0.5:1)*a*JL.GetValue(k2)*JR.GetValue(l2));
			do
			{
				//pick smallest
				r = 0;
				if( candidate[1].first < candidate[r].first ) r = 1;
				if( candidate[2].first < candidate[r].first ) r = 2;
				if( candidate[3].first < candidate[r].first ) r = 3;
				//all candidates are stub - exit
				if( candidate[r].first == stub_entry.first ) break;
				//record selected entry
				if( q > 0 && (output.GetIndex(q-1) == candidate[r].first) )
					output.GetValue(q-1) += candidate[r].second;
				else
				{
					output.GetIndex(q) = candidate[r].first;
					output.GetValue(q) = candidate[r].second;
					++q;
				}
				if( r == 0 ) //update left hessian index
				{
					if( ++i < HL.Size() )
						candidate[0] = make_entry(HL.GetIndex(i),b*HL.GetValue(i));
					else candidate[0] = stub_entry;
				}
				else if( r == 1 ) //update right hessian index
				{
					if( ++j < HR.Size() )
						candidate[1] = make_entry(HR.GetIndex(j),c*HR.GetValue(j));
					else candidate[1] = stub_entry;
				}
				else if( r == 2 ) //update jacobians indexes
				{
					if( ++k1 == JL.Size() )
					{
						++l1;
						if( l1 < JR.Size() )
						{
							while(kk1 < JL.Size() && JL.GetIndex(kk1) < JR.GetIndex(l1) ) kk1++;
							k1 = kk1;
						}
					}
					if( k1 < JL.Size() && l1 < JR.Size() )
						candidate[2] = make_entry(make_index(JL.GetIndex(k1),JR.GetIndex(l1)),(JL.GetIndex(k1)==JR.GetIndex(l1)?0.5:1)*a*JL.GetValue(k1)*JR.GetValue(l1));
					else
						candidate[2] = stub_entry;
				}
				else //update jacobians indexes
				{
					if( ++l2 == JR.Size() )
					{
						++k2;
						if( k2 < JL.Size() )
						{
							while(ll2 < JR.Size() && JL.GetIndex(k2) > JR.GetIndex(ll2) ) ll2++;
							l2 = ll2;
						}
					}
					if( k2 < JL.Size() && l2 < JR.Size() )
						candidate[3] = make_entry(make_index(JL.GetIndex(k2),JR.GetIndex(l2)),(JL.GetIndex(k2)==JR.GetIndex(l2)?0.5:1)*a*JL.GetValue(k2)*JR.GetValue(l2));
					else
						candidate[3] = stub_entry;
				}
			}
			while(true);
			output.Resize(q);
		}

		void HessianRow::MergeJacobianHessian(INMOST_DATA_REAL_TYPE a, const Row & JL, const Row & JR, INMOST_DATA_REAL_TYPE b, const HessianRow & H, HessianRow & output)
		{
			// merge three sorted arrays at once
			// one of the array is synthesized from JL and JR on the go
			static const entry stub_entry = make_entry(make_index(ENUMUNDEF,ENUMUNDEF),0.0);
			assert(JL.isSorted());
			assert(JR.isSorted());
			assert(H.isSorted());
			output.Resize(H.Size()+JL.Size()*JR.Size());
			INMOST_DATA_ENUM_TYPE i = 0, k1 = 0, l1 = 0, q = 0, k2 = 0, l2 = 0, r, ll2 = 0, kk1 = 0;
			while(kk1 < JL.Size() && JL.GetIndex(kk1) < JR.GetIndex(l1) ) kk1++;
			k1 = kk1;
			while(ll2 < JR.Size() && JL.GetIndex(k2) > JR.GetIndex(ll2) ) ll2++;
			l2 = ll2;
			entry candidate[3] = {stub_entry,stub_entry,stub_entry};
			if( i < H.Size() )
				candidate[0] = make_entry(H.GetIndex(i),b*H.GetValue(i));
			if( k1 < JL.Size() && l1 < JR.Size() )
				candidate[1] = make_entry(make_index(JL.GetIndex(k1),JR.GetIndex(l1)),(JL.GetIndex(k1)==JR.GetIndex(l1)?0.5:1)*a*JL.GetValue(k1)*JR.GetValue(l1));
			if( k2 < JL.Size() && l2 < JR.Size() )
				candidate[2] = make_entry(make_index(JL.GetIndex(k2),JR.GetIndex(l2)),(JL.GetIndex(k2)==JR.GetIndex(l2)?0.5:1)*a*JL.GetValue(k2)*JR.GetValue(l2));
			do
			{
				//pick smallest
				r = 0;
				if( candidate[1].first < candidate[r].first ) r = 1;
				if( candidate[2].first < candidate[r].first ) r = 2;
				//all candidates are stub - exit
				if( candidate[r].first == stub_entry.first ) break;
				//record selected entry
				if( q > 0 && (output.GetIndex(q-1) == candidate[r].first) )
				{
					output.GetValue(q-1) += candidate[r].second;
				}
				else
				{
					output.GetIndex(q) = candidate[r].first;
					output.GetValue(q) = candidate[r].second;
					++q;
				}
				if( r == 0 ) //update left hessian index
				{
					if( ++i < H.Size() )
						candidate[0] = make_entry(H.GetIndex(i),b*H.GetValue(i));
					else candidate[0] = stub_entry;
				}
				else if( r == 1 )
				{
					if( ++k1 == JL.Size() )
					{
						++l1;
						if( l1 < JR.Size() )
						{
							while(kk1 < JL.Size() && JL.GetIndex(kk1) < JR.GetIndex(l1) ) kk1++;
							k1 = kk1;
						}
					}
					if( k1 < JL.Size() && l1 < JR.Size() )
						candidate[1] = make_entry(make_index(JL.GetIndex(k1),JR.GetIndex(l1)),(JL.GetIndex(k1)==JR.GetIndex(l1)?0.5:1)*a*JL.GetValue(k1)*JR.GetValue(l1));
					else
						candidate[1] = stub_entry;
				}
				else //update jacobians indexes
				{
					if( ++l2 == JR.Size() )
					{
						++k2;
						if( k2 < JL.Size() )
						{
							while(ll2 < JR.Size() && JL.GetIndex(k2) > JR.GetIndex(ll2) ) ll2++;
							l2 = ll2;
						}
					}
					if( k2 < JL.Size() && l2 < JR.Size() )
						candidate[2] = make_entry(make_index(JL.GetIndex(k2),JR.GetIndex(l2)),(JL.GetIndex(k2)==JR.GetIndex(l2)?0.5:1)*a*JL.GetValue(k2)*JR.GetValue(l2));
					else
						candidate[2] = stub_entry;
				}
			}
			while(true);
			output.Resize(q);
		}

////////class Row
#if defined(USE_SOLVER)
		INMOST_DATA_REAL_TYPE   Row::RowVec(const Vector & x) const
		{
			INMOST_DATA_REAL_TYPE ret = 0;
			INMOST_DATA_ENUM_TYPE end = Size();
			for(INMOST_DATA_ENUM_TYPE i = 0; i < end; i++) ret = ret + x[GetIndex(i)]*GetValue(i);
			return ret;
		}
#endif //USE_SOLVER

		bool Row::CheckUnique()
		{
			for (INMOST_DATA_ENUM_TYPE k = 0; k < Size() - 1; ++k)
				for (INMOST_DATA_ENUM_TYPE l = k + 1; l < Size(); ++l)
					if (GetIndex(k) == GetIndex(l))
						return false;
			//has to be sorted
			//for (INMOST_DATA_ENUM_TYPE k = 0; k < Size() - 1; ++k)
			//	if (GetIndex(k) >= GetIndex(k + 1))
			//		return false;
			return true;
		}

		void Row::Unique()
		{
			INMOST_DATA_ENUM_TYPE k = 0, s = Size(), q = 0;
			while (++k < s)
			{
				if (GetIndex(q) == GetIndex(k))
					GetValue(q) += GetValue(k);
				else
				{
					++q;
					GetIndex(q) = GetIndex(k);
					GetValue(q) = GetValue(k);
				}
			}
			if (s) Resize(q + 1);
		}
		
		void Row::GetInterval(INMOST_DATA_ENUM_TYPE& beg, INMOST_DATA_ENUM_TYPE& end) const
		{
#if defined(ASSUME_SORTED)
			if (!Empty())
			{
				beg = std::min(beg, Begin()->first);
				end = std::max(end, rBegin()->first + 1);
			}
#else
			for (Sparse::Row::const_iterator it = Begin(); it != End(); ++it)
			{
				beg = std::min(beg, it->first);
				end = std::max(end, it->first + 1);
			}
#endif
		}

		bool Row::isSorted() const
		{
			for(INMOST_DATA_ENUM_TYPE k = 1; k < Size(); ++k)
				if( GetIndex(k) < GetIndex(k-1) ) return false;
			return true;
		}
		void Row::MergeSortedRows(INMOST_DATA_REAL_TYPE alpha, const Row & left, INMOST_DATA_REAL_TYPE beta, const Row & right, Row & output)
		{
			assert(left.isSorted());
			assert(right.isSorted());
			output.Resize(left.Size()+right.Size());
			INMOST_DATA_ENUM_TYPE i = 0, j = 0, q = 0;
			MergeSortedRows(alpha, left, i, beta, right, j, output, q);
			output.Resize(q);
		}

		void Row::MergeSortedRows(INMOST_DATA_REAL_TYPE alpha, const Row& left, INMOST_DATA_ENUM_TYPE &i, INMOST_DATA_REAL_TYPE beta, const Row& right, INMOST_DATA_ENUM_TYPE& j, Row& output, INMOST_DATA_ENUM_TYPE & q)
		{
			INMOST_DATA_REAL_TYPE v;
			while (i < left.Size() && j < right.Size())
			{
				if (left.GetIndex(i) < right.GetIndex(j))
				{
					v = alpha * left.GetValue(i);
					if (1.0 + v != 1.0)
					{
						output.GetIndex(q) = left.GetIndex(i);
						output.GetValue(q) = v;
						++q;
					}
					++i;
				}
				else if (left.GetIndex(i) == right.GetIndex(j))
				{
					v = alpha * left.GetValue(i) + beta * right.GetValue(j);
					if (1.0 + v != 1.0)
					{
						output.GetIndex(q) = left.GetIndex(i);
						output.GetValue(q) = v;
						++q;
					}
					++i;
					++j;
				}
				else //right is smaller
				{
					v = beta * right.GetValue(j);
					if (1.0 + v != 1.0)
					{
						output.GetIndex(q) = right.GetIndex(j);
						output.GetValue(q) = v;
						++q;
					}
					++j;
				}
			}
			while (i < left.Size())
			{
				output.GetIndex(q) = left.GetIndex(i);
				output.GetValue(q) = alpha * left.GetValue(i);
				++q;
				++i;
			}
			while (j < right.Size())
			{
				output.GetIndex(q) = right.GetIndex(j);
				output.GetValue(q) = beta * right.GetValue(j);
				++q;
				++j;
			}
		}
#endif //defined(USE_SOLVER) || defined(USE_AUTODIFF)


#if defined(USE_SOLVER)
////////class Vector
		Vector::Vector(std::string _name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm) :data(start,end)
		{
			comm = _comm;
			name = _name;
			is_parallel = false;
		}

		Vector::Vector(const Vector & other) : data(other.data)
		{
			comm = other.comm;
			name = other.name;
			is_parallel = other.is_parallel;
		}

		Vector & Vector::operator =(Vector const & other)
		{
			comm = other.comm;
			data = other.data;
			name = other.name;
			is_parallel = other.is_parallel;
			return *this;
		}

		Vector::~Vector()
		{
		}

		void     Vector::Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg, INMOST_DATA_ENUM_TYPE mend, std::string file_ord)
		{
			char str[16384];
			std::ifstream input(file.c_str());
			if( input.fail() ) throw -1;
			int state = 0, k;
			INMOST_DATA_ENUM_TYPE vec_size, vec_block, ind = 0;
			INMOST_DATA_REAL_TYPE val;
			int size = 1, rank = 0;
#if defined(USE_MPI)
			int flag = 0;
			MPI_Initialized(&flag);
			if( flag && mend == ENUMUNDEF && mbeg == ENUMUNDEF )
			{
				MPI_Comm_rank(GetCommunicator(),&rank);
				MPI_Comm_size(GetCommunicator(),&size);
			}
#endif
			INMOST_DATA_ENUM_TYPE * ord = NULL;
			if (file_ord != "")
			{
				std::ifstream input_ord;
				input_ord.open(file_ord.c_str(), std::ifstream::in);
				if( input_ord.fail() ) throw -2;
				int n;
				input_ord >> n;
				ord = (INMOST_DATA_ENUM_TYPE *)malloc(sizeof(INMOST_DATA_ENUM_TYPE)* n);
				assert(ord != NULL);
				for (int i=0; i<n; i++) input_ord >> ord[i];
				int nbl;
				input_ord >> nbl;
				if( nbl != size ) throw -3;
				int * ibl;
				ibl = (int *) malloc(sizeof(int) * (nbl+1));
				assert(ibl != NULL);
				for (int i=0; i<nbl+1; i++) input_ord >> ibl[i];
				if( mbeg == ENUMUNDEF ) mbeg = ibl[rank];
				if( mend == ENUMUNDEF ) mend = ibl[rank+1];
				for (int i=0; i<n; i++) ord[i] -= 1;
				free(ibl);
				input_ord.close();
				//std::cout<<"Vector::Load(): n="<<n<<" nbl="<<nbl<<" np="<<size<<" id="<<rank<<" mbeg="<<mbeg<<" mend="<<mend<<std::endl;//db
			}

			while( !input.getline(str,16384).eof() )
			{
				k = 0; while( isspace(str[k]) ) k++;
				if( str[k] == '%' || str[k] == '\0' ) continue;
				std::istringstream istr(str+k);
				switch(state)
				{
					case 0:
						istr >> vec_size; state = 1;
						vec_block = vec_size/size;
						if( mbeg == ENUMUNDEF ) mbeg = rank*vec_block;
						if( mend == ENUMUNDEF )
						{
							if( rank == size-1 ) mend = vec_size;
							else mend = mbeg+vec_block;
						}
						SetInterval(mbeg,mend);
						break;
					case 1:
						istr >> val;
						if( ord ) {
							if( ord[ind] >= mbeg && ord[ind] < mend ) data[ord[ind]] = val;
						} else {
							if( ind >= mbeg && ind < mend ) data[ind] = val;
						}
						ind++;
						break;
				}
			}
			input.close();
			if (file_ord != "") free(ord);
		}

		

		static void RecvBigBuffer(void* buffer, INMOST_DATA_BIG_ENUM_TYPE dsize, int orig, int tag, INMOST_MPI_Comm comm)
		{
#if defined(USE_MPI)
			INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0, wdsize = sizeof(char) * dsize;
			int t = 0;
			while (shift != wdsize)
			{
				chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX), wdsize - shift);
				MPI_Recv(static_cast<char*>(buffer) + shift, sizeof(char) * chunk, MPI_CHAR, orig, tag + t++, comm, MPI_STATUS_IGNORE);
				shift += chunk;
			}
#endif
		}

		static void SendBigBuffer(void* buffer, INMOST_DATA_BIG_ENUM_TYPE dsize, int dest, int tag, INMOST_MPI_Comm comm)
		{
#if defined(USE_MPI)
			INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0, wdsize = sizeof(char) * dsize;
			int t = 0;
			while (shift != wdsize)
			{
				chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX), wdsize - shift);
				MPI_Send(static_cast<char*>(buffer) + shift, sizeof(char) * chunk, MPI_CHAR, dest, tag + t++, comm);
				shift += chunk;
			}
#endif
		}


		void     Vector::SaveBinary(std::string file)
		{
			//std::cout << __FUNCTION__ << std::endl;
			int rank = 0, size = 1, compr = 0;
			INMOST_DATA_BIG_ENUM_TYPE vecsize = Size();
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> vecsizes(size, vecsize);
#if defined(USE_MPI)
			MPI_Comm_rank(GetCommunicator(), &rank);
			MPI_Comm_size(GetCommunicator(), &size);
			if (rank == 0) vecsizes.resize(size);
			MPI_Gather(&vecsize, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, &vecsizes[0], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
#endif
			size_t dsize = vecsize * sizeof(INMOST_DATA_REAL_TYPE), zdsize;
			void* buffer = static_cast<void*>(data.begin()), * zbuffer;
			if (zcompress(buffer, dsize, zbuffer, zdsize))
			{
				std::swap(buffer, zbuffer);
				std::swap(dsize, zdsize);
				compr = 1;
			}
			INMOST_DATA_BIG_ENUM_TYPE wdsize = dsize;
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> dsizes(size, wdsize);
#if defined(USE_MPI)
			MPI_Gather(&wdsize, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, &dsizes[0], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
#endif
			if (rank == 0)
			{
				//io_converter<INMOST_DATA_BIG_ENUM_TYPE, INMOST_DATA_REAL_TYPE> buconv;
				std::ofstream fout(file.c_str(), std::ios::binary);
				/*
				buconv.write_iByteOrder(fout);
				buconv.write_iByteSize(fout);
				//todo: use for data unpack
				buconv.write_fByteOrder(fout);
				buconv.write_fByteSize(fout);
				buconv.write_iValue(fout, compr);
				buconv.write_iValue(fout, size);
				for (int k = 0; k < size; ++k)
				{
					buconv.write_iValue(fout, vecsizes[k]);
					buconv.write_iValue(fout, dsizes[k]);
				}
				*/
				fout.write(reinterpret_cast<const char*>(&size), sizeof(int));
				fout.write(reinterpret_cast<const char*>(&vecsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * size);
				fout.write(reinterpret_cast<const char*>(&dsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * size);
				//Write my data
				fout.write(reinterpret_cast<const char*>(&compr), sizeof(int));
				fout.write(reinterpret_cast<const char*>(buffer), sizeof(char) * dsize);
				if (compr) free(buffer);
				//Get and write remote data
#if defined(USE_MPI)
				for (int it = 1; it < size; ++it)
				{
					buffer = malloc(sizeof(char) * dsizes[it]);
					MPI_Recv(&compr, 1, MPI_INT, it, it, GetCommunicator(),MPI_STATUS_IGNORE);
					RecvBigBuffer(buffer, dsizes[it], it, it, GetCommunicator());
					fout.write(reinterpret_cast<const char*>(&compr), sizeof(int));
					fout.write(reinterpret_cast<const char*>(buffer), sizeof(char) * dsizes[it]);
					free(buffer);
				}
#endif

			}
			else
			{
#if defined(USE_MPI)
				MPI_Send(&compr, 1, MPI_INT, 0, rank, GetCommunicator());
				SendBigBuffer(buffer, wdsize, 0, rank, GetCommunicator());
				if (compr) free(buffer);
#endif
			}
		}

		template<typename Type>
		static bool LoadBinaryPiece(void * buffer, size_t dsize, Type * data, INMOST_DATA_BIG_ENUM_TYPE shift, INMOST_DATA_BIG_ENUM_TYPE size, int compr)
		{
			if (compr)
			{
				size_t udsize = sizeof(Type) * size;
				if (!zuncompress(buffer, dsize, &data[shift], udsize))
				{
					std::cout << __FILE__ << ":" << __LINE__ << " cannot unpack" << std::endl;
					return false;
				}
			}
			else std::memcpy(&data[shift], buffer, dsize);
			return true;
		}

		void     Vector::LoadBinary(std::string file)
		{
			int rank = 0, size = 1, compr = 0, npart;
			INMOST_DATA_BIG_ENUM_TYPE vecsize, vecshift = 0;
#if defined(USE_MPI)
			MPI_Comm_rank(GetCommunicator(), &rank);
			MPI_Comm_size(GetCommunicator(), &size);
#endif
			if (rank == 0)
			{
				int rsize;
				//INMOST_DATA_BIG_ENUM_TYPE wrcompr, wrsize;
				//io_converter<INMOST_DATA_BIG_ENUM_TYPE, INMOST_DATA_REAL_TYPE> buconv;
				std::vector< INMOST_DATA_BIG_ENUM_TYPE> vecsizep, vecsizes, dsizes;
				std::vector<int> nparts(size, 0);
				std::ifstream fin(file.c_str(), std::ios::binary);
				/*
				buconv.read_iByteOrder(fin);
				buconv.read_iByteSize(fin);
				buconv.read_fByteOrder(fin);
				buconv.read_fByteSize(fin);
				buconv.read_iValue(fin, wrcompr);
				buconv.read_iValue(fin, wrsize);
				compr = wrcompr;
				rsize = wrsize;
				*/
				fin.read(reinterpret_cast<char*>(rsize), sizeof(int));
				//Read information on data pieces
				vecsizes.resize(rsize);
				dsizes.resize(rsize);
				fin.read(reinterpret_cast<char*>(&vecsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * rsize);
				fin.read(reinterpret_cast<char*>(&dsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * rsize);
				/*
				for (int k = 0; k < rsize; ++k)
				{
					buconv.read_iValue(fin, vecsizes[k]);
					buconv.read_iValue(fin, dsizes[k]);
				}
				*/
				//Distribute data pieces among processes
				int cnt = static_cast<int>(ceil(rsize / static_cast<double>(size))), tot = 0, part = 0;
				vecsizep.resize(size, 0);
				for (int k = 0; k < size; ++k)
				{
					npart = std::min(cnt, rsize - tot);
					for (int q = 0; q < npart; ++q)
						vecsizep[k] += vecsizes[part++];
					nparts[k] = npart;
					tot += npart;
				}
#if defined(USE_MPI)
				MPI_Scatter(&nparts[0], 1, MPI_INT, &npart, 1, MPI_INT, 0, GetCommunicator());
				MPI_Scatter(&vecsizep[0], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, &vecsize, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
				MPI_Exscan(&vecsize, &vecshift, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, MPI_SUM, GetCommunicator());
#endif
				SetInterval(vecshift, vecshift + vecsize);
				//read my own part
				part = 0; //current part number
				for (int k = 0; k < npart; ++k)
				{
					void * buffer = malloc(sizeof(char) * dsizes[part]);
					fin.read(reinterpret_cast<char*>(&compr), sizeof(int));
					fin.read(reinterpret_cast<char*>(buffer), sizeof(char) * dsizes[part]);
					LoadBinaryPiece(buffer, dsizes[part], &data[0], vecshift, vecsizes[part], compr);
					vecshift += vecsizes[part];
					free(buffer);
					part++;
				}
				//read and send others parts
#if defined(USE_MPI)
				for (int k = 1; k < size; ++k)
				{
					for (int q = 0; q < nparts[k]; ++q)
					{
						MPI_Send(&dsizes[part], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, k, k, GetCommunicator());
						MPI_Send(&vecsizes[part], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, k, k, GetCommunicator());
						fin.read(reinterpret_cast<char*>(&compr), sizeof(int));
						MPI_Send(&compr, 1, MPI_INT, k, k, GetCommunicator());
						void* buffer = malloc(sizeof(char) * dsizes[part]);
						fin.read(reinterpret_cast<char*>(buffer), sizeof(char) * dsizes[part]);
						SendBigBuffer(buffer, sizeof(char) * dsizes[part], k, k, GetCommunicator());
						free(buffer);
						part++;
					}
				}
#endif
			}
#if defined(USE_MPI)
			else
			{
				MPI_Scatter(NULL, 1, MPI_INT, &npart, 1, MPI_INT, 0, GetCommunicator());
				MPI_Scatter(NULL, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, &vecsize, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
				MPI_Exscan(&vecsize, &vecshift, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, MPI_SUM, GetCommunicator());
				SetInterval(vecshift, vecshift + vecsize);
				INMOST_DATA_BIG_ENUM_TYPE dsize, vecsizep;
				for (int k = 0; k < npart; ++k)
				{
					MPI_Recv(&dsize, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, rank, GetCommunicator(), MPI_STATUS_IGNORE);
					MPI_Recv(&vecsizep, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, rank, GetCommunicator(), MPI_STATUS_IGNORE);
					MPI_Recv(&compr, 1, MPI_INT, 0, rank, GetCommunicator(), MPI_STATUS_IGNORE);
					void* buffer = malloc(sizeof(char) * dsize);
					RecvBigBuffer(buffer, sizeof(char) * dsize, 0, rank, GetCommunicator());
					LoadBinaryPiece(buffer, dsize, &data[0], vecshift, vecsizep, compr);
					vecshift += vecsizep;
					free(buffer);
				}
			}
#endif
		}

		void     Vector::Save(std::string file)
		{
			INMOST_DATA_ENUM_TYPE vecsize = Size();

#if defined(USE_MPI)
			int rank = 0, size = 1;
			{
				MPI_Comm_rank(GetCommunicator(),&rank);
				MPI_Comm_size(GetCommunicator(),&size);
				INMOST_DATA_ENUM_TYPE temp = vecsize;
				MPI_Allreduce(&temp,&vecsize,1,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator());
			}
#endif
			std::stringstream rhs(std::ios::in | std::ios::out);
			rhs << std::scientific;
			rhs.precision(15);
			for(iterator it = Begin(); it != End(); ++it) rhs << *it << std::endl;
#if defined(USE_MPI) && defined(USE_MPI_FILE) // Use mpi files
			{
				int ierr;
				MPI_File fh;
				MPI_Status stat;
				ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()), MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_File_close(&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				std::vector<INMOST_DATA_BIG_ENUM_TYPE> datasizes(rank ? 1 : size), offsets(rank ? 1 : size);
				std::string buffer = rhs.str();
				INMOST_DATA_BIG_ENUM_TYPE datasize = buffer.size(), offset;
				ierr = MPI_Gather(&datasize,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&datasizes[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator());
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				if( rank == 0 )
				{
					MPI_Offset off;
					std::stringstream header;
					//header << "% vector " << name << std::endl;
					//header << "% is written by INMOST" << std::endl;
					//header << "% by MPI_File_* api" << std::endl;
					header << vecsize << std::endl;
					ierr = MPI_File_write(fh,const_cast<char *>(header.str().c_str()),static_cast<int>(header.str().size()),MPI_CHAR,&stat);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
					ierr = MPI_File_get_position(fh,&off);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
					offsets[0] = off;
					for(int k = 1; k < size; ++k)
						offsets[k] = offsets[k-1] + datasizes[k-1];
				}
				ierr = MPI_Scatter(&offsets[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&offset,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator());
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				if( datasize )
				{
					INMOST_DATA_BIG_ENUM_TYPE shift = 0, chunk;
					while( shift != datasize )
					{
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),datasize-shift);
						ierr = MPI_File_write_at(fh,offset+shift,buffer.c_str()+shift,static_cast<INMOST_MPI_SIZE>(chunk),MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
						shift += chunk;
					}
				}
				ierr = MPI_File_close(&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			}
#elif defined(USE_MPI) //USE_MPI alternative
			std::string senddata = rhs.str(), recvdata;
			INMOST_DATA_BIG_ENUM_TYPE sendsize = static_cast<int>(senddata.size());
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> recvsize(rank ? 1 : size);
			std::vector<MPI_Request> requests;
			int ierr = MPI_Gather(&sendsize,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&recvsize[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator());
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			int max_tag = 32767;
			int flag = 0;
			int * p_max_tag;
#if defined(USE_MPI2)
			MPI_Comm_get_attr(comm,MPI_TAG_UB,&p_max_tag,&flag);
#else //USE_MPI2
			MPI_Attr_get(comm,MPI_TAG_UB,&p_max_tag,&flag);
#endif //USE_MPI2
			if( flag ) max_tag = *p_max_tag;
			if( rank == 0 )
			{
				INMOST_DATA_BIG_ENUM_TYPE sizesum = recvsize[0];
				for(int k = 1; k < size; k++)
					sizesum += recvsize[k];
				recvdata.resize(sizesum);
				memcpy(&recvdata[0],&senddata[0],sizeof(char)*recvsize[0]);
				INMOST_DATA_BIG_ENUM_TYPE offset = recvsize[0];
				for(int k = 1; k < size; k++)
				{
					INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
					int it = 0; // for mpi tag
					while( shift != recvsize[k] )
					{
						MPI_Request req;
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),recvsize[k] - shift);
						ierr = MPI_Irecv(&senddata[offset+shift],chunk,MPI_CHAR,k, (k*1000 + it)%max_tag, GetCommunicator(), &req);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
						requests.push_back(req);
						shift += chunk;
						it++;
						//TODO: remove temporary check
						if( it >= 1000 )
							std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << recvsize[k] << std::endl;
					}
					offset += shift;
				}
			}
			else 
			{
				INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
				int it = 0; // for mpi tag
				while( shift != sendsize )
				{
					MPI_Request req;
					chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),sendsize - shift);
					ierr = MPI_Isend(&senddata[shift],chunk,MPI_CHAR, 0, (rank*1000 + it)%max_tag, GetCommunicator(), &req);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
					requests.push_back(req);
					shift += chunk;
					it++;
					//TODO: remove temporary check
					if( it >= 1000 )
						std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << sendsize << std::endl;
				}
			}
			if( !requests.empty() )
			{
				ierr = MPI_Waitall((int)requests.size(),&requests[0],MPI_STATUSES_IGNORE);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			}
			if( rank == 0 )
			{
				std::fstream output(file.c_str(),std::ios::out);
				output << vecsize << std::endl;
				output << recvdata;
			}
#else
			std::fstream output(file.c_str(),std::ios::out);
			//output << "% vector " << name << std::endl;
			//output << "% is written by INMOST" << std::endl;
			//output << "% by sequential write" << std::endl;
			output << vecsize << std::endl;
			output << rhs.rdbuf();
#endif
		}
		
		INMOST::Matrix<INMOST_DATA_REAL_TYPE> Vector::operator [](const INMOST::AbstractMatrix<INMOST_DATA_INTEGER_TYPE> & rows) const
		{
			INMOST::Matrix<INMOST_DATA_REAL_TYPE> ret(rows.Rows(),rows.Cols());
			for(INMOST_DATA_ENUM_TYPE i = 0; i < rows.Rows(); ++i)
				for(INMOST_DATA_ENUM_TYPE j = 0; j < rows.Cols(); ++j)
					ret(i,j) = data[rows.get(i,j)];
			return ret;
		}
#if defined(USE_AUTODIFF)
		INMOST::Matrix<value_reference> Vector::operator [](const INMOST::AbstractMatrix<INMOST_DATA_INTEGER_TYPE>& rows)
		{
			INMOST::Matrix<value_reference> ret(rows.Rows(), rows.Cols(),value_reference());
			for (INMOST_DATA_ENUM_TYPE i = 0; i < rows.Rows(); ++i)
				for (INMOST_DATA_ENUM_TYPE j = 0; j < rows.Cols(); ++j)
					new (&ret(i, j)) value_reference(data[rows.get(i, j)]);
			return ret;
		}
#endif //USE_AUTODIFF

////////class Matrix
		void Matrix::MatVec(INMOST_DATA_REAL_TYPE alpha, const Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & out) const //y = alpha*A*x + beta * y
		{
			INMOST_DATA_ENUM_TYPE mbeg, mend;
			INMOST_DATA_INTEGER_TYPE ind, imbeg, imend;
			if( out.Empty() )
			{
				INMOST_DATA_ENUM_TYPE vbeg,vend;
				GetInterval(vbeg,vend);
				out.SetInterval(vbeg,vend);
			}
			//CHECK SOMEHOW FOR DEBUG THAT PROVIDED VECTORS ARE OK
			//~ assert(GetFirstIndex() == out.GetFirstIndex());
			//~ assert(Size() == out.Size());
			GetInterval(mbeg,mend);
			imbeg = mbeg;
			imend = mend;
#if defined(USE_OMP)
#pragma omp for private(ind)
#endif
			for(ind = imbeg; ind < imend; ++ind) //iterate rows of matrix
				out[ind] = beta * out[ind] + alpha * (*this)[ind].RowVec(x);
			// outer procedure should update out vector, if needed
		}
		
		
		


		void Matrix::MatVecTranspose(INMOST_DATA_REAL_TYPE alpha, const Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & out) const //y = alpha*A*x + beta * y
		{
			INMOST_DATA_ENUM_TYPE mbeg, mend;
			INMOST_DATA_INTEGER_TYPE ind, imbeg, imend;
			if( out.Empty() )
			{
				INMOST_DATA_ENUM_TYPE vbeg,vend;
				GetInterval(vbeg,vend);
				out.SetInterval(vbeg,vend);
			}
			//CHECK SOMEHOW FOR DEBUG THAT PROVIDED VECTORS ARE OK
			//~ assert(GetFirstIndex() == out.GetFirstIndex());
			//~ assert(Size() == out.Size());
			GetInterval(mbeg,mend);
			imbeg = mbeg;
			imend = mend;
			if( beta ) for(Vector::iterator it = out.Begin(); it != out.End(); ++it) (*it) *= beta;
			else for(Vector::iterator it = out.Begin(); it != out.End(); ++it) (*it) = 0.0;
#if defined(USE_OMP)
#pragma omp for private(ind)
#endif
			for(ind = imbeg; ind < imend; ++ind)
			{
				for(Row::const_iterator it = (*this)[ind].Begin(); it != (*this)[ind].End(); ++it)
					out[it->first] += alpha * x[ind] * it->second;
			}
			// outer procedure should update out vector, if needed
		}


		Matrix::Matrix(std::string _name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm)
				:data(start,end)
		{
			is_parallel = false;
			comm = _comm;
			SetInterval(start,end);
			name = _name;
		}

		Matrix::Matrix(const Matrix & other) :data(other.data)
		{
			comm = other.comm;
			name = other.name;
		}


		Matrix & Matrix::operator =(Matrix const & other)
		{
			comm = other.comm;
			data = other.data;
			name = other.name;
			return *this;
		}

		Matrix::~Matrix() {}

		void      Matrix::MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size)
		{
			INMOST_DATA_ENUM_TYPE i = to + size, j = from + size;
			if( size > 0 && to != from )
				while( j != from ) data[--j].MoveRow(data[--i]);
		}

		void      HessianMatrix::MoveRows(INMOST_DATA_ENUM_TYPE from, INMOST_DATA_ENUM_TYPE to, INMOST_DATA_ENUM_TYPE size)
		{
			INMOST_DATA_ENUM_TYPE i = to + size, j = from + size;
			if( size > 0 && to != from )
				while( j != from ) data[--j].MoveRow(data[--i]);
		}


		void     Matrix::Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg, INMOST_DATA_ENUM_TYPE mend, std::string file_ord)
		{
			char str[16384];
			std::ifstream input(file.c_str());
			if( input.fail() ) throw -1;
			int state = 0, k;
			INMOST_DATA_ENUM_TYPE mat_size, max_lines, row, col, mat_block;
			INMOST_DATA_REAL_TYPE val;
			int size = 1, rank = 0;
#if defined(USE_MPI)
			int flag = 0;
			MPI_Initialized(&flag);
			if( flag && mend == ENUMUNDEF && mbeg == ENUMUNDEF )
			{
				MPI_Comm_rank(GetCommunicator(),&rank);
				MPI_Comm_size(GetCommunicator(),&size);
			}
#endif
			int * ord = NULL;
			if (file_ord != "")
			{
				std::ifstream input_ord;
				input_ord.open(file_ord.c_str(), std::ifstream::in);
				if( input_ord.fail() ) throw -2;
				int n;
				input_ord >> n; // check if( n == line )
				ord = (int *) malloc(sizeof(int) * n);
				assert(ord != NULL);
				for (int i=0; i<n; i++) input_ord >> ord[i];
				int nbl;
				input_ord >> nbl;
				if( nbl != size ) throw -3;
				int * ibl;
				ibl = (int *) malloc(sizeof(int) * (nbl+1));
				assert(ibl != NULL);
				for (int i=0; i<nbl+1; i++) input_ord >> ibl[i];
				if( mbeg == ENUMUNDEF ) mbeg = (unsigned int) ibl[rank];
				if( mend == ENUMUNDEF ) mend = (unsigned int) ibl[rank + 1];
				for (int i=0; i<n; i++) ord[i] -= 1;
				free(ibl);
				input_ord.close();
				//std::cout<<"Matrix::Load(): n="<<n<<" nbl="<<nbl<<" np="<<size<<" id="<<rank<<" mbeg="<<mbeg<<" mend="<<mend<<std::endl;//db
			}

			int line = 0;
			while( !input.getline(str,16384).eof() )
			{
				line++;
				k = 0; while( isspace(str[k]) ) k++;
				//TODO: check for %%MatrixMarket matrix coordinate real general
				if( str[k] == '%' || str[k] == '\0' ) continue;
				std::istringstream istr(str+k);
				switch(state)
				{
					case 0:
						istr >> mat_size >> mat_size >> max_lines; state = 1;
						mat_block = mat_size/size;
						if( mbeg == ENUMUNDEF ) mbeg = rank*mat_block;
						if( mend == ENUMUNDEF )
						{
							if( rank == size-1 ) mend = mat_size;
							else mend = mbeg+mat_block;
						}
						SetInterval(mbeg,mend);
						//~ std::cout << rank << " my interval " << mbeg << ":" << mend << std::endl;
						break;
					case 1:
						istr >> row >> col >> val;
						row--; col--;
						if( ord ) { row = ord[row]; col = ord[col]; }
						if( row >= mbeg && row < mend ) data[row][col] = val;
						break;
				}
			}
			int nonzero = 0;
			for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
			//~ std::cout << rank << " total nonzero " << max_lines << " my nonzero " << nonzero << std::endl;
			input.close();
			if (file_ord != "") free(ord);
		}

		
		void Matrix::SaveBinary(std::string file)
		{
			//std::cout << __FUNCTION__ << std::endl;
			int rank = 0, size = 1, compr[3] = { 0,0,0 }, rcompr[3] = { 0,0,0 };
			INMOST_DATA_BIG_ENUM_TYPE matsize = GetLastIndex() - GetFirstIndex(), nnzsize = Nonzeros();
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> matsizes(size, matsize), nnzsizes(size, nnzsize);
			//Gather data
			std::vector<INMOST_DATA_ENUM_TYPE> ia;
			std::vector<INMOST_DATA_ENUM_TYPE> ja;
			std::vector<INMOST_DATA_REAL_TYPE> va;
			ia.reserve(matsize + 1);
			ja.reserve(nnzsize);
			va.reserve(nnzsize);
			ia.push_back(0);
			nnzsize = 0;
			for (INMOST_DATA_ENUM_TYPE k = GetFirstIndex(); k < GetLastIndex(); ++k)
			{
				for (Row::iterator jt = (*this)[k].Begin(); jt != (*this)[k].End(); ++jt)
				{
					ja.push_back(jt->first);
					va.push_back(jt->second);
					nnzsize++;
				}
				ia.push_back(ja.size());
			}
			if (ia.size() != matsize + 1)
				std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
			if (ja.size() != nnzsize)
				std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
			if (va.size() != nnzsize)
				std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
			//Communicate sizes
#if defined(USE_MPI)
			MPI_Comm_rank(GetCommunicator(), &rank);
			MPI_Comm_size(GetCommunicator(), &size);
			if (rank == 0)
			{
				matsizes.resize(size);
				nnzsizes.resize(size);
			}
			MPI_Gather(&matsize, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, &matsizes[0], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
			MPI_Gather(&nnzsize, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, &nnzsizes[0], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
#endif
			//Compress data
			size_t dsize[3] =
			{
				ia.size() * sizeof(INMOST_DATA_ENUM_TYPE),
				ja.size() * sizeof(INMOST_DATA_ENUM_TYPE),
				va.size() * sizeof(INMOST_DATA_REAL_TYPE)
			}, zdsize;
			void* buffer[3] =
			{
				static_cast<void*>(&ia[0]),
				static_cast<void*>(&ja[0]),
				static_cast<void*>(&va[0])
			}, * zbuffer;
			for (int k = 0; k < 3; ++k)
			{
				if (zcompress(buffer[k], dsize[k], zbuffer, zdsize))
				{
					std::swap(buffer[k], zbuffer);
					std::swap(dsize[k], zdsize);
					compr[k] = 1;
				}
			}
			//Clear compressed memory
			if (compr[0]) { ia.clear();	std::vector<INMOST_DATA_ENUM_TYPE> empty; ia.swap(empty); }
			if (compr[1]) {	ja.clear();	std::vector<INMOST_DATA_ENUM_TYPE> empty; ja.swap(empty); }
			if (compr[2]) {	va.clear();	std::vector<INMOST_DATA_REAL_TYPE> empty; va.swap(empty); }
			//Exchange buffer sizes
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> dsizes[3];
			INMOST_DATA_BIG_ENUM_TYPE wdsize[3] = { dsize[0],dsize[1],dsize[2] };
			for (int k = 0; k < 3; ++k)
			{
				dsizes[k].resize(size, wdsize[k]);
#if defined(USE_MPI)
				MPI_Gather(&wdsize[k], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, &dsizes[k][0], 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
#endif
			}
			if (rank == 0)
			{
				std::ofstream fout(file.c_str(), std::ios::binary);
				fout.write(reinterpret_cast<const char*>(&size), sizeof(int));
				fout.write(reinterpret_cast<const char*>(&matsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * size);
				fout.write(reinterpret_cast<const char*>(&nnzsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * size);
				for (int k = 0; k < 3; ++k)
					fout.write(reinterpret_cast<const char*>(&dsizes[k][0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * size);
				//Write my data
				fout.write(reinterpret_cast<const char*>(compr), sizeof(int) * 3);
				for (int k = 0; k < 3; ++k)
				{
					fout.write(reinterpret_cast<const char*>(buffer[k]), sizeof(char) * dsizes[k][0]);
					if (compr[k]) free(buffer[k]);
				}
				//Get and write remote data
#if defined(USE_MPI)
				for (int it = 1; it < size; ++it)
				{
					//std::cout << rank << " recv from rank " << it << std::endl;
					MPI_Recv(rcompr, 3, MPI_INT, it, it, GetCommunicator(), MPI_STATUS_IGNORE);
					fout.write(reinterpret_cast<const char*>(rcompr), sizeof(int) * 3);
					for (int k = 0; k < 3; ++k)
					{
						buffer[k] = malloc(sizeof(char) * dsizes[k][it]);
						RecvBigBuffer(buffer[k], sizeof(char) * dsizes[k][it], it, it, GetCommunicator());
						fout.write(reinterpret_cast<const char*>(buffer[k]), sizeof(char) * dsizes[k][it]);
						free(buffer[k]);
					}
					//std::cout << rank << " end recv from rank " << it << std::endl;
				}
#endif
			}
#if defined(USE_MPI)
			else
			{
				//std::cout << "send from rank " << rank << std::endl;
				MPI_Send(compr, 3, MPI_INT, 0, rank, GetCommunicator());
				for (int k = 0; k < 3; ++k)
				{
					SendBigBuffer(buffer[k], sizeof(char) * wdsize[k], 0, rank, GetCommunicator());
					if (compr[k]) free(buffer[k]);
				}
				//std::cout << "end send from rank " << rank << std::endl;
			}
#endif
			//std::cout << "End " << __FUNCTION__ << " rank " << rank << std::endl;
		}

		void Matrix::SaveBinaryRaw(std::string file)
		{
			typedef unsigned idx_t;
			int rank = 0, size = 1, compr[3] = { 0,0,0 }, rcompr[3] = { 0,0,0 };
			INMOST_DATA_BIG_ENUM_TYPE matsize = GetLastIndex() - GetFirstIndex(), nnzsize = Nonzeros();
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> matsizes(size, matsize), nnzsizes(size, nnzsize);
			//Gather data
			std::vector<idx_t> ia;
			std::vector<idx_t> ja;
			std::vector<double> va;
			ia.reserve(matsize + 1);
			ja.reserve(nnzsize);
			va.reserve(nnzsize);
			ia.push_back(0);
			nnzsize = 0;
			for (INMOST_DATA_ENUM_TYPE k = GetFirstIndex(); k < GetLastIndex(); ++k)
			{
				for (Row::iterator jt = (*this)[k].Begin(); jt != (*this)[k].End(); ++jt)
				{
					ja.push_back(jt->first);
					va.push_back(jt->second);
					nnzsize++;
				}
				ia.push_back(ja.size());
			}
			if (ia.size() != matsize + 1)
				std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
			if (ja.size() != nnzsize)
				std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
			if (va.size() != nnzsize)
				std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
			std::ofstream output(file.c_str(), std::ios::binary);
			size_t header[3] = { ia.size(),ja.size(),va.size() };
			output.write(reinterpret_cast<const char*>(header), sizeof(size_t) * 3);
			output.write(reinterpret_cast<const char*>(&ia[0]), sizeof(idx_t) * ia.size());
			output.write(reinterpret_cast<const char*>(&ja[0]), sizeof(idx_t) * ja.size());
			output.write(reinterpret_cast<const char*>(&va[0]), sizeof(double) * va.size());
			output.close();
		}


		void     Matrix::LoadBinary(std::string file)
		{
			int rank = 0, size = 1, rsize = 1, compr[3] = { 0,0,0 }, npart = 1, spart = 0;
			INMOST_DATA_BIG_ENUM_TYPE matsize = 0, matshift = 0, sdsize[3];
			std::vector< INMOST_DATA_BIG_ENUM_TYPE> matsizes, nnzsizes;
			std::vector<INMOST_DATA_ENUM_TYPE> ia;
			std::vector<INMOST_DATA_ENUM_TYPE> ja;
			std::vector<INMOST_DATA_REAL_TYPE> va;
#if defined(USE_MPI)
			MPI_Comm_rank(GetCommunicator(), &rank);
			MPI_Comm_size(GetCommunicator(), &size);
#endif
			if (rank == 0)
			{
				std::vector< INMOST_DATA_BIG_ENUM_TYPE> dsizes[3];
				std::vector<int> nparts(size, 0);
				std::ifstream fin(file.c_str(), std::ios::binary);
				if (fin.fail()) throw BadFile;
				fin.read(reinterpret_cast<char*>(&rsize), sizeof(int));
				//Read information on data pieces
				matsizes.resize(rsize);
				nnzsizes.resize(rsize);
				fin.read(reinterpret_cast<char*>(&matsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * rsize);
				fin.read(reinterpret_cast<char*>(&nnzsizes[0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * rsize);
				for (int k = 0; k < 3; ++k)
				{
					dsizes[k].resize(rsize);
					fin.read(reinterpret_cast<char*>(&dsizes[k][0]), sizeof(INMOST_DATA_BIG_ENUM_TYPE) * rsize);
				}
				//Distribute data pieces among processes
				int cnt = static_cast<int>(ceil(rsize / static_cast<double>(size))), tot = 0, part = 0;
				for (int k = 0; k < size; ++k)
				{
					npart = std::min(cnt, rsize - tot);
					nparts[k] = npart;
					tot += npart;
				}
#if defined(USE_MPI)
				MPI_Scatter(&nparts[0], 1, MPI_INT, &npart, 1, MPI_INT, 0, GetCommunicator());
				MPI_Bcast(&rsize, 1, MPI_INT, 0, GetCommunicator());
				MPI_Bcast(&matsizes[0], rsize, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
				MPI_Bcast(&nnzsizes[0], rsize, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
				MPI_Exscan(&npart, &spart, 1, MPI_INT, MPI_SUM, GetCommunicator());
#endif
				for (int p = spart; p < spart + npart; ++p)
					matsize += matsizes[p];
#if defined(USE_MPI)
				MPI_Exscan(&matsize, &matshift, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, MPI_SUM, GetCommunicator());
#endif
				SetInterval(matshift, matshift + matsize);
				//read my own part
				part = 0; //current part number
				for (int p = 0; p < npart; ++p)
				{
					fin.read(reinterpret_cast<char*>(compr), sizeof(int) * 3);
					ia.resize(matsizes[part] + 1);
					ja.resize(nnzsizes[part]);
					va.resize(nnzsizes[part]);
					size_t dsize_max = std::max(dsizes[0][part], std::max(dsizes[1][part], dsizes[2][part]));
					void* buffer = malloc(sizeof(char) * dsize_max);
					fin.read(reinterpret_cast<char*>(buffer), sizeof(char) * dsizes[0][part]);
					LoadBinaryPiece(buffer, dsizes[0][part], &ia[0], 0, ia.size(), compr[0]);
					fin.read(reinterpret_cast<char*>(buffer), sizeof(char) * dsizes[1][part]);
					LoadBinaryPiece(buffer, dsizes[1][part], &ja[0], 0, ja.size(), compr[1]);
					fin.read(reinterpret_cast<char*>(buffer), sizeof(char) * dsizes[2][part]);
					LoadBinaryPiece(buffer, dsizes[2][part], &va[0], 0, va.size(), compr[2]);
					for (INMOST_DATA_BIG_ENUM_TYPE k = 0; k < matsizes[part]; ++k)
					{
						data[matshift + k].Clear();
						for (INMOST_DATA_BIG_ENUM_TYPE q = ia[k]; q < ia[k + 1]; ++q)
							data[matshift + k].Push(ja[q], va[q]);
					}
					matshift += matsizes[part];
					free(buffer);
					part++;
				}
				{ ia.clear(); std::vector<INMOST_DATA_ENUM_TYPE> empty; ia.swap(empty); }
				{ ja.clear(); std::vector<INMOST_DATA_ENUM_TYPE> empty; ja.swap(empty); }
				{ va.clear(); std::vector<INMOST_DATA_REAL_TYPE> empty; va.swap(empty); }
				//read and send others parts
#if defined(USE_MPI)
				for (int it = 1; it < size; ++it)
				{
					for (int q = 0; q < nparts[it]; ++q)
					{
						fin.read(reinterpret_cast<char*>(compr), sizeof(int) * 3);
						MPI_Send(compr, 3, MPI_INT, it, it, GetCommunicator());
						sdsize[0] = dsizes[0][part];
						sdsize[1] = dsizes[1][part];
						sdsize[2] = dsizes[2][part];
						MPI_Send(sdsize, 3, INMOST_MPI_DATA_BIG_ENUM_TYPE, it, it, GetCommunicator());
						INMOST_DATA_BIG_ENUM_TYPE dsize_max = std::max(sdsize[0], std::max(sdsize[1], sdsize[2]));
						void* buffer = malloc(sizeof(char) * dsize_max);
						for (int k = 0; k < 3; ++k)
						{
							fin.read(reinterpret_cast<char*>(buffer), sizeof(char) * dsizes[k][part]);
							SendBigBuffer(buffer, dsizes[k][part], it, it, GetCommunicator());
						}
						free(buffer);
						part++;
					}
				}
#endif
			}
#if defined(USE_MPI)
			else
			{
				MPI_Scatter(NULL, 1, MPI_INT, &npart, 1, MPI_INT, 0, GetCommunicator());
				MPI_Bcast(&rsize, 1, MPI_INT, 0, GetCommunicator());
				matsizes.resize(rsize);
				nnzsizes.resize(rsize);
				MPI_Bcast(&matsizes[0], rsize, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
				MPI_Bcast(&nnzsizes[0], rsize, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, GetCommunicator());
				MPI_Exscan(&npart, &spart, 1, MPI_INT, MPI_SUM, GetCommunicator());
				for (int p = spart; p < spart + npart; ++p)
					matsize += matsizes[p];
#if defined(USE_MPI)
				MPI_Exscan(&matsize, &matshift, 1, INMOST_MPI_DATA_BIG_ENUM_TYPE, MPI_SUM, GetCommunicator());
#endif
				SetInterval(matshift, matshift + matsize);
				for (int part = spart; part < spart + npart; ++part)
				{
					ia.resize(matsizes[part] + 1);
					ja.resize(nnzsizes[part]);
					va.resize(nnzsizes[part]);
					MPI_Recv(compr, 3, MPI_INT, 0, rank, GetCommunicator(), MPI_STATUS_IGNORE);
					MPI_Recv(sdsize, 3, INMOST_MPI_DATA_BIG_ENUM_TYPE, 0, rank, GetCommunicator(), MPI_STATUS_IGNORE);
					INMOST_DATA_BIG_ENUM_TYPE dsize_max = std::max(sdsize[0], std::max(sdsize[1], sdsize[2]));
					void* buffer = malloc(sizeof(char) * dsize_max);
					RecvBigBuffer(buffer, sdsize[0], 0, rank, GetCommunicator());
					LoadBinaryPiece(buffer, sdsize[0], &ia[0], 0, ia.size(), compr[0]);
					RecvBigBuffer(buffer, sdsize[1], 0, rank, GetCommunicator());
					LoadBinaryPiece(buffer, sdsize[1], &ja[0], 0, ja.size(), compr[1]);
					RecvBigBuffer(buffer, sdsize[2], 0, rank, GetCommunicator());
					LoadBinaryPiece(buffer, sdsize[2], &va[0], 0, va.size(), compr[2]);
					for (INMOST_DATA_BIG_ENUM_TYPE k = 0; k < matsizes[part]; ++k)
					{
						data[matshift + k].Clear();
						for (INMOST_DATA_BIG_ENUM_TYPE q = ia[k]; q < ia[k + 1]; ++q)
							data[matshift + k].Push(ja[q], va[q]);
					}
					matshift += matsizes[part];
					free(buffer);
				}
				{ ia.clear(); std::vector<INMOST_DATA_ENUM_TYPE> empty; ia.swap(empty); }
				{ ja.clear(); std::vector<INMOST_DATA_ENUM_TYPE> empty; ja.swap(empty); }
				{ va.clear(); std::vector<INMOST_DATA_REAL_TYPE> empty; va.swap(empty); }
			}
#endif
		}

		void     Matrix::Save(std::string file, const AnnotationService * text)
		{
			INMOST_DATA_ENUM_TYPE matsize = Size(), nonzero = 0, row = GetFirstIndex()+1;
			bool have_annotation = false;
			if( text && !text->Empty() )
			{
				if( text->GetFirstIndex() == GetFirstIndex() &&
					text->GetLastIndex() == GetLastIndex())
					have_annotation = true;
				else
				{
					std::cout << "Size of provided annotation (" << text->GetFirstIndex() << "," << text->GetLastIndex() << ")" << std::endl;
					std::cout << "differs from the size of the matrix (" << GetFirstIndex() << "," << GetLastIndex() << ")" << std::endl;
				}
			}
			for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
#if defined(USE_MPI)
			int rank = 0, size = 1;
			{
				MPI_Comm_rank(GetCommunicator(),&rank);
				MPI_Comm_size(GetCommunicator(),&size);
				INMOST_DATA_ENUM_TYPE temp_two[2] = {matsize,nonzero}, two[2];
				MPI_Allreduce(temp_two,two,2,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator());
				matsize = two[0];
				nonzero = two[1];
			}
#endif
			std::stringstream mtx(std::ios::in | std::ios::out);
			mtx << std::scientific;
			mtx.precision(15);
			for(INMOST_DATA_ENUM_TYPE k = GetFirstIndex(); k < GetLastIndex(); ++k)
			{
				if( have_annotation )
				{
					const std::string & str = text->GetAnnotation(k);
					if( !str.empty() ) mtx << "% " << str << "\n";
				}
				for(Row::iterator jt = (*this)[k].Begin(); jt != (*this)[k].End(); ++jt)
					mtx << row << " " << jt->first+1 << " " << jt->second << "\n";
				++row;
			}
#if defined(USE_MPI) && defined(USE_MPI_FILE) // USE_MPI2?
			{
				char estring[MPI_MAX_ERROR_STRING];
				int ierr, len;
				MPI_File fh;
				MPI_Status stat;
				std::string buffer = mtx.str();
				ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()), MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_File_close(&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_Barrier(GetCommunicator());
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
				if( ierr != MPI_SUCCESS )
				{
					MPI_Error_string(ierr,estring,&len);
					std::cout << estring << std::endl;
					MPI_Abort(GetCommunicator(),__LINE__);
				}
				std::vector<INMOST_DATA_BIG_ENUM_TYPE> datasizes(rank ? 1 : size), offsets(rank ? 1 : size);
				INMOST_DATA_BIG_ENUM_TYPE datasize = buffer.size(), offset;
				ierr = MPI_Gather(&datasize,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&datasizes[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator());
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				if( rank == 0 )
				{
					MPI_Offset off;
					std::stringstream header;
					header << "%%MatrixMarket matrix coordinate real general" << std::endl;
					header << "% matrix " << name << std::endl;
					header << "% is written by INMOST" << std::endl;
					header << "% by MPI_File_* api" << std::endl;
					header << matsize << " " << matsize << " " << nonzero << std::endl;
					//std::string header_data(header.str());
					ierr = MPI_File_write(fh,const_cast<char *>(header.str().c_str()),static_cast<int>(header.str().size()),MPI_CHAR,&stat);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
					ierr = MPI_File_get_position(fh,&off);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
					offsets[0] = off;
					for(int k = 1; k < size; ++k)
						offsets[k] = offsets[k-1] + datasizes[k-1];
				}
				ierr = MPI_Scatter(&offsets[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&offset,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator());
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				if( datasize )
				{
					INMOST_DATA_BIG_ENUM_TYPE shift = 0, chunk;
					while( shift != datasize )
					{
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),datasize-shift);
						ierr = MPI_File_write_at(fh,offset+shift,buffer.c_str()+shift,static_cast<INMOST_MPI_SIZE>(chunk),MPI_CHAR,&stat);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
						shift += chunk;
					}
				}
				ierr = MPI_File_close(&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			}
#elif defined(USE_MPI)//USE_MPI alternative
			std::string senddata = mtx.str(), recvdata;
			INMOST_DATA_BIG_ENUM_TYPE sendsize = senddata.size();
			std::vector<INMOST_DATA_BIG_ENUM_TYPE> recvsize(rank ? 1 : size);
			std::vector<MPI_Request> requests;
			int ierr = MPI_Gather(&sendsize,1,INMOST_MPI_DATA_BIG_ENUM_TYPE,&recvsize[0],1,INMOST_MPI_DATA_BIG_ENUM_TYPE,0,GetCommunicator());
			if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			int max_tag = 32767;
			int flag = 0;
			int * p_max_tag;
#if defined(USE_MPI2)
			MPI_Comm_get_attr(comm,MPI_TAG_UB,&p_max_tag,&flag);
#else //USE_MPI2
			MPI_Attr_get(comm,MPI_TAG_UB,&p_max_tag,&flag);
#endif //USE_MPI2
			if( flag ) max_tag = *p_max_tag;
			if( rank == 0 )
			{
				INMOST_DATA_BIG_ENUM_TYPE sizesum = recvsize[0];
				for(int k = 1; k < size; k++)
					sizesum += recvsize[k];
				recvdata.resize(sizesum);
				memcpy(&recvdata[0],&senddata[0],sizeof(char)*recvsize[0]);
				INMOST_DATA_BIG_ENUM_TYPE offset = recvsize[0];
				for(int k = 1; k < size; k++)
				{
					INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
					int it = 0; // for mpi tag
					while( shift != recvsize[k] )
					{
						MPI_Request req;
						chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),recvsize[k] - shift);
						ierr = MPI_Irecv(&senddata[offset+shift],chunk,MPI_CHAR,k, (k*1000 + it)%max_tag, GetCommunicator(), &req);
						if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
						requests.push_back(req);
						shift += chunk;
						it++;
						//TODO: remove temporary check
						if( it >= 1000 )
							std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << recvsize[k] << std::endl;
					}
					offset += shift;
				}
			}
			else 
			{
				INMOST_DATA_BIG_ENUM_TYPE chunk, shift = 0;
				int it = 0; // for mpi tag
				while( shift != sendsize )
				{
					MPI_Request req;
					chunk = std::min(static_cast<INMOST_DATA_BIG_ENUM_TYPE>(INT_MAX),sendsize - shift);
					ierr = MPI_Isend(&senddata[shift],chunk,MPI_CHAR, 0, (rank*1000 + it)%max_tag, GetCommunicator(), &req);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
					requests.push_back(req);
					shift += chunk;
					it++;
					//TODO: remove temporary check
					if( it >= 1000 )
						std::cout << __FILE__ << ":" << __LINE__ << " too many iterations!!! " << it << " datasize " << sendsize << std::endl;
				}
			}
			if( !requests.empty() )
			{
				ierr = MPI_Waitall((int)requests.size(),&requests[0],MPI_STATUSES_IGNORE);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			}
			
			if( rank == 0 )
			{
				std::fstream output(file.c_str(),std::ios::out);
				output << "%%MatrixMarket matrix coordinate real general" << std::endl;
				output << "% matrix " << name << std::endl;
				output << "% is written by INMOST" << std::endl;
				output << "% by MPI_Gather* api and sequential write" << std::endl;
				output << matsize << " " << matsize << " " << nonzero << std::endl;
				output << recvdata;
			}
#else
			std::fstream output(file.c_str(),std::ios::out);
			output << "%%MatrixMarket matrix coordinate real general" << std::endl;
			output << "% matrix " << name << std::endl;
			output << "% is written by INMOST" << std::endl;
			output << "% by sequential write " << std::endl;
			output << matsize << " " << matsize << " " << nonzero << std::endl;
			output << mtx.rdbuf();
#endif
		}

		void Matrix::Swap(Matrix & other)
		{
			data.swap(other.data);
			name.swap(other.name);
			std::swap(comm,other.comm);
			std::swap(is_parallel,other.is_parallel);
		}

#if defined(USE_OMP)
		bool LockService::HaveLocks() const
		{
			return !locks.empty();
		}
		bool LockService::Lock(INMOST_DATA_ENUM_TYPE row)
		{
			assert( !locks.empty() );
			omp_set_lock(&locks[row]);
			return true;
		}
		bool LockService::TestLock(INMOST_DATA_ENUM_TYPE row)
		{
			assert( !locks.empty() );
			if( omp_test_lock(&locks[row]) )
				return true;
			else
				return false;
		}
		bool LockService::UnLock(INMOST_DATA_ENUM_TYPE row)
		{
			assert( !locks.empty() );
			omp_unset_lock(&locks[row]);
			return true;
		}
		void LockService::SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end)
		{
			if( !locks.empty() ) DestroyLocks();
			if( end != beg )
			{
				locks.set_interval_beg(beg);
				locks.set_interval_end(end);
				for(INMOST_DATA_ENUM_TYPE k = beg; k < end; ++k)
					omp_init_lock(&locks[k]);
			}
		}
		void LockService::DestroyLocks()
		{
			INMOST_DATA_ENUM_TYPE kbeg,kend;
			kbeg = locks.get_interval_beg();
			kend = locks.get_interval_end();
			for(INMOST_DATA_ENUM_TYPE k = kbeg; k < kend; ++k)
				omp_destroy_lock(&locks[k]);
		}
		INMOST_DATA_ENUM_TYPE  LockService::GetFirstIndex() const {return locks.get_interval_beg();}
		INMOST_DATA_ENUM_TYPE  LockService::GetLastIndex() const {return locks.get_interval_end();}
		bool LockService::Empty() const {return locks.empty();}
		void LockService::GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = locks.get_interval_beg(); end = locks.get_interval_end();}
#else
		bool LockService::HaveLocks() const {return false;}
		bool LockService::Lock(INMOST_DATA_ENUM_TYPE row) {(void)row; return true;}
		bool LockService::TestLock(INMOST_DATA_ENUM_TYPE row) {(void)row; return true;}
		bool LockService::UnLock(INMOST_DATA_ENUM_TYPE row) {(void)row; return true;}
		void LockService::SetInterval(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end) {(void)beg; (void)end;}
		void LockService::DestroyLocks() {}
		INMOST_DATA_ENUM_TYPE  LockService::GetFirstIndex() const {return 0;}
		INMOST_DATA_ENUM_TYPE  LockService::GetLastIndex() const {return 0;}
		bool LockService::Empty() const {return true;}
		void LockService::GetInterval(INMOST_DATA_ENUM_TYPE & start, INMOST_DATA_ENUM_TYPE & end) const {start = 0; end = 0;}
#endif


////////class HessianMatrix

		void     HessianMatrix::Load(std::string file, INMOST_DATA_ENUM_TYPE mbeg, INMOST_DATA_ENUM_TYPE mend)
		{
			char str[16384];
			std::ifstream input(file.c_str());
			if( input.fail() ) throw -1;
			int state = 0, k;
			INMOST_DATA_ENUM_TYPE mat_size, max_lines, row, coli, colj, mat_block;
			INMOST_DATA_REAL_TYPE val;
			int size = 1, rank = 0;
#if defined(USE_MPI)
			int flag = 0;
			MPI_Initialized(&flag);
			if( flag && mend == ENUMUNDEF && mbeg == ENUMUNDEF )
			{
				MPI_Comm_rank(GetCommunicator(),&rank);
				MPI_Comm_size(GetCommunicator(),&size);
			}
#endif
			int line = 0;
			while( !input.getline(str,16384).eof() )
			{
				line++;
				k = 0; while( isspace(str[k]) ) k++;
				//TODO: check for %%MatrixMarket hessian coordinate real general
				if( str[k] == '%' || str[k] == '\0' ) continue;
				std::istringstream istr(str+k);
				switch(state)
				{
					case 0:
						istr >> mat_size >> mat_size >> max_lines; state = 1;
						mat_block = mat_size/size;
						if( mbeg == ENUMUNDEF ) mbeg = rank*mat_block;
						if( mend == ENUMUNDEF )
						{
							if( rank == size-1 ) mend = mat_size;
							else mend = mbeg+mat_block;
						}
						SetInterval(mbeg,mend);
						//~ std::cout << rank << " my interval " << mbeg << ":" << mend << std::endl;
						break;
					case 1:
						istr >> row >> coli >> colj >> val;
						row--; coli--; colj--;
						if( row >= mbeg && row < mend ) data[row][HessianRow::make_index(coli,colj)] = val;
						break;
				}
			}
			int nonzero = 0;
			for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
			//~ std::cout << rank << " total nonzero " << max_lines << " my nonzero " << nonzero << std::endl;
			input.close();
		}

		HessianMatrix::HessianMatrix(std::string _name, INMOST_DATA_ENUM_TYPE start, INMOST_DATA_ENUM_TYPE end, INMOST_MPI_Comm _comm)
				:data(start,end)
		{
			is_parallel = false;
			comm = _comm;
			SetInterval(start,end);
			name = _name;
		}

		HessianMatrix::HessianMatrix(const HessianMatrix & other) :data(other.data)
		{
			comm = other.comm;
			name = other.name;
		}


		HessianMatrix & HessianMatrix::operator =(HessianMatrix const & other)
		{
			comm = other.comm;
			data = other.data;
			name = other.name;
			return *this;
		}

		HessianMatrix::~HessianMatrix() {}

		void     HessianMatrix::Save(std::string file, const AnnotationService * text)
		{
			INMOST_DATA_ENUM_TYPE matsize = Size(), nonzero = 0, row = GetFirstIndex()+1;
			bool have_annotation = false;
			if( text && !text->Empty() )
			{
				if( text->GetFirstIndex() == GetFirstIndex() &&
					text->GetLastIndex() == GetLastIndex())
					have_annotation = true;
				else
				{
					std::cout << "Size of provided annotation (" << text->GetFirstIndex() << "," << text->GetLastIndex() << ")" << std::endl;
					std::cout << "differs from the size of the matrix (" << GetFirstIndex() << "," << GetLastIndex() << ")" << std::endl;
				}
			}
			for(iterator it = Begin(); it != End(); ++it) nonzero += it->Size();
#if defined(USE_MPI)
			int rank = 0, size = 1;
			{
				MPI_Comm_rank(GetCommunicator(),&rank);
				MPI_Comm_size(GetCommunicator(),&size);
				INMOST_DATA_ENUM_TYPE temp_two[2] = {matsize,nonzero}, two[2];
				MPI_Allreduce(temp_two,two,2,INMOST_MPI_DATA_ENUM_TYPE,MPI_SUM,GetCommunicator());
				matsize = two[0];
				nonzero = two[1];
			}
#endif
			std::stringstream mtx(std::ios::in | std::ios::out);
			mtx << std::scientific;
			mtx.precision(15);
			for(INMOST_DATA_ENUM_TYPE k = GetFirstIndex(); k < GetLastIndex(); ++k)
			{
				if( have_annotation )
				{
					const std::string & str = text->GetAnnotation(k);
					if( !str.empty() ) mtx << "% " << str << "\n";
				}
				for(HessianRow::iterator jt = (*this)[k].Begin(); jt != (*this)[k].End(); ++jt)
					mtx << row << " " << jt->first.first+1 << " " << jt->first.second+1 << " " << jt->second << "\n";
				++row;
			}
#if defined(USE_MPI) && defined(USE_MPI_FILE) // USE_MPI2?
			{
				int ierr;
				MPI_File fh;
				MPI_Status stat;
				ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()), MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_File_close(&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_File_open(GetCommunicator(),const_cast<char *>(file.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				if( rank == 0 )
				{
					std::stringstream header;
					header << "%%MatrixMarket hessian coordinate real general" << std::endl;
					header << "% matrix " << name << std::endl;
					header << "% is written by INMOST" << std::endl;
					header << "% by MPI_File_* api" << std::endl;
					header << matsize << " " << matsize << " " << nonzero << std::endl;
					//std::string header_data(header.str());
					ierr = MPI_File_write_shared(fh,const_cast<char *>(header.str().c_str()),static_cast<int>(header.str().size()),MPI_CHAR,&stat);
					if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				}
				ierr = MPI_File_write_ordered(fh,const_cast<char *>(mtx.str().c_str()),static_cast<int>(mtx.str().size()),MPI_CHAR,&stat);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
				ierr = MPI_File_close(&fh);
				if( ierr != MPI_SUCCESS ) MPI_Abort(GetCommunicator(),__LINE__);
			}
#elif defined(USE_MPI)//USE_MPI alternative
			std::string senddata = mtx.str(), recvdata;
			int sendsize = static_cast<int>(senddata.size());
			std::vector<int> recvsize(size), displ(size);
			MPI_Gather(&sendsize,1,MPI_INT,&recvsize[0],1,MPI_INT,0,GetCommunicator());
			if( rank == 0 )
			{
				int totsize = recvsize[0];

				displ[0] = 0;
				for(int i = 1; i < size; i++)
				{
					totsize += recvsize[i];
					displ[i] = displ[i-1]+recvsize[i-1];
				}
				recvdata.resize(totsize);
			}
			else recvdata.resize(1); //protect from dereferencing null
			MPI_Gatherv(&senddata[0],sendsize,MPI_CHAR,&recvdata[0],&recvsize[0],&displ[0],MPI_CHAR,0,GetCommunicator());
			if( rank == 0 )
			{
				std::fstream output(file.c_str(),std::ios::out);
				output << "%%MatrixMarket hessian coordinate real general" << std::endl;
				output << "% matrix " << name << std::endl;
				output << "% is written by INMOST" << std::endl;
				output << "% by MPI_Gather* api and sequential write" << std::endl;
				output << matsize << " " << matsize << " " << nonzero << std::endl;
				output << recvdata;
			}
#else
			std::fstream output(file.c_str(),std::ios::out);
			output << "%%MatrixMarket hessian coordinate real general" << std::endl;
			output << "% matrix " << name << std::endl;
			output << "% is written by INMOST" << std::endl;
			output << "% by sequential write " << std::endl;
			output << matsize << " " << matsize << " " << nonzero << std::endl;
			output << mtx.rdbuf();
#endif
		}

		void HessianMatrix::Swap(HessianMatrix & other)
		{
			data.swap(other.data);
			name.swap(other.name);
			std::swap(comm,other.comm);
			std::swap(is_parallel,other.is_parallel);
		}
		void HessianMatrix::MatVec(INMOST_DATA_REAL_TYPE alpha, const Matrix & U, INMOST_DATA_REAL_TYPE beta, Matrix & J) const //y = alpha*A*x + beta * y
		{
			INMOST_DATA_ENUM_TYPE mbeg, mend;
			INMOST_DATA_INTEGER_TYPE ind, imbeg, imend;
			if( J.Empty() )
			{
				INMOST_DATA_ENUM_TYPE vbeg,vend;
				GetInterval(vbeg,vend);
				J.SetInterval(vbeg,vend);
			}
			//CHECK SOMEHOW FOR DEBUG THAT PROVIDED VECTORS ARE OK
			//~ assert(GetFirstIndex() == out.GetFirstIndex());
			//~ assert(Size() == out.Size());
			GetInterval(mbeg,mend);
			imbeg = mbeg;
			imend = mend;
#if defined(USE_OMP)
#pragma omp for private(ind)
#endif
			for(ind = imbeg; ind < imend; ++ind) //iterate rows of matrix
				(*this)[ind].RowVec(alpha,U[ind],beta,J[ind]);
			// outer procedure should update J Matrix, if needed
		}


		void Sparse::Matrix::ZAXPBY(INMOST_DATA_REAL_TYPE alpha, const Sparse::Matrix& X, INMOST_DATA_REAL_TYPE beta, const Sparse::Matrix& Y, Sparse::Matrix& Z)
		{
			if (X.GetFirstIndex() != Y.GetFirstIndex())
			{
				std::cout << __FILE__ << ":" << __LINE__ << " X and Y first index mismatch!" << std::endl;
				throw Impossible;
			}
			if (X.GetLastIndex() != Y.GetLastIndex())
			{
				std::cout << __FILE__ << ":" << __LINE__ << " X and Y last index mismatch!" << std::endl;
				throw Impossible;
			}
			Z.SetInterval(X.GetFirstIndex(), X.GetLastIndex());
#if defined(USE_OMP)
#pragma omp parallel
#endif
			{
				Sparse::RowMerger merger;
#if defined(USE_OMP)
#pragma omp for
#endif
				for (INMOST_DATA_INTEGER_TYPE it = Z.GetFirstIndex(); it < static_cast<INMOST_DATA_INTEGER_TYPE>(Z.GetLastIndex()); ++it)
				{
					Z[it].Clear();
					merger.AddRow(alpha, X[it]);
					merger.AddRow(beta, Y[it]);
					merger.RetrieveRow(Z[it]);
					merger.Clear();
				}
			}
		}

#endif //USE_SOLVER
	}
}
