#define _CRT_SECURE_NO_WARNINGS
#include "inmost_sparse.h"
#include "inmost_dense.h"
#include <fstream>
#include <sstream>
#include <queue>
#if defined(USE_ZLIB)
#include "zlib.h"
#endif

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
		INMOST_DATA_REAL_TYPE   Row::RowVec(Vector & x) const
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

		static bool compress(const void* buffer, size_t dsize, void*& buffer_out, size_t& dsize_out)
		{
#if defined(USE_ZLIB)
			uLongf bufsize = compressBound((uLongf)dsize);
			void* zbuffer = malloc(bufsize);
			if (zbuffer)
			{
				int res = compress2(static_cast<Bytef*>(zbuffer), &bufsize, static_cast<const Bytef*>(buffer), (uLongf)dsize, 9);
				if (res != Z_OK)
					std::cout << __FILE__ << ":" << __LINE__ << " fail " << res << std::endl;
				else
				{
					buffer_out = zbuffer;
					dsize_out = bufsize;
					return true;
				}
			}
			else std::cout << __FILE__ << ":" << __LINE__ << " allocation of " << bufsize << " failed" << std::endl;
			return false;
#else
			return false;
#endif
		}


		void     Vector::SaveBinary(std::string file)
		{
			int rank = 0, size = 1, compr = 0;
			int vecsize = Size(), totsize = 0;
			std::vector<int> vecsizes(size, vecsize);
#if defined(USE_MPI)
			MPI_Comm_rank(GetCommunicator(), &rank);
			MPI_Comm_size(GetCommunicator(), &size);
			if (rank == 0)
				vecsizes.resize(size);
			MPI_Gather(&vecsize, 1, MPI_INT, &vecsizes[0], 1, MPI_INT, 0, GetCommunicator());
#endif
			size_t dsize = vecsize * sizeof(INMOST_DATA_REAL_TYPE), zdsize;
			void* buffer = static_cast<void*>(&data[0]), * zbuffer;
			if (compress(buffer, dsize, zbuffer, zdsize))
			{
				std::swap(buffer, zbuffer);
				std::swap(dsize, zdsize);
				compr = 1;
			}
			if (dsize > INT_MAX)
			{
				std::cout << __FILE__ << ":" << __LINE__;
				std::cout << " Current implementation cannot fit necessery size ";
				std::cout << dsize << ">" << INT_MAX << " on rank " << rank << std::endl;
				MPI_Abort(GetCommunicator(), __LINE__);
			}
			std::vector<int> dsizes(size, dsize), displ(size, 0);
#if defined(USE_MPI)
			MPI_Gather(&dsize, 1, MPI_INT, &dsizes[0], 1, MPI_INT, 0, GetCommunicator());
#endif
			size_t totbufsize = 1;
			if (rank == 0)
			{
				totbufsize = 0;
				for (int k = 0; k < size; ++k)
				{
					displ[k] = totbufsize;
					totbufsize += dsizes[k];
				}
			}
			std::vector<char> totbuffer(totbufsize);
#if defined(USE_MPI)
			MPI_Gatherv(buffer, dsize, MPI_CHAR, &totbuffer[0], &dsizes[0], &displ[0], MPI_INT, 0, GetCommunicator());
#else
			std::copy(buffer, buffer + dsize, &totbuffer[0]);
#endif
			if (compr) free(buffer); //temporary buffer
			if (rank == 0)
			{
				std::ofstream fout(file.c_str(), std::ios::binary);
				fout.write(reinterpret_cast<const char*>(&compr), sizeof(int));
				fout.write(reinterpret_cast<const char*>(&size), sizeof(int));
				fout.write(reinterpret_cast<const char*>(&vecsizes[0]), sizeof(int) * size);
				fout.write(reinterpret_cast<const char*>(&dsizes[0]), sizeof(int) * size);
				fout.write(reinterpret_cast<const char*>(&totbuffer[0]), sizeof(char) * totbufsize);
			}
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
		void Matrix::MatVec(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & out) const //y = alpha*A*x + beta * y
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
		
		
		


		void Matrix::MatVecTranspose(INMOST_DATA_REAL_TYPE alpha, Vector & x, INMOST_DATA_REAL_TYPE beta, Vector & out) const //y = alpha*A*x + beta * y
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
			int rank = 0, size = 1, compr[3] = { 0,0,0 };
			int matsize = GetLastIndex() - GetFirstIndex(), totsize = 0, nnzsize = Nonzeros();
			std::vector<int> matsizes(size, matsize), nnzsizes(size, nnzsize);
			std::vector<INMOST_DATA_ENUM_TYPE> ia;
			std::vector<INMOST_DATA_ENUM_TYPE> ja;
			std::vector<INMOST_DATA_ENUM_TYPE> va;
			ia.reserve(matsize);
			ja.reserve(nnzsize);
			va.reserve(nnzsize);
			ia.push_back(0);
			for (INMOST_DATA_ENUM_TYPE k = GetFirstIndex(); k < GetLastIndex(); ++k)
			{
				for (Row::iterator jt = (*this)[k].Begin(); jt != (*this)[k].End(); ++jt)
				{
					ja.push_back(jt->first);
					va.push_back(jt->second);
				}
				ia.push_back(ja.size());
			}
#if defined(USE_MPI)
			MPI_Comm_rank(GetCommunicator(), &rank);
			MPI_Comm_size(GetCommunicator(), &size);
			if (rank == 0)
				matsizes.resize(size);
			MPI_Gather(&matsize, 1, MPI_INT, &matsizes[0], 1, MPI_INT, 0, GetCommunicator());
			MPI_Gather(&nnzsize, 1, MPI_INT, &nnzsizes[0], 1, MPI_INT, 0, GetCommunicator());
#endif
			size_t dsize[3] =
			{
				matsize * sizeof(INMOST_DATA_ENUM_TYPE),
				nnzsize * sizeof(INMOST_DATA_ENUM_TYPE),
				nnzsize * sizeof(INMOST_DATA_REAL_TYPE)
			}, zdsize;
			void* buffer[3] =
			{
				static_cast<void*>(&ia[0]),
				static_cast<void*>(&ja[0]),
				static_cast<void*>(&va[0])
			}, * zbuffer;
			for (int k = 0; k < 3; ++k)
			{
				if (compress(buffer[k], dsize[k], zbuffer, zdsize))
				{
					std::swap(buffer[k], zbuffer);
					std::swap(dsize[k], zdsize);
					compr[0] = 1;
				}
				if (dsize[k] > INT_MAX)
				{
					std::cout << __FILE__ << ":" << __LINE__;
					std::cout << " Current implementation cannot fit necessery size ";
					std::cout << dsize[k] << ">" << INT_MAX << " k " << k << " on rank " << rank << std::endl;
					MPI_Abort(GetCommunicator(), __LINE__);
				}
			}
			if (compr[0])
			{
				ia.clear();
				ia.shrink_to_fit();
			}
			if (compr[1])
			{
				ja.clear();
				ja.shrink_to_fit();
			}
			if (compr[2])
			{
				va.clear();
				va.shrink_to_fit();
			}
			std::vector<int> dsizes[3], displ(size, 0);
			std::vector<char> totbuffer[3];
			for (int k = 0; k < 3; ++k)
			{
				dsizes[k].resize(size);
				dsizes[k][0] = dsize[k];
#if defined(USE_MPI)
				MPI_Gather(&dsize[k], 1, MPI_INT, &dsizes[k][0], 1, MPI_INT, 0, GetCommunicator());
#endif
				size_t totbufsize = 1;
				if (rank == 0)
				{
					totbufsize = 0;
					for (int q = 0; q < size; ++q)
					{
						displ[q] = totbufsize;
						totbufsize += dsizes[k][q];
					}
				}
				if (rank == 0)
					totbuffer[k].resize(totbufsize);
				else totbuffer[k].resize(1);
#if defined(USE_MPI)
				MPI_Gatherv(buffer[k], dsize[k], MPI_CHAR, &totbuffer[k][0], &dsizes[k][0], &displ[0], MPI_INT, 0, GetCommunicator());
#else
				std::copy(buffer[k], buffer[k] + dsize[k], &totbuffer[k][0]);
#endif
				if (compr[k]) free(buffer[k]); //temporary buffer
			}
			if (rank == 0)
			{
				std::ofstream fout(file.c_str(), std::ios::binary);
				fout.write(reinterpret_cast<const char*>(&compr), sizeof(int));
				fout.write(reinterpret_cast<const char*>(&size), sizeof(int));
				fout.write(reinterpret_cast<const char*>(&matsizes[0]), sizeof(int) * size);
				fout.write(reinterpret_cast<const char*>(&nnzsizes[0]), sizeof(int) * size);
				for (int k = 0; k < 3; ++k)
					fout.write(reinterpret_cast<const char*>(&dsizes[k][0]), sizeof(int) * size);
				for (int k = 0; k < 3; ++k)
					fout.write(reinterpret_cast<const char*>(&totbuffer[k][0]), sizeof(char) * totbuffer[k].size());
			}
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
