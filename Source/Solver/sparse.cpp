#define _CRT_SECURE_NO_WARNINGS
#include "inmost_sparse.h"
#include "inmost_dense.h"
#include <fstream>
#include <sstream>

namespace INMOST
{
	namespace Sparse
	{

#if defined(USE_SOLVER) || defined(USE_AUTODIFF)
		bool _hasRowEntryType = false;

		INMOST_MPI_Type RowEntryType = INMOST_MPI_DATATYPE_NULL;

		INMOST_MPI_Type GetRowEntryType() {return RowEntryType;}

		bool HaveRowEntryType() {return _hasRowEntryType;}

		void CreateRowEntryType()
		{
#if defined(USE_MPI)
			if( !HaveRowEntryType() )
			{
				int ierr;
				MPI_Datatype type[2] = { INMOST_MPI_DATA_ENUM_TYPE, INMOST_MPI_DATA_REAL_TYPE };
				int blocklen[2] = { 1, 1 };
				MPI_Aint disp[2];
				disp[0] = offsetof(Sparse::Row::entry,first);
				disp[1] = offsetof(Sparse::Row::entry,second);
				//disp[2] = sizeof(Sparse::Row::entry);
				ierr = MPI_Type_create_struct(2, blocklen, disp, type, &RowEntryType);
				if( ierr != MPI_SUCCESS )
				{
					std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Type_create_struct" << std::endl;
				}
				ierr = MPI_Type_commit(&RowEntryType);
				if( ierr != MPI_SUCCESS )
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
			if( HaveRowEntryType() )
			{
				MPI_Type_free(&RowEntryType);
				RowEntryType = INMOST_MPI_DATATYPE_NULL;
			}
#endif
			_hasRowEntryType = false;
		}

////////class RowMerger

		RowMerger::RowMerger() : Sorted(true), Nonzeros(0), IntervalBeg(0), IntervalEnd(0) {}

		INMOST_DATA_REAL_TYPE & RowMerger::operator[] (INMOST_DATA_ENUM_TYPE pos)
		{
			INMOST_DATA_ENUM_TYPE map = MapIndex(pos);
			if( map == ENUMUNDEF )
				return Nonlocal[pos];
			if( LinkedList[map+1].first != UNDEF ) return LinkedList[map+1].second;
			else
			{
				INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg(), next;
				if( Sorted )
				{
					next = index;
					while(next < map+1)
					{
						index = next;
						next = LinkedList[index].first;
					}
					assert(index < map+1);
					assert(map+1 < next);
					++Nonzeros;
					LinkedList[index].first = map+1;
					LinkedList[map+1].first = next;
					return LinkedList[map+1].second;
				}
				else
				{
					INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg();
					++Nonzeros;
					LinkedList[map+1].first = LinkedList[index].first;
					LinkedList[index].first = map+1;
					return LinkedList[map+1].second;
				}
			}
		}

		INMOST_DATA_REAL_TYPE RowMerger::operator[] (INMOST_DATA_ENUM_TYPE pos) const
		{
			INMOST_DATA_ENUM_TYPE map = MapIndex(pos);
			if( map == ENUMUNDEF )
			{
				std::map<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE>::const_iterator it = Nonlocal.find(pos);
				if( it != Nonlocal.end() )
					return it->second;
			}
			else
				if( LinkedList[map+1].first != UNDEF ) return LinkedList[map+1].second;
			throw -1;
		}

		RowMerger::RowMerger(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool Sorted)
				: Sorted(Sorted), Nonzeros(0), IntervalBeg(interval_begin), IntervalEnd(interval_end), LinkedList(interval_begin,interval_end+1,Row::make_entry(UNDEF,0.0))
		{
			LinkedList.begin()->first = EOL;
		}

		void RowMerger::Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end, bool _Sorted)
		{
			LinkedList.set_interval_beg(static_cast<INMOST_DATA_ENUM_TYPE>(interval_begin));// - NonlocalPre.size()));
			LinkedList.set_interval_end(static_cast<INMOST_DATA_ENUM_TYPE>(interval_end + 1));// + NonlocalPost.size()));
			IntervalBeg = interval_begin;
			IntervalEnd = interval_end;
			std::fill(LinkedList.begin(),LinkedList.end(),Row::make_entry(UNDEF,0.0));
			LinkedList.begin()->first = EOL;
			Nonzeros = 0;
			Sorted = _Sorted;
			Nonlocal.clear();
		}

#if defined(USE_SOLVER)
		void RowMerger::Resize(const Matrix & A, bool _Sorted)
		{
			INMOST_DATA_ENUM_TYPE mbeg, mend, k, l, ind;
			//retrieve interval of indices
			A.GetInterval(mbeg,mend);
			//gather non-local mapping from matrix
			std::set<INMOST_DATA_ENUM_TYPE> Pre, Post;
			for(k = mbeg; k < mend; ++k)
			{
				for(l = 0; l < A[k].Size(); ++l)
				{
					ind = A[k].GetIndex(l);
					if( ind < mbeg ) Pre.insert(ind);
					else if( ind >= mend ) Post.insert(ind);
				}
			}
			Resize(mbeg,mend,_Sorted);
			Nonlocal.clear();
		}

		RowMerger::RowMerger(const Matrix & A, bool Sorted) : Sorted(Sorted), Nonzeros(0)
		{
			Resize(A,Sorted);
		}
#endif

		RowMerger::~RowMerger() {}

		/*
		INMOST_DATA_ENUM_TYPE RowMerger::MapIndex(INMOST_DATA_ENUM_TYPE pos) const
		{
			if( pos < IntervalBeg || pos >= IntervalEnd )
				return ENUMUNDEF;
			return pos;
		}

		INMOST_DATA_ENUM_TYPE RowMerger::UnmapIndex(INMOST_DATA_ENUM_TYPE pos) const
		{
			return pos;
		}
		*/


		void RowMerger::Clear()
		{
			INMOST_DATA_ENUM_TYPE i = LinkedList.begin()->first, j;
			LinkedList.begin()->first = EOL;
			while( i != EOL )
			{
				j = LinkedList[i].first;
				LinkedList[i].first = UNDEF;
				LinkedList[i].second = 0.0;
				i = j;
			}
			Nonlocal.clear();
			Nonzeros = 0;
		}


		void RowMerger::PushRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow)
		{
			if( Sorted && PreSortRow ) std::sort(r.Begin(),r.End());
			PushRow(coef,r);
		}
		
		void RowMerger::Multiply(INMOST_DATA_REAL_TYPE coef)
		{
			INMOST_DATA_ENUM_TYPE i = LinkedList.begin()->first;
			while( i != EOL )
			{
				LinkedList[i].second *= coef;
				i = LinkedList[i].first;
			}
			if( !Nonlocal.empty() )
			{
				std::map< INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE >::iterator it;
				for(it = Nonlocal.begin(); it != Nonlocal.end(); ++it)
					it->second *= coef;
			}
		}
		

		void RowMerger::PushRow(INMOST_DATA_REAL_TYPE coef, const Row & r)
		{
			assert(Nonzeros == 0); //Linked list should be empty
			assert(LinkedList.begin()->first == EOL); //again check that list is empty
			if( coef )
			{
				INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg();
				Row::const_iterator it = r.Begin(), jt;
				while( it != r.End() )
				{
					INMOST_DATA_ENUM_TYPE pos = it->first;
					INMOST_DATA_ENUM_TYPE map = MapIndex(pos);
					if( map == ENUMUNDEF )
						Nonlocal[pos] = it->second*coef;
					else
					{
						LinkedList[index].first = map+1;
						LinkedList[map+1].first = EOL;
						LinkedList[map+1].second = it->second*coef;
						index = map+1;
						++Nonzeros;
					}
					jt = it;
					++it;
					assert(!Sorted || it == r.End() || jt->first < it->first);
				}
			}
		}

		void RowMerger::AddRow(INMOST_DATA_REAL_TYPE coef, Row & r, bool PreSortRow)
		{
			if( Sorted && PreSortRow ) std::sort(r.Begin(),r.End());
			AddRow(coef,r);
		}

		void RowMerger::AddRow(INMOST_DATA_REAL_TYPE coef, const Row & r)
		{
			if( coef )
			{
				INMOST_DATA_ENUM_TYPE index = LinkedList.get_interval_beg(), next;
				Row::const_iterator it = r.Begin(), itend = r.End(), jt;
				while( it != itend )
				{
					INMOST_DATA_ENUM_TYPE pos = it->first;
					INMOST_DATA_ENUM_TYPE map = MapIndex(pos);
					if( map == ENUMUNDEF )
						Nonlocal[pos] += coef*it->second;
					else
					{
						if( LinkedList[map+1].first != UNDEF )
							LinkedList[map+1].second += coef*it->second;
						else if( Sorted )
						{
							next = index;
							while(next < map+1)
							{
								index = next;
								next = LinkedList[index].first;
							}
							assert(index < map+1);
							assert(map+1 < next);
							LinkedList[index].first = map+1;
							LinkedList[map+1].first = next;
							LinkedList[map+1].second = coef*it->second;
							++Nonzeros;
						}
						else
						{
							LinkedList[map+1].first = LinkedList[index].first;
							LinkedList[map+1].second = coef*it->second;
							LinkedList[index].first = map+1;
							++Nonzeros;
						}
					}
					jt = it;
					++it;
					assert(!Sorted || it == itend || jt->first < it->first);
				}
			}
		}


		void RowMerger::RetrieveRow(Row & r)
		{
            r.Resize(static_cast<INMOST_DATA_ENUM_TYPE>(Nonzeros+Nonlocal.size()));
			INMOST_DATA_ENUM_TYPE i = LinkedList.begin()->first, k = 0;
			while( i != EOL )
			{
				if( LinkedList[i].second )
				{
					r.GetIndex(k) = UnmapIndex(i-1);
					r.GetValue(k) = LinkedList[i].second;
					++k;
				}
				i = LinkedList[i].first;
			}
			if( !Nonlocal.empty() )
			{
				std::map< INMOST_DATA_ENUM_TYPE, INMOST_DATA_REAL_TYPE >::iterator it;
				for(it = Nonlocal.begin(); it != Nonlocal.end(); ++it) if( it->second )
				{
					r.GetIndex(k) = it->first;
					r.GetValue(k) = it->second;
					++k;
				}
			}
			r.Resize(k);
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
					ret(i,j) = data[rows(i,j)];
			return ret;
		}
#if defined(USE_AUTODIFF)
		INMOST::Matrix<value_reference> Vector::operator [](const INMOST::AbstractMatrix<INMOST_DATA_INTEGER_TYPE>& rows)
		{
			INMOST::Matrix<value_reference> ret(rows.Rows(), rows.Cols());
			for (INMOST_DATA_ENUM_TYPE i = 0; i < rows.Rows(); ++i)
				for (INMOST_DATA_ENUM_TYPE j = 0; j < rows.Cols(); ++j)
					new (&ret(i, j)) value_reference(data[rows(i, j)]);
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


#endif //USE_SOLVER
	}
}
