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

		inline static  unsigned int flip(const unsigned int* fp)
		{
			unsigned int mask = -((int)(*fp >> 31)) | 0x80000000;
			return *fp ^ mask;
		}

		inline void RowMerger2::radix_sort()
		{
			temp.resize(inds.size());
			unsigned int i, size = static_cast<unsigned int>(inds.size());
			const unsigned int kHist = 2048;
			unsigned int  b0[kHist * 3];
			unsigned int* b1 = b0 + kHist;
			unsigned int* b2 = b1 + kHist;
			memset(b0, 0, sizeof(unsigned int) * kHist * 3);
			for (i = 0; i < size; i++)
			{
				unsigned int fi = flip((unsigned int*)&inds[i]);
				++b0[_m0(fi)]; ++b1[_m1(fi)]; ++b2[_m2(fi)];
			}
			{
				unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
				for (i = 0; i < kHist; i++)
				{
					b0[kHist - 1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist - 1];
					b1[kHist - 1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist - 1];
					b2[kHist - 1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist - 1];
				}
			}
			for (i = 0; i < size; i++) temp[++b0[_m0(flip((unsigned int*)&inds[i]))]] = inds[i];
			for (i = 0; i < size; i++) inds[++b1[_m1(flip((unsigned int*)&temp[i]))]] = temp[i];
			for (i = 0; i < size; i++) temp[++b2[_m2(flip((unsigned int*)&inds[i]))]] = inds[i];
			for (i = 0; i < size; i++) inds[i] = temp[i];
		}

		inline void RowMerger2::naive_sort() 
		{
			size_t i, j;
			for (i = 1; i < inds.size(); i++) 
			{
				INMOST_DATA_ENUM_TYPE tmp = inds[i];
				for (j = i; j >= 1 && tmp < inds[j - 1]; j--)
					inds[j] = inds[j - 1];
				inds[j] = tmp;
			}
		}

		void RowMerger2::clear()
		{
			//indset.clear();
			//bitset.clear();
			vals.clear();
			inds.clear();
		}
		void RowMerger2::set_vals()
		{
			//radix_sort();
			std::sort(inds.begin(), inds.end());
			//naive_sort();
			//inds.resize(indset.size());
			//std::copy(indset.begin(), indset.end(), inds.begin());
			vals.resize(inds.size(), 0.0);
		}
		void RowMerger2::set_bitset(INMOST_DATA_ENUM_TYPE beg, INMOST_DATA_ENUM_TYPE end)
		{
			if (end - beg > bitset.size()) bitset.resize(end - beg);
			std::fill(bitset.begin(), bitset.begin() + (end - beg), 0);
		}
		void RowMerger2::get_row(Sparse::Row& r)
		{
			INMOST_DATA_ENUM_TYPE k = 0, s = static_cast<INMOST_DATA_ENUM_TYPE>(inds.size());
			r.Resize(s);
			for (INMOST_DATA_ENUM_TYPE i = 0; i < r.Size(); ++i)
			{
				if (1.0 + vals[i] != 1.0)
				{
					r.GetIndex(k) = inds[i];
					r.GetValue(k) = vals[i];
					k++;
				}
			}
			r.Resize(k);
		}

		void RowMerger3::clear()
		{
			inds.clear();
			vals.clear();
		}

		void RowMerger3::set_vals()
		{
			vals.resize(inds.size());
			std::fill(vals.begin(), vals.end(), 0);
		}

		/*
		inline void RowMerger3::merge_indices(std::vector<INMOST_DATA_ENUM_TYPE>& inds, const Sparse::Row& r, std::vector<INMOST_DATA_ENUM_TYPE>& temp)
		{
			std::vector<INMOST_DATA_ENUM_TYPE>::iterator first1 = inds.begin(), last1 = inds.begin();
			Sparse:Row::const_iterator first2 = r.Begin(), last2 = r.End();
			temp.clear();
			while(first1 != last1)
			{
				if (first2 == last2)
				{
					temp.insert(temp.end(), first1, last1);
					break;
				}
					
				if (first2->first < *first1)
				{
					temp.push_back(first2->first);
					first2++;
				}
				else
				{
					temp.push_back(*first1);
					if (!(*first1 < first2->first))
						++first2;
					++first1;
				}
			}
			temp.insert(temp.end(), first2, last2);
			inds.swap(temp);
		}

		inline void RowMerger3::insert_index(std::vector<INMOST_DATA_ENUM_TYPE>& inds, INMOST_DATA_ENUM_TYPE ind)
		{
			std::vector<INMOST_DATA_ENUM_TYPE>::iterator it = std::lower_bound(inds.begin(), inds.end(), ind);
			if (*it != ind) inds.insert(it, ind);
		}

		inline void RowMerger3::set_indices(std::vector<INMOST_DATA_ENUM_TYPE>& inds, const Sparse::Row& r)
		{
			inds.resize(r.Size());
			for (unsigned k = 0; k < r.Size(); ++k)
				inds[k] = r.GetIndex(k);
		}
		*/


		void RowMerger3::get_row(Sparse::Row& r)
		{
			INMOST_DATA_ENUM_TYPE k = 0, s = static_cast<INMOST_DATA_ENUM_TYPE>(inds.size());
			r.Resize(s);
			for (INMOST_DATA_ENUM_TYPE i = 0; i < r.Size(); ++i)
			{
				if (1.0 + vals[i] != 1.0)
				{
					r.GetIndex(k) = inds[i];
					r.GetValue(k) = vals[i];
					k++;
				}
			}
			r.Resize(k);
		}

		void RowMerger4::clear()
		{
			inds.Clear();
		}

		void RowMerger4::get_row(Sparse::Row& r)
		{
			r.Swap(inds);
		}

		void RowMerger5::clear()
		{
#if defined(TEST_HASHTABLE)
			table.Clear();
#else
			table.clear();
#endif
			vals.clear();
		}

		void RowMerger5::add_value(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_REAL_TYPE val)
		{
#if defined(TEST_HASHTABLE)
			HashTable::iterator it = table.find(ind);
			if (it == table.end())
			{
				table.Insert(ind)->second = vals.size();
				vals.push_back(val);
			}
			else vals[it->second] += val;
#elif defined(TEST_JUDY1)
			uint64_t test = table.find(ind + 1);
			if (!table.success())
			{
				if (val)
				{
					table.insert(ind + 1, vals.size() + 1);
					vals.push_back(val);
				}
			}
			else vals[test - 1] += val;
#else
			if (1.0 + val != 1.0)
			{
				/*
				INMOST_DATA_ENUM_TYPE& key = table[ind];
				if (key == 0)
				{
					key = vals.size() + 1;
					vals.push_back(val);
				}
				else vals[key - 1] += val;
				*/
				table_t::iterator test = table.find(ind);
				if (test == table.end())
				{
					table[ind] = vals.size();
					vals.push_back(val);
				}
				else vals[test->second] += val;
				/*
				std::pair<table_t::iterator, bool> test
					= table.insert(std::make_pair(ind, static_cast<INMOST_DATA_ENUM_TYPE>(vals.size())));
				if (test.second)
					vals.push_back(val);
				else vals[test.first->second] += val;
				*/
			}
#endif
		}

		void RowMerger5::get_row(Sparse::Row& r)
		{
			INMOST_DATA_ENUM_TYPE k = 0, s = static_cast<INMOST_DATA_ENUM_TYPE>(vals.size());
			r.Resize(s);
#if defined(TEST_HASHTABLE)
			for (HashTable::iterator it = table.begin(); it != table.end(); ++it)
				if (1.0 + vals[it->second] != 1.0)
				{
					r.GetIndex(k) = it->first;
					r.GetValue(k) = vals[it->second];
					k++;
				}
#elif defined(TEST_JUDY1)
			const judyLArray< uint64_t, uint64_t>::pair * it = &table.begin();
			while (table.success())
			{
				if (1.0 + vals[it->value - 1] != 1.0)
				{
					r.GetIndex(k) = it->key - 1;
					r.GetValue(k) = vals[it->value - 1];
					k++;
				}
				it = &table.next();
			}
#else
			for(table_t::iterator it = table.begin(); it != table.end(); ++it)
				if (1.0 + vals[it->second] != 1.0)
				{
					r.GetIndex(k) = it->first;
					r.GetValue(k) = vals[it->second];
					k++;
				}
#endif
			r.Resize(k);
		}

		////////class RowMerger
				/*

				RowMerger::RowMerger() {}
				INMOST_DATA_REAL_TYPE RowMerger::operator[] (INMOST_DATA_ENUM_TYPE ind) const
				{
					container::const_iterator f = data.find(ind);
					if (f != data.end())
						return f->second;
					throw - 1;
				}


				//void RowMerger::Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end)
				void RowMerger::Resize(INMOST_DATA_ENUM_TYPE size)
				{
					assert(data.empty());
					data.reserve(size);
				}

				void RowMerger::Multiply(INMOST_DATA_REAL_TYPE coef)
				{
					for (container::iterator it = data.begin(); it != data.end(); ++it)
						it->second *= coef;
				}

				void RowMerger::PushRow(INMOST_DATA_REAL_TYPE coef, const Row& r)
				{
					assert(data.empty()); //Linked list should be empty
					if (coef)
					{
						for (Row::const_iterator it = r.Begin(); it != r.End(); ++it)
							data[it->first] = it->second * coef;
					}
				}

				void RowMerger::AddRow(INMOST_DATA_REAL_TYPE coef, const Row& r)
				{
					if (coef)
					{
						for (Row::const_iterator it = r.Begin(); it != r.End(); ++it)
							data[it->first] += coef * it->second;
					}
				}

				void RowMerger::RetriveRow(Row& r)
				{
					r.Resize(Size());
					INMOST_DATA_ENUM_TYPE k = 0;
					for (container::iterator it = data.begin(); it != data.end(); ++it)
					{
						if (it->second)
						{
							r.GetIndex(k) = it->first;
							r.GetValue(k) = it->second;
							k++;
						}
					}
					r.Resize(k);
				}
				*/

#if 0
		INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind) const
		{
			judyvalue key = ind;
			JudySlot* slot = (JudySlot*)judy_slot((Judy*)judy_array, (const unsigned char*)&key, JUDY_key_size);
			if (slot != NULL)
			{
				assert((*slot) - 1 < vals.size());
				return (*slot) - 1;
			}
			else throw Impossible;
		}
		INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind)
		{
			judyvalue key = ind;
			JudySlot* slot = (JudySlot*)judy_slot((Judy*)judy_array, (const unsigned char*)&key, JUDY_key_size);
			if (slot != NULL)
			{
				assert((*slot) - 1 < vals.size());
				return (*slot) - 1;
			}
			else
			{
				ins_pos(ind, Nonzeros);
				next.push_back(First);
				vals.push_back(0.0);
				First = ind;
				return Nonzeros++;
			}
		}
		void RowMerger::ins_pos(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE k)
		{
			judyvalue key = ind;
			JudySlot* slot = (JudySlot*)judy_cell((Judy*)judy_array, (const unsigned char*)&key, JUDY_key_size);
			assert(slot != NULL);
			*slot = static_cast<judyvalue>(k)+1;
		}
		RowMerger::RowMerger() : First(EOL), Nonzeros(0) 
		{
			judy_array = judy_open(JUDY_key_size, 1);
		}
		void RowMerger::Clear()
		{
			First = EOL;
			judyvalue key = 0;
			JudySlot* slot = NULL;
			while( NULL != (slot = judy_strt((Judy*)judy_array, (const unsigned char*)&key, 0)))
				judy_del((Judy*)judy_array);
			vals.clear();
			next.clear();
			Nonzeros = 0;
		}

#elif 0
		INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind) const
		{
			map_container::const_iterator f = pos.find(ind);
			if (f != pos.end())
				return f->second;
			else throw Impossible;
		}
		INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind)
		{
			map_container::iterator f = pos.find(ind);
			if (f != pos.end())
				return f->second;
			else
			{
				ins_pos(ind, Nonzeros);
				next.push_back(First);
				vals.push_back(0.0);
				First = ind;
				return Nonzeros++;
			}
		}
		void RowMerger::ins_pos(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE k)
		{
			pos[ind] = k;
		}
		RowMerger::RowMerger() : First(EOL), Nonzeros(0) {}
		void RowMerger::Clear()
		{
			First = EOL;
			pos.clear();
			vals.clear();
			next.clear();
			Nonzeros = 0;
		}
#elif defined(TEST_HASHTABLE)
	INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind) const
	{
		const map_container::Cell * f = pos.Lookup(ind);
		if (f != NULL)
			return f->second;
		else throw Impossible;
	}
	void RowMerger::ins_pos(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE val)
	{
		pos.Insert(ind)->second = val;
	}
	INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind)
	{
		map_container::Cell* f = pos.Lookup(ind);
		if (f != NULL)
			return f->second;
		else
		{
			ins_pos(ind, Nonzeros);
			vals.push_back(0.0);
			return Nonzeros++;
		}
	}
	RowMerger::RowMerger() : Nonzeros(0)
	{
		vals.reserve(1024);
		//pos.reserve(2048);
		//pos.max_load_factor(0.25);
	}
	void RowMerger::Clear()
	{
		pos.Clear();
		vals.clear();
		Nonzeros = 0;
	}
#else
	INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind) const
	{
		map_container::const_iterator f = pos.find(ind);
		if (f != pos.end())
			return f->second;
		else throw Impossible;
	}
	INMOST_DATA_ENUM_TYPE RowMerger::get_pos(INMOST_DATA_ENUM_TYPE ind)
	{
		map_container::iterator f = pos.find(ind);
		if (f != pos.end())
			return f->second;
		else
		{
			ins_pos(ind, Nonzeros);
			vals.push_back(0.0);
			return Nonzeros++;
		}
	}
	void RowMerger::ins_pos(INMOST_DATA_ENUM_TYPE ind, INMOST_DATA_ENUM_TYPE val)
	{
		pos[ind] = val;
	}
	RowMerger::RowMerger() : Nonzeros(0) 
	{
		vals.reserve(1024);
		pos.reserve(2048);
		//pos.max_load_factor(0.25);
	}
	void RowMerger::Clear()
	{
		pos.clear();
		vals.clear();
		Nonzeros = 0;
	}
#endif

		void RowMerger::AddRow(INMOST_DATA_REAL_TYPE coef, const Row& r)
		{
			if (coef)
			{
				for (Row::const_iterator it = r.Begin(); it != r.End(); ++it)
					vals[get_pos(it->first)] += coef * it->second;
			}
		}
		//void RowMerger::Resize(INMOST_DATA_ENUM_TYPE interval_begin, INMOST_DATA_ENUM_TYPE interval_end)
		void RowMerger::Resize(INMOST_DATA_ENUM_TYPE size)
		{
			assert(Nonzeros == 0);
			Nonzeros = 0;
			vals.clear();
		}

#if defined(USE_SOLVER)
		void RowMerger::Resize(const Matrix & A)
		{
			//INMOST_DATA_ENUM_TYPE mbeg, mend;
			//A.GetInterval(mbeg,mend);
			//Resize(mbeg,mend);
			Resize(A.Size());
		}

		RowMerger::RowMerger(const Matrix & A) : Nonzeros(0)
		{
			Resize(A);
		}
#endif

		RowMerger::~RowMerger() {}

	

		

		void RowMerger::Multiply(INMOST_DATA_REAL_TYPE coef)
		{
			for (size_t k = 0; k < vals.size(); ++k)
				vals[k] *= coef;
		}

		void RowMerger::PushRow(INMOST_DATA_REAL_TYPE coef, const Row & r)
		{
			assert(Nonzeros == 0); //Linked list should be empty
			if( coef )
			{
				vals.resize(r.Size());
				Nonzeros = r.Size();
				INMOST_DATA_ENUM_TYPE k = 0;
				for(Row::const_iterator it = r.Begin(); it != r.End(); ++it)
				{
					ins_pos(it->first, k);
					vals[k] = it->second * coef;
					++k;
				}
			}
		}

		
		


		void RowMerger::RetrieveRow(Row & r) const
		{
			INMOST_DATA_ENUM_TYPE k = 0;
            r.Resize(static_cast<INMOST_DATA_ENUM_TYPE>(Nonzeros));
			for (map_container::const_iterator it = pos.begin(); it != pos.end(); ++it)
			{
				if (1.0 + vals[it->second] != 1.0)
				{
					r.GetIndex(k) = it->first;
					r.GetValue(k) = vals[it->second];
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
		INMOST_DATA_REAL_TYPE   Row::RowVec(Vector & x) const
		{
			INMOST_DATA_REAL_TYPE ret = 0;
			INMOST_DATA_ENUM_TYPE end = Size();
			for(INMOST_DATA_ENUM_TYPE i = 0; i < end; i++) ret = ret + x[GetIndex(i)]*GetValue(i);
			return ret;
		}
#endif //USE_SOLVER
		
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

		void Row::GetIndices(INMOST_DATA_ENUM_TYPE shift, std::vector<Sparse::bit_type>& bitset, std::vector<INMOST_DATA_ENUM_TYPE>& inds) const
		{
			INMOST_DATA_ENUM_TYPE ind, sind;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k) 
			{
				ind = GetIndex(k);
				sind = ind - shift;
				if (!bitset[sind])
				{
					bitset[sind] = true;
					inds.push_back(ind);
				}
			}
		}

		void Row::GetIndices(std::set<INMOST_DATA_ENUM_TYPE>& indset) const
		{
			std::set<INMOST_DATA_ENUM_TYPE>::const_iterator hint = indset.cend();
			for (INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k)
				hint = indset.insert(hint, GetIndex(k));
		}

		void Row::GetPairs(INMOST_DATA_REAL_TYPE coef, Sparse::Row& inds, Sparse::Row& temp) const
		{
			INMOST_DATA_REAL_TYPE val;
			const_iterator first1 = inds.Begin(), last1 = inds.End();
			const_iterator first2 = Begin(), last2 = End();
			temp.Clear();
			temp.Reserve(inds.Size() + Size());
			while (first1 != last1 && first2 != last2)
			{
				if (first2->first < first1->first)
				{
					val = first2->second * coef;
					if(1. + val != 1.) temp.Push(first2->first, val);
					first2++;
				}
				else
				{
					if (!(first1->first < first2->first))
					{
						val = first2->second * coef + first1->second;
						if(1. + val != 1.)
							temp.Push(first1->first, val);
						++first2;
					}
					else if( 1. + first1->second != 1.)
						temp.Push(first1->first, first1->second);
					++first1;
				}
			}
			for (const_iterator it = first1; it < last1; ++it)
				temp.Push(it->first, it->second);
			for (const_iterator it = first2; it < last2; ++it)
				temp.Push(it->first, it->second * coef);
			inds.Swap(temp);
		}

		void Row::GetIndices(std::vector<INMOST_DATA_ENUM_TYPE>& inds, std::vector<INMOST_DATA_ENUM_TYPE>& temp) const
		{
			std::vector<INMOST_DATA_ENUM_TYPE>::const_iterator first1 = inds.begin(), last1 = inds.end();
			const_iterator first2 = Begin(), last2 = End();
			temp.clear();
			temp.reserve(inds.size() + Size());
			while (first1 != last1 && first2 != last2)
			{
				if (first2->first < *first1)
				{
					temp.push_back(first2->first);
					first2++;
				}
				else
				{
					temp.push_back(*first1);
					if (!(*first1 < first2->first))
						++first2;
					++first1;
				}
			}
			if (first1 != last1)
				temp.insert(temp.end(), first1, last1);
			for (const_iterator it = first2; it < last2; ++it)
				temp.push_back(it->first);
			//temp.insert(temp.end(), first2, last2);
			inds.swap(temp);
		}

		void Row::GetIndices(std::vector<INMOST_DATA_ENUM_TYPE>& inds) const
		{
			assert(inds.empty());
			inds.resize(Size());
			for (unsigned k = 0; k < Size(); ++k)
				inds[k] = GetIndex(k);
		}
		
		void Row::GetValues(INMOST_DATA_REAL_TYPE coef, const std::vector<INMOST_DATA_ENUM_TYPE>& inds, std::vector<INMOST_DATA_REAL_TYPE>& vals) const
		{
#if defined(ASSUME_SORTED)
			INMOST_DATA_ENUM_TYPE pos, ind;
			//const INMOST_DATA_ENUM_TYPE *beg = &inds[0], *beg0 = beg;
			std::vector<INMOST_DATA_ENUM_TYPE>::const_iterator beg = inds.begin();
			for (INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k)
			{
				ind = GetIndex(k);
				//beg = std::lower_bound(beg, inds.end(), ind);
				while (*beg != ind) beg++;
				pos = static_cast<INMOST_DATA_ENUM_TYPE>(beg - inds.cbegin());
				vals[pos] += GetValue(k) * coef;
			}
#else
			for (Sparse::Row::const_iterator it = Begin(); it != End(); ++it)
				vals[std::lower_bound(inds.begin(), inds.end(), it->first) - inds.begin()] += it->second * coef;
#endif
		}

		void Row::GetValues(INMOST_DATA_REAL_TYPE coef, const std::set<INMOST_DATA_ENUM_TYPE>& indset, std::vector<INMOST_DATA_REAL_TYPE>& vals) const
		{
			std::set< INMOST_DATA_ENUM_TYPE>::const_iterator find;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < Size(); ++k)
			{
				find = indset.find(GetIndex(k));
				assert(find != indset.end());
				vals[std::distance(indset.begin(), find)] += GetValue(k) * coef;
			}
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
					merger.PushRow(alpha, X[it]);
					merger.AddRow(beta, Y[it]);
					merger.RetrieveRow(Z[it]);
					merger.Clear();
				}
			}
		}

#endif //USE_SOLVER
	}
}
