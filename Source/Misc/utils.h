#ifndef INMOST_UTILS_H
#define INMOST_UTILS_H

#include "inmost_common.h"
#include <sstream>
#include <string>

namespace INMOST
{

    template<typename T>
    T from_string(std::string str)
	{
        T v;
        std::istringstream ss(str);
        ss >> v;
        return v;
    }

    template<typename T>
    std::string to_string(T value)
	{
        std::stringstream ss;
        ss << value;
        return ss.str();
    }

    std::string string_to_lower(const std::string &str);
	
	void MPIBarrier();

	int MPIGetRank();
	
	class BinaryHeap
	{
		INMOST_DATA_ENUM_TYPE size_max, size;
		std::vector<INMOST_DATA_ENUM_TYPE> heap;
		std::vector<INMOST_DATA_ENUM_TYPE> index;
		INMOST_DATA_REAL_TYPE * keys;
		
		void swap(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
		{
			INMOST_DATA_ENUM_TYPE t = heap[i];
			heap[i] = heap[j];
			heap[j] = t;
			index[heap[i]] = i;
			index[heap[j]] = j;
		}
		void bubbleUp(INMOST_DATA_ENUM_TYPE k)
		{
			while(k > 1 && keys[heap[k/2]] > keys[heap[k]])
			{
				swap(k, k/2);
				k = k/2;
			}
		}
		void bubbleDown(INMOST_DATA_ENUM_TYPE k)
		{
			INMOST_DATA_ENUM_TYPE j;
			while(2*k <= size)
			{
				j = 2*k;
				if(j < size && keys[heap[j]] > keys[heap[j+1]])
					j++;
				if(keys[heap[k]] <= keys[heap[j]])
					break;
				swap(k, j);
				k = j;
			}
		}
	public:
		BinaryHeap(INMOST_DATA_REAL_TYPE * pkeys, INMOST_DATA_ENUM_TYPE len)
		{
			size_max = len;
			keys = pkeys;
			size = 0;
			heap.resize(size_max+1);
			index.resize(size_max+1);
			for(INMOST_DATA_ENUM_TYPE i = 0; i <= size_max; i++)
				index[i] = ENUMUNDEF;
		}
		void PushHeap(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			size++;
			index[i] = size;
			heap[size] = i;
			keys[i] = key;
			bubbleUp(size);
		}
		INMOST_DATA_ENUM_TYPE PopHeap()
		{
			if( size == 0 ) return ENUMUNDEF;
			INMOST_DATA_ENUM_TYPE min = heap[1];
			swap(1, size--);
			bubbleDown(1);
			index[min] = ENUMUNDEF;
			heap[size+1] = ENUMUNDEF;
			return min;
		}
		void DecreaseKey(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			keys[i] = key;
			bubbleUp(index[i]);
		}
		void Clear()
		{
			while( size ) keys[PopHeap()] = std::numeric_limits<INMOST_DATA_REAL_TYPE>::max();
		}
		bool Contains(INMOST_DATA_ENUM_TYPE i) const
		{
			return index[i] != ENUMUNDEF;
		}
		INMOST_DATA_ENUM_TYPE Size() const {return size;}
		bool Empty() const {return size == 0;}
	};
	
	template<typename KeyType, typename Comparator>
	class BinaryHeapCustom
	{
		INMOST_DATA_ENUM_TYPE size_max, size;
		array<INMOST_DATA_ENUM_TYPE> heap;
		array<INMOST_DATA_ENUM_TYPE> index;
		array<KeyType> keys;
		
		void swap(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
		{
			INMOST_DATA_ENUM_TYPE t = heap[i];
			heap[i] = heap[j];
			heap[j] = t;
			index[heap[i]] = i;
			index[heap[j]] = j;
		}
		void bubbleUp(INMOST_DATA_ENUM_TYPE k)
		{
			while(k > 1 && !Comparator()(keys[heap[k/2]],keys[heap[k]]))
			{
				swap(k, k/2);
				k = k/2;
			}
		}
		void bubbleDown(INMOST_DATA_ENUM_TYPE k)
		{
			INMOST_DATA_ENUM_TYPE j;
			while(2*k <= size)
			{
				j = 2*k;
				if(j < size && !Comparator()(keys[heap[j]],keys[heap[j+1]]) )
					j++;
				if(Comparator()(keys[heap[k]],keys[heap[j]]) )
					break;
				swap(k, j);
				k = j;
			}
		}
	public:
		BinaryHeapCustom(INMOST_DATA_ENUM_TYPE len)
		{
			size_max = len;
			size = 0;
			keys.resize(size_max);
			heap.resize(size_max+1);
			index.resize(size_max+1);
			for(INMOST_DATA_ENUM_TYPE i = 0; i <= size_max; i++)
				index[i] = ENUMUNDEF;
		}
		void PushHeap(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			size++;
			index[i] = size;
			heap[size] = i;
			keys[i] = key;
			bubbleUp(size);
		}
		INMOST_DATA_ENUM_TYPE PopHeap()
		{
			if( size == 0 ) return ENUMUNDEF;
			INMOST_DATA_ENUM_TYPE min = heap[1];
			swap(1, size--);
			bubbleDown(1);
			index[min] = ENUMUNDEF;
			heap[size+1] = ENUMUNDEF;
			return min;
		}
		INMOST_DATA_ENUM_TYPE PeekHeap()
		{
			if( size == 0 ) return ENUMUNDEF;
			return heap[1];
		}
		void DecreaseKey(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			keys[i] = key;
			bubbleUp(index[i]);
		}
		void IncreaseKey(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			keys[i] = key;
			bubbleDown(index[i]);
		}
		void ChangeKey(INMOST_DATA_ENUM_TYPE i, INMOST_DATA_REAL_TYPE key)
		{
			keys[i] = key;
			if( !Comparator()(keys[i],key) )	
				bubbleUp(index[i]);
			else
				bubbleDown(index[i]);
		}
		KeyType & GetKey(INMOST_DATA_ENUM_TYPE i) {return keys[i];}
		const KeyType & GetKey(INMOST_DATA_ENUM_TYPE i) const {return keys[i];}
		void Clear()
		{
			while( size ) PopHeap(); 
		}
		bool Contains(INMOST_DATA_ENUM_TYPE i) const
		{
			return index[i] != ENUMUNDEF;
		}
		INMOST_DATA_ENUM_TYPE Size() const {return size;}
		bool Empty() const {return size == 0;}
	};
	
	//this structure is used to provide sorting routing for Reverse-Cuthill-McKee reordering
	struct RCM_Comparator
	{
		INMOST_DATA_ENUM_TYPE wbeg;
		std::vector<INMOST_DATA_ENUM_TYPE> & xadj;
	public:
		RCM_Comparator(INMOST_DATA_ENUM_TYPE _wbeg, std::vector<INMOST_DATA_ENUM_TYPE> & _xadj)
		: wbeg(_wbeg), xadj(_xadj) {}
		RCM_Comparator(const RCM_Comparator & b) : wbeg(b.wbeg), xadj(b.xadj) {}
		RCM_Comparator & operator = (RCM_Comparator const & b) {wbeg = b.wbeg; xadj = b.xadj; return *this;}
		bool operator () (INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
		{return xadj[i+1-wbeg] -xadj[i-wbeg] < xadj[j+1-wbeg] - xadj[j-wbeg];}
		
		bool operator () (const std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> & i, const std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> & j)
		{
			INMOST_DATA_ENUM_TYPE ci = xadj[i.first+1-wbeg] - xadj[i.first-wbeg];
			INMOST_DATA_ENUM_TYPE cj = xadj[j.first+1-wbeg] - xadj[j.first-wbeg];
			return (ci < cj) || (ci == cj && i.second < i.first);
		}
	};
	
	struct WRCM_Comparator
	{
		INMOST_DATA_ENUM_TYPE wbeg;
		std::vector<INMOST_DATA_REAL_TYPE> & wadj;
	public:
		WRCM_Comparator(INMOST_DATA_ENUM_TYPE _wbeg, std::vector<INMOST_DATA_REAL_TYPE> & _wadj)
		: wbeg(_wbeg), wadj(_wadj) {}
		WRCM_Comparator(const WRCM_Comparator & b) : wbeg(b.wbeg), wadj(b.wadj) {}
		WRCM_Comparator & operator = (WRCM_Comparator const & b) {wbeg = b.wbeg; wadj = b.wadj; return *this;}
		bool operator () (INMOST_DATA_ENUM_TYPE i, INMOST_DATA_ENUM_TYPE j)
		{return wadj[i-wbeg] < wadj[j-wbeg];}
		
		bool operator () (const std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> & i, const std::pair<INMOST_DATA_ENUM_TYPE,INMOST_DATA_REAL_TYPE> & j)
		{
			INMOST_DATA_REAL_TYPE ci = wadj[i.first-wbeg];
			INMOST_DATA_REAL_TYPE cj = wadj[j.first-wbeg];
			return (ci < cj) || (ci == cj && i.second < i.first);
		}
	};
	
	
}

#endif //INMOST_UTILS_H
