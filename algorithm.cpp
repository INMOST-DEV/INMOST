#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{	
#define CUTOFF 8
	
	static void shortsort(char *lo, char *hi, unsigned width, int (*comp)(const void *, const void *, void *),void (*swap)(void *, void *, void *), void * user_data) 
	{
		char *p, *max;
		
		while (hi > lo) 
		{
			max = lo;
			for (p = lo + width; p <= hi; p += width) if (comp(p, max,user_data) > 0) max = p;
			if( max != hi ) swap(max, hi,user_data);
			hi -= width;
		}
	}
		
	void sort(void *pbase, unsigned num, unsigned width, int (*comp)(const void *, const void *,void *), void (*swap)(void *, void *, void *), void * user_data)
	{
		if( num == 0 ) return;
		char *base = static_cast<char *>(pbase);
		char *lo, *hi;
		char *mid;
		char *l, *h;
		unsigned size;
		char *lostk[30], *histk[30];
		int stkptr;
		
		if (num < 2 || width == 0) return;
		stkptr = 0;
		
		lo = base;
		hi = (char *) base + width * (num - 1);
		
	recurse:
		size = static_cast<unsigned>((hi - lo) / width) + 1;
		
		if (size <= CUTOFF) 
		{
			shortsort(lo, hi, width, comp, swap, user_data);
		} 
		else 
		{
			mid = lo + (size / 2) * width;
			swap(mid, lo,user_data);
			
			l = lo;
			h = hi + width;
			
			while(true)
			{
				do { l += width; } while (l <= hi && comp(l, lo, user_data) <= 0);
				do { h -= width; } while (h > lo && comp(h, lo, user_data) >= 0);
				if (h <= l) break;
				swap(l, h,user_data);
			}
			
			if( lo != h ) swap(lo, h,user_data);
			
			if (h - 1 - lo >= hi - l) 
			{
				if (lo + width < h) 
				{
					lostk[stkptr] = lo;
					histk[stkptr] = h - width;
					++stkptr;
				}
				
				if (l < hi) 
				{
					lo = l;
					goto recurse;
				}
			} 
			else 
			{
				if (l < hi) 
				{
					lostk[stkptr] = l;
					histk[stkptr] = hi;
					++stkptr;
				}
				
				if (lo + width < h) 
				{
					hi = h - width;
					goto recurse;
				}
			}
		}
		
		--stkptr;
		if (stkptr >= 0) 
		{
			lo = lostk[stkptr];
			hi = histk[stkptr];
			goto recurse;
		}
	}
	int binary_search(void * pdata, unsigned num, unsigned width, int comp(const void *,void *), void * user_data)
	{
		if( num == 0 ) return -1;
		char * data = static_cast<char *>(pdata);
		unsigned imin = 0, imax = num-1;
		while (imin < imax)
		{
			int imid = (imax-imin)/2 + imin;
			if ( comp(data+imid*width,user_data) < 0 )
				imin = imid + 1;
			else
				imax = imid;
		}
		if (imax == imin && comp(data+imin*width,user_data) == 0)
			return imin;
		return -1;
	}	



	Mesh::Random::Random(unsigned int seed)
	{
		n = seed;
		a = 1103515245;
		c = 12345;
		m = 1u << (sizeof(unsigned int)*8-1);
	}
	Mesh::Random::Random(const Random & other)
	{
		n = other.n;
		a = other.a;
		c = other.c;
		m = other.m;
	}
	unsigned int Mesh::Random::Number()
	{
		n = (a*n + c)%m;
		return (n << 2) >> 18;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// HASH function 
	// detailed description: http://www.isthe.com/chongo/tech/comp/fnv/
	const size_t InitialFNV32 = 2166136261U;
	const size_t FNVMultiple32 = 16777619;
	size_t hashfunci32(const char * pos, int n)
	{
		size_t hash = InitialFNV32;
		for(int i = 0; i < n; i++)
			hash = (hash ^ pos[i]) * FNVMultiple32;
		return hash;
	}
	/*
	const uint64_t InitialFNV64 = 14695981039346656037U;
	const uint64_t FNVMultiple64 = 1099511628211;
	uint64_t hashfunci64(const char * pos, int n)
	{
		uint64_t hash = InitialFNV64;
		for(int i = 0; i < n; i++)
			hash = (hash ^ pos[i]) * FNVMultiple64;
		return hash;
	}
	*/
	///////////////////////////////////////////////////////////////////////////////////////////////
}
#endif
