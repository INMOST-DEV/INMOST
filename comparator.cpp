#include "inmost_mesh.h"
#if defined(USE_MESH)
namespace INMOST
{
	int Mesh::CentroidComparator::Compare(const real * a, const real * b)
	{
		real e = m->GetEpsilon();
		for(integer i = 0; i < m->GetDimensions(); i++) 
			if( fabs(a[i]-b[i]) > e ) 
			{
				if( a[i] < b[i] ) 
					return -1;
				else
					return 1;
			}
		return 0;
	}

	bool Mesh::CentroidComparator::operator () (HandleType a, HandleType b)
	{
		if( a == InvalidHandle() || b == InvalidHandle() ) return a > b;
		real ca[3] = {0,0,0}, cb[3] = {0,0,0};
		m->GetGeometricData(a,CENTROID,ca);
		m->GetGeometricData(b,CENTROID,ca);
		return Compare(ca,cb) <= 0;
	}

	bool Mesh::CentroidComparator::operator () (HandleType a, const real * cb)
	{
		if( a == InvalidHandle() ) return true;
		real ca[3] = {0,0,0}; 
		m->GetGeometricData(a,CENTROID,ca);
		return Compare(ca,cb) <= 0;
	}

	
	int Mesh::IerarhyComparator::CompareNodes(HandleType a, HandleType b)
	{
		real_array ca = m->RealArrayDF(a,m->CoordsTag());
		real_array cb = m->RealArrayDF(b,m->CoordsTag());
		real e = m->GetEpsilon();
		for(integer i = 0; i < m->GetDimensions(); ++i)
			if( fabs(ca[i]-cb[i]) > e )
			{
				if( ca[i] > cb[i] ) 
					return 1;
				else
					return -1;
			}
		return 0;
	}
	int Mesh::IerarhyComparator::CompareElements(HandleType a, HandleType b)
	{
		integer ia = GetHandleElementNum(a);
		integer ib = GetHandleElementNum(b);
		if( ia == ib )
		{
			if( ia == 0 ) return CompareNodes(a,b);
			Element::adj_type const & adja = m->LowConn(a);
			Element::adj_type const & adjb = m->LowConn(b);
			if( adja.size() == adjb.size() )
			{
				for(Element::adj_type::size_type k = 0; k < adja.size(); ++k)
				{
					int ret = CompareElements(adja[k],adjb[k]);
					if( ret != 0 ) return ret;
				}
				return 0;
			}
			else return static_cast<int>(adja.size()) - static_cast<int>(adjb.size());
		}
		else return ia - ib;
	}

	void Mesh::SortByGlobalID(HandleType * h, enumerator n)
	{
		if( sizeof(integer) == 4 && n > 64 )
		{
			Element::adj_type::size_type i, num = static_cast<Element::adj_type::size_type>(n);
			array<HandleType> tmp(num);
			HandleType * ptra = h;
			HandleType * ptrb = tmp.data();
			const unsigned int kHist = 256;
			unsigned int b0[kHist * 4];
			unsigned int *b1 = b0 + kHist;
			unsigned int *b2 = b0 + kHist*2;
			unsigned int *b3 = b0 + kHist*3;
			memset(b0,0,sizeof(unsigned int)*kHist*4);
			for (i = 0; i < num; i++) 
			{
				unsigned char * rec = reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]));
				++b0[rec[0]]; 
				++b1[rec[1]];
				++b2[rec[2]]; 
				++b3[rec[3]];
			}
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0, tsum;
			for (i = 0; i < kHist; i++) 
			{
				tsum = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = tsum;
				tsum = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = tsum;
				tsum = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = tsum;
				tsum = b3[i] + sum3; b3[i] = sum3 - 1; sum3 = tsum;
			}
			for (i = 0; i < num; i++) ptrb[++b0[reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]))[0]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b1[reinterpret_cast<unsigned char *>(&GlobalID(ptrb[i]))[1]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b2[reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]))[2]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b3[reinterpret_cast<unsigned char *>(&GlobalID(ptrb[i]))[3]]] = ptrb[i];
		}
		else if( sizeof(integer) == 8 && n > 128 )
		{
			Element::adj_type::size_type i, num = static_cast<Element::adj_type::size_type>(n);
			array<HandleType> tmp(num);
			HandleType * ptra = h;
			HandleType * ptrb = tmp.data();
			const unsigned int kHist = 256;
			unsigned int b0[kHist * 8];
			unsigned int *b1 = b0 + kHist;
			unsigned int *b2 = b0 + kHist*2;
			unsigned int *b3 = b0 + kHist*3;
			unsigned int *b4 = b0 + kHist*4;
			unsigned int *b5 = b0 + kHist*5;
			unsigned int *b6 = b0 + kHist*6;
			unsigned int *b7 = b0 + kHist*7;
			memset(b0,0,sizeof(unsigned int)*kHist*8);
			for (i = 0; i < num; i++) 
			{
				unsigned char * rec = reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]));
				++b0[rec[0]]; 
				++b1[rec[1]];
				++b2[rec[2]]; 
				++b3[rec[3]];
				++b4[rec[4]]; 
				++b5[rec[5]];
				++b6[rec[6]]; 
				++b7[rec[7]];
			}
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0, tsum;
			unsigned int sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0;
			for (i = 0; i < kHist; i++) 
			{
				tsum = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = tsum;
				tsum = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = tsum;
				tsum = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = tsum;
				tsum = b3[i] + sum3; b3[i] = sum3 - 1; sum3 = tsum;
				tsum = b4[i] + sum4; b4[i] = sum4 - 1; sum4 = tsum;
				tsum = b5[i] + sum5; b5[i] = sum5 - 1; sum5 = tsum;
				tsum = b6[i] + sum6; b6[i] = sum6 - 1; sum6 = tsum;
				tsum = b7[i] + sum7; b7[i] = sum7 - 1; sum7 = tsum;
			}
			for (i = 0; i < num; i++) ptrb[++b0[reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]))[0]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b1[reinterpret_cast<unsigned char *>(&GlobalID(ptrb[i]))[1]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b2[reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]))[2]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b3[reinterpret_cast<unsigned char *>(&GlobalID(ptrb[i]))[3]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b4[reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]))[4]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b5[reinterpret_cast<unsigned char *>(&GlobalID(ptrb[i]))[5]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b6[reinterpret_cast<unsigned char *>(&GlobalID(ptra[i]))[6]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b7[reinterpret_cast<unsigned char *>(&GlobalID(ptrb[i]))[7]]] = ptrb[i];
		}
		else std::sort(h,h+n,Mesh::GlobalIDComparator(this));
	}

	void Mesh::SortHandles(HandleType * h, enumerator n)
	{
		if( sizeof(HandleType) == 4 && n > 32 )
		{
			Element::adj_type::size_type i, num = static_cast<Element::adj_type::size_type>(n);
			array<HandleType> tmp(num);
			HandleType * ptra = h;
			HandleType * ptrb = tmp.data();
			//make so that invalid handles appear at he end
			for (i = 0; i < num; i++) if( ptra[i] == InvalidHandle() ) ptra[i] = ENUMUNDEF;
			const unsigned int kHist = 256;
			unsigned int b0[kHist * 4];
			unsigned int *b1 = b0 + kHist;
			unsigned int *b2 = b0 + kHist*2;
			unsigned int *b3 = b0 + kHist*3;
			unsigned char * l1 = reinterpret_cast<unsigned char *>(ptra);
			unsigned char * l2 = reinterpret_cast<unsigned char *>(ptrb);
			memset(b0,0,sizeof(unsigned int)*kHist*4);
			for (i = 0; i < num; i++) 
			{
				++b0[l1[4*i+0]]; 
				++b1[l1[4*i+1]];
				++b2[l1[4*i+2]]; 
				++b3[l1[4*i+3]];
			}
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0, tsum;
			for (i = 0; i < kHist; i++) 
			{
				tsum = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = tsum;
				tsum = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = tsum;
				tsum = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = tsum;
				tsum = b3[i] + sum3; b3[i] = sum3 - 1; sum3 = tsum;
			}
			for (i = 0; i < num; i++) ptrb[++b0[l1[4*i+0]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b1[l2[4*i+1]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b2[l1[4*i+2]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b3[l2[4*i+3]]] = ptrb[i];

			for (i = 0; i < num; i++) if( ptra[i] == ENUMUNDEF ) ptra[i] = InvalidHandle();
		}
		else if( sizeof(HandleType) == 8 && n > 64 )
		{
			Element::adj_type::size_type i, num = static_cast<Element::adj_type::size_type>(n);
			array<HandleType> tmp(num);
			HandleType * ptra = h;
			HandleType * ptrb = tmp.data();
			//make so that invalid handles appear at he end
			for (i = 0; i < num; i++) if( ptra[i] == InvalidHandle() ) ptra[i] = ENUMUNDEF;
			const unsigned int kHist = 256;
			unsigned int b0[kHist * 8];
			unsigned int *b1 = b0 + kHist;
			unsigned int *b2 = b0 + kHist*2;
			unsigned int *b3 = b0 + kHist*3;
			unsigned int *b4 = b0 + kHist*4;
			unsigned int *b5 = b0 + kHist*5;
			unsigned int *b6 = b0 + kHist*6;
			unsigned int *b7 = b0 + kHist*7;
			unsigned char * l1 = reinterpret_cast<unsigned char *>(ptra);
			unsigned char * l2 = reinterpret_cast<unsigned char *>(ptrb);
			memset(b0,0,sizeof(unsigned int)*kHist*8);
			for (i = 0; i < num; i++) 
			{
				++b0[l1[8*i+0]]; 
				++b1[l1[8*i+1]];
				++b2[l1[8*i+2]]; 
				++b3[l1[8*i+3]];
				++b4[l1[8*i+4]]; 
				++b5[l1[8*i+5]];
				++b6[l1[8*i+6]]; 
				++b7[l1[8*i+7]];
			}
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0, tsum;
			unsigned int sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0;
			for (i = 0; i < kHist; i++) 
			{
				tsum = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = tsum;
				tsum = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = tsum;
				tsum = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = tsum;
				tsum = b3[i] + sum3; b3[i] = sum3 - 1; sum3 = tsum;
				tsum = b4[i] + sum4; b4[i] = sum4 - 1; sum4 = tsum;
				tsum = b5[i] + sum5; b5[i] = sum5 - 1; sum5 = tsum;
				tsum = b6[i] + sum6; b6[i] = sum6 - 1; sum6 = tsum;
				tsum = b7[i] + sum7; b7[i] = sum7 - 1; sum7 = tsum;
			}
			for (i = 0; i < num; i++) ptrb[++b0[l1[8*i+0]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b1[l2[8*i+1]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b2[l1[8*i+2]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b3[l2[8*i+3]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b4[l1[8*i+4]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b5[l2[8*i+5]]] = ptrb[i];
			for (i = 0; i < num; i++) ptrb[++b6[l1[8*i+6]]] = ptra[i];
			for (i = 0; i < num; i++) ptra[++b7[l2[8*i+7]]] = ptrb[i];

			for (i = 0; i < num; i++) if( ptra[i] == ENUMUNDEF ) ptra[i] = InvalidHandle();
		}
		else std::sort(h,h+n);
	}
}
#endif