#include "inmost_mesh.h"
#if defined(USE_MESH)
#define _m0(x)	(x & 0x7FF)
#define _m1(x)	(x >> 11 & 0x7FF)
#define _m2(x)	(x >> 22 )

namespace INMOST
{
	template<typename bbox_type>
	inline int SearchKDTree::bbox_point(const Storage::real p[3], const bbox_type bbox[6], bool print)
	{
		for(int i = 0; i < 3; i++)
		{
			if( p[i] < bbox[i*2]-1.0e-3 || p[i] > bbox[i*2+1]+1.0e-3 )
			{
				if( print ) 
				{
					std::cout << "fail on " << i << " ";
					if( p[i] < bbox[i*2]-1.0e-3 ) std::cout << p[i] << " is less then " << bbox[i*2]-1.0e-3 << "(" << bbox[i*2] << ")";
					if( p[i] > bbox[i*2+1]+1.0e-3 ) std::cout << p[i] << " is greater then " << bbox[i*2+1]+1.0e-3 << "(" << bbox[i*2+1] << ")";
					std::cout << std::endl;
				}	
				return 0;
			}
		}
		if( print ) std::cout << "all ok!" << std::endl;
		return 1;
	}
	
	int SearchKDTree::cmpElements0(const void * a,const void * b)
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[0];
		float bd = eb->xyz[0];
		return (ad > bd) - (ad < bd);
	}
	
	int SearchKDTree::cmpElements1(const void * a,const void * b)
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[1];
		float bd = eb->xyz[1];
		return (ad > bd) - (ad < bd);
	}
	
	int SearchKDTree::cmpElements2(const void * a,const void * b)
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[2];
		float bd = eb->xyz[2];
		return (ad > bd) - (ad < bd);
	}
	
	inline unsigned int SearchKDTree::flip(const unsigned int * fp)
	{
		unsigned int mask = -((int)(*fp >> 31)) | 0x80000000;
		return *fp ^ mask;
	}
	
	void SearchKDTree::radix_sort(int dim, struct entry * temp)
	{
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0,0,sizeof(unsigned int)*kHist*3);
		for (i = 0; i < size; i++)
		{
			unsigned int fi = flip((unsigned int *)&set[i].xyz[dim]);
			++b0[_m0(fi)]; ++b1[_m1(fi)]; ++b2[_m2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++)
			{
				b0[kHist-1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist-1];
				b1[kHist-1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist-1];
				b2[kHist-1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist-1];
			}
		}
		for (i = 0; i < size; i++) temp[++b0[_m0(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[++b1[_m1(flip((unsigned int *)&temp[i].xyz[dim]))]] = temp[i];
		for (i = 0; i < size; i++) temp[++b2[_m2(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[i] = temp[i];
	}
	
	void SearchKDTree::kdtree_build(int dim, int & done, int total, struct entry * temp)
	{
		if( size > 1 )
		{
			if( size > 128 ) radix_sort(dim,temp); else
				switch(dim)
			{
				case 0: qsort(set,size,sizeof(entry),cmpElements0);break;
				case 1: qsort(set,size,sizeof(entry),cmpElements1);break;
				case 2: qsort(set,size,sizeof(entry),cmpElements2);break;
			}
			children = static_cast<SearchKDTree *>(malloc(sizeof(SearchKDTree)*2));//new kdtree[2];
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size/2;
			children[0].m = m;
			children[1].children = NULL;
			children[1].set = set+size/2;
			children[1].size = size - size/2;
			children[1].m = m;
			children[0].kdtree_build((dim+1)%3,done,total,temp);
			children[1].kdtree_build((dim+1)%3,done,total,temp);
			for(int k = 0; k < 3; k++)
			{
				bbox[0+2*k] = std::min<float>(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
				bbox[1+2*k] = std::max<float>(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
			}
		}
		else
		{
			assert(size == 1);
			{
				bool checkm = false, invm = false;
				if( m->HideMarker()  ) 
				{
					checkm = true;
					if( m->GetMarker(set[0].e,m->HideMarker()) ) invm = true;
				}
				Element::adj_type & nodes = m->HighConn(set[0].e); 
				bbox[0] = bbox[2] = bbox[4] = 1.0e20f;
				bbox[1] = bbox[3] = bbox[5] = -1.0e20f;
				if( checkm )
				{
					for(unsigned k = 0; k < nodes.size(); ++k)
					{
						bool hidn = m->GetMarker(nodes[k],m->HideMarker());
						bool newn = m->GetMarker(nodes[k],m->NewMarker());
						if( invm && newn ) continue;
						if( !invm && hidn ) continue;
						Storage::real_array coords = m->RealArrayDF(nodes[k],m->CoordsTag());;
						for(unsigned q = 0; q < coords.size(); q++)
						{
							bbox[q*2+0] = std::min<float>(bbox[q*2+0],(float)coords[q]);
							bbox[q*2+1] = std::max<float>(bbox[q*2+1],(float)coords[q]);
						}
					}
				}
				else for(unsigned k = 0; k < nodes.size(); ++k)
				{
					Storage::real_array coords = m->RealArrayDF(set[0].e,m->CoordsTag());;
					for(unsigned q = 0; q < coords.size(); q++)
					{
						bbox[q*2+0] = std::min<float>(bbox[q*2+0],(float)coords[q]);
						bbox[q*2+1] = std::max<float>(bbox[q*2+1],(float)coords[q]);
					}
				}
				for(int k = m->GetDimensions(); k < 3; ++k)
				{
					bbox[k*2+0] = -1.0e20f;
					bbox[k*2+1] = 1.0e20f;
				}
			}
		}
	}
	
	Cell SearchKDTree::SubSearchCell(const Storage::real p[3], bool print)
	{
		Cell ret = InvalidCell();
		if( size == 1 )
		{
			if( print ) std::cout << "test cell " << GetHandleID(set[0].e) << std::endl;
			if( m->HideMarker() && m->GetMarker(set[0].e,m->HideMarker()) ) 
			{
				m->SwapModification(false);
				if( cell_point(Cell(m,set[0].e),p) )
					ret = Cell(m,set[0].e);
				m->SwapModification(false);
			}
			else
			{
				if( cell_point(Cell(m,set[0].e),p) )
					ret = Cell(m,set[0].e);
			}
		}
		else
		{
			assert(size > 1);
			if( bbox_point(p,bbox) )
			{
				if( print )
				{
					std::cout << "point " << p[0] << " " << p[1] << " " << p[2] << " is in bbox ";
					std::cout << " x " << bbox[0] << ":" << bbox[1];
					std::cout << " y " << bbox[2] << ":" << bbox[3];
					std::cout << " z " << bbox[4] << ":" << bbox[5];
					std::cout << std::endl;
				}
				if( print ) std::cout << "try left child" << std::endl;
				ret = children[0].SubSearchCell(p,print);
				if( print ) std::cout << "ret " << (ret.isValid() ? ret.LocalID() : -1) << std::endl;
				if( !ret.isValid() )
				{
					if( print ) std::cout << "try right child" << std::endl;
					ret = children[1].SubSearchCell(p,print);
					if( print ) std::cout << "ret " << (ret.isValid() ? ret.LocalID() : -1) << std::endl;
				}
			}
			else if( print ) 
			{
				std::cout << "point " << p[0] << " " << p[1] << " " << p[2] << " is not in bbox ";
				std::cout << " x " << bbox[0] << ":" << bbox[1];
				std::cout << " y " << bbox[2] << ":" << bbox[3];
				std::cout << " z " << bbox[4] << ":" << bbox[5];
				std::cout << std::endl;
				bbox_point(p,bbox,true);
			}
		}
		return ret;
	}
	
	void SearchKDTree::clear_children() 
	{ 
		if( children ) 
		{
			children[0].clear_children(); 
			children[1].clear_children(); 
			free(children);
		}
	}
	
	SearchKDTree::SearchKDTree(Mesh * m) : m(m), children(NULL)
	{
		size = m->NumberOfCells();
		if( size )
		{
			set = new entry[size];
			INMOST_DATA_ENUM_TYPE k = 0;
			for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
			{
				set[k].e = it->GetHandle();
				Storage::real cnt[3] = {0,0,0};
				m->GetGeometricData(set[k].e,CENTROID,cnt);
				set[k].xyz[0] = (float)cnt[0];
				set[k].xyz[1] = (float)cnt[1];
				set[k].xyz[2] = (float)cnt[2];
				++k;
			}
			int done = 0, total = size;
			struct entry *  temp = new entry[size];
			kdtree_build(0,done,total,temp);
			delete [] temp;
			if( size > 1 )
			{
				for(int k = 0; k < 3; k++)
				{
					bbox[0+2*k] = std::min<float>(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
					bbox[1+2*k] = std::max<float>(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
				}
			}
		}
	}
	
	
	SearchKDTree::SearchKDTree(Mesh * m, HandleType * _set, unsigned set_size) : m(m), children(NULL)
	{
		size = set_size;
		if( size )
		{
			set = new entry[size];
			for(unsigned k = 0; k < set_size; ++k)
			{
				set[k].e = _set[k];
				Storage::real cnt[3] = {0,0,0};
				m->GetGeometricData(set[k].e,CENTROID,cnt);
				set[k].xyz[0] = (float)cnt[0];
				set[k].xyz[1] = (float)cnt[1];
				set[k].xyz[2] = (float)cnt[2];
			}
			int done = 0, total = size;
			struct entry *  temp = new entry[size];
			kdtree_build(0,done,total,temp);
			delete [] temp;
			if( size > 1 )
			{
				for(int k = 0; k < 3; k++)
				{
					bbox[0+2*k] = std::min<float>(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
					bbox[1+2*k] = std::max<float>(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
				}
			}
		}
		else set = NULL;
	}
	
	SearchKDTree::~SearchKDTree()
	{
		if( size )
		{
			delete [] set;
			clear_children();
		}
	}
	
	Cell SearchKDTree::SearchCell(const Storage::real * point, bool print)
	{
		return SubSearchCell(point, print);
	}
}
#endif
