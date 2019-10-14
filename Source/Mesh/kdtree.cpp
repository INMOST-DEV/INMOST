#include "inmost_mesh.h"
#if defined(USE_MESH)
#define _m0(x)	(x & 0x7FF)
#define _m1(x)	(x >> 11 & 0x7FF)
#define _m2(x)	(x >> 22 )

namespace INMOST
{
	inline static void crossproduct(const Storage::real vecin1[3],const Storage::real vecin2[3], Storage::real vecout[3])
	{
		vecout[0] = vecin1[1]*vecin2[2] - vecin1[2]*vecin2[1];
		vecout[1] = vecin1[2]*vecin2[0] - vecin1[0]*vecin2[2];
		vecout[2] = vecin1[0]*vecin2[1] - vecin1[1]*vecin2[0];
	}
	
	inline static Storage::real dotproduct(const  Storage::real * vecin1, const Storage::real * vecin2)
	{
		return vecin1[0]*vecin2[0]+vecin1[1]*vecin2[1]+vecin1[2]*vecin2[2];
	}
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
					Storage::real_array coords = m->RealArrayDF(nodes[k],m->CoordsTag());;
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
	
	inline int SearchKDTree::segment_tri(const Storage::real tri[3][3], const Storage::real p1[3], const Storage::real p2[3])
	{
		const Storage::real eps = 1.0e-7;
		Storage::real a[3],b[3],c[3],n[3], ray[3], d, k, m, l;
		Storage::real dot00,dot01,dot02,dot11,dot12,invdenom, uq,vq;
		ray[0] = p2[0]-p1[0];
		ray[1] = p2[1]-p1[1];
		ray[2] = p2[2]-p1[2];
		l = sqrt(dotproduct(ray,ray));
		ray[0] /= l;
		ray[1] /= l;
		ray[2] /= l;
		a[0] = tri[0][0] - tri[2][0];
		a[1] = tri[0][1] - tri[2][1];
		a[2] = tri[0][2] - tri[2][2];
		b[0] = tri[1][0] - tri[2][0];
		b[1] = tri[1][1] - tri[2][1];
		b[2] = tri[1][2] - tri[2][2];
		crossproduct(a,b,n);
		d = -dotproduct(n,tri[0]);
		m =  dotproduct(n,ray);
		if( fabs(m) < 1.0e-25 )
			return 0;
		k = -(d + dotproduct(n,p1))/m;
		if( k < 0 )
			return 0;
		if( k > l )
			return 0;
		c[0] = p1[0] + k*ray[0] - tri[2][0];
		c[1] = p1[1] + k*ray[1] - tri[2][1];
		c[2] = p1[2] + k*ray[2] - tri[2][2];
		dot00 = dotproduct(a,a);
		dot01 = dotproduct(a,b);
		dot02 = dotproduct(a,c);
		dot11 = dotproduct(b,b);
		dot12 = dotproduct(b,c);
		invdenom = (dot00*dot11 - dot01*dot01);
		uq = (dot11*dot02-dot01*dot12);
		vq = (dot00*dot12-dot01*dot02);
		if( fabs(invdenom) < 1.0e-25 && fabs(uq) > 0.0 && fabs(vq) > 0.0 )
			return 0;
		uq = uq/invdenom;
		vq = vq/invdenom;
		if( uq >= -eps && vq >= -eps && 1.0-(uq+vq) >= -eps )
			return 1;
		return 0;
	}
	
	inline int SearchKDTree::segment_bbox(const Storage::real p1[3], const Storage::real p2[3])
	{
		Storage::real tnear = -1.0e20, tfar = 1.0e20;
		Storage::real t1,t2,c;
		for(int i = 0; i < 3; i++)
		{
			if( fabs(p2[i]-p1[i]) < 1.0e-15 )
			{
				if( p1[i] < bbox[i*2] || p1[i] > bbox[i*2+1] )
					return 0;
			}
			else
			{
				t1 = (bbox[i*2+0] - p1[i])/(p2[i]-p1[i]);
				t2 = (bbox[i*2+1] - p1[i])/(p2[i]-p1[i]);
				if( t1 > t2 )
				{
					c = t1;
					t1 = t2;
					t2 = c;
				}
				if( t1 > tnear ) tnear = t1;
				if( t2 < tfar ) tfar = t2;
				if( tnear > tfar ) return 0;
				if( tfar < 0 ) return 0;
				if( tnear > 1 ) return 0;
			}
		}
		return 1;
	}
	
	inline bool SearchKDTree::segment_face(const Element & f, const Storage::real p1[3], const Storage::real p2[3])
	{
		Mesh * m = f->GetMeshLink();
		Storage::real tri[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
		m->GetGeometricData(f->GetHandle(),CENTROID,tri[2]);
		ElementArray<Node> nodes = f->getNodes();
		m->GetGeometricData(nodes[0]->GetHandle(),CENTROID,tri[1]);
		for(ElementArray<Node>::size_type k = 0; k < nodes.size(); ++k)
		{
			memcpy(tri[0],tri[1],sizeof(Storage::real)*3);
			//m->GetGeometricData(nodes[k]->GetHandle(),CENTROID,tri[0]);
			m->GetGeometricData(nodes[(k+1)%nodes.size()]->GetHandle(),CENTROID,tri[1]);
			if( segment_tri(tri,p1,p2) ) return true;
		}
		return false;
	}
	
	inline bool SearchKDTree::segment_cell(const Element & c, const Storage::real p1[3], const Storage::real p2[3])
	{
		ElementArray<Face> faces = c->getFaces();
		for(ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
			if( segment_face(it->self(),p1,p2) ) return true;
		return false;
	}
	
	void SearchKDTree::IntersectSegment(ElementArray<Cell> & cells, const Storage::real p1[3], const Storage::real p2[3])
	{
		ElementArray<Element> temp(m);
		MarkerType mrk = m->CreateMarker();
		sub_intersect_segment(temp, mrk, p1, p2);
		m->RemMarkerArray(temp.data(), static_cast<Storage::enumerator>(temp.size()), mrk);
		for (ElementArray<Element>::iterator it = temp.begin(); it != temp.end(); ++it)
		{
			if (it->GetElementType() == CELL && !it->GetMarker(mrk))
			{
				cells.push_back(it->getAsCell());
				it->SetMarker(mrk);
			}
			else if (it->GetElementType() == FACE)
			{
				ElementArray<Cell> f_cells = it->getCells();
				for (ElementArray<Cell>::iterator kt = f_cells.begin(); kt != f_cells.end(); ++kt) if (!kt->GetMarker(mrk))
				{
					cells.push_back(kt->self());
					kt->SetMarker(mrk);
				}
			}
		}
		m->RemMarkerArray(cells.data(), static_cast<Storage::enumerator>(cells.size()), mrk);
		m->ReleaseMarker(mrk);
	}
	
	
	void SearchKDTree::sub_intersect_segment(ElementArray<Element> & hits, MarkerType mrk, const Storage::real p1[3], const Storage::real p2[3])
	{
		if( size == 1 )
		{
			if( !m->GetMarker(set[0].e,mrk) )
			{
				Storage::integer edim = Element::GetGeometricDimension(m->GetGeometricType(set[0].e));
				if( edim == 2 )
				{
					if( segment_face(Element(m,set[0].e),p1,p2) )
					{
						hits.push_back(set[0].e);
						m->SetMarker(set[0].e,mrk);
					}
				}
				else if( edim == 3 )
				{
					if( segment_cell(Element(m,set[0].e),p1,p2) )
					{
						hits.push_back(set[0].e);
						m->SetMarker(set[0].e,mrk);
					}
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " kd-tree structure is not suitable to intersect edges with segments" << std::endl;
					exit(-1);
				}
			}
		}
		else
		{
			assert(size > 1);
			if( segment_bbox(p1,p2) )
			{
				children[0].sub_intersect_segment(hits,mrk,p1,p2);
				children[1].sub_intersect_segment(hits,mrk,p1,p2);
			}
		}
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
