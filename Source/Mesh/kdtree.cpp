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
	inline int SearchKDTree::bbox_point(const Storage::real p[3], const bbox_type bbox[6])
	{
		for(int i = 0; i < 3; i++)
		{
			if( p[i] < bbox[i*2]-1.0e-3 || p[i] > bbox[i*2+1]+1.0e-3 )
				return 0;
		}
		return 1;
	}

	template<typename bbox_type>
	inline int SearchKDTree::bbox_point_print(const Storage::real p[3], const bbox_type bbox[6], std::ostream & sout)
	{
		for (int i = 0; i < 3; i++)
		{
			if (p[i] < bbox[i * 2] - 1.0e-3 || p[i] > bbox[i * 2 + 1] + 1.0e-3)
			{
				sout << "fail on " << i << " ";
				if (p[i] < bbox[i * 2] - 1.0e-3) sout << p[i] << " is less than " << bbox[i * 2] - 1.0e-3 << "(" << bbox[i * 2] << ")";
				if (p[i] > bbox[i * 2 + 1] + 1.0e-3) sout << p[i] << " is greater then " << bbox[i * 2 + 1] + 1.0e-3 << "(" << bbox[i * 2 + 1] << ")";
				sout << std::endl;
				return 0;
			}
		}
		sout << "all ok!" << std::endl;
		return 1;
	}


	template<typename bbox_type>
	inline void SearchKDTree::bbox_closest_point(const Storage::real p[3], const bbox_type bbox[6], Storage::real pout[3])
	{
		for (int i = 0; i < 3; i++)
		{
			if (p[i] < bbox[i * 2] )
				pout[i] = bbox[i * 2];
			else if (p[i] > bbox[i * 2 + 1])
				pout[i] = bbox[i * 2 + 1];
			else
				pout[i] = p[i];
		}
		return;
	}

	template<typename bbox_type>
	inline int SearchKDTree::bbox_sphere(const Storage::real p[3], Storage::real r, const bbox_type bbox[6])
	{
		Storage::real pb[3], d;
		bbox_closest_point(p, bbox, pb);
		d = sqrt((pb[0] - p[0]) * (pb[0] - p[0]) + (pb[1] - p[1]) * (pb[1] - p[1]) + (pb[2] - p[2]) * (pb[2] - p[2]));
		return d <= r ? 1 : 0;
	}

	inline Storage::real SearchKDTree::segment_distance(const Storage::real a[3], const Storage::real b[3], const Storage::real p[3])
	{
		Storage::real pout[3];
		Storage::real v[3] = { b[0] - a[0],b[1] - a[1],b[2] - a[2] };
		Storage::real d = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
		Storage::real t = (v[0] * (p[0] - a[0]) + v[1] * (p[1] - a[1]) + v[2] * (p[2] - a[2])) / d;
		if (t >= 0.0 && t <= 1.0)
		{
			pout[0] = a[0] + t * v[0];
			pout[1] = a[1] + t * v[1];
			pout[2] = a[2] + t * v[2];
		}
		else if (t < 0.0)
		{
			pout[0] = a[0];
			pout[1] = a[1];
			pout[2] = a[2];
		}
		else
		{
			pout[0] = b[0];
			pout[1] = b[1];
			pout[2] = b[2];
		}
		return sqrt((pout[0] - p[0]) * (pout[0] - p[0]) + (pout[1] - p[1]) * (pout[1] - p[1]) + (pout[2] - p[2]) * (pout[2] - p[2]));
	}

	inline Storage::real SearchKDTree::triangle_distance(const Storage::real a[3], const Storage::real b[3], const Storage::real c[3], const Storage::real p[3])
	{
		Storage::real pout[3];
		Storage::real n[3], d, t, q1[3], q2[3], q3[3], q12, q23, q13;
		n[0] = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]);
		n[1] = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]);
		n[2] = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
		d = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (d)
		{
			n[0] /= d;
			n[1] /= d;
			n[2] /= d;
		}
		t = (a[0] - p[0]) * n[0] + (a[1] - p[1]) * n[1] + (a[2] - p[2]) * n[2];
		pout[0] = p[0] + t * n[0];
		pout[1] = p[1] + t * n[1];
		pout[2] = p[2] + t * n[2];
		q1[0] = (pout[1] - a[1]) * (pout[2] - b[2]) - (pout[1] - b[1]) * (pout[2] - a[2]);
		q1[1] = (pout[2] - a[2]) * (pout[0] - b[0]) - (pout[2] - b[2]) * (pout[0] - a[0]);
		q1[2] = (pout[0] - a[0]) * (pout[1] - b[1]) - (pout[0] - b[0]) * (pout[1] - a[1]);
		q2[0] = (pout[1] - b[1]) * (pout[2] - c[2]) - (pout[1] - c[1]) * (pout[2] - b[2]);
		q2[1] = (pout[2] - b[2]) * (pout[0] - c[0]) - (pout[2] - c[2]) * (pout[0] - b[0]);
		q2[2] = (pout[0] - b[0]) * (pout[1] - c[1]) - (pout[0] - c[0]) * (pout[1] - b[1]);
		q3[0] = (pout[1] - c[1]) * (pout[2] - a[2]) - (pout[1] - a[1]) * (pout[2] - c[2]);
		q3[1] = (pout[2] - c[2]) * (pout[0] - a[0]) - (pout[2] - a[2]) * (pout[0] - c[0]);
		q3[2] = (pout[0] - c[0]) * (pout[1] - a[1]) - (pout[0] - a[0]) * (pout[1] - c[1]);
		q12 = q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2];
		q23 = q2[0] * q3[0] + q2[1] * q3[1] + q2[2] * q3[2];
		q13 = q1[0] * q3[0] + q1[1] * q3[1] + q1[2] * q3[2];
		if (q12 >= 0 && q23 >= 0 && q13 >= 0)
			return fabs(t);
		Storage::real dmin = segment_distance(a, b, p);
		dmin = std::min(dmin, segment_distance(b, c, p));
		dmin = std::min(dmin, segment_distance(c, a, p));
		return dmin;
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
				std::vector<HandleType> nodes;
				if (GetHandleElementType(set[0].e) == CELL)
				{
					if (m->HighConnTag().isDefined(CELL))
					{
						Element::adj_type& cnodes = m->HighConn(set[0].e);
						nodes.insert(nodes.end(), cnodes.begin(), cnodes.end());
					}
					else
					{
						MarkerType mrk = m->CreatePrivateMarker();
						Element::adj_type& faces = m->LowConn(set[0].e);
						for (unsigned k = 0; k < faces.size(); ++k)
						{
							Element::adj_type& fedges = m->LowConn(faces[k]);
							for (unsigned q = 0; q < fedges.size(); ++q)
							{
								Element::adj_type& enodes = m->LowConn(fedges[q]);
								for (unsigned l = 0; l < enodes.size(); ++l) if (!m->GetPrivateMarker(enodes[l], mrk))
								{
									nodes.push_back(enodes[l]);
									m->SetPrivateMarker(enodes[l], mrk);
								}
							}
						}
						if (!nodes.empty()) m->RemPrivateMarkerArray(&nodes[0], (Storage::enumerator)nodes.size(), mrk);
						m->ReleasePrivateMarker(mrk);
					}
				}
				else if (GetHandleElementType(set[0].e) == FACE)
				{
					MarkerType mrk = m->CreatePrivateMarker();
					Element::adj_type& fedges = m->LowConn(set[0].e);
					for (unsigned q = 0; q < fedges.size(); ++q)
					{
						Element::adj_type& enodes = m->LowConn(fedges[q]);
						for (unsigned l = 0; l < enodes.size(); ++l) if (!m->GetPrivateMarker(enodes[l], mrk))
						{
							nodes.push_back(enodes[l]);
							m->SetPrivateMarker(enodes[l], mrk);
						}
					}
					if (!nodes.empty()) m->RemPrivateMarkerArray(&nodes[0], (Storage::enumerator)nodes.size(), mrk);
					m->ReleasePrivateMarker(mrk);
				}
				else
				{
					std::cout << "Unsupported element type in kd-tree" << std::endl;
					throw Impossible;
				}
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
	
	Cell SearchKDTree::SubSearchCell(const Storage::real p[3]) const
	{
		Cell ret = InvalidCell();
		if( size == 1 )
		{
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
				ret = children[0].SubSearchCell(p);
				if( !ret.isValid() )
					ret = children[1].SubSearchCell(p);
			}
		}
		return ret;
	}

	void PrintElement(Element c, std::ostream& sout)
	{
		sout << ElementTypeName(c.GetElementType()) << ":" << c.LocalID();
		sout << " " << Element::GeometricTypeName(c.GetGeometricType()) << std::endl;
		sout << " faces: " << c.nbAdjElements(FACE) << " edges: " << c.nbAdjElements(EDGE) << " nodes: " << c.nbAdjElements(NODE) << std::endl;
		if (c.GetElementType() == CELL)
		{
			Storage::real xc[3];
			c.Centroid(xc);
			sout << " closure: " << (c.getAsCell().Closure() ? "true" : "false");
			sout << " volume: " << c.getAsCell().Volume();
			sout << " x " << xc[0] << " " << xc[1] << " " << xc[2];
			sout << std::endl;

			ElementArray<Face> faces = c.getFaces();
			for (ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
			{
				Storage::real x[3], n[3];
				it->Centroid(x);
				it->UnitNormal(n);
				sout << "FACE:" << it->LocalID() << " nodes " << it->nbAdjElements(NODE);
				sout << " closure: " << (it->Closure() ? "true" : "false");
				sout << " flat: " << (it->Planarity() ? "true" : "false");
				sout << " bnd " << (it->Boundary() ? "true" : "false");
				sout << " ornt " << (it->CheckNormalOrientation() ? "true" : "false");
				sout << " x " << x[0] << " " << x[1] << " " << x[2] << " n " << n[0] << " " << n[1] << " " << n[2] << " " << " area " << it->Area() << std::endl;
			}
		}
		else if (c.GetElementType() == FACE)
		{
			Face it = c.getAsFace();
			Storage::real x[3], n[3];
			it->Centroid(x);
			it->UnitNormal(n);
			sout << " closure: " << (it->Closure() ? "true" : "false");
			sout << " flat: " << (it->Planarity() ? "true" : "false");
			sout << " bnd " << (it->Boundary() ? "true" : "false");
			sout << " ornt " << (it->CheckNormalOrientation() ? "true" : "false");
			sout << " x " << x[0] << " " << x[1] << " " << x[2] << " n " << n[0] << " " << n[1] << " " << n[2] << " " << " area " << it->Area() << std::endl;
		}
		ElementArray<Node> nodes = c.getNodes();
		for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); ++it)
			sout << "NODE:" << it->LocalID() << " x " << it->Coords()[0] << " " << it->Coords()[1] << " " << it->Coords()[2] << std::endl;
	}

	Cell SearchKDTree::SubSearchCellPrint(const Storage::real p[3], std::ostream & sout) const
	{
		Cell ret = InvalidCell();
		if (size == 1)
		{
			sout << "test cell " << GetHandleID(set[0].e) << std::endl;
			PrintElement(Cell(m, set[0].e), sout);
			if (m->HideMarker() && m->GetMarker(set[0].e, m->HideMarker()))
			{
				m->SwapModification(false);
				if (cell_point_print(Cell(m, set[0].e), p,sout))
					ret = Cell(m, set[0].e);
				m->SwapModification(false);
			}
			else
			{
				if (cell_point_print(Cell(m, set[0].e), p,sout))
					ret = Cell(m, set[0].e);
			}
		}
		else
		{
			assert(size > 1);
			if (bbox_point(p, bbox))
			{
				{
					sout << "point " << p[0] << " " << p[1] << " " << p[2] << " is in bbox ";
					sout << " x " << bbox[0] << ":" << bbox[1];
					sout << " y " << bbox[2] << ":" << bbox[3];
					sout << " z " << bbox[4] << ":" << bbox[5];
					sout << std::endl;
				}
				sout << "try left child" << std::endl;
				ret = children[0].SubSearchCellPrint(p, sout);
				sout << "ret " << (ret.isValid() ? ret.LocalID() : -1) << std::endl;
				if (!ret.isValid())
				{
					sout << "try right child" << std::endl;
					ret = children[1].SubSearchCellPrint(p, sout);
					sout << "ret " << (ret.isValid() ? ret.LocalID() : -1) << std::endl;
				}
			}
			else 
			{
				sout << "point " << p[0] << " " << p[1] << " " << p[2] << " is not in bbox ";
				sout << " x " << bbox[0] << ":" << bbox[1];
				sout << " y " << bbox[2] << ":" << bbox[3];
				sout << " z " << bbox[4] << ":" << bbox[5];
				sout << std::endl;
				bbox_point_print(p, bbox, sout);
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

	inline int SearchKDTree::sphere_tri(const Storage::real tri[3][3], const Storage::real p[3], Storage::real r) const
	{
		return triangle_distance(tri[0], tri[1], tri[2], p) <= r ? 1 : 0;
	}
	
	inline int SearchKDTree::segment_tri(const Storage::real tri[3][3], const Storage::real p1[3], const Storage::real p2[3]) const
	{
		const Storage::real eps = 1.0e-7;
		Storage::real a[3],b[3],c[3],n[3], ray[3], d, k, m, l;
		Storage::real dot00,dot01,dot02,dot11,dot12,invdenom, uq,vq;
		ray[0] = p2[0]-p1[0];
		ray[1] = p2[1]-p1[1];
		ray[2] = p2[2]-p1[2];
		/*
		l = sqrt(dotproduct(ray,ray));
		if (l)
		{
			ray[0] /= l;
			ray[1] /= l;
			ray[2] /= l;
		}
		*/
		a[0] = tri[0][0] - tri[2][0];
		a[1] = tri[0][1] - tri[2][1];
		a[2] = tri[0][2] - tri[2][2];
		b[0] = tri[1][0] - tri[2][0];
		b[1] = tri[1][1] - tri[2][1];
		b[2] = tri[1][2] - tri[2][2];
		crossproduct(a,b,n);
		l = sqrt(dotproduct(n, n));
		if (l)
		{
			n[0] /= l;
			n[1] /= l;
			n[2] /= l;
		}
		d = -dotproduct(n,tri[2]);
		m =  dotproduct(n,ray);
		if( fabs(m) < 1.0e-25 )
			return 0;
		k = -(d + dotproduct(n,p1))/m;
		if (k < -eps || k > 1.0 + eps)
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

	inline int SearchKDTree::segment_tri_print(const Storage::real tri[3][3], const Storage::real p1[3], const Storage::real p2[3], std::ostream& sout) const
	{
		const Storage::real eps = 1.0e-7;
		Storage::real a[3], b[3], c[3], n[3], ray[3], d, k, m, l;
		Storage::real dot00, dot01, dot02, dot11, dot12, invdenom, uq, vq;
		ray[0] = p2[0] - p1[0];
		ray[1] = p2[1] - p1[1];
		ray[2] = p2[2] - p1[2];
		/*
		l = sqrt(dotproduct(ray,ray));
		if (l)
		{
			ray[0] /= l;
			ray[1] /= l;
			ray[2] /= l;
		}
		*/
		sout << "r " << ray[0] << " " << ray[1] << " " << ray[2];
		a[0] = tri[0][0] - tri[2][0];
		a[1] = tri[0][1] - tri[2][1];
		a[2] = tri[0][2] - tri[2][2];
		b[0] = tri[1][0] - tri[2][0];
		b[1] = tri[1][1] - tri[2][1];
		b[2] = tri[1][2] - tri[2][2];
		crossproduct(a, b, n);
		l = sqrt(dotproduct(n, n));
		if (l)
		{
			n[0] /= l;
			n[1] /= l;
			n[2] /= l;
		}
		sout << "n " << n[0] << " " << n[1] << " " << n[2];
		d = -dotproduct(n, tri[2]);
		m = dotproduct(n, ray);
		sout << " d " << d << " m " << m;
		if (fabs(m) < 1.0e-25)
		{
			sout << std::endl;
			return 0;
		}
		k = -(d + dotproduct(n, p1)) / m;
		sout << " k " << k;
		if (k < -eps || k > 1.0 + eps)
		{
			sout << std::endl;
			return 0;
		}
		c[0] = p1[0] + k * ray[0] - tri[2][0];
		c[1] = p1[1] + k * ray[1] - tri[2][1];
		c[2] = p1[2] + k * ray[2] - tri[2][2];
		sout  << " c " << c[0] << " " << c[1] << " " << c[2];
		dot00 = dotproduct(a, a);
		dot01 = dotproduct(a, b);
		dot02 = dotproduct(a, c);
		dot11 = dotproduct(b, b);
		dot12 = dotproduct(b, c);
		invdenom = (dot00 * dot11 - dot01 * dot01);
		uq = (dot11 * dot02 - dot01 * dot12);
		vq = (dot00 * dot12 - dot01 * dot02);
		sout << " invdenom " << invdenom << " uq " << uq << " vq " << vq;
		if (fabs(invdenom) < 1.0e-25 && fabs(uq) > 0.0 && fabs(vq) > 0.0)
		{
			sout << std::endl;
			return 0;
		}
		uq = uq / invdenom;
		vq = vq / invdenom;
		sout << " coefs " << 1 - uq - vq << " " << uq << " " << vq;
		if (uq >= -eps && vq >= -eps && 1.0 - (uq + vq) >= -eps)
		{
			sout << " hit!" << std::endl;
			return 1;
		}
		sout << std::endl;
		return 0;
	}
	
	inline int SearchKDTree::segment_bbox(const Storage::real p1[3], const Storage::real p2[3]) const
	{
		Storage::real tnear = -1.0e20, tfar = 1.0e20;
		Storage::real t1,t2,c;
		for(int i = 0; i < 3; i++)
		{
			if( fabs(p2[i]-p1[i]) < 1.0e-15 )
			{
				if( p1[i] < bbox[i*2]-1.0e-3 || p1[i] > bbox[i*2+1]+1.0e-3 )
					return 0;
			}
			else
			{
				t1 = (bbox[i*2+0]-1.0e-3 - p1[i])/(p2[i]-p1[i]);
				t2 = (bbox[i*2+1]+1.0e-3 - p1[i])/(p2[i]-p1[i]);
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

	inline int SearchKDTree::sphere_bbox(const Storage::real p[3], Storage::real r) const
	{
		return bbox_sphere(p, r, bbox);
	}

	inline int SearchKDTree::ray_bbox(double pos[3], double ray[3], double closest) const
	{
		double tnear = -1.0e20, tfar = 1.0e20, t1, t2, c;
		for (int i = 0; i < 3; i++)
		{
			if (fabs(ray[i]) < 1.0e-15)
			{
				if (pos[i] < bbox[i * 2] || pos[i] > bbox[i * 2 + 1])
					return 0;
			}
			else
			{
				t1 = (bbox[i * 2 + 0] - pos[i]) / ray[i];
				t2 = (bbox[i * 2 + 1] - pos[i]) / ray[i];
				if (t1 > t2) { c = t1; t1 = t2; t2 = c; }
				if (t1 > tnear) tnear = t1;
				if (t2 < tfar) tfar = t2;
				if (tnear > closest) return 0;
				if (tnear > tfar) return 0;
				if (tfar < 0) return 0;
			}
		}
		return 1;
	}
	
	inline bool SearchKDTree::segment_face(const Element & f, const Storage::real p1[3], const Storage::real p2[3]) const
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

	inline bool SearchKDTree::segment_face_print(const Element& f, const Storage::real p1[3], const Storage::real p2[3], std::ostream & sout) const
	{
		Mesh* m = f->GetMeshLink();
		Storage::real tri[3][3] = { {0,0,0},{0,0,0},{0,0,0} }, nrm[3];
		f.getAsFace().UnitNormal(nrm);
		m->GetGeometricData(f->GetHandle(), CENTROID, tri[2]);
		sout << "segment_face FACE:" << f.LocalID();
		sout << " x " << tri[2][0] << " " << tri[2][1] << " " << tri[2][2];
		sout << " n " << nrm[0] << " " << nrm[1] << " " << nrm[2];
		sout << " bnd " << (f.getAsFace().Boundary() ? "true" : "false");
		sout << " ornt " << (f.getAsFace().CheckNormalOrientation() ? "true" : "false");
		sout << " nodes " << f.nbAdjElements(NODE);
		sout << std::endl;
		ElementArray<Node> nodes = f->getNodes();
		m->GetGeometricData(nodes[0]->GetHandle(), CENTROID, tri[1]);
		for (ElementArray<Node>::size_type k = 0; k < nodes.size(); ++k)
		{
			memcpy(tri[0], tri[1], sizeof(Storage::real) * 3);
			//m->GetGeometricData(nodes[k]->GetHandle(),CENTROID,tri[0]);
			m->GetGeometricData(nodes[(k + 1) % nodes.size()]->GetHandle(), CENTROID, tri[1]);
			if (segment_tri_print(tri, p1, p2,sout)) return true;
		}
		return false;
	}

	inline bool SearchKDTree::sphere_face(const Element& f, const Storage::real p[3], Storage::real r) const
	{
		Mesh* m = f->GetMeshLink();
		Storage::real tri[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
		m->GetGeometricData(f->GetHandle(), CENTROID, tri[2]);
		ElementArray<Node> nodes = f->getNodes();
		m->GetGeometricData(nodes[0]->GetHandle(), CENTROID, tri[1]);
		for (ElementArray<Node>::size_type k = 0; k < nodes.size(); ++k)
		{
			memcpy(tri[0], tri[1], sizeof(Storage::real) * 3);
			//m->GetGeometricData(nodes[k]->GetHandle(),CENTROID,tri[0]);
			m->GetGeometricData(nodes[(k + 1) % nodes.size()]->GetHandle(), CENTROID, tri[1]);
			if (sphere_tri(tri, p, r)) return true;
		}
		return false;
	}

	inline bool SearchKDTree::segment_cell(const Element & c, const Storage::real p1[3], const Storage::real p2[3]) const
	{
		ElementArray<Face> faces = c->getFaces();
		for(ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
			if( segment_face(it->self(),p1,p2) ) return true;
		return false;
	}

	inline bool SearchKDTree::sphere_cell(const Element& c, const Storage::real p[3], Storage::real r) const
	{
		ElementArray<Face> faces = c->getFaces();
		for (ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
			if (sphere_face(it->self(), p, r)) return true;
		return false;
	}
	
	void SearchKDTree::IntersectSegment(ElementArray<Cell> & cells, const Storage::real p1[3], const Storage::real p2[3]) const
	{
		ElementArray<Element> temp(m);
		MarkerType mrk = m->CreatePrivateMarker();
		sub_intersect_segment(temp, mrk, p1, p2);
		m->RemPrivateMarkerArray(temp.data(), static_cast<Storage::enumerator>(temp.size()), mrk);
		for (ElementArray<Element>::iterator it = temp.begin(); it != temp.end(); ++it)
		{
			if (it->GetElementType() == CELL && !it->GetPrivateMarker(mrk))
			{
				cells.push_back(it->getAsCell());
				it->SetPrivateMarker(mrk);
			}
			else if (it->GetElementType() == FACE)
			{
				ElementArray<Cell> f_cells = it->getCells();
				for (ElementArray<Cell>::iterator kt = f_cells.begin(); kt != f_cells.end(); ++kt) if (!kt->GetPrivateMarker(mrk))
				{
					cells.push_back(kt->self());
					kt->SetPrivateMarker(mrk);
				}
			}
		}
		m->RemPrivateMarkerArray(cells.data(), static_cast<Storage::enumerator>(cells.size()), mrk);
		m->ReleasePrivateMarker(mrk);
	}
	
	void SearchKDTree::IntersectSegment(ElementArray<Face>& faces, const Storage::real p1[3], const Storage::real p2[3]) const
	{
		ElementArray<Element> temp(m);
		MarkerType mrk = m->CreatePrivateMarker();
		sub_intersect_segment(temp, mrk, p1, p2);
		m->RemPrivateMarkerArray(temp.data(), static_cast<Storage::enumerator>(temp.size()), mrk);
		for (ElementArray<Element>::iterator it = temp.begin(); it != temp.end(); ++it)
		{
			if (it->GetElementType() == CELL)
			{
				ElementArray<Face> f = it->getFaces();
				for (ElementArray<Face>::iterator jt = f.begin(); jt != f.end(); ++jt) if( !jt->GetPrivateMarker(mrk) )
				{
					faces.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
			}
			else if (it->GetElementType() == FACE)
			{
				faces.push_back(it->getAsFace());
				it->SetPrivateMarker(mrk);
			}
		}
		m->RemPrivateMarkerArray(faces.data(), static_cast<Storage::enumerator>(faces.size()), mrk);
		m->ReleasePrivateMarker(mrk);
	}

	void SearchKDTree::IntersectSegmentPrint(ElementArray<Face>& faces, const Storage::real p1[3], const Storage::real p2[3], std::ostream & sout) const
	{
		ElementArray<Element> temp(m);
		MarkerType mrk = m->CreatePrivateMarker();
		sub_intersect_segment_print(temp, mrk, p1, p2, sout);
		sout << "found elements: " << temp.size() << std::endl;
		m->RemPrivateMarkerArray(temp.data(), static_cast<Storage::enumerator>(temp.size()), mrk);
		for (ElementArray<Element>::iterator it = temp.begin(); it != temp.end(); ++it)
		{
			sout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << " ";
			if (it->GetElementType() == CELL)
			{
				ElementArray<Face> f = it->getFaces();
				for (ElementArray<Face>::iterator jt = f.begin(); jt != f.end(); ++jt) if (!jt->GetPrivateMarker(mrk))
				{
					faces.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
			}
			else if (it->GetElementType() == FACE)
			{
				faces.push_back(it->getAsFace());
				it->SetPrivateMarker(mrk);
			}
		}
		sout << std::endl;
		sout << "Output faces: " << faces.size() << std::endl;
		for (ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
			sout << ElementTypeName(it->GetElementType()) << ":" << it->LocalID() << " ";
		sout << std::endl;
		m->RemPrivateMarkerArray(faces.data(), static_cast<Storage::enumerator>(faces.size()), mrk);
		m->ReleasePrivateMarker(mrk);
	}

	void SearchKDTree::IntersectSphere(ElementArray<Cell>& cells, const Storage::real p[3], Storage::real r) const
	{
		ElementArray<Element> temp(m);
		MarkerType mrk = m->CreatePrivateMarker();
		sub_intersect_sphere(temp, mrk, p, r);
		m->RemPrivateMarkerArray(temp.data(), static_cast<Storage::enumerator>(temp.size()), mrk);
		for (ElementArray<Element>::iterator it = temp.begin(); it != temp.end(); ++it)
		{
			if (it->GetElementType() == CELL && !it->GetPrivateMarker(mrk))
			{
				cells.push_back(it->getAsCell());
				it->SetPrivateMarker(mrk);
			}
			else if (it->GetElementType() == FACE)
			{
				ElementArray<Cell> f_cells = it->getCells();
				for (ElementArray<Cell>::iterator kt = f_cells.begin(); kt != f_cells.end(); ++kt) if (!kt->GetPrivateMarker(mrk))
				{
					cells.push_back(kt->self());
					kt->SetPrivateMarker(mrk);
				}
			}
		}
		m->RemPrivateMarkerArray(cells.data(), static_cast<Storage::enumerator>(cells.size()), mrk);
		m->ReleasePrivateMarker(mrk);
	}

	void SearchKDTree::IntersectSphere(ElementArray<Face>& faces, const Storage::real p[3], Storage::real r) const
	{
		ElementArray<Element> temp(m);
		MarkerType mrk = m->CreatePrivateMarker();
		sub_intersect_sphere(temp, mrk, p, r);
		m->RemPrivateMarkerArray(temp.data(), static_cast<Storage::enumerator>(temp.size()), mrk);
		for (ElementArray<Element>::iterator it = temp.begin(); it != temp.end(); ++it)
		{
			if (it->GetElementType() == CELL)
			{
				ElementArray<Face> f = it->getFaces();
				for (ElementArray<Face>::iterator jt = f.begin(); jt != f.end(); ++jt) if (!jt->GetPrivateMarker(mrk))
				{
					faces.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
			}
			else if (it->GetElementType() == FACE)
			{
				faces.push_back(it->getAsFace());
				it->SetPrivateMarker(mrk);
			}
		}
		m->RemPrivateMarkerArray(faces.data(), static_cast<Storage::enumerator>(faces.size()), mrk);
		m->ReleasePrivateMarker(mrk);
	}
	
	void SearchKDTree::sub_intersect_segment(ElementArray<Element> & hits, MarkerType mrk, const Storage::real p1[3], const Storage::real p2[3]) const
	{
		if( size == 1 )
		{
			if( !m->GetPrivateMarker(set[0].e,mrk) )
			{
				Storage::integer edim = Element::GetGeometricDimension(m->GetGeometricType(set[0].e));
				if( edim == 2 )
				{
					if( segment_face(Element(m,set[0].e),p1,p2) )
					{
						hits.push_back(set[0].e);
						m->SetPrivateMarker(set[0].e,mrk);
					}
				}
				else if( edim == 3 )
				{
					if( segment_cell(Element(m,set[0].e),p1,p2) )
					{
						hits.push_back(set[0].e);
						m->SetPrivateMarker(set[0].e,mrk);
					}
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " kd-tree structure is not implemented to intersect edges with segments" << std::endl;
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

	void SearchKDTree::sub_intersect_segment_print(ElementArray<Element>& hits, MarkerType mrk, const Storage::real p1[3], const Storage::real p2[3], std::ostream & sout) const
	{
		if (size == 1)
		{
			if (!m->GetPrivateMarker(set[0].e, mrk))
			{
				Storage::integer edim = Element::GetGeometricDimension(m->GetGeometricType(set[0].e));
				PrintElement(Element(m, set[0].e), sout);
				if (edim == 2)
				{
					if (segment_face_print(Element(m, set[0].e), p1, p2, sout))
					{
						hits.push_back(set[0].e);
						m->SetPrivateMarker(set[0].e, mrk);
					}
				}
				else if (edim == 3)
				{
					if (segment_cell(Element(m, set[0].e), p1, p2))
					{
						hits.push_back(set[0].e);
						m->SetPrivateMarker(set[0].e, mrk);
					}
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " kd-tree structure is not implemented to intersect edges with segments" << std::endl;
					exit(-1);
				}
			}
		}
		else
		{
			assert(size > 1);
			if (segment_bbox(p1, p2))
			{
				{
					sout << "segment ";
					sout << p1[0] << " " << p1[1] << " " << p1[2];
					sout << " <-> ";
					sout << p2[0] << " " << p2[1] << " " << p2[2];
					sout << " is in bbox ";
					sout << " x " << bbox[0] << ":" << bbox[1];
					sout << " y " << bbox[2] << ":" << bbox[3];
					sout << " z " << bbox[4] << ":" << bbox[5];
					sout << std::endl;
				}
				children[0].sub_intersect_segment_print(hits, mrk, p1, p2, sout);
				children[1].sub_intersect_segment_print(hits, mrk, p1, p2, sout);
			}
			else
			{
				sout << "segment ";
				sout << p1[0] << " " << p1[1] << " " << p1[2];
				sout << " <-> ";
				sout << p2[0] << " " << p2[1] << " " << p2[2];
				sout << " is not in bbox ";
				sout << " x " << bbox[0] << ":" << bbox[1];
				sout << " y " << bbox[2] << ":" << bbox[3];
				sout << " z " << bbox[4] << ":" << bbox[5];
				sout << std::endl;
			}
		}
	}

	void SearchKDTree::sub_intersect_sphere(ElementArray<Element>& hits, MarkerType mrk, const Storage::real p[3], Storage::real r) const
	{
		if (size == 1)
		{
			if (!m->GetPrivateMarker(set[0].e, mrk))
			{
				Storage::integer edim = Element::GetGeometricDimension(m->GetGeometricType(set[0].e));
				if (edim == 2)
				{
					if (sphere_face(Element(m, set[0].e), p, r))
					{
						hits.push_back(set[0].e);
						m->SetPrivateMarker(set[0].e, mrk);
					}
				}
				else if (edim == 3)
				{
					if (sphere_cell(Element(m, set[0].e), p, r))
					{
						hits.push_back(set[0].e);
						m->SetPrivateMarker(set[0].e, mrk);
					}
				}
				else
				{
					std::cout << __FILE__ << ":" << __LINE__ << " kd-tree structure is not implemented to intersect edges with spheres" << std::endl;
					exit(-1);
				}
			}
		}
		else
		{
			assert(size > 1);
			if (sphere_bbox(p, r))
			{
				children[0].sub_intersect_sphere(hits, mrk, p, r);
				children[1].sub_intersect_sphere(hits, mrk, p, r);
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
	
	Cell SearchKDTree::SearchCell(const Storage::real * point) const
	{
		return SubSearchCell(point);
	}

	Cell SearchKDTree::SearchCellPrint(const Storage::real* point, std::ostream & sout) const
	{
		return SubSearchCellPrint(point, sout);
	}
}
#endif

