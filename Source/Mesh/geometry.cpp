#include "inmost.h"
#if defined(USE_MESH)
#include <deque>
using namespace std;

const std::string normal_name = "PROTECTED_GEOM_UTIL_NORMAL";
const std::string measure_name = "PROTECTED_GEOM_UTIL_MEASURE";
const std::string centroid_name = "PROTECTED_GEOM_UTIL_CENTROID";
const std::string barycenter_name = "PROTECTED_GEOM_UTIL_BARYCENTER";

namespace INMOST
{
	typedef struct orient_face_t
	{
		Edge bridge;
		Node first;
		Face face;
		orient_face_t(Edge _bridge, Node _first, Face _face)
		:bridge(_bridge),first(_first),face(_face)
		{
		}
	} orient_face;

	__INLINE static void vec_diff(const Storage::real * vecin1, const Storage::real * vecin2, Storage::real * vecout, unsigned int size)
	{
		for(unsigned int i = 0; i < size; i++)
			vecout[i] = vecin1[i] - vecin2[i];
	}
	
	__INLINE static void vec_diff(const Storage::real_array & vecin1, const Storage::real * vecin2, Storage::real * vecout, unsigned int size)
	{
		for(unsigned int i = 0; i < size; i++)
			vecout[i] = vecin1[i] - vecin2[i];
	}
	
	__INLINE static void vec_diff(const Storage::real_array & vecin1,const Storage::real_array & vecin2, Storage::real * vecout, unsigned int size)
	{
		for(unsigned int i = 0; i < size; i++)
			vecout[i] = vecin1[i] - vecin2[i];
	}
	
	
	__INLINE static void vec_cross_product(const Storage::real * vecin1, const Storage::real * vecin2, Storage::real * vecout)
	{
		Storage::real temp[3];
		temp[0] = vecin1[1]*vecin2[2] - vecin1[2]*vecin2[1];
		temp[1] = vecin1[2]*vecin2[0] - vecin1[0]*vecin2[2];
		temp[2] = vecin1[0]*vecin2[1] - vecin1[1]*vecin2[0];
		vecout[0] = temp[0];
		vecout[1] = temp[1];
		vecout[2] = temp[2];
	}
	
	__INLINE static Storage::real vec_dot_product(const Storage::real * vecin1,const Storage::real * vecin2, unsigned int size)
	{
		Storage::real ret = 0;
		for(unsigned int i = 0; i < size; i++)
			ret += vecin1[i]*vecin2[i];
		return ret;
	}
	
	__INLINE static Storage::real vec_len2(const Storage::real * vecin, unsigned int size)
	{
		return vec_dot_product(vecin,vecin,size);
	}
	
	__INLINE static Storage::real vec_len(const Storage::real * vecin, unsigned int size)
	{
		return ::sqrt(vec_len2(vecin,size));
	}
	
	__INLINE static Storage::real det3d(Storage::real a, Storage::real b, Storage::real c,
	                         Storage::real d, Storage::real e, Storage::real f,
	                         Storage::real g, Storage::real h, Storage::real i ) 
	{
		return a*e*i - c*e*g + b*f*g - a*f*h + c*d*h - b*d*i;
	}
	
	__INLINE static Storage::real det3v(const Storage::real * x,const Storage::real * y,const Storage::real * z) 
	{
		return det3d(x[0], x[1], x[2],  y[0], y[1], y[2],  z[0], z[1], z[2]);
	}
	
	__INLINE static Storage::real det4v(const Storage::real * w, const Storage::real * x, const Storage::real * y, const Storage::real * z) 
	{
		return det3d(x[0]-w[0], x[1]-w[1], x[2]-w[2],  y[0]-w[0], y[1]-w[1], y[2]-w[2],  z[0]-w[0], z[1]-w[1], z[2]-w[2]);
	}	
	
	
	__INLINE static Storage::real vec_normalize(Storage::real * vecin, unsigned int size)
	{
		Storage::real len = 0;
		for(unsigned int i = 0; i < size; i++)
			len += vecin[i]*vecin[i];
		len = ::sqrt(len);
		for(unsigned int i = 0; i < size; i++)
			vecin[i] /= len;
		return len;
	}

	__INLINE static Storage::real triarea3d(const Storage::real p1[3], const Storage::real p2[3], const Storage::real p3[3])
	{
		Storage::real v12[3], v13[3], v23[3], l12, l13, l23, halfperim;
		vec_diff(p2,p1,v12,3);
		vec_diff(p3,p1,v13,3);
		vec_diff(p3,p2,v23,3);
		l12 = vec_len(v12,3);
		l13 = vec_len(v13,3);
		l23 = vec_len(v23,3);
		halfperim = 0.5*(l12+l13+l23);
		return sqrt(halfperim*(halfperim-l12)*(halfperim-l13)*(halfperim-l23));
	}
	
	
	ElementArray<Cell> Cell::NeighbouringCells() const
	{
		ElementArray<Cell> ret(GetMeshLink());
		ElementArray<Face> faces = getFaces();
		for(ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); f++)
		{
			Cell c = Neighbour(f->self());
			if( c.isValid() ) ret.push_back(c);
		}
		return ret;
	}
	
	
	Cell Cell::Neighbour(Face f) const
	{
		Cell b = f->BackCell();
		if( b == self() )
			return f->FrontCell();
		return b;
	}

	bool Face::Inside(const Storage::real * point) const
	{
		Mesh * mesh = GetMeshLink();
		real eps = mesh->GetEpsilon();
		integer dim = mesh->GetDimensions();
		integer mdim = GetElementDimension();
		
		if(mdim < 2)
		{
			// check whether point lies on edge
			real v1[3], v2[3], v12[3], r1[3], r2[3], nrm[3], area, h, len;
			ElementArray<Node> nodes = getNodes();
			assert(nodes.size() == 2);// 
			nodes[0].Centroid(v1);
			nodes[1].Centroid(v2);
			vec_diff(point,v1,r1,dim);
			vec_diff(point,v2,r2,dim);
			vec_diff(v2,v1,v12,dim);
			len = vec_len(v12,dim);
			assert(len > 1e-20);
			vec_cross_product(r1,v12,nrm);
			area = vec_len(nrm,dim);
			h = area / len;
			if( h > eps )	return false; // point does not lie on line
			return vec_dot_product(r1,r2,dim) <= 0.0; // point is between edge nodes
		}

		real nrm[3], cnt[3], v[3];
		UnitNormal(nrm);
		Centroid(cnt);
		vec_diff(cnt,point,v,dim);
		real d = vec_dot_product(nrm,v,dim);
		if(fabs(d) > eps)	return false;	// point is too far from the face plane

		// 2d algorithm from Cell::Inside 
		real data[9][3];
		for(int k = 0; k < 9; k++) for(int j = 0; j < 3; j++)	data[k][j] = 0;
		ElementArray<Node> nodes = getNodes();
		for(int i = 0; i < static_cast<int>(nodes.size()); i++)
		{
			int j = (i+1)%nodes.size();
			nodes[i].Centroid(data[1]);
			nodes[j].Centroid(data[2]);
			vec_diff(point,data[0],data[3],dim);
			vec_diff(point,data[1],data[4],dim);
			vec_diff(point,data[2],data[5],dim);
			vec_cross_product(data[3],data[4],data[6]);
			vec_cross_product(data[4],data[5],data[7]);
			vec_cross_product(data[5],data[3],data[8]);
			
			if( vec_dot_product(data[6],data[7],dim) >= 0 &&
				vec_dot_product(data[7],data[8],dim) >= 0 &&
				vec_dot_product(data[8],data[6],dim) >= 0 )
				return true; //inside one of the triangles
		}
		return false;
	}
	
	bool Cell::Inside(const Storage::real * point) const//check for 2d case
	{
		Mesh * mesh = GetMeshLink();
		integer dim = GetElementDimension();
		if( dim == 3 )
		{
			/*
			std::map<HandleType,real,16> hits;
			real ray[3];
			ray[0] = rand()/(real)RAND_MAX;
			ray[1] = rand()/(real)RAND_MAX;
			ray[2] = rand()/(real)RAND_MAX;
			CastRay(point,ray,hits);
			if( hits.size()%2 == 0 ) return false;
			return true;
			 */
			assert(mesh->GetDimensions() == 3);
			integer vp = 0;
			integer vm = 0;
			integer vz = 0;
			real eps = mesh->GetEpsilon();
			real c,d, fcnt[3];
			real_array v1,v2;
			ElementArray<Face> data = getFaces();
			Face cur = data[0];
			MarkerType mrk = mesh->CreatePrivateMarker();
			MarkerType rev = mesh->CreatePrivateMarker(); //reverse orientation
			data.SetPrivateMarker(mrk); //0-th face orientation is default
			cur->RemPrivateMarker(mrk);
			Node n1,n2; //to retrive edge
			bool reverse = false; //reverse orientation in considered face
			std::deque< orient_face > stack; //edge and first node and face for visiting
			//todo: can do faster by retriving edges and going over their nodes
			//should not use FindSharedAdjacency
			ElementArray<Edge> edges = cur->getEdges();
			do
			{
				//figure out starting node order
				if( edges[0]->getBeg() == edges[1]->getBeg() ||
				   edges[0]->getBeg() == edges[1]->getEnd() )
				{
					n1 = edges[0]->getEnd();
					n2 = edges[0]->getBeg();
				}
				else
				{
					n1 = edges[0]->getBeg();
					n2 = edges[0]->getEnd();
				}
				//schedule unvisited adjacent faces
				for(unsigned j = 0; j < edges.size(); j++)
				{
					//schedule face adjacent to considered edge
					ElementArray<Face> adjacent = edges[j]->getFaces(mrk);
					assert(adjacent.size() <= 1);
					if( !adjacent.empty() && adjacent[0].GetPrivateMarker(mrk))
					{
						adjacent.RemPrivateMarker(mrk);
						stack.push_back(orient_face(edges[j],reverse ? n2 : n1,adjacent[0]));
					}
					//update edge nodes
					n1 = n2; //current end is new begin
					//find new end
					if( n2 == edges[(j+1)%edges.size()]->getBeg() )
						n2 = edges[(j+1)%edges.size()]->getEnd();
					else
						n2 = edges[(j+1)%edges.size()]->getBeg();
				}
				if( stack.empty() ) break;
				//get entry from stack
				orient_face r = stack.front();
				//remove face from stack
				stack.pop_front();
				//retrive edges for new face
				edges = r.face->getEdges();
				reverse = false;
				//figure out starting node order
				if( edges[0]->getBeg() == edges[1]->getBeg() ||
				   edges[0]->getBeg() == edges[1]->getEnd() )
				{
					n1 = edges[0]->getEnd();
					n2 = edges[0]->getBeg();
				}
				else
				{
					n1 = edges[0]->getBeg();
					n2 = edges[0]->getEnd();
				}
				//find out common edge orientation
				for(unsigned j = 0; j < edges.size(); j++)
				{
					if( edges[j] == r.bridge ) //found the edge
					{
						//reverse ordering on this face
						if( r.first == n1 )
						{
							r.face->SetPrivateMarker(rev);
							reverse = true;
						}
						break;
					}
					//update edge nodes
					n1 = n2; //current end is new begin
					//find new end
					if( n2 == edges[(j+1)%edges.size()]->getBeg() )
						n2 = edges[(j+1)%edges.size()]->getEnd();
					else
						n2 = edges[(j+1)%edges.size()]->getBeg();
				}
			} while(true);
			data.RemPrivateMarker(mrk);
			mesh->ReleasePrivateMarker(mrk);
			for(ElementArray<Face>::size_type f = 0; f < data.size(); f++)
			{
				d = 0.0;
				data[f]->Centroid(fcnt);
				ElementArray<Node> nodes = data[f]->getNodes();
				v1 = nodes[0].Coords();
				for (ElementArray<Node>::size_type i=0; i<nodes.size(); i++)
				{
					v2 = nodes[(i+1)%nodes.size()].Coords();
					d += c = det4v(point, fcnt, v1.data(), v2.data());
					v1.swap(v2);
				}
				//if(!data[f]->FaceOrientedOutside(self()))
				if( data[f]->GetPrivateMarker(rev) )
					c = -1.0;
				else
					c = 1.0;
				if(c*d > eps)
					vp++;
				else if(c*d < -eps)
					vm++;
				else
					vz++;
			}
			data.RemPrivateMarker(rev);
			mesh->ReleasePrivateMarker(rev);
			if(vp*vm > 0) return false;
			else if( vz == 0 ) return true;
			else return true;
		}
		else
		{
			int mdim = mesh->GetDimensions();
			assert(mdim <= 3);
			Storage::real data[9][3];
			if( mdim < 3 )
			{
				memset(data,0,sizeof(Storage::real)*9*3);
			}
			Centroid(data[0]);
			ElementArray<Node> nodes = getNodes();
			for(int i = 0; i < static_cast<int>(nodes.size()); i++)
			{
				int j = (i+1)%nodes.size();
				nodes[i].Centroid(data[1]);
				nodes[j].Centroid(data[2]);
				vec_diff(point,data[0],data[3],mdim);
				vec_diff(point,data[1],data[4],mdim);
				vec_diff(point,data[2],data[5],mdim);
				vec_cross_product(data[3],data[4],data[6]);
				vec_cross_product(data[4],data[5],data[7]);
				vec_cross_product(data[5],data[3],data[8]);
				
				if( vec_dot_product(data[6],data[7],mdim) >= 0 &&
				   vec_dot_product(data[7],data[8],mdim) >= 0 &&
				   vec_dot_product(data[6],data[8],mdim) >= 0 )
					return true; //inside one of the triangles
			}
			return false;
		}
	}
	
	
	void Face::UnitNormal(real * nrm) const
	{
		Mesh * m = GetMeshLink();
		m->GetGeometricData(GetHandle(),NORMAL,nrm); 
		integer dim = m->GetDimensions();
		real    l   = ::sqrt(vec_dot_product(nrm,nrm,dim)); 
		if(::fabs(l) > m->GetEpsilon()) 
		{
			for(integer i = 0; i < dim; i++) 
				nrm[i] /= l; 
		}
	}
	
	void Face::OrientedNormal(Cell c, Storage::real * nrm) const
	{
		Normal(nrm); 
		if( !FaceOrientedOutside(c) )
		{
			integer dim = GetMeshLink()->GetDimensions();
			for(integer i = 0; i < dim; i++) 
				nrm[i] = -nrm[i];
		}
	}
	
	
	void Face::OrientedUnitNormal(Cell c, Storage::real * nrm) const
	{
		UnitNormal(nrm); 
		if( !FaceOrientedOutside(c) ) 
		{
			integer dim = GetMeshLink()->GetDimensions();
			for(integer i = 0; i < dim; i++) 
				nrm[i] = -nrm[i];
		}
	}
	
	
	
	bool Mesh::TestClosure(const HandleType * elements, integer num) const
	{
		integer i;
		std::map<HandleType,int> e_visit;
		std::map<HandleType,int>::iterator it;
		if( !HideMarker() )
		{
			for(i = 0; i < num; i++)
			{
				Element::adj_type const & lc = LowConn(elements[i]);
				for(Element::adj_type::size_type jt = 0; jt < lc.size(); jt++)
					e_visit[lc[jt]]++;
			}
		}
		else
		{
			for(i = 0; i < num; i++) if( !GetMarker(elements[i],HideMarker()) )
			{
				Element::adj_type const & lc = LowConn(elements[i]);
				for(Element::adj_type::size_type jt = 0; jt < lc.size(); jt++) 
				{
					if( !GetMarker(lc[jt],HideMarker()) ) 
						e_visit[lc[jt]]++;
				}
			}
		}
		for(it = e_visit.begin(); it != e_visit.end(); it++)
			if( it->second != 2 ) return false;
		return true;
	}
	
	Element::GeometricType Mesh::ComputeGeometricType(ElementType etype, const HandleType * lc, INMOST_DATA_ENUM_TYPE size)
	{
		Element::GeometricType ret = Element::Unset;
		INMOST_DATA_ENUM_TYPE s = 0;
		if( isMeshModified() )
			s = Mesh::Count(lc,size,HideMarker());
		else s = size;
		int dmax = -1, dmin = 4;
		if( s == 0 && etype != NODE) return ret;
		switch(etype)
		{
			case NODE: ret = Element::Vertex; break;
			case EDGE:
				if( s == 1 )
					ret = Element::Vertex;
				else if( s == 2 )
					ret = Element::Line;
				break;
			case FACE:
				for (INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k) if( !Hidden(lc[k]) )
				{
					int d = Element::GetGeometricDimension(GetGeometricType(lc[k]));
					if (dmax < d) dmax = d;
					if (dmin > d) dmin = d;
				}
				if (dmax != dmin)
				{
					ret = Element::MultiLine;
				}
				else if( dmax == 0 )
				{ 
					ret = Element::Line;
				}
				else
				{
					if( !GetTopologyCheck(NEED_TEST_CLOSURE) || TestClosure(lc,size) )
					{
						if( s == 3 )
							ret = Element::Tri;
						else if( s == 4 )
							ret = Element::Quad;
						else
							ret = Element::Polygon;
					}
					else ret = Element::MultiLine;
				}
				break;
			case CELL:
				for (INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k) if( !Hidden(lc[k]) )
				{
					int d = Element::GetGeometricDimension(GetGeometricType(lc[k]));
					if (dmax < d) dmax = d;
					if (dmin > d) dmin = d;
				}
				if (dmax != dmin)
				{
					ret = Element::MultiPolygon;
				}
				else if(  dmax == 1 )
				{
					if( !GetTopologyCheck(NEED_TEST_CLOSURE) || TestClosure(lc,size) )
					{
						if( s == 3 )
							ret = Element::Tri;
						else if( s == 4 )
							ret = Element::Quad;
						else
							ret = Element::Polygon;
					}
					else ret = Element::MultiLine;
				}
				else 
				{
					if( !GetTopologyCheck(NEED_TEST_CLOSURE) ||  TestClosure(lc,size) )
					{
						// for simple cells there is no more then one common edge
						// otherwise the cell should be treated as polyhedron
						bool check = true;
						MarkerType common = CreatePrivateMarker();
						for(INMOST_DATA_ENUM_TYPE i = 0; i < s; ++i)
						{
							const Element::adj_type & ilc = LowConn(lc[i]); // access edges of the i-th face
							for(INMOST_DATA_ENUM_TYPE r = 0; r < ilc.size(); ++r ) SetPrivateMarker(ilc[r],common);
							for(INMOST_DATA_ENUM_TYPE j = i+1; j < s; ++j )
							{
								const Element::adj_type & jlc = LowConn(lc[j]); // access edges of the j-th face
								int cnt = 0;
								for(INMOST_DATA_ENUM_TYPE r = 0; r < jlc.size(); ++r )
									if( GetPrivateMarker(jlc[r],common) ) cnt++;
								if( cnt > 1 ) check = false;
							}
							for(INMOST_DATA_ENUM_TYPE r = 0; r < ilc.size(); ++r ) RemPrivateMarker(ilc[r],common);
						}
						ReleasePrivateMarker(common);
						if( check )
						{
							//test c_faces closure, if no closure, set as MultiPolygon
							INMOST_DATA_ENUM_TYPE quads = 0,tris = 0,i;
							for(i = 0; i < size; i++) if( !Hidden(lc[i]) )
							{
								if( GetGeometricType(lc[i]) == Element::Tri )
									tris++;
								else if( GetGeometricType(lc[i]) == Element::Quad )
									quads++;
							}
							if( tris == 4 && s == 4 )
								ret = Element::Tet;
							else if( quads == 6 && s == 6 )
								ret = Element::Hex;
							else if( tris == 4 && quads == 1 && s == tris+quads)
								ret = Element::Pyramid;
							else if( quads == 3 && tris == 2 && s == tris+quads)
								ret = Element::Prism;
							else
								ret = Element::Polyhedron;
						}
						else ret = Element::Polyhedron;
					}
					else ret = Element::MultiPolygon;
				}
				break;
			case ESET: ret = Element::Set; break;
		}
		return ret;
	}

	Storage::real Edge::Length() const 
	{
		Storage::real ret; 
		GetMeshLink()->GetGeometricData(GetHandle(),MEASURE,&ret); 
		return ret;
	}

	Storage::real Face::Area() const 
	{
		real ret; 
		GetMeshLink()->GetGeometricData(GetHandle(),MEASURE,&ret); 
		return ret;
	}
	void Face::Normal(real * nrm) const 
	{
		GetMeshLink()->GetGeometricData(GetHandle(),NORMAL,nrm);
	}

	Storage::real Cell::Volume() const 
	{
		real ret; 
		GetMeshLink()->GetGeometricData(GetHandle(),MEASURE,&ret); 
		return ret;
	}

	void Element::ComputeGeometricType() const
	{
		GetMeshLink()->ComputeGeometricType(GetHandle());
	}
	
	void Mesh::ComputeGeometricType(HandleType h) 
	{
		SetGeometricType(h,Element::Unset);
		if( GetElementType() > NODE )
		{
			Element::adj_type const & lc = LowConn(h);
			if( !lc.empty() )
				SetGeometricType(h,ComputeGeometricType(GetHandleElementType(h),lc.data(),static_cast<integer>(lc.size())));
		}
		else SetGeometricType(h,Element::Vertex);
	}
	
	void Mesh::RecomputeGeometricData(HandleType e)
	{
		if (GetMarker(e, UpdateGeometryMarker())) //prevent recursive entry
			RemMarker(e, UpdateGeometryMarker());
		//static std::map<Element *, int> numfixes;
		GeometricData d ;
		for(d = CENTROID; d <= NORMAL; d++) // first compute centroids and normals 
		{
			if( HaveGeometricData(d,GetHandleElementType(e)) ) //compute centroid first
			{
				Tag t = GetGeometricTag(d);
				Storage::real * a = static_cast<Storage::real *>(MGetDenseLink(e,t));
				HideGeometricData(d,GetHandleElementType(e));
				GetGeometricData(e,d,a);
				ShowGeometricData(d,GetHandleElementType(e));
			}
		}


		if( HaveGeometricData(ORIENTATION,FACE) )
		{
			if( GetHandleElementType(e) == CELL ) //then correct the normal
			{
				Element::adj_type & lc = LowConn(e); //faces
				for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)
					if( !GetMarker(*it,HideMarker()) )
					{
						//Element::adj_type & hc = HighConn(e);
						//if( Mesh::Count(&hc[0],hc.size(),HideMarker()) == 1 )
							Face(this,*it)->FixNormalOrientation();
					}
			}
			else if( GetHandleElementType(e) == FACE )
				Face(this,e)->FixNormalOrientation();
		}

		for(d = MEASURE; d <= BARYCENTER; d++) // compute the rest
		{
			if( HaveGeometricData(d,GetHandleElementType(e)) )
			{
				Tag t = GetGeometricTag(d);
				Storage::real * a = static_cast<Storage::real *>(MGetDenseLink(e,t));
				HideGeometricData(d,GetHandleElementType(e));
				GetGeometricData(e,d,a);
				ShowGeometricData(d,GetHandleElementType(e));
			}
		}
		
		if( GetHandleElementType(e) == CELL ) //then correct the normal
		{
			
		}
	}
	
	
	void Mesh::RemoveGeometricData(GeomParam table)
	{
		for(GeomParam::iterator it = table.begin(); it != table.end(); ++it)
		{
			if( it->first == MEASURE    ) 
			{
				if(measure_tag.isValid())    
					measure_tag    = DeleteTag(measure_tag   ,it->second);
				for(ElementType etype = EDGE; etype <= CELL; etype = etype << 1) if( etype & it->second) HideGeometricData(MEASURE,etype);
			}
			if( it->first == CENTROID   ) 
			{
				if(centroid_tag.isValid())   
					centroid_tag   = DeleteTag(centroid_tag  ,it->second);
				for(ElementType etype = EDGE; etype <= CELL; etype = etype << 1) if( etype & it->second) HideGeometricData(CENTROID,etype);
			}
			if( it->first == BARYCENTER ) 
			{
				if(barycenter_tag.isValid()) 
					barycenter_tag = DeleteTag(barycenter_tag,it->second);
				for(ElementType etype = EDGE; etype <= CELL; etype = etype << 1) if( etype & it->second) HideGeometricData(BARYCENTER,etype);
			}
			if( it->first == NORMAL     ) 
			{
				if(normal_tag.isValid())
					normal_tag     = DeleteTag(normal_tag    ,it->second);
				for(ElementType etype = FACE; etype <= CELL; etype = etype << 1) if( etype & it->second) HideGeometricData(NORMAL,etype);
			}
			if( it->first == ORIENTATION) 
				if( FACE & it->second) HideGeometricData(ORIENTATION,FACE);
		}
	}
	
	void Mesh::RestoreGeometricTags()
	{
		for(GeometricData gtype = MEASURE; gtype <= NORMAL; gtype++)
		{
			bool restore = false;
			for(ElementType etype = EDGE; etype <= CELL && !restore; etype = NextElementType(etype))
				if( HaveGeometricData(gtype,etype) )
					restore = true;
			if( restore )
			{
				switch(gtype)
				{
				case MEASURE:       measure_tag = GetTag(measure_name);    break;
				case CENTROID:     centroid_tag = GetTag(centroid_name);   break;
				case BARYCENTER: barycenter_tag = GetTag(barycenter_name); break;
				case NORMAL:         normal_tag = GetTag(normal_name);     break;
				}
			}
		}
	}

	void Mesh::RepairGeometricTags()
	{
		if( HaveTag(measure_name) )
		{
			measure_tag = GetTag(measure_name);
			for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
				if( measure_tag.isDefined(etype) && !HaveGeometricData(MEASURE,etype) )
					ShowGeometricData(MEASURE,etype);
		}
		if( HaveTag(centroid_name) )
		{
			centroid_tag = GetTag(centroid_name);
			for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
				if( centroid_tag.isDefined(etype) && !HaveGeometricData(CENTROID,etype) )
					ShowGeometricData(CENTROID,etype);
		}
		if( HaveTag(barycenter_name) )
		{
			barycenter_tag = GetTag(barycenter_name);
			for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
				if( barycenter_tag.isDefined(etype) && !HaveGeometricData(BARYCENTER,etype) )
					ShowGeometricData(BARYCENTER,etype);
		}
		if( HaveTag(normal_name) )
		{
			normal_tag = GetTag(normal_name);
			for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
				if( normal_tag.isDefined(etype) && !HaveGeometricData(NORMAL,etype) )
					ShowGeometricData(NORMAL,etype);
		}
	}
	
	void Mesh::PrepareGeometricData(GeomParam table)
	{
		for(GeomParam::iterator it = table.begin(); it != table.end(); ++it)
		{
			GeometricData types = it->first;
			ElementType mask = it->second;
			if( types == ORIENTATION )
			{
				//std::cout << "ORIENTATION" << std::endl;
				if( mask & FACE )
				{
					if( HideMarker() )
					{
						MarkerType hm = HideMarker();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer e = 0; e < FaceLastLocalID(); ++e) 
						{
							if( isValidElement(FACE,e) )
							{
								HandleType h = ComposeHandle(FACE,e);
								if( !GetMarker(h,hm) )
									Face(this,h)->FixNormalOrientation();
							}
						}
					}
					else
					{
						//std::cout << "Fix orientation" << std::endl;
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer e = 0; e < FaceLastLocalID(); ++e) if( isValidElement(FACE,e) )
							Face(this,ComposeHandle(FACE,e))->FixNormalOrientation();
					}
				}
				ShowGeometricData(ORIENTATION,FACE);
			}
			if( types == MEASURE )
			{
				//std::cout << "MEASURE" << std::endl;
				for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
				{
					if( (mask & etype) && !HaveGeometricData(MEASURE,etype))
					{
						measure_tag = CreateTag(measure_name,DATA_REAL,etype,NONE,1);
						if( HideMarker() )
						{
							MarkerType hm = HideMarker();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
							{
								HandleType h = ComposeHandle(etype,e);
								if( !GetMarker(h,hm) ) GetGeometricData(h,MEASURE,static_cast<Storage::real *>(MGetDenseLink(h,measure_tag)));
							}
						}
						else
						{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
							{
								HandleType h = ComposeHandle(etype,e);
								GetGeometricData(h,MEASURE,static_cast<Storage::real *>(MGetDenseLink(h,measure_tag)));
							}
						}
						ShowGeometricData(MEASURE,etype);
					}
				}
			}
			if( types == CENTROID )
			{
				//std::cout << "CENTROID" << std::endl;
				for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
				{
					if( (mask & etype) && !HaveGeometricData(CENTROID,etype))
					{
						centroid_tag = CreateTag(centroid_name,DATA_REAL,etype,NONE,GetDimensions());
						if( HideMarker() )
						{
							MarkerType hm = HideMarker();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k) )
							{
								HandleType h = ComposeHandle(etype,k);
								if( !GetMarker(h,hm) ) GetGeometricData(h,CENTROID,static_cast<Storage::real *>(MGetDenseLink(h,centroid_tag)));
							}
						}
						else
						{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k) )
							{
								HandleType h = ComposeHandle(etype,k);
								GetGeometricData(h,CENTROID,static_cast<Storage::real *>(MGetDenseLink(h,centroid_tag)));
							}
						}
						ShowGeometricData(CENTROID,etype);
					}
				}
			}
			if( types == BARYCENTER )
			{
				//std::cout << "BARYCENTER" << std::endl;
				for(ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
				{
					if( (mask & etype) && !HaveGeometricData(BARYCENTER,etype))
					{
						barycenter_tag = CreateTag(barycenter_name,DATA_REAL,etype,NONE,GetDimensions());
						if( HideMarker() )
						{
							MarkerType hm = HideMarker();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
							{
								HandleType h = ComposeHandle(etype,e);
								if( !GetMarker(h,hm) ) GetGeometricData(h,BARYCENTER,static_cast<Storage::real *>(MGetDenseLink(h,barycenter_tag)));
							}
						}
						else
						{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
							{
								HandleType h = ComposeHandle(etype,e);
								GetGeometricData(h,BARYCENTER,static_cast<Storage::real *>(MGetDenseLink(h,barycenter_tag)));
							}
						}
						ShowGeometricData(BARYCENTER,etype);
					}
				}	
			}
			if( types == NORMAL )
			{
				//std::cout << "NORMAL" << std::endl;
				for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype))
				{
					if( (mask & etype) && !HaveGeometricData(NORMAL,etype))
					{
						normal_tag = CreateTag(normal_name,DATA_REAL,etype,NONE,GetDimensions());
						if( HideMarker() )
						{
							MarkerType hm = HideMarker();
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
							{
								HandleType h = ComposeHandle(etype,e);
								if( !GetMarker(h,hm) ) GetGeometricData(h,NORMAL,static_cast<Storage::real *>(MGetDenseLink(h,normal_tag)));
							}
						}
						else
						{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
							for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
							{
								HandleType h = ComposeHandle(etype,e);
								GetGeometricData(h,NORMAL,static_cast<Storage::real *>(MGetDenseLink(h,normal_tag)));
							}
						}
						ShowGeometricData(NORMAL,etype);
					}
				}
			}
		}
	}

	bool Cell::CheckConvexity() const {return GetMeshLink()->CheckConvexity(getFaces());}

	bool Mesh::CheckConvexity(const ElementArray<Face> & faces) const
	{
		if( !faces.empty() )
		{
			std::vector<Storage::real> x(faces.size()*6,0.0);
			Storage::real * n = &x[faces.size()*3];
			Storage::real eps = GetEpsilon();
			//precompute data for all faces
#if defined(USE_OMP)
#pragma omp critical (reorder_edges)
#endif
			{ //collect data safely so that nobody changes it
				for(ElementArray<Face>::size_type j = 0; j != faces.size(); ++j)
				{
					faces[j].Centroid(&x[j*3]);
					faces[j].UnitNormal(&n[j*3]);
				}
			}
			//run comparison
			for(ElementArray<Face>::size_type j = 0; j != faces.size(); ++j)
			{
				Storage::real dota = 0, dots = 0, dotv, v[3] = {0.,0.,0.};
				for(ElementArray<Face>::size_type m = 0; m != faces.size(); ++m)
				{
					vec_diff(&x[m*3],&x[j*3],v,3);
					dotv = vec_dot_product(&n[j*3],v,3);
					dots += dotv;
					dota += fabs(dotv);
				}
				if( fabs(fabs(dots) - dota) > eps )
					return false;
			}
		}
		return true;
	}
	
	void Mesh::FacesOrientation(ElementArray<Face> & faces, MarkerType rev)
	{
		//can copy orientation-independent algorithm from
		//incident_matrix.hpp: incident_matrix::compute_measure
		//assume mdim is of size 3 at most
		if( !faces.empty() )
		{
			if( CheckConvexity(faces) ) //simpler algorithm
			{
				std::vector<Storage::real> x(faces.size()*6);
				Storage::real * n = &x[faces.size()*3];
				Storage::real xc[3] = {0.,0.,0.}, v[3] = {0.,0.,0.};
				for(ElementArray<Face>::size_type j = 0; j != faces.size(); ++j)
				{
					faces[j].Centroid(&x[j*3]);
					faces[j].UnitNormal(&n[j*3]);
					xc[0] += x[j*3+0];
					xc[1] += x[j*3+1];
					xc[2] += x[j*3+2];
				}
				xc[0] /= (Storage::real)faces.size();
				xc[1] /= (Storage::real)faces.size();
				xc[2] /= (Storage::real)faces.size();
				for(ElementArray<Face>::size_type j = 0; j != faces.size(); ++j)
				{
					vec_diff(&x[j*3],xc,v,3);
					if( vec_dot_product(&n[j*3],v,3) < 0 )
					{
						if( isPrivate(rev) )
							faces[j].SetPrivateMarker(rev);
						else
							faces[j].SetMarker(rev);
					}
				}
			}
			else
			{
				//real was = *ret/3.0;
				Face cur = faces[0];
				Mesh * mesh = faces.GetMeshLink();
				//firstly, have to figure out orientation of each face
				//mark all faces, so that we can perform adjacency retrival
				MarkerType mrk = mesh->CreatePrivateMarker();
				//MarkerType rev = mesh->CreatePrivateMarker(); //reverse orientation
				faces.SetPrivateMarker(mrk); //0-th face orientation is default
				cur->RemPrivateMarker(mrk);
				Node n1,n2; //to retrive edge
				bool reverse = false; //reverse orientation in considered face
				std::deque< orient_face > stack; //edge and first node and face for visiting
				ElementArray<Edge> edges = cur->getEdges();
				do
				{
					//figure out starting node order
					if( edges[0]->getBeg() == edges[1]->getBeg() ||
						edges[0]->getBeg() == edges[1]->getEnd() )
					{
						n1 = edges[0]->getEnd();
						n2 = edges[0]->getBeg();
					}
					else
					{
						n1 = edges[0]->getBeg();
						n2 = edges[0]->getEnd();
					}
					//schedule unvisited adjacent faces
					for(unsigned j = 0; j < edges.size(); j++)
					{
						//schedule face adjacent to considered edge
						ElementArray<Face> adjacent = edges[j]->getFaces(mrk);
						assert(adjacent.size() <= 1);
						if( !adjacent.empty() )
						{
							adjacent.RemPrivateMarker(mrk);
							stack.push_back(orient_face(edges[j],reverse ? n2 : n1,adjacent[0]));
						}
						//update edge nodes
						n1 = n2; //current end is new begin
						//find new end
						if( n2 == edges[(j+1)%edges.size()]->getBeg() )
							n2 = edges[(j+1)%edges.size()]->getEnd();
							else
								n2 = edges[(j+1)%edges.size()]->getBeg();
								}
					if( stack.empty() ) break;
					//get entry from stack
					orient_face r = stack.front();
					//remove face from stack
					stack.pop_front();
					//retrive edges for new face
					edges = r.face->getEdges();
					reverse = false;
					//figure out starting node order
					if( edges[0]->getBeg() == edges[1]->getBeg() ||
						edges[0]->getBeg() == edges[1]->getEnd() )
					{
						n1 = edges[0]->getEnd();
						n2 = edges[0]->getBeg();
					}
					else
					{
						n1 = edges[0]->getBeg();
						n2 = edges[0]->getEnd();
					}
					//find out common edge orientation
					for(unsigned j = 0; j < edges.size(); j++)
					{
						if( edges[j] == r.bridge ) //found the edge
						{
							//reverse ordering on this face
							if( r.first == n1 )
							{
								if( isPrivate(rev) )
									r.face->SetPrivateMarker(rev);
								else
									r.face->SetMarker(rev);
								reverse = true;
							}
							break;
						}
						//update edge nodes
						n1 = n2; //current end is new begin
						//find new end
						if( n2 == edges[(j+1)%edges.size()]->getBeg() )
							n2 = edges[(j+1)%edges.size()]->getEnd();
						else
							n2 = edges[(j+1)%edges.size()]->getBeg();
					}
				} while(true);
				faces.RemPrivateMarker(mrk);
				mesh->ReleasePrivateMarker(mrk);
				//faces.RemPrivateMarker(rev);
			}
		}
	}
	
	void Mesh::GetGeometricData(HandleType e, GeometricData type, Storage::real * ret)
	{
		assert(e != InvalidHandle());
		assert(ret != NULL);
		assert(type == MEASURE ||
			   type == CENTROID ||
			   type == BARYCENTER ||
			   type == NORMAL);
		ElementType etype = GetHandleElementType(e);
		integer edim = Element::GetGeometricDimension(GetGeometricType(e));
		integer mdim = GetDimensions();
		switch(type)
		{
			case MEASURE:
			if( HaveGeometricData(MEASURE,etype) )
			{
				if (UpdateGeometryMarker() && GetMarker(e, UpdateGeometryMarker())) RecomputeGeometricData(e);
				*ret = static_cast<Storage::real *>(MGetDenseLink(e,measure_tag))[0];
				//~ if( isnan(*ret) || fabs(*ret) < 1e-15  ) throw -1;
			}
			else
			{
				switch(edim)
				{
					case 0: *ret = 0; break;
					case 1: //length of edge
					{
						ElementArray<Node> nodes = Element(this,e)->getNodes();
						if( nodes.size() > 1 )
						{
							real c[3] = {0,0,0};
							real_array v0 = nodes[0].Coords();
							real_array v1 = nodes[1].Coords();
							vec_diff(v0.data(),v1.data(),c,mdim);
							*ret = vec_len(c,mdim);
						}
						else *ret = 0;
						//~ if( isnan(*ret) || fabs(*ret) < 1e-15  ) throw -1;
						break;
					}
					case 2: //area of face
					{
						ElementArray<Node> nodes = Element(this,e)->getNodes();
						if( nodes.size() > 2 )
						{
							*ret = 0;
							real nt[3] = { 0,0,0 }, l1[3] = { 0,0,0 }, l2[3] = { 0,0,0 }, n0[3] = { 0,0,0 };// , ss, at;
							real_array v0, v1, v2;
							v0 = nodes[0].Coords();
							for(int i = 1; i < (int)nodes.size()-1; i++)
							{
								v1 = nodes[i].Coords();
								v2 = nodes[i+1].Coords();
								vec_diff(v1,v0,l1,mdim);
								vec_diff(v2,v0,l2,mdim);
								vec_cross_product(l1,l2,nt);
								for(int q = 0; q < 3; ++q)
									n0[q] += nt[q]*0.5;
							}
							*ret = vec_len(n0,3);
							//code below is not really needed
							/*
							for(int i = 1; i < (int)nodes.size()-1; i++)
							{
								v1 = nodes[i].Coords();
								v2 = nodes[i+1].Coords();
								vec_diff(v1,v0,l1,mdim);
								vec_diff(v2,v0,l2,mdim);
								vec_cross_product(l1,l2,nt);
								for(int q = 0; q < 3; ++q) 
									nt[q] *= 0.5;
								ss = vec_dot_product(n0,nt,3);
								if( ss ) ss /= fabs(ss);
								at = sqrt(vec_dot_product(nt,nt,3))*ss;
								*ret += at;
							}
							*/
							if( *ret != *ret ) std::cout << "area is nan" << std::endl;
							*ret = fabs(*ret);
						} else *ret = 0;
						//~ if( isnan(*ret) || fabs(*ret) < 1e-15  ) throw -1;
						break;
					}
					case 3: //volume of cell
					{
						//bool print = false;
//redo:
						Cell me = Cell(this,e);
						ElementArray<Face> faces = me->getFaces();
						bool ornt = true;//!HaveGeometricData(ORIENTATION,FACE);
						//bool ornt = !HaveGeometricData(ORIENTATION,FACE);
						//bool ornt = !CheckConvexity(faces);
						MarkerType rev = 0;
						if( ornt )
						{
							rev = CreatePrivateMarker();
							FacesOrientation(faces,rev);
						}
						real vol = 0, a, at, volp;
						real x[3] = {0,0,0}, n0[3] = {0,0,0}, s, ss;
						real l1[3] = {0,0,0}, l2[3] = {0,0,0};
						real nt[3] = {0,0,0};
						real vc[3] = {0,0,0};
						//me.Centroid(cx);
						//if( print ) std::cout << "cx: " << cx[0] << " " << cx[1] << " " << cx[2] << std::endl;
						for(unsigned j = 0; j < faces.size(); j++)
						{
							//compute normal to face
							ElementArray<Node> nodes = faces[j].getNodes();
							if( ornt )
								s = faces[j].GetPrivateMarker(rev) ? -1.0 : 1.0;
							else
								s = faces[j].FaceOrientedOutside(me) ? 1.0 : -1.0;
							x[0] = x[1] = x[2] = 0;
							n0[0] = n0[1] = n0[2] = 0;
							vc[0] = vc[1] = vc[2] = 0;
							a = 0;
							real_array v0 = nodes[0].Coords(), v1, v2;
							for(int i = 1; i < (int)nodes.size()-1; i++)
							{
								v1 = nodes[i].Coords();
								v2 = nodes[i+1].Coords();
								vec_diff(v1,v0,l1,mdim);
								vec_diff(v2,v0,l2,mdim);
								vec_cross_product(l1,l2,nt);
								for(int q = 0; q < 3; ++q)
									n0[q] += nt[q]*0.5;
							}
							for(int i = 1; i < (int)nodes.size()-1; i++)
							{
								v1 = nodes[i].Coords();
								v2 = nodes[i+1].Coords();
								vec_diff(v1,v0,l1,mdim);
								vec_diff(v2,v0,l2,mdim);
								vec_cross_product(l1,l2,nt);
								for(int q = 0; q < 3; ++q) 
									nt[q] *= 0.5;
								ss = vec_dot_product(n0,nt,3);
								if( ss ) 
									ss /= fabs(ss);
								at = sqrt(vec_dot_product(nt,nt,3))*ss;
								//same as faces[j].Barycenter(x)
								for(int q = 0; q < 3; ++q)
									x[q] += at*(v0[q]+v1[q]+v2[q])/3.0;
								a += at;
							}
							if( a )
							{
								for(int q = 0; q < 3; ++q) 
									x[q] = x[q]/a;
							}
							else Element(this,e).Centroid(x);
							volp = vec_dot_product(x,n0,3) / 3.0;
							vol += s*volp;
						}
						if( ornt )
						{
							if( vol < 0.0 ) vol = -vol;
							faces.RemPrivateMarker(rev);
							ReleasePrivateMarker(rev);
						}
						assert(vol > 0);
						*ret = vol;
						break;
					}
				}
			}
			//~ if( isnan(*ret) || fabs(*ret) < 1e-15  ) throw -1;
			break;
			case CENTROID:
				if(etype == NODE )
					memcpy(ret,MGetDenseLink(e,CoordsTag()),sizeof(real)*mdim);
				else if(HaveGeometricData(CENTROID,etype))
				{
					if (UpdateGeometryMarker() && GetMarker(e, UpdateGeometryMarker())) RecomputeGeometricData(e);
					memcpy(ret,MGetDenseLink(e,centroid_tag),sizeof(real)*mdim);
				}
				else
				{
					ElementArray<Node> nodes = Element(this,e).getNodes();
					memset(ret,0,sizeof(real)*mdim);
					for(unsigned k = 0; k < nodes.size(); ++k)
					{
						for(int q = 0; q < mdim; ++q)
							ret[q] += nodes[k].Coords()[q];
					}
					for(int q = 0; q < mdim; ++q)
						ret[q] /= (real)nodes.size();
				}
			break;
			case BARYCENTER:
			if(etype == NODE )
				memcpy(ret,MGetDenseLink(e,CoordsTag()),sizeof(real)*mdim);
			else if(HaveGeometricData(BARYCENTER,etype))
			{
				if (UpdateGeometryMarker() && GetMarker(e, UpdateGeometryMarker())) RecomputeGeometricData(e);
				memcpy(ret,MGetDenseLink(e,barycenter_tag),sizeof(real)*mdim);
			}
			else
			{
				memset(ret,0,sizeof(real)*mdim);
				if( edim == 1 )
				{
					ElementArray<Node> n = Element(this,e)->getNodes();
					if( n.size() == 2 )
					{
						real_array v0 = n[0].Coords();
						real_array v1 = n[1].Coords();
						for(integer j = 0; j < dim; j++)
							ret[j] = (v0[j] + v1[j])*0.5;
					}
					else if( n.size() == 1 )
					{
						real_array v0 = n[0].Coords();
						for(integer j = 0; j < dim; j++) ret[j] = v0[j];
					}
				}
				else if( edim == 2 )
				{
					ElementArray<Node> nodes = Element(this,e)->getNodes();
					if( nodes.size() > 2 )
					{
						*ret = 0;
						real nt[3] = {0,0,0}, l1[3] = {0,0,0}, l2[3] = {0,0,0};
						real c[3] = {0,0,0}, n0[3] = {0,0,0}, ss;
						real_array v0 = nodes[0].Coords(), v1, v2;
						for(int i = 1; i < (int)nodes.size()-1; i++)
						{
							v1 = nodes[i].Coords();
							v2 = nodes[i+1].Coords();
							vec_diff(v1,v0,l1,mdim);
							vec_diff(v2,v0,l2,mdim);
							vec_cross_product(l1,l2,nt);
							for(int q = 0; q < 3; ++q)
								n0[q] += nt[q]*0.5;
						}
						real a = 0, at;
						for(int i = 1; i < (int)nodes.size()-1; i++)
						{
							v1 = nodes[i].Coords();
							v2 = nodes[i+1].Coords();
							vec_diff(v1,v0,l1,mdim);
							vec_diff(v2,v0,l2,mdim);
							vec_cross_product(l1,l2,nt);
							for(int q = 0; q < 3; ++q) 
								nt[q] *= 0.5;
							ss = vec_dot_product(n0,nt,3);
							if( ss ) ss /= fabs(ss);
							at = sqrt(vec_dot_product(nt,nt,3))*ss;
							for(int q = 0; q < mdim; ++q)
								c[q] += at*(v0[q]+v1[q]+v2[q])/3.0;
							a += at;
						}
						if( a )
						{
							for(int q = 0; q < mdim; ++q) 
								ret[q] = c[q]/a;
						} 
						else Element(this,e).Centroid(ret);
					}
				}
				else if( edim == 3 )
				{
					Cell me = Cell(this,e);
					ElementArray<Face> faces = me->getFaces();
					bool ornt = true;//!HaveGeometricData(ORIENTATION,FACE);
					//bool ornt = !CheckConvexity(faces);
					//bool ornt = !HaveGeometricData(ORIENTATION,FACE);
					MarkerType rev = 0;
					if( ornt )
					{
						rev = CreatePrivateMarker();
						FacesOrientation(faces,rev);
					}
					real vol = 0, a, at, volp;
					real x[3] = {0,0,0}, nt[3] = {0,0,0}, s;
					real c[3] = {0,0,0};
					real n0[3] = {0,0,0}, ss;
					real l1[3] = {0,0,0}, l2[3] = {0,0,0};
					real vc[3] = {0,0,0};
					for(unsigned j = 0; j < faces.size(); j++)
					{
						//compute normal to face
						ElementArray<Node> nodes = faces[j].getNodes();
						if( ornt )
							s = faces[j].GetPrivateMarker(rev) ? -1.0 : 1.0;
						else
							s = faces[j].FaceOrientedOutside(me) ? 1.0 : -1.0;
						n0[0] = n0[1] = n0[2] = 0;
						x[0] = x[1] = x[2] = 0;
						vc[0] = vc[1] = vc[2] = 0;
						a = 0;
						real_array v0 = nodes[0].Coords(), v1, v2;
						for(int i = 1; i < (int)nodes.size()-1; i++)
						{
							v1 = nodes[i].Coords();
							v2 = nodes[i+1].Coords();
							vec_diff(v1,v0,l1,mdim);
							vec_diff(v2,v0,l2,mdim);
							vec_cross_product(l1,l2,nt);
							for(int q = 0; q < 3; ++q)
								n0[q] += nt[q]*0.5;
						}
						for(int i = 1; i < (int)nodes.size()-1; i++)
						{
							v1 = nodes[i].Coords();
							v2 = nodes[i+1].Coords();
							vec_diff(v1,v0,l1,mdim);
							vec_diff(v2,v0,l2,mdim);
							vec_cross_product(l1,l2,nt);
							for(int q = 0; q < 3; ++q)
								nt[q] *= 0.5;
							ss = vec_dot_product(n0,nt,3);
							if( ss ) 
								ss /= fabs(ss);
							at = sqrt(vec_dot_product(nt,nt,3))*ss;
							//same as faces[j].Centroid(x)
							for(int q = 0; q < 3; ++q)
							{
								x[q] += at*(v0[q]+v1[q]+v2[q])/3.0;
								c[q] += s*nt[q]*(pow(v0[q]+v1[q],2)+pow(v1[q]+v2[q],2)+pow(v2[q]+v0[q],2))/24.0;
							}
							a += at;
						}
						if( a )
						{
							for(int q = 0; q < 3; ++q) 
								x[q] = x[q]/a;
						} else faces[j].Centroid(x);
						volp = vec_dot_product(x,n0,3) / 3.0;
						//~ for(int q = 0; q < 3; ++q)
							//~ c[q] += 0.75*(x[q]-cx[q])*s*volp; // x - base barycenter, cx - apex
						vol += s*volp;
					}
					if( ornt )
					{
						if( vol < 0.0 )
						{
							vol = -vol;
							for(int q = 0; q < 3; ++q)
								c[q] = -c[q];
						}
						faces.RemPrivateMarker(rev);
						ReleasePrivateMarker(rev);
					}
					if( vol ) 
					{
						for(int q = 0; q < mdim; ++q)
							ret[q] = c[q]/vol;// + cx[q];
					}
					else me.Centroid(ret);//for(int q = 0; q < mdim; ++q) ret[q] = cx[q];
					//std::cout << ret[0] << " " << ret[1] << " " << ret[2] << std::endl;
				}
			}
			break;
			case NORMAL:
			{
				if( HaveGeometricData(NORMAL,etype) )
				{
					if (UpdateGeometryMarker() && GetMarker(e, UpdateGeometryMarker())) RecomputeGeometricData(e);
					memcpy(ret,MGetDenseLink(e,normal_tag),sizeof(real)*mdim);
				}
				else
				{
					memset(ret,0,sizeof(real)*mdim);
					if( edim == 2 )//&& mdim == 3)
					{
						ElementArray<Node> nodes = Element(this,e)->getNodes();
						real n[3] = {0,0,0}, l1[3] = {0,0,0}, l2[3] = {0,0,0};
						real nt[3] = {0,0,0};
						real_array v0 = nodes[0].Coords();
						for(int i = 1; i < (int)nodes.size()-1; i++)
						{
							real_array v1 = nodes[i].Coords();
							real_array v2 = nodes[i+1].Coords();
							vec_diff(v1,v0,l1,mdim);
							vec_diff(v2,v0,l2,mdim);
							vec_cross_product(l1,l2,nt);
							for(int q = 0; q < 3; ++q)
								n[q] += nt[q]*0.5;
						}
						for(int q = 0; q < mdim; ++q)
							ret[q] = n[q];
					}
					else if( edim == 1 )//&& mdim == 2 )
					{
						ElementArray<Node> nodes = Element(this,e)->getNodes();
						if( nodes.size() > 1 )
						{
							Storage::real_array a = nodes[0].Coords();
							Storage::real_array b = nodes[1].Coords();
							ret[0] = b[1] - a[1];
							ret[1] = a[0] - b[0];
							Storage::real l = ::sqrt(ret[0]*ret[0]+ret[1]*ret[1]);
							if( l )
							{
								ret[0] /= l;
								ret[1] /= l;
							}
							l = ::sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]));
							ret[0] *= l;
							ret[1] *= l;
						}
					}
				}
			}
			break;
		}
		//~ if( type == MEASURE )
		//~ {
			//~ if( isnan(*ret) || fabs(*ret) < 1e-15  ) throw -1;
		//~ }
	}
	


	bool Element::Planarity() const
	{
		Mesh * m = GetMeshLink();
		integer dim = m->GetDimensions();
		if( dim < 3 ) return true;
		ElementArray<Node> nodes = getNodes();
		if( nodes.size() <= 3 ) return true;
		/*
		ElementArray<Node>::size_type i, s = nodes.size();
		Storage::real c[3] = {0,0,0}, n[3] = {0,0,0}, v[3] = {0,0,0};
		getAsFace()->Normal(n);
		Centroid(c);
		for(i = 0; i < s; i++)
		{
			vec_diff(nodes[i].Coords().data(),c,v,3);
			if( ::fabs(vec_dot_product(v,n,3)) > m->GetEpsilon() ) return false;
		}
		*/
		INMOST_DATA_REAL_TYPE l1[3] = {0,0,0}, l2[3] = {0,0,0}, nt[3] = {0,0,0}, n0[3] = {0,0,0};
		real_array v0, v1, v2;
		v0 = nodes[0].Coords();
		for(int i = 1; i < (int)nodes.size()-1; i++)
		{
			v1 = nodes[i].Coords();
			v2 = nodes[i+1].Coords();
			vec_diff(v1,v0,l1,dim);
			vec_diff(v2,v0,l2,dim);
			vec_cross_product(l1,l2,nt);
			for(int q = 0; q < 3; ++q)
				n0[q] += nt[q]*0.5;
		}
		vec_normalize(n0,3);
		for(int i = 1; i < (int)nodes.size()-1; i++)
		{
			v1 = nodes[i].Coords();
			v2 = nodes[i+1].Coords();
			vec_diff(v1,v0,l1,dim);
			vec_diff(v2,v0,l2,dim);
			vec_cross_product(l1,l2,nt);
			vec_normalize(nt,3);
			if( fabs(fabs(vec_dot_product(n0,nt,3))-1) > m->GetEpsilon() ) 
				return false;
		}
		
		return true;
	}


	bool Cell::Closure() const
	{
		adj_type & lc = GetMeshLink()->LowConn(GetHandle());
		return lc.size() > 0 ? GetMeshLink()->TestClosure(lc.data(),static_cast<integer>(lc.size())) : false;
	}

	bool Face::Closure() const
	{
		adj_type & lc = GetMeshLink()->LowConn(GetHandle());
		return lc.size() > 0 ? GetMeshLink()->TestClosure(lc.data(),static_cast<integer>(lc.size())) : false;
	}



	

	

	bool Element::Boundary() const
	{
		switch(GetElementType())
		{
			case FACE:
				if( !getAsFace()->FrontCell().isValid() )
					return true;
				return false;
			case CELL:
			case EDGE:
			case NODE:
			{
				ElementArray<Element> faces = getAdjElements(FACE);
				for(ElementArray<Element>::iterator it = faces.begin(); it != faces.end(); it++)
					if( it->Boundary() ) return true;
				return false;
			}
			default: return false;
		}
		return false;
	}
	
	
	
	std::ostream & spaces(std::ostream & out, int print)
	{
		while(print) {out << " "; print--;}
		return out;
	}
	
	
	bool Face::CheckNormalOrientation() const
	{
		Mesh * mesh = GetMeshLink();
        //integer dim = mesh->GetDimensions();
		Cell c1 = BackCell();
		if( c1.isValid() )
		{
			if( c1.CheckConvexity() ) //simpler algorithm
			{
				Storage::real xc[3] = {0.,0.,0.};
				Storage::real xf[3] = {0.,0.,0.};
				Storage::real nf[3] = {0.,0.,0.};
				Storage::real v[3] = {0.,0.,0.};
				c1.Centroid(xc);
				Centroid(xf);
				UnitNormal(nf);
				vec_diff(xf,xc,v,3);
				if( vec_dot_product(nf,v,3) < 0 )
					return false;
			}
			else
			{
				MarkerType mrk = mesh->CreatePrivateMarker();
				MarkerType rev = mesh->CreatePrivateMarker(); //reverse orientation
				Storage::real measure = 0;
				ElementArray<Face> data = c1.getFaces();
				Face cur = *this;
				//firstly, have to figure out orientation of each face
				//mark all faces, so that we can perform adjacency retrival
#if defined(USE_OMP)
#pragma omp critical (reorder_edges)
#endif
				{
					data.SetPrivateMarker(mrk); //0-th face orientation is default
					cur->RemPrivateMarker(mrk);
					Node n1,n2; //to retrive edge
					bool reverse = false; //reverse orientation in considered face
					std::deque< orient_face > stack; //edge and first node and face for visiting
					ElementArray<Edge> edges;
					edges = cur->getEdges();
					do
					{
						//figure out starting node order
						if( edges[0]->getBeg() == edges[1]->getBeg() ||
							edges[0]->getBeg() == edges[1]->getEnd() )
						{
							n1 = edges[0]->getEnd();
							n2 = edges[0]->getBeg();
						}
						else
						{
							n1 = edges[0]->getBeg();
							n2 = edges[0]->getEnd();
						}
						//schedule unvisited adjacent faces
						for(unsigned j = 0; j < edges.size(); j++)
						{
							//schedule face adjacent to considered edge
							ElementArray<Face> adjacent = edges[j]->getFaces(mrk);
							assert(adjacent.size() <= 1);
							if( !adjacent.empty() )
							{
								adjacent.RemPrivateMarker(mrk);
								stack.push_back(orient_face(edges[j],reverse ? n2 : n1,adjacent[0]));
							}
							//update edge nodes
							n1 = n2; //current end is new begin
							//find new end
							if( n2 == edges[(j+1)%edges.size()]->getBeg() )
								n2 = edges[(j+1)%edges.size()]->getEnd();
							else
								n2 = edges[(j+1)%edges.size()]->getBeg();
						}
						if( stack.empty() ) break;
						//get entry from stack
						orient_face r = stack.front();
						//remove face from stack
						stack.pop_front();
						//retrive edges for new face
						edges = r.face->getEdges();
						reverse = false;
						//figure out starting node order
						if( edges[0]->getBeg() == edges[1]->getBeg() ||
						edges[0]->getBeg() == edges[1]->getEnd() )
						{
							n1 = edges[0]->getEnd();
							n2 = edges[0]->getBeg();
						}
						else
						{
							n1 = edges[0]->getBeg();
							n2 = edges[0]->getEnd();
						}
						//find out common edge orientation
						for(unsigned j = 0; j < edges.size(); j++)
						{
							if( edges[j] == r.bridge ) //found the edge
							{
								//reverse ordering on this face
								if( r.first == n1 )
								{
									r.face->SetPrivateMarker(rev);
									reverse = true;
								}
								break;
							}
							//update edge nodes
							n1 = n2; //current end is new begin
							//find new end
							if( n2 == edges[(j+1)%edges.size()]->getBeg() )
								n2 = edges[(j+1)%edges.size()]->getEnd();
							else
								n2 = edges[(j+1)%edges.size()]->getBeg();
						}
					} while(true);
					data.RemPrivateMarker(mrk);
					mesh->ReleasePrivateMarker(mrk);
					Storage::real nrm[3], cnt[3];
					for(unsigned j = 0; j < data.size(); j++)
					{
						//std::cout << (data[j].GetPrivateMarker(rev) ? 0:1);
						//compute normal to face
						data[j].Barycenter(cnt);
						data[j].Normal(nrm);
						measure += (data[j]->GetPrivateMarker(rev) ? -1.0 : 1.0)*vec_dot_product(cnt,nrm,3);
					}
				}
				//std::cout << "cur" << (cur->GetPrivateMarker(rev) ? 0:1) << " " << measure << " ";
				//bool have_rev = cur->GetPrivateMarker(rev);
				data.RemPrivateMarker(rev);
				mesh->ReleasePrivateMarker(rev);
				if( (measure < 0 ))// && !have_rev) || (measure > 0 && have_rev))
					return false;
			}
		}
		return true;
	}
	
	bool Face::FixNormalOrientation(bool allow_swap) const
	{
		if( !CheckNormalOrientation() )
		{
			if( allow_swap && FrontCell().isValid() )
				SwapCells();
			else
			{
#if defined(USE_OMP)
#pragma omp critical (reorder_edges)
#endif
				{
					ReorderEdges(); //this is not thread-safe with respect to CheckNormalOrientation
					if( GetMeshLink()->HaveGeometricData(NORMAL,FACE) )
					{
						real_array nrm = GetMeshLink()->RealArrayDF(GetHandle(),GetMeshLink()->GetGeometricTag(NORMAL));
						for(real_array::size_type it = 0; it < nrm.size(); ++it)
							nrm[it] = -nrm[it];
					}
				}
			}
			return true;
		}
		return false;
	}

	Storage::real meantri(Storage::real * v0, Storage::real * v1, Storage::real * v2, Storage::integer dim, Storage::real (*func)(Storage::real* x,Storage::real), Storage::real time)
	{
		Storage::real value = 0;
		static const Storage::real w[4] =   { -0.149570044467670, 0.175615257433204, 0.053347235608839 , 0.077113760890257};
		static const Storage::real a[4][3] =
		{
			{0.333333333333333,0.333333333333333,0.333333333333333},
			{0.479308067841923,0.260345966079038,0.260345966079038},
			{0.869739794195568,0.065130102902216,0.065130102902216},
			{0.638444188569809,0.312865496004875,0.048690315425316}
		};
		Storage::real XYG[13][3];
		for (Storage::integer i = 0 ; i < dim; i++)
			XYG[0][i] = 0.33333333333333333333*(v0[i]+v1[i]+v2[i]);
		 value += w[0] * func(XYG[0],time);
		for (int i = 0 ; i < 3 ; i++ )
		{
			for (Storage::integer j = 0 ; j < dim; j++)
				XYG[1+i][j] = v0[j] + (v1[j] - v0[j]) * a[1][i] + (v2[j] - v0[j])*a[1][(i+1)%3];
			value += w[1] * func(XYG[1+i],time);
		}
		for (int i = 0 ; i < 3 ; i++ )
		{
			for (Storage::integer j = 0 ; j < dim; j++)
				XYG[4+i][j] = v0[j] + (v1[j] - v0[j]) * a[2][i] + (v2[j] - v0[j])*a[2][(i+1)%3];
			value += w[2] * func(XYG[4+i],time);
		}
		for (int i = 0 ; i < 3 ; i++ )
		{
			for (Storage::integer j = 0 ; j < dim; j++)
			{
				XYG[7+2*i][j] = v0[j] + (v1[j] - v0[j]) * a[3][i] + (v2[j] - v0[j])*a[3][(i+1)%3];
				XYG[8+2*i][j] = v0[j] + (v1[j] - v0[j]) * a[3][(i+1)%3] + (v2[j] - v0[j])*a[3][i];
			}
			value += w[3] * func(XYG[7+2*i],time);
			value += w[3] * func(XYG[8+2*i],time);
		}
		return value;
	}

	Storage::real meantet(Storage::real * v0, Storage::real * v1, Storage::real * v2, Storage::real * v3, Storage::real (*func)(Storage::real* x,Storage::real),Storage::real time)
	{
		Storage::real value;
		static const Storage::real T5A = 0.25, W5A = 0.11851851851852;
		static const Storage::real T5B = 0.09197107805272, T5C = 0.72408676584183, W5B = 0.07193708377902;
		static const Storage::real T5D = 0.31979362782963, T5E = 0.04061911651111;
		static const Storage::real W5C = 0.06906820722627;
		static const Storage::real T5F = 0.05635083268963, T5G = 0.44364916731037;
		static const Storage::real W5D = 0.05291005291005;
		static const Storage::real w[15] = {W5A,W5B,W5B,W5B,W5B,W5C,W5C,W5C,W5C,W5D,W5D,W5D,W5D,W5D,W5D};
		Storage::real XYG[15][3];
		for (int i = 0; i < 3 ;i++)
		{
			XYG[0][i] = T5A * (v0[i] + v1[i] + v2[i] + v3[i]);
			XYG[1][i] = T5C * v0[i] + T5B * (v1[i]+v2[i]+v3[i]);
			XYG[2][i] = T5C * v1[i] + T5B * (v0[i]+v2[i]+v3[i]);
			XYG[3][i] = T5C * v2[i] + T5B * (v0[i]+v1[i]+v3[i]);
			XYG[4][i] = T5C * v3[i] + T5B * (v0[i]+v1[i]+v2[i]);

			XYG[5][i] = T5E * v0[i] + T5D * (v1[i]+v2[i]+v3[i]);
			XYG[6][i] = T5E * v1[i] + T5D * (v0[i]+v2[i]+v3[i]);
			XYG[7][i] = T5E * v2[i] + T5D * (v0[i]+v1[i]+v3[i]);
			XYG[8][i] = T5E * v3[i] + T5D * (v0[i]+v1[i]+v2[i]);

			XYG[9][i] = T5F * (v0[i] + v1[i]) + T5G * (v2[i] + v3[i]);
			XYG[10][i] = T5G * (v0[i] + v1[i]) + T5F * (v2[i] + v3[i]);
			XYG[11][i] = T5F * (v0[i] + v3[i]) + T5G * (v1[i] + v2[i]);
			XYG[12][i] = T5G * (v0[i] + v3[i]) + T5F * (v1[i] + v2[i]);
			XYG[13][i] = T5F * (v0[i] + v2[i]) + T5G * (v1[i] + v3[i]);
			XYG[14][i] = T5G * (v0[i] + v2[i]) + T5F * (v1[i] + v3[i]);
		}
		value = 0;
		for (int i = 0 ; i < 15 ; i++)
			value += w[i] * func(XYG[i],time);
		return value;
	}

	Storage::real Element::Mean(Storage::real (*func)(Storage::real* x,Storage::real),Storage::real time) const
	{
		Mesh * m = GetMeshLink();
		if( GetElementDimension() == 2 )
		{
			integer dim = m->GetDimensions();
			real val = 0, vol = 0, tvol,tval;
			real v1[3] = {0,0,0},v2[3] = {0,0,0}, product[3] = {0,0,0};
			ElementArray<Node> nodes = getNodes();
			real_array av0 = nodes.front().Coords();
			for(ElementArray<Node>::iterator it = ++nodes.begin(); it != nodes.end(); it++)
			{
				ElementArray<Node>::iterator jt = it++;
				if( it == nodes.end() ) break;
				real_array av1 = jt->Coords();
				real_array av2 = it->Coords();
				tval = meantri(av0.data(),av1.data(),av2.data(),dim,func,time);
				vec_diff(av1.data(),av0.data(),v1,dim);
				vec_diff(av2.data(),av0.data(),v2,dim);
				vec_cross_product(v1,v2,product);
				tvol = sqrt(vec_dot_product(product,product,dim))*0.5;
				val += tval*tvol;
				vol += tvol;
				it = jt;
			}
			return val / vol;
		}
		else if( GetElementDimension() == 3 )
		{
			
			ElementArray<Element> rfaces = getAdjElements(FACE);
			array<int> n(static_cast<array<int>::size_type>(rfaces.size()));
			array<real> v;
			int k = 0;
			for(ElementArray<Element>::iterator f = rfaces.begin(); f != rfaces.end(); f++)
			{
				ElementArray<Node> nodes = f->getNodes();
				int nn = n[k] = static_cast<int>(nodes.size());
				for(int i = 0; i < nn; i++)
				{
					real_array a = nodes[i].Coords();
					v.insert(v.end(),a.begin(),a.end());
				}
				k++;
			}
			int j = 0;
			real x[3] = {0,0,0}, y[3], d, vol = 0, c;
			for(array<int>::size_type i = 0; i < n.size(); i++)
			{
				y[0] = y[1] = y[2] = d = 0;
				for(array<int>::size_type k = 1; k < static_cast<array<int>::size_type>(n[i] - 1); k++)
				{
					d += c = det3v(&v[j*3],&v[(j+k)*3],&v[(j+k+1)*3]);
					y[0] += c * (v[j*3+0] + v[(j+k)*3+0] + v[(j+k+1)*3+0]);
					y[1] += c * (v[j*3+1] + v[(j+k)*3+1] + v[(j+k+1)*3+1]);
					y[2] += c * (v[j*3+2] + v[(j+k)*3+2] + v[(j+k+1)*3+2]);
				}
				int orient = rfaces[static_cast<ElementArray<Element>::size_type>(i)]->getAsFace()->FaceOrientedOutside(getAsCell())?1:-1;
				x[0] += orient*y[0];
				x[1] += orient*y[1];
				x[2] += orient*y[2];
				vol += orient*d;
				j += n[i];
			}
			x[0] /= 4.0 * vol;
			x[1] /= 4.0 * vol;
			x[2] /= 4.0 * vol;
			vol /= 6;
			j = 0;
			vol = 0;
			real tvol, tval, val = 0;
			real vv0[3], vv1[3], vv2[3], prod[3];
			for(array<int>::size_type i = 0; i < n.size(); i++)
			{
				for(array<int>::size_type k = 1; k < static_cast<array<int>::size_type>(n[i] - 1); k++)
				{
					vec_diff(&v[j*3],x,vv0,3);
					vec_diff(&v[(j+k)*3],x,vv1,3);
					vec_diff(&v[(j+k+1)*3],x,vv2,3);
					vec_cross_product(vv1,vv2,prod);
					tvol = vec_dot_product(vv0,prod,3)/6.0 * (rfaces[static_cast<ElementArray<Element>::size_type>(i)]->getAsFace()->FaceOrientedOutside(getAsCell())?1:-1);
					tval = meantet(x,&v[j*3],&v[(j+k)*3],&v[(j+k+1)*3],func,time);
					val += tval * tvol;
					vol += tvol;
				}
				j+=n[i];
			}
			return val / vol;
		}
		else if ( GetElementDimension() == 1 ) //Mean value over line.
		{
			ElementArray<Node> nodes = getNodes();
			integer dim = m->GetDimensions();
			real_array x1 = nodes[0].Coords();
			real_array x2 = nodes[1].Coords();
			real middle[3];
			for (integer i = 0 ; i < dim ; i++) middle[i] = (x1[i]+x2[i])*0.5;
			//Simpson formula
			return (func(x1.data(),time) + 4*func(middle,time) + func(x2.data(),time))/6.0;
		}
		return 0;
	}
	
	const Tag & Mesh::GetGeometricTag(GeometricData type) const
	{
		switch(type) 
		{
		case    MEASURE: return    measure_tag; 
		case   CENTROID: return   centroid_tag; 
		case BARYCENTER: return barycenter_tag;
		case     NORMAL: return     normal_tag;
		}
		assert(false);
		static const Tag t; //invalid tag
		return t;
	}
	
	ElementArray<Face> Mesh::GatherBoundaryFaces()
	{
		ElementArray<Face> ret(this);
		for(integer f = 0; f < FaceLastLocalID(); f++) if( isValidElement(FACE,f) )
		{
			Face ef = Face(this,ComposeHandle(FACE,f));
			if( ef->Boundary() ) ret.push_back(ef);
		}
		return ret;
	}

	ElementArray<Face> Mesh::GatherInteriorFaces()
	{
		ElementArray<Face> ret(this);
		if( GetMeshState() == Mesh::Serial )
		{
			for(integer f = 0; f < FaceLastLocalID(); f++) if( isValidElement(FACE,f) )
			{
				Face ef = Face(this,ComposeHandle(FACE,f));
				if( ef->nbAdjElements(CELL) == 2 ) ret.push_back(ef);
			}
		}
		else
		{
			for(integer f = 0; f < FaceLastLocalID(); f++) if( isValidElement(FACE,f) )
			{
				Face ef = Face(this,ComposeHandle(FACE,f));
				ElementArray<Cell> cells = ef->getCells();
				if( cells.size() == 2 && cells[0].GetStatus() != Element::Ghost && cells[1].GetStatus() != Element::Ghost)
					ret.push_back(ef);
			}
		}
		return ret;
	}

	Storage::integer Mesh::CountBoundaryFaces()
	{
		integer ret = 0;
		for(integer f = 0; f < FaceLastLocalID(); f++) if( isValidElement(FACE,f) )
		{
			Face ef = Face(this,ComposeHandle(FACE,f));
			if( ef->Boundary() ) ret++;
		}
		return ret;
	}

	Storage::integer Mesh::CountInteriorFaces()
	{
		integer ret = 0;
		if( GetMeshState() == Mesh::Serial )
		{
			for(integer f = 0; f < FaceLastLocalID(); f++) if( isValidElement(FACE,f) )
			{
				Face ef = Face(this,ComposeHandle(FACE,f));
				if( ef->nbAdjElements(CELL) == 2 ) ret++;
			}
		}
		else
		{
			for(integer f = 0; f < FaceLastLocalID(); f++) if( isValidElement(FACE,f) )
			{
				Face ef = Face(this,ComposeHandle(FACE,f));
				ElementArray<Cell> cells = ef->getCells();
				if( cells.size() == 2 && cells[0].GetStatus() != Element::Ghost && cells[1].GetStatus() != Element::Ghost)
					ret++;
			}
		}

		return ret;
	}
	
	

	
	void Element::CastRay(const real * pos, const real * dir, std::map<HandleType, real> & hits) const
	{
		Mesh * mm = GetMeshLink();
		unsigned dim = mm->GetDimensions();
		Storage::real eps = mm->GetEpsilon();
		switch(GetElementType())
		{
			case NODE:
			{
				Storage::real_array coords = getAsNode()->Coords();
				Storage::real vec[3],ndir[3], lvec, ldir;
				for(unsigned k = 0; k < dim; k++)
				{
					vec[k] = (coords[k]-pos[k]);
					ndir[k] = dir[k];
				}
				lvec = vec_normalize(vec,dim);
				ldir = vec_normalize(ndir,dim);
				if( vec_dot_product(vec,ndir,dim) >= 1.0 - eps )
					hits[GetHandle()] = lvec/ldir;
			}	 
			break;
			case EDGE:
			{
				throw NotImplemented;
			}
			break;
			case FACE:
			{
				HandleType shr_nodes[2];
				Storage::real tri[3][3];
				Centroid(tri[2]);
				ElementArray<Node> nodes = getNodes();
				memcpy(tri[0],nodes[0].Coords().data(),sizeof(Storage::real)*dim);
				for(ElementArray<Node>::size_type q = 0; q < nodes.size(); q++)
				{
					memcpy(tri[1],nodes[(q+1)%nodes.size()].Coords().data(),sizeof(Storage::real)*dim);
					Storage::real a[3],b[3],c[3],n[3], d, k, m, l;
					Storage::real dot00,dot01,dot02,dot11,dot12,invdenom, uq,vq;
					a[0] = tri[0][0] - tri[2][0];
					a[1] = tri[0][1] - tri[2][1];
					b[0] = tri[1][0] - tri[2][0];
					b[1] = tri[1][1] - tri[2][1];
					if( dim > 2 )
					{
						a[2] = tri[0][2] - tri[2][2];
						b[2] = tri[1][2] - tri[2][2];
					}
					else a[2] = b[2] = 0;
					vec_cross_product(a,b,n);
					d = -vec_dot_product(n,tri[2],dim);
					m =  vec_dot_product(n,dir,dim);
					if( !(::fabs(m) < 1.0e-25) )
					{
						l = vec_dot_product(n,pos,dim);
						k = -(d + l)/m;
						if( k >= 0 )
						{
							c[0] = pos[0] + k*dir[0] - tri[2][0];
							c[1] = pos[1] + k*dir[1] - tri[2][1];
							if( dim > 2 ) c[2] = pos[2] + k*dir[2] - tri[2][2]; else c[2] = 0;
							dot00 = vec_dot_product(a,a,dim);
							dot01 = vec_dot_product(a,b,dim);
							dot02 = vec_dot_product(a,c,dim);
							dot11 = vec_dot_product(b,b,dim);
							dot12 = vec_dot_product(b,c,dim);
							invdenom = (dot00*dot11 - dot01*dot01);
							uq = (dot11*dot02-dot01*dot12);
							vq = (dot00*dot12-dot01*dot02);
							//std::cout << "uq " << uq << " vq " << vq << " invdenom " << invdenom << std::endl;
							if( !(::fabs(invdenom) < 1.0e-25  ))//&& ::fabs(uq) >= 0.0 && ::fabs(vq) >= 0.0) )
							{
								uq = uq/invdenom;
								vq = vq/invdenom;
								if( uq >= -eps && vq >= -eps )
								{
									if( 1.0-(uq+vq) >= -eps )
									{
										if( 1.0-(uq+vq) <= eps )
										{
											shr_nodes[0] = nodes.at(q);
											shr_nodes[1] = nodes.at((q+1)%nodes.size());
											if( uq <= eps && vq >= -eps )
												hits[shr_nodes[0]] = k;
											else if ( uq >= -eps && vq <= eps )
												hits[shr_nodes[1]] = k;
											else
												hits[mm->FindSharedAdjacency(shr_nodes,2)]= k;
										} 
										else
										{
											hits[GetHandle()] = k;
										}
										break; //we shouldn't have more then one intersection
									}
								}
							}
						}
					}
					memcpy(tri[0],tri[1],sizeof(Storage::real)*dim);
				}
			}
			break;
			case CELL:
			{
				ElementArray<Face> faces = getFaces();
				for(ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
					it->CastRay(pos,dir,hits);
			}
			break;
			/*
			case ESET: //cast ray through set of elements, probably accelerate by tree
			throw NotImplemented;
			break;
			case MESH: //cast ray through all elements, probably accelerate by tree
			throw NotImplemented;
			break;
			*/
		}
	}

	/// Wachspress interpolation from nodes to arbitrary point inside 2d face
	/// 	x	point, should be inside face
	/// 	f	2d face
	/// 	nodes_stencil	array of pairs (handle of node, interpolation coefficient)
	///
	/// Example: 
	///		std::map<HandleType,double> nodes_stencil;
	///		wachspress_2d(x, f, nodes_stencil);
	///		double xvalue = 0.0;
	///		for(std::map<HandleType,double>::const_iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
	///		{
	///			xvalue += it->second * tag_value.Real(it->first);
	///		}
	void Mesh::WachspressInterpolation2D(const real * x, const Face & f, std::map<HandleType,real> & nodes_stencil) const
	{
		ElementArray<Node> nodes = f.getNodes();
		int n = (int)nodes.size();
		INMOST_DATA_REAL_TYPE xprev[3], xthis[3], xnext[3], dx[3];
		INMOST_DATA_REAL_TYPE weight_sum = 0.0;
		INMOST_DATA_REAL_TYPE areaface = f.Area();
		bool edge_interpolation = false;
		int edgenode1 = -1, edgenode2 = -1;
		for(int i = 0; i < n && !edge_interpolation; i++)
		{
			int prev = (i+n-1)%n, next = (i+1)%n;
			nodes[i].Centroid(xthis);
			vec_diff(x,xthis,dx,3);
			nodes[prev].Centroid(xprev);
			nodes[next].Centroid(xnext);
			INMOST_DATA_REAL_TYPE area1 = triarea3d(x,xprev,xthis);
			edge_interpolation = area1 < 1e-8 * areaface;
			if( !edge_interpolation )
			{
				INMOST_DATA_REAL_TYPE area2 = triarea3d(x,xthis,xnext);
				edge_interpolation = area2 < 1e-8 * areaface;
				if( !edge_interpolation )
				{
					INMOST_DATA_REAL_TYPE area = triarea3d(xthis,xprev,xnext);
					INMOST_DATA_REAL_TYPE weight = area / (area1 * area2);
					weight_sum += weight;
					nodes_stencil[nodes[i].GetHandle()] = weight;
				}
				else
				{
					edgenode1 = i;
					edgenode2 = next;
				}
			}
			else
			{
				edgenode1 = prev;
				edgenode2 = i;
			}
		}
		if( edge_interpolation )
		{
			if( !nodes_stencil.empty() )	nodes_stencil.clear();
			INMOST_DATA_REAL_TYPE x1[3],x2[3], v12[3], v[3], alpha;
			nodes[edgenode1].Centroid(x1);
			nodes[edgenode2].Centroid(x2);
			vec_diff(x2,x1,v12,3);
			vec_diff(x,x1,v,3);
			alpha = vec_dot_product(v12,v,3) / vec_dot_product(v12,v12,3);

			nodes_stencil[nodes[edgenode1].GetHandle()] = alpha;
			nodes_stencil[nodes[edgenode2].GetHandle()] = 1.0-alpha;
		}
		else if( weight_sum )
			for(std::map<HandleType,real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
				it->second /= weight_sum;
	}
	
	void Mesh::WachspressInterpolation3D(const real * x, const Cell & c, std::map<HandleType,real> & nodes_stencil) const
	{
		ElementArray<Node> nodes = c.getNodes();
		int n = (int)nodes.size();
		ElementArray<Face> cfaces = c.getFaces(), faces;
		Face swp;

		INMOST_DATA_REAL_TYPE weight_sum = 0;
		for(int i = 0; i < n; i++)
		{
			INMOST_DATA_REAL_TYPE xthis[3], v[3];
			nodes[i].Centroid(xthis);

			faces = cfaces;
			faces.Intersect(nodes[i].getFaces());
			int nf = (int)faces.size();
			INMOST_DATA_REAL_TYPE nrm_mean[3] = {0,0,0}, cross[3];
			INMOST_DATA_REAL_TYPE * nrm = new INMOST_DATA_REAL_TYPE[nf*3];	// oriented unit normal to face
			INMOST_DATA_REAL_TYPE * h = new INMOST_DATA_REAL_TYPE[nf];		// distance to face
			for(int j = 0; j < nf; j++)
			{
				faces[j].OrientedUnitNormal(c,nrm+3*j);
				for(int k = 0; k < 3; k++)	nrm_mean[k] += nrm[3*j+k];
			}
			for(int k = 0; k < 3; k++)	nrm_mean[k] /= nf;
			// sort faces counter-clockwise relative to mean normal
			for(int j = 0; j < nf-1; j++)
				for(int k = j+1; k < nf-1; k++)
				{
					vec_cross_product(nrm+3*j, nrm+3*k, cross);
					if( vec_dot_product(nrm_mean,cross,3) < 0 )
					{
						swp = faces[j];
						faces[j] = faces[k];
						faces[k] = swp;
						for(int l = 0; l < 3; l++)
						{
							cross[l] = nrm[3*j+l];
							nrm[3*j+l] = nrm[3*k+l];
							nrm[3*k+l] = cross[l];
						}
					}
				}

			vec_diff(xthis,x,v,3);
			for(int j = 0; j < nf; j++)
			{
				h[j] = vec_dot_product(nrm+3*j,v,3);
				if( h[j] < -1e-8 )
				{
					//std::cerr << __FILE__ << ":" << __LINE__ << " point is outside of cell" << std::endl;
					h[j] = fabs(h[j]);
				}
				else if(h[j] < 1e-8) // point on the face plane
				{
					INMOST_DATA_REAL_TYPE xf[3];
					for(int k = 0; k < 3; k++)	xf[k] = x[k] + h[j]*nrm[3*j+k];
					if( !nodes_stencil.empty() )	nodes_stencil.clear();
					WachspressInterpolation2D(xf,faces[j],nodes_stencil);
					return;
				}
			}

			INMOST_DATA_REAL_TYPE weight = 0.0;
			for(int j = 0; j < nf-2; j++)
			{
				INMOST_DATA_REAL_TYPE tet_weight = fabs(det3v(nrm+3*j, nrm+3*(j+1), nrm+3*(nf-1))) / (h[j]*h[j+1]*h[nf-1]);
				weight += tet_weight;
			}
			nodes_stencil[nodes[i].GetHandle()] = weight;
			weight_sum += weight;
			
			delete [] h;
			delete [] nrm;
		}

		if( weight_sum )
			for(std::map<HandleType,real>::iterator it = nodes_stencil.begin(); it != nodes_stencil.end(); ++it)
				it->second /= weight_sum;
	}


	
	void UnpackBoundary(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
			const Storage::integer * recv = static_cast<const Storage::integer *>(static_cast<const void *>(data));
			Storage::integer_array arr = element->IntegerArray(tag);
			for(INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k)
			{
				if( std::find(arr.begin(),arr.end(),recv[k]) == arr.end() )
					arr.push_back(recv[k]);
			}
		}
	}

	void Mesh::MarkBoundaryFaces(MarkerType boundary_marker) 
	{
		if( GetProcessorsNumber() > 1 )
		{
			Tag tag_bnd = CreateTag("CALC_BOUNDARY",DATA_INTEGER,FACE,NONE);
			//we need shared unique numbers on cells
			if( !HaveGlobalID(CELL) ) AssignGlobalID(CELL);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(integer it = 0; it < FaceLastLocalID(); ++it) if( isValidFace(it) )
			{
				Face face = FaceByLocalID(it);
				integer_array arr = face->IntegerArray(tag_bnd);
				ElementArray<Cell> adj = face->getCells();
				for(ElementArray<Cell>::iterator jt = adj.begin(); jt != adj.end(); jt++)
					arr.push_back(jt->GlobalID());
			}
			ReduceData(tag_bnd,FACE,0,UnpackBoundary);
			ExchangeData(tag_bnd,FACE,0);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(integer it = 0; it < FaceLastLocalID(); ++it) if( isValidFace(it) )
			{
				Face face = FaceByLocalID(it);
				assert(face->IntegerArray(tag_bnd).size() <= 2);
				if( face->IntegerArray(tag_bnd).size() == 1 )
					face->SetMarker(boundary_marker);
				else
					face->RemMarker(boundary_marker);
			}
			DeleteTag(tag_bnd);
		}
		else //nothing to worry about
		{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(integer it = 0; it < FaceLastLocalID(); ++it) if( isValidFace(it) )
			{
				Face face = FaceByLocalID(it);
				if( face->Boundary() ) 
					face->SetMarker(boundary_marker);
				else 
					face->RemMarker(boundary_marker);
			}

		}
	}

}
#endif
