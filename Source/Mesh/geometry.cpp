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
	
	bool Cell::Inside(const Storage::real * point) const//check for 2d case
	{
		Mesh * mesh = GetMeshLink();
		integer dim = GetElementDimension();
		if( dim == 3 )
		{
			/*
			tiny_map<HandleType,real,16> hits;
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
		tiny_map<HandleType,int,64> e_visit;
		tiny_map<HandleType,int,64>::iterator it;
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
	
	Element::GeometricType Mesh::ComputeGeometricType(ElementType etype, const HandleType * lc, INMOST_DATA_ENUM_TYPE s) const
	{
		Element::GeometricType ret = Element::Unset;
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

				if( Element::GetGeometricDimension(GetGeometricType(lc[0])) == 0 )
				{ 
					ret = Element::Line;
				}
				else
				{
					if( !GetTopologyCheck(NEED_TEST_CLOSURE) || TestClosure(lc,s) )
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
				if(  Element::GetGeometricDimension(GetGeometricType(lc[0])) == 1 )
				{
					if( !GetTopologyCheck(NEED_TEST_CLOSURE) || TestClosure(lc,s) )
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
					if( !GetTopologyCheck(NEED_TEST_CLOSURE) ||  TestClosure(lc,s) )
					{
						//test c_faces closure, if no closure, set as MultiPolygon
						INMOST_DATA_ENUM_TYPE quads = 0,tris = 0,i;
						for(i = 0; i < s; i++)
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
		Element::adj_type const & lc = LowConn(h);
		if( !lc.empty() )
			SetGeometricType(h,ComputeGeometricType(GetHandleElementType(h),lc.data(),static_cast<integer>(lc.size())));
	}
	
	void Mesh::RecomputeGeometricData(HandleType e)
	{
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


		if( GetHandleElementType(e) == CELL && HaveGeometricData(ORIENTATION,FACE)) //then correct the normal
		{
			Element::adj_type & lc = LowConn(e);
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)
				if( !GetMarker(*it,HideMarker()) && HighConn(*it).size() == 1 )
				{
					Face(this,*it)->FixNormalOrientation();
				}
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
		std::sort(&*table.begin(),&*table.end());
		for(GeomParam::iterator it = table.begin(); it != table.end(); ++it)
		{
			GeometricData types = it->first;
			ElementType mask = it->second;
			if( types == ORIENTATION )
			{
				//std::cout << "ORIENTATION" << std::endl;
				if( mask & FACE )
				{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for(integer e = 0; e < FaceLastLocalID(); ++e) 
					{
						if( isValidElement(FACE,e) )
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
						{
							HandleType h = ComposeHandle(etype,e);
							GetGeometricData(h,MEASURE,static_cast<Storage::real *>(MGetDenseLink(h,measure_tag)));
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k) )
						{
							HandleType h = ComposeHandle(etype,k);
							GetGeometricData(h,CENTROID,static_cast<Storage::real *>(MGetDenseLink(h,centroid_tag)));
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
						{
							HandleType h = ComposeHandle(etype,e);
							GetGeometricData(h,BARYCENTER,static_cast<Storage::real *>(MGetDenseLink(h,barycenter_tag)));
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
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer e = 0; e < LastLocalID(etype); ++e) if( isValidElement(etype,e) )
						{
							HandleType h = ComposeHandle(etype,e);
							GetGeometricData(h,NORMAL,static_cast<Storage::real *>(MGetDenseLink(h,normal_tag)));
						}
						ShowGeometricData(NORMAL,etype);
					}
				}
			}
		}
	}
	
	void Mesh::GetGeometricData(HandleType e, GeometricData type, Storage::real * ret)
	{
		assert(e != InvalidHandle());
		assert(ret != NULL);
		assert(type == MEASURE || type == CENTROID || type == BARYCENTER || type == NORMAL);
		ElementType etype = GetHandleElementType(e);
		integer edim = Element::GetGeometricDimension(GetGeometricType(e));
		integer mdim = GetDimensions();
		switch(type)
		{
			case MEASURE:
			if( HaveGeometricData(MEASURE,etype) )
			{
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
							Storage::real c[3];
							vec_diff(nodes[0]->Coords().data(),nodes[1]->Coords().data(),c,mdim);
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
							real x[3] = {0,0,0};
							Storage::real_array x0 = nodes[0].Coords();
							for(ElementArray<Node>::size_type i = 1; i < nodes.size()-1; i++)
							{
								Storage::real_array v1 = nodes[i].Coords();
								Storage::real_array v2 = nodes[i+1].Coords();
								if( mdim == 3 )
								{
									x[0] += (v1[1]-x0[1])*(v2[2]-x0[2]) - (v1[2]-x0[2])*(v2[1]-x0[1]);
									x[1] += (v1[2]-x0[2])*(v2[0]-x0[0]) - (v1[0]-x0[0])*(v2[2]-x0[2]);
								}
								x[2] += (v1[0]-x0[0])*(v2[1]-x0[1]) - (v1[1]-x0[1])*(v2[0]-x0[0]);
							}
							*ret = ::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*0.5;
						} else *ret = 0;
						//~ if( isnan(*ret) || fabs(*ret) < 1e-15  ) throw -1;
						break;
					}
					case 3: //volume of cell
					{
						
						Cell me = Cell(this,e);
						ElementArray<Face> faces = me->getFaces();
						*ret = 0;
						//can copy orientation-independent algorithm from
						//incident_matrix.hpp: incident_matrix::compute_measure
						//assume mdim is of size 3 at most
						dynarray<Storage::real,3> fcnt(mdim), fnrm(mdim), ccnt(mdim);
						/*
						me->Centroid(ccnt.data());
						for(ElementArray<Face>::size_type i = 0; i < faces.size(); i++)
						{
							faces[i]->Centroid(fcnt.data());
							for(int r = 0; r < mdim; ++r)
								fcnt[r] = fcnt[r]-ccnt[r];
							faces[i]->OrientedNormal(me,fnrm.data());
							*ret += vec_dot_product(fcnt.data(),fnrm.data(),mdim);
						}
						if( *ret < 0.0 ) //a robust algorithm that can handle unoriented cell
						 */
						if( !faces.empty() )
						{
							//real was = *ret/3.0;
							*ret = 0;
							Face cur = faces[0];
							Cell c1 = me;
							Mesh * mesh = c1.GetMeshLink();
							//firstly, have to figure out orientation of each face
							//mark all faces, so that we can perform adjacency retrival
							MarkerType mrk = mesh->CreatePrivateMarker();
							MarkerType rev = mesh->CreatePrivateMarker(); //reverse orientation
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
							faces.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
							c1->Centroid(ccnt.data());
							for(unsigned j = 0; j < faces.size(); j++)
							{
								//compute normal to face
								faces[j].Centroid(fcnt.data());
								faces[j].Normal(fnrm.data());
								for(int r = 0; r < mdim; ++r)
									fcnt[r] = fcnt[r]-ccnt[r];
								*ret += (faces[j]->GetPrivateMarker(rev) ? -1.0 : 1.0)*vec_dot_product(fcnt.data(),fnrm.data(),3);
							}
							*ret = fabs(*ret);
							faces.RemPrivateMarker(rev);
							mesh->ReleasePrivateMarker(rev);
							
							//std::cout << "volume is " << *ret/3.0 << " was " << was << " for " << me->LocalID() << std::endl;
						}
						*ret /= 3.0;
						
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
				memcpy(ret,MGetDenseLink(e,centroid_tag),sizeof(real)*mdim);
			}
			else
			{
				ElementArray<Node> nodes = Element(this,e)->getNodes();
				memset(ret,0,sizeof(real)*mdim);
				if(nodes.size() != 0)
				{
					Storage::real div = 1.0/nodes.size();
					for(ElementArray<Node>::size_type i = 0; i < nodes.size(); i++)
					{
						Storage::real_array c =nodes[i].Coords();
						for(integer j = 0; j < mdim; j++) ret[j] += c[j];
					}
					for(integer j = 0; j < mdim; j++) ret[j] *= div;
				}
			}
			break;
			case BARYCENTER:
			if( etype == NODE )
				memcpy(ret,MGetDenseLink(e,CoordsTag()),sizeof(real)*mdim);
			else if(HaveGeometricData(BARYCENTER,etype))
				memcpy(ret,MGetDenseLink(e,barycenter_tag),sizeof(real)*mdim);
			else
			{
				memset(ret,0,sizeof(real)*mdim);
				if( edim == 1 )
				{
					ElementArray<Node> n = Element(this,e)->getNodes();
					if( n.size() == 2 )
					{
						Storage::real_array a = n[0].Coords();
						Storage::real_array b = n[1].Coords();
						for(integer j = 0; j < dim; j++) 
							ret[j] = (a[j] + b[j])*0.5;
					}
					else if( n.size() == 1 )
					{
						Storage::real_array a = n[0].Coords();
						for(integer j = 0; j < dim; j++) ret[j] = a[j];
					}
				}
				else if( edim == 2 )
				{
					ElementArray<Node> nodes = Element(this,e)->getNodes();
					real s,d, x1[3] = {0,0,0},x2[3] = {0,0,0},x[3] = {0,0,0};
					//here we compute area of polygon
					//~ if( HaveGeometricData(MEASURE,etype) && HaveGeometricData(NORMAL,etype) )
					//~ {
						//~ s = e->RealDF(measure_tag);
						//~ memcpy(x,&e->RealDF(normal_tag),sizeof(Storage::real)*mdim);
					//~ }
					//~ else
					//~ {
					if( nodes.size() > 2 )
					{
						Storage::real_array x0 = nodes[0].Coords();
						for(ElementArray<Node>::size_type i = 1; i < nodes.size()-1; i++)
						{
							Storage::real_array v1 = nodes[i].Coords();
							Storage::real_array v2 = nodes[i+1].Coords();
							if( mdim == 3 )
							{
								x[0] += (v1[1]-x0[1])*(v2[2]-x0[2]) - (v1[2]-x0[2])*(v2[1]-x0[1]);
								x[1] += (v1[2]-x0[2])*(v2[0]-x0[0]) - (v1[0]-x0[0])*(v2[2]-x0[2]);
							}
							x[2] += (v1[0]-x0[0])*(v2[1]-x0[1]) - (v1[1]-x0[1])*(v2[0]-x0[0]);
						}
						s = ::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
						x[0] /= s; x[1] /= s; x[2] /= s; //here we obtain the unit normal
						//~ }
						//here we compute the center
						Storage::real_array v0 = nodes[0].Coords();
						for(ElementArray<Node>::size_type j = 1; j < nodes.size()-1; j++)
						{
							Storage::real_array v1 = nodes[j].Coords();
							Storage::real_array v2 = nodes[j+1].Coords();
							for(integer k = 0; k < mdim; k++)
							{
								x1[k] = v0[k] - v1[k];
								x2[k] = v0[k] - v2[k];
							}
							d = det3v(x1,x2,x); //here we use unit normal
							for(integer k = 0; k < mdim; k++) 
								ret[k] += d*(v0[k]+v1[k]+v2[k]);
						}
						for(integer k = 0; k < mdim; k++) ret[k] /= 3.0 * s;
					}
					else
					{
						for(ElementArray<Node>::size_type i = 0; i < nodes.size(); i++)
						{
							Storage::real_array c = nodes[i].Coords();
							for(integer k = 0; k < mdim; k++) ret[k] += c[k];
						}
						for(integer k = 0; k < mdim; k++) ret[k] /= static_cast<Storage::real>(nodes.size());
					}
				}
				else if( edim == 3 )
				{
					ElementArray<Face> faces = Element(this,e)->getFaces();
					real d,c,vol = 0, y[3];
					for(ElementArray<Face>::size_type i = 0; i < faces.size(); i++)
					{
						d = y[0] = y[1] = y[2] = 0;
						ElementArray<Node> nodes = faces[i].getNodes();
						Storage::real_array v0 = nodes[0].Coords();
						for(ElementArray<Node>::size_type j = 1; j < nodes.size()-1; j++)
						{
							Storage::real_array v1 = nodes[j].Coords();
							Storage::real_array v2 = nodes[j+1].Coords();
							c = det3v(v0.data(),v1.data(),v2.data());
							d += c;
							y[0] += c * (v0[0] + v1[0] + v2[0]);
							y[1] += c * (v0[1] + v1[1] + v2[1]);
							y[2] += c * (v0[2] + v1[2] + v2[2]);
						}
						c = faces[i].FaceOrientedOutside(Cell(this,e)) ? 1 : -1;
						ret[0] += c * y[0];
						ret[1] += c * y[1];
						ret[2] += c * y[2];
						vol += c*d;
					}
					ret[0] /= vol*4;
					ret[1] /= vol*4;
					ret[2] /= vol*4;
				}
			}
			break;
			case NORMAL:
			{
//				real sret[3];
//				bool cmp = false;
				if( HaveGeometricData(NORMAL,etype) )
				{
					memcpy(ret,MGetDenseLink(e,normal_tag),sizeof(real)*mdim);
//					cmp = true;
				}
				else
				{
					memset(ret,0,sizeof(real)*mdim);
					if( edim == 2 )//&& mdim == 3)
					{
						ElementArray<Node> nodes = Element(this,e)->getNodes();
						
						Storage::real_array x0 = nodes[0].Coords(), a = x0, b;
						for(ElementArray<Node>::size_type i = 0; i < nodes.size(); i++)
						{
							b = nodes[(i+1)%nodes.size()].Coords();
							ret[0] += (a[1]-x0[1])*(b[2]-x0[2]) - (a[2]-x0[2])*(b[1]-x0[1]);
							ret[1] += (a[2]-x0[2])*(b[0]-x0[0]) - (a[0]-x0[0])*(b[2]-x0[2]);
							ret[2] += (a[0]-x0[0])*(b[1]-x0[1]) - (a[1]-x0[1])*(b[0]-x0[0]);
							a.swap(b);
						}
						/*
						 for(unsigned i = 0; i < nodes.size(); i++)
						 {
						 Storage::real_array a = nodes[i].Coords();
						 Storage::real_array b = nodes[(i+1)%nodes.size()].Coords();
						 ret[0] += a[1]*b[2] - a[2]*b[1];
						 ret[1] += a[2]*b[0] - a[0]*b[2];
						 ret[2] += a[0]*b[1] - a[1]*b[0];
						 }
						 */
						ret[0] *= 0.5;
						ret[1] *= 0.5;
						ret[2] *= 0.5;
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
//					real err = 0;
//					for(int k = 0; k < 3; ++k)
//						err += (ret[0]-sret[0])*(ret[0]-sret[0]);
//					err = sqrt(err);
//					if( err > 1.0e-8 && cmp )
//					{
//						std::cout << "face " << GetHandleID(e) << " rank " << GetProcessorRank() <<  std::endl;
//						std::cout << "cn " << ret[0] << " " << ret[1] << " " << ret[2] << std::endl;
//						std::cout << "sn " << sret[0] << " " << sret[1] << " " << sret[2] << std::endl;
//					}
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
		ElementArray<Node> p = getNodes();
		if( p.size() <= 3 ) return true;
		ElementArray<Node>::size_type i, s = p.size();
		Storage::real v[2][3] = {{0,0,0},{0,0,0}};
		vec_diff(p[1].Coords().data(),p[0].Coords().data(),v[0],3);
		vec_diff(p[2].Coords().data(),p[0].Coords().data(),v[1],3);
		vec_cross_product(v[0],v[1],v[1]);
		for(i = 3; i < s; i++)
		{
			vec_diff(p[i].Coords().data(),p[0].Coords().data(),v[0],3);
			if( ::fabs(vec_dot_product(v[0],v[1],3)) > m->GetEpsilon() ) return false;
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
			Storage::real measure = 0;
			ElementArray<Face> data = c1.getFaces();
			Face cur = *this;
			//firstly, have to figure out orientation of each face
			//mark all faces, so that we can perform adjacency retrival
			MarkerType mrk = mesh->CreatePrivateMarker();
			MarkerType rev = mesh->CreatePrivateMarker(); //reverse orientation
			data.SetPrivateMarker(mrk); //0-th face orientation is default
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
			Storage::real nrm[3], cnt[3], ccnt[3];
			c1->Centroid(ccnt);
			for(unsigned j = 0; j < data.size(); j++)
			{
				//std::cout << (data[j].GetPrivateMarker(rev) ? 0:1);
				//compute normal to face
				data[j].Centroid(cnt);
				data[j].Normal(nrm);
				for(int r = 0; r < 3; ++r)
					cnt[r] = cnt[r]-ccnt[r];
				measure += (data[j]->GetPrivateMarker(rev) ? -1.0 : 1.0)*vec_dot_product(cnt,nrm,3);
			}
			//std::cout << "cur" << (cur->GetPrivateMarker(rev) ? 0:1) << " " << measure << " ";
            //bool have_rev = cur->GetPrivateMarker(rev);
			data.RemPrivateMarker(rev);
			mesh->ReleasePrivateMarker(rev);
			if( (measure < 0 ))// && !have_rev) || (measure > 0 && have_rev))
				return false;
		}
		return true;
	}
	
	bool Face::FixNormalOrientation() const
	{
		if( !CheckNormalOrientation() )
		{
			if( FrontCell().isValid() )
				SwapCells();
			else
			{
				ReorderEdges(); //this is not thread-safe with respect to CheckNormalOrientation
				if( GetMeshLink()->HaveGeometricData(NORMAL,FACE) )
				{
					real_array nrm = GetMeshLink()->RealArrayDF(GetHandle(),GetMeshLink()->GetGeometricTag(NORMAL));
					for(real_array::size_type it = 0; it < nrm.size(); ++it)
						nrm[it] = -nrm[it];
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
			real normal[3];
			real v1[3] = {0,0,0},v2[3] = {0,0,0}, product[3] = {0,0,0};
			m->GetGeometricData(GetHandle(),NORMAL,normal);
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
				tvol = vec_dot_product(product,normal,dim)*0.5;
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
	
	

	
	void Element::CastRay(const real * pos, const real * dir, tiny_map<HandleType, real, 16> & hits) const
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


	
	void UnpackBoundary(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		if( size )
		{
			bool flag = true;
			const Storage::integer * recv = static_cast<const Storage::integer *>(static_cast<const void *>(data));
			Storage::integer_array arr = element->IntegerArray(tag);
			for(Storage::integer_array::iterator it = arr.begin(); it != arr.end(); it++)
				if( *it == recv[0] )
				{
					flag = false;
					break;
				}
				if( flag ) 
					arr.push_back(recv[0]);
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
				if( face->IntegerArray(tag_bnd).size() == 1 )
					face->SetMarker(boundary_marker);
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
				if( face->Boundary() ) face->SetMarker(boundary_marker);
			}

		}
	}

}
#endif
