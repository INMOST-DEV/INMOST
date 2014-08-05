#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{

	const char * Element::GeometricTypeName(GeometricType t)
	{
		switch(t)
		{
			case Unset:        return "Unset";
			case Vertex:       return "Vertex";
			case Line:         return "Line";
			case MultiLine:    return "MultiLine";
			case Tri:          return "Tri";
			case Quad:         return "Quad";
			case Polygon:      return "Polygon";
			case MultiPolygon: return "MultiPolygon";
			case Tet:          return "Tet";
			case Hex:          return "Hex";
			case Prism:        return "Prism";
			case Pyramid:      return "Pyramid";
			case Polyhedron:   return "Polyhedron";
			case Set:          return "Set";
		};
		return "Unset";
	}
	
	
	unsigned int Element::GetGeometricDimension(GeometricType m_type)
	{
		switch(m_type)
		{
			case Vertex:		return 0;
			case Line:			return 1;
			case MultiLine:
			case Tri:
			case Quad: 			
			case Polygon:		return 2;
			case MultiPolygon:
			case Tet:
			case Hex:
			case Prism:
			case Pyramid:
			case Polyhedron:	return 3;
			default: return UINT_MAX;
		}
		return UINT_MAX;
	}
	
	const char * Element::StatusName(Status s)
	{
		switch(s)
		{
			case Owned:  return "Owned";
			case Shared: return "Shared";
			case Ghost:  return "Ghost";
		}
		return "Unknown";
	}
	
	Element::Element(Mesh * m, ElementType _etype)
	:Storage(m,_etype)
	{
		m_type = Unset;
#if !defined(NDEBUG)
		assert( !(_etype & MESH) );
		assert( !(_etype & ESET) );
#endif
	}

	Element::Element(const Element & other) :Storage(other) {throw NotImplemented;}
	
	Element::Element(Mesh *m, INMOST_DATA_ENUM_TYPE lid, const Element & other)
	:Storage(m,lid,other)
	{
		m_type = other.m_type;
		high_conn.clear();
		low_conn.clear();
		//high_conn and low_conn may be filled only from outside
	}
	
	Element::~Element()
	{
		//~ std::cout << " deleted " << Element::GeometricTypeName(GetGeometricType()) << std::endl;
		Disconnect(true);
		high_conn.clear();
		low_conn.clear();
	}
	
	
	
	
	
	
	
	INMOST_DATA_ENUM_TYPE Element::nbAdjElements(ElementType _etype)const 
	{
		assert( !(_etype & MESH) );
		assert( !(_etype & ESET) );
		
		dynarray<Element *,128> result;
		unsigned int conn[4] = {0,0,0,0};
		unsigned int myconn, i, ret = 0;
		
		for(ElementType e = NODE, i = 0; e <= CELL; i++, e = e << 1)
		{
			if( _etype & e ) conn[i] = 1;
			if( GetElementType() & e ) myconn = i;
		}
		if( !GetMeshLink()->HideMarker() )
		{
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					ret += 1;
				}
				else if( i == (myconn + 1 + 4)%4 )
				{
					ret += high_conn.size();
				}
				else if( i == (myconn - 1 + 4)%4 )
				{
					ret += low_conn.size();
				}
				else if( i == (myconn - 2 + 4)%4 )
				{
					if( (GetElementType() & NODE) || (GetElementType() & EDGE) )
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = high_conn.begin(); it != high_conn.end(); it++)
							for(adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++)
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
						ret += result.size();
						for(dynarray<Element *,128>::iterator it = result.begin(); it != result.end(); it++)
							(*it)->RemMarker(mrk);
						
						result.clear();
						m->ReleaseMarker(mrk);
					}
					else if( GetElementType() & FACE )
					{
						ret += low_conn.size();
					}
					else
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = low_conn.begin(); it != low_conn.end(); it++)
							for(adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++)
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
								
						ret += result.size();
						for(dynarray<Element *,128>::iterator it = result.begin(); it != result.end(); it++)
							(*it)->RemMarker(mrk);
						result.clear();
						m->ReleaseMarker(mrk);
					}
				}
			}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					if( !GetMarker(hm) ) ret ++;
				}
				else if( i == (myconn + 1 + 4)%4 )
				{
					for(adj_type::const_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
						if( !(*it)->GetMarker(hm) ) ret++;
				}
				else if( i == (myconn - 1 + 4)%4 )
				{
					for(adj_type::const_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
						if( !(*it)->GetMarker(hm) ) ret++;
				}
				else if( i == (myconn - 2 + 4)%4 )
				{
					if( (GetElementType() & NODE) || (GetElementType() & EDGE) )
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = high_conn.begin(); it != high_conn.end(); it++) if( !(*it)->GetMarker(hm) )
							for(adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) if( !(*jt)->GetMarker(hm) )
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
						ret += result.size();
						for(dynarray<Element *,128>::iterator it = result.begin(); it != result.end(); it++)
							(*it)->RemMarker(mrk);
						
						result.clear();
						m->ReleaseMarker(mrk);
					}
					else
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = low_conn.begin(); it != low_conn.end(); it++) if( !(*it)->GetMarker(hm) )
							for(adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) if( !(*jt)->GetMarker(hm) )
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
								
						ret += result.size();
						for(dynarray<Element *,128>::iterator it = result.begin(); it != result.end(); it++)
							(*it)->RemMarker(mrk);
						result.clear();
						m->ReleaseMarker(mrk);
					}
				}
			}
		}
		return ret;
	}
	
	adjacent<Element> Element::getAdjElements(ElementType _etype) const 
	{
		assert( !(_etype & MESH) );
		assert( !(_etype & ESET) );
		
		adjacent<Element> result;
		unsigned int conn[4] = {0,0,0,0};
		unsigned int myconn, i = 0;
		
		for(ElementType e = NODE; e <= CELL; e = e << 1)
		{
			if( _etype & e ) conn[i] = 1;
			if( GetElementType() & e ) myconn = i;
			i++;
		}
		if( !GetMeshLink()->HideMarker() )
		{
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					result.push_back(const_cast<Element *>(this));
				}
				else if( i == (myconn + 1 + 4)%4 )
				{
					result.insert(result.end(),high_conn.begin(),high_conn.end());
				}
				else if( i == (myconn - 1 + 4)%4 )
				{
					result.insert(result.end(),low_conn.begin(),low_conn.end());
				}
				else if( i == (myconn - 2 + 4)%4 )
				{
					if( (GetElementType() & NODE) || (GetElementType() & EDGE) )
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = high_conn.begin(); it != high_conn.end(); it++)
							for(adj_type::iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++)
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
						for(adjacent<Element>::iterator it = result.begin(); it != result.end(); it++)
							it->RemMarker(mrk);
						m->ReleaseMarker(mrk);
					}
					else
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = low_conn.begin(); it != low_conn.end(); it++)
							for(adj_type::iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++)
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
						for(adjacent<Element>::iterator it = result.begin(); it != result.end(); it++)
							it->RemMarker(mrk);
						m->ReleaseMarker(mrk);
					}
				}
			}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					if( !GetMarker(hm) )
						result.push_back(const_cast<Element *>(this));
				}
				else if( i == (myconn + 1 + 4)%4 )
				{
					for(adj_type::const_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
						if( !(*it)->GetMarker(hm) ) result.push_back(*it);
					
				}
				else if( i == (myconn - 1 + 4)%4 )
				{
					for(adj_type::const_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
						if( !(*it)->GetMarker(hm) ) result.push_back(*it);
				}
				else if( i == (myconn - 2 + 4)%4 )
				{
					if( (GetElementType() & NODE) || (GetElementType() & EDGE) )
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = high_conn.begin(); it != high_conn.end(); it++) if( !(*it)->GetMarker(hm) )
							for(adj_type::iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) if( !(*jt)->GetMarker(hm) )
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
						for(adjacent<Element>::iterator it = result.begin(); it != result.end(); it++)
							it->RemMarker(mrk);
						m->ReleaseMarker(mrk);
					}
					else
					{
						Mesh * m = GetMeshLink();
						MIDType mrk = m->CreateMarker();
						for(adj_type::const_iterator it = low_conn.begin(); it != low_conn.end(); it++) if( !(*it)->GetMarker(hm) )
							for(adj_type::iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) if( !(*jt)->GetMarker(hm) )
								if( !(*jt)->GetMarker(mrk) )
								{
									result.push_back(*jt);
									(*jt)->SetMarker(mrk);
								}
						for(adjacent<Element>::iterator it = result.begin(); it != result.end(); it++)
							it->RemMarker(mrk);
						m->ReleaseMarker(mrk);
					}
				}
			}
		}
		return result;
	}
	
	
	
	Element & Element::operator =(Element const & other)
	{
		Storage::operator =(other);
		high_conn.clear();
		low_conn.clear();
		m_type = other.m_type;
		high_conn = other.high_conn;
		low_conn = other.low_conn;
		return *this;
	}
	
	
	adjacent<Node> Element::getNodes()
	{
		adjacent<Node> ret;
		ret.push_back(this);
		return ret; //this operation will be virtualized by proper algorithm	
	}
	adjacent<Edge> Element::getEdges()
	{
		adjacent<Edge> ret;
		ret.push_back(this);
		return ret; //this operation will be virtualized by proper algorithm	
	}
	adjacent<Face> Element::getFaces()
	{
		adjacent<Face> ret;
		ret.push_back(this);
		return ret; //this operation will be virtualized by proper algorithm	
	}
	adjacent<Cell> Element::getCells()
	{
		adjacent<Cell> ret;
		ret.push_back(this);
		return ret; //this operation will be virtualized by proper algorithm	
	}
	
	
		
	void SwapElement(void * pa, void * pb, void * udata)
	{
		(void) udata;
		Element ** a = static_cast<Element **>(pa);
		Element ** b = static_cast<Element **>(pb);
				  
		Element * c = *b;
		int lid = (*b)->LocalID(), rid = (*a)->LocalID();
		(*a)->SwapData(lid);
		*b = *a;
		*a = c;
		(*b)->local_id = lid;
		(*a)->local_id = rid;
		
	}
	
	void SwapPointer(void * pa, void * pb, void *)
	{
		Element ** a = (Element **)pa;
		Element ** b = (Element **)pb;
		Element * temp = *a;
		*a = *b;
		*b = temp;
	}

	
	int CompareElement(const void * pa, const void * pb, void * udata)
	{
		(void) udata;
		return CompareElementsCCentroid(pa,pb);
	}
	
	int CompareElementPointer(const void * pa, const void * pb, void * udata)
	{
		(void) udata;
		return CompareElementsCPointer(pa,pb);
	}
	
	int CompareCoordSearch(const void * pa, void * udata)
	{
		Element ** e = (Element **)(pa);
		INMOST_DATA_REAL_TYPE * cb = static_cast<INMOST_DATA_REAL_TYPE *>(udata);
		Mesh * m = (*e)->GetMeshLink();
		Storage::real_array ca = (*e)->RealArray(m->CoordsTag());
		Storage::real epsilon = m->GetEpsilon();
		unsigned int dim = m->GetDimensions();
		for(unsigned i = 0; i < dim; i++)
			if( fabs(ca[i]-cb[i]) > epsilon )
			{
				if( ca[i] > cb[i] ) return 1;
				else return -1;
			}
		return 0;
	}
	
	int CompareGIDSearch(const void * pa, void * udata)
	{
		Element ** e = (Element **)(pa);
		return (*e)->GlobalID() - *(static_cast<Storage::integer *>(udata));
	}
	
	bool Element::CheckElementConnectivity()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			for(adj_iterator jt = high_conn.begin(); jt != high_conn.end(); jt++) //iterate over upper adjacent
			{
				bool found = false;
				for(adj_iterator kt = (*jt)->low_conn.begin(); kt != (*jt)->low_conn.end(); kt++) //search for the link to me
					if( *kt == this )
					{
						found = true;
						break;
					}
				if( !found ) return false;
			}
			for(adj_iterator jt = low_conn.begin(); jt != low_conn.end(); jt++) //iterate over lower adjacent
			{
				bool found = false;
				for(adj_iterator kt = (*jt)->high_conn.begin(); kt != (*jt)->high_conn.end(); kt++) //search for the link to me
					if( *kt == this )
					{
						found = true;
						break;
					}
				if( !found ) return false;
			}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(adj_iterator jt = high_conn.begin(); jt != high_conn.end(); jt++) if( !(*jt)->GetMarker(hm) )//iterate over upper adjacent
			{
				bool found = false;
				for(adj_iterator kt = (*jt)->low_conn.begin(); kt != (*jt)->low_conn.end(); kt++) if( !(*kt)->GetMarker(hm) )//search for the link to me
					if( *kt == this )
					{
						found = true;
						break;
					}
				if( !found ) return false;
			}
			for(adj_iterator jt = low_conn.begin(); jt != low_conn.end(); jt++) if( !(*jt)->GetMarker(hm) ) //iterate over lower adjacent
			{
				bool found = false;
				for(adj_iterator kt = (*jt)->high_conn.begin(); kt != (*jt)->high_conn.end(); kt++) if( !(*kt)->GetMarker(hm) ) //search for the link to me
					if( *kt == this )
					{
						found = true;
						break;
					}
				if( !found ) return false;
			}
		}
		return true;
	}
	
	bool Element::CheckConnectivity(Mesh * m)
	{
		for(Mesh::iteratorElement it = m->BeginElement(CELL | FACE | EDGE | NODE); it != m->EndElement(); it++)
			if( !it->CheckElementConnectivity() ) return false;
		return true;
	}
	
	
	
	adjacent<Element> Element::BridgeAdjacencies(ElementType Bridge, ElementType Dest, MIDType mask)
	{
		Mesh * m = GetMeshLink();
		MIDType mrk = m->CreateMarker();
		adjacent<Element> adjcells;
		adjacent<Element> adjfaces = getAdjElements(Bridge);
		adjacent<Element> my = Bridge & GetElementType() ? adjacent<Element>() : getAdjElements(Dest);
		for(adjacent<Element>::iterator it = my.begin(); it != my.end(); it++) it->SetMarker(mrk);
		for(adjacent<Element>::iterator it = adjfaces.begin(); it != adjfaces.end(); it++)
		{
			adjacent<Element> sub = it->getAdjElements(Dest);
			for(adjacent<Element>::iterator jt = sub.begin(); jt != sub.end(); jt++)
				if( !jt->GetMarker(mrk) && (mask == 0 || jt->GetMarker(mask)) )
				{
					adjcells.push_back(&*jt);
					jt->SetMarker(mrk);
				}
		}
		for(adjacent<Element>::iterator it = adjcells.begin(); it != adjcells.end(); it++) it->RemMarker(mrk);
		for(adjacent<Element>::iterator it = my.begin(); it != my.end(); it++) it->RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return adjcells;
	}


	Node * Element::getAsNode() { if (this != NULL && GetElementType() == NODE) return static_cast<Node *>(this); return NULL; }
	Edge * Element::getAsEdge() { if (this != NULL && GetElementType() == EDGE) return static_cast<Edge *>(this); return NULL; }
	Face * Element::getAsFace() { if( this != NULL && GetElementType() == FACE) return static_cast<Face *>(this); return NULL; }
	Cell * Element::getAsCell() { if (this != NULL && GetElementType() == CELL) return static_cast<Cell *>(this); return NULL; }
	Element::Status Element::GetStatus()              {Tag s = GetMeshLink()->SharedTag(); if( !s.isValid() ) return Owned; return BulkDF(s);}
	void   Element::SetStatus(Status status) {Tag s = GetMeshLink()->SharedTag(); if( !s.isValid() ) throw TagNotInitialized; BulkDF(s) = status;}
	Storage::integer & Element::GlobalID() {return IntegerDF(GetMeshLink()->GlobalIDTag());}
	void Element::Centroid  (Storage::real * cnt) {m_link->GetGeometricData(this,CENTROID,cnt);}
	void Element::Barycenter(Storage::real * cnt) {m_link->GetGeometricData(this,BARYCENTER,cnt);}
	
}
#endif
