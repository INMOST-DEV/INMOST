#include "inmost.h"
#if defined(USE_MESH)

namespace INMOST
{
	
	Cell::Cell(Mesh * m)
	:Element(m,CELL)
	{
	}

	Cell::Cell(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Cell & other)
	:Element(m,lid,other)
	{
	}
	
	Cell::Cell(const Cell & other)
	:Element(other)
	{
	}
	
	Cell & Cell::operator =(Cell const & other)
	{
        Element::operator =(other);
        return *this;
	}
	
	Cell::~Cell()
	{
	}
	
	bool Cell::CheckEdgeOrder()
	{
		return true;
	}
	
	bool Cell::FixEdgeOrder()
	{
		return false;
	}
	
	adjacent<Node> Cell::getNodes()
	{
		if( !GetMeshLink()->HideMarker() )
			return adjacent<Node>(high_conn.begin(),high_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Node> aret;
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}


	adjacent<Node> Cell::getNodes(MIDType mask, bool invert)
	{
		adjacent<Node> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( invert ^ (*it)->GetMarker(mask) ) aret.push_back((*it));
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( (invert ^ (*it)->GetMarker(mask)) && !(*it)->GetMarker(hm) ) aret.push_back((*it));
		}
		return aret;
	}

	adjacent<Edge> Cell::getEdges()
	{
		adjacent<Edge> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			if( GetElementDimension() == 2 ) // This cell is 2d face
			{
				aret.reserve(low_conn.size());
				Element * q = low_conn[0]; //edge 0
				aret.push_back(q->low_conn[0]); //node 0
				aret.push_back(q->low_conn[1]); //node 1
				Element * r = low_conn[1]; //edge 1
				if( aret.data()[0] == r->low_conn[0] || aret.data()[0] == r->low_conn[1] )
				{
					Edge * temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_iterator it = low_conn.begin()+1, iend = low_conn.end()-1;
						while(it != iend) //loop over edges
						{
							if( &aret.back() == (*it)->low_conn[0] ) aret.push_back((*it)->low_conn[1]);
							else aret.push_back((*it)->low_conn[0]);
							++it;
						}
			}
			else
			{
				Mesh * m = GetMeshLink();
				MIDType mrk = m->CreateMarker();
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) //faces
					for(Element::adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) //edges
						if( !(*jt)->GetMarker(mrk))
						{
							aret.push_back(*jt);
							(*jt)->SetMarker(mrk);
						}
				for(adjacent<Edge>::iterator it = aret.begin(); it != aret.end(); it++)
					it->RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 2 ) // This cell is 2d face
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1), k2;
				aret.reserve(low_conn.size());
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				Element * q = low_conn[i]; //edge 0
				k = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k,hm);
				aret.push_back(q->low_conn[k]); //node 0
				k = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k,hm);
				aret.push_back(q->low_conn[k]); //node 1
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				Element * r = low_conn[i]; //edge 1
				k1 = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k1,hm);
				k2 = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k1,hm);
				if( aret.data()[0] == r->low_conn[k1] || aret.data()[0] == r->low_conn[k2] )
				{
					Edge * temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_iterator it = low_conn.begin()+1, iend = low_conn.end()-1;
				while(it != iend) if( !(*it)->GetMarker(hm) ) //loop over edges
				{
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1); 
					k1 = Mesh::getNext(&(*it)->low_conn[0],(*it)->low_conn.size(),k1,hm);
					k2 = Mesh::getNext(&(*it)->low_conn[0],(*it)->low_conn.size(),k1,hm);
					if( &aret.back() == (*it)->low_conn[k1] ) aret.push_back((*it)->low_conn[k2]);
					else aret.push_back((*it)->low_conn[k1]);
					++it;
				}
			}
			else
			{
				Mesh * m = GetMeshLink();
				MIDType mrk = m->CreateMarker();
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) if( !(*it)->GetMarker(hm) ) //faces
					for(Element::adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) if( !(*jt)->GetMarker(hm) )//edges
						if( !(*jt)->GetMarker(mrk))
						{
							aret.push_back(*jt);
							(*jt)->SetMarker(mrk);
						}
				for(adjacent<Edge>::iterator it = aret.begin(); it != aret.end(); it++)
					it->RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		return aret;
	}


	adjacent<Edge> Cell::getEdges(MIDType mask, bool invert)
	{
		adjacent<Edge> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			if( GetElementDimension() == 2 ) // This cell is 2d face
			{
				aret.reserve(low_conn.size());
				Element * last, * first;
				Element * q = low_conn[0]; //edge 0
				if( invert ^ q->low_conn[0]->GetMarker(mask) ) aret.push_back(q->low_conn[0]); //node 0
				if( invert ^ q->low_conn[1]->GetMarker(mask) ) aret.push_back(q->low_conn[1]); //node 1
				first = q->low_conn[0];
				last  = q->low_conn[1];
				Element * r = low_conn[1]; //edge 1
				if( first == r->low_conn[0] || first == r->low_conn[1] )
				{
					last = first;
					if( aret.size() > 1 )
					{
						Edge * temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_iterator it = low_conn.begin()+1, iend = low_conn.end()-1;
						while(it != iend) //loop over edges
						{
							if( last == (*it)->low_conn[0] ) 
							{
								last = (*it)->low_conn[1];
								if( invert ^ last->GetMarker(mask) ) aret.push_back(last);
							}
							else 
							{
								last = (*it)->low_conn[0];
								if( invert ^ last->GetMarker(mask) ) aret.push_back(last);
							}
							++it;
						}
			}
			else
			{
				Mesh * m = GetMeshLink();
				MIDType mrk = m->CreateMarker();
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) //faces
					for(Element::adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) //edges
						if( (invert ^ (*jt)->GetMarker(mask)) && !(*jt)->GetMarker(mrk))
						{
							aret.push_back(*jt);
							(*jt)->SetMarker(mrk);
						}
				for(adjacent<Edge>::iterator it = aret.begin(); it != aret.end(); it++)
					it->RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 2 ) // This cell is 2d face
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1), k2;
				Element * last, * first;
				aret.reserve(low_conn.size());
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				Element * q = low_conn[i]; //edge 0
				k = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k,hm);
				if( invert ^ q->low_conn[k]->GetMarker(mask) ) aret.push_back(q->low_conn[k]); //node 0
				first = q->low_conn[k];
				k = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k,hm);
				if( invert ^ q->low_conn[k]->GetMarker(mask) ) aret.push_back(q->low_conn[k]); //node 1
				last = q->low_conn[k];
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				Element * r = low_conn[i]; //edge 1
				k1 = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k1,hm);
				k2 = Mesh::getNext(&q->low_conn[0],q->low_conn.size(),k1,hm);
				if( first == r->low_conn[k1] || first == r->low_conn[k2] )
				{
					last = first;
					if( aret.size() > 1 )
					{
						Edge * temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_iterator it = low_conn.begin()+1, iend = low_conn.end()-1;
				while(it != iend) if( !(*it)->GetMarker(hm) ) //loop over edges
				{
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1); 
					k1 = Mesh::getNext(&(*it)->low_conn[0],(*it)->low_conn.size(),k1,hm);
					k2 = Mesh::getNext(&(*it)->low_conn[0],(*it)->low_conn.size(),k1,hm);
					if( last == (*it)->low_conn[k1] ) 
					{
						last = (*it)->low_conn[k2];
						if( invert ^ last->GetMarker(mask) ) aret.push_back(last);
					}
					else 
					{
						last = (*it)->low_conn[k1];
						if( invert ^ last->GetMarker(mask) ) aret.push_back(last);
					}
					++it;
				}
			}
			else
			{
				Mesh * m = GetMeshLink();
				MIDType mrk = m->CreateMarker();
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) if( !(*it)->GetMarker(hm) ) //faces
					for(Element::adj_iterator jt = (*it)->low_conn.begin(); jt != (*it)->low_conn.end(); jt++) if( !(*jt)->GetMarker(hm) )//edges
						if( (invert ^ (*jt)->GetMarker(mask)) && !(*jt)->GetMarker(mrk))
						{
							aret.push_back(*jt);
							(*jt)->SetMarker(mrk);
						}
				for(adjacent<Edge>::iterator it = aret.begin(); it != aret.end(); it++)
					it->RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		return aret;
	}

	adjacent<Face> Cell::getFaces()
	{
		if( !GetMeshLink()->HideMarker() )
			return adjacent<Face>(low_conn.begin(),low_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Face> aret;
			for(adj_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}


	adjacent<Face> Cell::getFaces(MIDType mask, bool invert)
	{
		adjacent<Face> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			for(adj_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
				if( invert ^ (*it)->GetMarker(mask) ) aret.push_back((*it));
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(adj_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
				if( (invert ^ (*it)->GetMarker(mask)) && !(*it)->GetMarker(hm) ) aret.push_back((*it));
		}
		return aret;
	}
	
	
	
	
	Storage::real Cell::Volume() {Storage::real ret; m_link->GetGeometricData(this,MEASURE,&ret); return ret;}	
}
#endif
