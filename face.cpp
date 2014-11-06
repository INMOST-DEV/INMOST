#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{
	
	Face::Face(Mesh * m)
	:Element(m,FACE)
	{
	}

	Face::Face(Mesh * m,INMOST_DATA_ENUM_TYPE lid, const Face & other)
	:Element(m,lid,other)
	{
	}
	
	Face::Face(const Face & other)
	:Element(other)
	{
	}
	
	Face & Face::operator =(Face const & other)
	{
		Element::operator =(other);
		return *this;
	}
	
	Face::~Face()
	{
	}
	
	Cell * Face::BackCell() const
	{
		if( !GetMeshLink()->HideMarker() )
		{
			if( high_conn.size() > 0 ) 
				return high_conn[0]->getAsCell(); 
			return NULL;
		}
		else
		{
			if( !high_conn.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MIDType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(&high_conn[0],high_conn.size(),i,hm);
				if( i != high_conn.size() ) return high_conn[i]->getAsCell();
			}
			return NULL;
		}
	}
	Cell * Face::FrontCell() const
	{
		if( !GetMeshLink()->HideMarker() )
		{
			if( high_conn.size() > 1 ) 
				return high_conn[1]->getAsCell(); 
			return NULL;
		}
		else
		{
			if( !high_conn.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MIDType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(&high_conn[0],high_conn.size(),i,hm); //found first
				i = Mesh::getNext(&high_conn[0],high_conn.size(),i,hm); //found second
				if( i != high_conn.size() ) return high_conn[i]->getAsCell();
			}
			return NULL;
		}
	}

	Node * Face::getBeg() const 
	{
		assert(GetElementDimension()==1);
		if( !GetMeshLink()->HideMarker() )
		{
			if( low_conn.empty() )
				return NULL;
			return low_conn.front()->low_conn.front()->getAsNode();
		}
		else
		{
			if( !low_conn.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MIDType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				if( i != low_conn.size() ) 
				{
					INMOST_DATA_ENUM_TYPE j = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
					j = Mesh::getNext(&low_conn[i]->low_conn[0],low_conn[i]->low_conn.size(),j,hm);
					if( j != low_conn[i]->low_conn.size() ) return low_conn[i]->low_conn[j]->getAsNode();
				}
			}
			return NULL;
		}
	}
	Node * Face::getEnd() const 
	{
		assert(GetElementDimension()==1);
		if( !GetMeshLink()->HideMarker() )
		{
			if( low_conn.size() < 2 )
				return NULL;
			return low_conn.back()->low_conn.front()->getAsNode();
		}
		else
		{
			if( !low_conn.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MIDType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				if( i != low_conn.size() ) 
				{
					INMOST_DATA_ENUM_TYPE j = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
					j = Mesh::getNext(&low_conn[i]->low_conn[0],low_conn[i]->low_conn.size(),j,hm);
					if( j != low_conn[i]->low_conn.size() ) return low_conn[i]->low_conn[j]->getAsNode();
				}
			}
			return NULL;
		}
	}
	
	bool Face::FaceOrientedOutside(Cell * c) const
	{
		if( BackCell() == c ) 
			return true; 
		return false;
	}
	
	void Face::ReorderEdges()
	{
		for(unsigned int j = 0; j < low_conn.size()/2; j++) // reorder edges!
		{
			Element * t = low_conn[j];
			low_conn[j] = low_conn[low_conn.size()-1-j];
			low_conn[low_conn.size()-1-j] = t;
		}
	}
	
	bool Face::CheckEdgeOrder()
	{
		return true;
	}
	
	bool Face::FixEdgeOrder()
	{
		return false;
	}
	
	
	adjacent<Node> Face::getNodes()
	{
		adjacent<Node> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			if( GetElementDimension() == 1 ) // This face is 2d edge
			{
				aret.reserve(low_conn.size());
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) //iterate over edges that are of type Vertex
					aret.push_back((*it)->low_conn.front());
			}
			else
			{
				aret.reserve(low_conn.size());
				Element * q = low_conn[0]; //edge 0
				aret.push_back(q->low_conn[0]); //node 0
				aret.push_back(q->low_conn[1]); //node 1
				Element * r = low_conn[1]; //edge 1
				if( aret.data()[0] == r->low_conn[0] || aret.data()[0] == r->low_conn[1] )
				{
					Node * temp = aret.data()[0];
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
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 1 ) // This face is 2d edge
			{
				aret.reserve(low_conn.size());
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) //iterate over edges that are of type Vertex
					if( !(*it)->GetMarker(hm) ) 
					{
						INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
						i = Mesh::getNext(&(*it)->low_conn[0],low_conn.size(),i,hm);
						if( i < (*it)->low_conn.size() ) aret.push_back((*it)->low_conn[i]);
					}
			}
			else
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
					Node * temp = aret.data()[0];
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

		}
		return aret;
	}


	adjacent<Node> Face::getNodes(MIDType mask, bool invert)
	{
		adjacent<Node> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			if( GetElementDimension() == 1 ) // This face is 2d edge
			{
				aret.reserve(low_conn.size());
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) //iterate over edges that are of type Vertex
					if( invert ^ (*it)->low_conn.front()->GetMarker(mask) )
						aret.push_back((*it)->low_conn.front());
			}
			else
			{
				aret.reserve(low_conn.size());
				Element * q = low_conn[0], * first = q->low_conn[0], * last = q->low_conn[1]; //edge 0
				if( invert ^ q->low_conn[0]->GetMarker(mask) ) aret.push_back(q->low_conn[0]); //node 0
				if( invert ^ q->low_conn[1]->GetMarker(mask) ) aret.push_back(q->low_conn[1]); //node 1
				Element * r = low_conn[1]; //edge 1
				if( first == r->low_conn[0] || first == r->low_conn[1] ) 
				{
					last = first;
					if( aret.size() > 1 )
					{
						Node * temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_iterator it = low_conn.begin()+1, iend = low_conn.end()-1;
				while(it != iend) //loop over edges
				{
					if( last == (*it)->low_conn[0] ) 
					{
						if( invert ^ (*it)->low_conn[1]->GetMarker(mask) ) 
							aret.push_back((*it)->low_conn[1]);
						last = (*it)->low_conn[1];
					}
					else 
					{
						if( invert ^ (*it)->low_conn[0]->GetMarker(mask) )
							aret.push_back((*it)->low_conn[0]);
						last = (*it)->low_conn[0];
					}
					++it;
				}
			}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 1 ) // This face is 2d edge
			{
				aret.reserve(low_conn.size());
				for(Element::adj_iterator it = low_conn.begin(); it != low_conn.end(); it++) //iterate over edges that are of type Vertex
					if( !(*it)->GetMarker(hm) ) 
					{
						INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
						i = Mesh::getNext(&(*it)->low_conn[0],low_conn.size(),i,hm);
						if( i < (*it)->low_conn.size() && (invert ^ (*it)->low_conn[i]->GetMarker(mask)) ) aret.push_back((*it)->low_conn[i]);
					}
			}
			else
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1),
					k = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1), k2;
				aret.reserve(low_conn.size());
				i = Mesh::getNext(&low_conn[0],low_conn.size(),i,hm);
				Element * q = low_conn[i], * first, * last; //edge 0
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
						Node * temp = aret.data()[0];
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
						if( invert ^ (*it)->low_conn[k2]->GetMarker(mask) ) aret.push_back((*it)->low_conn[k2]);
						last = (*it)->low_conn[k2];
					}
					else 
					{
						if( invert ^ (*it)->low_conn[k1]->GetMarker(mask) ) aret.push_back((*it)->low_conn[k1]);
						last = (*it)->low_conn[k1];
					}
					++it;
				}
			}

		}
		return aret;
	}

	adjacent<Edge> Face::getEdges()
	{
		if( !GetMeshLink()->HideMarker() )
			return adjacent<Edge>(low_conn.begin(),low_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Edge> aret;
			for(adj_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}

	adjacent<Edge> Face::getEdges(MIDType mask, bool invert)
	{
		adjacent<Edge> aret;
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

	adjacent<Cell> Face::getCells()
	{
		if( !GetMeshLink()->HideMarker() )
			return adjacent<Cell>(high_conn.begin(),high_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Cell> aret;
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}

	adjacent<Cell> Face::getCells(MIDType mask, bool invert)
	{
		adjacent<Cell> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( invert ^ (*it)->GetMarker(mask)) aret.push_back((*it));
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( (invert ^ (*it)->GetMarker(mask)) && !(*it)->GetMarker(hm) ) aret.push_back((*it));
		}
		return aret;
	}

	Storage::real Face::Area() {Storage::real ret; m_link->GetGeometricData(this,MEASURE,&ret); return ret;}
	void Face::Normal(Storage::real * nrm) {m_link->GetGeometricData(this,NORMAL,nrm);}
}
#endif
