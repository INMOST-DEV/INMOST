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
	
	Cell * Face::BackCell()  
	{
		if( !GetMeshLink()->HideMarker() )
		{
			if( high_conn.size() > 0 ) 
				return high_conn[0]->getAsCell(); 
			return NULL;
		}
		else
		{
			INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
			MIDType hm = GetMeshLink()->HideMarker();
			i = Mesh::getNext(&high_conn[0],high_conn.size(),i,hm);
			if( i != high_conn.size() ) return high_conn[i]->getAsCell();
			return NULL;
		}
	}
	Cell * Face::FrontCell() 
	{
		if( !GetMeshLink()->HideMarker() )
		{
			if( high_conn.size() > 1 ) 
				return high_conn[1]->getAsCell(); 
			return NULL;
		}
		else
		{
			INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
			MIDType hm = GetMeshLink()->HideMarker();
			i = Mesh::getNext(&high_conn[0],high_conn.size(),i,hm); //found first
			i = Mesh::getNext(&high_conn[0],high_conn.size(),i,hm); //found second
			if( i != high_conn.size() ) return high_conn[i]->getAsCell();
			return NULL;
		}
	}
	
	bool Face::FaceOrientedOutside(Cell * c) 
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
					if( !(*it)->GetMarker(hm) ) aret.push_back((*it)->low_conn.front());
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
	Storage::real Face::Area() {Storage::real ret; m_link->GetGeometricData(this,MEASURE,&ret); return ret;}
	void Face::Normal(Storage::real * nrm) {m_link->GetGeometricData(this,NORMAL,nrm);}
}
#endif
