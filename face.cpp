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
			adj_type const & hc = HighConn();
			if( hc.size() > 0 ) 
				return hc[0]->getAsCell(); 
			return NULL;
		}
		else
		{
			adj_type const & hc = HighConn();
			if( !hc.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MarkerType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(hc.data(),hc.size(),i,hm);
				if( i != hc.size() ) return hc[i]->getAsCell();
			}
			return NULL;
		}
	}
	Cell * Face::FrontCell() const
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			if( hc.size() > 1 ) 
				return hc[1]->getAsCell(); 
			return NULL;
		}
		else
		{
			adj_type const & hc = HighConn();
			if( !hc.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MarkerType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(&hc[0],hc.size(),i,hm); //found first
				i = Mesh::getNext(&hc[0],hc.size(),i,hm); //found second
				if( i != hc.size() ) return hc[i]->getAsCell();
			}
			return NULL;
		}
	}

	Node * Face::getBeg() const 
	{
		assert(GetElementDimension()==1);
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & lc = LowConn();
			if( lc.empty() )
				return NULL;
			return lc.front()->LowConn().front()->getAsNode();
		}
		else
		{
			adj_type const & lc = LowConn();
			if( !lc.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MarkerType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				if( i != lc.size() ) 
				{
					adj_type const & llc = lc[i]->LowConn();
					INMOST_DATA_ENUM_TYPE j = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
					j = Mesh::getNext(llc.data(),llc.size(),j,hm);
					if( j != llc.size() ) return llc[j]->getAsNode();
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
			adj_type const & lc = LowConn();
			if( lc.size() < 2 )
				return NULL;
			return lc.back()->LowConn().front()->getAsNode();
		}
		else
		{
			adj_type const & lc = LowConn();
			if( !lc.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MarkerType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				if( i != lc.size() ) 
				{
					adj_type const & llc = lc[i]->LowConn();
					INMOST_DATA_ENUM_TYPE j = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
					j = Mesh::getNext(llc.data(),llc.size(),j,hm);
					if( j != llc.size() ) return llc[j]->getAsNode();
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
		adj_type & lc = LowConn();
		for(unsigned int j = 0; j < lc.size()/2; j++) // reorder edges!
		{
			Element * t = lc[j];
			lc[j] = lc[lc.size()-1-j];
			lc[lc.size()-1-j] = t;
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
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				for(adj_type::enumerator it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
					aret.push_back(lc[it]->LowConn().front());
			}
			else
			{
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				Element * q = lc[0]; //edge 0
				adj_type & qlc = q->LowConn();
				aret.push_back(qlc[0]); //node 0
				aret.push_back(qlc[1]); //node 1
				Element * r = lc[1]; //edge 1
				adj_type & rlc = r->LowConn();
				if( aret.data()[0] == rlc[0] || aret.data()[0] == rlc[1] )
				{
					Node * temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::enumerator it = 1, iend = lc.size()-1;
				while(it < iend) //loop over edges
				{
					adj_type & ilc = lc[it]->LowConn();
					if( &aret.back() == ilc[0] ) aret.push_back(ilc[1]);
					else aret.push_back(ilc[0]);
					++it;
				}
			}
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 1 ) // This face is 2d edge
			{
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				for(adj_type::enumerator it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
					if( !lc[it]->GetMarker(hm) ) 
					{
						adj_type & ilc = lc[it]->LowConn();
						INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
						i = Mesh::getNext(ilc.data(),ilc.size(),i,hm);
						if( i < ilc.size() ) aret.push_back(ilc[i]);
					}
			}
			else
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1),
					k = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1), k2;
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				Element * q = lc[i]; //edge 0
				adj_type & qlc = q->LowConn();
				k = Mesh::getNext(qlc.data(),qlc.size(),k,hm);
				aret.push_back(qlc[k]); //node 0
				k = Mesh::getNext(qlc.data(),qlc.size(),k,hm);
				aret.push_back(qlc[k]); //node 1
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				Element * r = lc[i]; //edge 1
				adj_type & rlc = r->LowConn();
				k1 = Mesh::getNext(rlc.data(),rlc.size(),k1,hm);
				k2 = Mesh::getNext(rlc.data(),rlc.size(),k1,hm);
				if( aret.data()[0] == rlc[k1] || aret.data()[0] == rlc[k2] )
				{
					Node * temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::enumerator it = 1, iend = lc.size()-1;
				while(it != iend) if( !lc[it]->GetMarker(hm) ) //loop over edges
				{
					adj_type & ilc = lc[it]->LowConn();
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1); 
					k1 = Mesh::getNext(ilc.data(),ilc.size(),k1,hm);
					k2 = Mesh::getNext(ilc.data(),ilc.size(),k1,hm);
					if( &aret.back() == ilc[k1] ) aret.push_back(ilc[k2]);
					else aret.push_back(ilc[k1]);
					++it;
				}
			}

		}
		return aret;
	}


	adjacent<Node> Face::getNodes(MarkerType mask, bool invert)
	{
		adjacent<Node> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			if( GetElementDimension() == 1 ) // This face is 2d edge
			{
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				for(adj_type::enumerator it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
				{
					Element * e = lc[it]->LowConn().front();
					if( invert ^ e->GetMarker(mask) ) aret.push_back(e);
				}
			}
			else
			{
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				Element * q = lc[0];
				adj_type & qlc = q->LowConn();
				Element * first = qlc[0], * last = qlc[1]; //edge 0
				if( invert ^ qlc[0]->GetMarker(mask) ) aret.push_back(qlc[0]); //node 0
				if( invert ^ qlc[1]->GetMarker(mask) ) aret.push_back(qlc[1]); //node 1
				Element * r = lc[1]; //edge 1
				adj_type & rlc = r->LowConn();
				if( first == rlc[0] || first == rlc[1] ) 
				{
					last = first;
					if( aret.size() > 1 )
					{
						Node * temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_type::enumerator it = 1, iend = lc.size()-1;
				while(it < iend) //loop over edges
				{
					adj_type & ilc = lc[it]->LowConn();
					if( last == ilc[0] ) 
					{
						if( invert ^ ilc[1]->GetMarker(mask) ) 
							aret.push_back(ilc[1]);
						last = ilc[1];
					}
					else 
					{
						if( invert ^ ilc[0]->GetMarker(mask) )
							aret.push_back(ilc[0]);
						last = ilc[0];
					}
					++it;
				}
			}
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 1 ) // This face is 2d edge
			{
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				for(adj_type::enumerator it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
					if( !lc[it]->GetMarker(hm) ) 
					{
						adj_type & ilc = lc[it]->LowConn();
						INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
						i = Mesh::getNext(ilc.data(),ilc.size(),i,hm);
						if( i < ilc.size() && (invert ^ ilc[i]->GetMarker(mask)) ) aret.push_back(ilc[i]);
					}
			}
			else
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1),
					k = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1), k2;
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				Element * q = lc[i], * first, * last; //edge 0
				adj_type & qlc = q->LowConn();
				k = Mesh::getNext(qlc.data(),qlc.size(),k,hm);
				if( invert ^ qlc[k]->GetMarker(mask) ) aret.push_back(qlc[k]); //node 0
				first = qlc[k];
				k = Mesh::getNext(qlc.data(),qlc.size(),k,hm);
				if( invert ^ qlc[k]->GetMarker(mask) ) aret.push_back(qlc[k]); //node 1
				last = qlc[k];
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				Element * r = lc[i]; //edge 1
				adj_type & rlc = r->LowConn();
				k1 = Mesh::getNext(rlc.data(),rlc.size(),k1,hm);
				k2 = Mesh::getNext(rlc.data(),rlc.size(),k1,hm);
				if( first == rlc[k1] || first == rlc[k2] )
				{
					last = first;
					if( aret.size() > 1 )
					{
						Node * temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_type::enumerator it = 1, iend = lc.size()-1;
				while(it < iend) if( !lc[it]->GetMarker(hm) ) //loop over edges
				{
					adj_type & ilc = lc[it]->LowConn();
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1); 
					k1 = Mesh::getNext(ilc.data(),ilc.size(),k1,hm);
					k2 = Mesh::getNext(ilc.data(),ilc.size(),k1,hm);
					if( last == ilc[k1] ) 
					{
						if( invert ^ ilc[k2]->GetMarker(mask) ) aret.push_back(ilc[k2]);
						last = ilc[k2];
					}
					else 
					{
						if( invert ^ ilc[k1]->GetMarker(mask) ) aret.push_back(ilc[k1]);
						last = ilc[k1];
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
		{
			adj_type & lc = LowConn();
			return adjacent<Edge>(lc.data(),lc.data()+lc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Edge> aret;
			adj_type & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( !lc[it]->GetMarker(hm) ) aret.push_back(lc[it]);
			return aret;
		}
	}

	adjacent<Edge> Face::getEdges(MarkerType mask, bool invert)
	{
		adjacent<Edge> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( invert ^ lc[it]->GetMarker(mask) ) aret.push_back(lc[it]);
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adj_type & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( (invert ^ lc[it]->GetMarker(mask)) && !lc[it]->GetMarker(hm) ) aret.push_back(lc[it]);
		}
		return aret;
	}

	adjacent<Cell> Face::getCells()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type & hc = HighConn();
			return adjacent<Cell>(hc.data(),hc.data()+hc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Cell> aret;
			adj_type & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( !hc[it]->GetMarker(hm) ) aret.push_back(hc[it]);
			return aret;
		}
	}

	adjacent<Cell> Face::getCells(MarkerType mask, bool invert)
	{
		adjacent<Cell> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( invert ^ hc[it]->GetMarker(mask)) aret.push_back(hc[it]);
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adj_type & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( (invert ^ hc[it]->GetMarker(mask)) && !hc[it]->GetMarker(hm) ) aret.push_back(hc[it]);
		}
		return aret;
	}

	Storage::real Face::Area() {Storage::real ret; m_link->GetGeometricData(this,MEASURE,&ret); return ret;}
	void Face::Normal(Storage::real * nrm) {m_link->GetGeometricData(this,NORMAL,nrm);}
}
#endif
