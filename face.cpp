#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{
	
	
	
	Cell Face::BackCell() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			if( hc.size() > 0 ) 
				return Cell(m,hc[0]); 
			return Cell(m,InvalidHandle());
		}
		else
		{
			adj_type const & hc = m->HighConn(GetHandle());
			if( !hc.empty() )
			{
				enumerator i = ENUMUNDEF;
				MarkerType hm = m->HideMarker();
				i = m->getNext(hc.data(),static_cast<enumerator>(hc.size()),i,hm);
				if( i != hc.size() ) return Cell(m,hc[i]);
			}
			return Cell(m,InvalidHandle());
		}
	}
	Cell Face::FrontCell() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			if( hc.size() > 1 ) 
				return Cell(m,hc[1]); 
			return Cell(m,InvalidHandle());
		}
		else
		{
			adj_type const & hc = m->HighConn(GetHandle());
			if( !hc.empty() )
			{
				enumerator i = ENUMUNDEF;
				MarkerType hm = m->HideMarker();
				i = m->getNext(hc.data(),static_cast<enumerator>(hc.size()),i,hm); //found first
				i = m->getNext(hc.data(),static_cast<enumerator>(hc.size()),i,hm); //found second
				if( i != hc.size() ) return Cell(m,hc[i]);
			}
			return Cell(m,InvalidHandle());
		}
	}

	Node Face::getBeg() const 
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		assert(GetElementDimension()==1);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & lc = m->LowConn(GetHandle());
			if( lc.empty() )
				return Node(m,InvalidHandle());
			return Node(m,m->LowConn(lc.front()).front());
		}
		else
		{
			adj_type const & lc = m->LowConn(GetHandle());
			if( !lc.empty() )
			{
				enumerator i = ENUMUNDEF;
				MarkerType hm = m->HideMarker();
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				if( i != static_cast<enumerator>(lc.size()) ) 
				{
					adj_type const & llc = m->LowConn(lc[i]);
					enumerator j = ENUMUNDEF;
					j = m->getNext(llc.data(),static_cast<enumerator>(llc.size()),j,hm);
					if( j != llc.size() ) return Node(m,llc[j]);
				}
			}
			return Node(m,InvalidHandle());
		}
	}
	Node Face::getEnd() const 
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		assert(GetElementDimension()==1);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & lc = m->LowConn(GetHandle());
			if( lc.size() < 2 )
				return Node(m,InvalidHandle());
			return Node(m,m->LowConn(lc.back()).front());
		}
		else
		{
			adj_type const & lc = m->LowConn(GetHandle());
			if( !lc.empty() )
			{
				enumerator i = ENUMUNDEF;
				MarkerType hm = m->HideMarker();
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				if( i != static_cast<enumerator>(lc.size()) ) 
				{
					adj_type const & llc = m->LowConn(lc[i]);
					enumerator j = ENUMUNDEF;
					j = m->getNext(llc.data(),static_cast<enumerator>(llc.size()),j,hm);
					if( j != static_cast<enumerator>(llc.size()) ) return Node(m,llc[j]);
				}
			}
			return Node(m,InvalidHandle());
		}
	}
	
	bool Face::FaceOrientedOutside(Cell c) const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		if( BackCell() == c ) 
			return true; 
		return false;
	}
	
	void Face::ReorderEdges() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		adj_type & lc = GetMeshLink()->LowConn(GetHandle());
		for(adj_type::size_type j = 0; j < lc.size()/2; j++) // reorder edges!
		{
			HandleType t = lc[j];
			lc[j] = lc[lc.size()-1-j];
			lc[lc.size()-1-j] = t;
		}
	}
	
	bool Face::CheckEdgeOrder() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		return true;
	}
	
	bool Face::FixEdgeOrder() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		return false;
	}
	
	
	ElementArray<Node> Face::getNodes() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		ElementArray<Node> aret(m);
		if( !m->HideMarker() )
		{
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 1 ) // This face is 2d edge
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				for(adj_type::size_type it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
					aret.push_back(m->LowConn(lc[it]).front());
			}
			else
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				HandleType q = lc[0]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				aret.push_back(qlc[0]); //node 0
				aret.push_back(qlc[1]); //node 1
				HandleType r = lc[1]; //edge 1
				adj_type & rlc = m->LowConn(r);
				if( aret.data()[0] == rlc[0] || aret.data()[0] == rlc[1] )
				{
					HandleType temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it < iend) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					if( aret.atback() == ilc[0] ) aret.push_back(ilc[1]);
					else aret.push_back(ilc[0]);
					++it;
				}
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 1 ) // This face is 2d edge
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				for(adj_type::size_type it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
					if( !m->GetMarker(lc[it],hm) ) 
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						enumerator i = ENUMUNDEF;
						i = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),i,hm);
						if( i < static_cast<enumerator>(ilc.size()) ) aret.push_back(ilc[i]);
					}
			}
			else
			{
				enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType q = lc[i]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				aret.push_back(qlc[k]); //node 0
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				aret.push_back(qlc[k]); //node 1
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType r = lc[i]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				k1 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				k2 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				if( aret.data()[0] == rlc[k1] || aret.data()[0] == rlc[k2] )
				{
					HandleType temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it != iend) if( !m->GetMarker(lc[it],hm) ) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					k1 = ENUMUNDEF; 
					k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
					k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
					if( aret.atback() == ilc[k1] ) 
						aret.push_back(ilc[k2]);
					else 
						aret.push_back(ilc[k1]);
					++it;
				}
			}

		}
		return aret;
	}


	ElementArray<Node> Face::getNodes(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		ElementArray<Node> aret(m);
		if( !m->HideMarker() )
		{
			if(  Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 1 ) // This face is 2d edge
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				for(adj_type::size_type it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
				{
					HandleType e = m->LowConn(lc[it]).front();
					if( invert ^ m->GetMarker(e,mask) ) aret.push_back(e);
				}
			}
			else
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				HandleType q = lc[0];
				adj_type const & qlc = m->LowConn(q);
				HandleType first = qlc[0], last = qlc[1]; //edge 0
				if( invert ^ m->GetMarker(qlc[0],mask) ) aret.push_back(qlc[0]); //node 0
				if( invert ^ m->GetMarker(qlc[1],mask) ) aret.push_back(qlc[1]); //node 1
				HandleType r = lc[1]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				if( first == rlc[0] || first == rlc[1] ) 
				{
					last = first;
					if( aret.size() > 1 )
					{
						HandleType temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it < iend) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					if( last == ilc[0] ) 
					{
						if( invert ^ m->GetMarker(ilc[1],mask) ) 
							aret.push_back(ilc[1]);
						last = ilc[1];
					}
					else 
					{
						if( invert ^ m->GetMarker(ilc[0],mask) )
							aret.push_back(ilc[0]);
						last = ilc[0];
					}
					++it;
				}
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 1 ) // This face is 2d edge
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				for(adj_type::size_type it = 0; it < lc.size(); it++) //iterate over edges that are of type Vertex
					if( !m->GetMarker(lc[it],hm) ) 
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						enumerator i = ENUMUNDEF;
						i = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),i,hm);
						if( i < static_cast<enumerator>(ilc.size()) && (invert ^ m->GetMarker(ilc[i],mask)) ) 
							aret.push_back(ilc[i]);
					}
			}
			else
			{
				enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType q = lc[i], first, last; //edge 0
				adj_type const & qlc = m->LowConn(q);
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				if( invert ^ m->GetMarker(qlc[k],mask) ) aret.push_back(qlc[k]); //node 0
				first = qlc[k];
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				if( invert ^ m->GetMarker(qlc[k],mask) ) aret.push_back(qlc[k]); //node 1
				last = qlc[k];
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType r = lc[i]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				k1 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				k2 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				if( first == rlc[k1] || first == rlc[k2] )
				{
					last = first;
					if( aret.size() > 1 )
					{
						HandleType temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it < iend) if( !m->GetMarker(lc[it],hm) ) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					k1 = ENUMUNDEF; 
					k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
					k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
					if( last == ilc[k1] ) 
					{
						if( invert ^ m->GetMarker(ilc[k2],mask) ) 
							aret.push_back(ilc[k2]);
						last = ilc[k2];
					}
					else 
					{
						if( invert ^ m->GetMarker(ilc[k1],mask) ) 
							aret.push_back(ilc[k1]);
						last = ilc[k1];
					}
					++it;
				}
			}

		}
		return aret;
	}

	ElementArray<Edge> Face::getEdges() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & lc = m->LowConn(GetHandle());
			return ElementArray<Edge>(m,lc.data(),lc.data()+lc.size());
		}
		else
		{
			MarkerType hm = m->HideMarker();
			ElementArray<Edge> aret(m);
			adj_type const & lc = m->LowConn(GetHandle());
			for(adj_type::size_type it = 0; it < lc.size(); ++it)
				if( !m->GetMarker(lc[it],hm) ) 
					aret.push_back(lc[it]);
			return aret;
		}
	}

	ElementArray<Edge> Face::getEdges(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		ElementArray<Edge> aret(m);
		if( !m->HideMarker() )
		{
			adj_type const & lc = m->LowConn(GetHandle());
			for(adj_type::size_type it = 0; it < lc.size(); ++it)
				if( invert ^ m->GetMarker(lc[it],mask) ) 
					aret.push_back(lc[it]);
		}
		else
		{
			MarkerType hm = m->HideMarker();
			adj_type const & lc = m->LowConn(GetHandle());
			for(adj_type::size_type it = 0; it < lc.size(); ++it)
				if( (invert ^ m->GetMarker(lc[it],mask)) && !m->GetMarker(lc[it],hm) ) 
					aret.push_back(lc[it]);
		}
		return aret;
	}

	ElementArray<Cell> Face::getCells() const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			return ElementArray<Cell>(m,hc.data(),hc.data()+hc.size());
		}
		else
		{
			MarkerType hm = m->HideMarker();
			ElementArray<Cell> aret(m);
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( !m->GetMarker(hc[it],hm) ) 
					aret.push_back(hc[it]);
			return aret;
		}
	}

	ElementArray<Cell> Face::getCells(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==FACE);
		Mesh * m = GetMeshLink();
		ElementArray<Cell> aret(m);
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( invert ^ m->GetMarker(hc[it],mask)) 
					aret.push_back(hc[it]);
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( (invert ^ m->GetMarker(hc[it],mask)) && !m->GetMarker(hc[it],hm) ) 
					aret.push_back(hc[it]);
		}
		return aret;
	}

	
}
#endif
