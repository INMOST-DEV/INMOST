#include "inmost.h"
#if defined(USE_MESH)

namespace INMOST
{
	
	
	bool Cell::CheckEdgeOrder() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		return true;
	}
	
	bool Cell::FixEdgeOrder() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		return false;
	}
	
	ElementArray<Node> Cell::getNodes() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			return ElementArray<Node>(m,hc.data(),hc.data()+hc.size());
		}
		else
		{
			MarkerType hm = m->HideMarker();
			ElementArray<Node> aret(m);
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( !m->GetMarker(hc[it],hm) ) aret.push_back(hc[it]);
			return aret;
		}
	}


	ElementArray<Node> Cell::getNodes(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		ElementArray<Node> aret(m);
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( invert ^ m->GetMarker(hc[it],mask) ) 
					aret.push_back(hc[it]);
		}
		else
		{
			MarkerType hm = m->HideMarker();
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( (invert ^ m->GetMarker(hc[it],mask)) && !m->GetMarker(hc[it],hm) ) 
					aret.push_back(hc[it]);
		}
		return aret;
	}

	ElementArray<Edge> Cell::getEdges() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		ElementArray<Edge> aret(m);
		if( !m->HideMarker() )
		{
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				HandleType q = lc[0]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				aret.push_back(qlc[0]); //node 0
				aret.push_back(qlc[1]); //node 1
				HandleType r = lc[1]; //edge 1
				adj_type const & rlc = m->LowConn(r);
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
					if( aret.atback() == ilc[0] ) 
						aret.push_back(ilc[1]);
					else 
						aret.push_back(ilc[0]);
					++it;
				}
			}
			else
			{
				MarkerType mrk = m->CreateMarker();
				adj_type const & lc = m->LowConn(GetHandle());
				for(adj_type::size_type it = 0; it < lc.size(); it++) //faces
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) //edges
						if( !m->GetMarker(ilc[jt],mrk))
						{
							aret.push_back(ilc[jt]);
							m->SetMarker(ilc[jt],mrk);
						}
				}
				for(ElementArray<Edge>::size_type it = 0; it < aret.size(); it++)
					m->RemMarker(aret.at(it),mrk);
				m->ReleaseMarker(mrk);
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
			{
				integer i = -1, k = -1, k1 = -1, k2;
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				i = m->getNext(lc.data(),static_cast<integer>(lc.size()),i,hm);
				HandleType q = lc[i]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				k = m->getNext(qlc.data(),static_cast<integer>(qlc.size()),k,hm);
				aret.push_back(qlc[k]); //node 0
				k = m->getNext(qlc.data(),static_cast<integer>(qlc.size()),k,hm);
				aret.push_back(qlc[k]); //node 1
				i = m->getNext(lc.data(),static_cast<integer>(lc.size()),i,hm);
				HandleType r = lc[i]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				k1 = m->getNext(rlc.data(),static_cast<integer>(rlc.size()),k1,hm);
				k2 = m->getNext(rlc.data(),static_cast<integer>(rlc.size()),k1,hm);
				if( aret.data()[0] == rlc[k1] || aret.data()[0] == rlc[k2] )
				{
					HandleType temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it < iend) if( !m->GetMarker(lc[it],hm) ) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					k1 = -1; 
					k1 = m->getNext(ilc.data(),static_cast<integer>(ilc.size()),k1,hm);
					k2 = m->getNext(ilc.data(),static_cast<integer>(ilc.size()),k1,hm);
					if( aret.atback() == ilc[k1] ) 
						aret.push_back(ilc[k2]);
					else 
						aret.push_back(ilc[k1]);
					++it;
				}
			}
			else
			{
				MarkerType mrk = m->CreateMarker();
				adj_type const & lc = m->LowConn(GetHandle());
				for(adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) ) //faces
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !m->GetMarker(ilc[jt],hm) )//edges
						if( !m->GetMarker(ilc[jt],mrk))
						{
							aret.push_back(ilc[jt]);
							m->SetMarker(ilc[jt],mrk);
						}
				}
				for(ElementArray<Edge>::size_type it = 0; it < aret.size(); it++)
					m->RemMarker(aret.at(it),mrk);
				m->ReleaseMarker(mrk);
			}
		}
		return aret;
	}


	ElementArray<Edge> Cell::getEdges(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		ElementArray<Edge> aret(m);
		if( !m->HideMarker() )
		{
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
			{
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				HandleType last, first;
				HandleType q = lc[0]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				if( invert ^ m->GetMarker(qlc[0],mask) ) aret.push_back(qlc[0]); //node 0
				if( invert ^ m->GetMarker(qlc[1],mask) ) aret.push_back(qlc[1]); //node 1
				first = qlc[0];
				last  = qlc[1];
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
					if( last == ilc[0] ) last = ilc[1];
					else last = ilc[0];
					if( invert ^ m->GetMarker(last,mask) ) 
						aret.push_back(last);
					++it;
				}
			}
			else
			{
				MarkerType mrk = m->CreateMarker();
				adj_type const & lc = m->LowConn(GetHandle());
				for(adj_type::size_type it = 0; it < lc.size(); it++) //faces
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					for(adj_type::size_type jt = 0; jt != ilc.size(); jt++) //edges
						if( (invert ^ m->GetMarker(ilc[jt],mask)) && !m->GetMarker(ilc[jt],mrk))
						{
							aret.push_back(ilc[jt]);
							m->SetMarker(ilc[jt],mrk);
						}
				}
				for(ElementArray<Edge>::size_type it = 0; it != aret.size(); it++)
					m->RemMarker(aret.at(it),mrk);
				m->ReleaseMarker(mrk);
			}
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
			{
				integer i = -1, k = -1, k1 = -1, k2;
				HandleType last, first;
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				i = m->getNext(lc.data(),static_cast<integer>(lc.size()),i,hm);
				HandleType q = lc[i]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				k = m->getNext(qlc.data(),static_cast<integer>(qlc.size()),k,hm);
				if( invert ^ m->GetMarker(qlc[k],mask) ) aret.push_back(qlc[k]); //node 0
				first = qlc[k];
				k = m->getNext(qlc.data(),static_cast<integer>(qlc.size()),k,hm);
				if( invert ^ m->GetMarker(qlc[k],mask) ) aret.push_back(qlc[k]); //node 1
				last = qlc[k];
				i = m->getNext(lc.data(),static_cast<integer>(lc.size()),i,hm);
				HandleType r = lc[i]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				k1 = m->getNext(rlc.data(),static_cast<integer>(rlc.size()),k1,hm);
				k2 = m->getNext(rlc.data(),static_cast<integer>(rlc.size()),k1,hm);
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
				while(it != iend) if( !m->GetMarker(lc[it],hm) ) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					k1 = -1; 
					k1 = m->getNext(ilc.data(),static_cast<integer>(ilc.size()),k1,hm);
					k2 = m->getNext(ilc.data(),static_cast<integer>(ilc.size()),k1,hm);
					if( last == ilc[k1] ) 
						last = ilc[k2];
					else last = ilc[k1];
					if( invert ^ m->GetMarker(last,mask) ) aret.push_back(last);
					++it;
				}
			}
			else
			{
				MarkerType mrk = m->CreateMarker();
				adj_type const & lc = m->LowConn(GetHandle());
				for(adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) ) //faces
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !m->GetMarker(ilc[jt],hm) )//edges
						if( (invert ^ m->GetMarker(ilc[jt],mask)) && !m->GetMarker(ilc[jt],mrk))
						{
							aret.push_back(ilc[jt]);
							m->SetMarker(ilc[jt],mrk);
						}
				}
				for(ElementArray<Edge>::size_type it = 0; it < aret.size(); it++) 
					m->RemMarker(aret.at(it),mrk);
				m->ReleaseMarker(mrk);
			}
		}
		return aret;
	}

	ElementArray<Face> Cell::getFaces() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & lc = m->LowConn(GetHandle());
			return ElementArray<Face>(m,lc.data(),lc.data()+lc.size());
		}
		else
		{
			MarkerType hm = m->HideMarker();
			ElementArray<Face> aret(m);
			adj_type & lc = m->LowConn(GetHandle());
			for(adj_type::size_type it = 0; it < lc.size(); ++it)
				if( !m->GetMarker(lc[it],hm) ) aret.push_back(lc[it]);
			return aret;
		}
	}


	ElementArray<Face> Cell::getFaces(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		ElementArray<Face> aret(m);
		if( !m->HideMarker() )
		{
			adj_type const & lc = m->LowConn(GetHandle());
			for(adj_type::size_type it = 0; it < lc.size(); ++it)
				if( invert ^ m->GetMarker(lc[it],mask) ) aret.push_back(lc[it]);
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
}
#endif
