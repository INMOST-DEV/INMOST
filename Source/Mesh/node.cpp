#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{

	
	ElementArray<Edge> Node::getEdges() const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			return ElementArray<Edge>(m,hc.data(),hc.data()+hc.size());
		}
		else
		{
			MarkerType hm = m->HideMarker();
			ElementArray<Edge> aret(m);
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( !m->GetMarker(hc[it],hm) ) aret.push_back(hc[it]);
			return aret;
		}
	}

	ElementArray<Edge> Node::getEdges(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Edge> aret(m);
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( invert ^ m->GetMarker(hc[it],mask) ) aret.push_back(hc[it]);
		}
		else
		{
			MarkerType hm = m->HideMarker();
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
				if( (invert ^ m->GetMarker(hc[it],mask)) && !m->GetMarker(hc[it],hm) ) aret.push_back(hc[it]);
		}
		return aret;
	}
	
	ElementArray<Face> Node::getFaces() const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Face> aret(m);
		MarkerType mrk = m->CreatePrivateMarker();
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); it++) //edges
			{
				adj_type const & ihc = m->HighConn(hc[it]);
				for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) //faces
					if( !m->GetPrivateMarker(ihc[jt],mrk))
					{
						aret.push_back(ihc[jt]);
						m->SetPrivateMarker(ihc[jt],mrk);
					}
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); it++) if( !m->GetMarker(hc[it],hm) )//edges
			{
				adj_type const & ihc = m->HighConn(hc[it]);
				for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !m->GetMarker(ihc[jt],hm) ) //faces
					if( !m->GetPrivateMarker(ihc[jt],mrk))
					{
						aret.push_back(ihc[jt]);
						m->SetPrivateMarker(ihc[jt],mrk);
					}
			}
		}
		for(ElementArray<Face>::size_type it = 0; it < aret.size(); it++) m->RemPrivateMarker(aret.at(it),mrk);
		m->ReleasePrivateMarker(mrk);
		return aret;
	}

	ElementArray<Face> Node::getFaces(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Face> aret(m);
		MarkerType mrk = m->CreatePrivateMarker();
		if( !m->HideMarker() )
		{
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); it++) //edges
			{
				adj_type const & ihc = m->HighConn(hc[it]);
				for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) //faces
					if( (invert ^ m->GetMarker(ihc[jt],mask)) && !m->GetPrivateMarker(ihc[jt],mrk))
					{
						aret.push_back(ihc[jt]);
						m->SetPrivateMarker(ihc[jt],mrk);
					}
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			adj_type const & hc = m->HighConn(GetHandle());
			for(adj_type::size_type it = 0; it < hc.size(); it++) if( !m->GetMarker(hc[it],hm) )//edges
			{
				adj_type const & ihc = m->HighConn(hc[it]);
				for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !m->GetMarker(ihc[jt],hm) ) //faces
					if( (invert ^ m->GetMarker(ihc[jt],mask)) && !m->GetPrivateMarker(ihc[jt],mrk))
					{
						aret.push_back(ihc[jt]);
						m->SetPrivateMarker(ihc[jt],mrk);
					}
			}
		}
		for(ElementArray<Face>::size_type it = 0; it < aret.size(); it++) m->RemMarker(aret.at(it),mrk);
		m->ReleasePrivateMarker(mrk);
		return aret;
	}
	
	ElementArray<Cell> Node::getCells() const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		if( !m->HideMarker() )
		{
			adj_type const & lc = m->LowConn(GetHandle());
			return ElementArray<Cell>(m,lc.data(),lc.data()+lc.size());
		}
		else
		{
			MarkerType hm = m->HideMarker();
			ElementArray<Cell> aret(m);
			adj_type const & lc = m->LowConn(GetHandle());
			for(adj_type::size_type it = 0; it < lc.size(); ++it)
				if( !m->GetMarker(lc[it],hm) ) aret.push_back(lc[it]);
			return aret;
		}
	}

	ElementArray<Cell> Node::getCells(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Cell> aret(m);
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
				if( (invert ^ m->GetMarker(lc[it],mask)) && !m->GetMarker(lc[it],hm) ) aret.push_back(lc[it]);
		}
		return aret;
	}

	Storage::real_array Node::Coords() const {return GetMeshLink()->RealArrayDF(GetHandle(),GetMeshLink()->CoordsTag());}
}

#endif
