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
			aret.reserve(32);
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
		aret.reserve(32);
		adj_type const & hc = m->HighConn(GetHandle());
		if( !m->HideMarker() )
		{
			if( isPrivate(mask) )
			{
				for(adj_type::size_type it = 0; it < hc.size(); ++it)
					if( invert ^ m->GetPrivateMarker(hc[it],mask) ) aret.push_back(hc[it]);
			}
			else
			{
				for(adj_type::size_type it = 0; it < hc.size(); ++it)
					if( invert ^ m->GetMarker(hc[it],mask) ) aret.push_back(hc[it]);
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			if( isPrivate(mask) )
			{
				for(adj_type::size_type it = 0; it < hc.size(); ++it)
					if( (invert ^ m->GetPrivateMarker(hc[it],mask)) && !m->GetMarker(hc[it],hm) ) aret.push_back(hc[it]);
			}
			else
			{
				for(adj_type::size_type it = 0; it < hc.size(); ++it)
					if( (invert ^ m->GetMarker(hc[it],mask)) && !m->GetMarker(hc[it],hm) ) aret.push_back(hc[it]);
			}
		}
		return aret;
	}
	
	ElementArray<Face> Node::getFaces() const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Face> aret(m);
		aret.reserve(32);
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
		aret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return aret;
	}

	ElementArray<Face> Node::getFaces(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Face> aret(m);
		aret.reserve(32);
		MarkerType mrk = m->CreatePrivateMarker();
		if( isPrivate(mask) )
		{
			if( !m->HideMarker() )
			{
				adj_type const & hc = m->HighConn(GetHandle());
				for(adj_type::size_type it = 0; it < hc.size(); it++) //edges
				{
					adj_type const & ihc = m->HighConn(hc[it]);
					for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) //faces
						if( (invert ^ m->GetPrivateMarker(ihc[jt],mask)) && !m->GetPrivateMarker(ihc[jt],mrk))
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
						if( (invert ^ m->GetPrivateMarker(ihc[jt],mask)) && !m->GetPrivateMarker(ihc[jt],mrk))
						{
							aret.push_back(ihc[jt]);
							m->SetPrivateMarker(ihc[jt],mrk);
						}
				}
			}
		}
		else
		{
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
		}
		aret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return aret;
	}
	
	ElementArray<Cell> Node::getCells() const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Cell> aret(m);
		aret.reserve(32);
		MarkerType mrk = m->CreatePrivateMarker();
		if( !m->HideMarker() )
		{
			adj_type const & edges = m->HighConn(GetHandle()); // all edges
			for(adj_type::const_iterator eit = edges.begin(); eit != edges.end(); ++eit)
			{
				adj_type const & faces = m->HighConn(*eit); // all faces of the edge
				for(adj_type::const_iterator fit = faces.begin(); fit != faces.end(); ++fit)
				{
					adj_type const & cells = m->HighConn(*fit);
					for(adj_type::const_iterator cit = cells.begin(); cit != cells.end(); ++cit)
					{
						if( !m->GetPrivateMarker(*cit,mrk) )
						{
							aret.push_back(*cit);
							m->SetPrivateMarker(*cit,mrk);
						}
					}
				}
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			adj_type const & edges = m->HighConn(GetHandle()); // all edges
			for(adj_type::const_iterator eit = edges.begin(); eit != edges.end(); ++eit) if( !m->GetMarker(*eit,hm) )
			{
				adj_type const & faces = m->HighConn(*eit); // all faces of the edge
				for(adj_type::const_iterator fit = faces.begin(); fit != faces.end(); ++fit) if( !m->GetMarker(*fit,hm) )
				{
					adj_type const & cells = m->HighConn(*fit);
					for(adj_type::const_iterator cit = cells.begin(); cit != cells.end(); ++cit) if( !m->GetMarker(*cit,hm) )
					{
						if( !m->GetPrivateMarker(*cit,mrk) )
						{
							aret.push_back(*cit);
							m->SetPrivateMarker(*cit,mrk);
						}
					}
				}
			}
		}
		aret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return aret;
	}

	ElementArray<Cell> Node::getCells(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==NODE);
		Mesh * m = GetMeshLink();
		ElementArray<Cell> aret(m);
		aret.reserve(32);
		MarkerType mrk = m->CreatePrivateMarker();
		if( !m->HideMarker() )
		{
			if( isPrivate(mask) )
			{
				adj_type const & edges = m->HighConn(GetHandle()); // all edges
				for(adj_type::const_iterator eit = edges.begin(); eit != edges.end(); ++eit)
				{
					adj_type const & faces = m->HighConn(*eit); // all faces of the edge
					for(adj_type::const_iterator fit = faces.begin(); fit != faces.end(); ++fit)
					{
						adj_type const & cells = m->HighConn(*fit);
						for(adj_type::const_iterator cit = cells.begin(); cit != cells.end(); ++cit)
						{
							if( (invert ^ m->GetPrivateMarker(*cit,mask)) && !m->GetPrivateMarker(*cit,mrk) )
							{
								aret.push_back(*cit);
								m->SetPrivateMarker(*cit,mrk);
							}
						}
					}
				}
			}
			else
			{
				adj_type const & edges = m->HighConn(GetHandle()); // all edges
				for(adj_type::const_iterator eit = edges.begin(); eit != edges.end(); ++eit)
				{
					adj_type const & faces = m->HighConn(*eit); // all faces of the edge
					for(adj_type::const_iterator fit = faces.begin(); fit != faces.end(); ++fit)
					{
						adj_type const & cells = m->HighConn(*fit);
						for(adj_type::const_iterator cit = cells.begin(); cit != cells.end(); ++cit)
						{
							if( (invert ^ m->GetMarker(*cit,mask)) && !m->GetPrivateMarker(*cit,mrk) )
							{
								aret.push_back(*cit);
								m->SetPrivateMarker(*cit,mrk);
							}
						}
					}
				}
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			if( isPrivate(mask) )
			{
				adj_type const & edges = m->HighConn(GetHandle()); // all edges
				for(adj_type::const_iterator eit = edges.begin(); eit != edges.end(); ++eit) if( !m->GetMarker(*eit,hm) )
				{
					adj_type const & faces = m->HighConn(*eit); // all faces of the edge
					for(adj_type::const_iterator fit = faces.begin(); fit != faces.end(); ++fit) if( !m->GetMarker(*fit,hm) )
					{
						adj_type const & cells = m->HighConn(*fit);
						for(adj_type::const_iterator cit = cells.begin(); cit != cells.end(); ++cit) if( !m->GetMarker(*cit,hm) )
						{
							if( (invert ^ m->GetPrivateMarker(*cit,mask)) && !m->GetPrivateMarker(*cit,mrk) )
							{
								aret.push_back(*cit);
								m->SetPrivateMarker(*cit,mrk);
							}
						}
					}
				}
			}
			else
			{
				adj_type const & edges = m->HighConn(GetHandle()); // all edges
				for(adj_type::const_iterator eit = edges.begin(); eit != edges.end(); ++eit) if( !m->GetMarker(*eit,hm) )
				{
					adj_type const & faces = m->HighConn(*eit); // all faces of the edge
					for(adj_type::const_iterator fit = faces.begin(); fit != faces.end(); ++fit) if( !m->GetMarker(*fit,hm) )
					{
						adj_type const & cells = m->HighConn(*fit);
						for(adj_type::const_iterator cit = cells.begin(); cit != cells.end(); ++cit) if( !m->GetMarker(*cit,hm) )
						{
							if( (invert ^ m->GetMarker(*cit,mask)) && !m->GetPrivateMarker(*cit,mrk) )
							{
								aret.push_back(*cit);
								m->SetPrivateMarker(*cit,mrk);
							}
						}
					}
				}
			}
		}
		aret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return aret;
	}

	Storage::real_array Node::Coords() const {return GetMeshLink()->RealArrayDF(GetHandle(),GetMeshLink()->CoordsTag());}
}

#endif
