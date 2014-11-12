#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{
	
	Edge::Edge(Mesh * m)
	:Element(m,EDGE)
	{
	}

	Edge::Edge(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Edge & other)
	:Element(m,lid,other)
	{
	}
	
	Edge::Edge(const Edge & other)
	:Element(other)
	{
	}
	
	Edge & Edge::operator =(Edge const & other)
	{
		Element::operator =(other);
		return *this;
	}
	
	Edge::~Edge()
	{
	}

	Node * Edge::getBeg() const 
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & lc = LowConn();
			if( lc.empty() )
				return NULL;
			return lc.front()->getAsNode();
		}
		else
		{
			adj_type const & lc = LowConn();
			if( !lc.empty() )
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1);
				MarkerType hm = GetMeshLink()->HideMarker();
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				if( i != lc.size() ) return lc[i]->getAsNode();
			}
			return NULL;
		}
	}
	Node * Edge::getEnd() const 
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & lc = LowConn();
			if( lc.size() < 2 )
				return NULL;
			return lc.back()->getAsNode();
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
				if( i != lc.size() ) return lc[i]->getAsNode();
			}
			return NULL;
		}
	}
	
	adjacent<Node> Edge::getNodes()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type & lc = LowConn();
			return adjacent<Node>(lc.data(),lc.data()+lc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Node> aret;
			adj_type const & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( !lc[it]->GetMarker(hm) ) aret.push_back(lc[it]);
			return aret;
		}
	}

	adjacent<Node> Edge::getNodes(MarkerType mask, bool invert)
	{
		adjacent<Node> aret;
		if( !GetMeshLink()->HideMarker() )
		{	
			adj_type const & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( invert ^ lc[it]->GetMarker(mask) ) aret.push_back(lc[it]);
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adj_type const & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( (invert ^ lc[it]->GetMarker(mask)) && !lc[it]->GetMarker(hm) ) aret.push_back(lc[it]);
		}
		return aret;
	}
	
	adjacent<Face> Edge::getFaces()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			return adjacent<Face>(hc.data(),hc.data()+hc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Face> aret;
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( !hc[it]->GetMarker(hm) ) aret.push_back(hc[it]);
			return aret;
		}
	}


	adjacent<Face> Edge::getFaces(MarkerType mask, bool invert)
	{
		adjacent<Face> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( (invert ^ hc[it]->GetMarker(mask)) ) aret.push_back(hc[it]);
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( (invert ^ hc[it]->GetMarker(mask)) && !hc[it]->GetMarker(hm) ) aret.push_back(hc[it]);
		}
		return aret;
	}

	adjacent<Cell> Edge::getCells()
	{
		adjacent<Cell> aret;
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreateMarker();
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); it++) //faces
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) //cels
					if( !ihc[jt]->GetMarker(mrk) )
					{
						aret.push_back(ihc[jt]);
						ihc[jt]->SetMarker(mrk);
					}
			}
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); it++) if( !hc[it]->GetMarker(hm) ) //faces
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) if( !ihc[jt]->GetMarker(hm) ) //cels
					if( !ihc[jt]->GetMarker(mrk) )
					{
						aret.push_back(ihc[jt]);
						ihc[jt]->SetMarker(mrk);
					}
			}
		}
		for(adjacent<Cell>::enumerator it = 0; it < aret.size(); it++) aret[it].RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return aret;
	}

	adjacent<Cell> Edge::getCells(MarkerType mask, bool invert)
	{
		adjacent<Cell> aret;
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreateMarker();
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); it++) //faces
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) //cels
					if( (invert ^ ihc[jt]->GetMarker(mask)) && !ihc[jt]->GetMarker(mrk) )
					{
						aret.push_back(ihc[jt]);
						ihc[jt]->SetMarker(mrk);
					}
			}
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); it++) if( !hc[it]->GetMarker(hm) ) //faces
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) if( !ihc[jt]->GetMarker(hm) ) //cels
					if( (invert ^ ihc[jt]->GetMarker(mask)) && !ihc[jt]->GetMarker(mrk) )
					{
						aret.push_back(ihc[jt]);
						ihc[jt]->SetMarker(mrk);
					}
			}
		}
		for(adjacent<Cell>::enumerator it = 0; it < aret.size(); it++) aret[it].RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return aret;
	}

	Storage::real Edge::Length() {Storage::real ret; m_link->GetGeometricData(this,MEASURE,&ret); return ret;}
}

#endif
