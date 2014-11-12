#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{

	Node::Node(Mesh * m) :Element(m,NODE)
	{
	}
	
	Node::Node(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Node & other) :Element(m,lid,other){}

	Node::Node(const Node & other) :Element(other){}

	
	Node & Node::operator =(Node const & other)
	{
		Element::operator =(other);
		return *this;
	}
	
	Node::~Node()
	{
		
	}
	
	adjacent<Edge> Node::getEdges()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			return adjacent<Edge>(hc.data(),hc.data()+hc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Edge> aret;
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( !hc[it]->GetMarker(hm) ) aret.push_back(hc[it]);
			return aret;
		}
	}

	adjacent<Edge> Node::getEdges(MarkerType mask, bool invert)
	{
		adjacent<Edge> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( invert ^ hc[it]->GetMarker(mask) ) aret.push_back(hc[it]);
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
	
	adjacent<Face> Node::getFaces()
	{
		adjacent<Face> aret;
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreateMarker();
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); it++) //edges
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) //faces
					if( !ihc[jt]->GetMarker(mrk))
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
			for(adj_type::enumerator it = 0; it < hc.size(); it++) if( !hc[it]->GetMarker(hm) )//edges
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) if( !ihc[jt]->GetMarker(hm) ) //faces
					if( !ihc[jt]->GetMarker(mrk))
					{
						aret.push_back(ihc[jt]);
						ihc[jt]->SetMarker(mrk);
					}
			}
		}
		for(adjacent<Face>::enumerator it = 0; it < aret.size(); it++) aret[it].RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return aret;
	}

	adjacent<Face> Node::getFaces(MarkerType mask, bool invert)
	{
		adjacent<Face> aret;
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreateMarker();
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); it++) //edges
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) //faces
					if( (invert ^ ihc[jt]->GetMarker(mask)) && !ihc[jt]->GetMarker(mrk))
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
			for(adj_type::enumerator it = 0; it < hc.size(); it++) if( !hc[it]->GetMarker(hm) )//edges
			{
				adj_type const & ihc = hc[it]->HighConn();
				for(adj_type::enumerator jt = 0; jt < ihc.size(); jt++) if( !ihc[jt]->GetMarker(hm) ) //faces
					if( (invert ^ ihc[jt]->GetMarker(mask)) && !ihc[jt]->GetMarker(mrk))
					{
						aret.push_back(ihc[jt]);
						ihc[jt]->SetMarker(mrk);
					}
			}
		}
		for(adjacent<Face>::enumerator it = 0; it < aret.size(); it++) aret[it].RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return aret;
	}
	
	adjacent<Cell> Node::getCells()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type const & lc = LowConn();
			return adjacent<Cell>(lc.data(),lc.data()+lc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Cell> aret;
			adj_type const & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( !lc[it]->GetMarker(hm) ) aret.push_back(lc[it]);
			return aret;
		}
	}

	adjacent<Cell> Node::getCells(MarkerType mask, bool invert)
	{
		adjacent<Cell> aret;
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

	Storage::real_array Node::Coords() {return RealArrayDF(GetMeshLink()->CoordsTag());}
}

#endif
