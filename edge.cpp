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
	
	adjacent<Node> Edge::getNodes()
	{
		if( !GetMeshLink()->HideMarker() )
			return adjacent<Node>(low_conn.begin(),low_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Node> aret;
			for(adj_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}
	
	adjacent<Face> Edge::getFaces()
	{
		if( !GetMeshLink()->HideMarker() )
			return adjacent<Face>(high_conn.begin(),high_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Face> aret;
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}
	adjacent<Cell> Edge::getCells()
	{
		adjacent<Cell> aret;
		Mesh * m = GetMeshLink();
		MIDType mrk = m->CreateMarker();
		if( !GetMeshLink()->HideMarker() )
		{
			for(Element::adj_iterator it = high_conn.begin(); it != high_conn.end(); it++) //faces
				for(Element::adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) //cels
					if( !(*jt)->GetMarker(mrk) )
					{
						aret.push_back(*jt);
						(*jt)->SetMarker(mrk);
					}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(Element::adj_iterator it = high_conn.begin(); it != high_conn.end(); it++) if( !(*it)->GetMarker(hm) ) //faces
				for(Element::adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) if( !(*jt)->GetMarker(hm) ) //cels
					if( !(*jt)->GetMarker(mrk) )
					{
						aret.push_back(*jt);
						(*jt)->SetMarker(mrk);
					}
		}
		for(adjacent<Cell>::iterator it = aret.begin(); it != aret.end(); it++)
			it->RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return aret;
	}
	Storage::real Edge::Length() {Storage::real ret; m_link->GetGeometricData(this,MEASURE,&ret); return ret;}
}

#endif
