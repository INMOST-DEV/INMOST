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
			return adjacent<Edge>(high_conn.begin(),high_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Edge> aret;
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}

	adjacent<Edge> Node::getEdges(MIDType mask, bool invert)
	{
		adjacent<Edge> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( invert ^ (*it)->GetMarker(mask) ) aret.push_back((*it));
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(adj_iterator it = high_conn.begin(); it != high_conn.end(); ++it)
				if( (invert ^ (*it)->GetMarker(mask)) && !(*it)->GetMarker(hm) ) aret.push_back((*it));
		}
		return aret;
	}
	
	adjacent<Face> Node::getFaces()
	{
		adjacent<Face> aret;
		Mesh * m = GetMeshLink();
		MIDType mrk = m->CreateMarker();
		if( !GetMeshLink()->HideMarker() )
		{
			for(Element::adj_iterator it = high_conn.begin(); it != high_conn.end(); it++) //edges
				for(Element::adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) //faces
					if( !(*jt)->GetMarker(mrk))
					{
						aret.push_back(*jt);
						(*jt)->SetMarker(mrk);
					}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(Element::adj_iterator it = high_conn.begin(); it != high_conn.end(); it++) if( !(*it)->GetMarker(hm) )//edges
				for(Element::adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) if( !(*jt)->GetMarker(hm) ) //faces
					if( !(*jt)->GetMarker(mrk))
					{
						aret.push_back(*jt);
						(*jt)->SetMarker(mrk);
					}
		}
		for(adjacent<Face>::iterator it = aret.begin(); it != aret.end(); it++)
			it->RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return aret;
	}

	adjacent<Face> Node::getFaces(MIDType mask, bool invert)
	{
		adjacent<Face> aret;
		Mesh * m = GetMeshLink();
		MIDType mrk = m->CreateMarker();
		if( !GetMeshLink()->HideMarker() )
		{
			for(Element::adj_iterator it = high_conn.begin(); it != high_conn.end(); it++) //edges
				for(Element::adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) //faces
					if( (invert ^ (*jt)->GetMarker(mask)) && !(*jt)->GetMarker(mrk))
					{
						aret.push_back(*jt);
						(*jt)->SetMarker(mrk);
					}
		}
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			for(Element::adj_iterator it = high_conn.begin(); it != high_conn.end(); it++) if( !(*it)->GetMarker(hm) )//edges
				for(Element::adj_iterator jt = (*it)->high_conn.begin(); jt != (*it)->high_conn.end(); jt++) if( !(*jt)->GetMarker(hm) ) //faces
					if( (invert ^ (*jt)->GetMarker(mask)) && !(*jt)->GetMarker(mrk))
					{
						aret.push_back(*jt);
						(*jt)->SetMarker(mrk);
					}
		}
		for(adjacent<Face>::iterator it = aret.begin(); it != aret.end(); it++)
			it->RemMarker(mrk);
		m->ReleaseMarker(mrk);
		return aret;
	}
	
	adjacent<Cell> Node::getCells()
	{
		if( !GetMeshLink()->HideMarker() )
			return adjacent<Cell>(low_conn.begin(),low_conn.end());
		else
		{
			MIDType hm = GetMeshLink()->HideMarker();
			adjacent<Cell> aret;
			for(adj_iterator it = low_conn.begin(); it != low_conn.end(); ++it)
				if( !(*it)->GetMarker(hm) ) aret.push_back((*it));
			return aret;
		}
	}

	adjacent<Cell> Node::getCells(MIDType mask, bool invert)
	{
		adjacent<Cell> aret;
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

	Storage::real_array Node::Coords() {return RealArrayDF(GetMeshLink()->CoordsTag());}
}

#endif
