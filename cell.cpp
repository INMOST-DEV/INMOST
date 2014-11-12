#include "inmost.h"
#if defined(USE_MESH)

namespace INMOST
{
	
	Cell::Cell(Mesh * m)
	:Element(m,CELL)
	{
	}

	Cell::Cell(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Cell & other)
	:Element(m,lid,other)
	{
	}
	
	Cell::Cell(const Cell & other)
	:Element(other)
	{
	}
	
	Cell & Cell::operator =(Cell const & other)
	{
        Element::operator =(other);
        return *this;
	}
	
	Cell::~Cell()
	{
	}
	
	bool Cell::CheckEdgeOrder()
	{
		return true;
	}
	
	bool Cell::FixEdgeOrder()
	{
		return false;
	}
	
	adjacent<Node> Cell::getNodes()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type & hc = HighConn();
			return adjacent<Node>(hc.data(),hc.data()+hc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Node> aret;
			adj_type & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( !hc[it]->GetMarker(hm) ) aret.push_back(hc[it]);
			return aret;
		}
	}


	adjacent<Node> Cell::getNodes(MarkerType mask, bool invert)
	{
		adjacent<Node> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type & hc = HighConn();
			for(adj_type::enumerator it = 0; it < hc.size(); ++it)
				if( invert ^ hc[it]->GetMarker(mask) ) aret.push_back(hc[it]);
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

	adjacent<Edge> Cell::getEdges()
	{
		adjacent<Edge> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			if( GetElementDimension() == 2 ) // This cell is 2d face
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
					Edge * temp = aret.data()[0];
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
			else
			{
				Mesh * m = GetMeshLink();
				MarkerType mrk = m->CreateMarker();
				adj_type & lc = LowConn();
				for(adj_type::enumerator it = 0; it < lc.size(); it++) //faces
				{
					adj_type & ilc = lc[it]->LowConn();
					for(adj_type::enumerator jt = 0; jt < ilc.size(); jt++) //edges
						if( !ilc[jt]->GetMarker(mrk))
						{
							aret.push_back(ilc[jt]);
							ilc[jt]->SetMarker(mrk);
						}
				}
				for(adjacent<Edge>::enumerator it = 0; it < aret.size(); it++)
					aret[it].RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 2 ) // This cell is 2d face
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
					Edge * temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::enumerator it = 1, iend = lc.size()-1;
				while(it < iend) if( !lc[it]->GetMarker(hm) ) //loop over edges
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
			else
			{
				Mesh * m = GetMeshLink();
				MarkerType mrk = m->CreateMarker();
				adj_type & lc = LowConn();
				for(adj_type::enumerator it = 0; it < lc.size(); it++) if( !lc[it]->GetMarker(hm) ) //faces
				{
					adj_type & ilc = lc[it]->LowConn();
					for(adj_type::enumerator jt = 0; jt < ilc.size(); jt++) if( !ilc[jt]->GetMarker(hm) )//edges
						if( !ilc[jt]->GetMarker(mrk))
						{
							aret.push_back(ilc[jt]);
							ilc[jt]->SetMarker(mrk);
						}
				}
				for(adjacent<Edge>::enumerator it = 0; it < aret.size(); it++)
					aret[it].RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		return aret;
	}


	adjacent<Edge> Cell::getEdges(MarkerType mask, bool invert)
	{
		adjacent<Edge> aret;
		if( !GetMeshLink()->HideMarker() )
		{
			if( GetElementDimension() == 2 ) // This cell is 2d face
			{
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				Element * last, * first;
				Element * q = lc[0]; //edge 0
				adj_type & qlc = q->LowConn();
				if( invert ^ qlc[0]->GetMarker(mask) ) aret.push_back(qlc[0]); //node 0
				if( invert ^ qlc[1]->GetMarker(mask) ) aret.push_back(qlc[1]); //node 1
				first = qlc[0];
				last  = qlc[1];
				Element * r = lc[1]; //edge 1
				adj_type & rlc = r->LowConn();
				if( first == rlc[0] || first == rlc[1] )
				{
					last = first;
					if( aret.size() > 1 )
					{
						Edge * temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_type::enumerator it = 1, iend = lc.size()-1;
				while(it < iend) //loop over edges
				{
					adj_type & ilc = lc[it]->LowConn();
					if( last == ilc[0] ) last = ilc[1];
					else last = ilc[0];
					if( invert ^ last->GetMarker(mask) ) aret.push_back(last);
					++it;
				}
			}
			else
			{
				Mesh * m = GetMeshLink();
				MarkerType mrk = m->CreateMarker();
				adj_type & lc = LowConn();
				for(adj_type::enumerator it = 0; it < lc.size(); it++) //faces
				{
					adj_type & ilc = lc[it]->LowConn();
					for(adj_type::enumerator jt = 0; jt != ilc.size(); jt++) //edges
						if( (invert ^ ilc[jt]->GetMarker(mask)) && !ilc[jt]->GetMarker(mrk))
						{
							aret.push_back(ilc[jt]);
							ilc[jt]->SetMarker(mrk);
						}
				}
				for(adjacent<Edge>::enumerator it = 0; it != aret.size(); it++)
					aret[it].RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			if( GetElementDimension() == 2 ) // This cell is 2d face
			{
				INMOST_DATA_ENUM_TYPE i = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k = static_cast<INMOST_DATA_ENUM_TYPE>(-1), 
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1), k2;
				Element * last, * first;
				adj_type & lc = LowConn();
				aret.reserve(lc.size());
				i = Mesh::getNext(lc.data(),lc.size(),i,hm);
				Element * q = lc[i]; //edge 0
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
						Edge * temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_type::enumerator it = 1, iend = lc.size()-1;
				while(it != iend) if( !lc[it]->GetMarker(hm) ) //loop over edges
				{
					adj_type & ilc = lc[it]->LowConn();
					k1 = static_cast<INMOST_DATA_ENUM_TYPE>(-1); 
					k1 = Mesh::getNext(ilc.data(),ilc.size(),k1,hm);
					k2 = Mesh::getNext(ilc.data(),ilc.size(),k1,hm);
					if( last == ilc[k1] ) last = ilc[k2];
					else last = ilc[k1];
					if( invert ^ last->GetMarker(mask) ) aret.push_back(last);
					++it;
				}
			}
			else
			{
				Mesh * m = GetMeshLink();
				MarkerType mrk = m->CreateMarker();
				adj_type & lc = LowConn();
				for(adj_type::enumerator it = 0; it < lc.size(); it++) if( !lc[it]->GetMarker(hm) ) //faces
				{
					adj_type & ilc = lc[it]->LowConn();
					for(adj_type::enumerator jt = 0; jt < ilc.size(); jt++) if( !ilc[jt]->GetMarker(hm) )//edges
						if( (invert ^ ilc[jt]->GetMarker(mask)) && !ilc[jt]->GetMarker(mrk))
						{
							aret.push_back(ilc[jt]);
							ilc[jt]->SetMarker(mrk);
						}
				}
				for(adjacent<Edge>::enumerator it = 0; it < aret.size(); it++) aret[it].RemMarker(mrk);
				m->ReleaseMarker(mrk);
			}
		}
		return aret;
	}

	adjacent<Face> Cell::getFaces()
	{
		if( !GetMeshLink()->HideMarker() )
		{
			adj_type & lc = LowConn();
			return adjacent<Face>(lc.data(),lc.data()+lc.size());
		}
		else
		{
			MarkerType hm = GetMeshLink()->HideMarker();
			adjacent<Face> aret;
			adj_type & lc = LowConn();
			for(adj_type::enumerator it = 0; it < lc.size(); ++it)
				if( !lc[it]->GetMarker(hm) ) aret.push_back(lc[it]);
			return aret;
		}
	}


	adjacent<Face> Cell::getFaces(MarkerType mask, bool invert)
	{
		adjacent<Face> aret;
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
	
	
	
	
	Storage::real Cell::Volume() {Storage::real ret; m_link->GetGeometricData(this,MEASURE,&ret); return ret;}	
}
#endif
