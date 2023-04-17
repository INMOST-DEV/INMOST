#include "inmost.h"
#include <iomanip>
#if defined(USE_MESH)
namespace INMOST
{

	const char * Element::GeometricTypeName(GeometricType t)
	{
		switch(t)
		{
			case Unset:        return "Unset";
			case Vertex:       return "Vertex";
			case Line:         return "Line";
			case MultiLine:    return "MultiLine";
			case Tri:          return "Tri";
			case Quad:         return "Quad";
			case Polygon:      return "Polygon";
			case MultiPolygon: return "MultiPolygon";
			case Tet:          return "Tet";
			case Hex:          return "Hex";
			case Prism:        return "Prism";
			case Pyramid:      return "Pyramid";
			case Polyhedron:   return "Polyhedron";
			case Set:          return "Set";
		};
		return "Unset";
	}
	
	
	unsigned int Element::GetGeometricDimension(GeometricType m_type)
	{
		switch(m_type)
		{
			case Vertex:		return 0;
			case Line:			return 1;
			case MultiLine:
			case Tri:
			case Quad: 			
			case Polygon:		return 2;
			case MultiPolygon:
			case Tet:
			case Hex:
			case Prism:
			case Pyramid:
			case Polyhedron:	return 3;
			default: return UINT_MAX;
		}
		return UINT_MAX;
	}
	
	const char * Element::StatusName(Status s)
	{
		switch(s)
		{
			case Owned:  return "Owned";
			case Shared: return "Shared";
			case Ghost:  return "Ghost";
		}
		return "Unknown";
	}
	
	
	
	
	
	Storage::enumerator Element::nbAdjElements(ElementType _etype)const 
	{
		if( GetHandleElementType(GetHandle()) == ESET ) return getAsSet()->nbAdjElements(_etype);
		assert( !(_etype & MESH) );
		assert( !(_etype & ESET) );
		Mesh * mesh = GetMeshLink();
		std::vector<HandleType> result;
		INMOST_DATA_INTEGER_TYPE conn[4] = {0,0,0,0};
		INMOST_DATA_INTEGER_TYPE myconn = -1, i = 0;
		enumerator ret = 0;
		
		for(ElementType e = NODE; e <= CELL; e = NextElementType(e))
		{
			if( _etype & e ) conn[i] = 1;
			if( GetElementType() & e ) myconn = i;
			i++;
		}
		if( !mesh->HideMarker() )
		{
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					ret += 1;
				}
				else if( i == myconn + 1 )
				{
					ret += static_cast<enumerator>(mesh->HighConn(GetHandle()).size());
				}
				else if( i == myconn - 1 )
				{
					ret += static_cast<enumerator>(mesh->LowConn(GetHandle()).size());
				}
				else if( i == myconn + 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & hc = mesh->HighConn(GetHandle());
					for(adj_type::size_type it = 0; it < hc.size(); it++)
					{
						adj_type const & ihc = mesh->HighConn(hc[it]);
						for(adj_type::size_type jt = 0; jt < ihc.size(); jt++)
							if( !mesh->GetPrivateMarker(ihc[jt],mrk) )
							{
								result.push_back(ihc[jt]);
								mesh->SetPrivateMarker(ihc[jt],mrk);
							}
					}
					ret += static_cast<enumerator>(result.size());
					if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0],(enumerator)result.size(),mrk);
					result.clear();
					mesh->ReleasePrivateMarker(mrk);
				}					
				else if( i == myconn - 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & lc = mesh->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++)
					{
						adj_type const & ilc = mesh->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt < ilc.size(); jt++)
							if( !mesh->GetPrivateMarker(ilc[jt],mrk) )
							{
								result.push_back(ilc[jt]);
								mesh->SetPrivateMarker(ilc[jt],mrk);
							}
					}		
					ret += static_cast<enumerator>(result.size());
					if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
					result.clear();
					mesh->ReleasePrivateMarker(mrk);
				}
				else if( i == myconn + 3 ) // these are cells of the node
				{
					if( mesh->LowConnTag().isDefined(NODE))
						ret += static_cast<enumerator>(mesh->LowConn(GetHandle()).size());
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& hc = mesh->HighConn(GetHandle());
						for (adj_type::size_type it = 0; it < hc.size(); it++)
						{
							adj_type const& ihc = mesh->HighConn(hc[it]);
							for (adj_type::size_type jt = 0; jt < ihc.size(); jt++)
							{
								adj_type const& jhc = mesh->HighConn(ihc[jt]);
								for (adj_type::size_type kt = 0; kt < jhc.size(); kt++)
									if (!mesh->GetPrivateMarker(jhc[kt], mrk))
									{
										result.push_back(jhc[kt]);
										mesh->SetPrivateMarker(jhc[kt], mrk);
									}
							}
						}
						ret += static_cast<enumerator>(result.size());
						if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
				}					
				else if( i == myconn - 3 ) // these are nodes of the cell
				{
					if( mesh->HighConnTag().isDefined(CELL))
						ret += static_cast<enumerator>(mesh->HighConn(GetHandle()).size());
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& lc = mesh->LowConn(GetHandle());
						for (adj_type::size_type it = 0; it < lc.size(); it++)
						{
							adj_type const& ilc = mesh->LowConn(lc[it]);
							for (adj_type::size_type jt = 0; jt < ilc.size(); jt++)
							{
								adj_type const& jlc = mesh->LowConn(ilc[jt]);
								for (adj_type::size_type kt = 0; kt < jlc.size(); kt++)
									if (!mesh->GetPrivateMarker(jlc[kt], mrk))
									{
										result.push_back(jlc[kt]);
										mesh->SetPrivateMarker(jlc[kt], mrk);
									}
							}
						}
						ret += static_cast<enumerator>(result.size());
						if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
				}
			}
		}
		else
		{
			MarkerType hm = mesh->HideMarker();
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					if( !GetMarker(hm) ) ret ++;
				}
				else if( i == myconn + 1 )
				{
					adj_type const & hc = mesh->HighConn(GetHandle());
					for(adj_type::size_type it = 0; it < hc.size(); ++it)
						if( !mesh->GetMarker(hc[it],hm) ) ret++;
				}
				else if( i == myconn - 1 )
				{
					adj_type const & lc = mesh->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); ++it)
						if( !mesh->GetMarker(lc[it],hm) ) ret++;
				}
				else if( i == myconn + 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & hc = mesh->HighConn(GetHandle());
					for(adj_type::size_type it = 0; it < hc.size(); it++) if( !mesh->GetMarker(hc[it],hm) )
					{
						adj_type const & ihc = mesh->HighConn(hc[it]);
						for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !mesh->GetMarker(ihc[jt],hm) )
						{
							if( !mesh->GetPrivateMarker(ihc[jt],mrk) )
							{
								result.push_back(ihc[jt]);
								mesh->SetPrivateMarker(ihc[jt],mrk);
							}
						}
					}
					ret += static_cast<enumerator>(result.size());
					if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
					result.clear();
					mesh->ReleasePrivateMarker(mrk);
				}
				else if( i == myconn - 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & lc = mesh->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++) if( !mesh->GetMarker(lc[it],hm) )
					{
						adj_type const & ilc = mesh->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !mesh->GetMarker(ilc[jt],hm) )
							if( !mesh->GetPrivateMarker(ilc[jt],mrk) )
							{
								result.push_back(ilc[jt]);
								mesh->SetPrivateMarker(ilc[jt],mrk);
							}
					}
					ret += static_cast<enumerator>(result.size());
					if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
					result.clear();
					mesh->ReleasePrivateMarker(mrk);
				}
				else if( i == myconn + 3 ) // these are cells of the node
				{
					if (mesh->LowConnTag().isDefined(NODE))
					{
						adj_type const& lc = mesh->LowConn(GetHandle());
						for (adj_type::size_type it = 0; it < lc.size(); ++it)
							if (!mesh->GetMarker(lc[it], hm)) ret++;
					}
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& hc = mesh->HighConn(GetHandle());
						for (adj_type::size_type it = 0; it < hc.size(); it++) if (!mesh->GetMarker(hc[it], hm))
						{
							adj_type const& ihc = mesh->HighConn(hc[it]);
							for (adj_type::size_type jt = 0; jt < ihc.size(); jt++) if (!mesh->GetMarker(ihc[jt], hm))
							{
								adj_type const& jhc = mesh->HighConn(ihc[jt]);
								for (adj_type::size_type kt = 0; kt < jhc.size(); kt++) if (!mesh->GetMarker(jhc[kt], hm))
									if (!mesh->GetPrivateMarker(jhc[kt], mrk))
									{
										result.push_back(jhc[kt]);
										mesh->SetPrivateMarker(jhc[kt], mrk);
									}
							}
						}
						ret += static_cast<enumerator>(result.size());
						if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
				}
				else if( i == myconn - 3 ) // these are nodes of the cell
				{
					if (mesh->HighConnTag().isDefined(CELL))
					{
						adj_type const& hc = mesh->HighConn(GetHandle());
						for (adj_type::size_type it = 0; it < hc.size(); ++it)
							if (!mesh->GetMarker(hc[it], hm)) ret++;
					}
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& lc = mesh->LowConn(GetHandle());
						for (adj_type::size_type it = 0; it < lc.size(); it++) if (!mesh->GetMarker(lc[it], hm))
						{
							adj_type const& ilc = mesh->LowConn(lc[it]);
							for (adj_type::size_type jt = 0; jt < ilc.size(); jt++) if (!mesh->GetMarker(ilc[jt], hm))
							{
								adj_type const& jlc = mesh->LowConn(ilc[jt]);
								for (adj_type::size_type kt = 0; kt < jlc.size(); kt++) if (!mesh->GetMarker(jlc[kt], hm))
									if (!mesh->GetPrivateMarker(jlc[kt], mrk))
									{
										result.push_back(jlc[kt]);
										mesh->SetPrivateMarker(jlc[kt], mrk);
									}
							}
						}
						ret += static_cast<enumerator>(result.size());
						if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
				}
			}
		}
		return ret;
	}


	Storage::enumerator Element::nbAdjElements(ElementType _etype, MarkerType mask, bool invert)const 
	{
		if( GetHandleElementType(GetHandle()) == ESET ) return getAsSet()->nbAdjElements(_etype,mask,invert);
		assert( !(_etype & MESH) );
		assert( !(_etype & ESET) );
		Mesh * mesh = GetMeshLink();
		std::vector<HandleType> result;
		INMOST_DATA_INTEGER_TYPE conn[4] = {0,0,0,0};
		INMOST_DATA_INTEGER_TYPE myconn = -1, i = 0;
		enumerator ret = 0;
		
		for(ElementType e = NODE; e <= CELL; e = NextElementType(e))
		{
			if( _etype & e ) conn[i] = 1;
			if( GetElementType() & e ) myconn = i;
			i++;
		}
		if( isPrivate(mask) )
		{
			if( !mesh->HideMarker() )
			{
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( invert ^ GetPrivateMarker(mask) ) ret += 1;
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it) 
							if( invert ^ mesh->GetPrivateMarker(hc[it],mask) ) ret++;
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it) 
							if( invert ^ mesh->GetPrivateMarker(lc[it],mask) ) ret++;
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++)
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								if( (invert ^ mesh->GetPrivateMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++)
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								if( (invert ^ mesh->GetPrivateMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 ) // these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if (invert ^ mesh->GetPrivateMarker(lc[it], mask)) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++)
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++)
									{
										if ((invert ^ mesh->GetPrivateMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 ) // these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if (invert ^ mesh->GetPrivateMarker(hc[it], mask)) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++)
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++)
									{
										if ((invert ^ mesh->GetPrivateMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
			else
			{
				MarkerType hm = mesh->HideMarker();
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( !GetMarker(hm) && (invert ^ GetPrivateMarker(mask)) ) ret ++;
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it)
							if( !mesh->GetMarker(hc[it],hm) && (invert ^ mesh->GetPrivateMarker(hc[it],mask)) ) ret++;
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it)
							if( !mesh->GetMarker(lc[it],hm) && (invert ^ mesh->GetPrivateMarker(lc[it],mask)) ) ret++;
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++) if( !mesh->GetMarker(hc[it],hm) )
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !mesh->GetMarker(ihc[jt],hm) )
								if( (invert ^ mesh->GetPrivateMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++) if( !mesh->GetMarker(lc[it],hm) )
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !mesh->GetMarker(ilc[jt],hm) )
								if( (invert ^ mesh->GetPrivateMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}		
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 ) // these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if (!mesh->GetMarker(lc[it], hm) && (invert ^ mesh->GetPrivateMarker(lc[it], mask))) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++) if (!mesh->GetMarker(hc[it], hm))
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++) if (!mesh->GetMarker(ihc[jt], hm))
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++) if (!mesh->GetMarker(jhc[kt], hm))
									{
										if ((invert ^ mesh->GetPrivateMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 ) // these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if (!mesh->GetMarker(hc[it], hm) && (invert ^ mesh->GetPrivateMarker(hc[it], mask))) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++) if (!mesh->GetMarker(lc[it], hm))
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++) if (!mesh->GetMarker(ilc[jt], hm))
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++) if (!mesh->GetMarker(jlc[kt], hm))
									{
										if ((invert ^ mesh->GetPrivateMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
		}
		else
		{
			if( !mesh->HideMarker() )
			{
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( invert ^ GetMarker(mask) ) ret += 1;
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it) 
							if( invert ^ mesh->GetMarker(hc[it],mask) ) ret++;
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it) 
							if( invert ^ mesh->GetMarker(lc[it],mask) ) ret++;
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++)
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								if( (invert ^ mesh->GetMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++)
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								if( (invert ^ mesh->GetMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 ) // these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if (invert ^ mesh->GetMarker(lc[it], mask)) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++)
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++)
									{
										if ((invert ^ mesh->GetMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 )// these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if (invert ^ mesh->GetMarker(hc[it], mask)) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++)
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++)
									{
										if ((invert ^ mesh->GetMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
			else
			{
				MarkerType hm = mesh->HideMarker();
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( !GetMarker(hm) && (invert ^ GetMarker(mask)) ) ret ++;
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it)
							if( !mesh->GetMarker(hc[it],hm) && (invert ^ mesh->GetMarker(hc[it],mask)) ) ret++;
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it)
							if( !mesh->GetMarker(lc[it],hm) && (invert ^ mesh->GetMarker(lc[it],mask)) ) ret++;
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++) if( !mesh->GetMarker(hc[it],hm) )
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !mesh->GetMarker(ihc[jt],hm) )
								if( (invert ^ mesh->GetMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++) if( !mesh->GetMarker(lc[it],hm) )
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !mesh->GetMarker(ilc[jt],hm) )
								if( (invert ^ mesh->GetMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}		
						ret += static_cast<enumerator>(result.size());
						if( !result.empty() ) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(),mrk);
						result.clear();
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 )// these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if (!mesh->GetMarker(lc[it], hm) && (invert ^ mesh->GetMarker(lc[it], mask))) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++) if (!mesh->GetMarker(hc[it], hm))
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++) if (!mesh->GetMarker(ihc[jt], hm))
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++) if (!mesh->GetMarker(jhc[kt], hm))
									{
										if ((invert ^ mesh->GetMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 )// these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if (!mesh->GetMarker(hc[it], hm) && (invert ^ mesh->GetMarker(hc[it], mask))) ret++;
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++) if (!mesh->GetMarker(lc[it], hm))
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++) if (!mesh->GetMarker(ilc[jt], hm))
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++) if (!mesh->GetMarker(jlc[kt], hm))
									{
										if ((invert ^ mesh->GetMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
									}
								}
							}
							ret += static_cast<enumerator>(result.size());
							if (!result.empty()) mesh->RemPrivateMarkerArray(&result[0], (enumerator)result.size(), mrk);
							result.clear();
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
		}
		return ret;
	}
	
	ElementArray<Element> Element::getAdjElements(ElementType _etype) const 
	{
		if( GetHandleElementType(GetHandle()) == ESET ) return getAsSet()->getAdjElements(_etype);
		assert( !(_etype & MESH) );
		assert( !(_etype & ESET) );
		INMOST::Mesh * mesh = GetMeshLink();
		ElementArray<Element> result(mesh);
		INMOST_DATA_INTEGER_TYPE conn[4] = {0,0,0,0};
		INMOST_DATA_INTEGER_TYPE myconn = -1, i = 0;
		
		for(ElementType e = NODE; e <= CELL; e = NextElementType(e))
		{
			if( _etype & e ) conn[i] = 1;
			if( GetElementType() & e ) myconn = i;
			i++;
		}
		if( !mesh->HideMarker() )
		{
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					result.push_back(GetHandle());
				}
				else if( i == myconn + 1 )
				{
					adj_type const & hc = mesh->HighConn(GetHandle());
					result.insert(result.end(),hc.begin(),hc.end());
				}
				else if( i == myconn - 1 )
				{
					adj_type const & lc = mesh->LowConn(GetHandle());
					result.insert(result.end(),lc.begin(),lc.end());
				}
				else if( i == myconn + 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & hc = mesh->HighConn(GetHandle());
					for(adj_type::size_type it = 0; it < hc.size(); it++)
					{
						adj_type const & ihc = mesh->HighConn(hc[it]);
						for(adj_type::size_type jt = 0; jt < ihc.size(); jt++)
							if( !mesh->GetPrivateMarker(ihc[jt],mrk) )
							{
								result.push_back(ihc[jt]);
								mesh->SetPrivateMarker(ihc[jt],mrk);
							}
					}
					result.RemPrivateMarker(mrk);
					mesh->ReleasePrivateMarker(mrk);
				}
				else if( i == myconn - 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & lc = mesh->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++)
					{
						adj_type const & ilc = mesh->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt < ilc.size(); jt++)
							if( !mesh->GetPrivateMarker(ilc[jt],mrk) )
							{
								result.push_back(ilc[jt]);
								mesh->SetPrivateMarker(ilc[jt],mrk);
							}
					}
					result.RemPrivateMarker(mrk);
					mesh->ReleasePrivateMarker(mrk);
				}
				else if( i == myconn + 3 ) // these are cells of the node
				{
					if (mesh->LowConnTag().isDefined(NODE))
					{
						adj_type const& lc = mesh->LowConn(GetHandle());
						result.insert(result.end(), lc.begin(), lc.end());
					}
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& hc = mesh->HighConn(GetHandle());
						for (adj_type::size_type it = 0; it < hc.size(); it++)
						{
							adj_type const& ihc = mesh->HighConn(hc[it]);
							for (adj_type::size_type jt = 0; jt < ihc.size(); jt++)
							{
								adj_type const& jhc = mesh->HighConn(ihc[jt]);
								for (adj_type::size_type kt = 0; kt < jhc.size(); kt++)
									if (!mesh->GetPrivateMarker(jhc[kt], mrk))
									{
										result.push_back(jhc[kt]);
										mesh->SetPrivateMarker(jhc[kt], mrk);
									}
							}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
				}
				else if( i == myconn - 3 ) // these are nodes of the cell
				{
					if (mesh->HighConnTag().isDefined(CELL))
					{
						adj_type const& hc = mesh->HighConn(GetHandle());
						result.insert(result.end(), hc.begin(), hc.end());
					}
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& lc = mesh->LowConn(GetHandle());
						for (adj_type::size_type it = 0; it < lc.size(); it++)
						{
							adj_type const& ilc = mesh->LowConn(lc[it]);
							for (adj_type::size_type jt = 0; jt < ilc.size(); jt++)
							{
								adj_type const& jlc = mesh->LowConn(ilc[jt]);
								for (adj_type::size_type kt = 0; kt < jlc.size(); kt++)
									if (!mesh->GetPrivateMarker(jlc[kt], mrk))
									{
										result.push_back(jlc[kt]);
										mesh->SetPrivateMarker(jlc[kt], mrk);
									}
							}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
				}
			}
		}
		else
		{
			MarkerType hm = mesh->HideMarker();
			for(i = 0; i < 4; i++) if( conn[i] )
			{
				if( i == myconn )
				{
					if( !GetMarker(hm) )
						result.push_back(GetHandle());
				}
				else if( i == myconn + 1 )
				{
					adj_type const & hc = mesh->HighConn(GetHandle());
					for(adj_type::size_type it = 0; it < hc.size(); ++it)
						if( !mesh->GetMarker(hc[it],hm) ) result.push_back(hc[it]);
					
				}
				else if( i == myconn - 1 )
				{
					adj_type const & lc = mesh->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); ++it)
						if( !mesh->GetMarker(lc[it],hm) ) result.push_back(lc[it]);
				}
				else if( i == myconn + 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & hc = mesh->HighConn(GetHandle());
					for(adj_type::size_type it = 0; it < hc.size(); it++) if( !mesh->GetMarker(hc[it],hm) )
					{
						adj_type const & ihc = mesh->HighConn(hc[it]);
						for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !mesh->GetMarker(ihc[jt],hm) )
							if( !mesh->GetPrivateMarker(ihc[jt],mrk) )
							{
								result.push_back(ihc[jt]);
								mesh->SetPrivateMarker(ihc[jt],mrk);
							}
					}
					result.RemPrivateMarker(mrk);
					mesh->ReleasePrivateMarker(mrk);
				}
				else if( i == myconn - 2 )
				{
					MarkerType mrk = mesh->CreatePrivateMarker();
					adj_type const & lc = mesh->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++) if( !mesh->GetMarker(lc[it],hm) )
					{
						adj_type const & ilc = mesh->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt != ilc.size(); jt++) if( !mesh->GetMarker(ilc[jt],hm) )
							if( !mesh->GetPrivateMarker(ilc[jt],mrk) )
							{
								result.push_back(ilc[jt]);
								mesh->SetPrivateMarker(ilc[jt],mrk);
							}
					}
					result.RemPrivateMarker(mrk);
					mesh->ReleasePrivateMarker(mrk);
				}
				else if( i == myconn + 3 )// these are cells of the node
				{
					if (mesh->LowConnTag().isDefined(NODE))
					{
						adj_type const& lc = mesh->LowConn(GetHandle());
						for (adj_type::size_type it = 0; it < lc.size(); ++it)
							if (!mesh->GetMarker(lc[it], hm)) result.push_back(lc[it]);
					}
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& hc = mesh->HighConn(GetHandle());
						for (adj_type::size_type it = 0; it < hc.size(); it++) if (!mesh->GetMarker(hc[it], hm))
						{
							adj_type const& ihc = mesh->HighConn(hc[it]);
							for (adj_type::size_type jt = 0; jt < ihc.size(); jt++) if (!mesh->GetMarker(ihc[jt], hm))
							{
								adj_type const& jhc = mesh->HighConn(ihc[jt]);
								for (adj_type::size_type kt = 0; kt < jhc.size(); kt++) if (!mesh->GetMarker(jhc[kt], hm))
								{
									if (!mesh->GetPrivateMarker(jhc[kt], mrk))
									{
										result.push_back(jhc[kt]);
										mesh->SetPrivateMarker(jhc[kt], mrk);
									}
								}
							}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
				}
				else if( i == myconn - 3 ) // these are nodes of the cell
				{
					if (mesh->HighConnTag().isDefined(CELL))
					{
						adj_type const& hc = mesh->HighConn(GetHandle());
						for (adj_type::size_type it = 0; it < hc.size(); ++it)
							if (!mesh->GetMarker(hc[it], hm)) result.push_back(hc[it]);
					}
					else
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const& lc = mesh->LowConn(GetHandle());
						for (adj_type::size_type it = 0; it < lc.size(); it++) if (!mesh->GetMarker(lc[it], hm))
						{
							adj_type const& ilc = mesh->LowConn(lc[it]);
							for (adj_type::size_type jt = 0; jt != ilc.size(); jt++) if (!mesh->GetMarker(ilc[jt], hm))
							{
								adj_type const& jlc = mesh->LowConn(ilc[jt]);
								for (adj_type::size_type kt = 0; kt < jlc.size(); kt++) if (!mesh->GetMarker(jlc[kt], hm))
								{
									if (!mesh->GetPrivateMarker(jlc[kt], mrk))
									{
										result.push_back(jlc[kt]);
										mesh->SetPrivateMarker(jlc[kt], mrk);
									}
								}
							}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
				}
			}
		}
		return result;
	}
	
	ElementArray<Element> Element::getAdjElements(ElementType _etype, MarkerType mask, bool invert) const 
	{
		if( GetHandleElementType(GetHandle()) == ESET ) return getAsSet()->getAdjElements(_etype,mask,invert);
		assert( !(_etype & MESH) );
		assert( !(_etype & ESET) );
		INMOST::Mesh * mesh = GetMeshLink();
		ElementArray<Element> result(mesh);
		INMOST_DATA_INTEGER_TYPE conn[4] = {0,0,0,0};
		INMOST_DATA_INTEGER_TYPE myconn = -1, i = 0;
		
		for(ElementType e = NODE; e <= CELL; e = NextElementType(e))
		{
			if( _etype & e ) conn[i] = 1;
			if( GetElementType() & e ) myconn = i;
			i++;
		}
		if( isPrivate(mask) )
		{
			if( !mesh->HideMarker() )
			{
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( invert ^ GetPrivateMarker(mask) ) result.push_back(GetHandle());
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it)
							if( invert ^ mesh->GetPrivateMarker(hc[it],mask) ) result.push_back(hc[it]);
					
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it)
							if( invert ^ mesh->GetPrivateMarker(lc[it],mask) ) result.push_back(lc[it]);
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++)
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								if( (invert ^ mesh->GetPrivateMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++)
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								if( (invert ^ mesh->GetPrivateMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 ) // these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if (invert ^ mesh->GetPrivateMarker(lc[it], mask)) result.push_back(lc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++)
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++)
										if ((invert ^ mesh->GetPrivateMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 ) // these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if (invert ^ mesh->GetPrivateMarker(hc[it], mask)) result.push_back(hc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++)
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++)
										if ((invert ^ mesh->GetPrivateMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
			else
			{
				MarkerType hm = mesh->HideMarker();
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( !GetMarker(hm) && (invert ^ GetPrivateMarker(mask)) )
							result.push_back(GetHandle());
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it)
							if( (invert ^ mesh->GetPrivateMarker(hc[it],mask)) && !mesh->GetMarker(hc[it],hm) ) result.push_back(hc[it]);
					
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it)
							if( (invert ^ mesh->GetPrivateMarker(lc[it],mask)) && !mesh->GetMarker(lc[it],hm) ) result.push_back(lc[it]);
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++) if( !mesh->GetMarker(hc[it],hm) )
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !mesh->GetMarker(ihc[jt],hm) )
								if( (invert ^ mesh->GetPrivateMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++) if( !mesh->GetMarker(lc[it],hm) )
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !mesh->GetMarker(ilc[jt],hm) )
								if( (invert ^ mesh->GetPrivateMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 ) // these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if ((invert ^ mesh->GetPrivateMarker(lc[it], mask)) && !mesh->GetMarker(lc[it], hm)) result.push_back(lc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++) if (!mesh->GetMarker(hc[it], hm))
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++) if (!mesh->GetMarker(ihc[jt], hm))
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++) if (!mesh->GetMarker(jhc[kt], hm))
									{
										if ((invert ^ mesh->GetPrivateMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
									}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 ) // these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if ((invert ^ mesh->GetPrivateMarker(hc[it], mask)) && !mesh->GetMarker(hc[it], hm)) result.push_back(hc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++) if (!mesh->GetMarker(lc[it], hm))
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++) if (!mesh->GetMarker(ilc[jt], hm))
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++) if (!mesh->GetMarker(jlc[kt], hm))
									{
										if ((invert ^ mesh->GetPrivateMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
									}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
		}
		else
		{
			if( !mesh->HideMarker() )
			{
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( invert ^ GetMarker(mask) ) result.push_back(GetHandle());
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it)
							if( invert ^ mesh->GetMarker(hc[it],mask) ) result.push_back(hc[it]);
					
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it)
							if( invert ^ mesh->GetMarker(lc[it],mask) ) result.push_back(lc[it]);
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++)
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								if( (invert ^ mesh->GetMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++)
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								if( (invert ^ mesh->GetMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 )// these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if (invert ^ mesh->GetMarker(lc[it], mask)) result.push_back(lc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++)
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++)
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++)
										if ((invert ^ mesh->GetMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 )// these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if (invert ^ mesh->GetMarker(hc[it], mask)) result.push_back(hc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++)
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++)
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++)
										if ((invert ^ mesh->GetMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
			else
			{
				MarkerType hm = mesh->HideMarker();
				for(i = 0; i < 4; i++) if( conn[i] )
				{
					if( i == myconn )
					{
						if( !GetMarker(hm) && (invert ^ GetMarker(mask)) )
							result.push_back(GetHandle());
					}
					else if( i == myconn + 1 )
					{
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); ++it)
							if( (invert ^ mesh->GetMarker(hc[it],mask)) && !mesh->GetMarker(hc[it],hm) ) result.push_back(hc[it]);
					
					}
					else if( i == myconn - 1 )
					{
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); ++it)
							if( (invert ^ mesh->GetMarker(lc[it],mask)) && !mesh->GetMarker(lc[it],hm) ) result.push_back(lc[it]);
					}
					else if( i == myconn + 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & hc = mesh->HighConn(GetHandle());
						for(adj_type::size_type it = 0; it < hc.size(); it++) if( !mesh->GetMarker(hc[it],hm) )
						{
							adj_type const & ihc = mesh->HighConn(hc[it]);
							for(adj_type::size_type jt = 0; jt < ihc.size(); jt++) if( !mesh->GetMarker(ihc[jt],hm) )
								if( (invert ^ mesh->GetMarker(ihc[jt],mask)) && !mesh->GetPrivateMarker(ihc[jt],mrk) )
								{
									result.push_back(ihc[jt]);
									mesh->SetPrivateMarker(ihc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn - 2 )
					{
						MarkerType mrk = mesh->CreatePrivateMarker();
						adj_type const & lc = mesh->LowConn(GetHandle());
						for(adj_type::size_type it = 0; it < lc.size(); it++) if( !mesh->GetMarker(lc[it],hm) )
						{
							adj_type const & ilc = mesh->LowConn(lc[it]);
							for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !mesh->GetMarker(ilc[jt],hm) )
								if( (invert ^ mesh->GetMarker(ilc[jt],mask)) && !mesh->GetPrivateMarker(ilc[jt],mrk) )
								{
									result.push_back(ilc[jt]);
									mesh->SetPrivateMarker(ilc[jt],mrk);
								}
						}
						result.RemPrivateMarker(mrk);
						mesh->ReleasePrivateMarker(mrk);
					}
					else if( i == myconn + 3 )// these are cells of the node
					{
						if (mesh->LowConnTag().isDefined(NODE))
						{
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); ++it)
								if ((invert ^ mesh->GetMarker(lc[it], mask)) && !mesh->GetMarker(lc[it], hm)) result.push_back(lc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); it++) if (!mesh->GetMarker(hc[it], hm))
							{
								adj_type const& ihc = mesh->HighConn(hc[it]);
								for (adj_type::size_type jt = 0; jt < ihc.size(); jt++) if (!mesh->GetMarker(ihc[jt], hm))
								{
									adj_type const& jhc = mesh->HighConn(ihc[jt]);
									for (adj_type::size_type kt = 0; kt < jhc.size(); kt++) if (!mesh->GetMarker(jhc[kt], hm))
									{
										if ((invert ^ mesh->GetMarker(jhc[kt], mask)) && !mesh->GetPrivateMarker(jhc[kt], mrk))
										{
											result.push_back(jhc[kt]);
											mesh->SetPrivateMarker(jhc[kt], mrk);
										}
									}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
					else if( i == myconn - 3 )// these are nodes of the cell
					{
						if (mesh->HighConnTag().isDefined(CELL))
						{
							adj_type const& hc = mesh->HighConn(GetHandle());
							for (adj_type::size_type it = 0; it < hc.size(); ++it)
								if ((invert ^ mesh->GetMarker(hc[it], mask)) && !mesh->GetMarker(hc[it], hm)) result.push_back(hc[it]);
						}
						else
						{
							MarkerType mrk = mesh->CreatePrivateMarker();
							adj_type const& lc = mesh->LowConn(GetHandle());
							for (adj_type::size_type it = 0; it < lc.size(); it++) if (!mesh->GetMarker(lc[it], hm))
							{
								adj_type const& ilc = mesh->LowConn(lc[it]);
								for (adj_type::size_type jt = 0; jt < ilc.size(); jt++) if (!mesh->GetMarker(ilc[jt], hm))
								{
									adj_type const& jlc = mesh->LowConn(ilc[jt]);
									for (adj_type::size_type kt = 0; kt < jlc.size(); kt++) if (!mesh->GetMarker(jlc[kt], hm))
									{
										if ((invert ^ mesh->GetMarker(jlc[kt], mask)) && !mesh->GetPrivateMarker(jlc[kt], mrk))
										{
											result.push_back(jlc[kt]);
											mesh->SetPrivateMarker(jlc[kt], mrk);
										}
									}
								}
							}
							result.RemPrivateMarker(mrk);
							mesh->ReleasePrivateMarker(mrk);
						}
					}
				}
			}
		}
		return result;
	}
	
	ElementArray<Node> Element::getNodes() const
	{
		switch(GetElementType())
		{
		case NODE: return ElementArray<Node>(GetMeshLink(),1,GetHandle());
		case EDGE: return getAsEdge()->getNodes();
		case FACE: return getAsFace()->getNodes();
		case CELL: return getAsCell()->getNodes();
		case ESET: return getAsSet()->getNodes();
		}
		return ElementArray<Node>(NULL);
	}
	ElementArray<Edge> Element::getEdges() const
	{
		switch(GetElementType())
		{
		case NODE: return getAsNode()->getEdges();
		case EDGE: return ElementArray<Edge>(GetMeshLink(),1,GetHandle());
		case FACE: return getAsFace()->getEdges();
		case CELL: return getAsCell()->getEdges();
		case ESET: return getAsSet()->getEdges();
		}
		return ElementArray<Edge>(NULL);
	}
	ElementArray<Face> Element::getFaces() const
	{
		switch(GetElementType())
		{
		case NODE: return getAsNode()->getFaces();
		case EDGE: return getAsEdge()->getFaces();
		case FACE: return ElementArray<Face>(GetMeshLink(),1,GetHandle());
		case CELL: return getAsCell()->getFaces();
		case ESET: return getAsSet()->getFaces();
		}
		return ElementArray<Face>(NULL);
	}
	ElementArray<Cell> Element::getCells() const
	{
		switch(GetElementType())
		{
		case NODE: return getAsNode()->getCells();
		case EDGE: return getAsEdge()->getCells();
		case FACE: return getAsFace()->getCells();
		case CELL: return ElementArray<Cell>(GetMeshLink(),1,GetHandle());
		case ESET: return getAsSet()->getCells();
		}
		return ElementArray<Cell>(NULL);
	}
	
	ElementArray<Node> Element::getNodes(MarkerType mask,bool invert) const
	{
		switch(GetElementType())
		{
		case NODE: return ElementArray<Node>(GetMeshLink(),(invert ^ (isPrivate(mask) ? GetPrivateMarker(mask) : GetMarker(mask)) ? 1 : 0),GetHandle());
		case EDGE: return getAsEdge()->getNodes(mask,invert);
		case FACE: return getAsFace()->getNodes(mask,invert);
		case CELL: return getAsCell()->getNodes(mask,invert);
		case ESET: return getAsSet()->getNodes(mask,invert);
		}
		return ElementArray<Node>(NULL);
	}
	ElementArray<Edge> Element::getEdges(MarkerType mask,bool invert) const
	{
		switch(GetElementType())
		{
		case NODE: return getAsNode()->getEdges(mask,invert);
		case EDGE: return ElementArray<Edge>(GetMeshLink(),(invert ^ (isPrivate(mask) ? GetPrivateMarker(mask) : GetMarker(mask)) ? 1 : 0),GetHandle());
		case FACE: return getAsFace()->getEdges(mask,invert);
		case CELL: return getAsCell()->getEdges(mask,invert);
		case ESET: return getAsSet()->getEdges(mask,invert);
		}
		return ElementArray<Edge>(NULL);	
	}
	ElementArray<Face> Element::getFaces(MarkerType mask,bool invert) const
	{
		switch(GetElementType())
		{
		case NODE: return getAsNode()->getFaces(mask,invert);
		case EDGE: return getAsEdge()->getFaces(mask,invert);
		case FACE: return ElementArray<Face>(GetMeshLink(),(invert ^ (isPrivate(mask) ? GetPrivateMarker(mask) : GetMarker(mask)) ? 1 : 0),GetHandle());
		case CELL: return getAsCell()->getFaces(mask,invert);
		case ESET: return getAsSet()->getFaces(mask,invert);
		}
		return ElementArray<Face>(NULL);
	}
	ElementArray<Cell> Element::getCells(MarkerType mask,bool invert) const
	{
		switch(GetElementType())
		{
		case NODE: return getAsNode()->getCells(mask,invert);
		case EDGE: return getAsEdge()->getCells(mask,invert);
		case FACE: return getAsFace()->getCells(mask,invert);
		case CELL: return ElementArray<Cell>(GetMeshLink(),(invert ^ (isPrivate(mask) ? GetPrivateMarker(mask) : GetMarker(mask)) ? 1 : 0),GetHandle());
		case ESET: return getAsSet()->getCells(mask,invert);
		}
		return ElementArray<Cell>(NULL);
	}
	
	bool Element::CheckElementConnectivity() const
	{
		Mesh * mesh = GetMeshLink();
		HandleType me = GetHandle();
		//MarkerType hm = mesh->HideMarker();
		//hidden element should maintain connections to hidden and unhidden elements
		//both hidden and unhidden elements should maintain connection back to hidden element
		//
		//for unhidden element
		
		if( GetElementType() < CELL || mesh->HighConnTag().isDefined(CELL) )
		{
			adj_type const & hc = mesh->HighConn(GetHandle());
			for(adj_type::size_type jt = 0; jt < hc.size(); jt++) //iterate over upper adjacent
			{
				if( hc[jt] == InvalidHandle() )
				{
					std::cout << "Invalid connection from ";
					std::cout << ElementTypeName(GetElementType()) << ":" << LocalID();
					std::cout << " in HighConn" << std::endl;
					return false;
				}
				int found = 0;
				adj_type const & ilc = mesh->LowConn(hc[jt]);
				for(adj_type::size_type kt = 0; kt < ilc.size(); kt++) //search for the link to me
					if( ilc[kt] == me ) found++;
				if( !found )
				{
					std::cout << "Not found connection from ";
					std::cout << ElementTypeName(GetHandleElementType(hc[jt])) << ":" << GetHandleID(hc[jt]);
					std::cout << " to ";
					std::cout << ElementTypeName(GetElementType()) << ":" << LocalID();
					std::cout << std::endl;
					return false;
				}
				else if( found > 1 )
				{
					std::cout << "Found " << found << " connections from ";
					std::cout << ElementTypeName(GetHandleElementType(hc[jt])) << ":" << GetHandleID(hc[jt]);
					std::cout << " to ";
					std::cout << ElementTypeName(GetElementType()) << ":" << LocalID();
					std::cout << std::endl;
					return false;
				}
			}
		}
		if( GetElementType() > NODE || mesh->LowConnTag().isDefined(NODE) )
		{
			adj_type const & lc = mesh->LowConn(GetHandle());
			for(adj_type::size_type jt = 0; jt < lc.size(); jt++) //iterate over lower adjacent
			{
				if( lc[jt] == InvalidHandle() )
				{
					std::cout << "Invalid connection from ";
					std::cout << ElementTypeName(GetElementType()) << ":" << LocalID();
					std::cout << " in LowConn" << std::endl;
					return false;
				}
				int found = 0;
				adj_type const & ihc = mesh->HighConn(lc[jt]);
				for(adj_type::size_type kt = 0; kt < ihc.size(); kt++) //search for the link to me
					if( ihc[kt] == me ) found++;
				if( !found )
				{
					std::cout << "Not found connection from ";
					std::cout << ElementTypeName(GetHandleElementType(lc[jt])) << ":" << GetHandleID(lc[jt]);
					std::cout << " to ";
					std::cout << ElementTypeName(GetElementType()) << ":" << LocalID();
					std::cout << std::endl;
					return false;
				}
				else if( found > 1 )
				{
					std::cout << "Found " << found << " connections from ";
					std::cout << ElementTypeName(GetHandleElementType(lc[jt])) << ":" << GetHandleID(lc[jt]);
					std::cout << " to ";
					std::cout << ElementTypeName(GetElementType()) << ":" << LocalID();
					std::cout << std::endl;
					return false;
				}
			}
		}
		return true;
	}

	void Element::PrintElementConnectivity() const
	{
		Mesh * mesh = GetMeshLink();
		HandleType me = GetHandle();
		std::cout << "Element " << ElementTypeName(GetElementType()) << " " << LocalID() << std::endl;
		if( GetElementType() < CELL || mesh->HighConnTag().isDefined(CELL) )
		{
			adj_type const & hc = mesh->HighConn(GetHandle());
			std::cout << "Upper adjacencies (" << hc.size() << "):" << std::endl;
			for(adj_type::size_type jt = 0; jt < hc.size(); jt++) //iterate over upper adjacent
			{
				std::cout << "[" << std::setw(3) << jt << "] " << ElementTypeName(GetHandleElementType(hc[jt])) << " " << GetHandleID(hc[jt]) << " lower ";
				bool found = false;
				adj_type const & ilc = mesh->LowConn(hc[jt]);
				std::cout << "(" << ilc.size() << "): ";
				for(adj_type::size_type kt = 0; kt < ilc.size(); kt++) //search for the link to me
				{
					std::cout << ElementTypeName(GetHandleElementType(ilc[kt])) << " " << GetHandleID(ilc[kt]) << " ";
					if( ilc[kt] == me ) found = true;
				}
				if( !found ) std::cout << " no me here! ";
				std::cout << std::endl;
			}
		}
		if( GetElementType() > NODE || mesh->LowConnTag().isDefined(NODE) )
		{
			adj_type const & lc = mesh->LowConn(GetHandle());
			std::cout << "Lower adjacencies (" << lc.size() << "):" << std::endl;
			for(adj_type::size_type jt = 0; jt < lc.size(); jt++) //iterate over lower adjacent
			{
				std::cout << "[" <<  std::setw(3) <<  jt << "] " << ElementTypeName(GetHandleElementType(lc[jt])) << " " << GetHandleID(lc[jt]) << " higher ";
				bool found = false;
				adj_type const & ihc = mesh->HighConn(lc[jt]);
				std::cout << "(" << ihc.size() << "): ";
				for(adj_type::size_type kt = 0; kt < ihc.size(); kt++) //search for the link to me
				{
					std::cout << ElementTypeName(GetHandleElementType(ihc[kt])) << " " << GetHandleID(ihc[kt]) << " ";
					if( ihc[kt] == me ) found = true;
				}
				if( !found ) std::cout << " no me here! ";
				std::cout << std::endl;
			}
		}
	}
	
	bool Element::CheckConnectivity(Mesh * m)
	{
		bool check = true;
		for(Mesh::iteratorElement it = m->BeginElement(CELL | FACE | EDGE | NODE); it != m->EndElement(); it++)
			if( !it->CheckElementConnectivity() ) check = false;
		return check;
	}
	
	
	
	ElementArray<Element> Element::BridgeAdjacencies(ElementType Bridge, ElementType Dest,  MarkerType bridge_mask, bool bridge_invert, MarkerType target_mask, bool target_invert) const
	{
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		ElementArray<Element> adjcells(m);
		ElementArray<Element> adjfaces = getAdjElements(Bridge);
		ElementArray<Element> my = Bridge & GetElementType() ? ElementArray<Element>(m) : getAdjElements(Dest);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->SetPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = adjfaces.begin(); it != adjfaces.end(); it++) if(bridge_mask == 0 || (bridge_invert ^ (isPrivate(bridge_mask) ? it->GetPrivateMarker(bridge_mask) : it->GetMarker(bridge_mask))))
		{
			ElementArray<Element> sub = it->getAdjElements(Dest);
			for(ElementArray<Element>::iterator jt = sub.begin(); jt != sub.end(); jt++)
				if( !jt->GetPrivateMarker(mrk) && (target_mask == 0 || (target_invert ^ (isPrivate(target_mask) ? jt->GetPrivateMarker(target_mask) : jt->GetMarker(target_mask)))) )
				{
					adjcells.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
		}
		for(ElementArray<Element>::iterator it = adjcells.begin(); it != adjcells.end(); it++) it->RemPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return adjcells;
	}


	ElementArray<Node> Element::BridgeAdjacencies2Node(ElementType Bridge, MarkerType bridge_mask, bool bridge_invert, MarkerType target_mask, bool target_invert) const
	{
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		ElementArray<Node> adjcells(m);
		ElementArray<Element> adjfaces = getAdjElements(Bridge);
		ElementArray<Element> my = Bridge & GetElementType() ? ElementArray<Element>(m) : getAdjElements(NODE);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->SetPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = adjfaces.begin(); it != adjfaces.end(); it++) if(bridge_mask == 0 || (bridge_invert ^ (isPrivate(bridge_mask) ? it->GetPrivateMarker(bridge_mask) : it->GetMarker(bridge_mask))))
		{
			ElementArray<Node> sub = it->getNodes();
			for(ElementArray<Node>::iterator jt = sub.begin(); jt != sub.end(); jt++)
				if( !jt->GetPrivateMarker(mrk) && (target_mask == 0 || (target_invert ^ (isPrivate(target_mask) ? jt->GetPrivateMarker(target_mask) : jt->GetMarker(target_mask)))) )
				{
					adjcells.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
		}
		for(ElementArray<Node>::iterator it = adjcells.begin(); it != adjcells.end(); it++) it->RemPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return adjcells;
	}

	ElementArray<Edge> Element::BridgeAdjacencies2Edge(ElementType Bridge,  MarkerType bridge_mask, bool bridge_invert, MarkerType target_mask, bool target_invert) const
	{
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		ElementArray<Edge> adjcells(m);
		ElementArray<Element> adjfaces = getAdjElements(Bridge);
		ElementArray<Element> my = Bridge & GetElementType() ? ElementArray<Element>(m) : getAdjElements(EDGE);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->SetPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = adjfaces.begin(); it != adjfaces.end(); it++) if(bridge_mask == 0 || (bridge_invert ^ (isPrivate(bridge_mask) ? it->GetPrivateMarker(bridge_mask) : it->GetMarker(bridge_mask))))
		{
			ElementArray<Edge> sub = it->getEdges();
			for(ElementArray<Edge>::iterator jt = sub.begin(); jt != sub.end(); jt++)
				if( !jt->GetPrivateMarker(mrk) && (target_mask == 0 || (target_invert ^ (isPrivate(target_mask) ? jt->GetPrivateMarker(target_mask) : jt->GetMarker(target_mask)))) )
				{
					adjcells.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
		}
		for(ElementArray<Edge>::iterator it = adjcells.begin(); it != adjcells.end(); it++) it->RemPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return adjcells;
	}

	ElementArray<Face> Element::BridgeAdjacencies2Face(ElementType Bridge, MarkerType bridge_mask, bool bridge_invert, MarkerType target_mask, bool target_invert) const
	{
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		ElementArray<Face> adjcells(m);
		ElementArray<Element> adjfaces = getAdjElements(Bridge);
		ElementArray<Element> my = Bridge & GetElementType() ? ElementArray<Element>(m) : getAdjElements(FACE);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->SetPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = adjfaces.begin(); it != adjfaces.end(); it++) if(bridge_mask == 0 || (bridge_invert ^ (isPrivate(bridge_mask) ? it->GetPrivateMarker(bridge_mask) : it->GetMarker(bridge_mask))))
		{
			ElementArray<Face> sub = it->getFaces();
			for(ElementArray<Face>::iterator jt = sub.begin(); jt != sub.end(); jt++)
				if( !jt->GetPrivateMarker(mrk) && (target_mask == 0 || (target_invert ^ (isPrivate(target_mask) ? jt->GetPrivateMarker(target_mask) : jt->GetMarker(target_mask)))) )
				{
					adjcells.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
		}
		for(ElementArray<Face>::iterator it = adjcells.begin(); it != adjcells.end(); it++) it->RemPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return adjcells;
	}

	ElementArray<Cell> Element::BridgeAdjacencies2Cell(ElementType Bridge,  MarkerType bridge_mask, bool bridge_invert, MarkerType target_mask, bool target_invert) const
	{
		Mesh * m = GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		ElementArray<Cell> adjcells(m);
		ElementArray<Element> adjfaces = getAdjElements(Bridge);
		ElementArray<Element> my = Bridge & GetElementType() ? ElementArray<Element>(m) : getAdjElements(CELL);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->SetPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = adjfaces.begin(); it != adjfaces.end(); it++) if(bridge_mask == 0 || (bridge_invert ^ (isPrivate(bridge_mask) ? it->GetPrivateMarker(bridge_mask) : it->GetMarker(bridge_mask))))
		{
			ElementArray<Cell> sub = it->getCells();
			for(ElementArray<Cell>::iterator jt = sub.begin(); jt != sub.end(); jt++)
				if( !jt->GetPrivateMarker(mrk) && (target_mask == 0 || (target_invert ^ (isPrivate(target_mask) ? jt->GetPrivateMarker(target_mask) : jt->GetMarker(target_mask)))) )
				{
					adjcells.push_back(*jt);
					jt->SetPrivateMarker(mrk);
				}
		}
		for(ElementArray<Cell>::iterator it = adjcells.begin(); it != adjcells.end(); it++) it->RemPrivateMarker(mrk);
		for(ElementArray<Element>::iterator it = my.begin(); it != my.end(); it++) it->RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return adjcells;
	}


	Element::Status Element::GetStatus() const
	{
		assert(isValid());
		return GetMeshLink()->GetStatus(GetHandle());
	}
	void Element::SetStatus(Status status) const
	{
		GetMeshLink()->SetStatus(GetHandle(),status);
	}
	Storage::integer & Element::GlobalID() const
	{
		return GetMeshLink()->GlobalID(GetHandle());
	}
	void Element::Centroid  (Storage::real * cnt) const
	{
		GetMeshLink()->GetGeometricData(GetHandle(),CENTROID,cnt);
	}
	void Element::Barycenter(Storage::real * cnt) const
	{
		GetMeshLink()->GetGeometricData(GetHandle(),BARYCENTER,cnt);
	}
	/*
	Element::adj_type &       Element::HighConn() 
	{
		assert(isValid()); 
		return GetMeshLink()->HighConn(GetHandle());
	}
	Element::adj_type &       Element::LowConn () 
	{
		assert(isValid()); 
		return GetMeshLink()->LowConn(GetHandle());
	}
	Element::adj_type const & Element::HighConn() const 
	{
		assert(isValid()); 
		return GetMeshLink()->HighConn(GetHandle());
	}
	Element::adj_type const & Element::LowConn () const 
	{
		assert(isValid()); 
		return GetMeshLink()->LowConn(GetHandle());
	}
	*/
	Element::GeometricType Element::GetGeometricType() const 
	{
		assert(isValid()); 
		return GetMeshLink()->GetGeometricType(GetHandle());
	}
	void Element::SetGeometricType(GeometricType t)
	{
		assert(isValid()); 
		GetMeshLink()->SetGeometricType(GetHandle(),t); 
	}

	template<typename InputContainer>
	void UnionSendTo(const Element * e, InputContainer cont)
	{
		Storage::integer_array set_procs = e->IntegerArray(e->GetMeshLink()->ProcessorsTag());
		Storage::integer_array sendto = e->IntegerArray(e->GetMeshLink()->SendtoTag());
		std::sort(sendto.begin(),sendto.end());
		std::vector<Storage::integer> tmp1(cont.size()),tmp2;
		tmp1.resize(std::set_difference(cont.begin(),cont.end(),set_procs.begin(),set_procs.end(),tmp1.begin())-tmp1.begin());
		tmp2.resize(tmp1.size()+sendto.size());
		tmp2.resize(std::set_union(tmp1.begin(),tmp1.end(),sendto.begin(),sendto.end(),tmp2.begin())-tmp2.begin());
		sendto.replace(sendto.begin(),sendto.end(),tmp2.begin(),tmp2.end());
	}

	void Element::SendTo(std::set<Storage::integer> & procs) const
	{
		if( GetMeshLink()->GetMeshState() != Mesh::Serial ) UnionSendTo(this,procs);
	}	

	void Element::SendTo(std::vector<Storage::integer> & procs) const
	{
		if( GetMeshLink()->GetMeshState() != Mesh::Serial ) UnionSendTo(this,procs);
	}

	void Element::SendTo(Storage::integer_array procs) const
	{
		if( GetMeshLink()->GetMeshState() != Mesh::Serial ) UnionSendTo(this,procs);
	}

}
#endif
