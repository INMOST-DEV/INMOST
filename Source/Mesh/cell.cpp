#include "inmost.h"
#if defined(USE_MESH)

namespace INMOST
{
	
	
	bool Cell::CheckEdgeOrder() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 )
		{
			if( !m->HideMarker() )
			{
				HandleType last, first;
				adj_type const & lc = m->LowConn(GetHandle());
				if( lc.size() < 3 ) return false;
				HandleType q = lc[0]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				first = qlc[0]; //node 0
				last = qlc[1]; //node 1
				HandleType r = lc[1]; //edge 1
				adj_type & rlc = m->LowConn(r);
				if( first == rlc[0] || first == rlc[1] )
				{
					HandleType temp = first;
					first = last;
					last = temp;
				}
				else if ( !(last == rlc[0] || last == rlc[1]) )
					return false;
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it < iend) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					if( last == ilc[0] ) last = ilc[1];
					else if( last == ilc[1] ) last = ilc[0];
					else return false;
					++it;
				}
			}
			else
			{
				MarkerType hm = m->HideMarker();
				HandleType first, last;
				enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
				adj_type const & lc = m->LowConn(GetHandle());
				if( m->Count(lc.data(),lc.size(),hm) < 3 ) return false;
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType q = lc[i]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				first = qlc[k]; //node 0
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				last = qlc[k]; //node 1
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType r = lc[i]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				k1 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				k2 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				if( first == rlc[k1] || first == rlc[k2] )
				{
					HandleType temp = first;
					first = last;
					last = temp;
				}
				else if ( !(last == rlc[k1] || last == rlc[k2]) )
					return false;
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it != iend)
				{
					if( !m->GetMarker(lc[it],hm) ) //loop over edges
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						k1 = ENUMUNDEF; 
						k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
						k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
						if( last == ilc[k1] ) last = ilc[k2];
						else if( last == ilc[k2] ) last = ilc[k1];
						else return false;
					}
					++it;
				}
			}
		}
		return true;
	}
	
	bool Cell::FixEdgeOrder() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 )
		{
			if( !m->HideMarker() )
			{
				HandleType last, first;
				adj_type & lc = m->LowConn(GetHandle());
				if( lc.size() < 3 ) return false;
				HandleType q = lc[0]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				if( qlc.size() != 2 ) return false;
				first = qlc[0]; //node 0
				last = qlc[1]; //node 1
				HandleType r = lc[1]; //edge 1
				adj_type & rlc = m->LowConn(r);
				if( rlc.size() != 2 ) return false;
				if( first == rlc[0] || first == rlc[1] )
				{
					HandleType temp = first;
					first = last;
					last = temp;
				}
				else if ( !(last == rlc[0] || last == rlc[1]) )
				{
					adj_type::size_type jt = 2, jend = lc.size();
					while(jt < jend)
					{
						adj_type const & ilc = m->LowConn(lc[jt]);
						if( ilc.size() != 2 ) return false;
						if( first == ilc[0] || first == ilc[1] )
						{
							HandleType temp = lc[1];
							lc[1] = lc[jt];
							lc[jt] = temp;
							temp = first;
							first = last;
							last = temp;
							break;
						}
						else if( last == ilc[0] || last == ilc[1] )
						{
							HandleType temp = lc[1];
							lc[1] = lc[jt];
							lc[jt] = temp;
							break;
						}
						++jt;
					}
					if( jt == jend ) return false; //no matching edge
				}
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it < iend) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					if( ilc.size() != 2 ) return false;
					if( last == ilc[0] ) 
					{
						last = ilc[1];
						++it;
					}
					else if( last == ilc[1] ) 
					{
						last = ilc[0];
						++it;
					}
					else //search for the connected edge and swap with current
					{
						adj_type::size_type jt = it+1, jend = lc.size();
						while(jt < jend)
						{
							adj_type const & ilc = m->LowConn(lc[jt]);
							if( ilc.size() != 2 ) return false;
							if( last == ilc[0] || last == ilc[1] )
							{
								HandleType temp = lc[it];
								lc[it] = lc[jt];
								lc[jt] = temp;
								break;
							}
						}
						if( jt == jend ) return false; //no matching edge
					}
				}
				//check that the loop is closed
				adj_type const & ilc = m->LowConn(lc[iend]);
				if( ilc.size() != 2 ) return false;
				if( !( (ilc[0] == last && ilc[1] == first) || (ilc[0] == first && ilc[1] == last) ) )
					return false;
			}
			else
			{
				HandleType first, last;
				enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
				MarkerType hm = m->HideMarker();
				adj_type & lc = m->LowConn(GetHandle());
				if( m->Count(lc.data(),lc.size(),hm) < 3 ) return false;
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType q = lc[i]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				if(m->Count(qlc.data(),qlc.size(),hm)!=2) return false;
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				first = qlc[k]; //node 0
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				last = qlc[k]; //node 1
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				HandleType r = lc[i]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				if(m->Count(rlc.data(),rlc.size(),hm)!=2) return false;
				k1 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				k2 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				if( first == rlc[k1] || first == rlc[k2] )
				{
					HandleType temp = first;
					first = last;
					last = temp;
				}
				else if( !(last == rlc[k1] || last == rlc[k2]) )
				{
					adj_type::size_type jt = i+1, jend = lc.size();
					while(jt < jend)
					{
						if( !m->GetMarker(lc[jt],hm) ) //loop over edges
						{
							adj_type const & ilc = m->LowConn(lc[jt]);
							if( m->Count(ilc.data(),ilc.size(),hm) != 2 ) return false;
							k1 = ENUMUNDEF; 
							k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							if( first == ilc[k1] || first == ilc[k2] )
							{
								HandleType temp = lc[i];
								lc[i] = lc[jt];
								lc[jt] = temp;
								temp = first;
								first = last;
								last = temp;
								break;
							}
							else if( last == ilc[k1] || last == ilc[k2] )
							{
								HandleType temp = lc[i];
								lc[i] = lc[jt];
								lc[jt] = temp;
								break;
							}
						}
						++jt;
					}
					if( jt == jend ) return false; //no matching edge
				}
				adj_type::size_type it = i, iend = lc.size()-1;
				while(it != iend) if( !m->GetMarker(lc[it],hm) ) //loop over edges
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					if(m->Count(ilc.data(),ilc.size(),hm)!=2) return false;
					k1 = ENUMUNDEF; 
					k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
					k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
					if( last == ilc[k1] ) 
					{
						last = ilc[k2];
						++it;
					}
					else if( last == ilc[k2] )
					{
						last = ilc[k1];
						++it;
					}
					else//search for the connected edge and swap with current
					{
						adj_type::size_type jt = it+1, jend = lc.size();
						while(jt < jend) if( !m->GetMarker(lc[jt],hm) ) //loop over edges
						{
							adj_type const & ilc = m->LowConn(lc[jt]);
							if(m->Count(ilc.data(),ilc.size(),hm)!=2) return false;
							k1 = ENUMUNDEF; 
							k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							if( last == ilc[k1] || last == ilc[k2] )
							{
								HandleType temp = lc[it];
								lc[it] = lc[jt];
								lc[jt] = temp;
								break;
							}
						}
						if( jt == jend ) return false; //no matching edge
					}
				} else ++it;
				//check that the loop is closed
				adj_type const & ilc = m->LowConn(lc[iend]);
				if(m->Count(ilc.data(),ilc.size(),hm)!=2) return false;
				k1 = ENUMUNDEF; 
				k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
				k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
				if( !( (ilc[k1] == last && ilc[k2] == first) || (ilc[k1] == first && ilc[k2] == last) ) )
					return false;
			}
		}
		return true;
	}
	
	ElementArray<Node> Cell::getNodes() const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		ElementArray<Node> ret(m);
		MarkerType mrk = m->CreatePrivateMarker();
		HandleType hc = GetHandle();
		if( !m->HideMarker() )
		{
			if (Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2) // This cell is 2d face
			{
				ElementArray<Edge> aret;
				adj_type const& lc = m->LowConn(GetHandle()); // set of edges
				assert(lc.size() > 2); // it should be at least triangle
				aret.reserve(lc.size());
				HandleType q = lc[0]; //edge 0
				adj_type const& qlc = m->LowConn(q);
				assert(qlc.size() == 2);
				aret.push_back(qlc[0]); //node 0
				aret.push_back(qlc[1]); //node 1
				HandleType r = lc[1]; //edge 1
				adj_type& rlc = m->LowConn(r);
				assert(rlc.size() == 2);
				if (aret.data()[0] == rlc[0] || aret.data()[0] == rlc[1])
				{
					HandleType temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::size_type it = 1, iend = lc.size() - 1;
				while (it < iend) //loop over edges
				{
					adj_type const& ilc = m->LowConn(lc[it]);
					assert(ilc.size() == 2);
					if (aret.atback() == ilc[0]) aret.push_back(ilc[1]);
					else aret.push_back(ilc[0]);
					++it;
				}
				ret.reserve(aret.size());
				for (ElementArray<Edge>::iterator it = aret.begin(); it != aret.end(); ++it)
					ret.push_back(m->LowConn(*it)[0]);
			}
			else
			{
				Element::adj_type& lc = m->LowConn(hc);
				for (Element::adj_type::size_type i = 0; i < lc.size(); i++) //iterate over faces
				{
					Element::adj_type& ilc = m->LowConn(lc[i]);
					for (Element::adj_type::size_type j = 0; j < ilc.size(); j++) //iterate over face edges
					{
						Element::adj_type& jlc = m->LowConn(ilc[j]);
						for (Element::adj_type::size_type k = 0; k < jlc.size(); k++) //iterator over edge nodes
						{
							if (!m->GetPrivateMarker(jlc[k], mrk))
							{
								m->SetPrivateMarker(jlc[k], mrk);
								ret.push_back(jlc[k]);
							}
						}
					}
				}
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			if (Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2) // This cell is 2d face
			{
				ElementArray<Edge> aret;
				enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
				adj_type const& lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				i = m->getNext(lc.data(), static_cast<enumerator>(lc.size()), i, hm);
				HandleType q = lc[i]; //edge 0
				adj_type const& qlc = m->LowConn(q);
				assert(m->Count(qlc.data(), qlc.size(), hm) == 2);
				k = m->getNext(qlc.data(), static_cast<enumerator>(qlc.size()), k, hm);
				aret.push_back(qlc[k]); //node 0
				k = m->getNext(qlc.data(), static_cast<enumerator>(qlc.size()), k, hm);
				aret.push_back(qlc[k]); //node 1
				i = m->getNext(lc.data(), static_cast<enumerator>(lc.size()), i, hm);
				HandleType r = lc[i]; //edge 1
				adj_type const& rlc = m->LowConn(r);
				assert(m->Count(rlc.data(), rlc.size(), hm) == 2);
				k1 = m->getNext(rlc.data(), static_cast<enumerator>(rlc.size()), k1, hm);
				k2 = m->getNext(rlc.data(), static_cast<enumerator>(rlc.size()), k1, hm);
				if (aret.data()[0] == rlc[k1] || aret.data()[0] == rlc[k2])
				{
					HandleType temp = aret.data()[0];
					aret.data()[0] = aret.data()[1];
					aret.data()[1] = temp;
				}
				adj_type::size_type it = i, iend = lc.size() - 1;

				while (m->GetMarker(lc[iend], hm) && iend > 0) iend--;

				while (it != iend)
				{
					if (!m->GetMarker(lc[it], hm)) //loop over edges
					{
						adj_type const& ilc = m->LowConn(lc[it]);
						assert(m->Count(ilc.data(), ilc.size(), hm) == 2);
						k1 = ENUMUNDEF;
						k1 = m->getNext(ilc.data(), static_cast<enumerator>(ilc.size()), k1, hm);
						k2 = m->getNext(ilc.data(), static_cast<enumerator>(ilc.size()), k1, hm);
						if (aret.atback() == ilc[k1])
							aret.push_back(ilc[k2]);
						else
							aret.push_back(ilc[k1]);
					}
					++it;
				}
				ret.reserve(aret.size());
				for (ElementArray<Edge>::iterator it = aret.begin(); it != aret.end(); ++it)
				{
					Element::adj_type adj = m->LowConn(*it);
					assert(m->Count(adj.data(), static_cast<enumerator>(adj.size()),hm) == 1);
					ret.push_back(adj[m->getNext(adj.data(), static_cast<enumerator>(adj.size()),ENUMUNDEF,hm)]);
				}
			}
			else
			{
				Element::adj_type& lc = m->LowConn(hc);
				for (Element::adj_type::size_type i = 0; i < lc.size(); i++) if (!m->GetMarker(lc[i], hm))  //iterate over faces
				{
					Element::adj_type& ilc = m->LowConn(lc[i]);
					for (Element::adj_type::size_type j = 0; j < ilc.size(); j++) if (!m->GetMarker(ilc[j], hm)) //iterate over face edges
					{
						Element::adj_type& jlc = m->LowConn(ilc[j]);
						for (Element::adj_type::size_type k = 0; k < jlc.size(); k++) if (!m->GetMarker(jlc[k], hm)) //iterator over edge nodes
						{
							if (!m->GetPrivateMarker(jlc[k], mrk))
							{
								m->SetPrivateMarker(jlc[k], mrk);
								ret.push_back(jlc[k]);
							}
						}
					}
				}
			}
		}
		ret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return ret;
	}


	ElementArray<Node> Cell::getNodes(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		ElementArray<Node> ret(m);
		MarkerType mrk = m->CreatePrivateMarker();
		HandleType hc = GetHandle();
		if( isPrivate(mask) )
		{
			if( !m->HideMarker() )
			{
				if (Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2) // This cell is 2d face
				{
					ElementArray<Edge> aret;
					adj_type const& lc = m->LowConn(GetHandle());
					aret.reserve(lc.size());
					HandleType q = lc[0];
					adj_type const& qlc = m->LowConn(q);
					HandleType first = qlc[0], last = qlc[1]; //edge 0
					if (invert ^ m->GetPrivateMarker(qlc[0], mask)) aret.push_back(qlc[0]); //node 0
					if (invert ^ m->GetPrivateMarker(qlc[1], mask)) aret.push_back(qlc[1]); //node 1
					HandleType r = lc[1]; //edge 1
					adj_type const& rlc = m->LowConn(r);
					if (first == rlc[0] || first == rlc[1])
					{
						last = first;
						if (aret.size() > 1)
						{
							HandleType temp = aret.data()[0];
							aret.data()[0] = aret.data()[1];
							aret.data()[1] = temp;
						}
					}
					adj_type::size_type it = 1, iend = lc.size() - 1;
					while (it < iend) //loop over edges
					{
						adj_type const& ilc = m->LowConn(lc[it]);
						if (last == ilc[0])
						{
							if (invert ^ m->GetPrivateMarker(ilc[1], mask))
								aret.push_back(ilc[1]);
							last = ilc[1];
						}
						else
						{
							if (invert ^ m->GetPrivateMarker(ilc[0], mask))
								aret.push_back(ilc[0]);
							last = ilc[0];
						}
						++it;
					}
					ret.reserve(aret.size());
					for (ElementArray<Edge>::iterator it = aret.begin(); it != aret.end(); ++it)
						ret.push_back(m->LowConn(*it)[0]);
				}
				else
				{
					Element::adj_type& lc = m->LowConn(hc);
					for (Element::adj_type::size_type i = 0; i < lc.size(); i++) //iterate over faces
					{
						Element::adj_type& ilc = m->LowConn(lc[i]);
						for (Element::adj_type::size_type j = 0; j < ilc.size(); j++) //iterate over face edges
						{
							Element::adj_type& jlc = m->LowConn(ilc[j]);
							for (Element::adj_type::size_type k = 0; k < jlc.size(); k++) //iterator over edge nodes
							{
								if ((invert ^ m->GetPrivateMarker(jlc[k], mask)) && !m->GetPrivateMarker(jlc[k], mrk))
								{
									m->SetPrivateMarker(jlc[k], mrk);
									ret.push_back(jlc[k]);
								}
							}
						}
					}
				}
			}
			else
			{
				MarkerType hm = m->HideMarker();
				if (Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2) // This cell is 2d face
				{
					ElementArray<Edge> aret;
					enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
					adj_type const& lc = m->LowConn(GetHandle());
					aret.reserve(lc.size());
					i = m->getNext(lc.data(), static_cast<enumerator>(lc.size()), i, hm);
					HandleType q = lc[i], first, last; //edge 0
					adj_type const& qlc = m->LowConn(q);
					k = m->getNext(qlc.data(), static_cast<enumerator>(qlc.size()), k, hm);
					if (invert ^ m->GetPrivateMarker(qlc[k], mask)) aret.push_back(qlc[k]); //node 0
					first = qlc[k];
					k = m->getNext(qlc.data(), static_cast<enumerator>(qlc.size()), k, hm);
					if (invert ^ m->GetPrivateMarker(qlc[k], mask)) aret.push_back(qlc[k]); //node 1
					last = qlc[k];
					i = m->getNext(lc.data(), static_cast<enumerator>(lc.size()), i, hm);
					HandleType r = lc[i]; //edge 1
					adj_type const& rlc = m->LowConn(r);
					k1 = m->getNext(rlc.data(), static_cast<enumerator>(rlc.size()), k1, hm);
					k2 = m->getNext(rlc.data(), static_cast<enumerator>(rlc.size()), k1, hm);
					if (first == rlc[k1] || first == rlc[k2])
					{
						last = first;
						if (aret.size() > 1)
						{
							HandleType temp = aret.data()[0];
							aret.data()[0] = aret.data()[1];
							aret.data()[1] = temp;
						}
					}
					adj_type::size_type it = 1, iend = lc.size() - 1;
					while (it < iend) if (!m->GetMarker(lc[it], hm)) //loop over edges
					{
						adj_type const& ilc = m->LowConn(lc[it]);
						k1 = ENUMUNDEF;
						k1 = m->getNext(ilc.data(), static_cast<enumerator>(ilc.size()), k1, hm);
						k2 = m->getNext(ilc.data(), static_cast<enumerator>(ilc.size()), k1, hm);
						if (last == ilc[k1])
						{
							if (invert ^ m->GetPrivateMarker(ilc[k2], mask))
								aret.push_back(ilc[k2]);
							last = ilc[k2];
						}
						else
						{
							if (invert ^ m->GetPrivateMarker(ilc[k1], mask))
								aret.push_back(ilc[k1]);
							last = ilc[k1];
						}
						++it;
					}
					else ++it;
					ret.reserve(aret.size());
					for (ElementArray<Edge>::iterator it = aret.begin(); it != aret.end(); ++it)
					{
						Element::adj_type adj = m->LowConn(*it);
						assert(m->Count(adj.data(), static_cast<enumerator>(adj.size()), hm) == 1);
						ret.push_back(adj[m->getNext(adj.data(), static_cast<enumerator>(adj.size()), ENUMUNDEF, hm)]);
					}
				}
				else
				{
					Element::adj_type& lc = m->LowConn(hc);
					for (Element::adj_type::size_type i = 0; i < lc.size(); i++) if (!m->GetMarker(lc[i], hm))  //iterate over faces
					{
						Element::adj_type& ilc = m->LowConn(lc[i]);
						for (Element::adj_type::size_type j = 0; j < ilc.size(); j++) if (!m->GetMarker(ilc[j], hm)) //iterate over face edges
						{
							Element::adj_type& jlc = m->LowConn(ilc[j]);
							for (Element::adj_type::size_type k = 0; k < jlc.size(); k++) if (!m->GetMarker(jlc[k], hm)) //iterator over edge nodes
							{
								if ((invert ^ m->GetPrivateMarker(jlc[k], mask)) && !m->GetPrivateMarker(jlc[k], mrk))
								{
									m->SetPrivateMarker(jlc[k], mrk);
									ret.push_back(jlc[k]);
								}
							}
						}
					}
				}
			}
		}
		else
		{
			if( !m->HideMarker() )
			{
				if (Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2) // This cell is 2d face
				{
					ElementArray<Edge> aret;
					adj_type const& lc = m->LowConn(GetHandle());
					aret.reserve(lc.size());
					HandleType q = lc[0];
					adj_type const& qlc = m->LowConn(q);
					HandleType first = qlc[0], last = qlc[1]; //edge 0
					if (invert ^ m->GetMarker(qlc[0], mask)) aret.push_back(qlc[0]); //node 0
					if (invert ^ m->GetMarker(qlc[1], mask)) aret.push_back(qlc[1]); //node 1
					HandleType r = lc[1]; //edge 1
					adj_type const& rlc = m->LowConn(r);
					if (first == rlc[0] || first == rlc[1])
					{
						last = first;
						if (aret.size() > 1)
						{
							HandleType temp = aret.data()[0];
							aret.data()[0] = aret.data()[1];
							aret.data()[1] = temp;
						}
					}
					adj_type::size_type it = 1, iend = lc.size() - 1;
					while (it < iend) //loop over edges
					{
						adj_type const& ilc = m->LowConn(lc[it]);
						if (last == ilc[0])
						{
							if (invert ^ m->GetMarker(ilc[1], mask))
								aret.push_back(ilc[1]);
							last = ilc[1];
						}
						else
						{
							if (invert ^ m->GetMarker(ilc[0], mask))
								aret.push_back(ilc[0]);
							last = ilc[0];
						}
						++it;
					}
					ret.reserve(aret.size());
					for (ElementArray<Edge>::iterator it = aret.begin(); it != aret.end(); ++it)
						ret.push_back(m->LowConn(*it)[0]);
				}
				else
				{
					Element::adj_type& lc = m->LowConn(hc);
					for (Element::adj_type::size_type i = 0; i < lc.size(); i++) //iterate over faces
					{
						Element::adj_type& ilc = m->LowConn(lc[i]);
						for (Element::adj_type::size_type j = 0; j < ilc.size(); j++) //iterate over face edges
						{
							Element::adj_type& jlc = m->LowConn(ilc[j]);
							for (Element::adj_type::size_type k = 0; k < jlc.size(); k++) //iterator over edge nodes
							{
								if ((invert ^ m->GetMarker(jlc[k], mask)) && !m->GetPrivateMarker(jlc[k], mrk))
								{
									m->SetPrivateMarker(jlc[k], mrk);
									ret.push_back(jlc[k]);
								}
							}
						}
					}
				}
			}
			else
			{
				MarkerType hm = m->HideMarker();
				if (Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2) // This cell is 2d face
				{
					ElementArray<Edge> aret;
					enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
					adj_type const& lc = m->LowConn(GetHandle());
					aret.reserve(lc.size());
					i = m->getNext(lc.data(), static_cast<enumerator>(lc.size()), i, hm);
					HandleType q = lc[i], first, last; //edge 0
					adj_type const& qlc = m->LowConn(q);
					k = m->getNext(qlc.data(), static_cast<enumerator>(qlc.size()), k, hm);
					if (invert ^ m->GetMarker(qlc[k], mask)) aret.push_back(qlc[k]); //node 0
					first = qlc[k];
					k = m->getNext(qlc.data(), static_cast<enumerator>(qlc.size()), k, hm);
					if (invert ^ m->GetMarker(qlc[k], mask)) aret.push_back(qlc[k]); //node 1
					last = qlc[k];
					i = m->getNext(lc.data(), static_cast<enumerator>(lc.size()), i, hm);
					HandleType r = lc[i]; //edge 1
					adj_type const& rlc = m->LowConn(r);
					k1 = m->getNext(rlc.data(), static_cast<enumerator>(rlc.size()), k1, hm);
					k2 = m->getNext(rlc.data(), static_cast<enumerator>(rlc.size()), k1, hm);
					if (first == rlc[k1] || first == rlc[k2])
					{
						last = first;
						if (aret.size() > 1)
						{
							HandleType temp = aret.data()[0];
							aret.data()[0] = aret.data()[1];
							aret.data()[1] = temp;
						}
					}
					adj_type::size_type it = 1, iend = lc.size() - 1;
					while (m->GetMarker(lc[iend], hm) && iend > 0) iend--;
					while (it < iend)
					{
						if (!m->GetMarker(lc[it], hm)) //loop over edges
						{
							adj_type const& ilc = m->LowConn(lc[it]);
							k1 = ENUMUNDEF;
							k1 = m->getNext(ilc.data(), static_cast<enumerator>(ilc.size()), k1, hm);
							k2 = m->getNext(ilc.data(), static_cast<enumerator>(ilc.size()), k1, hm);
							if (last == ilc[k1])
							{
								if (invert ^ m->GetMarker(ilc[k2], mask))
									aret.push_back(ilc[k2]);
								last = ilc[k2];
							}
							else
							{
								if (invert ^ m->GetMarker(ilc[k1], mask))
									aret.push_back(ilc[k1]);
								last = ilc[k1];
							}
						}
						++it;
					}
					ret.reserve(aret.size());
					for (ElementArray<Edge>::iterator it = aret.begin(); it != aret.end(); ++it)
					{
						Element::adj_type adj = m->LowConn(*it);
						assert(m->Count(adj.data(), static_cast<enumerator>(adj.size()), hm) == 1);
						ret.push_back(adj[m->getNext(adj.data(), static_cast<enumerator>(adj.size()), ENUMUNDEF, hm)]);
					}
				}
				else
				{
					Element::adj_type& lc = m->LowConn(hc);
					for (Element::adj_type::size_type i = 0; i < lc.size(); i++) if (!m->GetMarker(lc[i], hm))  //iterate over faces
					{
						Element::adj_type& ilc = m->LowConn(lc[i]);
						for (Element::adj_type::size_type j = 0; j < ilc.size(); j++) if (!m->GetMarker(ilc[j], hm)) //iterate over face edges
						{
							Element::adj_type& jlc = m->LowConn(ilc[j]);
							for (Element::adj_type::size_type k = 0; k < jlc.size(); k++) if (!m->GetMarker(jlc[k], hm)) //iterator over edge nodes
							{
								if ((invert ^ m->GetMarker(jlc[k], mask)) && !m->GetPrivateMarker(jlc[k], mrk))
								{
									m->SetPrivateMarker(jlc[k], mrk);
									ret.push_back(jlc[k]);
								}
							}
						}
					}
				}
			}
		}
		ret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return ret;
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
				if( lc.size() < 1 ) return aret;
				HandleType q = lc[0]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				if( qlc.size() > 0 ) aret.push_back(qlc[0]); //node 0
				if( qlc.size() > 1 ) aret.push_back(qlc[1]); //node 1
				if( lc.size() < 2 ) return aret;
				HandleType r = lc[1]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				if( aret.size() == 2 && rlc.size() == 2 )
				{
					if( aret.data()[0] == rlc[0] || aret.data()[0] == rlc[1] )
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
					if( aret.atback() == ilc[0] ) 
						aret.push_back(ilc[1]);
					else 
						aret.push_back(ilc[0]);
					++it;
				}
			}
			else
			{
				MarkerType mrk = m->CreatePrivateMarker();
				adj_type const & lc = m->LowConn(GetHandle());
				for(adj_type::size_type it = 0; it < lc.size(); it++) //faces
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) //edges
						if( !m->GetPrivateMarker(ilc[jt],mrk))
						{
							aret.push_back(ilc[jt]);
							m->SetPrivateMarker(ilc[jt],mrk);
						}
				}
				for(ElementArray<Edge>::size_type it = 0; it < aret.size(); it++)
					m->RemPrivateMarker(aret.at(it),mrk);
				m->ReleasePrivateMarker(mrk);
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
			{
				enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
				adj_type const & lc = m->LowConn(GetHandle());
				aret.reserve(lc.size());
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				if( i == lc.size() ) return aret;
				HandleType q = lc[i]; //edge 0
				adj_type const & qlc = m->LowConn(q);
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				if( k != qlc.size() ) aret.push_back(qlc[k]); //node 0
				k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
				if( k != qlc.size() ) aret.push_back(qlc[k]); //node 1
				i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
				if( i == lc.size() ) return aret;
				HandleType r = lc[i]; //edge 1
				adj_type const & rlc = m->LowConn(r);
				k1 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				k2 = m->getNext(rlc.data(),static_cast<enumerator>(rlc.size()),k1,hm);
				if( k1 != rlc.size() && k2 != rlc.size() && aret.size() == 2 )
				{
					if( aret.data()[0] == rlc[k1] || aret.data()[0] == rlc[k2] )
					{
						HandleType temp = aret.data()[0];
						aret.data()[0] = aret.data()[1];
						aret.data()[1] = temp;
					}
				}
				adj_type::size_type it = 1, iend = lc.size()-1;
				while(it < iend)
				{
					if( !m->GetMarker(lc[it],hm) ) //loop over edges
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						k1 = ENUMUNDEF; 
						k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
						k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
						if( aret.atback() == ilc[k1] ) 
							aret.push_back(ilc[k2]);
						else 
							aret.push_back(ilc[k1]);
					}
					++it;
				}
			}
			else
			{
				MarkerType mrk = m->CreatePrivateMarker();
				adj_type const & lc = m->LowConn(GetHandle());
				for(adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) ) //faces
				{
					adj_type const & ilc = m->LowConn(lc[it]);
					for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !m->GetMarker(ilc[jt],hm) )//edges
						if( !m->GetPrivateMarker(ilc[jt],mrk))
						{
							aret.push_back(ilc[jt]);
							m->SetPrivateMarker(ilc[jt],mrk);
						}
				}
				for(ElementArray<Edge>::size_type it = 0; it < aret.size(); it++)
					m->RemPrivateMarker(aret.at(it),mrk);
				m->ReleasePrivateMarker(mrk);
			}
		}
		return aret;
	}


	ElementArray<Edge> Cell::getEdges(MarkerType mask, bool invert) const
	{
		assert(GetHandleElementType(GetHandle())==CELL);
		Mesh * m = GetMeshLink();
		ElementArray<Edge> aret(m);
		if( isPrivate(mask) )
		{
			if( !m->HideMarker() )
			{
				if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
				{
					adj_type const & lc = m->LowConn(GetHandle());
					aret.reserve(lc.size());
					HandleType last, first;
					HandleType q = lc[0]; //edge 0
					adj_type const & qlc = m->LowConn(q);
					if( invert ^ m->GetPrivateMarker(qlc[0],mask) ) aret.push_back(qlc[0]); //node 0
					if( invert ^ m->GetPrivateMarker(qlc[1],mask) ) aret.push_back(qlc[1]); //node 1
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
						if( invert ^ m->GetPrivateMarker(last,mask) ) 
							aret.push_back(last);
						++it;
					}
				}
				else
				{
					MarkerType mrk = m->CreatePrivateMarker();
					adj_type const & lc = m->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++) //faces
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt != ilc.size(); jt++) //edges
							if( (invert ^ m->GetPrivateMarker(ilc[jt],mask)) && !m->GetPrivateMarker(ilc[jt],mrk))
							{
								aret.push_back(ilc[jt]);
								m->SetPrivateMarker(ilc[jt],mrk);
							}
					}
					for(ElementArray<Edge>::size_type it = 0; it != aret.size(); it++)
						m->RemPrivateMarker(aret.at(it),mrk);
					m->ReleasePrivateMarker(mrk);
				}
			}
			else
			{
				MarkerType hm = GetMeshLink()->HideMarker();
				if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
				{
					enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
					HandleType last, first;
					adj_type const & lc = m->LowConn(GetHandle());
					aret.reserve(lc.size());
					i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
					HandleType q = lc[i]; //edge 0
					adj_type const & qlc = m->LowConn(q);
					k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
					if( invert ^ m->GetPrivateMarker(qlc[k],mask) ) aret.push_back(qlc[k]); //node 0
					first = qlc[k];
					k = m->getNext(qlc.data(),static_cast<enumerator>(qlc.size()),k,hm);
					if( invert ^ m->GetPrivateMarker(qlc[k],mask) ) aret.push_back(qlc[k]); //node 1
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
					while(it != iend)
					{
						if( !m->GetMarker(lc[it],hm) ) //loop over edges
						{
							adj_type const & ilc = m->LowConn(lc[it]);
							k1 = ENUMUNDEF; 
							k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							if( last == ilc[k1] ) 
								last = ilc[k2];
							else last = ilc[k1];
							if( invert ^ m->GetPrivateMarker(last,mask) ) aret.push_back(last);
						}
						++it;
					}
				}
				else
				{
					MarkerType mrk = m->CreatePrivateMarker();
					adj_type const & lc = m->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) ) //faces
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !m->GetMarker(ilc[jt],hm) )//edges
							if( (invert ^ m->GetPrivateMarker(ilc[jt],mask)) && !m->GetPrivateMarker(ilc[jt],mrk))
							{
								aret.push_back(ilc[jt]);
								m->SetPrivateMarker(ilc[jt],mrk);
							}
					}
					for(ElementArray<Edge>::size_type it = 0; it < aret.size(); it++) 
						m->RemPrivateMarker(aret.at(it),mrk);
					m->ReleasePrivateMarker(mrk);
				}
			}
		}
		else
		{
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
					MarkerType mrk = m->CreatePrivateMarker();
					adj_type const & lc = m->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++) //faces
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt != ilc.size(); jt++) //edges
							if( (invert ^ m->GetMarker(ilc[jt],mask)) && !m->GetPrivateMarker(ilc[jt],mrk))
							{
								aret.push_back(ilc[jt]);
								m->SetPrivateMarker(ilc[jt],mrk);
							}
					}
					for(ElementArray<Edge>::size_type it = 0; it != aret.size(); it++)
						m->RemPrivateMarker(aret.at(it),mrk);
					m->ReleasePrivateMarker(mrk);
				}
			}
			else
			{
				MarkerType hm = GetMeshLink()->HideMarker();
				if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 ) // This cell is 2d face
				{
					enumerator i = ENUMUNDEF, k = ENUMUNDEF, k1 = ENUMUNDEF, k2;
					HandleType last, first;
					adj_type const & lc = m->LowConn(GetHandle());
					aret.reserve(lc.size());
					i = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),i,hm);
					HandleType q = lc[i]; //edge 0
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
					while(it != iend) 
					{
						if( !m->GetMarker(lc[it],hm) ) //loop over edges
						{
							adj_type const & ilc = m->LowConn(lc[it]);
							k1 = ENUMUNDEF; 
							k1 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							k2 = m->getNext(ilc.data(),static_cast<enumerator>(ilc.size()),k1,hm);
							if( last == ilc[k1] ) 
								last = ilc[k2];
							else last = ilc[k1];
							if( invert ^ m->GetMarker(last,mask) ) aret.push_back(last);
						}
						++it;
					}
				}
				else
				{
					MarkerType mrk = m->CreatePrivateMarker();
					adj_type const & lc = m->LowConn(GetHandle());
					for(adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) ) //faces
					{
						adj_type const & ilc = m->LowConn(lc[it]);
						for(adj_type::size_type jt = 0; jt < ilc.size(); jt++) if( !m->GetMarker(ilc[jt],hm) )//edges
							if( (invert ^ m->GetMarker(ilc[jt],mask)) && !m->GetPrivateMarker(ilc[jt],mrk))
							{
								aret.push_back(ilc[jt]);
								m->SetPrivateMarker(ilc[jt],mrk);
							}
					}
					for(ElementArray<Edge>::size_type it = 0; it < aret.size(); it++) 
						m->RemPrivateMarker(aret.at(it),mrk);
					m->ReleasePrivateMarker(mrk);
				}
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
			if( isPrivate(mask) )
			{
				for(adj_type::size_type it = 0; it < lc.size(); ++it)
					if( invert ^ m->GetPrivateMarker(lc[it],mask) ) aret.push_back(lc[it]);
			}
			else
			{
				for(adj_type::size_type it = 0; it < lc.size(); ++it)
					if( invert ^ m->GetMarker(lc[it],mask) ) aret.push_back(lc[it]);
			}
		}
		else
		{
			MarkerType hm = m->HideMarker();
			adj_type const & lc = m->LowConn(GetHandle());
			if( isPrivate(mask) )
			{
				for(adj_type::size_type it = 0; it < lc.size(); ++it)
					if( (invert ^ m->GetPrivateMarker(lc[it],mask)) && !m->GetMarker(lc[it],hm) ) 
						aret.push_back(lc[it]);
			}
			else
			{
				for(adj_type::size_type it = 0; it < lc.size(); ++it)
					if( (invert ^ m->GetMarker(lc[it],mask)) && !m->GetMarker(lc[it],hm) ) 
						aret.push_back(lc[it]);
			}
		}
		return aret;
	}
}
#endif
