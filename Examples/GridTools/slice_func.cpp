#include "inmost.h"
#include "slice_func.h"
using namespace INMOST;



Storage::real Slice::Search(Storage::real r0, Storage::real r1, Storage::real c0[3], Storage::real c1[3], Storage::real p[3], bool binary) const
{
	Storage::real rp = 1.0e20, rp_min = 1.0e20, p_min[3] = {0,0,0};
	int iters = 0;
	do
	{
		Storage::real m1 = r0/(r0-r1), m0;
		if( m1 < 0.0 || m1 > 1.0 || binary ) m1 = 0.5;
		m0 = 1.0-m1;
		p[0] = m1*c1[0] + m0*c0[0];
		p[1] = m1*c1[1] + m0*c0[1];
		p[2] = m1*c1[2] + m0*c0[2];
		rp = LevelFunction(p);
		if( rp_min > rp )
		{
			p_min[0] = p[0];
			p_min[1] = p[1];
			p_min[2] = p[2];
			rp_min = rp;
		}
		if( rp*r0 < 0.0 )
		{
			c1[0] = p[0];
			c1[1] = p[1];
			c1[2] = p[2];
			r1 = rp;
		}
		else if( rp*r1 < 0.0 )
		{
			c0[0] = p[0];
			c0[1] = p[1];
			c0[2] = p[2];
			r0 = rp;
		}
		iters++;
	} 
	while( fabs(rp) > epsf && iters < maxits );
	
	if( rp > rp_min )
	{
		p[0] = p_min[0];
		p[1] = p_min[1];
		p[2] = p_min[2];
		rp = rp_min;
	}
	
	return rp;
}


Storage::real Slice::SearchZero(Storage::real r0, Storage::real r1, Storage::real c0[3], Storage::real c1[3], Storage::real p[3]) const
{
	(void)r1;
	int iters = 0;
	Storage::real rp = 1.0e+20;
	//std::cout << "r0 " << r0 << " r1 " << r1;
	do
	{
		p[0] = 0.5*c1[0] + 0.5*c0[0];
		p[1] = 0.5*c1[1] + 0.5*c0[1];
		p[2] = 0.5*c1[2] + 0.5*c0[2];
		rp = LevelFunction(p);
		//std::cout << " rp " << rp;
		if( fabs(rp) < epsf )
		{
			if( fabs(r0) < epsf )
			{
				c0[0] = p[0];
				c0[1] = p[1];
				c0[2] = p[2];
				r0 = rp;
				//std::cout << "<(l)";
			}
			else
			{
				c1[0] = p[0];
				c1[1] = p[1];
				c1[2] = p[2];
				r1 = rp;
				//std::cout << "<(r)";
			}
		}
		else
		{
			if( fabs(r0) > epsf )
			{
				c0[0] = p[0];
				c0[1] = p[1];
				c0[2] = p[2];
				r0 = rp;
				//std::cout << ">(l)";
			}
			else
			{
				c1[0] = p[0];
				c1[1] = p[1];
				c1[2] = p[2];
				r1 = rp;
				//std::cout << ">(r)";
			}
		}
		iters++;
	} while( iters < maxits_zero );
	//std::cout << std::endl;
	
	return rp;
}


void Slice::SliceMesh(Mesh & m, bool remove_material_zero)
{
	//m.SetTopologyCheck(NEED_TEST_CLOSURE|PROHIBIT_MULTILINE|PROHIBIT_MULTIPOLYGON|GRID_CONFORMITY|DEGENERATE_EDGE|DEGENERATE_FACE|DEGENERATE_CELL | FACE_EDGES_ORDER | MARK_ON_ERROR | ADJACENT_DUPLICATE);
	//m->RemTopologyCheck(THROW_EXCEPTION);
	TagInteger material = m.CreateTag("MATERIAL",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag sliced = m.CreateTag("SLICED",DATA_BULK,FACE|EDGE|NODE,FACE|EDGE|NODE,1);
	TagReal level = m.CreateTag("level",DATA_REAL,NODE,NONE,1);

	
	MarkerType original = m.CreateMarker();
	std::cout << "marker original " << original << std::endl;
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int k = 0; k < m.NodeLastLocalID(); ++k) if(m.isValidNode(k))
	{
		Node it = m.NodeByLocalID(k);
		it->SetMarker(original);
		material[it] = 3;
		level[it] = LevelFunction(it->Coords().data());
	}
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it) it->SetMarker(original);

	MarkerType slice = m.CreateMarker();
	std::cout << "marker slice " << slice << std::endl;
	MarkerType mrk = m.CreateMarker();
	std::cout << "marker mrk " << mrk << std::endl;
	int nslice = 0, nmark = 0;
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it) if( !it->GetMarker(mrk) )
	{
		Storage::real p[3],pc0[3],pc1[3];
		Storage::real_array c0 = it->getBeg()->Coords();
		Storage::real_array c1 = it->getEnd()->Coords();
		Storage::real r0 = LevelFunction(c0.data());
		Storage::real r1 = LevelFunction(c1.data());
		int m0 = material[it->getBeg()];
		int m1 = material[it->getEnd()];
		//std::cout << "m0 " << m0 << " m1 " << m1;
		bool cut = true;
		if (m0 == 3)
		{
			if( fabs(r0) < epsf )
			{
				m0 = 2;
				it->getBeg()->SetMarker(slice);
				cut = false;
				nmark++;
			}
			else if ( r0 < -epsf )
				m0 = 0;
			else
				m0 = 1;
		}
		if (m1 == 3)
		{
			if( fabs(r1) < epsf )
			{
				m1 = 2;
				it->getEnd()->SetMarker(slice);
				cut = false;
				nmark++;
			}
			else if( r1 < -epsf )
				m1 = 0;
			else
				m1 = 1;
		}
		material[it->getBeg()] = m0;
		material[it->getEnd()] = m1;
		//std::cout << " e " << it->LocalID() << " r0 " << r0 << " r1 " << r1 << " m0 " << m0 << " m1 " << m1;
		if( (cut || ((fabs(r0) < epsf) ^ (fabs(r1) < epsf))) && m0 != m1 )
		//if( (r0*r1 < -1.0e-12) || (fabs(r0*r1) < 1.0e-12 && ((fabs(r0) < 1.0e-6) ^ (fabs(r1) < 1.0e-6))) )
		{
			int alg = -1;
			(void)alg;
			pc0[0] = c0[0], pc0[1] = c0[1], pc0[2] = c0[2];
			pc1[0] = c1[0], pc1[1] = c1[1], pc1[2] = c1[2];
			if((fabs(r0) < epsf) ^ (fabs(r1) < epsf))
			{
				SearchZero(r0,r1,pc0,pc1,p);
				alg = 1;
			}
			else
			{
				Storage::real rp = Search(r0,r1,pc0,pc1,p);
				alg = 2;
				if( fabs(rp) > epsf ) //cannot find intersection
				{
					rp = Search(r0,r1,pc0,pc1,p,true);
					//p[0] = c0[0]*0.5+c1[0]*0.5;
					//p[1] = c0[1]*0.5+c1[1]*0.5;
					//p[2] = c0[2]*0.5+c1[2]*0.5;
					alg = 3;
				}
			}
			//p[0] = (r0*c1[0] - r1*c0[0])/(r0-r1);
			//p[1] = (r0*c1[1] - r1*c0[1])/(r0-r1);
			//p[2] = (r0*c1[2] - r1*c0[2])/(r0-r1);
			//pstd::cout << " p " << p[0] << " " << p[1] << " " << p[2];
			//distance to the corners
			
			Storage::real l0 = 0, l1 = 0, l;
			for(int r = 0; r < 3; ++r)
			{
				l0 += (p[r]-c0[r])*(p[r]-c0[r]);
				l1 += (p[r]-c1[r])*(p[r]-c1[r]);
			}
			l0 = sqrt(l0);
			l1 = sqrt(l1);
			l = l0+l1;
			//std::cout << " l0 " << l0 << " l1 " << l1;
			if( l0 < epsl*l )
			{
				material[it->getBeg()] = 2;
				it->getBeg()->SetMarker(slice);
				nmark++;
				//std::cout << " left " << func(c0[0],c0[1],c0[2],type) << " alg " << alg;
			}
			else if( l1 < epsl*l )
			{
				material[it->getEnd()] = 2;
				it->getEnd()->SetMarker(slice);
				nmark++;
				//std::cout << " right " << func(c1[0],c1[1],c1[2],type) << " alg " << alg;
			}
			else
			{
				//std::cout << " new " << func(p[0],p[1],p[2],type) << " alg " << alg;
				Node n = m.CreateNode(p);
				level[n] = LevelFunction(p);
				material[n] = 2;
				n.Bulk(sliced) = 1;
				n.SetMarker(slice);
				bool was_sliced = it->HaveData(sliced) ? true : false;
				ElementArray<Edge> ret = Edge::SplitEdge(it->self(),ElementArray<Node>(&m,1,n.GetHandle()),0);
				ret.SetMarker(mrk);
				if( was_sliced ) for(size_t q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
				nslice++;
			}
		}
		//std::cout << std::endl;
	}
	
	

	std::cout << "sliced edges: " << nslice << " marked nodes: " << nmark << std::endl;

	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;
	
	{
		int mat[3] = {0,0,0};
		int tot = 0;
		for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
		{
			mat[material[*it]]++;
			tot++;
		}
		std::cout << "node materials, 0: " << mat[0] << " 1: " << mat[1] << " 2: " << mat[2] << std::endl;
	}
	
	
	
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
	{
		int mat[3] = {0,0,0};
		mat[material[it->getBeg()]]++;
		mat[material[it->getEnd()]]++;
		if( !(mat[0] == 0 || mat[1] == 0) )
		{
			material[*it] = 2;
			std::cout << "oops, materials for edge nodes were not split, 0: " << mat[0] << ", 1: " << mat[1] << ", 2: " << mat[2] << std::endl;
                        
		}
		else if( mat[0] != 0 ) material[*it] = 0;
		else if( mat[1] != 0 ) material[*it] = 1;
		else material[*it] = 2;
	}
	
	
	nslice = 0;
	MarkerType unique = m.CreateMarker();
	std::cout << "marker unique " << unique << std::endl;
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) if( !it->GetMarker(mrk) )
	{
		ElementArray<Node> nodes = it->getNodes(slice); //those nodes should be ordered so that each pair forms an edge
		if( nodes.size() > 1 ) // if there is 1, then only one vertex touches the plane
		{
			//int nsliced = 0;
			//for(int q = 0; q < (int)nodes.size(); ++q) if( !nodes[q].GetMarker(original) ) nsliced++;
			//if there is more then two, then original face is non-convex
			//if( nsliced && nodes.size() > 2 ) std::cout << "Looks like face " << it->LocalID() << " is nonconvex, there is " << it->nbAdjElements(NODE) << " nodes, out of them " << nsliced << " new cuts on face" << " slice " << slice << " original " << original << std::endl;
			//else
			if( nodes.size() == 2 )
			{
				Edge e = m.CreateEdge(nodes).first;
				material[e] = 2; //on slice
				e.Bulk(sliced) = 1;
				e.SetMarker(slice);
				bool was_sliced = it->HaveData(sliced) ? true : false;
				ElementArray<Face> ret = Face::SplitFace(it->self(),ElementArray<Edge>(&m,1,e.GetHandle()),0);
				ret.SetMarker(mrk);
				if( was_sliced ) for(size_t q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
				nslice++;
			}
			else if( it->nbAdjElements(NODE,slice) != it->nbAdjElements(NODE) ) //not entire face is sliced
			{
				//std::cout << "Face " << it->LocalID();
				//std::cout << " sliced nodes " << nodes.size() << "/";
				nodes = it->getNodes();
				//std::cout << nodes.size();// << std::endl;

				
				Storage::real c0[3],c1[3],pc0[3],pc1[3],p[3];
				it->Centroid(c0);
				Storage::real r0 = LevelFunction(c0);
				int m0,m1;
				if( fabs(r0) < epsf ) m0 = 2;
				else if( r0 < -epsf ) m0 = 0;
				else m0 = 1;
				
				//std::cout << " centernode m " << m0 << " r " << r0;
				Node centernode = InvalidNode();
				ElementArray<Edge> split_edges(&m);
				ElementArray<Node> cutnodes(&m,nodes.size()), edge_nodes(&m,2);
				
				//calculate nodes that cut along the edges connecting centernode
				for(int q = 0; q < (int)nodes.size(); ++q) if( !nodes[q].GetMarker(slice) )
				{
					nodes[q].Centroid(c1);
					Storage::real r1 = LevelFunction(c1);
					m1 = material[nodes[q]];
					bool cut = !(m0 == 2 || m1 == 2);
					//std::cout << " n" << q << " m " << m1 << " r1 " << r1;
					//if( (r0*r1 < -1.0e-12) || (fabs(r0*r1) < 1.0e-12 && ((fabs(r0) < 1.0e-6) ^ (fabs(r1) < 1.0e-6))))
					if( (cut || ((fabs(r0) < epsf) ^ (fabs(r1) < epsf))) && m0 != m1 )
					{
						pc0[0] = c0[0], pc0[1] = c0[1], pc0[2] = c0[2];
						pc1[0] = c1[0], pc1[1] = c1[1], pc1[2] = c1[2];
						if( ((fabs(r0) < epsf) ^ (fabs(r1) < epsf)) )
							SearchZero(r0,r1,pc0,pc1,p);
						else
						{						
							Storage::real rp = Search(r0,r1,pc0,pc1,p);
							if( fabs(rp) > epsf )
							{
								//std::cout << "inaccurate search " << rp;
								rp = Search(r0,r1,pc0,pc1,p,true);
								//std::cout << " binary " << rp << std::endl;
								//p[0] = c0[0]*0.5+c1[0]*0.5;
								//p[1] = c0[1]*0.5+c1[1]*0.5;
								//p[2] = c0[2]*0.5+c1[2]*0.5;
								
							}
						}
						//distance to the center node
						Storage::real l0 = 0, l1 = 0, l;
						for(int r = 0; r < 3; ++r)
						{
							l0 += (p[r]-c0[r])*(p[r]-c0[r]);
							l1 += (p[r]-c1[r])*(p[r]-c1[r]);
						}
						l0 = sqrt(l0);
						l1 = sqrt(l1);
						l = l0+l1;
						if( l0 < epsl*l ) //edge goes through centernode
						{
							if( !centernode.isValid() )
							{
								centernode = m.CreateNode(c0);
								level[centernode] = r0;
							}
							cutnodes[q] = centernode;
							//std::cout << " (c) ";
						}
						else if( l1 > epsl*l )
						{
							cutnodes[q] = m.CreateNode(p);
							level[cutnodes[q]] = LevelFunction(p);
							//std::cout << " (n) ";
						}
						else
						{
							//std::cout << " (o) ";
							material[nodes[q]] = 2;
						}
					}
					else if( m0 == 2 )
					{
						if( !centernode.isValid() )
						{
							centernode = m.CreateNode(c0);
							level[centernode] = r0;
						}
						cutnodes[q] = centernode;
						//std::cout << " (C) ";
					}
				}
				//std::cout << std::endl;
				
				for(int q = 0; q < (int)cutnodes.size(); ++q) if( cutnodes[q].isValid() )
				{
					//std::cout << "New cut node " << cutnodes[q].LocalID() << " at " << q << " on line with " << nodes[q].LocalID() << " r=" << func(cutnodes[q].Coords()[0],cutnodes[q].Coords()[1],cutnodes[q].Coords()[2],type) << std::endl;
					Node n = cutnodes[q];
					material[n] = 2;
					n.Bulk(sliced) = 1;
					n.SetMarker(slice);
				}
				
				
				//go over triangles and figure out the path of the cut
				for(int q = 0; q < (int)nodes.size(); ++q)
				{
					int i1 = q;
					int i2 = (q+1)%nodes.size();
					Node n1 = nodes[i1];
					Node n2 = nodes[i2];
					bool s1 = n1->GetMarker(slice);
					bool s2 = n2->GetMarker(slice);
					if( s1 && s2 )
					{
						//cut passing through edge
						//skip this
					}
					else if( s2 )
					{
						//check cut on oposite edge
						if( cutnodes[i1].isValid() )
						{
							edge_nodes[0] = n2;
							edge_nodes[1] = cutnodes[i1];
							std::pair<Edge,bool> en = m.CreateEdge(edge_nodes);
							//if( en.second )
							{
								Edge e = en.first;
								material[e] = 2; //on slice
								e.Bulk(sliced) = 1;
								e.SetMarker(slice);
								split_edges.push_back(e);
							}
						}
					}
					else if( s1 )
					{
						//check cut on oposite edge
						if( cutnodes[i2].isValid() )
						{
							edge_nodes[0] = n1;
							edge_nodes[1] = cutnodes[i2];
							std::pair<Edge,bool> en = m.CreateEdge(edge_nodes);
							//if( en.second )
							{
								Edge e = en.first;
								material[e] = 2; //on slice
								e.Bulk(sliced) = 1;
								e.SetMarker(slice);
								split_edges.push_back(e);
							}
						}
					}
					else if( cutnodes[i1].isValid() && cutnodes[i2].isValid() && cutnodes[i1] != cutnodes[i2])
					{
						edge_nodes[0] = cutnodes[i1];
						edge_nodes[1] = cutnodes[i2];
						std::pair<Edge,bool> en = m.CreateEdge(edge_nodes);
						//if( en.second )
						{
							Edge e = en.first;
							material[e] = 2; //on slice
							e.Bulk(sliced) = 1;
							e.SetMarker(slice);
							split_edges.push_back(e);
						}
					}
				}
				//split face with multiple edges
				if( !split_edges.empty() )
				{
					
					int k = 0;
					for(size_t q = 0; q < split_edges.size(); ++q)
					{
						if( !split_edges[q].GetMarker(unique) )
						{
							split_edges[q].SetMarker(unique);
							split_edges[k++] = split_edges[q];
						}
					}
					split_edges.RemMarker(unique);
					split_edges.resize(k);
					//std::cout << "Split edges " << split_edges.size() << std::endl;
					//for(int q = 0; q < (int)split_edges.size(); ++q)
					//	std::cout << split_edges[q].getBeg().LocalID() << "<->" << split_edges[q].getEnd().LocalID() << std::endl;
					bool was_sliced = it->HaveData(sliced) ? true : false;
					//std::cout << "Split face " << it->LocalID();
					ElementArray<Face> ret = Face::SplitFace(it->self(),split_edges,0);
					ret.SetMarker(mrk);
					//std::cout << " new faces: " << ret.size() << ":";
					for(ElementArray<Face>::size_type q = 0; q < ret.size(); ++q)
					{
						int mat[3] = {0,0,0};
						ElementArray<Edge> fe = ret[q].getEdges();
						for( ElementArray<Edge>::size_type l = 0; l < fe.size(); ++l) mat[material[fe[l]]]++;
						//std::cout << " " << ret[q].LocalID() << "[" << mat[0] << "," << mat[1] << "," << mat[2] << "]";
					}
					//std::cout << std::endl;
					if( was_sliced ) for(INMOST_DATA_ENUM_TYPE q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
					nslice++;
				}
				//else std::cout << "No split edges " << std::endl;
			}
			//else std::cout << "No new edges on face " << it->LocalID() << std::endl;
		}
		//else std::cout << "Only one adjacent slice node, face " << it->LocalID() << std::endl;
	}
	m.ReleaseMarker(unique);
	std::cout << "release marker unique " << unique << std::endl;
	

	nmark = 0;
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
		if( !it->GetMarker(slice) && it->getBeg()->GetMarker(slice) && it->getEnd()->GetMarker(slice) )
		{
			if( material[*it] != 2 )
				std::cout << "Edge " << it->LocalID() << " supposed to get material 2, but have " << material[*it] << " nodes " << material[it->getBeg()] << " " << material[it->getEnd()] << std::endl;
			it->SetMarker(slice);
			nmark++;
		}

	std::cout << "sliced faces: " << nslice << " marked edges: " << nmark << std::endl;

	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;
	
	{
		int mat[3] = {0,0,0};
		int tot = 0;
		for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
		{
			mat[material[*it]]++;
			tot++;
		}
		std::cout << "edge materials, 0: " << mat[0] << " 1: " << mat[1] << " 2: " << mat[2] << std::endl;
	}
	
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
	{
		int mat[3] = {0,0,0};
		ElementArray<Edge> edges = it->getEdges();
		for(int k = 0; k < (int)edges.size(); ++k)
			mat[material[edges[k]]]++;
		if( !(mat[0] == 0 || mat[1] == 0) )
		{
			material[*it] = 2;
			std::cout << "oops, materials for face " << it->LocalID() << " edges were not split, 0: " << mat[0] << " ,1: " << mat[1] << " ,2: " << mat[2] <<
			" slice edges " << it->getEdges(slice).size() << "/" << it->getEdges().size() << " nodes " << it->getNodes(slice).size() <<  "/" << it->getNodes().size() << " mrk " << (it->GetMarker(mrk) ? "yes":"no") << std::endl;
		}
		else if( mat[0] != 0 ) material[*it] = 0;
		else if( mat[1] != 0 ) material[*it] = 1;
		else material[*it] = 2;
	}
	
	nslice = 0;
	TagInteger indx = m.CreateTag("TEMP_INDX",DATA_INTEGER,NODE|EDGE,NONE,1);
	MarkerType visited = m.CreateMarker();
	std::cout << "marker visited " << visited << std::endl;
	MarkerType cmrk = m.CreateMarker();
	std::cout << "marker cmrk " << cmrk << std::endl;
	MarkerType isolate = m.CreateMarker();
	std::cout << "marker isolate " << isolate << std::endl;
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( !it->GetMarker(mrk) )
	{
		ElementArray<Edge> edges = it->getEdges(slice);
		if( !edges.empty() ) //these should form a triangle
		{
			//check edges form a simple loops, each node should be visited twice
			ElementArray<Face> split_faces(&m);
			std::map<Node,int> visit_count;
			for(ElementArray<Edge>::iterator jt = edges.begin(); jt != edges.end(); ++jt)
			{
				visit_count[jt->getBeg()]++;
				visit_count[jt->getEnd()]++;
			}
			bool simple = true;
			for(std::map<Node,int>::iterator jt = visit_count.begin(); jt != visit_count.end(); ++jt)
				simple &= (jt->second == 2);
			
			//std::cout << "Cell " << it->LocalID() << " slice edges " << edges.size() << "/" << it->nbAdjElements(EDGE) << " " << (simple?"simple":"complex") << std::endl;
			
			if( simple )
			{
				//gather loop by loop and create faces
				//int loop_cnt = 0;
				//order edges
				ElementArray<Edge> order_edges(&m);
				
				Node last;
				ElementArray<Edge>::size_type nvisited = 0;
				bool found = false;
				for(ElementArray<Edge>::size_type k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(original) )
				{
					nvisited++;
					last = InvalidNode();
					order_edges.push_back(edges[k]);
					order_edges.SetMarker(visited);
					found = true;
					break;
				}
				
				if( !found ) continue;
				
				
				while(nvisited != edges.size() )
				{
					found = false;
					for(ElementArray<Edge>::size_type k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(visited) )
					{
						bool match = false;
						if( last.isValid() )
						{
							if( edges[k]->getBeg() == last )
							{
								match = true;
								last = edges[k]->getEnd();
							}
							else if( edges[k]->getEnd() == last )
							{
								match = true;
								last = edges[k]->getBeg();
							}
						}
						else
						{
							if( edges[k]->getBeg() == order_edges.back()->getBeg() || edges[k]->getBeg() == order_edges.back()->getEnd() )
							{
								match = true;
								last = edges[k]->getEnd();
							}
							else if( edges[k]->getEnd() == order_edges.back()->getBeg() || edges[k]->getEnd() == order_edges.back()->getEnd() )
							{
								match = true;
								last = edges[k]->getBeg();
							}
						}
						if( match )
						{
							nvisited++;
							order_edges.push_back(edges[k]);
							order_edges.back().SetMarker(visited);
							found = true;
						}
					}
					if( !found || nvisited == edges.size() )
					{
						if( !order_edges.empty() )
						{
							//std::cout << "New loop " << ++loop_cnt << ": " << order_edges.size() << " total edges " << edges.size() << std::endl;
							//for(int k = 0; k < order_edges.size(); ++k)
							//	std::cout << order_edges[k].getBeg().LocalID() << "<->" << order_edges[k].getEnd().LocalID() << std::endl;
							
							//std::cout << "All:" << std::endl;
							//for(int k = 0; k < edges.size(); ++k)
							//	std::cout << edges[k].getBeg().LocalID() << "<->" << edges[k].getEnd().LocalID() << " " << (edges[k].GetMarker(visited) ? "visited":"") << " " << (edges[k].GetMarker(original) ? "original":"")<< std::endl;
							if( order_edges.size() > 2 )
							{
								std::pair<Face,bool> f = m.CreateFace(order_edges);
								material[f.first] = 2;
								f.first->Bulk(sliced) = 1;
								split_faces.push_back(f.first);
							}
							order_edges.clear();
							found = false;
							for(ElementArray<Edge>::size_type k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(visited) && !edges[k]->GetMarker(original))
							{
								nvisited++;
								last = InvalidNode();
								order_edges.push_back(edges[k]);
								order_edges.back().SetMarker(visited);
								found = true;
								break;
							}
							if( !found ) break;
						}
						else break;
					}
				}
				
				//if( !order_edges.empty() ) std::cout << "size: " << order_edges.size() << std::endl;
				
				edges.RemMarker(visited);
				
			}
			else
			{
				//edges do not breakup into simple loops, have to run complicated algorithm
				//split into pyramids joining faces
				ElementArray<Face> cfaces = it->getFaces();
				ElementArray<Edge> cedges = it->getEdges();
				ElementArray<Node> cnodes = it->getNodes();
				ElementArray<Node> cutnodes(&m,cnodes.size());
				ElementArray<Edge> cutedges(&m,cedges.size());
				ElementArray<Node> edge_nodes(&m,2);
				
				for(Storage::integer k = 0; k < (Storage::integer)cnodes.size(); ++k) indx[cnodes[k]] = k;
				for(Storage::integer k = 0; k < (Storage::integer)cedges.size(); ++k) indx[cedges[k]] = k;
				
				Storage::real c0[3],c1[3],pc0[3],pc1[3],p[3];
				it->Centroid(c0);
				Storage::real r0 = LevelFunction(c0);

				int m0,m1;
				if( fabs(r0) < epsf ) m0 = 2;
				else if( r0 < -epsf ) m0 = 0;
				else m0 = 1;
				
				
				//std::cout << "Number of cut edges: " << edges.size() << std::endl;
				//for(std::map<Node,int>::iterator jt = visit_count.begin(); jt != visit_count.end(); ++jt)
				//std::cout << "NODE:" << jt->first.LocalID() << " visited " << jt->second << " times" << std::endl;
				//std::cout << "Centernode m " << m0 << " r " << r0 << std::endl;
				
				Node centernode = InvalidNode();
				
				
				//calculate nodes that cut along the edges connecting centernode
				for(ElementArray<Node>::size_type q = 0; q < cnodes.size(); ++q)
				{
					if( !cnodes[q].GetMarker(slice) )
					{
						cnodes[q].Centroid(c1);
						Storage::real r1 = LevelFunction(c1);
						m1 = material[cnodes[q]];
						bool cut = !(m0 == 2 || m1 == 2);
						//std::cout << "NODE:" << cnodes[q].LocalID() << " m " << m1 << " r " << r1  << " " << (cnodes[q].GetMarker(slice)?"":"not ") << "sliced" << std::endl;
						if( (cut || ((fabs(r0) < epsf) ^ (fabs(r1) < epsf))) && m0 != m1 )
						//if( (r0*r1 < -1.0e-12) || (fabs(r0*r1) < 1.0e-12 && ((fabs(r0) < 1.0e-6) ^ (fabs(r1) < 1.0e-6))))
						{
							pc0[0] = c0[0], pc0[1] = c0[1], pc0[2] = c0[2];
							pc1[0] = c1[0], pc1[1] = c1[1], pc1[2] = c1[2];
							if((fabs(r0) < epsf) ^ (fabs(r1) < epsf))
								SearchZero(r0,r1,pc0,pc1,p);
							else
							{						
								Storage::real rp = Search(r0,r1,pc0,pc1,p);
								if( fabs(rp) > epsf )
								{
									//std::cout << "inaccurate search " << rp;
									rp = Search(r0,r1,pc0,pc1,p,true);
									//std::cout << " binary " << rp << std::endl;
									//p[0] = c0[0]*0.5+c1[0]*0.5;
									//p[1] = c0[1]*0.5+c1[1]*0.5;
									//p[2] = c0[2]*0.5+c1[2]*0.5;
									
								}
							}
							//distance to the center node
							Storage::real l0 = 0, l1 = 0, l;
							for(int r = 0; r < 3; ++r)
							{
								l0 += (p[r]-c0[r])*(p[r]-c0[r]);
								l1 += (p[r]-c1[r])*(p[r]-c1[r]);
							}
							l0 = sqrt(l0);
							l1 = sqrt(l1);
							l = l0+l1;
							//std::cout << "l0 " << l0 << " l1 " << l1 << std::endl;
							if( l0 < epsl*l ) //edge goes through centernode
							{
								if( !centernode.isValid() )
								{
									centernode = m.CreateNode(c0);
									level[centernode] = r0;
								}
								cutnodes[q] = centernode;
								//std::cout << "selected centernode " << cutnodes[q].LocalID() << std::endl;
							}
							else if( l1 > epsl*l )
							{
								cutnodes[q] = m.CreateNode(p);
								level[cutnodes[q]] = LevelFunction(p);
								//std::cout << "created new node " << cutnodes[q].LocalID() << std::endl;
							}
							else if( !cnodes[q].GetMarker(slice))
							{
								material[cnodes[q]] = 2;
								cnodes[q].SetMarker(slice);
								ElementArray<Edge> nedges = cnodes[q]->getEdges();
								for(INMOST_DATA_ENUM_TYPE r = 0; r < nedges.size();++r)
									if( material[nedges[r]->getBeg()] == 2 && material[nedges[r]->getEnd()] == 2)
										material[nedges[r]] = 2;
								//std::cout << "use old node " << std::endl;
							}
						}
						else if( m0 == 2 )
						{
							if( !centernode.isValid() )
							{
								centernode = m.CreateNode(c0);
								level[centernode] = r0;
							}
							cutnodes[q] = centernode;
							//std::cout << "selected centernode(2) " << cutnodes[q].LocalID() << std::endl;
						}
					}
					//else cutnodes[q] = cnodes[q];
				}
				
				for(int q = 0; q < (int)cutnodes.size(); ++q) if( cutnodes[q].isValid() )
				{
					//std::cout << "New cut node " << cutnodes[q].LocalID() << " at " << q << " on line with " << cnodes[q].LocalID() << " r=" << LevelFunction(cutnodes[q].Coords().data()) << std::endl;
					Node n = cutnodes[q];
					material[n] = 2;
					n.Bulk(sliced) = 1;
					n.SetMarker(slice);
				}
				
				for(int q = 0; q < (int)cnodes.size(); ++q) if( cnodes[q].GetMarker(slice) ) cutnodes[q] = cnodes[q];
				
				//now find cutedges
				
				for(int q = 0; q < (int)cedges.size(); ++q )
				{
					if( !cedges[q].GetMarker(slice ) )
					{
						Node n1 = cutnodes[indx[cedges[q].getBeg()]];
						Node n2 = cutnodes[indx[cedges[q].getEnd()]];
						if( n1.isValid() && n2.isValid() && n1 != n2 )
						{
							edge_nodes[0] = n1;
							edge_nodes[1] = n2;
							cutedges[q] = m.CreateEdge(edge_nodes).first;
						}
					}
				}
				
				ElementArray<Edge> alledges(&m); //to draw
				
				for(int q = 0; q < (int)cedges.size(); ++q )
				{
					//std::cout << "Cell edge " << cedges[q].LocalID() << " m " << material[cedges[q]]  << " " << (cedges[q].GetMarker(slice) ? "slice":"orig") << " at " << q << " " << cedges[q].getBeg().LocalID() << "<->" << cedges[q].getEnd().LocalID() << " m " << material[cedges[q].getBeg()] << "<->" << material[cedges[q].getEnd()] << std::endl;
					alledges.push_back(cedges[q]);
				}
				
				for(int q = 0; q < (int)cutedges.size(); ++q) if( cutedges[q].isValid() )
				{
					//std::cout << "New cut edge " << cutedges[q].LocalID() << " at " << q << " on plane with " << cedges[q].LocalID() << " " << cutedges[q].getBeg().LocalID() << "<->" << cutedges[q].getEnd().LocalID() << " original " << cedges[q].getBeg().LocalID() << "<->" << cedges[q].getEnd().LocalID() << std::endl;
					Edge n = cutedges[q];
					material[n] = 2;
					n.Bulk(sliced) = 1;
					n.SetMarker(slice);
					alledges.push_back(n);
				}
				
				for(int q = 0; q < (int)cedges.size(); ++q) if( cedges[q].GetMarker(slice) ) cutedges[q] = cedges[q];
				
				ElementArray<Edge> split_edges(&m);
				//run over pyramids and collect faces, they should be already ordered into loops
				//although we still have to check they form closed loop
				std::map<Edge,int> vstcnt;
				for(int q = 0; q < (int)cfaces.size(); ++q)
				{
					split_edges.clear();
					ElementArray<Edge> fedges = cfaces[q].getEdges();
					for(int r = 0; r < (int)fedges.size(); ++r)
					{
						if( cutedges[indx[fedges[r]]].isValid() )
							split_edges.push_back(cutedges[indx[fedges[r]]]);
					}
					
					if( split_edges.size() > 2 ) //at least a triangle
					{
						visit_count.clear();
						for(int r = 0; r < (int)split_edges.size(); ++r)
						{
							visit_count[split_edges[r].getBeg()]++;
							visit_count[split_edges[r].getEnd()]++;
						}
						bool simple = true;
						for(std::map<Node,int>::iterator jt = visit_count.begin(); jt != visit_count.end(); ++jt)
							simple &= (jt->second == 2);
						
						if( simple )
						{
							for(INMOST_DATA_ENUM_TYPE r = 0; r < split_edges.size(); ++r) vstcnt[split_edges[r]]++;
							std::pair<Face,bool> f = m.CreateFace(split_edges);
							//std::cout << "Created face " << f.first.LocalID() << " with " << split_edges.size() << " edges: ";
							//for(int r = 0; r < split_edges.size(); ++r) std::cout << "EDGE:" << split_edges[r].LocalID()  << " ";
							//std::cout << " original FACE:" << cfaces[q].LocalID() << std::endl;
							material[f.first] = 2;
							f.first.Bulk(sliced) = 1;
							if( f.second ) split_faces.push_back(f.first);
							//else std::cout << "already existed" << std::endl;
						}
						else
						{
							//std::cout << "Face not created, node visits:";
							//for(std::map<Node,int>::iterator jt = visit_count.begin(); jt != visit_count.end(); ++jt)
							//	std::cout << " NODE:" << jt->first.LocalID() << " " << jt->second;
							//std::cout << " edges: ";
							//for(int r = 0; r < split_edges.size(); ++r) std::cout << "EDGE:" << split_edges[r].LocalID()  << " ";
							//std::cout << " original FACE:" << cfaces[q].LocalID() << std::endl;
						}
					}
					else if( !split_edges.empty() )
					{
						//std::cout << "Face not created, edges:";
						//for(int r = 0; r < (int)split_edges.size(); ++r)
						//	std::cout << " EDGE:" << split_edges[r].LocalID();
						//std::cout << " original FACE:" << cfaces[q].LocalID() << std::endl;
					}
				}
				
				//resulting surface may touch the element with one or more edge
				//making it impossible to separate in conformal way
				//
				// in this case we have to isolate this edge with additional faces
				//
				// check if the edge is counted twice and appears on original element
				visit_count.clear();
				cedges.SetMarker(cmrk);
				
				bool isolate_alogirthm = false;
				for(std::map<Edge,int>::iterator jt = vstcnt.begin(); jt != vstcnt.end(); ++jt)
					if( jt->first.GetMarker(cmrk) && jt->second > 1 )
					{
						//std::cout << "edge " << jt->first.LocalID() << " m " << material[jt->first] << " " << jt->first.getBeg().LocalID() << "<->" << jt->first.getEnd().LocalID() << " material " << material[jt->first.getBeg()] << "<->" << material[jt->first.getEnd()] << " visit count " << jt->second << std::endl;
						isolate_alogirthm = true;
						jt->first.SetMarker(isolate);
					}
				
				if( isolate_alogirthm )
				//if( false )
				{
					//std::cout << "Isolation algorithm" << std::endl;
					//on pyramids, build faces that start from cutedges and end up with edge on original edge
					for(int q = 0; q < (int)cfaces.size(); ++q)
					{
						split_edges.clear();
						
						ElementArray<Edge> fedges = cfaces[q].getEdges();
						bool have_isolate = false;
						for(int r = 0; r < (int)fedges.size(); ++r) have_isolate |= fedges[r].GetMarker(isolate);
						if( have_isolate )
						{
							for(int r = 0; r < (int)fedges.size(); ++r) if( !fedges[r].GetMarker(isolate) && !fedges[r].GetMarker(slice) && cutedges[indx[fedges[r]]].isValid() )
							{
								split_edges.clear();
								int cnt = 0;
								Edge e;
								Face f;
								int mat = material[fedges[r]];
								e = fedges[r];
								Node e1n1 = e.getBeg();
								Node e1n2 = e.getEnd();
								Node e2n1 = cutnodes[indx[e.getBeg()]];
								Node e2n2 = cutnodes[indx[e.getEnd()]];
							
								//std::cout << "Cell edge: " << e.LocalID() << " " << e.getBeg().LocalID() << "<->" << e.getEnd().LocalID() << std::endl;
								//std::cout << "Cut edge: " << e.LocalID() << " " << cutnodes[indx[e.getBeg()]].LocalID() << "<->" << cutnodes[indx[e.getEnd()]].LocalID() << std::endl;
								split_edges.push_back(e);
								
								if( e2n1.isValid() && e2n1 != e1n1 )
								{
									edge_nodes[0] = e1n1;
									edge_nodes[1] = e2n1;
									e = m.CreateEdge(edge_nodes).first;
									material[e] = mat;
									e.Bulk(sliced) = 1;
									//std::cout << "New edge: " << e.LocalID() << " " << e.getBeg().LocalID() << "<->" << e.getEnd().LocalID() << std::endl;
									split_edges.push_back(e);
									cnt++;
									alledges.push_back(e);
								}
								
								e = cutedges[indx[fedges[r]]];
								//std::cout << "Cut edge: " << e.LocalID() << " " << e.getBeg().LocalID() << "<->" << e.getEnd().LocalID() << std::endl;
								split_edges.push_back(e);
								
								if( e2n2.isValid() && e2n2 != e1n2 )
								{
									edge_nodes[0] = e1n2;
									edge_nodes[1] = e2n2;
									e = m.CreateEdge(edge_nodes).first;
									material[e] = mat;
									e.Bulk(sliced) = 1;
									//std::cout << "New edge: " << e.LocalID() << " " << e.getBeg().LocalID() << "<->" << e.getEnd().LocalID() << std::endl;
									split_edges.push_back(e);
									cnt++;
									alledges.push_back(e);
								}
								
								if( cnt )
								{
									f = m.CreateFace(split_edges).first;
									//std::cout << "New face " << f.LocalID();
									//std::cout << " [" << Element::GeometricTypeName(f.GetGeometricType()) << "]";
									//std::cout << " nodes (";
									ElementArray<Node> fnodes = f.getNodes();
									//for(int l = 0; l < fnodes.size(); ++l)
									//	std::cout << fnodes[l].LocalID() << " ";
									//std::cout << ") " << std::endl;
									material[f] = mat;
									f.Bulk(sliced) = 1;
									split_faces.push_back(f);
								}
							}
						}
					}
					cedges.RemMarker(isolate);
				}
				
				cedges.RemMarker(cmrk);
				/*
				for(INMOST_DATA_ENUM_TYPE q = 0; q < alledges.size(); ++q)
				{
					Storage::real_array a = alledges[q].getBeg().Coords();
					Storage::real_array b = alledges[q].getEnd().Coords();
					int aid = alledges[q].getBeg().LocalID();
					int bid = alledges[q].getEnd().LocalID();
					//std::cout << alledges[q].LocalID() << " ";
					//std::cout << "(" << aid << "," << a[0] << "," << a[1] << "," << a[2] << ")";
					//std::cout << " <-> ";
					//std::cout << "(" << bid << "," << b[0] << "," << b[1] << "," << b[2] << ")";
					//std::cout << std::endl;
				}
				*/
			}
			
			
			
			if( !split_faces.empty() )
			{
				int lid = it->LocalID();
				ElementArray<Cell> ret = Cell::SplitCell(it->self(),split_faces,0);
				
				
				//if(!simple )
				if( false )
				{
					std::cout << (simple?"simple":"complex") << " algorithm, split cell " << lid << " with " << split_faces.size() << " faces ";
					for(INMOST_DATA_ENUM_TYPE q = 0; q < split_faces.size(); ++q)
					{
						std::cout << split_faces[q].LocalID();
						std::cout << " [" << Element::GeometricTypeName(split_faces[q].GetGeometricType()) << "]";
						std::cout << " nodes (";
						ElementArray<Node> fnodes = split_faces[q].getNodes();
						for(INMOST_DATA_ENUM_TYPE l = 0; l < fnodes.size(); ++l)
							std::cout << fnodes[l].LocalID() << " ";
						std::cout << ") ";
					}
					std::cout << " result in " << ret.size() << " cells:";
					for(INMOST_DATA_ENUM_TYPE q = 0; q < ret.size(); ++q)
					{
						std::cout << " " << ret[q].LocalID();
						std::cout << " [" << Element::GeometricTypeName(ret[q].GetGeometricType()) << "]";
						std::cout << " nodes (";
						ElementArray<Node> cnodes = ret[q].getNodes();
						for(INMOST_DATA_ENUM_TYPE l = 0; l < cnodes.size(); ++l)
							std::cout << cnodes[l].LocalID() << " ";
						std::cout << ") faces (";
						ElementArray<Face> cfaces = ret[q].getFaces();
						for(INMOST_DATA_ENUM_TYPE l = 0; l < cfaces.size(); ++l)
						{
							std::cout << cfaces[l].LocalID() << " [" << Element::GeometricTypeName(cfaces[l].GetGeometricType()) << "] n{";
							ElementArray<Node> fcnodes = cfaces[l].getNodes();
							for(INMOST_DATA_ENUM_TYPE r = 0; r < fcnodes.size(); ++r)
								std::cout << fcnodes[r].LocalID() << " ";
							std::cout << "} ";
						}
						std::cout << ")";
					}
					std::cout << std::endl;
				}
				//if( ret.size() == 1 ) scanf("%*c");
				
				ret.SetMarker(mrk);
				nslice++;
				
				if( false )//ret.size() == 1 && !simple )
				{
					ElementArray<Edge> cedges = ret[0].getEdges();
					
					
					cedges.SetMarker(cmrk);
					
					std::cout << "Original cut edges: " << std::endl;
					for(INMOST_DATA_ENUM_TYPE k = 0; k < edges.size(); ++k) std::cout << edges[k].LocalID() << " ";
					std::cout << std::endl;
					
					std::cout << "Original cell edges: " << std::endl;
					for(INMOST_DATA_ENUM_TYPE k = 0; k < cedges.size(); ++k) std::cout << cedges[k].LocalID() << " ";
					std::cout << std::endl;
					
					std::cout << "Cut faces edges: " << std::endl;
					std::map<Edge,int> vstcnt;
					for(INMOST_DATA_ENUM_TYPE q = 0; q < split_faces.size(); ++q)
					{
						ElementArray<Edge> sedges = split_faces[q].getEdges();
						std::cout << "Face " << split_faces[q].LocalID() << ":";
						for(INMOST_DATA_ENUM_TYPE r = 0; r < sedges.size(); ++r)
						{
							vstcnt[sedges[r]]++;
							std::cout << " " << sedges[r].LocalID() << "(" << (sedges[r].GetMarker(cmrk)?"s":"i") << ")";
						}
						std::cout << std::endl;
					}
					std::cout  << "Visit counts:";
					for(std::map<Edge,int>::iterator jt = vstcnt.begin(); jt != vstcnt.end(); ++jt)
						std::cout << " EDGE:" << jt->first.LocalID() << "(" << (jt->first.GetMarker(cmrk)?"s":"i") <<  ") " << jt->second;
					std::cout <<std::endl;
					
					std::cout << "Cell faces edges: " << std::endl;
					ElementArray<Face> cfaces = ret[0].getFaces();
					for(INMOST_DATA_ENUM_TYPE k = 0; k < cfaces.size(); ++k)
					{
						ElementArray<Edge> sedges = cfaces[k].getEdges();
						std::cout << "Face " << cfaces[k].LocalID() << ":";
						for(INMOST_DATA_ENUM_TYPE r = 0; r < sedges.size(); ++r)
						{
							vstcnt[sedges[r]]++;
							std::cout << " " << sedges[r].LocalID() << "(" << (sedges[r].GetMarker(cmrk)?"s":"i") << ")";
						}
						std::cout << std::endl;
					}
					std::cout  << "Visit counts including cell's faces:";
					for(std::map<Edge,int>::iterator jt = vstcnt.begin(); jt != vstcnt.end(); ++jt)
						std::cout << " EDGE:" << jt->first.LocalID() << "(" << (jt->first.GetMarker(cmrk)?"s":"i") <<  ") " << jt->second;
					std::cout <<std::endl;
					
					
					
					std::cout << "Cell:" << std::endl;
					for(INMOST_DATA_ENUM_TYPE k = 0; k < cedges.size(); ++k)
					{
						Storage::real_array c1 = cedges[k].getBeg().Coords();
						Storage::real_array c2 = cedges[k].getEnd().Coords();
						std::cout << "(" << c1[0] << "," << c1[1] << "," << c1[2] << ")<->(" << c2[0] << "," << c2[1] << "," << c2[2] << ")" << std::endl;
					}
					std::cout << "Faces:" << std::endl;
					for(INMOST_DATA_ENUM_TYPE q = 0; q < split_faces.size(); ++q)
					{
						ElementArray<Edge> sedges = split_faces[q].getEdges();
						for(INMOST_DATA_ENUM_TYPE r = 0; r < sedges.size(); ++r)
						{
							Storage::real_array c1 = sedges[r].getBeg().Coords();
							Storage::real_array c2 = sedges[r].getEnd().Coords();
							std::cout << "(" << c1[0] << "," << c1[1] << "," << c1[2] << ")<->(" << c2[0] << "," << c2[1] << "," << c2[2] << ")" << std::endl;
						}
					}
					
					cedges.RemMarker(cmrk);
					
				}
			}
		}
	}
	m.ReleaseMarker(isolate);
	std::cout << "release marker isolate " << isolate << std::endl;
	m.ReleaseMarker(cmrk);
	std::cout << "release marker cmrk " << cmrk << std::endl;
	m.ReleaseMarker(visited);
	std::cout << "release marker visited " << visited << std::endl;
	m.DeleteTag(indx);

	std::cout << "sliced cells: " << nslice << std::endl;
	
	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;
	
	{
		int mat[3] = {0,0,0};
		int tot = 0;
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
		{
			mat[material[*it]]++;
			tot++;
		}
		std::cout << "face materials, 0: " << mat[0] << " 1: " << mat[1] << " 2: " << mat[2] << std::endl;
	}
	
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	{
		int mat[3] = {0,0,0};
		ElementArray<Face> faces = it->getFaces();
		for(int k = 0; k < (int)faces.size(); ++k)
			mat[material[faces[k]]]++;
		if( !(mat[0] == 0 || mat[1] == 0) )
		{
			material[*it] = 2;
			std::cout << "oops, materials for cell " << it->LocalID() << " faces were not split, 0: " << mat[0] << " ,1: " << mat[1] << " ,2: " << mat[2] << " slice edges " << it->getEdges(slice).size() << std::endl;
		}
		else if( mat[0] != 0 ) material[*it] = 0;
		else if( mat[1] != 0 ) material[*it] = 1;
		else
		{
			//double cnt[3];
			//it->Centroid(cnt);
			//double v = func(cnt[0],cnt[1],cnt[2],type);
			//it->Integer(material) = (v <= 0.0 ? 0 : 1);
			material[*it] = 2;
			//std::cout << "oops cannot determine material for cell, all faces have type 2, set to " << material[*it] << std::endl;
		}
	}
	
	
	{
		int mat[3] = {0,0,0};
		int tot = 0;
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			mat[material[*it]]++;
			tot++;
		}
		std::cout << "cell materials, 0: " << mat[0] << " 1: " << mat[1] << " 2: " << mat[2] << std::endl;
	}
	
	
	if( remove_material_zero )
	{
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
			if( material[*it] == 0 )//|| it->Integer(collapse) )
				it->Delete();
		
		for(ElementType etype = FACE; etype >= NODE; etype = PrevElementType(etype) )
			for(Mesh::iteratorElement it = m.BeginElement(etype); it != m.EndElement(); ++it)
				if( it->nbAdjElements(CELL) == 0 ) it->Delete();
	}

	m.ReleaseMarker(slice,NODE|EDGE);
	std::cout << "release marker slice " << slice << std::endl;
	m.ReleaseMarker(original,NODE|EDGE);
	std::cout << "release marker original " << original << std::endl;
	m.ReleaseMarker(mrk,EDGE|FACE|NODE|CELL);
	std::cout << "release marker mrk " << mrk << std::endl;
        
    // m.DeleteTag(material);
    // m.DeleteTag(level);
        
	//m.Save(grid_out);
	//return 0;
}
