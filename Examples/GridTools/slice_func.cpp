#include "inmost.h"

using namespace INMOST;

double func(double x, double y, double z, int n)
{
	if( n == 0 )
		return sqrt(x*x+y*y)-0.25;
	else if( n == 1 )
		return -(sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))-0.5);
	else if( n == 2 )
	{
		if( x > 0.5 )
			return -(sqrt((x-0.7)*(x-0.7)+(y-0.2)*(y-0.2)+(z-0.4)*(z-0.4))-0.4);
		else
			return sqrt((x-0.3)*(x-0.3)+(y-0.4)*(y-0.4)+(z-0.5)*(z-0.5))-0.3;
	}
	else if( n == 3 )
		return y-(4*x*x*x*x-2*x*x+0.5)+z*z*z;
	else if( n == 4 )
		return y-(10*x*x*x*x-8*x*x*x-5*x*x+0.2+4*x)+z*z*z;
	else if( n == 5 )
	{
		
		//double Lx = 0.2, Rx = 0.4, Ly = 0.1, Ry = 0.3;
		/*
		double Lx = 0., Rx = 2, Ly = 0.0, Ry = 4.5;
		if (x > Rx){
			if (y > Ry) return -sqrt( (x-Rx)*(x-Rx) + (y-Ry)*(y-Ry) );
			else return -(x-Rx);
		}
		if (x < Lx){
			if (y < Ly)  return -sqrt( (x-Lx)*(x-Lx) + (y-Ly)*(y-Ly) );
			else return (x - Lx);
		}
		if (y > Ry) return Ry - y;
		if (y < Ly) return y - Ly;
		return fmin( fmin(x-Lx,Rx-x), fmin(y-Ly, Ry-y) );
		*/
		
		double Lx = 0., Rx = 2, Ly = 0.0, Ry = 4.5;
	    if (x > Rx)
	    {
	        if (y > Ry) return -sqrt( (x-Rx)*(x-Rx) + (y-Ry)*(y-Ry) );
	        if (y < Ly) return -sqrt( (x-Rx)*(x-Rx) + (y-Ly)*(y-Ly) );
	        return -(x-Rx);
	    }
	    if (x < Lx)
	    {
	        if (y < Ly)  return -sqrt( (x-Lx)*(x-Lx) + (y-Ly)*(y-Ly) );
	        if (y > Ry)  return -sqrt( (x-Lx)*(x-Lx) + (y-Ry)*(y-Ry) );
	        else return (x - Lx);
	    }
	    if (y > Ry) return Ry - y;
	    if (y < Ly) return y - Ly;
	    return fmin( fmin(x-Lx,Rx-x), fmin(y-Ly, Ry-y) );
	    
	}
	else if( n == 6 )
	{
		return std::min(std::min(sqrt((x-0.48)*(x-0.48) + (z-1)*(z-1))-0.25,sqrt((y-0.53)*(y-0.53) + (z-3)*(z-3))-0.24),sqrt((x-0.6)*(x-0.6) + (z-5)*(z-5))-0.3);
	}
	return 1;
}

double search(double r0, double r1, double c0[3], double c1[3], double p[3],int type, bool binary = false)
{
	double rp = 1.0e20, rp_min = 1.0e20, p_min[3];
	int iters = 0;
	do
	{
		double m1 = r0/(r0-r1), m0;
		if( m1 < 0.0 || m1 > 1.0 || binary ) m1 = 0.5;
		m0 = 1.0-m1;
		p[0] = m1*c1[0] + m0*c0[0];
		p[1] = m1*c1[1] + m0*c0[1];
		p[2] = m1*c1[2] + m0*c0[2];
		rp = func(p[0],p[1],p[2],type);
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
	} while( fabs(rp) > 1.0e-9 && iters < 100 );
	
	if( rp > rp_min )
	{
		p[0] = p_min[0];
		p[1] = p_min[1];
		p[2] = p_min[2];
		rp = rp_min;
	}
	
	return rp;
}


double search_zero(double r0, double r1, double c0[3], double c1[3], double p[3],int type)
{
	int iters = 0;
	double rp = 1.0e+20;
	do
	{
		p[0] = 0.5*c1[0] + 0.5*c0[0];
		p[1] = 0.5*c1[1] + 0.5*c0[1];
		p[2] = 0.5*c1[2] + 0.5*c0[2];
		rp = func(p[0],p[1],p[2],type);
		if( fabs(rp) < 1.0e-8 )
		{
			if( fabs(r0) < 1.0e-8 )
			{
				c0[0] = p[0];
				c0[1] = p[1];
				c0[2] = p[2];
				r0 = rp;
			}
			else
			{
				c1[0] = p[0];
				c1[1] = p[1];
				c1[2] = p[2];
				r1 = rp;
			}
		}
		else
		{
			if( fabs(r0) > 1.0e-8 )
			{
				c0[0] = p[0];
				c0[1] = p[1];
				c0[2] = p[2];
				r0 = rp;
			}
			else
			{
				c1[0] = p[0];
				c1[1] = p[1];
				c1[2] = p[2];
				r1 = rp;
			}
		}
		iters++;
	} while( iters < 20 );
	
	return rp;
}

int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh  [mesh_out=grid.pmf] [type=perforated_strip]" << std::endl;
		return -1;
	}
	
	std::string grid_out = "grid.pmf";

	int type = 0;
	if( argc > 2 ) grid_out = std::string(argv[2]);
	if( argc > 3 ) type = atoi(argv[3]);
	
	

	Mesh m;
	m.Load(argv[1]);
	m.SetTopologyCheck(NEED_TEST_CLOSURE|PROHIBIT_MULTILINE|PROHIBIT_MULTIPOLYGON|GRID_CONFORMITY|DEGENERATE_EDGE|DEGENERATE_FACE|DEGENERATE_CELL | FACE_EDGES_ORDER | MARK_ON_ERROR);
	//m.RemTopologyCheck(THROW_EXCEPTION);
	TagInteger material = m.CreateTag("MATERIAL",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag sliced = m.CreateTag("SLICED",DATA_BULK,FACE|EDGE|NODE,FACE|EDGE|NODE,1);

	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	
	MarkerType original = m.CreateMarker();
	for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it) it->SetMarker(original);
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it) it->SetMarker(original);

	MarkerType slice = m.CreateMarker();
	MarkerType mrk = m.CreateMarker();
	int nslice = 0, nmark = 0;
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it) if( !it->GetMarker(mrk) )
	{
		double p[3],pc0[3],pc1[3];
		Storage::real_array c0 = it->getBeg()->Coords();
		Storage::real_array c1 = it->getEnd()->Coords();
		double r0 = func(c0[0],c0[1],c0[2],type);
		double r1 = func(c1[0],c1[1],c1[2],type);
		material[it->getBeg()] = (r0 <= 0)? 0 : 1;
		material[it->getEnd()] = (r1 <= 0)? 0 : 1;
		//std::cout << "r0 " << r0 << " r1 " << r1 << std::endl;
		if( (r0*r1 < -1.0e-12) || (fabs(r0*r1) < 1.0e-12 && ((fabs(r0) < 1.0e-12) ^ (fabs(r1) < 1.0e-12))) )
		{
			pc0[0] = c0[0], pc0[1] = c0[1], pc0[2] = c0[2];
			pc1[0] = c1[0], pc1[1] = c1[1], pc1[2] = c1[2];
			if((fabs(r0) < 1.0e-12) ^ (fabs(r1) < 1.0e-12))
			{
				search_zero(r0,r1,pc0,pc1,p,type);
			}
			else
			{
				double rp = search(r0,r1,pc0,pc1,p,type);
				if( rp > 1.0e-3 ) //cannot find intersection
				{
					rp = search(r0,r1,pc0,pc1,p,type,true);
					p[0] = c0[0]*0.5+c1[0]*0.5;
					p[1] = c0[1]*0.5+c1[1]*0.5;
					p[2] = c0[2]*0.5+c1[2]*0.5;
				}
			}
			//p[0] = (r0*c1[0] - r1*c0[0])/(r0-r1);
			//p[1] = (r0*c1[1] - r1*c0[1])/(r0-r1);
			//p[2] = (r0*c1[2] - r1*c0[2])/(r0-r1);
			//std::cout << "p " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			//distance to the corners
			
			double l0 = 0, l1 = 0, l;
			for(int r = 0; r < 3; ++r)
			{
				l0 += (p[r]-c0[r])*(p[r]-c0[r]);
				l1 += (p[r]-c1[r])*(p[r]-c1[r]);
			}
			l0 = sqrt(l0);
			l1 = sqrt(l1);
			l = l0+l1;
			if( l0 < 1.0e-5*l )
			{
				material[it->getBeg()] = 2;
				it->getBeg()->SetMarker(slice);
				nmark++;
			}
			else if( l1 < 1.0e-5*l )
			{
				material[it->getEnd()] = 2;
				it->getEnd()->SetMarker(slice);
				nmark++;
			}
			else
			{
				Node n = m.CreateNode(p);
				material[n] = 2;
				n.Bulk(sliced) = 1;
				n.SetMarker(slice);
				bool was_sliced = it->HaveData(sliced) ? true : false;
				ElementArray<Edge> ret = Edge::SplitEdge(it->self(),ElementArray<Node>(&m,1,n.GetHandle()),0);
				ret.SetMarker(mrk);
				if( was_sliced ) for(int q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
				nslice++;
			}
		}
		else
		{
			if( fabs(r0) < 1.0e-6 )
			{
				material[it->getBeg()] = 2;
				it->getBeg()->SetMarker(slice);
				nmark++;
			}
			if( fabs(r1) < 1.0e-6 )
			{
				material[it->getEnd()] = 2;
				it->getEnd()->SetMarker(slice);
				nmark++;
			}
		}
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
			std::cout << "oops, materials for edge nodes were not split, 0: " << mat[0] << " ,1: " << mat[1] << " ,2: " << mat[2] << std::endl;
		}
		else if( mat[0] != 0 ) material[*it] = 0;
		else if( mat[1] != 0 ) material[*it] = 1;
		else material[*it] = 2;
	}
	
	
	nslice = 0;
	MarkerType unique = m.CreateMarker();
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
				if( was_sliced ) for(int q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
				nslice++;
			}
			else if( it->nbAdjElements(NODE,slice) != it->nbAdjElements(NODE) ) //not entire face is sliced
			{
				
				//std::cout << "Sliced nodes " << nodes.size() << " total nodes: ";
				nodes = it->getNodes();
				//std::cout << nodes.size() << std::endl;
				for(int q = 0; q < (int)nodes.size(); ++q)
				{
					Storage::real_array c0 = nodes[q].Coords();
					double r0 = func(c0[0],c0[1],c0[2],type);
					//std::cout << "NODE:" << nodes[q].LocalID() << " " << (nodes[q]->GetMarker(slice)?"":"not ") << "slice " << (nodes[q]->GetMarker(original)?"":"not ") << "original r="<<r0 << " material " << material[nodes[q]] << std::endl;
				}
				ElementArray<Edge> fedges = it->getEdges();
				for(int q = 0; q < (int)fedges.size(); ++q)
				{
					//std::cout << "EDGE:" << fedges[q].LocalID() << " " << fedges[q].getBeg().LocalID() << "<->" << fedges[q].getEnd().LocalID() << " " << (fedges[q]->GetMarker(slice)?"":"not ") << "slice material " << material[fedges[q]] << std::endl;
				}

				
				double c0[3],c1[3],pc0[3],pc1[3],p[3];
				it->Centroid(c0);
				double r0 = func(c0[0],c0[1],c0[2],type);
				int material0 = (r0 <= 0)? 0 : 1;
				bool s0 = false;
				if( fabs(r0) < 1.0e-6 ) s0 = true;
				
				//std::cout << "Centernode " << (s0?"sliced":"") << " r=" << r0 << std::endl;
				Node centernode = InvalidNode();
				ElementArray<Edge> split_edges(&m);
				ElementArray<Node> cutnodes(&m,nodes.size()), edge_nodes(&m,2);
				
				//calculate nodes that cut along the edges connecting centernode
				for(int q = 0; q < (int)nodes.size(); ++q) if( !nodes[q].GetMarker(slice) )
				{
					nodes[q].Centroid(c1);
					double r1 = func(c1[0],c1[1],c1[2],type);
					//std::cout << "NODE:" << nodes[q].LocalID() << " r0 " << r0 << " r1 " << r1 << " r0*r1 " << r0*r1 << std::endl;
					if( (r0*r1 < -1.0e-12) || (fabs(r0*r1) < 1.0e-12 && ((fabs(r0) < 1.0e-12) ^ (fabs(r1) < 1.0e-12))))
					{
						pc0[0] = c0[0], pc0[1] = c0[1], pc0[2] = c0[2];
						pc1[0] = c1[0], pc1[1] = c1[1], pc1[2] = c1[2];
						if((fabs(r0) < 1.0e-12) ^ (fabs(r1) < 1.0e-12))
						{
							search_zero(r0,r1,pc0,pc1,p,type);
						}
						else
						{						
							double rp = search(r0,r1,pc0,pc1,p,type);
							if( rp > 1.0e-3 )
							{
								//std::cout << "inaccurate search " << rp;
								rp = search(r0,r1,pc0,pc1,p,type,true);
								//std::cout << " binary " << rp << std::endl;
								p[0] = c0[0]*0.5+c1[0]*0.5;
								p[1] = c0[1]*0.5+c1[1]*0.5;
								p[2] = c0[2]*0.5+c1[2]*0.5;
								
							}
						}
						//distance to the center node
						double l0 = 0, l1 = 0, l;
						for(int r = 0; r < 3; ++r)
						{
							l0 += (p[r]-c0[r])*(p[r]-c0[r]);
							l1 += (p[r]-c1[r])*(p[r]-c1[r]);
						}
						l0 = sqrt(l0);
						l1 = sqrt(l1);
						l = l0+l1;
						if( l0 < 1.0e-3*l ) //edge goes through centernode
						{
							if( !centernode.isValid() )
								centernode = m.CreateNode(c0);
							cutnodes[q] = centernode;
							//std::cout << "selected centernode " << std::endl;
						}
						else if( l1 > 1.0e-3*l )
						{
							cutnodes[q] = m.CreateNode(p);
							//std::cout << "created new node " << std::endl;
						}
						else
						{
							material[nodes[q]] = 2;
							nodes[q].SetMarker(slice);
						}
						//else std::cout << "use old node " << std::endl;
					}
					else if( s0 )
					{
						if( !centernode.isValid() ) centernode = m.CreateNode(c0);
						cutnodes[q] = centernode;
					}
				}
				
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
					else if( cutnodes[i1].isValid() && cutnodes[i2].isValid() )
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
					for(int q = 0; q < split_edges.size(); ++q)
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
					ElementArray<Face> ret = Face::SplitFace(it->self(),split_edges,0);
					ret.SetMarker(mrk);
					//std::cout << "New faces: " << ret.size() << ":";
					for(int q = 0; q < ret.size(); ++q)
					{
						int mat[3] = {0,0,0};
						ElementArray<Edge> fe = ret[q].getEdges();
						for( int l = 0; l < fe.size(); ++l) mat[material[fe[l]]]++;
						//std::cout << " " << ret[q].LocalID() << "[" << mat[0] << "," << mat[1] << "," << mat[2] << "]";
					}
					//std::cout << std::endl;
					if( was_sliced ) for(int q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
					nslice++;
				}
				//else std::cout << "No split edges " << std::endl;
			}
			//else std::cout << "No new edges on face " << it->LocalID() << std::endl;
		}
		//else std::cout << "Only one adjacent slice node, face " << it->LocalID() << std::endl;
	}
	m.ReleaseMarker(unique);
	

	nmark = 0;
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
		if( !it->GetMarker(slice) && it->getBeg()->GetMarker(slice) && it->getEnd()->GetMarker(slice) )
		{
			if( material[*it] != 2 ) std::cout << "Edge supposed to get material 2, but have " << material[*it] << std::endl;
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
			std::cout << "oops, materials for face " << it->LocalID() << " edges were not split, 0: " << mat[0] << " ,1: " << mat[1] << " ,2: " << mat[2] << std::endl;
		}
		else if( mat[0] != 0 ) material[*it] = 0;
		else if( mat[1] != 0 ) material[*it] = 1;
		else material[*it] = 2;
	}
	
	nslice = 0;
	TagInteger indx = m.CreateTag("TEMP_INDX",DATA_INTEGER,NODE|EDGE,NONE,1);
	MarkerType visited = m.CreateMarker();
	MarkerType cmrk = m.CreateMarker();
	MarkerType isolate = m.CreateMarker();
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
			
			if( simple )
			{
				//gather loop by loop and create faces
				int loop_cnt = 0;
				//order edges
				ElementArray<Edge> order_edges(&m);
				
				Node last;
				int nvisited = 0;
				bool found = false;
				for(int k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(original) )
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
					for(int k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(visited) )
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
							for(int k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(visited) && !edges[k]->GetMarker(original))
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
				
				for(int k = 0; k < (int)cnodes.size(); ++k) indx[cnodes[k]] = k;
				for(int k = 0; k < (int)cedges.size(); ++k) indx[cedges[k]] = k;
				
				double c0[3],c1[3],pc0[3],pc1[3],p[3];
				it->Centroid(c0);
				double r0 = func(c0[0],c0[1],c0[2],type);
				int material0 = (r0 <= 0)? 0 : 1;
				bool s0 = false;
				if( fabs(r0) < 1.0e-6 ) s0 = true;
				
				//std::cout << "Number of cut edges: " << edges.size() << std::endl;
				//for(std::map<Node,int>::iterator jt = visit_count.begin(); jt != visit_count.end(); ++jt)
				//	std::cout << "NODE:" << jt->first.LocalID() << " visited " << jt->second << " times" << std::endl;
				//std::cout << "Centernode " << (s0?"sliced":"") << " r=" << r0 << std::endl;
				Node centernode = InvalidNode();
				
				
				//calculate nodes that cut along the edges connecting centernode
				for(int q = 0; q < (int)cnodes.size(); ++q)
				{
					if( !cnodes[q].GetMarker(slice) )
					{
						cnodes[q].Centroid(c1);
						double r1 = func(c1[0],c1[1],c1[2],type);
						//std::cout << "NODE:" << cnodes[q].LocalID() << " r0 " << r0 << " r1 " << r1 << " r0*r1 " << r0*r1 << " " << (cnodes[q].GetMarker(slice)?"":"not ") << "sliced" << std::endl;
						
						if( (r0*r1 < -1.0e-12) || (fabs(r0*r1) < 1.0e-12 && ((fabs(r0) < 1.0e-12) ^ (fabs(r1) < 1.0e-12))))
						{
							pc0[0] = c0[0], pc0[1] = c0[1], pc0[2] = c0[2];
							pc1[0] = c1[0], pc1[1] = c1[1], pc1[2] = c1[2];
							if((fabs(r0) < 1.0e-12) ^ (fabs(r1) < 1.0e-12))
							{
								search_zero(r0,r1,pc0,pc1,p,type);
							}
							else
							{						
								double rp = search(r0,r1,pc0,pc1,p,type);
								if( rp > 1.0e-3 )
								{
									//std::cout << "inaccurate search " << rp;
									rp = search(r0,r1,pc0,pc1,p,type,true);
									//std::cout << " binary " << rp << std::endl;
									p[0] = c0[0]*0.5+c1[0]*0.5;
									p[1] = c0[1]*0.5+c1[1]*0.5;
									p[2] = c0[2]*0.5+c1[2]*0.5;
									
								}
							}
							//distance to the center node
							double l0 = 0, l1 = 0, l;
							for(int r = 0; r < 3; ++r)
							{
								l0 += (p[r]-c0[r])*(p[r]-c0[r]);
								l1 += (p[r]-c1[r])*(p[r]-c1[r]);
							}
							l0 = sqrt(l0);
							l1 = sqrt(l1);
							l = l0+l1;
							if( l0 < 1.0e-3*l ) //edge goes through centernode
							{
								if( !centernode.isValid() )
									centernode = m.CreateNode(c0);
								cutnodes[q] = centernode;
								//std::cout << "selected centernode " << cutnodes[q].LocalID() << std::endl;
							}
							else if( l1 > 1.0e-3*l )
							{
								cutnodes[q] = m.CreateNode(p);
								//std::cout << "created new node " << cutnodes[q].LocalID() << std::endl;
							}
							else
							{
								material[cnodes[q]] = 2;
								cnodes[q].SetMarker(slice);
								//std::cout << "use old node " << std::endl;
							}
						}
						else if( s0 )
						{
							if( !centernode.isValid() ) centernode = m.CreateNode(c0);
							cutnodes[q] = centernode;
						}
					}
					//else cutnodes[q] = cnodes[q];
				}
				
				for(int q = 0; q < (int)cutnodes.size(); ++q) if( cutnodes[q].isValid() )
				{
					//std::cout << "New cut node " << cutnodes[q].LocalID() << " at " << q << " on line with " << cnodes[q].LocalID() << " r=" << func(cutnodes[q].Coords()[0],cutnodes[q].Coords()[1],cutnodes[q].Coords()[2],type) << std::endl;
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
				
				for(int q = 0; q < (int)cutedges.size(); ++q) if( cutedges[q].isValid() )
				{
					//std::cout << "New cut edge " << cutedges[q].LocalID() << " at " << q << " on plane with " << cedges[q].LocalID() << " " << cutedges[q].getBeg().LocalID() << "<->" << cutedges[q].getEnd().LocalID() << " original " << cedges[q].getBeg().LocalID() << "<->" << cedges[q].getEnd().LocalID() << std::endl;
					Edge n = cutedges[q];
					material[n] = 2;
					n.Bulk(sliced) = 1;
					n.SetMarker(slice);
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
							for(int r = 0; r < split_edges.size(); ++r) vstcnt[split_edges[r]]++;
							Face f = m.CreateFace(split_edges).first;
							//std::cout << "Created face " << f.LocalID() << " with " << split_edges.size() << " edges: ";
							//for(int r = 0; r < split_edges.size(); ++r) std::cout << "EDGE:" << split_edges[r].LocalID()  << " ";
							//std::cout << " original FACE:" << cfaces[q].LocalID() << std::endl;
							material[f] = 2;
							f.Bulk(sliced) = 1;
							split_faces.push_back(f);
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
								}
								
								if( cnt )
								{
									f = m.CreateFace(split_edges).first;
									//std::cout << "New face " << f.LocalID() << std::endl;
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
			}
			
			if( !split_faces.empty() )
			{
				int lid = it->LocalID();
				ElementArray<Cell> ret = Cell::SplitCell(it->self(),split_faces,0);
				
				
				if(!simple )
				{
					//std::cout << (simple?"simple":"complex") << " algorithm, split cell " << lid << " with " << split_faces.size() << " faces ";
					//std::cout << " result in " << ret.size() << " cells:";
					//for(int q = 0; q < (int)ret.size(); ++q) std::cout << " " << ret[q].LocalID();
					//std::cout << std::endl;
				}
				
				ret.SetMarker(mrk);
				nslice++;
				
				if( false )//ret.size() == 1 && !simple )
				{
					ElementArray<Edge> cedges = ret[0].getEdges();
					
					
					cedges.SetMarker(cmrk);
					
					std::cout << "Original cut edges: " << std::endl;
					for(int k = 0; k < edges.size(); ++k) std::cout << edges[k].LocalID() << " ";
					std::cout << std::endl;
					
					std::cout << "Original cell edges: " << std::endl;
					for(int k = 0; k < cedges.size(); ++k) std::cout << cedges[k].LocalID() << " ";
					std::cout << std::endl;
					
					std::cout << "Cut faces edges: " << std::endl;
					std::map<Edge,int> vstcnt;
					for(int q = 0; q < split_faces.size(); ++q)
					{
						ElementArray<Edge> sedges = split_faces[q].getEdges();
						std::cout << "Face " << split_faces[q].LocalID() << ":";
						for(int r = 0; r < sedges.size(); ++r)
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
					for(int k = 0; k < cfaces.size(); ++k)
					{
						ElementArray<Edge> sedges = cfaces[k].getEdges();
						std::cout << "Face " << cfaces[k].LocalID() << ":";
						for(int r = 0; r < sedges.size(); ++r)
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
					for(int k = 0; k < cedges.size(); ++k)
					{
						Storage::real_array c1 = cedges[k].getBeg().Coords();
						Storage::real_array c2 = cedges[k].getEnd().Coords();
						std::cout << "(" << c1[0] << "," << c1[1] << "," << c1[2] << ")<->(" << c2[0] << "," << c2[1] << "," << c2[2] << ")" << std::endl;
					}
					std::cout << "Faces:" << std::endl;
					for(int q = 0; q < split_faces.size(); ++q)
					{
						ElementArray<Edge> sedges = split_faces[q].getEdges();
						for(int r = 0; r < sedges.size(); ++r)
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
	m.ReleaseMarker(cmrk);
	m.ReleaseMarker(visited);
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
			std::cout << "oops cannot determine material for cell, all faces have type 2, set to " << material[*it];
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
	
	Tag collapse = m.CreateTag("COLLAPSE",DATA_INTEGER,CELL,NONE,1);
	Tag el_mat = m.CreateTag("ELLIPSE_MATRIX",DATA_REAL,CELL,CELL,9);
	Tag el_cnt = m.CreateTag("ELLIPSE_CENTER",DATA_REAL,CELL,CELL,3);
	Tag level = m.CreateTag("LEVEL",DATA_REAL,CELL|FACE|EDGE|NODE,NONE,1);
	
	int ncollapsed = 0;
	for(Mesh::iteratorElement it = m.BeginElement(NODE|EDGE|FACE|CELL); it != m.EndElement(); ++it)
	{
		double cnt[3];
		it->Centroid(cnt);
		it->Real(level) = func(cnt[0],cnt[1],cnt[2],type);
	}
	
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	{
		double vv, v, nvv, cnt[3];
		v = it->Volume();
		vv = v;
		nvv = 1;
		ElementArray<Cell> cells = it->NeighbouringCells();
		for(ElementArray<Cell>::iterator jt = cells.begin(); jt != cells.end(); ++jt)
		{
			vv += jt->Volume();
			nvv++;
		}
		vv /= nvv;
		if( v < vv*0.01 )
		{
			if( false )
			{
				std::cout << "looks like cell " << it->LocalID() << " is too small, volume " << it->Volume() << std::endl;
				for(ElementArray<Cell>::iterator jt = cells.begin(); jt != cells.end(); ++jt)
				{
					std::cout << "adjacent cell " << jt->LocalID() << " volume " << jt->Volume() << std::endl;
				}
			}
			it->Integer(collapse) = 1;
			
			/*
			it->Centroid(cnt);
			vMatrix A(3,3), b(3,1);
			A(0,0) = unknown(1,0);
			A(1,1) = unknown(1,1);
			A(2,2) = unknown(1,2);
			A(1,2) = A(2,1) = unknown(0,3);
			A(1,3) = A(3,1) = unknown(0,4);
			A(2,3) = A(3,2) = unknown(0,5);
			b(0,0) = unknown(cnt[0],6);
			b(1,0) = unknown(cnt[1],7);
			b(2,0) = unknown(cnt[2],8);
			
			ElementArray<Node> nodes = it->getNodes();
			*/
			
			
			//find biggest face
			/*
			ElementArray<Face> faces = it->getFaces();
			Cell n = InvalidCell();
			double area = 0;
			for(int k = 0; k < (int)faces.size(); ++k)
			{
				Cell nn = it->Neighbour(faces[k]);
				if( nn.isValid() )
				{
					if( faces[k].Area() > area )
					{
						n = nn;
						area = faces[k].Area();
					}
				}
			}
			//merge cells
			if( n.isValid() )
			{
				ElementArray<Cell> unite(&m,2);
				unite[0] = it->self();
				unite[1] = n;
				Cell::UniteCells(unite,0);
			}
			 */
		}
		
	}
	/*
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		if( material[*it] == 0 || it->Integer(collapse) )
			it->Delete();
	
	for(ElementType etype = FACE; etype >= NODE; etype = PrevElementType(etype) )
		for(Mesh::iteratorElement it = m.BeginElement(etype); it != m.EndElement(); ++it)
			if( it->nbAdjElements(CELL) == 0 ) it->Delete();
	*/

	m.ReleaseMarker(slice,NODE|EDGE);
	m.ReleaseMarker(original,NODE|EDGE);
	m.ReleaseMarker(mrk,EDGE|FACE|NODE);
	m.Save(grid_out);
	return 0;
}
