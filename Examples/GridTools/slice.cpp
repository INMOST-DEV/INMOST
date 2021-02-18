#include "inmost.h"

using namespace INMOST;

int main(int argc, char ** argv)
{
	double nx = 2.0/7.0, ny = 6.0/7.0, nz = 3.0/7.0;
	double px = 0.5, py = 0.5, pz = 0.5;

	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid.pmf] [nx=0] [ny=0] [nz=1] [px=0.5] [py=0.5] [pz=0.5]" << std::endl;
		return -1;
	}

	std::string grid_out = "grid.pmf";

	if( argc > 2 ) grid_out = std::string(argv[2]);
	if( argc > 3 ) nx = atof(argv[3]);
	if( argc > 4 ) ny = atof(argv[4]);
	if( argc > 5 ) nz = atof(argv[5]);
	if( argc > 6 ) px = atof(argv[6]);
	if( argc > 7 ) py = atof(argv[7]);
	if( argc > 8 ) pz = atof(argv[8]);

	double d = nx*px+ny*py+nz*pz;

	Mesh m;
	m.Load(argv[1]);
	m.SetTopologyCheck(NEED_TEST_CLOSURE|PROHIBIT_MULTILINE|PROHIBIT_MULTIPOLYGON|GRID_CONFORMITY|DEGENERATE_EDGE|DEGENERATE_FACE|DEGENERATE_CELL | FACE_EDGES_ORDER);
	//m.RemTopologyCheck(THROW_EXCEPTION);
	Tag sliced = m.CreateTag("SLICED",DATA_BULK,FACE|EDGE|NODE,FACE|EDGE|NODE,1);

	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;

	MarkerType slice = m.CreateMarker();
	int nslice = 0, nmark = 0;
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
	{
		double p[3];
		Storage::real_array c0 = it->getBeg()->Coords();
		Storage::real_array c1 = it->getEnd()->Coords();
		double r0 = c0[0]*nx+c0[1]*ny+c0[2]*nz - d;
		double r1 = c1[0]*nx+c1[1]*ny+c1[2]*nz - d;
		//std::cout << "r0 " << r0 << " r1 " << r1 << std::endl;
		if( r0*r1 < -1.0e-12 )
		{
			p[0] = (r0*c1[0] - r1*c0[0])/(r0-r1);
			p[1] = (r0*c1[1] - r1*c0[1])/(r0-r1);
			p[2] = (r0*c1[2] - r1*c0[2])/(r0-r1);
			//std::cout << "p " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			Node n = m.CreateNode(p);
			n.Bulk(sliced) = 1;
			n.SetMarker(slice);
			bool was_sliced = it->HaveData(sliced) ? true : false;
			ElementArray<Edge> ret = Edge::SplitEdge(it->self(),ElementArray<Node>(&m,1,n.GetHandle()),0);
			if( was_sliced ) for(INMOST_DATA_ENUM_TYPE q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
			nslice++;
		}
		else
		{
			if( fabs(r0) < 1.0e-6 )
			{
				it->getBeg()->SetMarker(slice);
				nmark++;
			}
			if( fabs(r1) < 1.0e-6 )
			{
				it->getEnd()->SetMarker(slice);
				nmark++;
			}
		}
	}

	std::cout << "sliced edges: " << nslice << " marked nodes: " << nmark << std::endl;

	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;
	
	nslice = 0;
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
	{
		ElementArray<Node> nodes = it->getNodes(slice); //those nodes should be ordered so that each pair forms an edge
		if( nodes.size() > 1 ) // if there is 1, then only one vertex touches the plane
		{
			//if there is more then two, then original face is non-convex
			if( nodes.size() > 2 ) std::cout << "Looks like face " << it->LocalID() << " is nonconvex" << std::endl;
			else
			{
				Edge e = m.CreateEdge(nodes).first;
				e.Bulk(sliced) = 1;
				e.SetMarker(slice);
				bool was_sliced = it->HaveData(sliced) ? true : false;
				ElementArray<Face> ret = Face::SplitFace(it->self(),ElementArray<Edge>(&m,1,e.GetHandle()),0);
				if( was_sliced ) for(INMOST_DATA_ENUM_TYPE q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
				nslice++;
			}
		}
		//else std::cout << "Only one adjacent slice node, face " << it->LocalID() << std::endl;
	}

	nmark = 0;
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
		if( !it->GetMarker(slice) && it->getBeg()->GetMarker(slice) && it->getEnd()->GetMarker(slice) )
		{
			it->SetMarker(slice);
			nmark++;
		}

	std::cout << "sliced faces: " << nslice << " marked edges: " << nmark << std::endl;

	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;
		
	
	nslice = 0;
	MarkerType visited = m.CreateMarker();
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	{
		ElementArray<Edge> edges = it->getEdges(slice);
		if( edges.size() >= 3 ) //these should form a triangle
		{
			//order edges
			ElementArray<Edge> order_edges(&m);
			order_edges.push_back(edges[0]);
			order_edges.SetMarker(visited);
			while(order_edges.size() != edges.size() )
			{
				for(INMOST_DATA_ENUM_TYPE k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(visited) )
				{
					if( edges[k]->getBeg() == order_edges.back()->getBeg() || edges[k]->getBeg() == order_edges.back()->getEnd() ||
						edges[k]->getEnd() == order_edges.back()->getBeg() || edges[k]->getEnd() == order_edges.back()->getEnd() )
					{
						order_edges.push_back(edges[k]);
						order_edges.back().SetMarker(visited);
					}
				}
			}
			edges.RemMarker(visited);
			Face f = m.CreateFace(order_edges).first;
			f.Bulk(sliced) = 1;
			Cell::SplitCell(it->self(),ElementArray<Face>(&m,1,f.GetHandle()),0);
			nslice++;
		}
	}
	m.ReleaseMarker(visited);

	std::cout << "sliced cells: " << nslice << std::endl;
	
	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;

	Tag material = m.CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);

	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	{
		double cnt[3];
		it->Centroid(cnt);
		double v = cnt[0]*nx+cnt[1]*ny+cnt[2]*nz-d;
		if( v < 0.0 )
			it->Integer(material) = 0;
		else
			it->Integer(material) = 1;
	}
	

	m.ReleaseMarker(slice,NODE|EDGE);

	m.Save(grid_out);
	return 0;
}
