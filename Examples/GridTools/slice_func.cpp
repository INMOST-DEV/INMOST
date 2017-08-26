#include "inmost.h"

using namespace INMOST;


double func(double x, double y, double z, int n)
{
	if( n == 0 )
		return sqrt(x*x+y*y)-0.25;
	else if( n == 1 )
		return -(sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))-0.5);
	return 1;
}

void search(double r0, double r1, double c0[3], double c1[3], double p[3],int type)
{
	double rp;
	do
	{
		p[0] = (r0*c1[0] - r1*c0[0])/(r0-r1);
		p[1] = (r0*c1[1] - r1*c0[1])/(r0-r1);
		p[2] = (r0*c1[2] - r1*c0[2])/(r0-r1);
		rp = func(p[0],p[1],p[2],type);
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
	} while( fabs(rp) > 1.0e-9 );
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
	m.SetTopologyCheck(NEED_TEST_CLOSURE|PROHIBIT_MULTILINE|PROHIBIT_MULTIPOLYGON|GRID_CONFORMITY|DEGENERATE_EDGE|DEGENERATE_FACE|DEGENERATE_CELL | FACE_EDGES_ORDER);
	//m.RemTopologyCheck(THROW_EXCEPTION);
	Tag material = m.CreateTag("MATERIAL",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
	Tag sliced = m.CreateTag("SLICED",DATA_BULK,FACE|EDGE|NODE,FACE|EDGE|NODE,1);

	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	
	MarkerType original = m.CreateMarker();
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
		it->getBeg()->Integer(material) = (r0 <= 0)? 0 : 1;
		it->getEnd()->Integer(material) = (r1 <= 0)? 0 : 1;
		//std::cout << "r0 " << r0 << " r1 " << r1 << std::endl;
		if( r0*r1 < -1.0e-12 )
		{
			pc0[0] = c0[0], pc0[1] = c0[1], pc0[2] = c0[2];
			pc1[0] = c1[0], pc1[1] = c1[1], pc1[2] = c1[2];
			search(r0,r1,pc0,pc1,p,type);
			//p[0] = (r0*c1[0] - r1*c0[0])/(r0-r1);
			//p[1] = (r0*c1[1] - r1*c0[1])/(r0-r1);
			//p[2] = (r0*c1[2] - r1*c0[2])/(r0-r1);
			//std::cout << "p " << p[0] << " " << p[1] << " " << p[2] << std::endl;
			Node n = m.CreateNode(p);
			n->Integer(material) = 2;
			n.Bulk(sliced) = 1;
			n.SetMarker(slice);
			bool was_sliced = it->HaveData(sliced) ? true : false;
			ElementArray<Edge> ret = Edge::SplitEdge(it->self(),ElementArray<Node>(&m,1,n.GetHandle()),0);
			ret.SetMarker(mrk);
			if( was_sliced ) for(int q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
			nslice++;
		}
		else
		{
			if( fabs(r0) < 1.0e-6 )
			{
				it->getBeg()->Integer(material) = 2;
				it->getBeg()->SetMarker(slice);
				nmark++;
			}
			if( fabs(r1) < 1.0e-6 )
			{
				it->getEnd()->Integer(material) = 2;
				it->getEnd()->SetMarker(slice);
				nmark++;
			}
		}
	}
	
	m.ReleaseMarker(mrk,EDGE);

	std::cout << "sliced edges: " << nslice << " marked nodes: " << nmark << std::endl;

	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;
	
	{
		int mat[3] = {0,0,0};
		int tot = 0;
		for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
		{
			mat[it->Integer(material)]++;
			tot++;
		}
		std::cout << "node materials, 0: " << mat[0] << " 1: " << mat[1] << " 2: " << mat[2] << std::endl;
	}
	
	
	
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
	{
		int mat[3] = {0,0,0};
		mat[it->getBeg()->Integer(material)]++;
		mat[it->getEnd()->Integer(material)]++;
		if( !(mat[0] == 0 || mat[1] == 0) )
		{
			it->Integer(material) = 2;
			std::cout << "oops, materials for edge nodes were not split, 0: " << mat[0] << " ,1: " << mat[1] << " ,2: " << mat[2] << std::endl;
		}
		else if( mat[0] != 0 ) it->Integer(material) = 0;
		else if( mat[1] != 0 ) it->Integer(material) = 1;
		else it->Integer(material) = 2;
	}
	
	
	nslice = 0;
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
	{
		ElementArray<Node> nodes = it->getNodes(slice); //those nodes should be ordered so that each pair forms an edge
		if( nodes.size() > 1 ) // if there is 1, then only one vertex touches the plane
		{
			int nsliced = 0;
			for(int q = 0; q < (int)nodes.size(); ++q) if( !nodes[q].GetMarker(original) ) nsliced++;
			//if there is more then two, then original face is non-convex
			if( nsliced && nodes.size() > 2 ) std::cout << "Looks like face " << it->LocalID() << " is nonconvex, there is " << nodes.size() << " nodes, out of them " << nsliced << " new cuts on edge" << std::endl;
			else if( nodes.size() == 2 )
			{
				Edge e = m.CreateEdge(nodes).first;
				e.Integer(material) = 2; //on slice
				e.Bulk(sliced) = 1;
				e.SetMarker(slice);
				bool was_sliced = it->HaveData(sliced) ? true : false;
				ElementArray<Face> ret = Face::SplitFace(it->self(),ElementArray<Edge>(&m,1,e.GetHandle()),0);
				if( was_sliced ) for(int q = 0; q < ret.size(); ++q) ret[q]->Bulk(sliced) = 1;
				nslice++;
			}
			//else std::cout << "No new edges on face " << it->LocalID() << std::endl;
		}
		//else std::cout << "Only one adjacent slice node, face " << it->LocalID() << std::endl;
	}
	
	m.ReleaseMarker(original,NODE);

	nmark = 0;
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it)
		if( !it->GetMarker(slice) && it->getBeg()->GetMarker(slice) && it->getEnd()->GetMarker(slice) )
		{
			if( it->Integer(material) != 2 ) std::cout << "Edge supposed to get material 2, but have " << it->Integer(material) << std::endl;
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
			mat[it->Integer(material)]++;
			tot++;
		}
		std::cout << "edge materials, 0: " << mat[0] << " 1: " << mat[1] << " 2: " << mat[2] << std::endl;
	}
	
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
	{
		int mat[3] = {0,0,0};
		ElementArray<Edge> edges = it->getEdges();
		for(int k = 0; k < (int)edges.size(); ++k)
			mat[edges[k]->Integer(material)]++;
		if( !(mat[0] == 0 || mat[1] == 0) )
		{
			it->Integer(material) = 2;
			std::cout << "oops, materials for face edges were not split, 0: " << mat[0] << " ,1: " << mat[1] << " ,2: " << mat[2] << std::endl;
		}
		else if( mat[0] != 0 ) it->Integer(material) = 0;
		else if( mat[1] != 0 ) it->Integer(material) = 1;
		else it->Integer(material) = 2;
	}
	
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
				for(int k = 0; k < edges.size(); ++k) if( !edges[k]->GetMarker(visited) )
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
			std::pair<Face,bool> f = m.CreateFace(order_edges);
			f.first->Integer(material) = 2;
			f.first->Bulk(sliced) = 1;
			ElementArray<Cell> ret = Cell::SplitCell(it->self(),ElementArray<Face>(&m,1,f.first.GetHandle()),0);
			nslice++;
		}
	}
	m.ReleaseMarker(visited);

	std::cout << "sliced cells: " << nslice << std::endl;
	
	if( !Element::CheckConnectivity(&m) )
		std::cout << "Connectivity is broken" << std::endl;
	
	{
		int mat[3] = {0,0,0};
		int tot = 0;
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
		{
			mat[it->Integer(material)]++;
			tot++;
		}
		std::cout << "face materials, 0: " << mat[0] << " 1: " << mat[1] << " 2: " << mat[2] << std::endl;
	}
	
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	{
		int mat[3] = {0,0,0};
		ElementArray<Face> faces = it->getFaces();
		for(int k = 0; k < (int)faces.size(); ++k)
			mat[faces[k]->Integer(material)]++;
		if( !(mat[0] == 0 || mat[1] == 0) )
		{
			it->Integer(material) = 2;
			std::cout << "oops, materials for cell " << it->LocalID() << " faces were not split, 0: " << mat[0] << " ,1: " << mat[1] << " ,2: " << mat[2] << " slice edges " << it->getEdges(slice).size() << std::endl;
		}
		else if( mat[0] != 0 ) it->Integer(material) = 0;
		else if( mat[1] != 0 ) it->Integer(material) = 1;
		else
		{
			//double cnt[3];
			//it->Centroid(cnt);
			//double v = func(cnt[0],cnt[1],cnt[2],type);
			//it->Integer(material) = (v <= 0.0 ? 0 : 1);
			it->Integer(material) = 2;
			std::cout << "oops cannot determine material for cell, all faces have type 2, set to " << it->Integer(material);
		}
	}
	
	
	{
		int mat[3] = {0,0,0};
		int tot = 0;
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			mat[it->Integer(material)]++;
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
		if( v < vv*0.15 )
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
	
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		if( it->Integer(material) == 0 )
			it->Delete();
	
	for(ElementType etype = FACE; etype >= NODE; etype = PrevElementType(etype) )
		for(Mesh::iteratorElement it = m.BeginElement(etype); it != m.EndElement(); ++it)
			if( it->nbAdjElements(CELL) == 0 ) it->Delete();
	

	m.ReleaseMarker(slice,NODE|EDGE);

	m.Save(grid_out);
	return 0;
}