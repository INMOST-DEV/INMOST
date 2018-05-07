#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 2 ) 
	{
		printf("Usage: %s input_mesh [output_mesh]\n",argv[0]);
		return -1;
	}

	Mesh A("orig"),B("voro");
	A.Load(argv[1]);
	
	Mesh::GeomParam table;
	table[ORIENTATION] = FACE;
	table[BARYCENTER] = FACE | CELL | EDGE;
	A.PrepareGeometricData(table);
	
	
	//create nodes of the voronoi mesh
	Tag cell2node = A.CreateTag("CELL2NODE",DATA_REMOTE_REFERENCE,CELL|FACE|EDGE|NODE,FACE|EDGE|NODE,1);
	for(Mesh::iteratorCell it = A.BeginCell(); it != A.EndCell(); ++it)
	{
		real cnt[3];
		it->Barycenter(cnt);
		it->RemoteReference(cell2node) = RemoteHandleType(&B,B.CreateNode(cnt).GetHandle());
	}
	for(Mesh::iteratorFace it = A.BeginFace(); it != A.EndFace(); ++it) if( it->Boundary() )
	{
		real cnt[3];
		it->Barycenter(cnt);
		it->RemoteReference(cell2node) = RemoteHandleType(&B,B.CreateNode(cnt).GetHandle());
	}
	//add corner nodes
	int corners = 0;
	MarkerType corner = A.CreateMarker();
	for(Mesh::iteratorEdge it = A.BeginEdge(); it != A.EndEdge(); ++it) if( it->Boundary() )
	{
		real nrm[2][3], dot, cnt[3];
		int nbndfaces = 0;
		ElementArray<Face> faces = it->getFaces();
		for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end(); ++kt) if( kt->Boundary() )
			kt->UnitNormal(nrm[nbndfaces++]);
		assert(nbndfaces==2);
		dot = 0;
		for(int k = 0; k < 3; ++k) dot += nrm[0][k]*nrm[1][k];
		if( dot < 0.95 )
		{
			corners++;
			it->SetMarker(corner);
			it->Barycenter(cnt);
			it->RemoteReference(cell2node) = RemoteHandleType(&B,B.CreateNode(cnt).GetHandle());
		}
	}
	for(Mesh::iteratorNode it = A.BeginNode(); it != A.EndNode(); ++it) if( it->Boundary() )
	{
		real cnt[3];
		int ncorners = it->nbAdjElements(EDGE,corner);
		if( ncorners > 2 )
		{
			corners++;
			it->SetMarker(corner);
			it->Barycenter(cnt);
			it->RemoteReference(cell2node) = RemoteHandleType(&B,B.CreateNode(cnt).GetHandle());
		}
		else if( ncorners == 2 )
		{
			ElementArray<Edge> edges = it->getEdges(corner);
			Node n1 = edges[0]->getBeg() == it->self() ? edges[0]->getEnd() : edges[0]->getBeg();
			Node n2 = edges[1]->getBeg() == it->self() ? edges[1]->getEnd() : edges[1]->getBeg();
			real cnt1[3], cnt2[3], r1[3],r2[3], l, dot;
			it->Barycenter(cnt);
			n1->Barycenter(cnt1);
			n2->Barycenter(cnt2);
			r1[0] = cnt[0] - cnt1[0];
			r1[1] = cnt[1] - cnt1[1];
			r1[2] = cnt[2] - cnt1[2];
			l = sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
			r1[0] /= l;
			r1[1] /= l;
			r1[2] /= l;
			r2[0] = cnt2[0] - cnt1[0];
			r2[1] = cnt2[1] - cnt1[1];
			r2[2] = cnt2[2] - cnt1[2];
			l = sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
			r2[0] /= l;
			r2[1] /= l;
			r2[2] /= l;
			dot = 0;
			for(int k = 0; k < 3; ++k) dot += r1[k]*r2[k];
			if( dot < 0.95 )
			{
				corners++;
				it->SetMarker(corner);
				it->Barycenter(cnt);
				it->RemoteReference(cell2node) = RemoteHandleType(&B,B.CreateNode(cnt).GetHandle());
			}
		}
	}
	std::cout << "corner nodes: " << corners << std::endl;
	//create edges of the voronoi mesh
	Tag face2edge = A.CreateTag("FACE2EDGE",DATA_REMOTE_REFERENCE,FACE|EDGE,EDGE,1);
	Tag face2edge_corner = A.CreateTag("FACE2EDGE_CORNER",DATA_REMOTE_REFERENCE,EDGE,EDGE,2);
	Tag face2edge_corner_node = A.CreateTag("FACE2EDGE_CORNER_NODE",DATA_REMOTE_REFERENCE,NODE,NODE);
	ElementArray<Node> edgenodes(&B,2);
	for(Mesh::iteratorFace it = A.BeginFace(); it != A.EndFace(); ++it)
	{
		Cell c1 = it->BackCell();
		Cell c2 = it->FrontCell();
		if( c2.isValid() )
		{
			edgenodes[0] = MakeElement(c1.RemoteReference(cell2node)).getAsNode();
			edgenodes[1] = MakeElement(c2.RemoteReference(cell2node)).getAsNode();
		}
		else
		{
			edgenodes[0] = MakeElement(c1.RemoteReference(cell2node)).getAsNode();
			edgenodes[1] = MakeElement(it->RemoteReference(cell2node)).getAsNode();
		}
		it->RemoteReference(face2edge) = RemoteHandleType(&B,B.CreateEdge(edgenodes).first.GetHandle());
	}
	int process_corners = 0;
	for(Mesh::iteratorEdge it = A.BeginEdge(); it != A.EndEdge(); ++it) if( it->Boundary() )
	{
		ElementArray<Face> faces = it->getFaces();
		ElementArray<Face> bndfaces(&A);
		for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end(); ++kt) if( kt->Boundary() )
			bndfaces.push_back(kt->self());
		assert(bndfaces.size() == 2);
		
		if( it->GetMarker(corner) )
		{
			process_corners++;
			edgenodes[0] = MakeElement(bndfaces[0].RemoteReference(cell2node)).getAsNode();
			edgenodes[1] = MakeElement(it->RemoteReference(cell2node)).getAsNode();
			
			it->RemoteReferenceArray(face2edge_corner).at(0) = RemoteHandleType(&B,B.CreateEdge(edgenodes).first.GetHandle());
			
			edgenodes[0] = MakeElement(bndfaces[1].RemoteReference(cell2node)).getAsNode();
			edgenodes[1] = MakeElement(it->RemoteReference(cell2node)).getAsNode();
			
			it->RemoteReferenceArray(face2edge_corner).at(1) = RemoteHandleType(&B,B.CreateEdge(edgenodes).first.GetHandle());
		}
		else
		{
			edgenodes[0] = MakeElement(bndfaces[0].RemoteReference(cell2node)).getAsNode();
			edgenodes[1] = MakeElement(bndfaces[1].RemoteReference(cell2node)).getAsNode();
		
			it->RemoteReference(face2edge) = RemoteHandleType(&B,B.CreateEdge(edgenodes).first.GetHandle());
		}
	}
	//add edges connecting corners
	for(Mesh::iteratorNode it = A.BeginNode(); it != A.EndNode(); ++it) if( it->Boundary() )
	{
		ElementArray<Edge> edges = it->getEdges(corner);
		if( it->GetMarker(corner) ) //connect each corner node
		{
			edgenodes[0] = MakeElement(it->RemoteReference(cell2node)).getAsNode();
			for(ElementArray<Edge>::iterator kt = edges.begin(); kt != edges.end(); ++kt)
			{
				process_corners++;
				edgenodes[1] = MakeElement(kt->RemoteReference(cell2node)).getAsNode();
				it->RemoteReferenceArray(face2edge_corner_node).push_back(&B,B.CreateEdge(edgenodes).first.GetHandle());
			}
		}
		else if( edges.size() == 2 ) //connect two corner nodes
		{
			process_corners++;
			edgenodes[0] = MakeElement(edges[0]->RemoteReference(cell2node)).getAsNode();
			edgenodes[1] = MakeElement(edges[1]->RemoteReference(cell2node)).getAsNode();
			
			it->RemoteReferenceArray(face2edge_corner_node).push_back(&B,B.CreateEdge(edgenodes).first.GetHandle());
		}
	}
	std::cout << "corners processed: " << process_corners << std::endl;
	//create faces of the voronoi mesh
	Tag edge2face = A.CreateTag("EDGE2FACE",DATA_REMOTE_REFERENCE,EDGE|NODE,NODE,1);
	ElementArray<Edge> faceedges(&B);
	process_corners = 0;
	for(Mesh::iteratorEdge it = A.BeginEdge(); it != A.EndEdge(); ++it)
	{
		ElementArray<Face> faces = it->getFaces();
		ElementArray<Node> nodes = it->getNodes();
		for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end(); ++kt)
			faceedges.push_back(MakeElement(kt->RemoteReference(face2edge)));
		
		if( it->GetMarker(corner) )
		{
			process_corners++;
			faceedges.push_back(it->RemoteReferenceArray(face2edge_corner)[0]);
			faceedges.push_back(it->RemoteReferenceArray(face2edge_corner)[1]);
		}
		else if( it->Boundary() )
			faceedges.push_back(MakeElement(it->RemoteReference(face2edge)));
		
		Face f = B.CreateFace(faceedges).first;
		f.FixEdgeOrder();
		it->RemoteReference(edge2face) = RemoteHandleType(&B,B.CreateFace(faceedges).first.GetHandle());
		faceedges.clear();
	}
	//construct boundary faces without corner nodes
	for(Mesh::iteratorNode it = A.BeginNode(); it != A.EndNode(); ++it) if( it->Boundary() )
	{
		ElementArray<Edge> edges = it->getEdges();
		for(ElementArray<Edge>::iterator kt = edges.begin(); kt != edges.end(); ++kt) if( kt->Boundary() )
		{
			if( kt->GetMarker(corner) )
			{
				process_corners++;
				faceedges.push_back(kt->RemoteReferenceArray(face2edge_corner)[0]);
				faceedges.push_back(kt->RemoteReferenceArray(face2edge_corner)[1]);
			}
			else
				faceedges.push_back(MakeElement(kt->RemoteReference(face2edge)));
		}
		Face f = B.CreateFace(faceedges).first;
		f.FixEdgeOrder();
		it->RemoteReference(edge2face) = RemoteHandleType(&B,B.CreateFace(faceedges).first.GetHandle());
		faceedges.clear();
	}
	std::cout << "corners processed: " << process_corners << std::endl;
	//construct cells
	ElementArray<Face> cellfaces(&B);
	for(Mesh::iteratorNode it = A.BeginNode(); it != A.EndNode(); ++it)
	{
		ElementArray<Edge> edges = it->getEdges();
		for(ElementArray<Edge>::iterator kt = edges.begin(); kt != edges.end(); ++kt)
			cellfaces.push_back(MakeElement(kt->RemoteReference(edge2face)));
		if( it->Boundary() )
			cellfaces.push_back(MakeElement(it->RemoteReference(edge2face)));
		B.CreateCell(cellfaces);
		cellfaces.clear();
	}
	
	//split corner faces by edges connected to corner nodes
	process_corners= 0;
	for(Mesh::iteratorNode it = A.BeginNode(); it != A.EndNode(); ++it) if( it->HaveData(face2edge_corner_node) )
	{
		process_corners++;
		Storage::remote_reference_array arr = it->RemoteReferenceArray(face2edge_corner_node);
		ElementArray<Edge> edges(&B);
		for(Storage::remote_reference_array::iterator kt = arr.begin(); kt != arr.end(); ++kt)
			edges.push_back(kt->GetHandle());
		Face f = MakeElement(it->RemoteReference(edge2face)).getAsFace();
		Face::SplitFace(f,edges,0);
	}
	std::cout << "corners splitted: " << process_corners << std::endl;

	if( A.HaveTag("GRIDNAME") )
	{
		Storage::bulk_array nameA = A.self().BulkArray(A.GetTag("GRIDNAME"));
		Storage::bulk_array nameB = B.self().BulkArray(B.CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
		std::string ins = "_dual";
		nameB.insert(nameB.end(),nameA.begin(),nameA.end());
		nameB.insert(nameB.end(),ins.c_str(),ins.c_str()+ins.size());
	}
	
	if( argc > 2 )
	{
		std::cout << "Save to " << argv[2] << std::endl;
		B.Save(argv[2]);
	}
	else
	{
		std::cout << "Save to out.vtk" << std::endl;
		B.Save("out.vtk");
	}

	return 0;
}
