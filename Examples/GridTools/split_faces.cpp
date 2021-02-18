#include "inmost.h"
#include <stdio.h>
#include <deque>
//#include "tetgen/tetgen.h"
using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;



int main(int argc, char ** argv)
{
	if( argc < 2 ) 
	{
		std::cout << "Usage: " << argv[0] << " input_mesh [output_mesh]\n";
		return -1;
	}
	
	//~ int special_triangle_split = 1;
	//~ if( argc > 3 ) special_triangle_split = atoi(argv[3]);

	Mesh m;
	m.SetFileOption("ECL_CURVILINEAR","FALSE");
	m.SetFileOption("ECL_SPLIT_GLUED","TRUE");
	m.Load(argv[1]);
	double cnt[3];
	ElementArray<Node> edge_nodes(&m,2);
	MarkerType new_node = m.CreateMarker();
	m.BeginModification();
	for(Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it) if( !it->New() )
	{
		for(int k = 0; k < 3; ++k)
			cnt[k] = (it->getBeg().Coords()[k]+it->getEnd().Coords()[k])*0.5;
		Node n = m.CreateNode(cnt);
		n.SetMarker(new_node);
		Edge::SplitEdge(it->self(),ElementArray<Node>(&m,1,n.GetHandle()),0);
	}
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) if( !it->New() )
	{
		ElementArray<Node> new_nodes = it->getNodes(new_node);
		ElementArray<Edge> new_edges(&m,new_nodes.size());
		if( new_nodes.size() == 3 )
		{
			for(int k = 0; k < 3; ++k)
			{
				edge_nodes[0] = new_nodes[k];
				edge_nodes[1] = new_nodes[(k+1)%3];
				new_edges[k] = m.CreateEdge(edge_nodes).first;
			}
		}
		else
		{
			cnt[0] = cnt[1] = cnt[2] = 0;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < new_nodes.size(); ++k)
			{
				cnt[0] += new_nodes[k].Coords()[0];
				cnt[1] += new_nodes[k].Coords()[1];
				cnt[2] += new_nodes[k].Coords()[2];
			}
			cnt[0] /= new_nodes.size();
			cnt[1] /= new_nodes.size();
			cnt[2] /= new_nodes.size();
			Node n = m.CreateNode(cnt);
			edge_nodes[0] = n;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < new_nodes.size(); ++k)
			{
				edge_nodes[1] = new_nodes[k];
				new_edges[k] = m.CreateEdge(edge_nodes).first;
			}
		}
		Face::SplitFace(it->self(),new_edges,0);
	}
	
	m.EndModification();
	m.ReleaseMarker(new_node,NODE);
							  
							  
	if( argc > 2 )
		m.Save(argv[2]);
	else
		m.Save("out.vtk");

	return 0;
}
