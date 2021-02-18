#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		printf("Usage: %s input_mesh edge_number [output_mesh]\n",argv[0]);
		return -1;
	}
	
	Mesh A("");
	A.Load(argv[1]);
	
	Mesh::GeomParam table;
	table[ORIENTATION] = FACE;
	table[BARYCENTER] = FACE | CELL | EDGE;
	A.PrepareGeometricData(table);
	
	int edge_num = atoi(argv[2]);
	
	if( edge_num >= A.EdgeLastLocalID() )
	{
		std::cout << "input edge is " << edge_num << " but only ";
		std::cout << A.EdgeLastLocalID() << " edges in the mesh" << std::endl;;
		return -1;
	}

	Edge::CollapseEdge(Edge(&A,ComposeHandle(EDGE,edge_num)));
	
	if( A.HaveTag("GRIDNAME") )
	{
		Storage::bulk_array nameA = A.self().BulkArray(A.GetTag("GRIDNAME"));
		std::string ins = "_collapse_degenerate";
		nameA.insert(nameA.end(),ins.c_str(),ins.c_str()+ins.size());
	}
	
	if( argc > 3 )
	{
		std::cout << "Save to " << argv[3] << std::endl;
		A.Save(argv[3]);
	}
	else
	{
		std::cout << "Save to out.vtk" << std::endl;
		A.Save("out.vtk");
	}

	return 0;
}
