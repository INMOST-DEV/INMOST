#include "inmost.h"

using namespace INMOST;

//todo: want to separate all faces into those that share edges
//see todo: inside text


int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [tag=MATERIAL] [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}
	
	Mesh::Initialize(&argc,&argv);
	
	std::string tag_name = "MATERIAL";
	if( argc > 2) tag_name = std::string(argv[2]);
	
	std::string grid_out = "grid.vtk";
	if (argc > 3) grid_out = std::string(argv[3]);


	Mesh m;
	m.SetCommunicator(INMOST_MPI_COMM_WORLD);
	m.Load(argv[1]);
	
	if( !m.HaveTag(tag_name) )
	{
		std::cout << "mesh does not have tag named " << tag_name << std::endl;
		exit(-1);
	}
	
	TagInteger mat = m.GetTag(tag_name);

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	std::map< int,ElementArray<Cell> > mat_cells;
	
	for(int q = 0; q < m.CellLastLocalID(); ++q) if( m.isValidCell(q) )
	{
		Cell n = m.CellByLocalID(q);
		mat_cells[ mat[n] ].push_back(n);
	}
	
	//Unite cells
	//todo: split cells that are not connected by face
	std::cout << "Unite cells" << std::endl;
	int nunited = 0;
	for(std::map< int,ElementArray<Cell> >::iterator it = mat_cells.begin();
		it != mat_cells.end(); ++it)
	{
		if( it->second.size() > 1 )
		{
			mat[Cell::UniteCells(it->second,0)] = it->first;
			nunited++;
		}
	}
	std::cout << "united: " << nunited << std::endl;
	/*
	std::cout << "Unite faces" << std::endl;
	
	std::map< std::pair<int,int>, ElementArray<Face> > faces;
	
	nunited = 0;
	
	
	
	for(int q = 0; q < m.FaceLastLocalID(); ++q) if( m.isValidFace(q) )
	{
		Face n = m.FaceByLocalID(q);
		std::pair<int,int> f_cells;
		int bc = n.BackCell().LocalID();;
		int fc = n.FrontCell().isValid() ? n.FrontCell().LocalID() : -1;
		f_cells.first = std::min(bc,fc);
		f_cells.second = std::max(bc,fc);
		faces[f_cells].push_back(n);
	}
	
	std::cout << "computed faces" << std::endl;
	
	m.BeginModification();
	
	 
	 //todo: split faces that are not connected by edge
	for(std::map< std::pair<int,int>, ElementArray<Face> >::iterator it = faces.begin();
		it != faces.end(); ++it)
	{
		if( it->second.size() > 1 )
		{
			std::cout << it->first.first << "<->" << it->first.second << " faces " << it->second.size() << std::endl;
			Face::UniteFaces(it->second,0);
			nunited++;
		}
	}
	
	m.EndModification();
	
	std::cout << "united: " << nunited << std::endl;
	*/
	
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;


	m.Save(grid_out);
	
	Mesh::Finalize();
	return 0;
}
