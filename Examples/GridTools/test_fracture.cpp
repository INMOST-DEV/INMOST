#include "inmost.h"
#include "fracture.h"
using namespace INMOST;

// cell-> node, edge, face
// face-> node, edge
// edge-> node
// based on eigenvalues in bounding ellipse


int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [fill_fracture=true] [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}

	std::string grid_out = "grid.vtk";
	
	bool fill_fracture = true;
	if( argc > 2 ) fill_fracture = atoi(argv[2]) ? true : false;
	if (argc > 3) grid_out = std::string(argv[3]);


	Mesh m;
	m.Load(argv[1]);
	
	//Mesh::GeomParam table;
	//table[MEASURE] = CELL;
	//m.PrepareGeometricData(table);
	
	TagReal aperture = m.CreateTag("aperture",DATA_REAL,FACE,FACE,1);
	TagBulk sliced;
	if( m.HaveTag("SLICED") )
		sliced = m.GetTag("SLICED");
	else
	{
		std::cout << "No tag named SLICED on the grid" << std::endl;
		return -1;
	}
	
	for(int k = 0; k < m.FaceLastLocalID(); ++k) if( m.isValidFace(k) )
	{
		Face f = m.FaceByLocalID(k);
		if( f.HaveData(sliced) )
			aperture[f] = 0.001;
	}
	

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;
	double tt = Timer();
	
	Fracture f(m);
	f.Open(aperture,fill_fracture);
	
	
	std::cout << "Time to open fracture:" << Timer() - tt << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	m.Save(grid_out);
	return 0;
}
