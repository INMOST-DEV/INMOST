#include "inmost.h"

using namespace INMOST;



int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}

	std::string grid_out = "grid.vtk";
	
	if (argc > 2) grid_out = std::string(argv[2]);


	Mesh m;
	m.Load(argv[1]);
	
	Mesh::GeomParam table;
	table[MEASURE] = CELL;
	m.PrepareGeometricData(table);
	

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;
	double tt = Timer();
	
	MarkerType collapse = m.CreateMarker();
	
	int ncells = 0;
	
#pragma omp parallel for reduction(+:ncells)
	for(int k = 0; k < m.CellLastLocalID(); ++k) if( m.isValidCell(k) )
	{
		Cell c = m.CellByLocalID(k);
		ElementArray<Cell> around = c.BridgeAdjacencies2Cell(FACE);
		double vol0 = c.Volume(), vol = vol0;
		for(unsigned l = 0; l < around.size(); ++l)
			vol += around[l].Volume();
		if( vol0 / vol < 0.03 )
		{
			c.SetMarker(collapse);
			ncells++;
		}
	}
	std::cout << "collapse cells: " << ncells << std::endl;
	
	int nedges = 0;
	
#pragma omp parallel for reduction(+:nedges)	
	for(int k = 0; k < m.EdgeLastLocalID(); ++k) if( m.isValidEdge(k) )
	{
		Edge e = m.EdgeByLocalID(k);
		ElementArray<Cell> around = e.getCells();
		for(unsigned l = 0; l < around.size(); ++l)
			if( around[l].GetMarker(collapse) )
			{
				e.SetMarker(collapse);
				nedges++;
				break;
			}
	}
	
	std::cout << "collapse edges: " << nedges << std::endl;
	
	int ncollapsed = 0;
	for(int k = 0; k < m.EdgeLastLocalID(); ++k) if( m.isValidEdge(k) )
	{
		Edge e = m.EdgeByLocalID(k);
		if( e.GetMarker(collapse) ) 
		{
			Edge::CollapseEdge(e,0);
			ncollapsed++;
			if( ncollapsed % 100 == 0 )
				m.Save("collapsed.vtk");
			if( !Element::CheckConnectivity(&m) )
			{
				std::cout << "connectivity is bad" << std::endl;
				throw -1;
			}
		}
	}
	
	
	
	std::cout << "Time to fix tiny:" << Timer() - tt << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	m.Save(grid_out);
	return 0;
}
