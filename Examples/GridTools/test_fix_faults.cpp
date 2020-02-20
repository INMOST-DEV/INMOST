#include "inmost.h"
#include "fix_faults.h"


int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}

	std::string grid_out = "grid.vtk";
	if (argc > 2) grid_out = std::string(argv[2]);


	{
		Mesh m;
		m.Load(argv[1]);
		FixFaults fix(m);	
		fix.FixMeshFaults();
		m.Save(grid_out);
	}
	return 0;
}
