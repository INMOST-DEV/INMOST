#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 2 ) 
	{
		printf("Usage: %s input_mesh1 [input_mesh2] ... [output_mesh]\n",argv[0]);
		return -1;
	}

	Mesh m;
	for(int k = 1 ; k < argc; ++k)
	{
		std::cout << "load " << argv[k] << std::endl;
		m.Load(argv[k]);
	}
	
	double xmax = -1.0e20, xmin = 1.0e20;
	double ymax = -1.0e20, ymin = 1.0e20;
	double zmax = -1.0e20, zmin = 1.0e20;
	
	for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
	{
		double x = n->Coords()[0], y = n->Coords()[1], z = n->Coords()[2];
		if( x > xmax ) xmax = x;
		if( x < xmin ) xmin = x;
		if( y > ymax ) ymax = y;
		if( y < ymin ) ymin = y;
		if( z > zmax ) zmax = z;
		if( z < zmin ) zmin = z;
	}
	
	std::cout << "x: " << xmin << ":" << xmax << std::endl;
	std::cout << "y: " << ymin << ":" << ymax << std::endl;
	std::cout << "z: " << zmin << ":" << zmax << std::endl;
	
	std::cout << "Save to out.pmf" << std::endl;
	m.Save("out.pmf");
	
	return 0;
}
