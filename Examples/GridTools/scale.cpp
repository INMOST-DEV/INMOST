#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 2 ) 
	{
		printf("Usage: %s input_mesh [scale_x=1] [scale_y=1] [scale_z=1] [output_mesh]\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.Load(argv[1]);
	double sx = argc > 2 ? atof(argv[2]) : 1;
	double sy = argc > 3 ? atof(argv[3]) : 1;
	double sz = argc > 4 ? atof(argv[4]) : 1;
	std::string save_file = argc > 5 ? std::string(argv[5]) : "out.vtk";
	
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
	
	for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
	{
		n->Coords()[0] *= sx;
		n->Coords()[1] *= sy;
		n->Coords()[2] *= sz;
	}
	
	xmax = -1.0e20, xmin = 1.0e20;
	ymax = -1.0e20, ymin = 1.0e20;
	zmax = -1.0e20, zmin = 1.0e20;
	
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
	
	std::cout << "Save to " << save_file << std::endl;
	m.Save(save_file);
	
	return 0;
}
