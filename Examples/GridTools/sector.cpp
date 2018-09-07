#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 2 ) 
	{
		printf("Usage: %s input_mesh [rotation_angle=0 degrees] [refine_boundary=1] [output_mesh]\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.Load(argv[1]);
	
	const double pi = 3.1415926536;
	double theta = 0, ct, st;
	int refine = 1;
	if( argc > 2 ) theta = atof(argv[2])/180.0*pi;
	if( argc > 3 ) refine = atoi(argv[3]);
	ct = cos(theta);
	st = sin(theta);
	
	double xmin = 1.0e20, xmax = -1.0e20;
	double ymin = 1.0e20, ymax = -1.0e20;
	
	for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
	{
		double x = n->Coords()[0];
		double y = n->Coords()[1];
		if( x > xmax ) xmax = x;
		if( x < xmin ) xmin = x;
		if( y > ymax ) ymax = y;
		if( y < ymin ) ymin = y;
	}
	
	std::cout << "x: " << xmin << ":" << xmax << std::endl;
	std::cout << "y: " << ymin << ":" << ymax << std::endl;
	
	if( refine )
	{
		for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
		{
			double a = (n->Coords()[0]-xmin)/(xmax-xmin);
			a = sin(3.14159265359*a/2.0);
			n->Coords()[0] = (xmax-xmin)*a + xmin;
		}
	}
	
	for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
	{
		double x = n->Coords()[0], mx = x;
		double y = n->Coords()[1], my = y;
		double alpha = pi/4*(2*y-1);
		mx += 1 - x*(1-cos(alpha)) + x*cos(2*alpha)*0.2;
		my += x*sin(alpha);
		n->Coords()[0] = (mx-0.5)*ct - (my-0.5)*st + 0.5;
		n->Coords()[1] = (mx-0.5)*st + (my-0.5)*ct + 0.5;
	}
	
	xmin = 1.0e20, xmax = -1.0e20;
	ymin = 1.0e20, ymax = -1.0e20;
	
	for(Mesh::iteratorNode n = m.BeginNode(); n != m.EndNode(); ++n)
	{
		double x = n->Coords()[0];
		double y = n->Coords()[1];
		if( x > xmax ) xmax = x;
		if( x < xmin ) xmin = x;
		if( y > ymax ) ymax = y;
		if( y < ymin ) ymin = y;
	}
	
	std::cout << "x: " << xmin << ":" << xmax << std::endl;
	std::cout << "y: " << ymin << ":" << ymax << std::endl;
	
	
	
	if( argc > 4 )
	{
		std::cout << "Save to " << argv[4] << std::endl;
		m.Save(argv[4]);
	}
	else
	{
		std::cout << "Save to out.vtk" << std::endl;
		m.Save("out.vtk");
	}

	return 0;
}
