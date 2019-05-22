#include "amesh.h"

using namespace INMOST;

int main(int argc, char ** argv)
{
	Mesh::Initialize(&argc,&argv);
	
	if( argc > 1 )
	{
		Mesh mm;
		mm.Load(argv[1]);
		TagInteger indicator = mm.CreateTag("INDICATOR",DATA_INTEGER,CELL,NONE,1);
		AdaptiveMesh m(mm);
		
		int max_levels = 2;
		if( argc > 2 ) max_levels = atoi(argv[2]);
		
		int numref;
		do
		{
			numref = 0;
			for(Mesh::iteratorCell it = mm.BeginCell(); it != mm.EndCell(); ++it)
				if( m.GetLevel(it->self()) < max_levels )
				{
					double x[3];
					it->Centroid(x);
					if( x[0] > 0.3 && x[0] < 0.7 && x[1] > 0.3 && x[1] < 0.7 && x[2] > 0.3 && x[2] < 0.7)
					{
						indicator[it->self()] = 1;
						numref++;
					}
				}
			if( numref )
			{
				if( !m.Refine(indicator) ) break;
				for(Mesh::iteratorCell it = mm.BeginCell(); it != mm.EndCell(); ++it) indicator[it->self()] = 0;
			}
		}
		while(numref);
		std::string file = "out.vtk";
		if( argc > 3 ) file = std::string(argv[3]);
		mm.Save(file);
	}
	else std::cout << "Usage: " << argv[0] << " mesh_file [max_levels=2] [mesh_out=out.vtk]" << std::endl;
}
