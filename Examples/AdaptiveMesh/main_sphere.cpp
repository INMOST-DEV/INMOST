#include "amesh.h"

using namespace INMOST;

int main(int argc, char ** argv)
{
	Mesh::Initialize(&argc,&argv);
	
	if( argc > 1 )
	{
		Mesh m;
		m.Load(argv[1]);
		TagInteger indicator = m.CreateTag("INDICATOR",DATA_INTEGER,CELL,NONE,1);
		AdaptiveMesh am(m);
		
		int max_levels = 2;
		if( argc > 2 ) max_levels = atoi(argv[2]);
		
		int numref;
		do
		{
			numref = 0;
			for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
				if( am.GetLevel(it->self()) < max_levels )
				{
					INMOST_DATA_REAL_TYPE x[3];
					it->Centroid(x);
					if( sqrt((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5)+(x[2]-0.5)*(x[2]-0.5)) < 0.25 )
					{
						indicator[it->self()] = 1;
						numref++;
					}
				}
			if( numref )
			{
				if( !am.Refine(indicator) ) break;
				for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) indicator[it->self()] = 0;
			}
		}
		while(numref);
		std::string file = "out.vtk";
		if( argc > 3 ) file = std::string(argv[3]);
		m.Save(file);
	}
	else std::cout << "Usage: " << argv[0] << " mesh_file [max_levels=2] [mesh_out=out.vtk]" << std::endl;
}
