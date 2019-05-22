#include "save_mesh.h"


std::string SaveMesh(INMOST::Mesh * m, std::string prefix, OutputMeshFormat format)
{
	std::string file = prefix;
	switch(format)
	{
		case VTK:
			if( m->GetProcessorsNumber() == 1 )
				prefix += ".vtk";
			else
				prefix += ".pvtk";
			break;
		case GMV:
			prefix += ".gmv";
			break;
		case PMF:
			prefix += ".pmf";
			break;
		case XML:
			prefix += ".xml";
			break;
	}
	m->Save(prefix);
	return prefix;
}

std::string SaveMesh(INMOST::Mesh * m, std::string prefix, int counter, OutputMeshFormat format)
{
	std::stringstream prefix_counter;
	prefix_counter << prefix << std::setw(5) << std::setfill('0') << counter;
	return SaveMesh(m,prefix_counter.str(),format);
}
