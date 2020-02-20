#include "inmost.h"
#include "set_layers.h"

std::string mesh_name = "_TRANSFORMED";

using namespace INMOST;




int main(int argc, char *argv[]) 
{
	
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file [output=out.pmf] [nlayers=10] [layer_coef=0.1] [height_coef=1]\n",argv[0]);
		return -1;
	}
	
	std::string infile(argv[1]);
	std::string outfile = "out.pmf";
	int nlayers = 10;
	double layer_coef = 0.1;
	double height_coef = 1;
	if( argc > 2 ) outfile = std::string(argv[2]);
	if( argc > 3 ) nlayers = atoi(argv[3]);
	if( argc > 4 ) layer_coef = atof(argv[4]);
	if( argc > 5 ) height_coef = atof(argv[5]);
	  
	
	Mesh mesh; 
	mesh.Load(infile);
	SetLayers sl(mesh);
	
	sl.RandomLayersZ(nlayers);
	sl.DeformLayers(layer_coef);
	sl.DeformHeight(height_coef);
	
	Storage::bulk_array name = mesh->self()->BulkArray(mesh.CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
	name.insert(name.end(),mesh_name.begin(),mesh_name.end());
	  
	printf("I'm ready!\n");
	mesh.Save(outfile);
	return 0;
}
