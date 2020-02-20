#include "inmost.h"
#include "cube_transform.h"

std::string mesh_name = "_TRANSFORMED";

using namespace INMOST;




int main(int argc, char *argv[]) 
{
	
	if( argc < 3 )
	{
		printf("Usage: %s mapping_file mesh_file [output=out.pmf] [dims=3]\n",argv[0]);
		printf("mapping file contains 8 (for 3D) 4 (for 2D) 2 (for 1D) coordinates mapping each corner of the cube:\n");
		for(int k = 0; k < 8; ++k)
			printf("x%d y%d z%d\n",k,k,k);
		printf("      (6)----------(7)\n");
		printf("      /|           /|\n");
		printf("     / |          / |\n");
		printf("    /  |         /  |\n");
		printf("  (4)----------(5)  |\n");
		printf("   |   |        |   |\n");
		printf("   |   |        |   |\n");
		printf("   |  (2)-------|--(3)\n");
		printf("   |  /         |  /\n");
		printf("   | /          | /\n");
		printf("   |/           |/\n");
		printf("z (0)----------(1)\n");
		printf("^ y\n");
		printf("|/\n");
		printf("*-->x\n");
		return -1;
	}
	std::string mapfile(argv[1]);
	std::string infile(argv[2]);
	std::string outfile = "out.pmf";
	int dims = 3;
	if( argc > 3 ) outfile = std::string(argv[3]);
	if( argc > 4 ) dims = atoi(argv[4]);
	  
	
	Mesh mesh; 
	mesh.Load(infile);
	CubeTransform ct(mesh);
	
	ct.LoadMapping(mapfile);
	ct.PrintMapping();
	ct.Transform();
	
	Storage::bulk_array name = mesh->self()->BulkArray(mesh.CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
	name.insert(name.end(),mesh_name.begin(),mesh_name.end());
	
	  
	printf("I'm ready!\n");
	mesh.Save(outfile);
	return 0;
}
