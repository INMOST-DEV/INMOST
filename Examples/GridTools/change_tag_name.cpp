#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 4 )
	{
		printf("Usage: %s input_mesh tag_name new_tag_name [output_mesh=out.pmf]\n",argv[0]);
		return -1;
	}
	
	std::string inpmesh;
	std::string outmesh = "out.pmf";
	std::string tagname;
	std::string newtagname;
	if( argc > 1 ) inpmesh = std::string(argv[1]);
	if( argc > 2 ) tagname = std::string(argv[2]);
	if( argc > 3 ) newtagname = std::string(argv[3]);
	if( argc > 4 ) outmesh = std::string(argv[4]);
	Mesh m;
	std::cout << "load " << inpmesh << std::endl;
	m.Load(inpmesh);
	std::cout << "rename " << tagname << " -> " << newtagname << std::endl;
	if( m.RenameTag(tagname,newtagname) )
		std::cout << "success!" << std::endl;
	else
		std::cout << "failure" << std::endl;
	std::cout << "save " << outmesh << std::endl;
	m.Save(outmesh);
	
	
	return 0;
}
