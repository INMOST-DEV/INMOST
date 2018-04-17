#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		printf("Usage: %s input_mesh1 [input_mesh2] ... output_format\n",argv[0]);
		return -1;
	}
	
	std::string output_fmt = argv[argc-1];

	std::cout << "Output format: " << output_fmt << std::endl;
	
	
	for(int k = 1 ; k < argc-1; ++k)
	{
		Mesh m;
		std::string load_file = std::string(argv[k]);
		std::string save_file = load_file.substr(0,load_file.find_last_of(".")) + "." + output_fmt;
		
		std::cout << load_file << " -> " << save_file << std::endl;
		m.Load(load_file);
		m.Save(save_file);
	}
	
	return 0;
}
