#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


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
	  
	std::fstream map(mapfile.c_str(),std::ios::in);
	double xyz[8][3] =  {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
	for(int k = 0; k < (1<<dims); ++k)
		for(int d = 0; d < dims; ++d)
			map >> xyz[k][d];
	map.close();
	
	for(int k = (1<<dims); k < 8; k+=(1<<dims))
	{
		for(int q = 0; q < (1<<dims); ++q)
			for(int d = 0; d < dims; ++d)
				xyz[k+q][d] = xyz[k-(1<<dims)+q][d];
	}
	
	for(int k = 0; k < 8; ++k)
		printf("%g %g %g\n",xyz[k][0],xyz[k][1],xyz[k][2]);
	
	Mesh mesh; 
	mesh.Load(infile);
	
	double cmin[3], cmax[3];
	for(int d = 0; d < 3; ++d)
	{
		cmin[d]=1.0e+20;
		cmax[d]=-1.0e+20;
	}
	
	for(Mesh::iteratorNode n = mesh.BeginNode(); n != mesh.EndNode(); ++n)
	{
		for(int d = 0; d < 3; ++d)
		{
			if( cmin[d] > n->Coords()[d] ) cmin[d] = n->Coords()[d];
			if( cmax[d] < n->Coords()[d] ) cmax[d] = n->Coords()[d];
		}
	}
	printf("before\n");
	for(int d = 0; d < 3; ++d)
		printf("%d %g:%g\n",d,cmin[d],cmax[d]);
		
	  
	for(Mesh::iteratorNode n = mesh.BeginNode(); n != mesh.EndNode(); ++n)
	{
		double c[3] = {0,0,0};
		for(int d = 0; d < 3; ++d)
			c[d] = (n->Coords()[d]-cmin[d])/(cmax[d]-cmin[d]);
		for(int d = 0; d < dims; ++d)
			n->Coords()[d] = ((xyz[0][d]*(1-c[0]) + xyz[1][d]*c[0])*(1-c[1])+
			                  (xyz[2][d]*(1-c[0]) + xyz[3][d]*c[0])*c[1])*(1-c[2])+
			                 ((xyz[4][d]*(1-c[0]) + xyz[5][d]*c[0])*(1-c[1])+
			                  (xyz[6][d]*(1-c[0]) + xyz[7][d]*c[0])*c[1])*c[2];
	} 
	
	for(int d = 0; d < 3; ++d)
	{
		cmin[d]=1.0e+20;
		cmax[d]=-1.0e+20;
	}
	
	for(Mesh::iteratorNode n = mesh.BeginNode(); n != mesh.EndNode(); ++n)
	{
		for(int d = 0; d < 3; ++d)
		{
			if( cmin[d] > n->Coords()[d] ) cmin[d] = n->Coords()[d];
			if( cmax[d] < n->Coords()[d] ) cmax[d] = n->Coords()[d];
		}
	}
	
	printf("after\n");
	for(int d = 0; d < 3; ++d)
		printf("%d %g:%g\n",d,cmin[d],cmax[d]);
		
	
	Storage::bulk_array name = mesh->self()->BulkArray(mesh.CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
	name.insert(name.end(),mesh_name.begin(),mesh_name.end());
	
	  
	printf("I'm ready!\n");
	mesh.Save(outfile);
	  
	
	return 0;
}
