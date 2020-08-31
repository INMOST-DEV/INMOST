#include "inmost.h"
#include "cube_transform.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


using namespace INMOST;

void CubeTransform::PrintMapping()
{
	for(int k = 0; k < 8; ++k)
		printf("%g %g %g\n",xyz[k][0],xyz[k][1],xyz[k][2]);
}


void CubeTransform::LoadMapping(std::string mapfile, int dims)
{
	std::fstream map(mapfile.c_str(),std::ios::in);
	for(int k = 0; k < (1<<dims); ++k)
		for(int d = 0; d < dims; ++d)
			map >> xyz[k][d];
	map.close();
	//complete definition to other dims
	for(int k = (1<<dims); k < 8; k+=(1<<dims))
	{
		for(int q = 0; q < (1<<dims); ++q)
			for(int d = 0; d < dims; ++d)
				xyz[k+q][d] = xyz[k-(1<<dims)+q][d];
	}
	
	
}


void CubeTransform::SetMapping(double * map, int dims)
{
	int q = 0;
	for(int k = 0; k < (1<<dims); ++k)
		for(int d = 0; d < dims; ++d)
			xyz[k][d] = map[q++];
	//complete definition to other dims		
	for(int k = (1<<dims); k < 8; k+=(1<<dims))
	{
		for(int q = 0; q < (1<<dims); ++q)
			for(int d = 0; d < dims; ++d)
				xyz[k+q][d] = xyz[k-(1<<dims)+q][d];
	}
}

void CubeTransform::Transform()
{
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
		for(int d = 0; d < 3; ++d)
			n->Coords()[d] = ((xyz[0][d]*(1-c[0]) + xyz[1][d]*c[0])*(1-c[1])+
			                  (xyz[2][d]*(1-c[0]) + xyz[3][d]*c[0])*c[1])*(1-c[2])+
			                 ((xyz[4][d]*(1-c[0]) + xyz[5][d]*c[0])*(1-c[1])+
			                  (xyz[6][d]*(1-c[0]) + xyz[7][d]*c[0])*c[1])*c[2];
	} 
	
	for(int d = 0; d < 3; ++d)
	{
		cmin[d]=+1.0e+20;
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
		
}	

