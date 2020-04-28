#ifndef _CUBE_TRANSFORM_H
#define _CUBE_TRANSFORM_H

#include "inmost.h"
using namespace INMOST;


class CubeTransform
{
	Mesh & mesh;
	double xyz[8][3];
public:
	CubeTransform(Mesh & m) : mesh(m),
	xyz({{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}})
	{
		xyz[0][0] = 0; xyz[0][1] = 0; xyz[0][2] = 0;
		xyz[1][0] = 1; xyz[1][1] = 0; xyz[1][2] = 0;
		xyz[2][0] = 0; xyz[2][1] = 1; xyz[2][2] = 0;
		xyz[3][0] = 1; xyz[3][1] = 1; xyz[3][2] = 0;
		xyz[4][0] = 0; xyz[4][1] = 0; xyz[4][2] = 1;
		xyz[5][0] = 1; xyz[5][1] = 0; xyz[5][2] = 1;
		xyz[6][0] = 0; xyz[6][1] = 1; xyz[6][2] = 1;
		xyz[7][0] = 1; xyz[7][1] = 1; xyz[7][2] = 1;
	}
	CubeTransform(const CubeTransform & b) : mesh(b.mesh) {}
	CubeTransform & operator = (CubeTransform const & b) {mesh = b.mesh; return *this;}
	~CubeTransform() {}
	
	/// Set mapping.
	/// @param map Array map should contain (2^dims)*dims integers.
	/// @param dims Number of dimensions.
	void SetMapping(double * map, int dims = 3);
	/// Load mapping from file.
	/// @param file File name.
	/// @param dims Number of dimensions.
	void LoadMapping(std::string file, int dims = 3);
	/// Output mapping to std::cout.
	void PrintMapping();
	///Perform transformation.
	void Transform();
};


#endif //_CUBE_TRANSFORM_H
