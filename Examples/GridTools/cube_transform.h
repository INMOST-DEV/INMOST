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
	{}
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
