#ifndef __FRACTURE_H
#define __FRACTURE_H

#include "inmost.h"

using namespace INMOST;

class Fracture
{
	Mesh * m;
	Tag fracture_aperture;
	Tag fracture_volume;
	Tag fracture_area;
	MarkerType fracture_marker;
	MarkerType multiple_fracture_joints;
	MarkerType matrix_fracture;
public:
	MarkerType isFracture() const {return fracture_marker;}
	MarkerType MultipleJoints() const {return multiple_fracture_joints;}
	MarkerType MatrixFracture() const {return matrix_fracture;}
	///Default constructor, defines empty set.
	Fracture(Mesh & _m) : m(&_m) {}
	///Copy constructor.
	Fracture(const Fracture & b) : m(b.m) {}
	///Assignment operator.
	Fracture & operator = (Fracture const & b) {m = b.m; return * this;}
	
	
	INMOST_DATA_REAL_TYPE Volume(Cell c) const;
	INMOST_DATA_REAL_TYPE Area(Face f) const;
	void FaceCenter(Face f, INMOST_DATA_REAL_TYPE cnt[3]) const;
	void Open(Tag aperture);
};

#endif //__FRACTURE_H