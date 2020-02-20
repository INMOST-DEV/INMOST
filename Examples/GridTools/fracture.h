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
	/// Marker indicates that the cell is the fracture cell.
	/// For face indicates that the face appears between two
	/// fracture cells.
	MarkerType isFracture() const {return fracture_marker;}
	/// Marker indicates that there additional adjacent fracture cells
	/// and multiple fracture cells meet at one of the edges of the face.
	MarkerType MultipleJoints() const {return multiple_fracture_joints;}
	/// Marker indicates that face appears between matrix and fracture.
	MarkerType MatrixFracture() const {return matrix_fracture;}
	///Default constructor, defines empty set.
	Fracture(Mesh & _m) : m(&_m) {}
	///Copy constructor.
	Fracture(const Fracture & b) : m(b.m) {}
	///Assignment operator.
	Fracture & operator = (Fracture const & b) {m = b.m; return * this;}
	/// Retrive cell volume that accounts for fracture aperture.
	INMOST_DATA_REAL_TYPE Volume(Cell c) const;
	/// Retrive face area that accounts for fracture aperture.
	INMOST_DATA_REAL_TYPE Area(Face f) const;
	/// Retrive face center that accounts for fracture aperture.
	void FaceCenter(Face f, INMOST_DATA_REAL_TYPE cnt[3]) const;
	
	/// @param aperture Defines the opening of the fracture. Geometry does not
	/// necesserely follows this parameter, however the functions for Volume and
	/// Area are affected by this parameter.
	/// @param fill_fracture Tells whether it is needed to connect the mesh with cells 
	/// in the fracture region. If parameter is false the mesh is disconnected.
	/// @param gap_multiplier This parameter averages the placement of the new mesh node
	/// position between the original node of the fracture face and averaged centers of cells
	/// that appear on the same side of fracture. If the parameter is 1, then the new node
	/// is positioned at the same location as original node.
	void Open(Tag aperture, bool fill_fracture, double gap_multiplier = 0.95);
};

#endif //__FRACTURE_H
