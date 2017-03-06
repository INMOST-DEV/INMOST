#ifndef _AGRID_H
#define _AGRID_H
#include "inmost.h"


class AdaptiveMesh : public Mesh
{
	MarkerType hanging_node;
	
	void SplitFace();
public:
	AdaptiveMesh() : Mesh() {}
	~AdaptiveMesh() {}
	/// Indicator must be 1 on cells to be refined
	/// and 0 on all other cells
	void Refine(const TagInteger & indicator);
	void Coarse(const TagInteger & indicator);
};

#endif //_AGRID_H
