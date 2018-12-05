#ifndef _AMESH_H
#define _AMESH_H
#include "inmost.h"

namespace INMOST
{
	class AdaptiveMesh
	{
		Mesh * m;
		Model * model;
		ElementSet root; //< Root set that links all the other sets for coarsements
		TagInteger level; //< Refinement level of the cell
		TagReference parent_set; //<Link to the set that contains an element.
		TagReferenceArray hanging_nodes; //< Link to current hanging nodes of the cell.
		/// Prepare sets for coarsements.
		/// Do not do this in constructor, since mesh may contain no cells.
		void PrepareSet();
	public:
		Storage::integer GetLevel(const Storage & e) {return level[e];}
		AdaptiveMesh(Mesh & m);
		~AdaptiveMesh();
		/// Indicator must be 1 on cells to be refined
		/// and 0 on all other cells
		bool Refine(TagInteger & indicator);
		bool Coarse(TagInteger & indicator);
		/// Delete all data related to mesh refinement-coarsement.
		void ClearData();
		void SetModel(Model * mm) {model = mm;}
	};
}

#endif //_AMESH_H
