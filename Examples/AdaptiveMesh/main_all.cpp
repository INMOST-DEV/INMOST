#include "amesh.h"

using namespace INMOST;

class RefineSets : public AdaptiveMeshCallback
{
	TagReferenceArray element_set;
	Mesh& m;
public:
	RefineSets(Mesh& m) : m(m) {}

	void NewNode(Cell c, Node n, Storage::reference_array cell_hanging_nodes) {}
	void NewNode(Face f, Node n, Storage::reference_array face_hanging_nodes) {}
	void NewNode(Edge e, Node n) {}

	void NewEdge(Cell c, Edge e) {}
	void NewEdge(Face f, Edge e) {}

	void NewFace(Cell c, Face f)
	{

	}

	void CellRefinement(Cell old_cell, ElementArray<Cell>& new_cells, ElementSet new_cell_set)
	{
		Storage::reference_array sets = element_set[old_cell];
		for (Storage::reference_array::iterator it = sets.begin(); it != sets.end(); ++it)
			it->getAsSet().AddElements(new_cells);
	}
	void FaceRefinement(Face old_face, ElementArray<Face>& new_faces)
	{
		Storage::reference_array sets = element_set[old_face];
		for (Storage::reference_array::iterator it = sets.begin(); it != sets.end(); ++it)
			it->getAsSet().AddElements(new_faces);
	}
	void EdgeRefinement(Edge old_edge, ElementArray<Edge>& new_edges) {}

	void CellCoarsening(ElementArray<Cell>& old_cells, Cell new_cell, ElementSet old_cells_set) {}
	void FaceCoarsening(ElementArray<Face>& old_faces, Face new_face) {}
	void EdgeCoarsening(ElementArray<Edge>& old_edges, Edge new_edge) {}

	void BeginRefinement()
	{
		element_set = m.CreateTag("element_sets", DATA_REFERENCE, FACE | CELL, FACE | CELL);
		for (Mesh::iteratorSet set = m.BeginSet(); set != m.EndSet(); ++set)
		{
			for (ElementSet::iterator it = set->Begin(); it != set->End(); ++it)
				element_set[*it].push_back(*set);
		}
	}
	void Refinement() {}
	void EndRefinement()
	{
		element_set = m.DeleteTag(element_set);
	}

	void BeginCoarsening() {}
	void Coarsening() {}
	void EndCoarsening() {}

	//Set indicator to 1 if cell is to be refined and to 0 if not.
	//Indicator is exchanged after all callbacks.
	void RefineIndicator(AdaptiveMesh& am, TagInteger tag_I) {}
	//Set indicator to 0 if cell is not allowed to be coarsened and to 1 otherwise.
	//Indicator is exchanged after all callbacks.
	void CoarseIndicator(AdaptiveMesh& am, TagInteger tag_I) {}
};

int main(int argc, char ** argv)
{
	Mesh::Initialize(&argc,&argv);
	
	if( argc > 1 )
	{
		Mesh mm;
		mm.Load(argv[1]);
		TagInteger indicator = mm.CreateTag("INDICATOR",DATA_INTEGER,CELL,NONE,1);
		
		int max_levels = 1;
		if( argc > 2 ) max_levels = atoi(argv[2]);

		bool make2d = false;
		if (argc > 3) make2d = atoi(argv[3]);
		
		AdaptiveMesh m(mm, false, make2d);
		RefineSets refsets(mm);
		m.AddCallback(&refsets);
		int numref;
		do
		{
			numref = 0;
			for(Mesh::iteratorCell it = mm.BeginCell(); it != mm.EndCell(); ++it)
				if( m.GetLevel(it->self()) < max_levels )
				{
					indicator[it->self()] = 1;
					numref++;
				}
			if( numref )
			{
				if( !m.Refine(indicator) ) break;
				for(Mesh::iteratorCell it = mm.BeginCell(); it != mm.EndCell(); ++it) indicator[it->self()] = 0;
			}
		}
		while(numref);
		std::string file = "out.pmf";
		if( argc > 4 ) file = std::string(argv[4]);
		mm.Save(file);
	}
	else std::cout << "Usage: " << argv[0] << " mesh_file [max_levels=1] [make2d=0] [mesh_out=out.pmf]" << std::endl;
}
