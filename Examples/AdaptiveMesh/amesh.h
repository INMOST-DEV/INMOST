#ifndef _AMESH_H
#define _AMESH_H
#include "inmost.h"

namespace INMOST
{
	struct AdaptiveMeshCallback
	{
		virtual void NewNode(Cell& c, Node& n, Storage::reference_array cell_hanging_nodes) = 0;
		virtual void NewNode(Face& f, Node& n, Storage::reference_array face_hanging_nodes) = 0;
		virtual void NewNode(Edge& e, Node& n) = 0;

		virtual void NewEdge(Cell& c, Edge& e) = 0;
		virtual void NewEdge(Face& f, Edge& e) = 0;

		virtual void NewFace(Cell& c, Face& f) = 0;

		virtual void CellRefinement(Cell& old_cell, ElementArray<Cell>& new_cells, ElementSet& new_cell_set) = 0;
		virtual void FaceRefinement(Face& old_face, ElementArray<Face>& new_faces) = 0;
		virtual void EdgeRefinement(Edge& old_edge, ElementArray<Edge>& new_edges) = 0;

		virtual void CellCoarsening(ElementArray<Cell>& old_cells, Cell& new_cell, ElementSet& old_cells_set) = 0;
		virtual void FaceCoarsening(ElementArray<Face>& old_faces, Face& new_face) = 0;
		virtual void EdgeCoarsening(ElementArray<Edge>& old_edges, Edge& new_edge) = 0;

		virtual void Adaptation() const = 0;
		virtual void BeginAdaptation() = 0;
		virtual void EndAdaptation() = 0;
	};


	class AdaptiveMesh
	{
		Mesh * m;
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//		Model * model;
//#endif
		ElementSet root; //< Root set that links all the other sets for coarsements
		//TagInteger tag_status;
		TagInteger set_id;
		//TagInteger tag_an; 
        int rank;
        int size;
		/// Prepare sets for coarsements.
		/// Do not do this in constructor, since mesh may contain no cells.
		void PrepareSet();
		void CheckClosure(std::string file, int line);
        //void PrintSetLocal(std::string offset, ElementSet it, std::stringstream& ss);
        //void SynchronizeIndicated(TagInteger& indicator);
		bool skip_tri;
		std::vector<AdaptiveMeshCallback*> callbacks;
	public:
		void ReportSets(std::fstream & fout);
		void CheckParentSet(std::string file, int line);//, TagInteger indicator);
		TagReference parent_set; //<Link to the set that contains an element.
		TagReferenceArray hanging_nodes; //< Link to current hanging nodes of the cell.
		TagReferenceArray tri_hanging_edges; //< Link of the cell to hanging edges of triangles.
		TagInteger level; //< Refinement level of the cell
		//TagReferenceArray ref_tag; //<Link to the set that contains an element.
		Storage::integer GetLevel(const Storage & e) {return level[e];}
        AdaptiveMesh(Mesh & m, bool skip_tri = false);
		~AdaptiveMesh();
		/// Indicator must be 1 on cells to be refined
		/// and 0 on all other cells
		bool Refine(TagInteger & indicator);
		bool Coarse(TagInteger & indicator);
		/// Delete all data related to mesh refinement-coarsement.
		void ClearData();
		void PrintSet(std::ostream & fout, ElementSet set);
//#if defined(USE_AUTODIFF) && defined(USE_SOLVER)
//		void SetModel(Model * mm) {model = mm;}
//#endif
		void AddCallback(AdaptiveMeshCallback* callback) { callbacks.push_back(callback); }
		//the work on each cell is supposed to be proportional to the number of cells refined
		//this number is equal to number of original nodes
		void ComputeWeightRefine(TagInteger indicator, TagReal weight);
		//the work on each cell is supposed to be proportional to the number of cells united
		void ComputeWeightCoarse(TagInteger indicator, TagReal weight);
        //void Test();
        //void PrintMesh(std::ostream& os, int cell = 0, int face = 0, int edge = 0, int node = 0);
        //void PrintSet();
        //void UpdateStatus();
        //void test_sets();
	};
}

#endif //_AMESH_H
