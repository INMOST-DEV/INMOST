#include "amesh.h"

/// todo:
/// 1. coarsment
/// 2. strategy for faces/edges with faults
/// 3. geom model support
/// 4. make class abstract virtual for user implementation of refinement and coarsment indicators
/// see in code todo:
namespace INMOST
{
	void CleanupSets(ElementSet set)
	{
		ElementSet::iterator it = set.Begin();
		while(it != set.End())
		{
			if( it->isValid() ) ++it;
			else it = set.Erase(it);
		}
		for(ElementSet child = set.GetChild(); child.isValid(); child = child.GetSibling())
			CleanupSets(child);
	}
	
	void ReduceMax(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element->Integer(tag) = std::max(element->Integer(tag),*((const INMOST_DATA_INTEGER_TYPE *)data));
	}
	
	void ReduceMin(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
	{
		(void) size;
		element->Integer(tag) = std::min(element->Integer(tag),*((const INMOST_DATA_INTEGER_TYPE *)data));
	}
	
	
	void AdaptiveMesh::ClearData()
	{
		level         = m->DeleteTag(level);
		hanging_nodes = m->DeleteTag(hanging_nodes);
		parent_set    = m->DeleteTag(parent_set);
		root.DeleteSetTree();
	}
	
	void AdaptiveMesh::PrepareSet()
	{
		//retrive set for coarsening, initialize set if is not present
		if( !root.isValid() )
		{
			root = m->GetSet("ROOT_SET");
			if( root == InvalidElement() )
			{
				root = m->CreateSetUnique("ROOT_SET").first;
				level[root] = 0;
				for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				{
					root.PutElement(it->self());
					parent_set[it->self()] = root.GetHandle();
				}
			}
		}
		if( !m->HaveGlobalID(CELL) ) m->AssignGlobalID(CELL); //for unique set names
	}
	
	AdaptiveMesh::AdaptiveMesh(Mesh & _m) : m(&_m)
	{
		model = NULL;
		//create a tag that stores maximal refinement level of each element
		level = m->CreateTag("REFINEMENT_LEVEL",DATA_INTEGER,CELL|FACE|EDGE|NODE|ESET,NONE,1);
		//create a tag that stores links to all the hanging nodes of the cell
		hanging_nodes = m->CreateTag("HANGING_NODES",DATA_REFERENCE,CELL|FACE,NONE);
		//create a tag that stores links to sets
		parent_set = m->CreateTag("PARENT_SET",DATA_REFERENCE,CELL,NONE,1);
	}
	
	AdaptiveMesh::~AdaptiveMesh()
	{
		//do not delete tags, user may want to repetitively use this class
		//as extension of class mesh in limited code span
	}
	
	bool AdaptiveMesh::Refine(TagInteger & indicator)
	{
		//std::cout << __FUNCTION__ << std::endl;
		static int call_counter = 0;
		int ret = 0; //return number of refined cells
		//initialize tree structure
		PrepareSet();
		int schedule_counter = 1; //indicates order in which refinement will be scheduled
		int scheduled = 1; //indicates that at least one element was scheduled on current sweep
		//0. Extend indicator for edges and faces
		indicator = m->CreateTag(indicator.GetTagName(),DATA_INTEGER,FACE|EDGE,NONE,1);
		while(scheduled)
		{
			//1.Communicate indicator - it may be not synced
			m->ExchangeData(indicator,CELL,0);
			//2.Propogate indicator down to the faces,edges
			//  select schedule for them
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] == schedule_counter )
				{
					ElementArray<Element> adj = c.getAdjElements(FACE|EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if( level[adj[kt]] == level[c] ) //do not schedule finer or coarser elements
							indicator[adj[kt]] = schedule_counter; //refine together with current cell
					}
				}
			}
			//3.Communicate indicator on faces and edges
			m->ExchangeData(indicator,FACE|EDGE,0);
			//4.Check for each cell if there is
			//  any hanging node with adjacent in a need to refine,
			//  schedule for refinement earlier.
			scheduled = 0;
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				//already scheduled cells may be required to be refined first
				//if( indicator[c] == 0 ) //some optimization
				{
					bool scheduled_c = false;
					//any finer level edge is scheduled to be refined first
					ElementArray<Edge> edges = c->getEdges();
					for(ElementArray<Edge>::size_type kt = 0; kt < edges.size() && !scheduled_c; ++kt)
					{
						//if a finer edge is scheduled
						//then this cell should be refined first
						if( indicator[edges[kt]] != 0 &&
							level[edges[kt]] > level[c] &&
							indicator[edges[kt]] >= indicator[c] )
						{
							indicator[c] = schedule_counter+1;
							scheduled++;
							scheduled_c = true;
						}
					}
				}
			}
			//5.Go back to 1 until no new elements scheduled
			scheduled = m->Integrate(scheduled);
			if( scheduled ) schedule_counter++;
		}
		//6.Refine
		m->BeginModification();
		while(schedule_counter)
		{
			Storage::real xyz[3] = {0,0,0};
			//7.split all edges of the current schedule
			for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
			{
				Edge e = m->EdgeByLocalID(it);
				if( !e.Hidden() && indicator[e] == schedule_counter )
				{
					//remember adjacent faces that should get information about new hanging node
					ElementArray<Face> edge_faces = e.getFaces();
					//location on the center of the edge
					for(Storage::integer d = 0; d < m->GetDimensions(); ++d)
						xyz[d] = (e.getBeg().Coords()[d]+e.getEnd().Coords()[d])*0.5;
					//todo: request transformation of node location according to geometrical model
					//create middle node
					Node n = m->CreateNode(xyz);
					//set increased level for new node
					level[n] = level[e.getBeg()] = level[e.getEnd()] = level[e]+1;
					//for each face provide link to a new hanging node
					for(ElementArray<Face>::size_type kt = 0; kt < edge_faces.size(); ++kt)
						hanging_nodes[edge_faces[kt]].push_back(n);
					//split the edge by the middle node
					ElementArray<Edge> new_edges = Edge::SplitEdge(e,ElementArray<Node>(m,1,n.GetHandle()),0);
					//set increased level for new edges
					level[new_edges[0]] = level[new_edges[1]] = level[e]+1;
				}
			}
			//8.split all faces of the current schedule, using hanging nodes on edges
			for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
			{
				Face f = m->FaceByLocalID(it);
				if( !f.Hidden() && indicator[f] == schedule_counter )
				{
					//connect face center to hanging nodes of the face
					Storage::reference_array face_hanging_nodes = hanging_nodes[f];
					//remember adjacent cells that should get information about new hanging node
					//and new hanging edges
					ElementArray<Cell> face_cells = f.getCells();
					//create node at face center
					//f->Centroid(xyz);
					for(int d = 0; d < 3; ++d) xyz[d] = 0.0;
					for(Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
						for(int d = 0; d < 3; ++d) xyz[d] += face_hanging_nodes[kt].getAsNode().Coords()[d];
					for(int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)face_hanging_nodes.size();
					//todo: request transformation of node location according to geometrical model
					//create middle node
					Node n = m->CreateNode(xyz);
					//set increased level for the new node
					level[n] = level[f]+1;
					//for each cell provide link to new hanging node
					for(ElementArray<Face>::size_type kt = 0; kt < face_cells.size(); ++kt)
						hanging_nodes[face_cells[kt]].push_back(n);
					ElementArray<Node> edge_nodes(m,2); //to create new edges
					ElementArray<Edge> hanging_edges(m,face_hanging_nodes.size());
					edge_nodes[0] = n;
					for(Storage::reference_array::size_type kt = 0; kt < face_hanging_nodes.size(); ++kt)
					{
						edge_nodes[1] = face_hanging_nodes[kt].getAsNode();
						hanging_edges[kt] = m->CreateEdge(edge_nodes).first;
						//set increased level for new edges
						level[hanging_edges[kt]] = level[f]+1;
					}
					//split the face by these edges
					ElementArray<Face> new_faces = Face::SplitFace(f,hanging_edges,0);
					//set increased level to new faces
					for(ElementArray<Face>::size_type kt = 0; kt < new_faces.size(); ++kt)
						level[new_faces[kt]] = level[f]+1;
				}
			}
			//this tag helps recreate internal face
			TagReferenceArray internal_face_edges = m->CreateTag("INTERNAL_FACE_EDGES",DATA_REFERENCE,NODE,NODE,4);
			//this marker helps detect edges of current cell only
			MarkerType mark_cell_edges = m->CreateMarker();
			//this marker helps detect nodes hanging on edges of unrefined cell
			MarkerType mark_hanging_nodes = m->CreateMarker();
			//9.split all cells of the current schedule
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( !c.Hidden() && indicator[c] == schedule_counter )
				{
					Storage::reference_array cell_hanging_nodes = hanging_nodes[c]; //nodes to be connected
					//create node at cell center
					for(int d = 0; d < 3; ++d) xyz[d] = 0.0;
					for(Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
						for(int d = 0; d < 3; ++d) xyz[d] += cell_hanging_nodes[kt].getAsNode().Coords()[d];
					for(int d = 0; d < 3; ++d) xyz[d] /= (Storage::real)cell_hanging_nodes.size();
					//c->Centroid(xyz);
					//todo: request transformation of node location according to geometrical model
					//create middle node
					Node n = m->CreateNode(xyz);
					//set increased level for the new node
					level[n] = level[c]+1;
					//retrive all edges of current face to mark them
					ElementArray<Edge> cell_edges = c.getEdges();
					//mark all edges so that we can retive them later
					cell_edges.SetMarker(mark_cell_edges);
					//connect face center to centers of faces by edges
					ElementArray<Node> edge_nodes(m,2);
					ElementArray<Edge> edges_to_faces(m,cell_hanging_nodes.size());
					edge_nodes[0] = n;
					for(Storage::reference_array::size_type kt = 0; kt < cell_hanging_nodes.size(); ++kt)
					{
						assert(cell_hanging_nodes[kt].isValid());
						//todo: unmark hanging node on edge if no more cells depend on it
						edge_nodes[1] = cell_hanging_nodes[kt].getAsNode();
						edges_to_faces[kt] = m->CreateEdge(edge_nodes).first;
						//set increased level for new edges
						level[edges_to_faces[kt]] = level[c]+1;
						//for each node other then the hanging node of the face
						//(this is hanging node on the edge)
						//we record a pair of edges to reconstruct internal faces
						ElementArray<Edge> hanging_edges = cell_hanging_nodes[kt].getEdges(mark_cell_edges,0);
						for(ElementArray<Edge>::size_type lt = 0; lt < hanging_edges.size(); ++lt)
						{
							//get hanging node on the edge
							assert(hanging_edges[lt].getBeg() == cell_hanging_nodes[kt] || hanging_edges[lt].getEnd() == cell_hanging_nodes[kt]);
							Node v = hanging_edges[lt].getBeg() == cell_hanging_nodes[kt]? hanging_edges[lt].getEnd() : hanging_edges[lt].getBeg();
							//mark so that we can collect all of them
							v.SetMarker(mark_hanging_nodes);
							//fill the edges
							Storage::reference_array face_edges = internal_face_edges[v];
							//fill first two in forward order
							//this way we make a closed loop
							assert(face_edges[0] == InvalidElement() || face_edges[2] == InvalidElement());
							if( face_edges[0] == InvalidElement() )
							{
								face_edges[0] = edges_to_faces[kt];
								face_edges[1] = hanging_edges[lt];
							}
							else //last two in reverse
							{
								assert(face_edges[2] ==InvalidElement());
								face_edges[2] = hanging_edges[lt];
								face_edges[3] = edges_to_faces[kt];
							}
						}
					}
					//remove marker from cell edges
					cell_edges.RemMarker(mark_cell_edges);
					//now we have to create internal faces
					ElementArray<Node> edge_hanging_nodes = c.getNodes(mark_hanging_nodes,0);
					ElementArray<Face> internal_faces(m,edge_hanging_nodes.size());
					//unmark hanging nodes on edges
					edge_hanging_nodes.RemMarker(mark_hanging_nodes);
					for(ElementArray<Node>::size_type kt = 0; kt < edge_hanging_nodes.size(); ++kt)
					{
						//create a face based on collected edges
						Storage::reference_array face_edges = internal_face_edges[edge_hanging_nodes[kt]];
						assert(face_edges[0].isValid());
						assert(face_edges[1].isValid());
						assert(face_edges[2].isValid());
						assert(face_edges[3].isValid());
						internal_faces[kt] = m->CreateFace(ElementArray<Edge>(m,face_edges.begin(),face_edges.end())).first;
						//set increased level
						level[internal_faces[kt]] = level[c]+1;
						//clean up structure, so that other cells can use it
						edge_hanging_nodes[kt].DelData(internal_face_edges);
					}
					//split the cell
					ElementArray<Cell> new_cells = Cell::SplitCell(c,internal_faces,0);
					//retrive parent set
					ElementSet parent(m,parent_set[c]);
					//create set corresponding to old coarse cell
					std::stringstream set_name;
					set_name << parent.GetName() << "_C" << c.GlobalID(); //rand may be unsafe
					ElementSet cell_set = m->CreateSetUnique(set_name.str()).first;
					level[cell_set] = level[c]+1;
					//set up increased level for the new cells
					for(ElementArray<Cell>::size_type kt = 0; kt < new_cells.size(); ++kt)
					{
						level[new_cells[kt]] = level[c]+1;
						cell_set.PutElement(new_cells[kt]);
						parent_set[new_cells[kt]] = cell_set.GetHandle();
					}
					parent.AddChild(cell_set);
					//increment number of refined cells
					ret++;
				}
			}
			m->ReleaseMarker(mark_hanging_nodes);
			m->ReleaseMarker(mark_cell_edges);
			m->DeleteTag(internal_face_edges);
			//10.jump to later schedule, and go to 7.
			schedule_counter--;
		}
		//free created tag
		m->DeleteTag(indicator,FACE|EDGE);
		//11. Restore parallel connectivity, global ids
		//ResolveShared(true);
		//ResolveModification();
		//12. Let the user update their data
		//todo: call back function for New() cells
		if( model ) model->Adaptation(*m);
		//13. Delete old elements of the mesh
		m->ApplyModification();
		//14. Done
		m->EndModification();
		//ExchangeData(hanging_nodes,CELL | FACE,0);
		//reorder element's data to free up space
		m->ReorderEmpty(CELL|FACE|EDGE|NODE);
		//return number of refined cells
		call_counter++;
		return ret != 0;
	}
	
	
	
	
	
	
	bool AdaptiveMesh::Coarse(TagInteger & indicator)
	{
		//std::cout << __FUNCTION__ << std::endl;
		static int call_counter = 0;
		//return number of coarsened cells
		int ret = 0;
		//initialize tree structure
		PrepareSet();
		int schedule_counter = 1; //indicates order in which refinement will be scheduled
		int scheduled = 1, unscheduled = 0; //indicates that at least one element was scheduled on current sweep
		//TagInteger coarsened = CreateTag("COARSENED",DATA_INTEGER,CELL,NONE,1);
		TagInteger coarse_indicator = m->CreateTag("COARSE_INDICATOR",DATA_INTEGER,EDGE,NONE,1); //used to find on fine cells indicator on coarse cells
		//0. Extend indicator for sets, edges and faces
		indicator = m->CreateTag(indicator.GetTagName(),DATA_INTEGER,FACE|EDGE,NONE,1);
		while(scheduled || unscheduled)
		{
			// rules
			// a) If there is adjacent finer edge that is not marked for coarsening
			// then this cell should not be coarsened
			// b) If there is adjacent coarser cell, then this cell should be coarsened
			// first
			//0.Communicate indicator - it may be not synced
			m->ExchangeData(indicator,CELL,0);
			//1. Mark each adjacent face/edge for coarsement schedule
			// problem: should mark so that if every adjacent cell is coarsened
			// then adjacent face/edge are also coarsened
			for(ElementType etype = EDGE; etype <= FACE; etype = NextElementType(etype))
			{
				//for(Storage::integer it = 0; it < LastLocalID(etype); ++it) if( isValidElement(etype,it) )
				//	indicator[ElementByLocalID(etype,it)] = 0;
				for(Storage::integer it = 0; it < m->LastLocalID(etype); ++it) if( m->isValidElement(etype,it) )
				{
					Element e = m->ElementByLocalID(etype,it);
					ElementArray<Cell> adj = e.getCells();
					indicator[e] = INT_MAX;
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
						if( level[e] == level[adj[kt]]) indicator[e] = std::min(indicator[e],indicator[adj[kt]]);
					assert(indicator[e] != INT_MAX);
				}
			}
			//2.Communicate indicator on faces and edges
			m->ReduceData(indicator,FACE|EDGE,0,ReduceMin);
			//3.If there is adjacent finer edge that are not marked for coarsening
			// then this cell should not be coarsened
			unscheduled = scheduled = 0;
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementArray<Edge> edges = c.getEdges();
					for(ElementArray<Edge>::size_type kt = 0; kt < edges.size(); ++kt)
					{
						if( level[edges[kt]] > level[c] && indicator[edges[kt]] == 0 )
						{
							indicator[c] = 0;
							unscheduled++;
						}
					}
				}
			}
			//4. Propogate coarsement info over set tree to detect valid coarsenings.
			// go down over sets, if set does not have children and all of the cells
			// of the set are marked for coarsening, then mark the set for coarsement
			// otherwise unmark.
			// Unmark all cells that are not to be coarsened
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementSet parent(m,parent_set[c]);
					//intermediate cell may not be coarsened
					//root set may not have coarsening cells
					if( parent.HaveChild() || !parent.HaveParent() )
					{
						indicator[c] = 0;
						unscheduled++;
					}
					else
					{
						Storage::integer schedule_first = 0;
						bool check = true;
						//check that all elements of the set are to be coarsened
						for(ElementSet::iterator it = parent.Begin(); it != parent.End(); ++it)
						{
							check &= (indicator[it->self()] != 0);
							schedule_first = std::max(schedule_first,indicator[it->self()]);
						}
						if(!check)
						{
							indicator[c] = 0;
							unscheduled++;
						}
						else if( indicator[c] != schedule_first )
						{
							indicator[c] = schedule_first;
							unscheduled++;
						}
					}
				}
			}
			//5.If there is an adjacent coarser element to be refined, then
			//   this one should be scheduled to be refined first
			//a) clean up coarse indicator tag
			for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
				coarse_indicator[m->EdgeByLocalID(it)] = 0;
			//b) each cell mark it's finer edges with cell's schedule
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementArray<Element> adj = c.getAdjElements(EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if( level[adj[kt]] > level[c] ) //only finer edges
							coarse_indicator[adj[kt]] = std::max(coarse_indicator[adj[kt]],indicator[c]);
					}
				}
			}
			//c) data reduction to get maximum over mesh partition
			m->ReduceData(coarse_indicator,EDGE,0,ReduceMax);
			//d) look from cells if any edge is coarsened earlier
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( indicator[c] )
				{
					ElementArray<Element> adj = c.getAdjElements(EDGE);
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
					{
						if( level[c] == level[adj[kt]] && //do not look from coarse cell onto finer edge
						    indicator[c] <= coarse_indicator[adj[kt]])
						{
							indicator[c] = coarse_indicator[adj[kt]]+1;
							scheduled++;
						}
					}
				}
			}
			//5.Go back to 1 until no new elements scheduled
			scheduled = m->Integrate(scheduled);
			unscheduled = m->Integrate(unscheduled);
			if( scheduled ) schedule_counter++;
		}
		//cleanup
		coarse_indicator = m->DeleteTag(coarse_indicator);
		//Make schedule which elements should be refined earlier.
		m->BeginModification();
		while(schedule_counter)
		{
			//unite cells
			//should find and set hanging nodes on faces
			//find single node at the center, all other nodes,
			//adjacent over edge are hanging nodes
			for(Storage::integer it = 0; it < m->CellLastLocalID(); ++it) if( m->isValidCell(it) )
			{
				Cell c = m->CellByLocalID(it);
				if( !c.Hidden() && indicator[c] == schedule_counter )
				{
					//this set contains all the cells to be united
					ElementSet parent(m,parent_set[c]);
					ElementArray<Cell> unite_cells(m,parent.Size());
					//unmark indicator to prevent coarsement with next element
					Storage::integer kt = 0;
					for(ElementSet::iterator jt = parent.Begin(); jt != parent.End(); ++jt)
					{
						unite_cells[kt++] = jt->getAsCell();
						indicator[jt->self()] = 0; //prevent algorithm from visiting again
					}
					//find a node common to all the cells
					ElementArray<Node> center_node = unite_cells[0].getNodes();
					for(kt = 1; kt < unite_cells.size(); ++kt)
						center_node.Intersect(unite_cells[kt].getNodes());
					//only one should be found
					assert(center_node.size() == 1);
					ElementArray<Node> hanging = center_node[0].BridgeAdjacencies2Node(EDGE);
					Cell v = Cell::UniteCells(unite_cells,0);
					//connect hanging nodes to the cell
					assert(hanging_nodes[v].size() == 0);
					for(ElementArray<Node>::size_type kt = 0; kt < hanging.size(); ++kt)
						hanging_nodes[v].push_back(hanging[kt]);
					//set new parent
					parent_set[v] = parent.GetParent().GetHandle();
					//add cell to parent set
					ElementSet(m,parent_set[v]).PutElement(v);
					//set level for new cell
					level[v] = level[c]-1;
					//delete set that contained cells
					//tree structure should be resolved on ApplyModification
					parent.DeleteSet();
					//increment number of coarsened cells
					ret++;
				}
			}
			//unite faces
			//should find and set hanging nodes on edges
			//find single node at the center, all other nodes,
			//adjacent over edge of the face are hanging nodes
			int numcoarsened = 0;
			for(Storage::integer it = 0; it < m->FaceLastLocalID(); ++it) if( m->isValidFace(it) )
			{
				Face f = m->FaceByLocalID(it);
				if( !f.Hidden() && indicator[f] == schedule_counter )
				{
					//one (or both) of the adjacent cells were coarsened and has lower level
					bool visited = false;
					ElementArray<Cell> cells = f.getCells();
					for(ElementArray<Cell>::size_type kt = 0; kt < cells.size(); ++kt)
					{
						assert(level[cells[kt]] < level[f]);
					}
					for(ElementArray<Cell>::size_type kt = 0; kt < cells.size(); ++kt)
					{
						if( level[cells[kt]] < level[f] )
						{
							//cell has one hanging node in common with current face
							ElementArray<Node> nodes = f.getNodes();
							Storage::reference_array search_hanging = hanging_nodes[cells[kt]];
							nodes.Intersect(search_hanging.data(),search_hanging.size());
							assert(nodes.size() == 1);
							//faces that hanging node shares with the cell are
							//those to be united
							ElementArray<Face> unite_faces = cells[kt].getFaces();
							unite_faces.Intersect(nodes[0].getFaces());
							//unmark faces to prevent visit
							for(ElementArray<Face>::size_type lt = 0; lt < unite_faces.size(); ++lt)
								indicator[unite_faces[lt]] = 0;
							//nodes connected by edges to hanging node and
							//common to the cell are hanging nodes on edges
							ElementArray<Node> hanging = cells[kt].getNodes();
							hanging.Intersect(nodes[0].BridgeAdjacencies(EDGE,NODE));
							//unite faces
							Face v = Face::UniteFaces(unite_faces,0);
							//set level for new face
							level[v] = level[f]-1;
							//connect new face to hanging nodes
							for(ElementArray<Node>::size_type lt = 0; lt < hanging.size(); ++lt)
								hanging_nodes[v].push_back(hanging[lt]);
							visited = true;
							numcoarsened++;
							break; //no need to visit the other cell
						}
					}
					assert(visited);
				}
			}
			//unite edges
			for(Storage::integer it = 0; it < m->EdgeLastLocalID(); ++it) if( m->isValidEdge(it) )
			{
				Edge e = m->EdgeByLocalID(it);
				if( !e.Hidden() && indicator[e] == schedule_counter )
				{
					//at least one face must have lower level
					bool visited = false;
					ElementArray<Face> faces = e.getFaces();
					for(ElementArray<Face>::size_type kt = 0; kt < faces.size(); ++kt)
					{
						if( level[faces[kt]] < level[e] )
						{
							//face has one hanging node in common with current edge
							ElementArray<Node> nodes = e.getNodes();
							Storage::reference_array search_hanging = hanging_nodes[faces[kt]];
							nodes.Intersect(search_hanging.data(),search_hanging.size());
							assert(nodes.size() == 1);
							//edges that hanging node shares with the face are those to
							//be united
							ElementArray<Edge> unite_edges = faces[kt].getEdges();
							unite_edges.Intersect(nodes[0].getEdges());
							//unmark edges to prevent visit
							for(ElementArray<Edge>::size_type lt = 0; lt < unite_edges.size(); ++lt)
								indicator[unite_edges[lt]] = 0;
							//unite edges
							Edge v = Edge::UniteEdges(unite_edges,0);
							//set level for new edge
							level[v] = level[e]-1;
							visited = true;
							break; //no need to visit any other face
						}
					}
					assert(visited);
				}
			}
			//jump to later schedule
			schedule_counter--;
		}
		//free created tag
		m->DeleteTag(indicator,FACE|EDGE);
		//todo:
		//ResolveShared(true);
		//ResolveModification();
		//todo:
		//let the user update their data
		if( model ) model->Adaptation(*m);
		m->ApplyModification();
		//done
		m->EndModification();
		//ExchangeData(hanging_nodes,CELL | FACE,0);
		//cleanup null links to hanging nodes
		for(ElementType etype = FACE; etype <= CELL; etype = NextElementType(etype))
		{
			for(Storage::integer it = 0; it < m->LastLocalID(etype); ++it) if( m->isValidElement(etype,it) )
			{
				Storage::reference_array arr = hanging_nodes[m->ElementByLocalID(etype,it)];
				Storage::reference_array::size_type jt = 0;
				for(Storage::reference_array::size_type kt = 0; kt < arr.size(); ++kt)
					if( arr[kt] != InvalidElement() ) arr[jt++] = arr[kt];
				arr.resize(jt);
			}
		}
		//cleanup null links in sets
		CleanupSets(root);
		//reorder element's data to free up space
		m->ReorderEmpty(CELL|FACE|EDGE|NODE|ESET);
		
		call_counter++;
		return ret != 0;
	}
}
