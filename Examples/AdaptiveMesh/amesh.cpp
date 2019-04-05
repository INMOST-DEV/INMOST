#include "amesh.h"
#include <iomanip>
#include <set>

//using namespace std;

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


    void AdaptiveMesh::PrintSetLocal(std::string offset, ElementSet it, std::stringstream& ss)
    {
        std::stringstream ss1;
        ss1 << offset << rank << ": Set : " << std::setw(5) << it.GetName() << " ";
        for (int i = ss1.str().length(); i < 23; i++) ss1 << " ";
        ss << ss1.str();
        ss << std::setw(6);
        if (it.GetStatus() == Element::Shared) ss << "shared";
        else if (it.GetStatus() == Element::Ghost) ss << "ghost";
        else if (it.GetStatus() == Element::Owned) ss << "owned";
        else ss << "none";

        ss << " tag_owner (" << it.IntegerDF(m->OwnerTag()) << ")";

        //ss << "   level (" << level[it.self()] << ")  ";
        ss << " tag_processors (";
        std::stringstream ss2;
        Storage::integer_array arr = it.IntegerArrayDV(m->ProcessorsTag());
        for (int i = 0; i < arr.size(); i++)
            ss2 << arr[i] << " ";
        ss << std::setw(5) << ss2.str() <<")";
        

        ElementSet::iterator p = it.Begin();
        ss << "     | Refs: ";
        int first = 0;
        while(p != it.End())
        {
        //    if (first++ == 0) ss << endl << offset << "   ";
            std::string type = "unknw";
            if (p->GetElementType() == CELL) type = "cell";
            if (p->GetElementType() == FACE) type = "face";
            if (p->GetElementType() == EDGE) type = "edge";
            if (p->GetElementType() == NODE) type = "node";
            ss << type << "-" << std::setw(2) << p->GlobalID() << " ";
            p++;
        }
        ss << std::endl;

        for(ElementSet child = it.GetChild(); child.isValid(); child = child.GetSibling())
        {
            PrintSetLocal(offset + "   ",child,ss);
        }
    }


    void AdaptiveMesh::PrintSet()
    {
        std::stringstream ss;
        for(Mesh::iteratorSet it = m->BeginSet(); it != m->EndSet(); ++it) 
        {
            //if (it->HaveParent()) continue;
            PrintSetLocal("",ElementSet(m,*it),ss);
        }
        std::cout << ss.str() << std::endl;
    }

    void PrintRefs(std::ostream& os, Storage::reference_array refs)
    {
        for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
        {
            std::string type = "unknw";
            if (refs[i].GetElementType() == CELL) type = "cell";
            if (refs[i].GetElementType() == FACE) type = "face";
            if (refs[i].GetElementType() == EDGE) type = "edge";
            if (refs[i].GetElementType() == NODE) type = "node";
            os << "(" << type << "," << refs[i]->GlobalID() << ") ";
			Storage::real xyz[3] = {0,0,0};
            refs[i]->Centroid(xyz);
            os << "(" << xyz[0] << "," << xyz[1] << "," << xyz[2] <<")" <<  std::endl;
        }


    }

    void AdaptiveMesh::PrintMesh(std::ostream& os, int cell, int face, int edge, int node)
    {
        if (cell + face + edge + node == 0) return;
        std::stringstream ss;
        ss << "================= " << rank << " =====================" << std::endl;
        if (cell)
        {
            ss << "Cells: " << m->NumberOfCells() <<  std::endl;
            for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
            {
                ss << rank << ": " << it->GlobalID() << " - " << it->LocalID() << " - ";            
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";
                
                /*
                ss << " (";
                Storage::integer_array arr = it->IntegerArrayDV(tag_processors);
                for (int i = 0; i < arr.size(); i++) ss << arr[i] << " ";
                ss << ") ";
*/
                /*
                Storage::reference_array refs = hanging_nodes[it->self()];
                if (refs.size() > 0)
                {
                    ss << std::endl << "   Hanging nodes: ";
                    PrintRefs(ss,refs);
                }
                ss << std::endl;
                ss << "   ParentSet: " << ElementSet(this,parent_set[*it]).GetName();
                */
/*
                orage::reference_array refs = ref_tag[it->self()];
                PrintRefs(refs);
                */
                    ss << std::endl;
            }
        }

        if (face)
        {
            ss << "Faces: " << m->NumberOfFaces() <<  std::endl;
            for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
            {
                ss << rank << ": " << std::setw(2) << it->LocalID() << " " << std::setw(2) << it->GlobalID() << " - " ;            
                ss << std::setw(6);
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";


                double xyz[3];
                it->Centroid(xyz);
                ss << "   (" << std::setw(5) << xyz[0] << " " << std::setw(5) << xyz[1] << " " << std::setw(5) << xyz[2] << ")";

                ss << "  " << m->GetMarker(*it,m->NewMarker());

                ss << " nc(" << it->getNodes().size() << ": "; 
                ElementArray<Node> nodes = it->getNodes();
                for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
                    ss << std::setw(2) << node->GlobalID() << " ";
                ss << ")";

                /*
                Storage::reference_array refs = ref_tag[it->self()];
                if (refs.size() > 0) ss << ". Ref: ";
                for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
                {
                    string type = "unknw";
                    if (refs[i].GetElementType() == CELL) type = "cell";
                    if (refs[i].GetElementType() == FACE) type = "face";
                    if (refs[i].GetElementType() == EDGE) type = "edge";
                    if (refs[i].GetElementType() == NODE) type = "node";
                    ss << "(" << type << "," << refs[i]->GlobalID() << ") ";
                }

                ss << std::endl;

                {
                    Storage::reference_array refs = hanging_nodes[it->self()];
                    if (refs.size() > 0)
                    {
                        ss << std::endl << "   Hanging nodes: ";
                        PrintRefs(ss,refs);
                    }
                    ss << std::endl;
                }
                */
                ss << std::endl;
            }
        }

        if (edge)
        {
            ss << "Edges: " << m->NumberOfEdges() <<  std::endl;
            for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it) 
            {
                ss << rank << ": " << it->GlobalID() << " - " ;            
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";
				/*
                Storage::reference_array refs = ref_tag[it->self()];
                if (refs.size() > 0) ss << ". Ref: ";
                for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
                {
                    std::string type = "unknw";
                    if (refs[i].GetElementType() == CELL) type = "cell";
                    if (refs[i].GetElementType() == FACE) type = "face";
                    if (refs[i].GetElementType() == EDGE) type = "edge";
                    if (refs[i].GetElementType() == NODE) type = "node";
                    ss << "(" << type << "," << refs[i]->GlobalID() << ") ";
                }
                */
                ss << std::endl;
            }
        }

        if (node)
        {
            ss << "Nodes:" << std::endl;
            for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) 
            {
                ss << rank << ": " << std::setw(2) << it->GlobalID() << " - " ;            
                ss << std::setw(6);
                if (it->GetStatus() == Element::Shared) ss << "shared";
                else if (it->GetStatus() == Element::Ghost) ss << "ghost";
                else ss << "none";

                {
					/*
                    Storage::reference_array refs = ref_tag[it->self()];
                    if (refs.size() > 0) ss << ". Ref: ";
                    for(Storage::reference_array::size_type i = 0; i < refs.size(); ++i)
                    {
                        std::string type = "unknw";
                        if (refs[i].GetElementType() == CELL) type = "cell";
						if (refs[i].GetElementType() == FACE) type = "face";
						if (refs[i].GetElementType() == EDGE) type = "edge";
                        if (refs[i].GetElementType() == NODE) type = "node";
                        ss << "(" << type << "," << refs[i]->GlobalID() << ") ";
                        
                    }
					*/
                    ss << "  " << m->GetMarker(*it,m->NewMarker());

                    ss << "(" << 
                        std::setw(3) << it->RealArray(m->CoordsTag())[0] << " " << 
                        std::setw(3) << it->RealArray(m->CoordsTag())[1] << " " << 
                        std::setw(3) << it->RealArray(m->CoordsTag())[2] << ")";

                }
                ss << std::endl;
            }
        }

        ss << "=========================================" << std::endl;
        os << ss.str() << std::endl;
    }

    void AdaptiveMesh::UpdateStatus()
    {

        for(ElementType mask = CELL; mask >= NODE; mask = PrevElementType(mask))
        {
            for(Mesh::iteratorElement it = m->BeginElement(mask); it != m->EndElement(); it++)
            {
                int stat = 0;
                if (it->GetStatus() == Element::Shared) stat = 1;
                else if (it->GetStatus() == Element::Ghost)  stat = 2;

                tag_status[it->self()] = stat;
            }
        }
    }
    
    void AdaptiveMesh::PrintSet(ElementSet set, std::string offset)
    {
        std::cout << offset << "Set: " << std::endl;

		ElementSet::iterator it = set.Begin();
		while(it != set.End())
        {
		    ElementSet parent(m,parent_set[*it]);
            std::cout << offset << it->GlobalID() << " - " << level[*it] << " : ";
		    
            ElementSet::iterator p = parent.Begin();
    		while(p != parent.End())
            {
                std::cout << p->GlobalID() << " ";
                p++;
            }
            std::cout << std::endl;
            it++;
        }

		for(ElementSet child = set.GetChild(); child.isValid(); child = child.GetSibling())
        {
            PrintSet(child,offset + "  ");
        }

    }

    void AdaptiveMesh::SynchronizeSet(ElementSet set)
    {
    #ifdef USE_MPI
        int size = m->GetProcessorsNumber();
        int rank = m->GetProcessorRank();
        for (int i = 0; i < size; i++)
        {
            set.IntegerArray(m->SendtoTag()).push_back(i);
            m->ExchangeMarked();
        }
    #endif
    }

    void AdaptiveMesh::Test()
    {
        std::cout << rank << ": ================" << std::endl;
        PrintSet(root,"");
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
		tag_status = m->CreateTag("TAG_STATUS",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		//tag_an = m->CreateTag("TAG_AN",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		//ref_tag = m->CreateTag("REF",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
		//create a tag that stores links to all the hanging nodes of the cell
		hanging_nodes = m->CreateTag("HANGING_NODES",DATA_REFERENCE,CELL|FACE,NONE);
		//create a tag that stores links to sets
		parent_set = m->CreateTag("PARENT_SET",DATA_REFERENCE,CELL,NONE,1);
	    size = m->GetProcessorsNumber();
    	rank = m->GetProcessorRank();
	}
	
	AdaptiveMesh::~AdaptiveMesh()
	{
		//do not delete tags, user may want to repetitively use this class
		//as extension of class mesh in limited code span
	}
	
	bool AdaptiveMesh::Refine(TagInteger & indicator)
	{
        std::cout << rank << " enter " << __FUNCTION__ << std::endl;
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
		m->ExchangeData(indicator,CELL | FACE | EDGE,0);
		m->Save("indicator.pmf");
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


        MarkerType marker_new = m->CreateMarker();
        for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it) 
        {
            if (it->GetMarker(m->NewMarker()) == false) continue;
            it->SetMarker(marker_new);

        }

        for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
        {
            if (it->GetMarker(m->NewMarker()) == false) continue;
            it->SetMarker(marker_new);
        }

        for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it) 
        {
            if (it->GetMarker(m->NewMarker()) == false) continue;
            it->SetMarker(marker_new);
        }

        for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) 
        {
            if (it->GetMarker(m->NewMarker()) == false) continue;
            it->SetMarker(marker_new);
        }


		//11. Restore parallel connectivity, global ids
        
        //if (call_counter == 0)
		m->ResolveModification();
        //ExchangeGhost(3,NODE); // Construct Ghost cells in 2 layers connected via nodes
		//12. Let the user update their data
		//todo: call back function for New() cells
		if( model ) model->Adaptation(*m);
		//13. Delete old elements of the mesh
		m->ApplyModification();
		//14. Done
        //cout << rank << ": Before end " << std::endl;
		m->EndModification();
		//ExchangeData(hanging_nodes,CELL | FACE,0);
        //m->ResolveSets();

		m->CheckCentroids();

        //m->BeginModification();
        //    m->ExchangeGhost(1,NODE,marker_new); // Construct Ghost cells in 2 layers connected via nodes
        //    m->ReleaseMarker(marker_new,CELL|FACE|EDGE|NODE);
		//m->ApplyModification();
    	//m->EndModification();
    	//PrintSet();
    	//m->ExchangeData(parent_set,CELL,0);
    	
    			//restore face orientation
		//BUG: bad orientation not fixed automatically
		int nfixed = 0;
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
			if( !it->CheckNormalOrientation() )
			{
				it->FixNormalOrientation();
				nfixed++;
			}
			//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
		if( nfixed ) std::cout << "fixed " << nfixed << " faces" << std::endl;


        //ExchangeData(hanging_nodes,CELL | FACE,0);
        //cout << rank << ": After end " << std::endl;
		//reorder element's data to free up space
		m->ReorderEmpty(CELL|FACE|EDGE|NODE);
		//return number of refined cells
		call_counter++;
        std::cout << rank << " exit " << __FUNCTION__ << std::endl;
		return ret != 0;
	}

    void OperationMin(const Tag & tag, const Element & element, const INMOST_DATA_BULK_TYPE * data, INMOST_DATA_ENUM_TYPE size)
    {
        int value  =  (int)*((int*)data);


        if (value < element->Integer(tag))
        {
            element->Integer(tag) = value;
        }
        (void)size;
    }
	
    void AdaptiveMesh::SynchronizeIndicated(TagInteger& indicator)
    {
        if (m->GetProcessorsNumber() == 1) return;
        int rank = m->GetProcessorRank();

        // Check all sets. All elements in sets must be indicated. At first we check indicator in local processor, and second integrate data
        TagInteger tag_indicated = m->CreateTag("INDICATED",DATA_INTEGER,ESET,NONE,1);
        for(Mesh::iteratorSet it = m->BeginSet(); it != m->EndSet(); ++it) 
        {
            ElementSet set = ElementSet(m,*it);
            set.Integer(tag_indicated) = 1;
            ElementSet::iterator p = set.Begin();
            while(p != set.End())
            {
                if (indicator[*p] == 0)
                {
                    tag_indicated[set] = 0;
                    //std::cout << rank << ": Set " << set.GetName() << " not all indicated" << endl;
                    break;
                }

                p++;
            }
        }

        m->ReduceData(tag_indicated,ESET,0,OperationMin);
        m->ExchangeData(tag_indicated,ESET,0);

        for(Mesh::iteratorSet it = m->BeginSet(); it != m->EndSet(); ++it) 
            if (it->Integer(tag_indicated) == 0)
            {
                ElementSet::iterator p = it->Begin();
                while(p != it->End())
                {
                    p->Integer(indicator) = 0;
                    p++;
                }
            }
        /*
                stringstream ss;
        for(Mesh::iteratorSet it = BeginSet(); it != EndSet(); ++it) 
        {
            ss << it->GetName() << " - " << it->Integer(tag_indicated) << endl;;
        }
        cout << rank << " Sets: \n" << ss.str() << endl;
        */
    }

	bool AdaptiveMesh::Coarse(TagInteger & indicator)
	{
        std::cout << rank << " enter " << __FUNCTION__ << std::endl;
        SynchronizeIndicated(indicator);
        return false;

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
					//if (indicator[e] == INT_MAX) cout << rank << ": " << ElementTypeName(e.GetElementType()) << e.GlobalID() << endl;
					//assert(indicator[e] != INT_MAX);
				}
			}
			//2.Communicate indicator on faces and edges
			m->ReduceData(indicator,FACE|EDGE,0,ReduceMin);
			m->ExchangeData(indicator,FACE|EDGE,0);
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
                //if (!isValidElement(c.GetHandle())) continue;
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
			m->ExchangeData(coarse_indicator,EDGE,0);
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
			m->ExchangeData(indicator,CELL|FACE|EDGE,0);
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
			/*
			for(Storage::integer it = 0; it < NodeLastLocalID(); ++it) if( isValidNode(it) )
			{
				Node e = NodeByLocalID(it);
				if( !e.Hidden() )
				{
					int my_level = -1;
					ElementArray<Edge> edges = e.getEdges();
					for(ElementArray<Edge>::iterator kt = edges.begin(); kt != edges.end(); ++kt)
						my_level = std::max(my_level,level[*kt]);
					level[e] = my_level;
				}
			}
			*/
			//jump to later schedule
			schedule_counter--;
		}
		//free created tag
		m->DeleteTag(indicator,FACE|EDGE);
		//todo:
		m->ResolveModification();
		//todo:
		//let the user update their data
		if( model ) model->Adaptation(*m);
		m->ApplyModification();
		//done
		m->EndModification();
		
		m->CheckCentroids();
		//CheckCentroids();
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
		
		
		//restore face orientation
		//BUG: bad orientation not fixed automatically
		int nfixed = 0;
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) 
			if( !it->CheckNormalOrientation() )
			{
				it->FixNormalOrientation();
				nfixed++;
			}
			//std::cout << "Face " << it->LocalID() << " oriented incorrectly " << std::endl;
		if( nfixed ) std::cout << "fixed " << nfixed << " faces" << std::endl;
		
		//reorder element's data to free up space
		m->ReorderEmpty(CELL|FACE|EDGE|NODE|ESET);
		
		call_counter++;
        std::cout << rank << " exit " << __FUNCTION__ << std::endl;
		return ret != 0;
	}
}
