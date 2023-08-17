#include "inmost.h"
#include "incident_matrix.hpp"
#include <queue>
#if defined(USE_MESH)

//TODO:
// incident_matrix class should measure for minimal volume,
// possibly check and update from projects/OctreeCutcell/octgrid.cpp
#if defined(USE_PARALLEL_WRITE_TIME)
__INLINE std::string NameSlash(std::string input)
{
	for(size_t l = input.size(); l > 0; --l)
		if( input[l-1] == '/' || input[l-1] == '\\' )
			return std::string(input.c_str() + l);
	return input;
}
#define REPORT_MPI(x) {WriteTab(out_time) << "<MPI><![CDATA[" << #x << "]]></MPI>" << std::endl; x;}
#define REPORT_STR(x) {WriteTab(out_time) << "<TEXT><![CDATA[" << x << "]]></TEXT>" << std::endl;}
#define REPORT_VAL(str,x) {WriteTab(out_time) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]></CONTENT> <CODE><![CDATA[" << #x << "]]></CODE></VALUE>" << std::endl;}
#define ENTER_FUNC() double all_time = Timer(); {WriteTab(out_time) << "<FUNCTION name=\"" << __FUNCTION__ << "\" id=\"func" << func_id++ << "\">" << std::endl; Enter();}
#define ENTER_BLOCK() { double btime = Timer(); WriteTab(out_time) << "<FUNCTION name=\"" << __FUNCTION__ << ":" << NameSlash(__FILE__) << ":" << __LINE__ << "\" id=\"func" << GetFuncID()++ << "\">" << std::endl; Enter();
#define EXIT_BLOCK() WriteTab(out_time) << "<TIME>" << Timer() - btime << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#define EXIT_FUNC() {WriteTab(out_time) << "<TIME>" << Timer() - all_time << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#define EXIT_FUNC_DIE() {WriteTab(out_time) << "<TIME>" << -1 << "</TIME>" << std::endl; Exit(); WriteTab(out_time) << "</FUNCTION>" << std::endl;}
#else
#define REPORT_MPI(x) x
#define REPORT_STR(x) {}
#define REPORT_VAL(str,x) {}
#define ENTER_FUNC() {}
#define ENTER_BLOCK()
#define EXIT_BLOCK()
#define EXIT_FUNC() {}
#define EXIT_FUNC_DIE()  {}
#endif

const bool optimize_split_geom = true;
const bool optimize_unite_geom = true;

namespace INMOST
{

	


	
	
	void Element::Disconnect(bool del_upper) const
	{
		
		std::vector<adj_type::size_type> del;
		Mesh * m = GetMeshLink();
		//Reorder face edges, so that they always apear in right direction
		if( GetElementType() == CELL ) //This is a cell
			getAsCell()->SwapBackCell();
		if( GetElementType() < CELL || m->HighConnTag().isDefined(CELL) )
		{
			adj_type & hc = m->HighConn(GetHandle()); // node->edges, edge->faces, face->cells
			adj_type::size_type i = 0;
			for( adj_type::size_type it = 0; it < hc.size(); it++)
			{
				int flag = 0;
				adj_type & ilc = m->LowConn(hc[it]); //edge->nodes, face->edges, cell->faces
				adj_type::iterator jt = ilc.begin();
				while (jt != ilc.end())
				{
					if(*jt == GetHandle())
					{
						flag = 1;
						jt = ilc.erase(jt);
					}
					else ++jt;
				}
				if( flag ) del.push_back(i);
				i++;
			}
			if( del_upper )
			{
				if (GetElementType() < CELL)
				{
					if (Hidden())
					{
						for (std::vector<adj_type::size_type>::iterator it = del.begin(); it != del.end(); it++)
							if (m->GetMarker(hc[*it], m->HideMarker())) m->Delete(hc[*it]);
					}
					else
					{
						for (std::vector<adj_type::size_type>::iterator it = del.begin(); it != del.end(); it++)
							m->Delete(hc[*it]);
					}
				}
			}
			hc.clear();
		}
		if( GetElementType() > NODE || m->LowConnTag().isDefined(NODE) )
		{
			adj_type & lc = m->LowConn(GetHandle()); // edge->nodes, face->edges, cell->faces
			for(adj_type::size_type it = 0; it < lc.size(); it++)
			{
				adj_type & ihc = m->HighConn(lc[it]); // node->edges, edge->faces, face->cells
				adj_type::iterator jt = ihc.begin();
				while (jt != ihc.end())
				{
					if(*jt == GetHandle()) jt = ihc.erase(jt);
					else ++jt;
				}
			}
			lc.clear();
		}
	}
	
	Cell Cell::UniteCells(const ElementArray<Cell> & unite, MarkerType del_protect)
	{
		Mesh * m = unite.GetMeshLink();
		if( unite.empty() ) return Cell(m,InvalidHandle());
		std::map<HandleType, int> face_visit; // we check all edges of faces, inner edges are visited several times, outer - once
		ElementArray<Face> faces(m); 
		std::vector<HandleType> inner_faces;
		bool doexit = false;
		MarkerType hm = m->HideMarker();
		//Compute which faces of the provided cells lay on boundary (1 visit) and
		//inside of the set of cells (2 visits)
		for(ElementArray<Cell>::size_type j = 0; j < unite.size(); j++)
		{
			if( unite[j]->GetMarker(del_protect) ) doexit = true;
			//access faces
			adj_type const & lc = m->LowConn(unite.at(j));
			for(adj_type::size_type it = 0; it < lc.size(); it++)
				if( !m->GetMarker(lc[it],hm) ) face_visit[lc[it]]++;
		}
		if( doexit ) return Cell(m,InvalidHandle());
		MarkerType visited = m->CreatePrivateMarker(), rem = m->CreatePrivateMarker();
		std::vector<HandleType> edges;
		std::vector<HandleType> nodes;
		//gather boundary faces into set that will be used to create new cell
		//mark internal faces to be deleted. For internal faces find out
		//all internal edges that should be deleted as well.
		for(std::map<HandleType,int>::iterator it = face_visit.begin(); it != face_visit.end(); it++)
		{
			if( it->second == 1 ) //boundary faces, use for new cell
				faces.push_back(it->first);
			else //internal faces
			{
				if( m->GetMarker(it->first,del_protect) ) doexit = true;
				//mark face to be deleted
				m->SetPrivateMarker(it->first,rem);
				//access edges of the face, gather into array
				adj_type const & lc = m->LowConn(it->first);
				for(adj_type::size_type jt = 0; jt < lc.size(); jt++) if( !m->GetMarker(lc[jt],hm) )
					if( !m->GetPrivateMarker(lc[jt],visited) )
					{
						m->SetPrivateMarker(lc[jt],visited);
						edges.push_back(lc[jt]);
					}
				//gather internal faces
				inner_faces.push_back(it->first);
			}
		}
		//for edges of internal faces gather their nodes,
		//for each edge check if all it's faces are to be deleted,
		//then the edge should be deleted, otherwise keep the edge
		for(size_t i = 0; i < edges.size(); i++) 
		{
			m->RemPrivateMarker(edges[i],visited);
			//access nodes of the edge, gather into array
			adj_type const & lc = m->LowConn(edges[i]);
			for(adj_type::size_type jt = 0; jt < lc.size(); jt++) if( !m->GetMarker(lc[jt],hm) )
			{
				if( !m->GetPrivateMarker(lc[jt],visited) )
				{
					m->SetPrivateMarker(lc[jt],visited);
					nodes.push_back(lc[jt]);
				}
			}
			//access faces of the edge, check is there any
			//face that would not be deleted
			int nonzero = 0;
			adj_type const & hc = m->HighConn(edges[i]); //faces
			for(adj_type::size_type jt = 0; jt < hc.size(); jt++) if( !m->GetMarker(hc[jt],hm) )//iterate over faces
				if( !m->GetPrivateMarker(hc[jt],rem) ) nonzero++; // check if it is not deleted
			if( nonzero == 0 ) //all faces should be deleted, edge to remove
			{
				//mark edge to be deleted
				m->SetPrivateMarker(edges[i],rem);
				if( m->GetMarker(edges[i],del_protect) ) doexit = true;
			}
		}
		//for nodes of internal faces check is there any edges
		//that should not be deleted
		for(size_t i = 0; i < nodes.size(); i++) 
		{
			m->RemPrivateMarker(nodes[i],visited);
			int nonzero = 0;
			//acces edges of the node, check is there any
			//edge that would not be deleted
			adj_type const & hc = m->HighConn(nodes[i]); //edges
			for(adj_type::size_type jt = 0; jt < hc.size(); jt++) if( !m->GetMarker(hc[jt],hm) )
				if( !m->GetPrivateMarker(hc[jt],rem) ) nonzero++;
			if( nonzero == 0 ) //all edges should be deleted, node to remove
			{
				//mark node to be deleted
				m->SetPrivateMarker(nodes[i],rem);
				if( m->GetMarker(nodes[i],del_protect) ) doexit = true;
			}
		}
		m->ReleasePrivateMarker(visited);
		if( doexit )
		{
			m->RemPrivateMarkerArray(inner_faces.data(),(enumerator)inner_faces.size(),rem);
			m->RemPrivateMarkerArray(edges.data(), (enumerator)edges.size(), rem);
			m->RemPrivateMarkerArray(nodes.data(), (enumerator)nodes.size(), rem);
			m->ReleasePrivateMarker(rem);
			return Cell(m,InvalidHandle());
		}
		//delete cells
		for(ElementArray<Cell>::size_type j = 0; j < unite.size(); j++)
		{
			if( m->GetMarker(unite.at(j),del_protect) )
				std::cout << __FUNCTION__ << " deleted protected cells, united " << unite.size() << " cells " << std::endl;
			if( m->HideMarker() )
				m->Hide(unite.at(j));
			else
				m->Delete(unite.at(j));
		}
		//delete inner faces
		for(size_t j = 0; j < inner_faces.size(); j++)
		{
			//std::cout << "delete face " << GetHandleID(inner_faces[j]) << std::endl;
			if( m->GetPrivateMarker(inner_faces[j],rem) )
			{
				m->RemPrivateMarker(inner_faces[j],rem);
				if( m->GetMarker(inner_faces[j],del_protect) )
					std::cout << __FUNCTION__ << " deleted protected faces, united " << unite.size() << " cells " << std::endl;
				if( m->HideMarker() )
					m->Hide(inner_faces[j]);
				else
					m->Delete(inner_faces[j]);
			}
		}
		//delete unused edges
		for(size_t j = 0; j < edges.size(); j++)
		{
			if( m->GetPrivateMarker(edges[j],rem) )
			{
				m->RemPrivateMarker(edges[j],rem);
				if( m->GetMarker(edges[j],del_protect) )
					std::cout << __FUNCTION__ << " deleted protected edge, united " << unite.size() << " cells " << std::endl;
				if( m->HideMarker() )
					m->Hide(edges[j]);
				else
					m->Delete(edges[j]);
			}
		}
		//delete unused nodes
		for(size_t j = 0; j < nodes.size(); j++)
		{

			if( m->GetPrivateMarker(nodes[j],rem) ) //there are no edges that use this edge
			{
				m->RemPrivateMarker(nodes[j],rem);
				// there must be no cells that use this node
				assert(!m->LowConnTag().isDefined(NODE) || m->LowConn(nodes[j]).empty() || m->Count(m->LowConn(nodes[j]).data(), static_cast<integer>(m->LowConn(nodes[j]).size()), hm) == 0);
				if( m->GetMarker(nodes[j],del_protect) )
					std::cout << __FUNCTION__ << " deleted protected node, united " << unite.size() << " cells " << std::endl;
				if( m->HideMarker() )
					m->Hide(nodes[j]);
				else
					m->Delete(nodes[j]);
			}
		}
		m->ReleasePrivateMarker(rem);
		//reconstruct cell by outer faces
		return m->CreateCell(faces).first;
	}
	
	
	bool Cell::TestUniteCells(const ElementArray<Cell> & unite, MarkerType del_protect)
	{
		if( unite.empty() ) return false;
		Mesh * m = unite.GetMeshLink();
		std::map<HandleType,int> face_visit; // we check all edges of faces, inner edges are visited several times, outer - once
		bool doexit = false;
		MarkerType hm = m->HideMarker();			
		for(ElementArray<Cell>::size_type j = 0; j < unite.size(); j++)
		{
			if( m->GetMarker(unite.at(j),del_protect) ) doexit = true; 
			adj_type const & lc = m->LowConn(unite.at(j));
			for(adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) )
				face_visit[lc[it]]++;
		}
		if( doexit ) return false;
		MarkerType rem = m->CreatePrivateMarker();
		std::vector<HandleType> edges;
		std::vector<HandleType> nodes;		
		for(std::map<HandleType,int>::iterator it = face_visit.begin(); it != face_visit.end(); it++)
		{
			if( it->second != 1 )
			{
				if( m->GetMarker(it->first,del_protect) ) doexit = true;
				m->SetPrivateMarker(it->first,rem);
				adj_type const & lc = m->LowConn(it->first);
				for(adj_type::size_type jt = 0; jt < lc.size(); jt++) if( !m->GetMarker(lc[jt],hm) )
					if( !m->GetPrivateMarker(lc[jt],rem) )
					{
						m->SetPrivateMarker(lc[jt],rem);
						edges.push_back(lc[jt]);
					}
			}
		}
		for(size_t i = 0; i < edges.size(); i++) 
		{
			m->RemPrivateMarker(edges[i],rem);
			adj_type const & lc = m->LowConn(edges[i]);
			for(adj_type::size_type jt = 0; jt < lc.size(); jt++) if( !m->GetMarker(lc[jt],hm) )
			{
				if( !m->GetPrivateMarker(lc[jt],rem) )
				{
					m->SetPrivateMarker(lc[jt],rem);
					nodes.push_back(lc[jt]);
				}
			}
			int nonzero = 0;
			//access faces of the edge, check is there any
			//that would not be deleted
			adj_type const & hc = m->HighConn(edges[i]); //faces
			for(adj_type::size_type jt = 0; jt < hc.size(); jt++) if( !m->GetMarker(hc[jt],hm) ) //iterate over faces
				if( !m->GetPrivateMarker(hc[jt],rem) ) nonzero++; // check if it is not deleted
			//all faces should be deleted, edge to remove
			if( nonzero == 0 && m->GetMarker(edges[i],del_protect) )
				doexit = true;
		}
		for(size_t i = 0; i < nodes.size(); i++) 
		{
			int nonzero = 0;
			//access edges of the node, check is there
			//any that would not be deleted
			adj_type const & hc = m->HighConn(nodes[i]); //edges
			for(adj_type::size_type jt = 0; jt < hc.size(); jt++) if( !m->GetMarker(hc[jt],hm) )
				if( !m->GetPrivateMarker(hc[jt],rem) ) nonzero++;
			//all edges should be deleted, node to remove
			if( nonzero == 0 && m->GetMarker(edges[i],del_protect) )
				doexit = true;
		}
		for(std::map<HandleType,int>::iterator it = face_visit.begin(); it != face_visit.end(); it++) m->RemPrivateMarker(it->first,rem);
		m->RemPrivateMarkerArray(edges.data(), (enumerator)edges.size(),rem);
		m->RemPrivateMarkerArray(nodes.data(), (enumerator)nodes.size(), rem);
		m->ReleasePrivateMarker(rem);
		if( doexit ) return false;
		return true;
	}
	
	
	Face Face::UniteFaces(const ElementArray<Face> & unite, MarkerType del_protect)
	{
		Mesh * m = const_cast<Mesh *>(unite.GetMeshLink());
		assert(m != NULL);
		if( unite.empty() ) return Face(m,InvalidHandle());
		MarkerType hm = m->HideMarker();
		bool doexit = false, dothrow = false;
		(void)dothrow;
		for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
			if( m->GetMarker(unite.at(j),del_protect) ) doexit = true;
		if( doexit ) return Face(m,InvalidHandle());
		std::vector<HandleType> cells;
		MarkerType edge_set = m->CreatePrivateMarker();
		MarkerType rem = m->CreatePrivateMarker();

		bool plane = optimize_unite_geom && Face::SamePlane(unite);
		
		//gather cells adjacent to united faces
		for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
		{
			adj_type const & hc = m->HighConn(unite.at(j)); //cells (unite is face)
			for(adj_type::size_type it = 0; it < hc.size(); it++) if( !m->GetMarker(hc[it],hm) )
			{
				if( !m->GetPrivateMarker(hc[it],edge_set) )
				{
					cells.push_back(hc[it]);
					m->SetPrivateMarker(hc[it],edge_set);
				}
			}
		}
		m->RemPrivateMarkerArray(cells.data(), (enumerator)cells.size(), edge_set);
		assert(cells.size() <= 2);
		//check is there a topological problem
		//new face should be adjacent to no more than two cells
		if( m->GetTopologyCheck(TRIPLE_SHARED_FACE) && cells.size() > 2 )
		{
			m->SetTopologyError(TRIPLE_SHARED_FACE);
			if( m->GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(TRIPLE_SHARED_FACE) << std::endl;
			if( m->GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
		}
		
		std::vector<HandleType> nodes;
		std::map<HandleType, int> edge_visit;
		ElementArray<Edge> edges(m);
		//compute how many times each edge is visited
		for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
		{
			//access edges of the face
			adj_type const & lc = m->LowConn(unite.at(j));
			for(adj_type::size_type it = 0; it < lc.size(); it++)
				if( !m->GetMarker(lc[it],hm) ) edge_visit[lc[it]]++;
		}
		HandleType first = InvalidHandle(); //first edge
		HandleType prev = InvalidHandle(); //prev node
		//Mark all the edges on boundary to recreate the face,
		//assemble set of nodes
		int expect_q = 0; //used to test for consistency of the loop of edges
		for(std::map<HandleType,int>::iterator it = edge_visit.begin(); it != edge_visit.end(); it++)
		{
			if( it->second == 1 )
			{
				expect_q++;
				m->SetPrivateMarker(it->first,edge_set);
				if( first == InvalidHandle() ) first = it->first;
			}
			else if( it->second == 2 )
			{
				//edge is protected
				if( m->GetMarker(it->first,del_protect) ) doexit = true;
				//mark edge to be deleted
				m->SetPrivateMarker(it->first,rem);
				//access nodes of the edge
				adj_type const & lc = m->LowConn(it->first);
				for(adj_type::size_type jt = 0; jt < lc.size(); jt++) if( !m->GetMarker(lc[jt],hm) )
				{
					if( !m->GetPrivateMarker(lc[jt],edge_set) )
					{
						m->SetPrivateMarker(lc[jt],edge_set);
						nodes.push_back(lc[jt]);
					}
				}
			}
			else 
			{
				doexit = true; //the faces we unite would not appear as one surface
				dothrow = true;
			}
		}
		//Find out set of nodes to be deleted
		for(size_t j = 0; j < nodes.size(); j++) 
		{
			m->RemPrivateMarker(nodes[j],edge_set);
			int nonzero = 0;
			//access edges of the nodes, find out whether all
			//of them are deleted
			adj_type const & hc = m->HighConn(nodes[j]); // edge
			for(adj_type::size_type it = 0; it < hc.size(); it++) if( !m->GetMarker(hc[it],hm) ) //iterate through edges of the node
				if( !m->GetPrivateMarker(hc[it],rem) ) nonzero++; // check if edge should not be deleted
			if( nonzero == 0 ) //all edges are deleted but the node is protected
			{
				m->SetPrivateMarker(nodes[j],rem);
				if( m->GetMarker(nodes[j],del_protect) ) doexit = true;
			}
		}
			
		if( doexit )
		{
			for(std::map<HandleType,int>::iterator it = edge_visit.begin(); it != edge_visit.end(); it++)
			{
				m->RemPrivateMarker(it->first,rem);
				m->RemPrivateMarker(it->first,edge_set);
			}
			m->RemPrivateMarkerArray(nodes.data(), (enumerator)nodes.size(), rem);
			m->ReleasePrivateMarker(edge_set);
			m->ReleasePrivateMarker(rem);
			assert( !dothrow ); //report the situation, because user need to debug the input
			return Face(m,InvalidHandle());
		}
		//Order edges on the boundary of united faces into loop
		edges.push_back(first);
		edges.back()->RemPrivateMarker(edge_set);
		bool done = false;
		int q = 1;
		while( !done )
		{
			enumerator k1 = ENUMUNDEF,k2;
			//access nodes of the last pushed edge
			adj_type const & lc = m->LowConn(edges.atback()); // nodes
			//find out first node
			k1 = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),k1,hm);
			//find out second node
			k2 = m->getNext(lc.data(),static_cast<enumerator>(lc.size()),k1,hm);
			//find out which of the nodes should connect to the next edge
			//at the beginning we can select any of the two
			if( lc[k1] != prev ) 
				prev = lc[k1]; 
			else 
				prev = lc[k2];
			//find out the next edge connected to the previous node
			bool found = false; //detect that edge was found, otherwise there is no loop
			(void)found;
			adj_type const & hc = m->HighConn(prev); // edge (prev is node)
			for(adj_type::size_type it = 0; it < hc.size(); ++it) if( !m->GetMarker(hc[it],hm) )
			{
				//we came back to the first edge, the loop is closed
				if( hc[it] == first && q != 1)
				{
					found = done = true;
					break;
				}
				if( m->GetPrivateMarker(hc[it],edge_set) )
				{
					q++;
					found = true;
					edges.push_back(hc[it]);
					m->RemPrivateMarker(hc[it],edge_set);
					break;
				}
			}
			assert(found); //there is no loop
			//if( !found ) throw Failure;
		}
		m->ReleasePrivateMarker(edge_set);
		
		//number of edges collected matches number of edges expected
		assert(expect_q == q);
		

		
		
		if( !m->HideMarker() ) //we can't hide elements
		{
			unite.SetPrivateMarker(rem);
			//untie every face from the cell
			for(size_t it = 0; it < cells.size(); it++)
			{
				adj_type & lc = m->LowConn(cells[it]);
				adj_type::iterator jt = lc.begin();
				while( jt != lc.end()) //don't need to check is it hidden
					if( m->GetPrivateMarker(*jt,rem) )
						jt = lc.erase(jt);
					else jt++;
			}
			unite.RemPrivateMarker(rem);
		}
		
		//delete old faces
		if( m->HideMarker() )
			for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
				m->Hide(unite.at(j));
		else
		{
			// delete all faces
			for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
			{
				//untie cells from face
				m->HighConn(unite.at(j)).clear(); // cells (unite is face)
				if( m->GetMarker(unite.at(j),del_protect) )
					std::cout << __FUNCTION__ << " deleted protected face, united " << unite.size() << " faces " << std::endl;
				m->Delete(unite.at(j)); 
			}
		}

		for(std::map<HandleType,int>::iterator it = edge_visit.begin(); it != edge_visit.end(); it++)
			if( it->second != 1 )
			{
				m->RemPrivateMarker(it->first,rem);
				assert( m->HighConn(it->first).empty() || m->Count(m->HighConn(it->first).data(),static_cast<integer>(m->HighConn(it->first).size()),hm) == 0 ); //it's connected to face somewhere
				if( m->GetMarker(it->first,del_protect) )
					std::cout << __FUNCTION__ << " deleted protected edge, united " << unite.size() << " faces " << std::endl;
				if( m->HideMarker() )
					m->Hide(it->first);
				else
					m->Delete(it->first);
			}
			
		

		for(size_t it = 0; it < nodes.size(); it++) //delete nodes inside the face
		{
			if( m->GetPrivateMarker(nodes[it],rem) )
			{
				assert( m->HighConn(nodes[it]).empty() || m->Count(m->HighConn(nodes[it]).data(),static_cast<integer>(m->HighConn(nodes[it]).size()),hm) == 0 ); //it's connected to edge somewhere
				m->RemPrivateMarker(nodes[it],rem);
				if( !m->HideMarker() )
				{
					if (m->LowConnTag().isDefined(NODE))
					{
						adj_type& lc = m->LowConn(nodes[it]);
						adj_type::iterator jt = lc.begin();
						while (jt != lc.end()) // iterate through cells of the node
						{
							adj_type& ihc = m->HighConn(*jt);
							adj_type::iterator qt = ihc.begin(); //iterate through nodes of the cell
							while (qt != ihc.end())
							{
								if (*qt == nodes[it])
									qt = ihc.erase(qt); //remove links from the cell to the node
								else ++qt;
							}
							++jt;
						}
						lc.clear(); // remove links to cells
					}
					if( m->GetMarker(nodes[it],del_protect) )
						std::cout << __FUNCTION__ << " deleted protected node, united " << unite.size() << " faces " << std::endl;
					m->Destroy(nodes[it]);
				}
				else m->Hide(nodes[it]);
			}
		}
		
		m->ReleasePrivateMarker(rem);
				
		Face ret = m->CreateFace(edges).first;
		
		
		assert( ret->GetGeometricType() != MultiLine );
		
		
		adj_type & hc = m->HighConn(ret->GetHandle()); // cells of face
		for(size_t it = 0; it < cells.size(); it++)  //tie new face to old cells
		{
			hc.push_back(cells[it]); // connect new face to cells
			m->LowConn(cells[it]).push_back(ret.GetHandle()); // connect cells to new face
			
		}
		
		//determine orientation for connected cells
		m->RecomputeGeometricData(ret.GetHandle(), ORIENTATION);
		
		for (size_t it = 0; it < cells.size(); it++)  //tie new face to old cells
		{
			assert(m->Count(m->LowConn(cells[it]).data(), m->LowConn(cells[it]).size(), hm) >= 4);
			//compute geometric data
			m->ComputeGeometricType(cells[it]);
			assert(m->GetGeometricType(cells[it]) != MultiPolygon);
		}
		for (size_t it = 0; it < cells.size(); it++) m->RecomputeGeometricData(cells[it], CENTROID);
		if (!plane)
		{
			for (size_t it = 0; it < cells.size(); it++) m->RecomputeGeometricData(cells[it], NORMAL);
			//for (size_t it = 0; it < cells.size(); it++) m->RecomputeGeometricData(cells[it], ORIENTATION);
			for (size_t it = 0; it < cells.size(); it++) m->RecomputeGeometricData(cells[it], MEASURE);
			for (size_t it = 0; it < cells.size(); it++) m->RecomputeGeometricData(cells[it], BARYCENTER);
		}
		for (size_t it = 0; it < cells.size(); it++) m->EndTopologyCheck(cells[it], 0);
		
		return ret;
	}
	
	
	bool Face::TestUniteFaces(const ElementArray<Face> & unite, MarkerType del_protect)
	{
		Mesh * m = const_cast<Mesh *>(unite.GetMeshLink());
		if( unite.size() == 0 ) return false;
		bool doexit = false, dothrow = false;
		(void)dothrow;
		for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
			if( m->GetMarker(unite.at(j),del_protect) ) doexit = true;
		MarkerType hm = m->HideMarker();
		if( doexit ) return false;
		std::vector<HandleType> cells;
		MarkerType rem = m->CreatePrivateMarker();
		for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
		{
			adj_type const & hc = m->HighConn(unite.at(j)); // cells (unite is face)
			for(adj_type::size_type it = 0; it < hc.size(); it++) if( !m->GetMarker(hc[it],hm) )
				if( !m->GetPrivateMarker(hc[it],rem) )
				{
					cells.push_back(hc[it]);
					m->SetPrivateMarker(hc[it],rem);
				}
		}
		m->RemPrivateMarkerArray(cells.data(), (enumerator)cells.size(), rem);
		std::vector<HandleType> nodes;
		std::map<HandleType, int> edge_visit;
		for(ElementArray<Face>::size_type j = 0; j < unite.size(); j++)
		{
			adj_type const & lc = m->LowConn(unite.at(j));
			for(adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) )
				edge_visit[lc[it]]++;
		}
		for(std::map<HandleType,int>::iterator it = edge_visit.begin(); it != edge_visit.end(); it++)
		{
			if( it->second == 2 )
			{
				if( m->GetMarker(it->first,del_protect) ) doexit = true; // edge is protected
				m->SetMarker(it->first,rem);
				adj_type const & lc = m->LowConn(it->first);
				for(adj_type::size_type jt = 0; jt < lc.size(); jt++) if( !m->GetMarker(lc[jt],hm) )
				{
					if( !m->GetPrivateMarker(lc[jt],rem) )
					{
						m->SetPrivateMarker(lc[jt],rem);
						nodes.push_back(lc[jt]);
					}
				}
			}
			else if( it->second != 1 )
			{
				doexit = true;
				dothrow = true;
			}
		}
		for(size_t j = 0; j < nodes.size(); j++) 
		{
			m->RemPrivateMarker(nodes[j],rem);
			int nonzero = 0;
			adj_type const & hc = m->HighConn(nodes[j]); // edges
			for(adj_type::size_type it = 0; it < hc.size(); it++) if( !m->GetMarker(hc[it],hm) )//iterate through edges of the node
				if( !m->GetPrivateMarker(hc[it],rem) ) nonzero++; // check if edge should not be deleted
			//all edges are deleted but the node is protected
			if( nonzero == 0 && m->GetMarker(nodes[j],del_protect)) doexit = true;
		}
		for(std::map<HandleType,int>::iterator it = edge_visit.begin(); it != edge_visit.end(); it++)
			if( it->second != 1 ) m->RemPrivateMarker(it->first,rem);
		m->ReleasePrivateMarker(rem);
		if( doexit )
		{
			assert(!dothrow);
			//if( dothrow ) throw Impossible; //something went the way it shouldn't, possibly bad input
			return false;
		}
		return true;
	}

	Edge Edge::UniteEdges(const ElementArray<Edge> & edges, MarkerType del_protect)
	{
		Mesh * m = const_cast<Mesh *>(edges.GetMeshLink());
		if( edges.size() == 0 ) return Edge(m,InvalidHandle());
		
		bool doexit = false, dothrow = false;
		bool line = optimize_unite_geom && Edge::SameLine(edges);
		(void)dothrow;
		MarkerType hm = m->HideMarker();
		MarkerType rem = m->CreatePrivateMarker();
		std::vector<HandleType> cells;
		std::vector<HandleType> faces;
		std::map<HandleType,int> nodes;
		ElementArray<Node> build_nodes(m);
		for(ElementArray<Edge>::size_type it = 0; it < edges.size(); ++it)
		{
			if( m->GetMarker(edges.at(it),del_protect) )
				doexit = true;
			adj_type const & hc = m->HighConn(edges.at(it)); // faces
			for(adj_type::size_type jt = 0; jt != hc.size(); ++jt) if( !m->GetMarker(hc[jt],hm) )
			{
				if( !m->GetPrivateMarker(hc[jt],rem) )
				{
					faces.push_back(hc[jt]);
					m->SetPrivateMarker(hc[jt],rem);
				}
			}
			adj_type const & lc = m->LowConn(edges.at(it));
			for(adj_type::size_type jt = 0; jt < lc.size(); ++jt) if( !m->GetMarker(lc[jt],hm) )
				nodes[lc[jt]]++;
		}
		
		for(size_t it = 0; it < faces.size(); ++it)
		{
			m->RemPrivateMarker(faces[it],rem);
			adj_type const & hc = m->HighConn(faces[it]); // cells
			for(adj_type::size_type jt = 0; jt < hc.size(); ++jt)
			{
				if( !m->GetPrivateMarker(hc[jt],rem) )
				{
					m->SetPrivateMarker(hc[jt],rem);
					cells.push_back(hc[jt]);
				}
			}
		}
		
		if( !cells.empty() ) m->RemPrivateMarkerArray(&cells[0], (enumerator)cells.size(), rem);
		
		

		if( doexit )
		{
			m->ReleasePrivateMarker(rem);
			return Edge(m,InvalidHandle());
		}
		
		for(std::map<HandleType,int>::iterator it = nodes.begin(); it != nodes.end(); it++)
		{
			
			if( it->second == 1 )
				build_nodes.push_back(it->first);
			else if( it->second == 2 )
			{
				if( m->GetMarker(it->first,del_protect) )
					doexit = true;
			}
			else 
			{
				doexit = true; //edges have some loop
				dothrow = true;
			}
		}
		
		
		assert(build_nodes.size() == 2);

		if( doexit )
		{
			m->ReleasePrivateMarker(rem);
			assert(!dothrow);
			//if( dothrow ) throw Impossible; // bad input
			return Edge(m,InvalidHandle());
		}


		std::vector<adj_type::size_type> insert_pos; //position where we insert new edge

		for(ElementArray<Edge>::size_type it = 0; it < edges.size(); ++it)
			m->SetPrivateMarker(edges.at(it),rem);

		for(size_t it = 0; it < faces.size(); ++it)
		{
			adj_type const & lc = m->LowConn(faces[it]); //edges of face
			bool found_rem = false;
			for(adj_type::size_type k = 0; k < lc.size(); k++) //insert new edge to the first position where we delete old edge
			{
				if( m->GetPrivateMarker(lc[k],rem) )
				{
					//all united edges should appear in consecutive order in deleted face
					/*
					adj_type::size_type sum = 0, j = k;
					while( j < lc.size() && !m->GetMarker(lc[j],hm) && m->GetMarker(lc[j],rem) )
					{
						sum++; j++;
					}
					if( sum != static_cast<adj_type::size_type>(edges.size()) ) //not all edges belong to current face - face will be broken!
					{
						doexit = true;
						dothrow = true;
					}
					 */
					insert_pos.push_back(k);
					found_rem = true;
					break;
				}
			}
			if( !found_rem )
			{
				doexit = true;
				dothrow = true;
			}
		}

		if( doexit )
		{
			m->RemPrivateMarkerArray(edges.data(), (enumerator)edges.size(), rem);
			m->ReleasePrivateMarker(rem);
			assert(!dothrow);
			//if( dothrow ) throw Impossible; // bad input
			return Edge(m,InvalidHandle());
		}

		if( !m->HideMarker() ) //disconnect if cannot hide
		{
			for(size_t it = 0; it < faces.size(); ++it)
			{
				adj_type & lc = m->LowConn(faces[it]);
				adj_iterator jt = lc.begin(); //iterate over edges of faces
				while( jt != lc.end())
				{
					if( m->GetPrivateMarker(*jt,rem) )
						jt = lc.erase(jt);
					else ++jt;
				}
			}
		}


		for(ElementArray<Edge>::size_type it = 0; it < edges.size(); ++it)//delete edges
		{
			m->RemPrivateMarker(edges.at(it),rem);
			if( !m->Hide(edges.at(it)) ) //cannot hide
			{
				m->HighConn(edges.at(it)).clear(); //remove connection from edge to faces
				m->Destroy(edges.at(it));
			}
		}

		m->ReleasePrivateMarker(rem);

		for(std::map<HandleType,int>::iterator it = nodes.begin(); it != nodes.end(); it++)
		{
			if( it->second != 1 )
			{
				if (!m->Hide(it->first)) //cannot hide, we have to untie cells from nodes
				{
					if (m->LowConnTag().isDefined(NODE))
					{
						adj_type& lc = m->LowConn(it->first);
						for (adj_type::size_type qt = 0; qt < lc.size(); ++qt) //iterate through cells of the node
						{
							adj_type& hc = m->HighConn(lc[qt]);
							adj_iterator jt = hc.begin();
							while (jt != hc.end()) // iterate through nodes of the cell
							{
								if (*jt == it->first)
								{
									jt = hc.erase(jt); //erase link to node
								}
								else ++jt;
							}
						}
					}
					m->Destroy(it->first);
				}
			}
		}
			
		Edge e = m->CreateEdge(build_nodes).first;
		adj_type & ehc = m->HighConn(e->GetHandle());

		for (size_t k = 0; k < insert_pos.size(); k++)
		{
			adj_type& lc = m->LowConn(faces[k]);
			std::vector<HandleType> tlc(lc.begin(), lc.end());
			lc.insert(lc.begin() + insert_pos[k], e->GetHandle());
			ehc.push_back(faces[k]);
		}
		bool nonplanar = false;
		MarkerType upd = 0;
		for (size_t k = 0; k < faces.size(); k++) m->ComputeGeometricType(faces[k]);
		for (size_t k = 0; k < faces.size(); k++) m->RecomputeGeometricData(faces[k], CENTROID);
		if (line)
		{
			upd = m->CreatePrivateMarker();
			for (size_t it = 0; it < faces.size(); ++it)
			{
				if (!Face(m, faces[it]).Planarity())
				{
					m->SetPrivateMarker(faces[it], upd);
					nonplanar = true;
				}
			}
			if (!nonplanar)
			{
				m->ReleasePrivateMarker(upd);
				upd = 0;
			}
		}
		if (!line || nonplanar )
		{
			for (size_t k = 0; k < faces.size(); k++) if (!upd || m->GetPrivateMarker(faces[k], upd)) m->RecomputeGeometricData(faces[k], NORMAL);
			for (size_t k = 0; k < faces.size(); k++) if (!upd || m->GetPrivateMarker(faces[k], upd)) m->RecomputeGeometricData(faces[k], ORIENTATION);
			for (size_t k = 0; k < faces.size(); k++) if (!upd || m->GetPrivateMarker(faces[k], upd)) m->RecomputeGeometricData(faces[k], MEASURE);
			for (size_t k = 0; k < faces.size(); k++) if (!upd || m->GetPrivateMarker(faces[k], upd)) m->RecomputeGeometricData(faces[k], BARYCENTER);
		}
		for (size_t k = 0; k < faces.size(); k++) m->EndTopologyCheck(faces[k],0);

		for (size_t it = 0; it < cells.size(); ++it) m->ComputeGeometricType(cells[it]);
		for (size_t it = 0; it < cells.size(); ++it) m->RecomputeGeometricData(cells[it], CENTROID);
		if (upd)
		{
			for (size_t it = 0; it < faces.size(); ++it) if (m->GetPrivateMarker(faces[it], upd))
			{
				Element::adj_type const& hc = m->HighConn(faces[it]);
				for (Element::adj_type::size_type jt = 0; jt < hc.size(); ++jt)
					m->SetPrivateMarker(hc[jt], upd);
			}
		}
		if (!line || nonplanar)
		{
			for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], NORMAL);
			//for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], ORIENTATION);
			for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], MEASURE);
			for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], BARYCENTER);
		}
		for (size_t it = 0; it < cells.size(); ++it) m->EndTopologyCheck(cells[it],0);
		if (upd)
		{
			m->RemPrivateMarkerArray(&faces[0], (enumerator)faces.size(), upd);
			m->RemPrivateMarkerArray(&cells[0], (enumerator)cells.size(), upd);
			m->ReleasePrivateMarker(upd);
		}
		return e;
	}
	bool Edge::TestUniteEdges(const ElementArray<Edge> & edges, MarkerType del_protect)
	{
		if( edges.empty() ) return false;
		Mesh * m = edges.GetMeshLink();
		bool doexit = false, dothrow = false;
		(void)dothrow;
		MarkerType hm = m->HideMarker();
		MarkerType rem = m->CreatePrivateMarker();
		std::vector<HandleType> faces;
		std::map<HandleType,int> nodes;

		for(ElementArray<Edge>::size_type it = 0; it < edges.size(); ++it)
		{
			if( m->GetMarker(edges.at(it),del_protect) )
				doexit = true;
			adj_type const & hc = m->HighConn(edges.at(it)); //faces
			for(adj_type::size_type jt = 0; jt < hc.size(); ++jt) if( !m->GetMarker(hc[jt],hm) )
			{
				if( !m->GetPrivateMarker(hc[jt],rem) )
				{
					faces.push_back(hc[jt]);
					m->SetPrivateMarker(hc[jt],rem);
				}
			}
			adj_type const & lc = m->LowConn(edges.at(it));
			for(adj_type::size_type jt = 0; jt < lc.size(); ++jt) if( !m->GetMarker(lc[jt],hm) )
				nodes[lc[jt]]++;
		}
		
		for(size_t it = 0; it < faces.size(); ++it)
			m->RemPrivateMarker(faces[it],rem);

		if( doexit )
		{
			m->ReleasePrivateMarker(rem);
			return false;
		}
		
		for(std::map<HandleType,int>::iterator it = nodes.begin(); it != nodes.end(); it++)
		{
			if( it->second == 1 )
			{
//				build_nodes.push_back(it->first);
			}
			else if( it->second == 2 )
			{
				if( m->GetMarker(it->first,del_protect) )
					doexit = true;
			}
			else 
			{
				doexit = true;
				dothrow = true;
			}
		}

		if( doexit )
		{
			m->ReleasePrivateMarker(rem);
			assert(!dothrow); //inner loop in deleted edges
			//if( dothrow ) throw Impossible; // bad input
			return false;
		}



		for(ElementArray<Edge>::size_type it = 0; it < edges.size(); ++it) m->SetPrivateMarker(edges.at(it),rem);

		for(size_t it = 0; it < faces.size(); ++it)
		{
			adj_type const & lc = m->LowConn(faces[it]);
			bool found_rem = false;
			for(adj_type::size_type k = 0; k < lc.size(); k++) //insert new edge to the first position where we delete old edge
			{
				if( m->GetPrivateMarker(lc[k],rem) )
				{
					/*
					//all united edges should appear in consecutive order in deleted face
					adj_type::size_type sum = 0, j = k;
					while( !m->GetMarker(lc[j],hm) && m->GetMarker(lc[j],rem) )
					{
						sum++; j++;
					}
					if( sum != static_cast<adj_type::size_type>(edges.size()) )
					{
						doexit = true;
						dothrow = true;
					}
					 */
					found_rem = true;
					break;
				}
			}
			if( !found_rem )
			{
				doexit = true;
				dothrow = true;
			}
		}
		for(ElementArray<Edge>::size_type it = 0; it < edges.size(); ++it) m->RemPrivateMarker(edges.at(it),rem);

		if( doexit )
		{
			m->ReleasePrivateMarker(rem);
			assert(!dothrow);
			//if( dothrow ) throw Impossible; // bad input
			return false;
		}

		return true;
	}
	
	
	Node Edge::CollapseEdge(Edge e)
	{
		Node n = InvalidNode();
		Mesh & m = *e.GetMeshLink();
		//e.PrintElementConnectivity();
		{
			INMOST_DATA_REAL_TYPE a = 0.5, cnt[3] = {0,0,0};
			if( e.getBeg().Boundary() ) a -= 0.5;
			if( e.getEnd().Boundary() ) a += 0.5;
			for(INMOST_DATA_INTEGER_TYPE k = 0; k < m.GetDimensions(); ++k)
				cnt[k] = (1-a)*e.getBeg().Coords()[k] + a*e.getEnd().Coords()[k];
			n = m.CreateNode(cnt);
		}
		TagInteger index = m.CreateTag("element_index",DATA_INTEGER,CELL|FACE|EDGE,NONE,1);
		ElementArray<Node> nodes = e.getNodes();
		ElementArray<Edge> edges(&m), new_edges(&m);
		ElementArray<Face> faces(&m), new_faces(&m);
		ElementArray<Cell> cells(&m), new_cells(&m);
		for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		{
			edges.Unite(kt->getEdges());
			faces.Unite(kt->getFaces());
			cells.Unite(kt->getCells());
		}
		for(ElementArray<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
		{
			ElementArray<Node> edge_nodes = it->getNodes();
			if( edge_nodes[0] == nodes[0] || edge_nodes[0] == nodes[1] ) edge_nodes[0] = n;
			if( edge_nodes[1] == nodes[0] || edge_nodes[1] == nodes[1] ) edge_nodes[1] = n;
			if( edge_nodes[0] != edge_nodes[1] )
			{
				Edge ne = m.CreateEdge(edge_nodes).first;
				new_edges.push_back(ne);
				index[*it] = (INMOST_DATA_INTEGER_TYPE)new_edges.size();
			}
			else index[*it] = -1;
		}
		for(ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
		{
			ElementArray<Edge> face_edges = it->getEdges();
			ElementArray<Edge>::iterator jt = face_edges.begin();
			while(jt != face_edges.end())
			{
				if( index[*jt] == -1 )
					jt = face_edges.erase(jt);
				else
				{
					if( index[*jt] > 0 ) *jt = new_edges[index[*jt]-1].GetHandle();
					jt++;
				}
			}
			if( face_edges.size() >= 3 )
			{
				Face nf = m.CreateFace(face_edges).first;
				new_faces.push_back(nf);
				index[*it] = (INMOST_DATA_INTEGER_TYPE)new_faces.size();
			}
			else index[*it] = -1;
		}
		//transfer edge data
		for(ElementArray<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
		{
			if( index[*it] > 0 )
				Mesh::CopyData(new_edges[index[*it]-1],it->self());
			index[*it] = 0;
		}
		for(ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
		{
			ElementArray<Face> cell_faces = it->getFaces();
			ElementArray<Face>::iterator jt = cell_faces.begin();
			while(jt != cell_faces.end())
			{
				if( index[*jt] == -1 )
					jt = cell_faces.erase(jt);
				else
				{
					if( index[*jt] > 0 ) *jt = new_faces[index[*jt]-1].GetHandle();
					jt++;
				}
			}
			if( cell_faces.size() >= 4 )
			{
				ElementArray<Node> check_nodes(&m);
				for(jt = cell_faces.begin(); jt != cell_faces.end(); ++jt)
					check_nodes.Unite(jt->getNodes());
				bool flattened = false;
				for(jt = cell_faces.begin(); jt != cell_faces.end(); ++jt)
					if( check_nodes.size() == jt->nbAdjElements(NODE) )
						flattened = true;
				if( !flattened )
				{
					Cell nc = m.CreateCell(cell_faces).first;
					new_cells.push_back(nc);
					index[*it] = (INMOST_DATA_INTEGER_TYPE)new_cells.size();
				}
				else index[*it] = -1;
			}
			else index[*it] = -1;
		}
		//transfer face data
		for(ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it)
		{
			if( index[*it] > 0 )
				Mesh::CopyData(new_faces[index[*it]-1],it->self());
			index[*it] = 0;
		}
		//transfer cell data
		for(ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
		{
			if( index[*it] > 0 )
				Mesh::CopyData(new_cells[index[*it]-1],it->self());
			index[*it] = 0;
		}
		//delete old elements
		for(ElementArray<Cell>::iterator it = cells.begin(); it != cells.end(); ++it) it->Delete();
		for(ElementArray<Face>::iterator it = faces.begin(); it != faces.end(); ++it) it->Delete();
		for(ElementArray<Edge>::iterator it = edges.begin(); it != edges.end(); ++it) it->Delete();
		for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++) kt->Delete();
		m.DeleteTag(index);
		return n;
	}


	ElementArray<Edge> Edge::SplitEdge(Edge e, const ElementArray<Node> & nodes, MarkerType del_protect)
	{
		Mesh * m = e->GetMeshLink();
		ElementArray<Edge> ret(m);
		std::vector<HandleType> faces;
		std::vector<HandleType> cells;
		HandleType n[2];
		if( e->GetMarker(del_protect) || nodes.empty() ) return ret;
		ret.reserve(nodes.size()+1);
		MarkerType hm = m->HideMarker();
		MarkerType dup = m->CreatePrivateMarker();
		adj_type & hc = m->HighConn(e->GetHandle()); // faces

		//for geometry recomputation, compute that all nodes are on the line of initial segment
		bool line = optimize_split_geom && e.SameLine(nodes);
		
		for(adj_type::size_type it = 0; it < hc.size(); ++it) if( !m->GetMarker(hc[it],hm) )
		{
			faces.push_back(hc[it]);
			adj_type const & ihc = m->HighConn(hc[it]); //cells
			for(adj_type::size_type jt = 0; jt < ihc.size(); ++jt) if( !m->GetMarker(ihc[jt],hm) )
			{
				if( !m->GetPrivateMarker(ihc[jt],dup) )
				{
					cells.push_back(ihc[jt]);
					m->SetPrivateMarker(ihc[jt],dup);
				}
			}
		}

		for(size_t it = 0; it < cells.size(); ++it)
			m->RemPrivateMarker(cells[it],dup);

		

		m->ReleasePrivateMarker(dup);

		int k = 0;

		adj_type const & lc = m->LowConn(e->GetHandle());

		for(adj_type::size_type it = 0; it < lc.size(); ++it)
			if( !m->GetMarker(lc[it],hm) )
				n[k++] = lc[it];

		assert( k == 2 );

		std::vector<adj_type::size_type> insert_pos; //position where we insert new edges

		for(size_t it = 0; it < faces.size(); ++it)
		{
			bool found = false;
			adj_type const ilc = m->LowConn(faces[it]);
			for(adj_type::size_type jt = 0; jt < ilc.size(); ++jt) if( !m->GetMarker(ilc[jt],hm) )
			{
				if( ilc[jt] == e->GetHandle() )
				{
					insert_pos.push_back(jt);
					found = true;
					break;
				}
			}
			if( !found )
			{
				std::cout << __FILE__ << ":" << __LINE__ << " cannot find edge " << e->GetHandle() << " in face " << faces[it] << " connections " << std::endl;
			}
			assert(found);
		}
		assert(insert_pos.size() == faces.size());
		

		if( !e->Hide() ) //we cannot hide the edge, should disconnect faces
		{
			for(adj_type::size_type it = 0; it < hc.size(); ++it)
			{
				adj_type & ilc = m->LowConn(hc[it]);
				adj_iterator jt = ilc.begin();
				while(jt != ilc.end())
					if( *jt == e->GetHandle() ) jt = ilc.erase(jt);
					else ++jt;
			}
			hc.clear();
			m->Destroy(e->GetHandle());
		}
		
		{
			ElementArray<Node> build_nodes(m,2);
			build_nodes.at(0) = n[0];
			build_nodes.at(1) = nodes.at(0);
			ret.push_back(m->CreateEdge(build_nodes).first->GetHandle());
			assert(!m->GetMarker(ret.back().GetHandle(),hm));
			
			for(ElementArray<Node>::size_type k = 0; k < nodes.size()-1; k++)
			{
				build_nodes.at(0) = nodes.at(k);
				build_nodes.at(1) = nodes.at(k+1);
				ret.push_back(m->CreateEdge(build_nodes).first->GetHandle());
				assert(!m->GetMarker(ret.back().GetHandle(),hm));
			}

			build_nodes.at(0) = nodes.at(nodes.size()-1);
			build_nodes.at(1) = n[1];
			ret.push_back(m->CreateEdge(build_nodes).first->GetHandle());
			assert(!m->GetMarker(ret.back().GetHandle(),hm));
		}

		//connect new edges to faces
		for (size_t it = 0; it < faces.size(); ++it)
		{
			adj_type& lc = m->LowConn(faces[it]);
			//check that that one of the nodes of previous edge match n[0],
			//otherwise we have to insert in reverse
			const adj_type& phc = m->LowConn(lc[(insert_pos[it] + lc.size() - 1) % lc.size()]);
			if (phc[0] == n[0] || phc[1] == n[0])
				lc.insert(lc.begin() + insert_pos[it], ret.begin(), ret.end());
			else
				lc.insert(lc.begin() + insert_pos[it], ret.rbegin(), ret.rend());
		}

		//inform edges that they are connected to faces
		for (ElementArray<Edge>::iterator kt = ret.begin(); kt != ret.end(); ++kt)
		{
			adj_type& hc = m->HighConn(kt->GetHandle()); //faces
			hc.insert(hc.end(), faces.begin(), faces.end());
		}

		//inform nodes that they are connected to cells of this edge
		if (m->LowConnTag().isDefined(NODE))
		{
			for (ElementArray<Node>::const_iterator kt = nodes.begin(); kt != nodes.end(); ++kt)
			{
				adj_type& lc = m->LowConn(kt->GetHandle()); //cells of the node
				lc.insert(lc.end(), cells.begin(), cells.end());
			}
		}

		//inform cells that they are connected to nodes
		if (m->HighConnTag().isDefined(CELL))
		{
			for (size_t it = 0; it < cells.size(); ++it)
			{
				adj_type& hc = m->HighConn(cells[it]); //nodes of the cell
				hc.insert(hc.end(), nodes.begin(), nodes.end());
			}
		}

		
		bool nonplanar = false;
		MarkerType upd = 0;
		for (size_t it = 0; it < faces.size(); ++it) m->ComputeGeometricType(faces[it]);
		for (size_t it = 0; it < faces.size(); ++it) m->RecomputeGeometricData(faces[it], CENTROID);
		if (line)
		{
			upd = m->CreatePrivateMarker();
			for (size_t it = 0; it < faces.size(); ++it)
			{
				if (!Face(m, faces[it]).Planarity())
				{
					m->SetPrivateMarker(faces[it], upd);
					nonplanar = true;
				}
			}
			if (!nonplanar)
			{
				m->ReleasePrivateMarker(upd);
				upd = 0;
			}
		}
		if (!line || nonplanar)
		{
			for (size_t it = 0; it < faces.size(); ++it) if (!upd || m->GetPrivateMarker(faces[it], upd)) m->RecomputeGeometricData(faces[it], NORMAL);
			for (size_t it = 0; it < faces.size(); ++it) if (!upd || m->GetPrivateMarker(faces[it], upd)) m->RecomputeGeometricData(faces[it], ORIENTATION);
			for (size_t it = 0; it < faces.size(); ++it) if (!upd || m->GetPrivateMarker(faces[it], upd)) m->RecomputeGeometricData(faces[it], MEASURE);
			for (size_t it = 0; it < faces.size(); ++it) if (!upd || m->GetPrivateMarker(faces[it], upd)) m->RecomputeGeometricData(faces[it], BARYCENTER);
		}
		for (size_t it = 0; it < faces.size(); ++it) m->EndTopologyCheck(faces[it],0);

		

		for (size_t it = 0; it < cells.size(); ++it) m->ComputeGeometricType(cells[it]);
		for (size_t it = 0; it < cells.size(); ++it) m->RecomputeGeometricData(cells[it], CENTROID);
		if (upd)
		{
			for (size_t it = 0; it < faces.size(); ++it) if (m->GetPrivateMarker(faces[it], upd))
			{
				Element::adj_type const& hc = m->HighConn(faces[it]);
				for (Element::adj_type::size_type jt = 0; jt < hc.size(); ++jt)
					m->SetPrivateMarker(hc[jt], upd);
			}
		}
		if (!line || nonplanar)
		{
			for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], NORMAL);
			//for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], ORIENTATION);
			for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], MEASURE);
			for (size_t it = 0; it < cells.size(); ++it) if (!upd || m->GetPrivateMarker(cells[it], upd)) m->RecomputeGeometricData(cells[it], BARYCENTER);
		}
		for (size_t it = 0; it < cells.size(); ++it) m->EndTopologyCheck(cells[it],0);
		if (upd)
		{
			m->RemPrivateMarkerArray(&faces[0], (enumerator)faces.size(), upd);
			m->RemPrivateMarkerArray(&cells[0], (enumerator)cells.size(), upd);
			m->ReleasePrivateMarker(upd);
		}
		return ret;
	}
	
	bool Edge::TestSplitEdge(Edge e, const ElementArray<Node> & nodes, MarkerType del_protect)
	{
		return !nodes.empty() && !e->GetMarker(del_protect);
	}

	
	ElementArray<Face> Face::SplitFace(Face face, const ElementArray<Edge> & edges, MarkerType del_protect)
	{
		//bool did_restart = false;
		//bool report = false;
		//restart_algorithm:
		Mesh * m = face->GetMeshLink();
		ElementArray<Edge> loop(m);
		ElementArray<Face> ret(m);
		ElementArray<Node> cnodes(m);
		std::vector<HandleType> temp;
		if( edges.empty() || face->GetMarker(del_protect) ) return ret;
		MarkerType hm = m->HideMarker();
		std::vector<HandleType> cells;
		
		//if( report ) std::cout << "Marker for hidden elements: " << hm << std::endl;

		bool plane = optimize_split_geom && face.Planarity() && face.SamePlane(edges);

		//collect nodes that should be connected to cell
		//these are nodes not present at face
		if (m->HighConnTag().isDefined(CELL) || m->LowConnTag().isDefined(NODE))
		{
			MarkerType mrk = m->CreatePrivateMarker();
			ElementArray<Node> fnodes = face.getNodes();
			fnodes.SetPrivateMarker(mrk);
			for (ElementArray<Edge>::const_iterator it = edges.begin(); it != edges.end(); ++it)
			{
				Element::adj_type const& lc = m->LowConn(*it); //edge nodes
				for (Element::adj_type::const_iterator ilc = lc.begin(); ilc != lc.end(); ++ilc) if (!m->GetMarker(*ilc, hm))
				{
					if (!m->GetPrivateMarker(*ilc, mrk))
					{
						cnodes.push_back(*ilc);
						m->SetPrivateMarker(*ilc, mrk);
					}
				}
				//ElementArray<Node> enodes = it->getNodes(mrk, true); //nodes does not have marker
				//cnodes.insert(cnodes.end(), enodes.begin(), enodes.end());
				//enodes.SetPrivateMarker(mrk);
			}
			fnodes.RemPrivateMarker(mrk);
			cnodes.RemPrivateMarker(mrk);
			m->ReleasePrivateMarker(mrk);
		}

		adj_type & hc = m->HighConn(face->GetHandle()); // cells
		//if( report ) std::cout << "Adjacent cells for face: ";
		for(adj_type::size_type it = 0; it < hc.size(); ++it) 
		{
			if( !m->GetMarker(hc[it],hm) )
			{
				cells.push_back(hc[it]);
				//if( report ) std::cout << GetHandleID(hc[it]) << " ";
			}
		}
		//if( report ) std::cout << std::endl;

		//assert(cells.size() == 2);

		
		MarkerType outer = m->CreatePrivateMarker();

		int ninner = 0;
		adj_type & lc = m->LowConn(face->GetHandle());
		//if( report ) std::cout << "Edges for face: ";
		for(adj_type::size_type it = 0; it < lc.size(); ++it) if( !m->GetMarker(lc[it],hm) )
		{
			m->SetPrivateMarker(lc[it],outer);
			//if( report ) std::cout << GetHandleID(lc[it]) << " ";
		}
		//if( report ) std::cout << std::endl;
		
		
		//if( report ) std::cout << "Split edges: ";
		for(ElementArray<Edge>::size_type it = 0; it < edges.size(); ++it)
			if( !m->GetPrivateMarker(edges[it].GetHandle(),outer) )
			{
				temp.push_back(edges[it].GetHandle());
				//if( report ) std::cout << GetHandleID(edges[it].GetHandle()) << " ";
				ninner++;
			}
		//if( report ) std::cout << std::endl;

		for(adj_type::size_type it = 0; it < lc.size(); ++it) if( !m->GetMarker(lc[it],hm) )
		{
			temp.push_back(lc[it]);
			m->RemPrivateMarker(lc[it],outer);
		}
		m->ReleasePrivateMarker(outer);

		//all provided edges coincide with boundary edges of current face, no action required
		if( ninner == 0 ) 
		{
			ret.push_back(face);
			return ret;
		}
		
		/*
		if( !did_restart )
		{
			incident_matrix<Edge> mat(m,temp.begin(),temp.end(),ninner);
			do
			{
				mat.find_shortest_loop(loop);
				if (!loop.empty())
				{
					if( loop.size() > 2 )
					{
						//ret.push_back(m->CreateFace(loop).first);
						//adj_type & hc = m->HighConn(ret.back()->GetHandle());
						//hc.replace(hc.begin(),hc.end(),cells.begin(),cells.end());
					}
					else report = true;
				}
				else break;
			} while(true);
			
			
			if(report )//!mat.all_visited())
			{
				mat.print_matrix();
				incident_matrix<Edge> mat(m,temp.begin(),temp.end(),ninner,GridCoords(),true);
				do
				{
					mat.find_shortest_loop(loop);
					if( loop.empty() ) break;
				} while( true );
				did_restart = true;
				goto restart_algorithm;
			}
		}
		*/

		
		if( !face->Hide() )
		{
			for(size_t k = 0; k < cells.size(); k++)
			{
				adj_type & ilc = m->LowConn(cells[k]);
				for(adj_type::size_type it = 0; it < ilc.size(); ++it)
				{
					if( ilc[it] == face->GetHandle() )
					{
						ilc.erase(ilc.begin()+it);
						break;
					}
				}
			}
			hc.clear();
			m->Destroy(face->GetHandle());
		}
		
		incident_matrix<Edge> mat(m,temp.begin(),temp.end(),ninner);
		
		//mat.print_matrix();
		//std::cin.get();

		
		do
		{
			mat.find_shortest_loop(loop);
			if (!loop.empty())
			{
				if( loop.size() > 2 )
				{
					std::pair<Face,bool> ff = m->CreateFace(loop);
					ret.push_back(ff.first); //collect new faces
					adj_type & hc = m->HighConn(ret.back()->GetHandle()); //cells
					if( ff.second ) //connect brand new face to cells
					{
						//std::cout << "hello1 from " << ff.first.LocalID() << std::endl;
						hc.replace(hc.begin(),hc.end(),cells.begin(),cells.end());
						//m->RecomputeGeometricData(ff.first.GetHandle(), ORIENTATION);
					}
					else //face already existed before and may be connected to some cells
					{
						//std::cout << "hello2 from " << ff.first.LocalID() << std::endl;
						//std::cout << "cells: ";
						//for(int k = 0; k < (int)cells.size(); ++k) std::cout << GetHandleID(cells[k]) << " ";
						//std::cout << std::endl;
						//std::cout << "conns: ";
						//for(int k = 0; k < (int)hc.size(); ++k) std::cout << GetHandleID(hc[k]) << " ";
						//std::cout << std::endl;
						for(int k = 0; k < (int)cells.size(); ++k)
						{
							bool add = true;
							for(int j = 0; j < (int)hc.size() && add; ++j)
								if( hc[j] == cells[k] ) add = false;
							if (add) hc.push_back(cells[k]);
						} //FIXME: what if more then 2 connections? have to fire topology exception, unless there are hidden cells
					}
				}
			}
			else break;
		} while(true);



		for (size_t it = 0; it < cells.size(); ++it)
		{
			//connect cell faces 
			adj_type& lc = m->LowConn(cells[it]); //cell faces
			lc.insert(lc.end(), ret.begin(), ret.end());
			//connect cell nodes
			if (m->HighConnTag().isDefined(CELL))
			{
				adj_type& hc = m->HighConn(cells[it]);
				hc.insert(hc.end(), cnodes.begin(), cnodes.end());
				//do we need to check duplicates?
			}
			//connect nodes to cell
			if (m->LowConnTag().isDefined(NODE))
			{
				for (ElementArray<Node>::iterator jt = cnodes.begin(); jt != cnodes.end(); ++jt)
				{
					adj_type& nlc = m->LowConn(*jt); //node cells
					//nlc.push_back(cells[it]);
					//do we need to check duplicates?
					bool have_cell = false;
					for (adj_type::iterator kt = nlc.begin(); kt != nlc.end() && !have_cell; ++kt)
						if (*kt == cells[it]) have_cell = true;
					if (!have_cell)
						nlc.push_back(cells[it]);
				}
			}
		}

		for(unsigned it = 0; it < ret.size(); ++it)
			m->RecomputeGeometricData(ret[it].GetHandle(), ORIENTATION);

		for (size_t it = 0; it < cells.size(); ++it) m->ComputeGeometricType(cells[it]);
		for (size_t it = 0; it < cells.size(); ++it) m->RecomputeGeometricData(cells[it], CENTROID);
		if (!plane)
		{
			for (size_t it = 0; it < cells.size(); ++it) m->RecomputeGeometricData(cells[it], NORMAL);
			//for (size_t it = 0; it < cells.size(); ++it) m->RecomputeGeometricData(cells[it], ORIENTATION);
			for (size_t it = 0; it < cells.size(); ++it) m->RecomputeGeometricData(cells[it], MEASURE);
			for (size_t it = 0; it < cells.size(); ++it) m->RecomputeGeometricData(cells[it], BARYCENTER);
		}
		for (size_t it = 0; it < cells.size(); ++it) m->EndTopologyCheck(cells[it],0);

		return ret;
	}
	bool Face::TestSplitFace(Face face, const ElementArray<Edge> & edges, MarkerType del_protect)
	{
		return !edges.empty() && !face->GetMarker(del_protect);
	}


	ElementArray<Cell> Cell::SplitCell(Cell cell, const ElementArray<Face> & faces, MarkerType del_protect)
	{
		Mesh * m = cell->GetMeshLink();
		ElementArray<Cell> ret(m);
		ElementArray<Face> loop(m);
		std::vector<HandleType> temp;
		if( faces.empty() || cell->GetMarker(del_protect) ) return ret;
		MarkerType hm = m->HideMarker();



		MarkerType outer = m->CreatePrivateMarker();
		int ninner = 0;
		adj_type & lc = m->LowConn(cell->GetHandle());
		for(adj_type::size_type it = 0; it < lc.size(); ++it) if( !m->GetMarker(lc[it],hm) )
			m->SetPrivateMarker(lc[it],outer);

		for(ElementArray<Face>::size_type it = 0; it < faces.size(); ++it)
			if( !m->GetPrivateMarker(faces[it].GetHandle(),outer) )
			{
				temp.push_back(faces[it].GetHandle());
				ninner++;
			}

		for(adj_type::size_type it = 0; it < lc.size(); ++it) if( !m->GetMarker(lc[it],hm) )
		{
			temp.push_back(lc[it]);
			m->RemPrivateMarker(lc[it],outer);
		}
		m->ReleasePrivateMarker(outer);

		//all provided faces coincide with boundary faces of current cell, no action required
		if( ninner == 0 )
		{
			ret.push_back(cell);
			return ret;
		}
			
		incident_matrix<Face> mat(m,temp.begin(),temp.end(),ninner);
		//mat.print_matrix();
		cell->Delete();
		
		bool report = false;
		do
		{
			mat.find_shortest_loop(loop);
			if (!loop.empty() ) 
			{
				if( loop.size() > 3 )
					ret.push_back(m->CreateCell(loop).first);
				else
					report = true;
			}
			else break;
		} while( true );
		
		//debug
		
		if( false ) if(report || !mat.all_visited())
		{
			mat.print_matrix();
			incident_matrix<Face> mat(m,temp.begin(),temp.end(),ninner,GridCoords(),true);
			do
			{
				mat.find_shortest_loop(loop);
				if( loop.empty() ) break;
			} while( true );
		}
		 

		return ret;
	}
	bool Cell::TestSplitCell(Cell cell, const ElementArray<Face> & faces, MarkerType del_protect)
	{
		return !faces.empty() && !cell->GetMarker(del_protect);
	}
	
	void Mesh::BeginModification()
	{
		ENTER_FUNC();
		hide_element = CreateMarker();
		new_element = CreateMarker();
		//update_geometry = 0;
		//update_geometry = CreateMarker();
		EXIT_FUNC();
	}
	
	void Mesh::SwapModification(bool recompute_geometry)
	{
		//~ ENTER_FUNC();
		MarkerType temp = hide_element;
		hide_element = new_element;
		new_element = temp;

		integer tmp[6];
		memcpy(tmp,hidden_count,sizeof(integer)*6);
		memcpy(hidden_count,hidden_count_zero,sizeof(integer)*6);
		memcpy(hidden_count_zero,tmp,sizeof(integer)*6);

		if( recompute_geometry ) //TODO ????????
		{
			for(ElementType etype = EDGE; etype <= CELL; etype = etype << 1)
			{
				for(integer it = 0; it < LastLocalID(etype); ++it) if( isValidElement(etype,it) )
				{
					HandleType h = ComposeHandle(etype,it);
					if( GetMarker(h,NewMarker()) )//&& !GetMarker(h,UpdateGeometryMarker()) ) //element is new and geometry was already recomputed
					{
						ComputeGeometricType(h);
						RecomputeGeometricData(h);
					}
				}
			}
		}
		//~ EXIT_FUNC();
	}

	void Mesh::RemoveLinksToDeletedElements(MarkerType mrk)
	{
		ENTER_FUNC();
		ENTER_BLOCK();
		for (Mesh::iteratorTag it = BeginTag(); it != EndTag(); ++it)
		{
			Tag t = *it;
			if (t.GetDataType() == DATA_REFERENCE)
			{
				if (t == HighConnTag() || t == LowConnTag()) continue;
				for (ElementType etype = NODE; etype <= MESH; etype = etype << 1)
					if (t.isDefined(etype))
					{
						if (t.isSparse(etype))
						{
							if (isPrivate(mrk))
							{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
								for (integer jt = 0; jt < LastLocalID(etype); ++jt) if (isValidElement(etype, jt))
								{
									HandleType h = ComposeHandle(etype, jt);
									if (HaveData(h, t))
									{
										reference_array arr = ReferenceArray(h, t);
										for (reference_array::size_type qt = 0; qt < arr.size(); ++qt)
											if (arr.at(qt) != InvalidHandle() && GetPrivateMarker(arr.at(qt), mrk)) arr.at(qt) = InvalidHandle();
									}
								}
							}
							else
							{
#if defined(USE_OMP)
#pragma omp parallel for
#endif

								for (integer jt = 0; jt < LastLocalID(etype); ++jt) if (isValidElement(etype, jt))
								{
									HandleType h = ComposeHandle(etype, jt);
									if (HaveData(h, t))
									{
										reference_array arr = ReferenceArray(h, t);
										for (reference_array::size_type qt = 0; qt < arr.size(); ++qt)
											if (arr.at(qt) != InvalidHandle() && GetMarker(arr.at(qt),mrk)) arr.at(qt) = InvalidHandle();
									}
								}
							}
						}
						else
						{
							if (isPrivate(mrk))
							{
#if defined(USE_OMP)
#pragma omp parallel for
#endif

								for (integer jt = 0; jt < LastLocalID(etype); ++jt) if (isValidElement(etype, jt))
								{
									HandleType h = ComposeHandle(etype, jt);
									reference_array arr = ReferenceArray(h, t);
									for (reference_array::size_type qt = 0; qt < arr.size(); ++qt)
										if (arr.at(qt) != InvalidHandle() && GetPrivateMarker(arr.at(qt), mrk)) arr.at(qt) = InvalidHandle();
								}
							}
							else
							{
#if defined(USE_OMP)
#pragma omp parallel for
#endif

								for (integer jt = 0; jt < LastLocalID(etype); ++jt) if (isValidElement(etype, jt))
								{
									HandleType h = ComposeHandle(etype, jt);
									reference_array arr = ReferenceArray(h, t);
									for (reference_array::size_type qt = 0; qt < arr.size(); ++qt)
										if (arr.at(qt) != InvalidHandle() && GetMarker(arr.at(qt),mrk)) arr.at(qt) = InvalidHandle();
								}
							}
						}
					}
			}
		}
		EXIT_BLOCK();
		ENTER_BLOCK();
//somehow openmp does not work here, maybe non-private markers
//#if defined(USE_OMP)
//#pragma omp parallel for schedule(dynamic,1)
//#endif
		for (integer jt = 0; jt < EsetLastLocalID(); ++jt) if (isValidElementSet(jt))
		//for (Mesh::iteratorSet it = BeginSet(); it != EndSet(); it++)
		{
			ElementSet it = EsetByLocalID(jt);
			//std::cout << "set name: " << it->GetName() << " size " << it->Size() << " id " << it->LocalID() << std::endl;
			while (it->HaveChild() && it->GetChild()->GetMarker(mrk))
				it->RemChild(it->GetChild());
			while (it->HaveSibling() && it->GetSibling()->GetMarker(mrk))
				it->RemSibling(it->GetSibling());
			if (it->HaveParent() && it->GetParent()->GetMarker(mrk))
				it->GetParent()->RemChild(it->self());

			ElementSet::iterator jt = it->Begin();
			//int q = 0;
			while (jt != it->End())
			{
				//++q;
				//std::cout << "check element " << ElementTypeName(jt->GetElementType()) << " num " << jt->LocalID() << " handle " << jt->GetHandle() << std::endl;
				if (jt->GetMarker(mrk) && jt->GetElementType() != MESH)
				{
					//std::cout << "erase element " << ElementTypeName(jt->GetElementType()) << " num " << jt->LocalID() << " handle " << jt->GetHandle() << std::endl;
					jt = it->Erase(jt);
				}
				else ++jt;
			}
			//std::cout << "size " << it->Size() << " traversed " << q << std::endl;
			//it->Subtract(erase); //old approach
		}
		EXIT_BLOCK();
		EXIT_FUNC();
	}

	//from parallel.cpp
	void DeleteUnpack(const Tag& tag, const Element& e, const INMOST_DATA_BULK_TYPE* data, INMOST_DATA_ENUM_TYPE size);
	
	void Mesh::ApplyModification()
	{
		ENTER_FUNC();
		ENTER_BLOCK();
		temp_hide_element = hide_element;
		hide_element = 0;
		RemoveLinksToDeletedElements(temp_hide_element);
		hide_element = temp_hide_element;
		/*
		for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); ++it)
		{
			Tag t = *it;
			if( t.GetDataType() == DATA_REFERENCE )
			{
				if( t == HighConnTag() || t  == LowConnTag() ) continue;
				for(ElementType etype = NODE; etype <= CELL; etype = etype << 1)
					if( t.isDefined(etype) )
					{
						if( t.isSparse(etype) )
						{
							for(integer jt = 0; jt < LastLocalID(etype); ++jt) if( isValidElement(etype,jt) )
							{
								HandleType h = ComposeHandle(etype,jt);
								if( HaveData(h,t) )
								{
									reference_array arr = ReferenceArray(h,t);
									for(reference_array::size_type qt = 0; qt < arr.size(); ++qt)
										if( arr.at(qt) != InvalidHandle() && Hidden(arr.at(qt)) ) arr.at(qt) = InvalidHandle();
								}
							}
						}
						else
						{
							for(integer jt = 0; jt < LastLocalID(etype); ++jt) if( isValidElement(etype,jt) )
							{
								HandleType h = ComposeHandle(etype,jt);
								reference_array arr = ReferenceArray(h,t);
								for(reference_array::size_type qt = 0; qt < arr.size(); ++qt)
									if( arr.at(qt) != InvalidHandle() && Hidden(arr.at(qt)) ) arr.at(qt) = InvalidHandle();
							}
						}
					}
			}
		}
		*/
		EXIT_BLOCK();
		//need to gather the set of deleted elements
		/*//old approach
		ElementSet erase = CreateSet("TEMPORARY_ERASE_SET").first;
		for(ElementType etype = NODE; etype <= ESET; etype = etype << 1)
		{
			for(integer it = 0; it < LastLocalID(etype); ++it)
				if( isValidElement(etype,it) && Hidden(ComposeHandle(etype,it)) )
					erase->PutElement(ComposeHandle(etype,it)); 
		}
		//the formed above set will be automatically sorted by handles, it does not harm to explicitly indicate that,
		//in case any algorithm further will use this fact (Subtract currently would not use this fact)
		//but may use by retrieving lower_bound/higher_bound O(log(n)) operations to narrow performed operations
		erase->BulkDF(SetComparatorTag()) = ElementSet::HANDLE_COMPARATOR;
		*/

		ENTER_BLOCK();
		//for(integer jt = 0; jt < LastLocalID(ESET); ++jt) if( isValidElementSet(jt) )
		/*
		for(Mesh::iteratorSet it = BeginSet(); it != EndSet(); it++)
		{
			//ElementSet it = EsetByLocalID(jt);
			//std::cout << "set name: " << it->GetName() << " size " << it->Size() << " id " << it->LocalID() << std::endl;
			if( it->HaveParent() && it->GetParent()->Hidden() )
				it->GetParent()->RemChild(it->self());
			while( it->HaveChild() && it->GetChild()->Hidden() )
				it->RemChild(it->GetChild());
			while( it->HaveSibling() && it->GetSibling()->Hidden() )
				it->RemSibling(it->GetSibling());
			temp_hide_element = hide_element;
			hide_element = 0;
			ElementSet::iterator jt = it->Begin();
			//int q = 0;
			while(jt != it->End() )
			{
				//++q;
				//std::cout << "check element " << ElementTypeName(jt->GetElementType()) << " num " << jt->LocalID() << " handle " << jt->GetHandle() << std::endl;
				if( jt->GetMarker(temp_hide_element) ) 
				{
					//std::cout << "erase element " << ElementTypeName(jt->GetElementType()) << " num " << jt->LocalID() << " handle " << jt->GetHandle() << std::endl;
					jt = it->Erase(jt);
				}
				else ++jt;
			}
			//std::cout << "size " << it->Size() << " traversed " << q << std::endl;
			hide_element = temp_hide_element;
			//it->Subtract(erase); //old approach
		}
		*/
		EXIT_BLOCK();
		//ENTER_BLOCK();
		//SwapModification(false);
		//SwapModification(false);
		//EXIT_BLOCK();
		/*
		ENTER_BLOCK();
#if defined(USE_PARALLEL_STORAGE)
		for(parallel_storage::iterator it = shared_elements.begin(); it != shared_elements.end(); it++)
			for(int i = 0; i < 5; i++)
			{
				unsigned k = 0, l;
				for(l = 0; l < it->second[i].size(); ++l)
				{
					if( !GetMarker(it->second[i][l],hide_element) )
						it->second[i][k++] = it->second[i][l];
				}
				it->second[i].resize(k);
			}
		for(parallel_storage::iterator it = ghost_elements.begin(); it != ghost_elements.end(); it++)
			for(int i = 0; i < 5; i++)
			{
				unsigned k = 0, l;
				for(l = 0; l < it->second[i].size(); ++l)
				{
					if( !GetMarker(it->second[i][l],hide_element) )
						it->second[i][k++] = it->second[i][l];
				}
				it->second[i].resize(k);
			}
#endif
		EXIT_BLOCK();
		*/
#if 0//defined(USE_MPI)
		ENTER_BLOCK();
		SwapModification(false);
		//RecomputeParallelStorage(ESET | CELL | FACE | EDGE | NODE);
		CheckGhostSharedCount(__FILE__, __LINE__);
		MarkerType del = NewMarker();
		int mpirank = GetProcessorRank();
		Tag tag_delete = CreateTag("TEMPORARY_DELETE_GHOST_ELEMENTS_TAG", DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE, ESET | CELL | FACE | EDGE | NODE);
#if defined(USE_PARALLEL_STORAGE)
		proc_elements del_shared, del_ghost;
#endif //USE_PARALLEL_STORAGE

		for (ElementType mask = ESET; mask >= NODE; mask = PrevElementType(mask))
		{
			INMOST_DATA_ENUM_TYPE cnt = 0;
			ENTER_BLOCK();
			for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				if (it->GetMarker(del) && (GetHandleElementType(*it) & mask) && GetStatus(*it) == Element::Ghost)
				{
					Integer(*it, tag_delete) = mpirank;
					cnt++;
				}
			EXIT_BLOCK();
			if (mask & (ESET | CELL))
			{
				cnt = Integrate(cnt);
				if (!cnt) continue;
			}
			ENTER_BLOCK();
			if (mask & (FACE | EDGE | NODE))
			{
				for(Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
				{
					if (it->Hidden()) continue;
					Element::Status estat = it->GetStatus();
					if (estat == Element::Owned || estat == Element::Shared) continue;
					if (it->nbAdjElements(NextElementType(mask)) == 0) it->Integer(tag_delete) = mpirank;
				}
			}
			EXIT_BLOCK();
			ReduceData(tag_delete, mask, 0, DeleteUnpack);
			ExchangeData(tag_delete, mask, 0);
			ENTER_BLOCK();
			for (Mesh::iteratorElement it = BeginElement(mask); it != EndElement(); ++it)
			{
				if (it->Hidden()) continue;
				Element::Status estat = it->GetStatus();
				if (estat == Element::Owned) continue;
				if (it->HaveData(tag_delete))
				{
					Storage::integer_array del_procs = it->IntegerArray(tag_delete);
					std::sort(del_procs.begin(), del_procs.end());

					if (estat == Element::Ghost && std::binary_search(del_procs.begin(), del_procs.end(), mpirank))
					{
#if defined(USE_PARALLEL_STORAGE)
						del_ghost[it->IntegerDF(tag_owner)].push_back(it->GetHandle());
#endif //USE_PARALLEL_STORAGE
					}
					else
					{
						Storage::integer_array procs = it->IntegerArrayDV(tag_processors);
						std::vector<Storage::integer> result(procs.size());
#if defined(USE_PARALLEL_STORAGE)
						if (estat == Element::Shared)
						{
							for (Storage::integer_array::iterator vit = del_procs.begin(); vit != del_procs.end(); vit++)
								del_shared[*vit].push_back(it->GetHandle());
						}
#endif	//USE_PARALLEL_STORAGE
						std::vector<Storage::integer>::iterator end = std::set_difference(procs.begin(), procs.end(), del_procs.begin(), del_procs.end(), result.begin());
						result.resize(end - result.begin());
						procs.clear();
						procs.insert(procs.begin(), result.begin(), result.end());

						if (procs.size() == 1 && procs[0] == mpirank)
							it->SetStatus(Element::Owned);
					}
				}
			}
			EXIT_BLOCK();
			ENTER_BLOCK();
#if defined(USE_PARALLEL_STORAGE)
			for (proc_elements::iterator it = del_ghost.begin(); it != del_ghost.end(); it++)
			{
				//std::cout << GetProcessorRank() << " ghost delete size " << it->second.size() << std::endl;
				element_set& ref = ghost_elements[it->first][ElementNum(mask)];
				if (!it->second.empty())
				{
					if (HaveGlobalID(mask))
						std::sort(it->second.begin(), it->second.end(), GlobalIDComparator(this));
					else
						std::sort(it->second.begin(), it->second.end(), CentroidComparator(this));
				}
				element_set result(ref.size());
				element_set::iterator end;
				if (HaveGlobalID(mask))
					end = std::set_difference(ref.begin(), ref.end(), it->second.begin(), it->second.end(), result.begin(), GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(), ref.end(), it->second.begin(), it->second.end(), result.begin(), CentroidComparator(this));
				result.resize(end - result.begin());
				ref.swap(result);
			}
			del_ghost.clear();
			for (proc_elements::iterator it = del_shared.begin(); it != del_shared.end(); it++)
			{
				//std::cout << GetProcessorRank() << " shared delete size " << it->second.size() << std::endl;
				element_set& ref = shared_elements[it->first][ElementNum(mask)];
				if (!it->second.empty())
				{
					if (HaveGlobalID(mask))
						std::sort(it->second.begin(), it->second.end(), GlobalIDComparator(this));
					else
						std::sort(it->second.begin(), it->second.end(), CentroidComparator(this));
				}
				element_set result(ref.size());
				element_set::iterator end;
				if (HaveGlobalID(mask))
					end = std::set_difference(ref.begin(), ref.end(), it->second.begin(), it->second.end(), result.begin(), GlobalIDComparator(this));
				else
					end = std::set_difference(ref.begin(), ref.end(), it->second.begin(), it->second.end(), result.begin(), CentroidComparator(this));
				result.resize(end - result.begin());
				ref.swap(result);
			}
			del_shared.clear();
#endif //USE_PARALLEL_STORAGE
			EXIT_BLOCK();
		}
		DeleteTag(tag_delete);

		SwapModification(false);
		EXIT_BLOCK();
#endif //USE_MPI
		/*
		ENTER_BLOCK();
		if (UpdateGeometryMarker())
		{
			for (ElementType etype = EDGE; etype <= CELL; etype = NextElementType(etype))
			{
				int updated = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:updated)
#endif
				for (integer it = 0; it < LastLocalID(etype); ++it) if (isValidElement(etype, it))
				{
					Element e = ElementByLocalID(etype, it);
					if (e.GetMarker(UpdateGeometryMarker()))
					{
						updated++;
						RecomputeGeometricData(e.GetHandle());
					}
				}
				std::cout << "Geometry updated for " << updated << " element type " << ElementTypeName(etype) << std::endl;
			}
			ReleaseMarker(UpdateGeometryMarker());
			update_geometry = 0;
		}
		EXIT_BLOCK();
		*/
		//Destroy(erase);//old approach
		EXIT_FUNC();
	}
	void Mesh::ResolveModification()
	{
		if( GetProcessorsNumber() == 1 ) return;
		ENTER_FUNC();
		//ReportParallelStorage();
		//CheckCentroids(__FILE__,__LINE__);
		/*
		int h = 0, n = 0, hn = 0;
		for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
		for(int k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k) )
		{
			Element it = ElementByLocalID(etype,k);
			if( it->New() ) n++;
			if( it->Hidden() ) h++;
			if( it->New() && it->Hidden() ) hn++;
		}
		std::cout << GetProcessorRank() << " before resolve shared new " << n << " hidden " << h << " both " << hn << std::endl;
		*/
		ENTER_BLOCK();
		CheckSetLinks(__FILE__,__LINE__);
		ResolveSets();
		CheckSetLinks(__FILE__,__LINE__);
		ResolveShared(true);
		//ResolveShared();
		CheckSetLinks(__FILE__,__LINE__);
		EXIT_BLOCK();
		//ReportParallelStorage();
		//RecomputeParallelStorage(ESET | CELL | FACE | EDGE | NODE);
		CheckCentroids(__FILE__,__LINE__);
		/*
		h = 0, n = 0, hn = 0;
		for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
		for(int k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k) )
		{
			Element it = ElementByLocalID(etype,k);
			if( it->New() ) n++;
			if( it->Hidden() ) h++;
			if( it->New() && it->Hidden() ) hn++;
		}
		std::cout << GetProcessorRank() << " before exchange ghost new " << n << " hidden " << h << " both " << hn << std::endl;
		*/
		CheckGhostSharedCount(__FILE__, __LINE__);
		//std::cout << "layers " << Integer(GetHandle(),tag_layers) << " bridge " << ElementTypeName(ElementType(Integer(GetHandle(),tag_bridge))) << std::endl;
		ENTER_BLOCK();
		if( Integer(GetHandle(),tag_layers) )
			ExchangeGhost(Integer(GetHandle(),tag_layers),Integer(GetHandle(),tag_bridge));//,NewMarker()); //TODO!!!!
		EXIT_BLOCK();

		CheckGhostSharedCount(__FILE__, __LINE__);
		//ReportParallelStorage();
		//CheckCentroids(__FILE__,__LINE__);
		//Save("after_exchange_ghost.pvtk");
		/*
		h = 0, n = 0, hn = 0;
		for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
		for(int k = 0; k < LastLocalID(etype); ++k) if( isValidElement(etype,k) )
		{
			Element it = ElementByLocalID(etype,k);
			if( it->New() ) n++;
			if( it->Hidden() ) h++;
			if( it->New() && it->Hidden() ) hn++;
		}
		std::cout << GetProcessorRank() << " after exchange ghost new " << n << " hidden " << h << " both " << hn << std::endl;
		*/
		
		//ReportParallelStorage();
		//CheckCentroids(__FILE__,__LINE__);
		//exit(-1);
		EXIT_FUNC();
	}
	
	void Mesh::EndModification()
	{
		ENTER_FUNC();
		//ApplyModification();
		//temp_hide_element = hide_element;
		//hide_element = 0;
		/*
		for(ElementType etype = FACE; etype >= NODE; etype = PrevElementType(etype))
		{
			for(integer it = 0; it < LastLocalID(etype); ++it) if( isValidElement(etype,it) )
			{
				//all upper elements are deleted
				if( ElementByLocalID(etype,it).nbAdjElements(NextElementType(etype),hide_element) ==
				    ElementByLocalID(etype,it).nbAdjElements(NextElementType(etype)) )
					SetMarker(ComposeHandle(etype,it),hide_element);
			}
		}
		 */
		MarkerType nm = new_element;
		ENTER_BLOCK();
		new_element = 0;
		for(ElementType etype = ESET; etype >= NODE; etype = PrevElementType(etype))
		{
			for(integer it = 0; it < LastLocalID(etype); ++it) if( isValidElement(etype,it) )
			{
				HandleType h = ComposeHandle(etype,it);
				RemMarker(h,nm);
				if( GetMarker(h,hide_element) )
					Destroy(h);
			}
		}
		new_element = nm;
		EXIT_BLOCK();
		/*
		for(ElementType etype = FACE; etype >= NODE; etype = PrevElementType(etype))
		{
			for(integer it = 0; it < LastLocalID(etype); ++it) if( isValidElement(etype,it) )
			{
				if( ElementByLocalID(etype,it).nbAdjElements(NextElementType(etype)) == 0 )
					Destroy(ComposeHandle(etype,it));
			}
		}
		 */
		//RecomputeParallelStorage(ESET|CELL|FACE|EDGE|NODE);
		//ExchangeGhost(Integer(GetHandle(),tag_layers),Integer(GetHandle(),tag_bridge),NewMarker()); //TODO!!!!
		memset(hidden_count,0,sizeof(integer)*6);
		memset(hidden_count_zero,0,sizeof(integer)*6);
		ReleaseMarker(hide_element);
		ReleaseMarker(new_element);
		hide_element = 0;
		new_element = 0;
		//This should be done in ResolveModification
		ElementType have_global_id = NONE;
		ENTER_BLOCK();
		if( GlobalIDTag().isValid() )
		{
			for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
				if( GlobalIDTag().isDefined(etype) ) have_global_id |= etype;
		}
		EXIT_BLOCK();
		if( have_global_id ) AssignGlobalID(have_global_id);
		EXIT_FUNC();
	}


	void Mesh::Destroy(HandleType h)
	{
		//if( h == 1073743552 || h == 1073742933)
		//	std::cout << GetProcessorRank() << " destroy " << ElementTypeName(GetHandleElementType(h)) << " id " << GetHandleID(h) << " handle " << h << std::endl;// << (GetMarker(h,temp_hide_element) ? " hidden " : " not hidden ") << std::endl;
		assert(isValidHandleRange(h));
		assert(isValidHandle(h));
		ElementType htype = GetHandleElementType(h);
		//if( Hidden(h) ) Show(h);
		if( htype & (NODE|EDGE|FACE|CELL) )
			ElementByHandle(h).Disconnect(true);
		else if( htype & ESET )
		{
			ElementSet eset(this,h);
			set_search.erase(eset->GetName());
			if( eset->HaveParent() )
				eset->GetParent()->RemChild(eset);
			while( eset->HaveChild() )
				eset->RemChild(eset->GetChild());
		}
		for(iteratorTag t = BeginTag(); t != EndTag(); ++t) 
			if( t->isDefined(htype) ) DelData(h,*t);
		UntieElement(GetHandleElementNum(h),GetHandleID(h));
	}

	bool Mesh::Hide(HandleType h)
	{
		if(HideMarker()) 
		{
			if( !Hidden(h) )
			{
				//TODO if we hide a back cell for a face we have to swap it's normal!
				if( GetHandleElementType(h) == CELL ) Cell(this,h)->SwapBackCell();
				hidden_count[GetHandleElementNum(h)]++;
				SetMarker(h,HideMarker()); 
			}
			return true;
		} 
		return false;
	}
	
	bool Mesh::Show(HandleType h)
	{
		if(HideMarker()) 
		{
			if( Hidden(h) )
			{
				hidden_count[GetHandleElementNum(h)]--;
				RemMarker(h,HideMarker()); 
			}
			return true;
		} 
		return false;
	}
	
	bool Mesh::Delete(HandleType h) 
	{
		if(/*!New(h) &&*/ Hide(h))
		{
			//mark all elements that rely on this that they should be deleted
			if( GetHandleElementType(h) < CELL )
			{
				Element::adj_type & hc = HighConn(h);
				for(Element::adj_type::size_type it = 0; it < hc.size(); ++it)
					if( !Hidden(hc[it]) ) Delete(hc[it]);
			}
			return false;
		}
		else 
		{
			Destroy(h);
			return true;
		}
	}
	
	bool Mesh::Hidden(HandleType h) const
	{
		assert(h != InvalidHandle());
		return GetMarker(h,HideMarker());
	}
	
	bool Mesh::New(HandleType h) const
	{
		assert(h != InvalidHandle());
		return GetMarker(h,NewMarker());
	}
	
	bool Element::Hide() const
	{
		return GetMeshLink()->Hide(GetHandle());
	}
	
	bool Element::Show() const
	{
		return GetMeshLink()->Show(GetHandle());
	}
	
	bool Element::Delete() 
	{
		bool ret = GetMeshLink()->Delete(GetHandle());
		if( ret ) handle = InvalidHandle();
		return ret;
	}
	
	bool Element::Hidden() const
	{
		return GetMeshLink()->Hidden(GetHandle());
	}
	
	bool Element::New() const
	{
		return GetMeshLink()->New(GetHandle());
	}
	
	
	Storage::enumerator Mesh::getNext(const HandleType * arr, enumerator size, enumerator k, MarkerType marker) const
	{
		k++;
		while(k < size && GetMarker(arr[k],marker)) k++;
		return k;
	}
	
	Storage::enumerator Mesh::Count(const HandleType * arr, enumerator size, MarkerType marker) const
	{
		enumerator ret = 0, k = 0;
		while(k < size) 
		{
			if( !GetMarker(arr[k],marker) ) ret++;
			k++;
		}
		return ret;
	}

	//void Element::UpdateGeometricData() const
	//{
	//	GetMeshLink()->RecomputeGeometricData(GetHandle());
	//}
	
	void Cell::SwapBackCell() const
	{
		Mesh * m = GetMeshLink();
		
		//if( m->HaveGeometricData(ORIENTATION,FACE) )
		{
			//retrieve faces
			MarkerType hm = m->HideMarker();
			adj_type & lc = m->LowConn(GetHandle());
			for( adj_type::size_type it = 0; it < lc.size(); it++) if( !m->GetMarker(lc[it],hm) ) //iterator over unhidden faces
			{
				adj_type & hc = m->HighConn(lc[it]); //cells of face
				if( !hc.empty() && m->Count(hc.data(),static_cast<enumerator>(hc.size()),hm) == 2 ) //there are two cells connected to face
				{
					enumerator k1 = ENUMUNDEF, k2;
					k1 = m->getNext(hc.data(),static_cast<enumerator>(hc.size()),k1,hm);
					if( hc[k1] == GetHandle() ) //the first cell is current
					{
						k2 = m->getNext(hc.data(),static_cast<enumerator>(hc.size()),k1,hm);
						std::swap(hc[k1],hc[k2]);
						//hc[k2] = GetHandle(); //cannot use the cell because virtualization table is already destroyed and FixNormalOrientation will do bad things
						//hc.resize(1); //just remove element, we will do this anyway later
						if (m->HaveGeometricData(ORIENTATION, FACE))
							Face(m, lc[it])->FixNormalOrientation(false); //restore orientation
					}
				}
				//~ if( !Face(m,lc[it])->CheckNormalOrientation() )
				//~ {
					//~ std::cout << __FILE__ << ":" << __LINE__ << " on " << GetMeshLink()->GetProcessorRank() << " face " << lc[it] << " of cell " << GetHandle() << " bad orientation, count " << m->Count(hc.data(),static_cast<enumerator>(hc.size()),hm) << " hc " << hc.size() << std::endl;
				//~ }
			}
		}
	}
	
	void Element::Disconnect(const HandleType * adjacent, INMOST_DATA_ENUM_TYPE num) const
	{
		Mesh * m = GetMeshLink();
		INMOST_DATA_ENUM_TYPE k = 0;
		std::vector<HandleType> arr[4];
		MarkerType mrk = m->CreatePrivateMarker(), mod = m->CreatePrivateMarker();
		for(k = 0; k < num; k++) m->SetPrivateMarker(adjacent[k], mrk);
		adj_iterator it;
		// disconnects nodes from current edge, edges from current face, faces from current cell
		if( GetElementType() > NODE ) //cannot disconnect cells from current node
		{
			adj_type & lc = m->LowConn(GetHandle());
			it = lc.begin();
			//look into nodes (edge), edges (face), faces (cell)
			while(it != lc.end() ) 
			{
				if( m->GetPrivateMarker(*it,mrk) ) //element should be disconnected
				{
					adj_type & hc = m->HighConn(*it); // edges (node), faces (edge), cells (face)
					if( GetElementType() == CELL ) //update some geometric data
					{
						MarkerType hm = m->HideMarker();
						if( !hc.empty() && m->Count(hc.data(),static_cast<enumerator>(hc.size()),hm) == 2 )
						{
							enumerator k1 = ENUMUNDEF, k2;
							k1 = m->getNext(hc.data(),static_cast<enumerator>(hc.size()),k1,hm);
							if( hc[k1] == GetHandle() ) //the first cell is current
							{
								k2 = m->getNext(hc.data(),static_cast<enumerator>(hc.size()),k1,hm);
								hc[k1] = hc[k2];
								hc[k2] = GetHandle();
								Face(m,*it)->FixNormalOrientation(false); //restore orientation
							}
						}
					}
					//look into edges of node (edge) faces of edge (face) cells of face (cell)
					for(adj_iterator jt = hc.begin(); jt != hc.end(); ++jt)
					{
						if( *jt == GetHandle() ) //remove me
						{
							hc.erase(jt);
							//don't update me
							//only lower adjacencies affect geometric data of element
							//this may be subject to change
							break;
						}
					}
					//update my data
					if( !m->GetPrivateMarker(GetHandle(),mod) ) 
					{
						m->SetPrivateMarker(GetHandle(),mod);
						arr[GetElementNum()].push_back(GetHandle());
					}
					//disconnect element
					it = lc.erase(it);
				}
				else it++;
			}
		}
		// disconnect edges from current node, faces from current edge, cells from current face
		if( GetElementType() < CELL ) //Cell cannot be disconnected from nodes
		{
			adj_type & hc = m->HighConn(GetHandle());
			it = hc.begin();
			//go over edges (node), faces (edge), cells (face)
			while(it != hc.end() ) 
			{
				if( m->GetPrivateMarker(*it,mrk) ) //should be disconnected
				{
					adj_type & lc = m->LowConn(*it);
					//look into nodes of edge (node), edges of face (edge), faces of cell (face)
					for(adj_iterator jt = lc.begin(); jt != lc.end(); ++jt)
					{
						if( *jt == GetHandle() ) //remove me
						{
							lc.erase(jt);
							// lower adjacencies of edge (node), face (edge), cell (face) where modified
							// update it's geometric data
							if( !m->GetPrivateMarker(*it,mod) ) 
							{
								m->SetPrivateMarker(*it,mod);
								arr[GetHandleElementNum(*it)].push_back(*it);
							}
							break;
						}
					}
					/* //don't update me
					//only lower adjacencies affect geometric data of element
					//this may be subject to change
					if( !GetMarker(mod) )
					{
						SetMarker(mod);
						arr[GetElementNum()].push_back(this);
					}
					*/
					//disconnect element
					it = hc.erase(it);
				}
				else it++;
			}
		}
		for(k = 0; k < num; k++) m->RemPrivateMarker(adjacent[k],mrk);
		//if element placed below on hierarchy was modified, then all upper elements
		//should be also modified, start from lowest elements in hierarchy and go up
		for(ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype) )
		{
			int el_num = ElementNum(etype);
			for (size_t it = 0; it < arr[el_num].size(); it++)
			{
				assert(GetHandleElementType(arr[el_num][it]) == etype);
				if (etype < CELL) //check for upper adjacencies of current element
				{
					//check all upper adjacencies that may be affected by modification of me
					adj_type& hc = m->HighConn(arr[el_num][it]);
					for (adj_type::size_type jt = 0; jt < hc.size(); ++jt)
					{
						if (!m->GetPrivateMarker(hc[jt], mod))
						{
							m->SetPrivateMarker(hc[jt], mod);
							arr[GetHandleElementNum(hc[jt])].push_back(hc[jt]);
						}
					}
				}
				else if( m->HighConnTag().isDefined(CELL) && m->LowConnTag().isDefined(NODE) ) //update nodes for current cell
				{
					//mark remaining connected nodes
					ElementArray<Node> cnodes(m);
					MarkerType keep = m->CreatePrivateMarker();
					adj_type& lc = m->LowConn(arr[el_num][it]); //cell faces
					for (adj_type::iterator ilc = lc.begin(); ilc != lc.end(); ++ilc) //loop remaining faces
					{
						adj_type& elc = m->LowConn(*ilc);//edges of the face
						for (adj_type::iterator jlc = elc.begin(); jlc != elc.end(); ++jlc) //loop edges
						{
							adj_type& nlc = m->LowConn(*jlc);//nodes of the edge
							for (adj_type::iterator klc = nlc.begin(); klc != nlc.end(); ++klc) if(!m->GetPrivateMarker(*klc,keep))
							{
								cnodes.push_back(*klc);
								m->SetPrivateMarker(*klc, keep);
							}
						}
					}
					//remove me from disconnected nodes
					adj_type& hc = m->HighConn(arr[el_num][it]);
					adj_iterator jt = hc.begin();
					while(jt != hc.end()) //iterate over my nodes
					{
						if (!m->GetPrivateMarker(*jt, keep))
						{
							adj_type& lc = m->LowConn(*jt);
							adj_iterator kt = lc.begin();
							while (kt != lc.end()) //iterate over nodes's cells
							{
								if ((*kt) == arr[el_num][it]) // if nodes's cell is equal to modified cell, then remove connection
								{
									kt = lc.erase(kt);
									break;
								}
								else kt++;
							}
							jt = hc.erase(jt);
						}
						else jt++;
					}
					//check that we have an extra node
					cnodes.RemPrivateMarker(keep);
					m->SetPrivateMarkerArray(hc.data(), hc.size(), keep);
					for (ElementArray<Node>::iterator jt = cnodes.begin(); jt != cnodes.end(); ++jt) if (!jt->GetPrivateMarker(keep))
					{
						hc.push_back(*jt);
						//connect cell to node (if not)
						adj_type& nlc = m->LowConn(*jt); 
						bool have_cell = false;
						for (adj_type::iterator kt = nlc.begin(); kt != nlc.end() && !have_cell; ++kt)
							if (*kt == arr[el_num][it]) have_cell = true;
						if (!have_cell)
							nlc.push_back(arr[el_num][it]);
						jt->SetPrivateMarker(keep);
					}
					m->RemPrivateMarkerArray(hc.data(), hc.size(), keep);
					m->ReleasePrivateMarker(keep);
				}
			}
			if (etype > NODE)
			{
				for (size_t it = 0; it < arr[el_num].size(); it++) m->ComputeGeometricType(arr[el_num][it]);
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], CENTROID);
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], NORMAL);
				if (etype >= FACE)
				{
					for (size_t it = 0; it < arr[el_num].size(); it++)
						m->RecomputeGeometricData(arr[el_num][it], ORIENTATION);
				}
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], MEASURE);
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], BARYCENTER);
			}
			if (!arr[el_num].empty()) m->RemPrivateMarkerArray(&arr[el_num][0], (enumerator)arr[el_num].size(), mod);
		}
		m->ReleasePrivateMarker(mod);
		m->ReleasePrivateMarker(mrk);
	}

	void Element::Connect(const HandleType * adjacent, INMOST_DATA_ENUM_TYPE num) const
	{
		assert( !(GetElementType() == EDGE && GetMeshLink()->LowConn(GetHandle()).size() > 2) ); // cannot add another node to edge
		Mesh * m = GetMeshLink();
		std::vector<HandleType> arr[4];
		MarkerType mod = m->CreatePrivateMarker();
		for(INMOST_DATA_ENUM_TYPE k = 0; k < num; k++)
		{
			assert( GetHandleElementType(adjacent[k]) != CELL );
			assert( GetHandleElementType(adjacent[k]) == PrevElementType(GetElementType()) ); //only lower dimension elements can be connected
			assert( !(GetHandleElementType(adjacent[k]) == FACE && m->HighConn(adjacent[k]).size() > 1) ); // face already connected to two cells

			m->HighConn(adjacent[k]).push_back(GetHandle());
			//this may be dead slow
			//LowConn().push_back(adjacent[k]);
			if (!m->GetPrivateMarker(adjacent[k], mod))
			{
				m->SetPrivateMarker(adjacent[k], mod);
				arr[GetHandleElementNum(adjacent[k])].push_back(adjacent[k]);
			}
		}
		adj_type & lc = m->LowConn(GetHandle());
		//TODO:
		//for face have to find positions where to attach edges
		lc.insert(lc.end(),adjacent,adjacent+num);
		if (!m->GetPrivateMarker(GetHandle(), mod))
		{
			m->SetPrivateMarker(GetHandle(), mod);
			arr[GetElementNum()].push_back(GetHandle());
		}
		//perform fix instead of algorithm above
		if( Element::GetGeometricDimension(m->GetGeometricType(GetHandle())) == 2 )
		{
			if( GetElementType() == CELL ) 
			{
				getAsCell()->FixEdgeOrder();
			}
			else if( GetElementType() == FACE ) 
			{
				if( getAsFace()->FixEdgeOrder() )
					getAsFace()->FixNormalOrientation();
			}
		}
		ComputeGeometricType();
		//if (m->UpdateGeometryMarker())
		//{
		//	m->SetMarker(GetHandle(), m->UpdateGeometryMarker());
		//	if (GetElementType() == FACE)
		//		Face(m, GetHandle()).FixNormalOrientation();
		//}
		//else 
		m->RecomputeGeometricData(GetHandle());
		//if element placed below on hierarchy was modified, then all upper elements
		//should be also modified, start from lowest elements in hierarchy and go up
		for (ElementType etype = NODE; etype <= CELL; etype = NextElementType(etype))
		{
			int el_num = ElementNum(etype);
			for (size_t it = 0; it < arr[el_num].size(); it++)
			{
				assert(GetHandleElementType(arr[el_num][it]) == etype);
				if (etype < CELL) //check for upper adjacencies of current element
				{
					//check all upper adjacencies that may be affected by modification of me
					adj_type& hc = m->HighConn(arr[el_num][it]);
					for (adj_type::size_type jt = 0; jt < hc.size(); ++jt)
					{
						if (!m->GetPrivateMarker(hc[jt], mod))
						{
							m->SetPrivateMarker(hc[jt], mod);
							arr[GetHandleElementNum(hc[jt])].push_back(hc[jt]);
						}
					}
				}
				else if (m->HighConnTag().isDefined(CELL) && m->LowConnTag().isDefined(NODE)) //update nodes for current cell
				{
					//mark remaining connected nodes
					ElementArray<Node> cnodes(m);
					MarkerType keep = m->CreatePrivateMarker();
					adj_type& lc = m->LowConn(arr[el_num][it]); //cell faces
					for (adj_type::iterator ilc = lc.begin(); ilc != lc.end(); ++ilc) //loop remaining faces
					{
						adj_type& elc = m->LowConn(*ilc);//edges of the face
						for (adj_type::iterator jlc = elc.begin(); jlc != elc.end(); ++jlc) //loop edges
						{
							adj_type& nlc = m->LowConn(*jlc);//nodes of the edge
							for (adj_type::iterator klc = nlc.begin(); klc != nlc.end(); ++klc) if (!m->GetPrivateMarker(*klc, keep))
							{
								cnodes.push_back(*klc);
								m->SetPrivateMarker(*klc, keep);
							}
						}
					}
					//remove me from disconnected nodes
					adj_type& hc = m->HighConn(arr[el_num][it]);
					adj_iterator jt = hc.begin();
					while (jt != hc.end()) //iterate over my nodes
					{
						if (!m->GetPrivateMarker(*jt, keep))
						{
							adj_type& lc = m->LowConn(*jt);
							adj_iterator kt = lc.begin();
							while (kt != lc.end()) //iterate over nodes's cells
							{
								if ((*kt) == arr[el_num][it]) // if nodes's cell is equal to modified cell, then remove connection
								{
									kt = lc.erase(kt);
									break;
								}
								else kt++;
							}
							jt = hc.erase(jt);
						}
						else jt++;
					}
					//check that we have an extra node
					cnodes.RemPrivateMarker(keep);
					m->SetPrivateMarkerArray(hc.data(), hc.size(), keep);
					for (ElementArray<Node>::iterator jt = cnodes.begin(); jt != cnodes.end(); ++jt) if (!jt->GetPrivateMarker(keep))
					{
						hc.push_back(*jt);
						//connect cell to node (if not)
						adj_type& nlc = m->LowConn(*jt);
						bool have_cell = false;
						for (adj_type::iterator kt = nlc.begin(); kt != nlc.end() && !have_cell; ++kt)
							if (*kt == arr[el_num][it]) have_cell = true;
						if (!have_cell)
							nlc.push_back(arr[el_num][it]);
						jt->SetPrivateMarker(keep);
					}
					m->RemPrivateMarkerArray(hc.data(), hc.size(), keep);
					m->ReleasePrivateMarker(keep);
				}
			}
			if (etype > NODE)
			{
				for (size_t it = 0; it < arr[el_num].size(); it++) m->ComputeGeometricType(arr[el_num][it]);
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], CENTROID);
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], NORMAL);
				if (etype >= FACE)
				{
					for (size_t it = 0; it < arr[el_num].size(); it++)
						m->RecomputeGeometricData(arr[el_num][it], ORIENTATION);
				}
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], MEASURE);
				for (size_t it = 0; it < arr[el_num].size(); it++) m->RecomputeGeometricData(arr[el_num][it], BARYCENTER);
			}
			if (!arr[el_num].empty())m->RemPrivateMarkerArray(&arr[el_num][0], (enumerator)arr[el_num].size(), mod);
		}
		m->ReleasePrivateMarker(mod);
	}
}




#endif
