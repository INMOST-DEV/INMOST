#include "inmost.h"
#include <sstream>
#if defined(USE_MESH)
#define WAITNL 	{char c;scanf("%c",&c);}

//#define HEAVY_DUPLICATE_CHECK //depricated, left for future debug

#include <unordered_map>

namespace INMOST
{
	int CompareElementsCPointer(const void * pa, const void * pb);
	
	//~ void Mesh::GetElementsBySparseTag(Tag t, ElementType types, std::vector<Element *> & ret) const
	//~ {
		//~ for(ElementType mask = NODE; mask <= CELL; mask = mask << 1)
			//~ if( mask & types && t.isDefined(mask) && t.isSparse(mask) )
			//~ {
				//~ TagManager::sparse_sub_type & s = t.GetTagManager()->GetSparseData(t.GetPosition(mask));
				//~ for(TagManager::sparse_sub_type::iterator it = s.begin(); it != s.end(); ++it)
					//~ ret.push_back(reinterpret_cast<Element *>(it->first));
			//~ }
	//~ }
	
	Mesh::Mesh()
	:TagManager(),Storage(this,MESH)
	//for chunk_array
	,
	cells(1),empty_cells(1),
	faces(1),empty_faces(1),
	edges(1),empty_edges(1),
	nodes(1),empty_nodes(1),
	sets(1), empty_sets(1)
	{
		
		
		
		dim = 3;
		have_global_id = NONE;
		checkset = DEFAULT_CHECK;
		errorset = 0;
		new_element = hide_element = 0;

		memset(remember,0,sizeof(remember));
		
		//~ nodes.reserve(2048);
		//~ edges.reserve(4096);
		//~ faces.reserve(4096);
		//~ cells.reserve(1024);
		
		//tag_global_id = CreateTag("GLOBAL_ID",DATA_INTEGER, ESET | CELL | FACE | EDGE | NODE,NONE,1);
		tag_coords        = CreateTag("COORD",DATA_REAL, NODE,NONE,dim);
		tag_topologyerror = CreateTag("TOPOLOGY_ERROR_TAG",DATA_INTEGER,CELL | FACE | EDGE,CELL | FACE | EDGE,1);
#if defined(NEW_CONNECTIONS)
		tag_high_conn     = CreateTag("HIGH_CONN",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
		tag_low_conn      = CreateTag("LOW_CONN",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
#endif
#if defined(NEW_MARKERS)
		tag_markers       = CreateTag("MARKERS",DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFields);
#endif
		tag_geom_type     = CreateTag("GEOM_TYPE",DATA_BULK,CELL|FACE|EDGE|NODE,NONE,1);
		//tag_shared_elems = CreateTag("SHARED_ELEMS_TAG",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);

		epsilon = 1.0e-8;
		m_state = Mesh::Serial;

#if defined(USE_MPI)
		{
			int test;
			MPI_Initialized(&test);
			if( test == 0 ) MPI_Init(NULL,NULL);
			comm = INMOST_MPI_COMM_WORLD;
		}
#endif

#if defined(USE_PARALLEL_WRITE_TIME)
		num_exchanges = 0;
		std::stringstream temp;
		temp << "time_" << GetProcessorRank() << ".xml";
		out_time.open(temp.str().c_str(),std::ios::out);
		out_time << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl;
		out_time << "<?xml-stylesheet type=\"text/xsl\" href=\"style.xsl\"?>" << std::endl;
		out_time << "<Debug>" << std::endl;
		tab = 1;
		func_id = 0;
#endif
	}
	
	Mesh::Mesh(const Mesh & other)
	:TagManager(this,other),Storage(this,0,other),
	empty_cells(other.empty_cells),empty_faces(other.empty_faces),
	empty_edges(other.empty_edges),empty_nodes(other.empty_nodes),
	empty_sets(other.empty_sets)
	,
	cells(1),
	faces(1),
	edges(1),
	nodes(1),
	sets(1) 
	{
		INMOST_DATA_ENUM_TYPE i;
#if defined(NEW_MARKERS)
		Storage::bulk marker_space[MarkerFields];
		other.GetMarkerSpace(marker_space);
		SetMarkerSpace(marker_space);
#else
		markers = other.GetMarkerSpace();
#endif
		m_state = other.m_state;
		
		checkset = other.checkset;
		errorset = other.errorset;

		new_element = other.new_element;
		hide_element = other.hide_element;
		
		memcpy(remember,other.remember,sizeof(remember));
		
		dim = other.dim;

		
		tag_coords        = CreateTag("COORD",DATA_REAL, NODE,NONE,dim);
		tag_topologyerror = CreateTag("TOPOLOGY_ERROR_TAG",DATA_INTEGER,CELL | FACE | EDGE,CELL | FACE | EDGE,1);
#if defined(NEW_CONNECTIONS)
		tag_high_conn     = CreateTag("HIGH_CONN",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
		tag_low_conn      = CreateTag("LOW_CONN",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
#endif
#if defined(NEW_MARKERS)
		tag_markers       = CreateTag("MARKERS",DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFields);
#endif
		tag_geom_type     = CreateTag("GEOM_TYPE",DATA_BULK,CELL|FACE|EDGE|NODE,NONE,1);
		//tag_shared_elems = CreateTag("SHARED_ELEMS_TAG",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		
		if( m_state == Mesh::Parallel ) SetCommunicator(other.comm);
		else comm = INMOST_MPI_COMM_WORLD;
		
		RestoreGeometricTags();

		//~ nodes.reserve(other.nodes.size());
		//~ edges.reserve(other.edges.size());
		//~ faces.reserve(other.faces.size());
		//~ cells.reserve(other.cells.size());
		//~ sets.reserve(other.sets.size());
		for(ElementType etype = NODE; etype <= MESH; etype = etype << 1 )
			ReallocateData(etype);

		/*
		for(iteratorTag t = BeginTag(); t != EndTag(); t++)
			for(ElementType etype = NODE; etype <= MESH; etype = etype << 1 )
				if( t->isDefined(etype) && !t->isSparse(etype) )
					t->AllocateData(etype);
		*/
		
		for(nodes_container::const_iterator it = other.nodes.begin(); it != other.nodes.end(); it++) if( (*it) != NULL ) nodes.push_back(new Node(this,(*it)->LocalID(),*(*it))); else nodes.push_back(NULL);
		for(edges_container::const_iterator it = other.edges.begin(); it != other.edges.end(); it++) if( (*it) != NULL ) edges.push_back(new Edge(this,(*it)->LocalID(),*(*it))); else edges.push_back(NULL);
		for(faces_container::const_iterator it = other.faces.begin(); it != other.faces.end(); it++) if( (*it) != NULL ) faces.push_back(new Face(this,(*it)->LocalID(),*(*it))); else faces.push_back(NULL);
		for(cells_container::const_iterator it = other.cells.begin(); it != other.cells.end(); it++) if( (*it) != NULL ) cells.push_back(new Cell(this,(*it)->LocalID(),*(*it))); else cells.push_back(NULL);
		for(sets_container::const_iterator it = other.sets.begin(); it != other.sets.end(); it++) if( (*it) != NULL ) sets.push_back(new ElementSet(this,(*it)->LocalID(),*(*it))); else sets.push_back(NULL);
		i = 0;
		for(nodes_container::const_iterator it = other.nodes.begin(); it != other.nodes.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) nodes[i]->LowConn().push_back(cells[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) nodes[i]->HighConn().push_back(edges[(*jt)->LocalID()]);
		}
		i = 0;
		for(edges_container::const_iterator it = other.edges.begin(); it != other.edges.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) edges[i]->LowConn().push_back(nodes[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) edges[i]->HighConn().push_back(faces[(*jt)->LocalID()]);
		}
		i = 0;
		for(faces_container::const_iterator it = other.faces.begin(); it != other.faces.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) faces[i]->LowConn().push_back(edges[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) faces[i]->HighConn().push_back(cells[(*jt)->LocalID()]);
		}
		i = 0;
		for(cells_container::const_iterator it = other.cells.begin(); it != other.cells.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) cells[i]->LowConn().push_back(faces[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) cells[i]->HighConn().push_back(nodes[(*jt)->LocalID()]);
		}
		i = 0;
		for(sets_container::const_iterator it = other.sets.begin(); it != other.sets.end(); it++, i++) if( (*it) != NULL )
		{
			for(ElementSet::iterator jt = (*it)->begin(); jt != (*it)->end(); jt++)
				switch(jt->GetElementType())
				{
					case NODE: sets[i]->Insert(nodes[jt->LocalID()]); break;
					case EDGE: sets[i]->Insert(edges[jt->LocalID()]); break;
					case FACE: sets[i]->Insert(faces[jt->LocalID()]); break;
					case CELL: sets[i]->Insert(cells[jt->LocalID()]); break;
				}
		}
		epsilon = other.epsilon;
		AssignGlobalID(other.have_global_id);
	}
	
	Mesh & Mesh::operator =(Mesh const & other)
	{
		INMOST_DATA_ENUM_TYPE i;
		if( this == &other ) return *this; //don't do anything
		//first delete everithing
		
		for(sets_container::iterator it = sets.begin(); it != sets.end(); it++) if( *it != NULL ) delete *it;
		for(cells_container::iterator it = cells.begin(); it != cells.end(); it++) if( *it != NULL ) delete *it;
		for(faces_container::iterator it = faces.begin(); it != faces.end(); it++) if( *it != NULL ) delete *it;
		for(edges_container::iterator it = edges.begin(); it != edges.end(); it++) if( *it != NULL ) delete *it;
		for(nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++) if( *it != NULL ) delete *it;
#if defined(USE_MPI)
#if defined(USE_MPI2)
		if( m_state == Mesh::Parallel )
		{
			MPI_Free_mem(shared_space);
			MPI_Win_free(&window);
		}
#endif
#endif
		nodes.clear();
		edges.clear();
		faces.clear();
		cells.clear();
		sets.clear();
		empty_cells.clear();
		empty_faces.clear();
		empty_edges.clear();
		empty_nodes.clear();
		empty_sets.clear();

		
		TagManager::assign(this,other);
		Storage::assign(this,0,other);

		m_state = other.m_state;
		dim = other.dim;
		checkset = other.checkset;
		errorset = other.errorset;
		memcpy(remember,other.remember,sizeof(remember));

		
		tag_coords        = CreateTag("COORD",DATA_REAL, NODE,NONE,dim);
		tag_topologyerror = CreateTag("TOPOLOGY_ERROR_TAG",DATA_INTEGER,CELL | FACE | EDGE,CELL | FACE | EDGE,1);
#if defined(NEW_CONNECTIONS)
		tag_high_conn     = CreateTag("HIGH_CONN",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
		tag_low_conn      = CreateTag("LOW_CONN",DATA_REFERENCE,CELL|FACE|EDGE|NODE,NONE);
#endif
#if defined(NEW_MARKERS)
		tag_markers       = CreateTag("MARKERS",DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFields);
#endif
		tag_geom_type     = CreateTag("GEOM_TYPE",DATA_BULK,CELL|FACE|EDGE|NODE,NONE,1);
		//tag_shared_elems = CreateTag("SHARED_ELEMS_TAG",DATA_INTEGER,CELL|FACE|EDGE|NODE,NONE,1);
		
		
		if( m_state == Mesh::Parallel ) SetCommunicator(other.comm);
		else comm = INMOST_MPI_COMM_WORLD;

		RestoreGeometricTags();

		
		empty_cells = other.empty_cells;
		empty_faces = other.empty_faces;
		empty_edges = other.empty_edges;
		empty_nodes = other.empty_nodes;
		empty_sets = other.empty_sets; 

		//~ nodes.reserve(other.nodes.size());
		//~ edges.reserve(other.edges.size());
		//~ faces.reserve(other.faces.size());
		//~ cells.reserve(other.cells.size());
		//~ sets.reserve(other.sets.size());

		for(ElementType etype = NODE; etype <= MESH; etype = etype << 1 )
			ReallocateData(etype);
		/*
		for(iteratorTag t = BeginTag(); t != EndTag(); t++)
			for(ElementType etype = NODE; etype <= MESH; etype = etype << 1 )
				if( t->isDefined(etype) && !t->isSparse(etype) )
					t->AllocateData(etype);
		*/
		for(nodes_container::const_iterator it = other.nodes.begin(); it != other.nodes.end(); it++) if( (*it) != NULL ) nodes.push_back(new Node(this,(*it)->LocalID(),*(*it))); else nodes.push_back(NULL);
		for(edges_container::const_iterator it = other.edges.begin(); it != other.edges.end(); it++) if( (*it) != NULL ) edges.push_back(new Edge(this,(*it)->LocalID(),*(*it))); else edges.push_back(NULL);
		for(faces_container::const_iterator it = other.faces.begin(); it != other.faces.end(); it++) if( (*it) != NULL ) faces.push_back(new Face(this,(*it)->LocalID(),*(*it))); else faces.push_back(NULL);
		for(cells_container::const_iterator it = other.cells.begin(); it != other.cells.end(); it++) if( (*it) != NULL ) cells.push_back(new Cell(this,(*it)->LocalID(),*(*it))); else cells.push_back(NULL);
		for(sets_container::const_iterator it = other.sets.begin(); it != other.sets.end(); it++) if( (*it) != NULL ) sets.push_back(new ElementSet(this,(*it)->LocalID(),*(*it))); else sets.push_back(NULL);
		i = 0;
		for(nodes_container::const_iterator it = other.nodes.begin(); it != other.nodes.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) nodes[i]->LowConn().push_back(cells[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) nodes[i]->HighConn().push_back(edges[(*jt)->LocalID()]);
		}
		i = 0;
		for(edges_container::const_iterator it = other.edges.begin(); it != other.edges.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) edges[i]->LowConn().push_back(nodes[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) edges[i]->HighConn().push_back(faces[(*jt)->LocalID()]);
		}
		i = 0;
		for(faces_container::const_iterator it = other.faces.begin(); it != other.faces.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) faces[i]->LowConn().push_back(edges[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) faces[i]->HighConn().push_back(cells[(*jt)->LocalID()]);
		}
		i = 0;
		for(cells_container::const_iterator it = other.cells.begin(); it != other.cells.end(); it++, i++) if( (*it) != NULL )
		{
			for(Element::adj_iterator jt = (*it)->LowConn().begin(); jt != (*it)->LowConn().end(); jt++) cells[i]->LowConn().push_back(faces[(*jt)->LocalID()]);
			for(Element::adj_iterator jt = (*it)->HighConn().begin(); jt != (*it)->HighConn().end(); jt++) cells[i]->HighConn().push_back(nodes[(*jt)->LocalID()]);
		}
		i = 0;
		for(sets_container::const_iterator it = other.sets.begin(); it != other.sets.end(); it++, i++) if( (*it) != NULL )
		{
			for(ElementSet::iterator jt = (*it)->begin(); jt != (*it)->end(); jt++)
				switch(jt->GetElementType())
				{
					case NODE: sets[i]->Insert(nodes[jt->LocalID()]); break;
					case EDGE: sets[i]->Insert(edges[jt->LocalID()]); break;
					case FACE: sets[i]->Insert(faces[jt->LocalID()]); break;
					case CELL: sets[i]->Insert(cells[jt->LocalID()]); break;
				}
		}
		epsilon = other.epsilon;
		AssignGlobalID(other.have_global_id);
		return *this;
	}
	
	Mesh::~Mesh()
	{
		//Mesh::iteratorTag t = BeginTag();
		//while(t != EndTag()) DeleteTag(*t);
		for(sets_container::iterator it = sets.begin(); it != sets.end(); it++) if( *it != NULL ) delete *it;
		for(cells_container::iterator it = cells.begin(); it != cells.end(); it++) if( *it != NULL ) 
		{
			(*it)->HighConn().clear();
			(*it)->LowConn().clear();
			delete *it;
		}
		for(faces_container::iterator it = faces.begin(); it != faces.end(); it++) if( *it != NULL ) 
		{
			(*it)->HighConn().clear();
			(*it)->LowConn().clear();
			delete *it;
		}
		for(edges_container::iterator it = edges.begin(); it != edges.end(); it++) if( *it != NULL ) 
		{
			(*it)->HighConn().clear();
			(*it)->LowConn().clear();
			delete *it;
		}
		for(nodes_container::iterator it = nodes.begin(); it != nodes.end(); it++) if( *it != NULL ) 
		{
			(*it)->HighConn().clear();
			(*it)->LowConn().clear();
			delete *it;
		}
#if defined(USE_MPI)
#if defined(USE_MPI2)
		if( m_state == Mesh::Parallel )
		{
			MPI_Free_mem(shared_space);
			MPI_Win_free(&window);
			//~ MPI_Comm_free(&comm);
		}
#endif
#endif
#if defined(USE_PARALLEL_WRITE_TIME)
		out_time << "</Debug>" << std::endl;
		out_time.close();
#endif
	}
	
	
	
	Tag Mesh::CreateTag(std::string name, DataType dtype, ElementType etype,ElementType sparse, INMOST_DATA_ENUM_TYPE size)
	{
		Tag ret = TagManager::CreateTag(this,name,dtype,etype,sparse,size);
		return ret;
	}
	Tag Mesh::DeleteTag(Tag tag, ElementType type_mask)
	{
		for(ElementType mask = NODE; mask <= MESH; mask = mask << 1 ) if( mask & type_mask )
		{
			if( tag.isDefined(mask) && (tag.isSparse(mask) || tag.GetSize() == ENUMUNDEF) )
				for(iteratorElement it = BeginElement(mask); it != EndElement(); it++)
					it->DelData(tag);
		}
		tag = TagManager::DeleteTag(tag,type_mask);
		return tag;
	}
	
	
	
	Element * Mesh::FindSharedAdjacency(Element * const * arr, unsigned s) 
	{
		if( s == 0 ) return NULL;
		if( !HideMarker() )
		{
			/*
			dynarray<Element *, 64> inter(arr[0]->HighConn().begin(),arr[0]->HighConn().end());
			MarkerType mrk = CreateMarker();
			for(unsigned i = 1; i < s; i++)
			{
				for(Element::adj_iterator jt = arr[i]->HighConn().begin(); jt != arr[i]->HighConn().end(); ++jt) (*jt)->SetMarker(mrk);
				{
					size_t m = 0, n = 0;
					while( m < inter.size() ) 
					{
						if( inter[m]->GetMarker(mrk) )
							inter[n++] = inter[m];
						m++;
					}
					inter.resize(n);
				}
				for(Element::adj_iterator jt = arr[i]->HighConn().begin(); jt != arr[i]->HighConn().end(); ++jt) (*jt)->RemMarker(mrk);
				if( inter.empty() )
				{
					ReleaseMarker(mrk);
					return NULL;
				}
			}	
			ReleaseMarker(mrk);
			return inter[0];
			*/
			/*
			if( arr[0]->GetElementType() == NODE )// this is edge
			{
				//there are 2 nodes
				if( s == 2 )
				{
					//check that any of edges of current node links to the other node
					Element::adj_type const & hc = arr[0]->HighConn();
					int i, iend = hc.size(), ss = s;
					for(i = 0; i < iend; ++i)
					{
						Element::adj_type const & lc = hc[i]->LowConn();
						//assert(lc.size() == 2);
						if( (lc[0] == arr[1] || lc[1] == arr[1]) ) return hc[i];
					}
				}
				else if( s == 1 )
				{
					Element::adj_type const & hc = arr[0]->HighConn();
					if( !hc.empty() ) return hc.front();
				}
				return NULL;
			}
			else
			*/
			{
				int flag0, flag1, i , ss = s;
				dynarray<Element::adj_type const *, 64> hcarr(ss);
				for(i = 0; i < ss; i++) hcarr[i] = &arr[i]->HighConn();
				int it, iend = hcarr[0]->size(), jt, jend;
				for(it = 0; it < iend; ++it)
				{
					flag0 = 0;
					for(i = 1; i < ss; ++i)
					{
						jend = hcarr[i]->size();
						flag1 = 0;
						for(jt = 0; jt < jend; ++jt)
						{
							if( hcarr[0]->at(it) == hcarr[i]->at(jt) )
							{
								flag0++;
								flag1 = 1;
								break;
							}
						}
						if( flag1 == 0 ) break;
					}
					if( flag0 == ss-1 ) return hcarr[0]->at(it);
				}
			}
			
			/*
			{
				//typedef std::map<Element *,int> set_type;
				//typedef tiny_map<Element *,int,64> set_type;
				//typedef small_hash<unsigned,int,64> set_type;
				typedef std::unordered_map<unsigned,int> set_type;
				set_type visits;
				for(unsigned i = 0; i < s; ++i)
				{
					unsigned jend = arr[i]->HighConn().size();
					for(unsigned jt = 0; jt < jend; ++jt)
						visits[reinterpret_cast<unsigned>(arr[i]->HighConn()[jt])]++;
				}
				for(set_type::iterator it = visits.begin(); it != visits.end(); ++it)
					if( it->second == s ) return reinterpret_cast<Element *>(it->first);
			}
			*/
			/*
			{
				//typedef std::map<Element *,int> set_type;
				//typedef tiny_map<Element *,int,64> set_type;
				//typedef small_hash<unsigned,int,64> set_type;
				Element * ret = NULL;
				unsigned q = s;
				for(unsigned i = 0; i < s && !ret; ++i)
				{
					unsigned jend = arr[i]->HighConn().size();
					for(unsigned jt = 0; jt < jend && !ret; ++jt)
					{
						if( ++arr[i]->HighConn()[jt]->RealDF(tag_shared_elems) == s)
						{
							ret = arr[i]->HighConn()[jt];
							q = i+1;
						}
					}
				}
				for(unsigned i = 0; i < q; ++i)
				{
					unsigned jend = arr[i]->HighConn().size();
					for(unsigned jt = 0; jt < jend; ++jt)
						arr[i]->HighConn()[jt]->RealDF(tag_shared_elems) = 0.0;
				}
				return ret;
			}
			*/
		}
		else
		{
			unsigned ss = Mesh::Count(arr,s,HideMarker());
			Element::adj_type & const hc0 = arr[0]->HighConn();
			unsigned iend = hc0.size();
			for(unsigned it = 0; it < iend; ++it) if( !hc0[it]->Hidden() )
			{
				unsigned int flag0 = 0;
				for(unsigned i = 1; i < s; i++)
				{
					bool flag1 = false;
					Element::adj_type & const hci = arr[i]->HighConn();
					unsigned jend = hci.size();
					for(unsigned jt = 0; jt < jend; ++jt) if( !hci[jt]->Hidden() )
					{
						if( hc0[it] == hci[jt] )
						{
							flag0++;
							flag1 = true;
							break;
						}
					}
					if( !flag1 )
						break;
				}
				if( flag0 == ss-1 ) return hc0[it];
			}
		}
		return NULL;
	}
	
	
	TopologyCheck Mesh::BeginTopologyCheck(ElementType etype, Element ** adj, INMOST_DATA_ENUM_TYPE s)
	{
		INMOST_DATA_ENUM_TYPE i,j,d = ENUMUNDEF;
		TopologyCheck chk = 0;
		if( GetTopologyCheck(ADJACENT_HIDDEN) )
		{
			for(i = 0; i < s; i++)
				if( adj[i]->GetMarker(HideMarker()) )
				{
					chk |= ADJACENT_HIDDEN;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_HIDDEN) << std::endl;
				}
		}
		if( GetTopologyCheck(ADJACENT_ALIEN) )
		{
			for(i = 0; i < s; i++)
				if( adj[i]->GetMeshLink() != this )
				{
					chk |= ADJACENT_ALIEN;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_ALIEN) << std::endl;
				}
		}
		if( GetTopologyCheck(ADJACENT_DUPLICATE) )
		{
			bool have_dup = false;
			MarkerType dup = CreateMarker();
			for(i = 0; i < s; i++)
			{
				if( adj[i]->GetMarker(dup) ) 
				{
					have_dup = true; // duplication of element
					break;
				}
				else adj[i]->SetMarker(dup);
			}
			for(i = 0; i < s; i++) adj[i]->RemMarker(dup);
			ReleaseMarker(dup);	
			
			if( have_dup )
			{
				chk |= ADJACENT_DUPLICATE;
				if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_DUPLICATE) << std::endl;
			}
		}
		
		if( GetTopologyCheck(ADJACENT_DIMENSION | DEGENERATE_CELL | DEGENERATE_FACE | DEGENERATE_EDGE | DISABLE_2D) )
		{
			bool happen = false;
			for(i = 0; i < s; i++)
			{
				j = adj[i]->GetElementDimension();
				if( j == ENUMUNDEF ) 
				{
					chk |= ADJACENT_DIMENSION;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_DIMENSION) << std::endl;
					happen = true;
					break;
				}
				if( d == ENUMUNDEF ) d = j;
				if( d != j ) 
				{
					chk |= ADJACENT_DIMENSION;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_DIMENSION) << std::endl;
					happen = true;
					break;
				}
			}
			if(!happen) switch(etype) //if happen = true, we cannot be sure about value in d
			{
				case EDGE: 
					if( GetTopologyCheck(DISABLE_2D) )
					{
						if( s == 1 && d == 0 ) 
						{
							chk |= DISABLE_2D;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DISABLE_2D) << std::endl;
						}
					}
					if( GetTopologyCheck(ADJACENT_DIMENSION) )
					{
						if( d > 0 ) 
						{
							chk |= ADJACENT_DIMENSION;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_DIMENSION) << std::endl;
						}
					}
					if( GetTopologyCheck(DEGENERATE_EDGE) )
					{
						if( s > 2 )
						{
							chk |= DEGENERATE_EDGE;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DEGENERATE_EDGE) << std::endl;
						}
					}
				break;
				case FACE:
					if( GetTopologyCheck(DISABLE_2D) )
					{
						if( s == 2 && d == 0 ) 
						{
							chk |= DISABLE_2D;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DISABLE_2D) << std::endl;
						}
					}
					if( GetTopologyCheck(ADJACENT_DIMENSION) )
					{
						if( d > 1 ) 
						{
							chk |= ADJACENT_DIMENSION;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_DIMENSION) << std::endl;
						}
					}
					if( GetTopologyCheck(DEGENERATE_EDGE) )
					{
						if(d == 0 && s != 2) //This should be Edge (2D)
						{
							chk |= DEGENERATE_EDGE;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DEGENERATE_EDGE) << std::endl;
						}
					}
					if( GetTopologyCheck() & DEGENERATE_FACE )
					{
						if(d == 1 && s < 3) // Should not be less then 3 Line or Curve
						{
							chk |= DEGENERATE_FACE;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DEGENERATE_FACE) << std::endl;
						}					
					}
				break;
				case CELL:
					if( GetTopologyCheck(DISABLE_2D) )
					{
						if( d == 1 ) 
						{
							chk |= DISABLE_2D;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DISABLE_2D) << std::endl;
						}
					}
					if( GetTopologyCheck(ADJACENT_DIMENSION) )
					{
						if( d > 2 || d == 0 ) 
						{
							chk |= ADJACENT_DIMENSION;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_DIMENSION) << std::endl;
						}
					}
					if( GetTopologyCheck(DEGENERATE_FACE) )
					{
						if( d == 1 && s < 3 ) //Should not be less then 3 Line
						{
							chk |= DEGENERATE_FACE;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DEGENERATE_FACE) << std::endl;
						} 
					}
					if( GetTopologyCheck(DEGENERATE_CELL) )
					{
						if( d == 2 && s < 4 ) //Should not be less then 4 Tri, Quad, Polygon
						{
							chk |= DEGENERATE_CELL;
							if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(DEGENERATE_CELL) << std::endl;
						} 
					}
				break;
			}
		}
		if( etype == CELL )
		{
			if( GetTopologyCheck(TRIPLE_SHARED_FACE) )
			{
				bool check = false;
				for(i = 0; i < s; i++)
				{
					if( adj[i]->nbAdjElements(CELL) == 2 )
					{
						check = true;
						break;
					}
				}
				if( check )
				{
					bool happen = true;
					if( (GetTopologyCheck(DUPLICATE_CELL)) && (FindSharedAdjacency(adj,s) != NULL) ) happen = false;
					if( happen )
					{
						chk |= TRIPLE_SHARED_FACE;
						if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(TRIPLE_SHARED_FACE) << std::endl;
					}
				}
			}
		}
		else if( etype == FACE )
		{
			if( GetTopologyCheck(INTERLEAVED_FACES) )
			{
				adjacent<Node> allnodes;
				for(i = 0; i < s; i++) allnodes.unite(adj[i]->getNodes());
				adjacent<Face> faces = allnodes[0].getFaces();
				for(i = 1; i < allnodes.size() && !faces.empty(); i++) faces.intersect(allnodes[i].getFaces());
				for(i = 0; i < faces.size(); i++)
				{
					if( faces[i].nbAdjElements(NODE) != allnodes.size() )
					{
						chk |= INTERLEAVED_FACES;
						if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(INTERLEAVED_FACES) << std::endl;
					}
				}
			}
		}
		errorset |= chk;
		return chk;
	}
	
	TopologyCheck Mesh::EndTopologyCheck(Element * e)
	{
		TopologyCheck chk = 0;
		switch(e->GetGeometricDimension(e->GetGeometricType()))
		{
		case 3:
			if( GetTopologyCheck(FLATTENED_CELL) )
			{
				adjacent<Face> faces = e->getFaces();
				INMOST_DATA_ENUM_TYPE num = e->nbAdjElements(NODE);
				for(INMOST_DATA_ENUM_TYPE i = 0; i < faces.size(); i++)
				{
					if( faces[i].nbAdjElements(NODE) == num )
					{
						chk |= FLATTENED_CELL;
						if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(FLATTENED_CELL) << std::endl;
					}
				}
			}
			if( GetTopologyCheck(FACE_ORIENTATION) )
			{
				adjacent<Face> faces = e->getFaces();
				for(INMOST_DATA_ENUM_TYPE i = 0; i < faces.size(); i++)
				{
					if( faces[i].nbAdjElements(CELL) == 1 && !faces[i].CheckNormalOrientation() )
					{
						chk |= FACE_ORIENTATION;
						if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(FACE_ORIENTATION) << std::endl;
						break;
					}
				}
			}
			break;
		case 2:
			if( GetTopologyCheck(FACE_PLANARITY) )
			{
				if( (e->GetElementType() == FACE && !e->getAsFace()->Planarity()) || (e->GetElementType() == CELL && !e->getAsCell()->Planarity()) )
				{
					chk |= FACE_PLANARITY;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(FACE_PLANARITY) << std::endl;
				}
			}
			if( GetTopologyCheck(FACE_EDGES_ORDER) )
			{
				if( (e->GetElementType() == FACE && !e->getAsFace()->CheckEdgeOrder()) || (e->GetElementType() == CELL && !e->getAsCell()->CheckEdgeOrder()) )
				{
					chk |= FACE_EDGES_ORDER;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(FACE_EDGES_ORDER) << std::endl;
				}
			}
			break;
		}
		if( GetTopologyCheck(PROHIBIT_MULTILINE) )
		{
			if( e->GetGeometricType() == Element::MultiLine )
			{
				chk |= PROHIBIT_MULTILINE;
				if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(PROHIBIT_MULTILINE) << std::endl;
			}
		}
		if( GetTopologyCheck(PROHIBIT_MULTIPOLYGON) )
		{
			if( e->GetGeometricType() == Element::MultiPolygon )
			{
				chk |= PROHIBIT_MULTIPOLYGON;
				if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(PROHIBIT_MULTIPOLYGON) << std::endl;
			}
		}
		if( GetTopologyCheck(PROHIBIT_POLYGON) )
		{
			if( e->GetGeometricType() == Element::Polygon )
			{
				chk |= PROHIBIT_POLYGON;
				if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(PROHIBIT_POLYGON) << std::endl;
			}
		}
		if( GetTopologyCheck(PROHIBIT_POLYHEDRON) )
		{
			if( e->GetGeometricType() == Element::Polyhedron )
			{
				chk |= PROHIBIT_POLYHEDRON;
				if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(PROHIBIT_POLYHEDRON) << std::endl;
			}
		}
		errorset |= chk;
		return chk;
	}
	
	Node * Mesh::CreateNode(Storage::real * coords)
	{
		INMOST_DATA_ENUM_TYPE i = 0;
		Node * e = new Node(this);
		e->SetGeometricType(Element::Vertex);
		Storage::real_array v = e->RealArray(tag_coords);
		for(i = 0; i < dim; i++) v[i] = coords[i];
		e->SetMarker(NewMarker());
		return e;
	}
	std::pair<Edge *,bool> Mesh::CreateEdge(Node ** e_nodes, INMOST_DATA_ENUM_TYPE s)
	{
		INMOST_DATA_ENUM_TYPE i;
		Edge * e = NULL;
		if( s > 0 )
		{
			TopologyCheck chk = 0;
			chk |= BeginTopologyCheck(EDGE,reinterpret_cast<Element **>(&e_nodes[0]),s);
			if( chk != 0 )
			{
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
			if( GetTopologyCheck() & DUPLICATE_EDGE )
			{
				Element * test = FindSharedAdjacency(reinterpret_cast<Element **>(&e_nodes[0]),s);
				if (test != NULL) return std::pair<Edge *, bool>(static_cast<Edge *>(test),false);
			}		
			e = new Edge(this);
			for(i = 0; i < s; i++)
			{
				e_nodes[i]->HighConn().push_back(e);
				//e->LowConn().push_back(e_nodes[i]);
			}
			Element::adj_type & lc = e->LowConn();
			lc.insert(lc.end(),reinterpret_cast<Element **>(e_nodes),reinterpret_cast<Element **>(e_nodes+s));
			e->ComputeGeometricType();
			e->SetMarker(NewMarker());
			RecomputeGeometricData(e);
			chk |= EndTopologyCheck(e);
			if( chk != 0 )
			{
				if( GetTopologyCheck(MARK_ON_ERROR) ) e->Integer(TopologyErrorTag()) = chk;
				if( GetTopologyCheck(DELETE_ON_ERROR) ) { delete e; e = NULL;}
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
		}
		return std::pair<Edge *, bool>(e,true);
	}
	
	
	std::pair<Face *,bool> Mesh::CreateFace(Node ** f_nodes, INMOST_DATA_ENUM_TYPE s)
	{
		dynarray<Edge *,64> f_edges(s);
		Node * e_nodes[2];
		if( s == 2 ) //This is an edge for 2d!
		{
			e_nodes[0] = f_nodes[0];
			f_edges[0] = CreateEdge(e_nodes,1).first;
			e_nodes[0] = f_nodes[1];
			f_edges[1] = CreateEdge(e_nodes,1).first;
		}
		else
		{
			for(unsigned int i = 0; i < s; i++)
			{
				e_nodes[0] = f_nodes[i];
				e_nodes[1] = f_nodes[(i+1)%s];
				f_edges[i] = CreateEdge(e_nodes,2).first;
			}
		}
		return CreateFace(&f_edges[0],s);
	}
	
	std::pair<Face *,bool> Mesh::CreateFace(Edge ** f_edges, INMOST_DATA_ENUM_TYPE s)
	{
		INMOST_DATA_ENUM_TYPE i;
		Face * e = NULL;
		if( s > 0 )
		{
			TopologyCheck chk = 0;
			chk |= BeginTopologyCheck(FACE,reinterpret_cast<Element **>(&f_edges[0]),s);
			if( chk != 0 )
			{
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
			if( GetTopologyCheck(DUPLICATE_FACE) )
			{
				Element * test = FindSharedAdjacency(reinterpret_cast<Element **>(&f_edges[0]),s);
				if (test != NULL) return std::pair<Face *, bool>(static_cast<Face *>(test),false);
			}
			e = new Face(this);
			for(i = 0; i < s; i++)
			{
				f_edges[i]->HighConn().push_back(e);
				//e->LowConn().push_back(f_edges[i]);
			}
			Element::adj_type & lc = e->LowConn();
			lc.insert(lc.end(),reinterpret_cast<Element **>(f_edges),reinterpret_cast<Element **>(f_edges+s));
			e->ComputeGeometricType();
			e->SetMarker(NewMarker());
			RecomputeGeometricData(e);
			chk |= EndTopologyCheck(e);
			if( chk != 0 )
			{
				if( GetTopologyCheck(MARK_ON_ERROR)   ) e->Integer(TopologyErrorTag()) = chk;
				if( GetTopologyCheck(DELETE_ON_ERROR) ) { delete e; e = NULL;}
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
		}
		return std::pair<Face *, bool>(e,true);
	}
	
	
	std::pair<Cell *,bool> Mesh::CreateCell(Node ** c_f_nodes, const INMOST_DATA_ENUM_TYPE * c_f_sizes, INMOST_DATA_ENUM_TYPE s, Node ** suggest_nodes_order, INMOST_DATA_ENUM_TYPE sn)
	{
		INMOST_DATA_ENUM_TYPE i,j = 0;
		dynarray<Face *,64> c_faces;
		c_faces.resize(s);
		for(i = 0; i < s; i++)
		{
			c_faces[i] = CreateFace(&c_f_nodes[j], c_f_sizes[i]).first;
			j += c_f_sizes[i];
		}
		return CreateCell(&c_faces[0],s,suggest_nodes_order,sn);
	}
	

	std::pair<Cell *, bool> Mesh::CreateCell(Node ** c_f_nodes, const INMOST_DATA_ENUM_TYPE * c_f_nodeinds, const INMOST_DATA_ENUM_TYPE * c_f_numnodes, INMOST_DATA_ENUM_TYPE s, Node ** suggest_nodes_order, INMOST_DATA_ENUM_TYPE sn)
	{
		INMOST_DATA_ENUM_TYPE i,k,j = 0; //, sn = 0;
		dynarray<Node *,64> temp;
		dynarray<Face *,64> c_faces;
		c_faces.resize(s);
		for(i = 0; i < s; i++)
		{
			for(k = j; k < j+c_f_numnodes[i]; k++)
			{
				//if( c_f_nodeinds[k]+1 > sn ) sn = c_f_nodeinds[k]+1;
				temp.push_back(c_f_nodes[c_f_nodeinds[k]]);
			}
			c_faces[i] = CreateFace(&temp[0], c_f_numnodes[i]).first;
			j += c_f_numnodes[i];
			temp.clear();
		}
		return CreateCell(&c_faces[0],s,suggest_nodes_order,sn);
	}
	
	
	void Mesh::RestoreCellNodes(Cell * c, dynarray<Node *,64> & ret)
	{
		ret.clear();
		switch(c->GetGeometricType())
		{
			case Element::Vertex:
			case Element::Line:
			{
				throw CorruptedIerarchy;
			}
			case Element::MultiLine:
			case Element::Tri:
			case Element::Quad:
			case Element::Polygon:
			{
				adjacent<Edge> edges = c->getEdges();
				ret.reserve(edges.size());
				for(adjacent<Edge>::iterator it = edges.begin(); it != edges.end(); it++)
				{
					adjacent<Node> nodes = it->getNodes();
					if( nodes.size() != 1 ) throw CorruptedIerarchy;
					ret.push_back(nodes.data()[0]);
				}
				break;
			}
				/*
				 6 7
				 4 5
				 2 3
				 0 1
				 */
			case Element::Hex:
			{
				MarkerType mrk = CreateMarker();
				MarkerType cemrk = CreateMarker();
				MarkerType femrk = CreateMarker();
				//printf("%lx %lx %lx\n",mrk,cemrk,femrk);
				adjacent<Face> faces = c->getFaces();
				Face * face = &faces[0];
				ret.reserve(8);
				adjacent<Node> verts = face->getNodes();
				if( face->BackCell() == c )
					for(adjacent<Node>::iterator it = verts.begin(); it != verts.end(); it++)
					{
						ret.push_back(&*it);
						it->SetMarker(mrk);
					}
				else
					for(adjacent<Node>::reverse_iterator it = verts.rbegin(); it != verts.rend(); it++)
					{
						ret.push_back(&*it);
						it->SetMarker(mrk);
					}
				adjacent<Edge> c_edges = c->getEdges();
				for(adjacent<Edge>::iterator it = c_edges.begin(); it != c_edges.end(); it++)
					it->SetMarker(cemrk);
				adjacent<Edge> f_edges = face->getEdges();
				for(adjacent<Edge>::iterator it = f_edges.begin(); it != f_edges.end(); it++)
					it->SetMarker(femrk);
				for(unsigned int k = 0; k < 4; k++)
				{
					adjacent<Edge> v_edges = ret[k]->getEdges();
					for(adjacent<Edge>::iterator it = v_edges.begin(); it != v_edges.end(); it++)
					{
						if( it->GetMarker(cemrk) && !it->GetMarker(femrk) )
						{
							adjacent<Node> nn = it->getNodes();
							if( nn[0].GetMarker(mrk) )
								ret.push_back(&nn[1]);
							else
								ret.push_back(&nn[0]);
							break;
						}
					}
					ret[k]->RemMarker(mrk);
				}
				for(adjacent<Edge>::iterator it = c_edges.begin(); it != c_edges.end(); it++) it->RemMarker(cemrk);
				for(adjacent<Edge>::iterator it = f_edges.begin(); it != f_edges.end(); it++) it->RemMarker(femrk);
				ReleaseMarker(mrk);
				ReleaseMarker(cemrk);
				ReleaseMarker(femrk);
				break;
			}
				/*
				 5
				 3 4
				 2 
				 0 1
				 */
			case Element::Prism:
			{
				MarkerType mrk = CreateMarker();
				MarkerType cemrk = CreateMarker();
				MarkerType femrk = CreateMarker();
				ret.reserve(6);
				Face * face;
				adjacent<Face> faces = c->getFaces();
				for(unsigned int i = 0 ; i < faces.size(); i++) //iterate over faces
					if( faces[i].nbAdjElements(EDGE) == 3 ) //number of edges in i-th face
					{
						face = &faces[i];
						break;
					}
				adjacent<Node> verts = face->getNodes();
				if( face->BackCell() == c )
					for(adjacent<Node>::iterator it = verts.begin(); it != verts.end(); it++)
					{
						ret.push_back(&*it);
						it->SetMarker(mrk);
					}
				else
					for(adjacent<Node>::reverse_iterator it = verts.rbegin(); it != verts.rend(); it++)
					{
						ret.push_back(&*it);
						it->SetMarker(mrk);
					}
				adjacent<Edge> c_edges = c->getEdges();
				for(adjacent<Edge>::iterator it = c_edges.begin(); it != c_edges.end(); it++) it->SetMarker(cemrk);
				adjacent<Edge> f_edges = face->getEdges();
				for(adjacent<Edge>::iterator it = f_edges.begin(); it != f_edges.end(); it++) it->SetMarker(femrk);
				for(unsigned int k = 0; k < 3; k++)
				{
					adjacent<Edge> v_edges = ret[k]->getEdges();
					for(adjacent<Edge>::iterator it = v_edges.begin(); it != v_edges.end(); it++)
						if( it->GetMarker(cemrk) && !it->GetMarker(femrk) )
						{
							adjacent<Node> nn = it->getNodes();
							if( nn[0].GetMarker(mrk) )
								ret.push_back(&nn[1]);
							else
								ret.push_back(&nn[0]);
							break;
						}
					ret[k]->RemMarker(mrk);
				}
				for(adjacent<Edge>::iterator it = c_edges.begin(); it != c_edges.end(); it++) it->RemMarker(cemrk);
				for(adjacent<Edge>::iterator it = f_edges.begin(); it != f_edges.end(); it++) it->RemMarker(femrk);
				ReleaseMarker(mrk);
				ReleaseMarker(cemrk);
				ReleaseMarker(femrk);
				break;
			}
				/*
				  4
				 2 3
				 0 1
				 */
			case Element::Pyramid:
			{
				ret.reserve(5);
				Face * quad, * triangle;
				MarkerType mrk = CreateMarker();
				adjacent<Face> faces = c->getFaces();
				for(unsigned int i = 0; i < faces.size(); i++) //go over faces
					if( faces[i].nbAdjElements(EDGE) == 4 ) //check if number of edges = 4
					{
						quad = &faces[i];
						break;
					}
				for(unsigned int i = 0; i < faces.size(); i++) //go over faces
					if( faces[i].nbAdjElements(EDGE) == 3 ) //check if number of edges = 3
					{
						triangle = &faces[i];
						break;
					}
				adjacent<Node> base_nodes = quad->getNodes();
				if( quad->BackCell() == c )
					for(adjacent<Node>::iterator it = base_nodes.begin(); it != base_nodes.end(); it++)
					{
						ret.push_back(&*it);
						it->SetMarker(mrk);
					}
				else
					for(adjacent<Node>::reverse_iterator it = base_nodes.rbegin(); it != base_nodes.rend(); it++)
					{
						ret.push_back(&*it);
						it->SetMarker(mrk);
					}
				adjacent<Node> tri_nodes = triangle->getNodes();
				for(adjacent<Node>::iterator it = tri_nodes.begin(); it != tri_nodes.end(); it++)
				{
					if( !it->GetMarker(mrk) )
					{
						ret.push_back(&*it);
						break;
					}
				}
				for(dynarray<Node *,64>::iterator it = ret.begin(); it != ret.end(); it++)
					(*it)->RemMarker(mrk);
				ReleaseMarker(mrk);
				break;
			}
			default: //Tet, MultiPolygon, Polyhedron
			{
				MarkerType mrk = CreateMarker();
				if( !HideMarker() )
				{
					Element::adj_type & lc = c->LowConn();
					for(Element::adj_type::enumerator i = 0; i < lc.size(); i++) //iterate over faces
					{
						Element::adj_type & ilc = lc[i]->LowConn();
						for(Element::adj_type::enumerator j = 0; j < ilc.size(); j++) //iterate over face edges
						{
							Element::adj_type & jlc = ilc[j]->LowConn();
							for(Element::adj_type::enumerator k = 0; k < jlc.size(); k++) //iterator over edge nodes
							{
								if( !jlc[k]->GetMarker(mrk) )
								{
									jlc[k]->SetMarker(mrk);
									ret.push_back(static_cast<Node *>(jlc[k]));
								}
							}
						}
					}
				}
				else
				{
					Element::adj_type & lc = c->LowConn();
					for(Element::adj_type::enumerator i = 0; i < lc.size(); i++) if( !lc[i]->Hidden() )  //iterate over faces
					{
						Element::adj_type & ilc = lc[i]->LowConn();
						for(Element::adj_type::enumerator j = 0; j < ilc.size(); j++) if( !ilc[j]->Hidden() ) //iterate over face edges
						{
							Element::adj_type & jlc = ilc[j]->LowConn();
							for(Element::adj_type::enumerator k = 0; k < jlc.size(); k++) if( !jlc[k]->Hidden() ) //iterator over edge nodes
							{
								if( !jlc[k]->GetMarker(mrk) )
								{
									jlc[k]->SetMarker(mrk);
									ret.push_back(static_cast<Node *>(jlc[k]));
								}
							}
						}
					}
				}
				for(dynarray<Node *,64>::enumerator it = 0; it < ret.size(); it++) ret[it]->RemMarker(mrk);
				ReleaseMarker(mrk);
				break;
			}
		}

	}
	
	std::pair<Cell *,bool> Mesh::CreateCell(Face ** c_faces, INMOST_DATA_ENUM_TYPE s, Node ** c_nodes, INMOST_DATA_ENUM_TYPE sn)
	{
		INMOST_DATA_ENUM_TYPE i;
		Cell * e = NULL;
		if( s > 0 )
		{
			TopologyCheck chk = 0;
			chk |= BeginTopologyCheck(CELL,reinterpret_cast<Element **>(&c_faces[0]),s);
			if( chk != 0 )
			{
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
			if( GetTopologyCheck(DUPLICATE_CELL) )
			{
				Element * test = FindSharedAdjacency(reinterpret_cast<Element **>(&c_faces[0]),s);
				if (test != NULL) return std::pair<Cell *, bool>(static_cast<Cell *>(test),false);
			}
			e = new Cell(this);		
			for(i = 0; i < s; i++)
			{
				c_faces[i]->HighConn().push_back(e);
				//e->LowConn().push_back(c_faces[i]);
			}
			Element::adj_type & lc = e->LowConn();
			lc.insert(lc.begin(),reinterpret_cast<Element **>(c_faces),reinterpret_cast<Element **>(c_faces+s));
			e->ComputeGeometricType();		
			{
				dynarray<Node *,64> temp_nodes;
				if( sn == 0 ) 
				{
					RestoreCellNodes(e,temp_nodes);
					c_nodes = &temp_nodes[0];
					sn = temp_nodes.size();
				}
				for(INMOST_DATA_ENUM_TYPE k = 0; k < sn; k++)
				{
					c_nodes[k]->LowConn().push_back(static_cast<Element *>(e));
					//e->HighConn().push_back(static_cast<Element *>(c_nodes[k]));
				}
				Element::adj_type & hc = e->HighConn();
				hc.insert(hc.begin(),reinterpret_cast<Element **>(c_nodes),reinterpret_cast<Element **>(c_nodes+sn));
				e->SetMarker(NewMarker());
				RecomputeGeometricData(e);
				chk |= EndTopologyCheck(e);
				if( chk != 0 )
				{
					if( GetTopologyCheck(MARK_ON_ERROR)   ) e->Integer(TopologyErrorTag()) = chk;
					if( GetTopologyCheck(DELETE_ON_ERROR) ) { delete e; e = NULL;}
					if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
				}
			}
		}
		return std::pair<Cell *, bool>(e,true);
	}
	
	
	ElementSet * Mesh::CreateSet()
	{
		ElementSet * e = new ElementSet(this,false);
		return e;
	}
	
	ElementSet * Mesh::CreateOrderedSet()
	{
		ElementSet * e = new ElementSet(this,true);
		return e;
	}
	
	
	void Mesh::SetDimensions(unsigned int dims)
	{
		if( dim == dims ) return;
		if( dims < 2 || dims > 3 ) throw DimensionIsNotSupportedByGeometry;

		if( NumberOfNodes() > 0)
		{
			Storage::real * temp = new Storage::real[dims*NumberOfNodes()];
			Storage::integer j = 0;
			for(Storage::integer k = 0; k < MaxLocalIDNODE(); ++k)
			{
				Node * n = NodeByLocalID(k);
				if( n != NULL )
				{
					Storage::real_array c = n->Coords();
					for(Storage::integer i = 0; i < dims; i++)
					{
						temp[j+i] = c[i];
					}
					j+=dims;
				}
			}
			DeleteTag(tag_coords);
			tag_coords = CreateTag("COORD",DATA_REAL,NODE,NONE,dims);
			j = 0;
			for(Storage::integer k = 0; k < MaxLocalIDNODE(); ++k)
			{
				Node * n = NodeByLocalID(k);
				if( n != NULL )
				{
					Storage::real_array c = n->Coords();
					for(Storage::integer i = 0; i < dims; i++)
					{
						c[i] = temp[j+i];
					}
					j+=dims;
				}
			}
			dim = dims;
		}
		else
		{
			dim = dims;
			DeleteTag(tag_coords);
			tag_coords = CreateTag("COORD",DATA_REAL, NODE,NONE,dim);
		}
	}
	
	
	bool LessElements(Element * a, Element * b)
	{
		return USE_COMPARE(a,b) < 0;
	}
	int CompareElementsC(const void * pa, const void * pb)
	{
		Element ** a = (Element **)pa;
		Element ** b = (Element **)pb;
		return USE_COMPARE(*a,*b);
	}
	int CompareElements(Element *a, Element * b)
	{
		return USE_COMPARE(a,b);
	}
	
	int CompareElementsCCentroid(const void * pa, const void * pb)
	{
		Element ** a = (Element **)pa;
		Element ** b = (Element **)pb;
		return CompareElementsCentroid(*a,*b);
	}
	
	int CompareElementsCGID(const void * pa, const void * pb)
	{
		Element ** a = (Element **)pa;
		Element ** b = (Element **)pb;
		return CompareElementsGID(*a,*b);
	}

	
	int CompareElementsCUnique(const void * pa, const void * pb)
	{
		Element ** a = (Element **)pa;
		Element ** b = (Element **)pb;
		return CompareElementsUnique(*a,*b);
	}
	
	int CompareElementsCPointer(const void * pa, const void * pb)
	{
		Element ** a = (Element **)pa;
		Element ** b = (Element **)pb;
		return CompareElementsPointer(*a,*b);
	}
	
	bool LessElementsGID(Element * a, Element * b)
	{
		return CompareElementsGID(a,b) < 0;
	}
	
	bool LessElementsCentroid(Element * a, Element * b)
	{
		return CompareElementsCentroid(a,b) < 0;
	}
	
	bool LessElementsUnique(Element * a, Element * b)
	{
		return CompareElementsUnique(a,b) < 0;
	}

	bool LessElementsPointer(Element * a, Element * b)
	{
		return CompareElementsPointer(a,b) < 0;
	}

	
	int CompareElementsUnique(Element * a,Element * b)
	{
		Mesh * ma = a->GetMeshLink(), * mb = b->GetMeshLink();
		assert(ma == mb);
		if( a == b ) return 0;
		if( a == NULL ) return 1;
		if( b == NULL ) return -1;
		if( a->GetGeometricType() == b->GetGeometricType() )
		{
			if( a->GetElementType() != b->GetElementType() )
			{
				if( a->GetElementType() < b->GetElementType() )
					return CompareElementsUnique(a,b->LowConn()[0]);
				else
					return CompareElementsUnique(a->LowConn()[0],b);
			}
			switch( a->GetGeometricType() )
			{
				case Element::Vertex:
				{
					if( a->GetElementType() > NODE )
						return CompareElementsUnique(a->LowConn()[0],b->LowConn()[0]);
					Storage::real_array 
						ca = a->RealArray(ma->CoordsTag()),
						cb = b->RealArray(mb->CoordsTag());
					Storage::real e = std::max(ma->GetEpsilon(),mb->GetEpsilon());
					if( ca.size() != cb.size() ) return ca.size()-cb.size();
					for(unsigned int i = 0; i < ca.size(); i++)
						if( fabs(ca[i]-cb[i]) > e )
						{
							if( ca[i] > cb[i] ) return 1;
							else return -1;
						}
					return 0;
				}
				case Element::Line:
				{
					if( a->GetElementType() > EDGE )
						return CompareElementsUnique(a->LowConn()[0],b->LowConn()[0]);
					Element::adj_type ca, cb;
					if( CompareElementsUnique(a->LowConn().front(),a->LowConn().back()) <= 0 )
						ca.insert(ca.begin(),a->LowConn().begin(),a->LowConn().end());
					else
						ca.insert(ca.begin(),a->LowConn().rbegin(),a->LowConn().rend());
					if( CompareElementsUnique(b->LowConn().front(),b->LowConn().back()) <= 0 )
						cb.insert(cb.begin(),b->LowConn().begin(),b->LowConn().end());
					else
						cb.insert(cb.begin(),b->LowConn().rbegin(),b->LowConn().rend());
					for(unsigned int i = 0; i < ca.size(); i++)
					{
						int test = CompareElementsUnique(ca[i],cb[i]);
						if( test != 0 ) return test;
					}
					/*
					 for(int i = 0; i < ca.size(); i++)
					 if( ca[i] != cb[i] )
					 {
					 if( ca[i] > cb[i] ) return 1;
					 else return -1;
					 }
					 */
					return 0;
				}
				case Element::MultiLine:
				case Element::Tri:
				case Element::Quad:
				case Element::Polygon:
				{
					if( a->GetElementType() > FACE )
						return CompareElementsUnique(a->LowConn()[0],b->LowConn()[0]);
					if( a->LowConn().size() != b->LowConn().size() )
						return a->LowConn().size() - b->LowConn().size();
					Element::adj_type ca(a->LowConn()), cb(b->LowConn());
#if defined(USE_QSORT)
					qsort(&ca[0],ca.size(),sizeof(Element *),CompareElementsCUnique);
					qsort(&cb[0],cb.size(),sizeof(Element *),CompareElementsCUnique);
#else
					std::sort(ca.begin(),ca.end(),LessElementsUnique);
					std::sort(cb.begin(),cb.end(),LessElementsUnique);
#endif
					for(unsigned int i = 0; i < ca.size(); i++)
					{
						int test = CompareElementsUnique(ca[i],cb[i]);
						if( test != 0 ) return test;
					}
					return 0;
				}
				case Element::MultiPolygon:
				case Element::Tet:
				case Element::Hex:
				case Element::Prism:
				case Element::Pyramid:
				case Element::Polyhedron:
				{
					if( a->LowConn().size() != b->LowConn().size() )
						return a->LowConn().size() - b->LowConn().size();
					Element::adj_type ca(a->LowConn()), cb(b->LowConn());
#if defined(USE_QSORT)
					qsort(&ca[0],ca.size(),sizeof(Element *),CompareElementsCUnique);
					qsort(&cb[0],cb.size(),sizeof(Element *),CompareElementsCUnique);
#else
					std::sort(ca.begin(),ca.end(),LessElementsUnique);
					std::sort(cb.begin(),cb.end(),LessElementsUnique);
#endif
					for(unsigned int i = 0; i < ca.size(); i++)
					{
						int test = CompareElementsUnique(ca[i],cb[i]);
						if( test != 0 ) return test;
					}
					return 0;
				}
				default: throw Impossible;
			}
			throw Impossible;
		}
		else return a->GetGeometricType() - b->GetGeometricType();
	}
	
	int CompareElementsCentroid(Element * a,Element * b)
	{
		Mesh * ma = a->GetMeshLink(), * mb = b->GetMeshLink();
		assert(ma == mb);
		if( a == b ) return 0;
		if( a == NULL ) return 1;
		if( b == NULL ) return -1;
		if( a->GetGeometricType() == b->GetGeometricType() )
		{
			if( a->GetElementType() != b->GetElementType() )
			{
				if( a->GetElementType() < b->GetElementType() )
					return CompareElementsCentroid(a,b->LowConn()[0]);
				else
					return CompareElementsCentroid(a->LowConn()[0],b);
			}
			Storage::real ca[3] = {0,0,0};
			Storage::real cb[3] = {0,0,0};
			a->Centroid(ca);
			b->Centroid(cb);
			Storage::real e = std::max(ma->GetEpsilon(),mb->GetEpsilon());
			for(unsigned int i = 0; i < 3; i++)
				if( fabs(ca[i]-cb[i]) > e )
				{
					if( ca[i] > cb[i] ) return 1;
					else return -1;
				}
			return 0;
		}
		else return a->GetGeometricType() - b->GetGeometricType();
	}
	
	int CompareElementsGID(Element * a,Element * b)
	{
		if( a == b ) return 0;
		if( a == NULL ) return 1;
		if( b == NULL ) return -1;
		if( a->GetGeometricType() == b->GetGeometricType() )
			return a->GlobalID()-b->GlobalID();
		else return a->GetGeometricType() - b->GetGeometricType();
	}
	
	int CompareElementsPointer(Element * a, Element * b)
	{
		//printf("%ld\n",a-b);
		if( a > b ) return 1;
		if( a < b ) return -1;
		return 0;
	}
	
	int CompareElementsHybrid(Element * a,Element * b)
	{
		int test = CompareElementsCentroid(a,b);
		if( test == 0 ) test = CompareElementsUnique(a,b);
		return test;
	}
	
	
	void Mesh::UntieElement(Storage * e)
	{
#if defined(USE_OMP)
#pragma omp critical [storage_interraction]
#endif
		{
			INMOST_DATA_ENUM_TYPE pos = e->LocalID();
			switch(e->GetElementType())
			{
				case ESET:
					if( pos >= sets.size() || sets[pos] != e ) return;
					empty_sets.push_back(pos);
					sets[pos] = NULL;
					break;
				case CELL:
					if( pos >= cells.size() || cells[pos] != e ) return;
					empty_cells.push_back(pos);
					cells[pos] = NULL;
					break;
				case FACE:
					if( pos >= faces.size() || faces[pos] != e ) return;
					empty_faces.push_back(pos);
					faces[pos] = NULL;
					break;
				case EDGE:
					if( pos >= edges.size() || edges[pos] != e ) return;
					empty_edges.push_back(pos);
					edges[pos] = NULL;
					break;
				case NODE:
					if( pos >= nodes.size() || nodes[pos] != e ) return;
					empty_nodes.push_back(pos);
					nodes[pos] = NULL;
					break;
				default:
					throw Impossible;
			}
		}
		// This should be done through BeginModification - ApplyModification - EndModification
		//if( e->GetElementType() & (CELL | FACE | EDGE | NODE) )
		//	for(unsigned int i = 0; i < sets.size(); i++) if( sets[i] != NULL ) sets[i]->Erase(static_cast<Element *>(e));
	}
	
	
	void Mesh::TieElement(Storage * e)
	{
#if defined(USE_OMP)
#pragma omp critical [storage_interraction]
#endif
		{
			ElementType etype = e->GetElementType();
			INMOST_DATA_ENUM_TYPE old_size = GetArrayCapacity(etype);
			switch(e->GetElementType())
			{
				case ESET:
					if( !empty_sets.empty() && !isMeshModified() )
					{
						e->local_id = empty_sets.back();
						sets[empty_sets.back()] = static_cast<ElementSet *>(e);
						empty_sets.pop_back();
					}
					else
					{
						e->local_id = sets.size();
						sets.push_back(static_cast<ElementSet *>(e));
					}
					break;
				case NODE:
					if( !empty_nodes.empty() && !isMeshModified() )
					{
						e->local_id = empty_nodes.back();
						nodes[empty_nodes.back()] = static_cast<Node *>(e);
						empty_nodes.pop_back();
					}
					else 
					{
						e->local_id = nodes.size();
						nodes.push_back(static_cast<Node *>(e));
					}
					break;
				case EDGE:
					if( !empty_edges.empty() && !isMeshModified() )
					{
						e->local_id = empty_edges.back();
						edges[empty_edges.back()] = static_cast<Edge *>(e);
						empty_edges.pop_back();
					}
					else 
					{
						e->local_id = edges.size();
						edges.push_back(static_cast<Edge *>(e));
					}
					break;
				case FACE:
					if( !empty_faces.empty() && !isMeshModified() )
					{
						e->local_id = empty_faces.back();
						faces[empty_faces.back()] = static_cast<Face *>(e);
						empty_faces.pop_back();
					}
					else 
					{
						e->local_id = faces.size();
						faces.push_back(static_cast<Face *>(e));
					}
					break;
				case CELL:
					if( !empty_cells.empty() && !isMeshModified() )
					{
						e->local_id = empty_cells.back();
						cells[empty_cells.back()] = static_cast<Cell *>(e);
						empty_cells.pop_back();
					}
					else 
					{
						e->local_id = cells.size();
						cells.push_back(static_cast<Cell *>(e));
					}
					break;
				case MESH:
					{
						e->local_id = 0;
					}
					break;
			}
			if( GetArrayCapacity(etype) != old_size )
			{
				ReallocateData(etype);
				/*
				for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); ++t)
					if( t->isDefined(etype) && !t->isSparse(etype) )
						t->AllocateData(etype);
				*/
			}
		}
	}
	
	
	void Mesh::MoveStorage(Storage * e, int new_local_id)
	{
		e->MoveData(new_local_id);
		switch(e->GetElementType())
		{
			case CELL: cells[e->local_id] = NULL; cells[new_local_id] = static_cast<Cell *>(e); break;
			case FACE: faces[e->local_id] = NULL; faces[new_local_id] = static_cast<Face *>(e); break;
			case EDGE: edges[e->local_id] = NULL; edges[new_local_id] = static_cast<Edge *>(e); break;
			case NODE: nodes[e->local_id] = NULL; nodes[new_local_id] = static_cast<Node *>(e); break;
			case ESET: sets[e->local_id] = NULL;  sets[new_local_id] = static_cast<ElementSet *>(e); break;
		}
		e->local_id = new_local_id;
	}
	
#define DOWNSIZE_FACTOR 2

	int CompareElementsForReorderApply(const void * a, const void * b, void * udata)
	{
		return (*(Storage **)a)->Integer(*(Tag *)udata) - (*(Storage **)b)->Integer(*(Tag *)udata);

	}

	
	void Mesh::ReorderApply(Tag index, ElementType mask)
	{
		ReorderEmpty(mask);
		if( mask & ESET )
		{
			sort(&sets[0],sets.size(),sizeof(ElementSet *),CompareElementsForReorderApply,SwapElement,&index);
		}
		if( mask & CELL )
		{
			sort(&cells[0],cells.size(),sizeof(Cell *),CompareElementsForReorderApply,SwapElement,&index);
		}
		if( mask & FACE )
		{
			sort(&faces[0],faces.size(),sizeof(Face *),CompareElementsForReorderApply,SwapElement,&index);
		}
		if( mask & EDGE )
		{
			sort(&edges[0],edges.size(),sizeof(Edge *),CompareElementsForReorderApply,SwapElement,&index);
		}
		if( mask & NODE )
		{
			sort(&nodes[0],nodes.size(),sizeof(Node *),CompareElementsForReorderApply,SwapElement,&index);
		}
		return;

		throw NotImplemented;
		/*
		if( mask & CELL )
		{
			cells_container new_cells(cells.size(),static_cast<Cell *>(NULL));
			for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); ++it)
				new_cells[it->Integer(index)] = &*it;
			cells.swap(new_cells);
			empty_cells.clear();
			for(unsigned k = 0; k < cells.size(); ++k)
				if( cells[k] == NULL ) empty_cells.push_back(k);
		}
		
		if( mask & FACE )
		{
			faces_container new_faces(faces.size(),static_cast<Face *>(NULL));
			for(Mesh::iteratorFace it = BeginFace(); it != EndFace(); ++it)
				new_faces[it->Integer(index)] = &*it;
			faces.swap(new_faces);
			empty_faces.clear();
			for(unsigned k = 0; k < faces.size(); ++k)
				if( faces[k] == NULL ) empty_faces.push_back(k);
		}
		
		if( mask & EDGE )
		{
			edges_container new_edges(edges.size(),static_cast<Edge *>(NULL));
			for(Mesh::iteratorEdge it = BeginEdge(); it != EndEdge(); ++it)
				new_edges[it->Integer(index)] = &*it;
			edges.swap(new_edges);
			empty_edges.clear();
			for(unsigned k = 0; k < edges.size(); ++k)
				if( edges[k] == NULL ) empty_edges.push_back(k);
		}
		if( mask & NODE )
		{
			nodes_container new_nodes(nodes.size(),static_cast<Node *>(NULL));
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
				new_nodes[it->Integer(index)] = &*it;
			nodes.swap(new_nodes);
			empty_nodes.clear();
			for(unsigned k = 0; k < nodes.size(); ++k)
				if( nodes[k] == NULL ) empty_nodes.push_back(k);
		}
		*/
	}
	
	void Mesh::ReorderEmpty(ElementType etype)
	{
		//~ Tag temp = CreateTag("REORDER_EMPTY_OLD_POSITION",DATA_REFERENCE,etype,NONE,1);
		//~ for(Mesh::iteratorElement it = BeginElement(etype); it != EndElement(); ++it) 
			//~ it->ReferenceDF(temp) = &*it;
		if( etype & CELL )
		{
			int cend = cells.size()-1;
			while( cend >= 0 && !empty_cells.empty())
			{
				if( cells[cend] != NULL )
				{
					if( static_cast<int>(empty_cells.back()) < cend )
					{
						MoveStorage(cells[cend],empty_cells.back());
						cend--;
					}
					empty_cells.pop_back();
				}
				else cend--;
			}
			double factor = static_cast<double>(cells.capacity())/static_cast<double>(cend);
			if( cend >= 0 )
			{
				//std::cout << " resize cells " << cells.size() << " to " << cend+1 << std::endl;
				cells.resize(cend+1);
				empty_cells.clear();
				if( factor > DOWNSIZE_FACTOR )
				{
					//~ cells_container(cells).swap(cells);
					//~ empty_container(empty_cells).swap(empty_cells);
					
					//Downsize data
					//std::cout << "Downsize cell data" << std::endl;
					ReallocateData(CELL);

				}
			}
			else
			{
				cells.clear();
				empty_cells.clear();
			}
			//printf("factor %g cells size %ld capacity %ld\n",factor,cells.size(),cells.capacity());
		}
		if( etype & FACE )
		{
			int cend = faces.size()-1;
			while( cend >= 0 && !empty_faces.empty())
			{
				if( faces[cend] != NULL )
				{
					if( static_cast<int>(empty_faces.back()) < cend )
					{
						MoveStorage(faces[cend],empty_faces.back());
						cend--;
					}
					empty_faces.pop_back();
				}
				else cend--;
			}
			double factor = static_cast<double>(faces.capacity())/static_cast<double>(cend);
			if( cend >= 0 )
			{
				//std::cout << " resize faces " << cells.size() << " to " << cend+1 << std::endl;
				faces.resize(cend+1);
				empty_faces.clear();
				if( factor > DOWNSIZE_FACTOR )
				{
					//~ faces_container(faces).swap(faces);
					//~ empty_container(empty_faces).swap(empty_faces);
					//Downsize data
					//std::cout << "Downsize face data" << std::endl;
					ReallocateData(FACE);
				}
			}
			else
			{
				faces.clear();
				empty_faces.clear();
			}
			//printf("factor %g faces size %ld capacity %ld\n",factor,faces.size(),faces.capacity());
		}
		if( etype & EDGE )
		{
			int cend = edges.size()-1;
			while( cend >= 0 && !empty_edges.empty())
			{
				if( edges[cend] != NULL )
				{
					if( static_cast<int>(empty_edges.back()) < cend )
					{
						MoveStorage(edges[cend],empty_edges.back());
						cend--;
					}
					empty_edges.pop_back();
				}
				else cend--;
			}
			double factor = static_cast<double>(edges.capacity())/static_cast<double>(cend);
			if( cend >= 0 )
			{
				//std::cout << " resize edges " << cells.size() << " to " << cend+1 << std::endl;
				edges.resize(cend+1);
				empty_edges.clear();
				if( factor > DOWNSIZE_FACTOR )
				{
					//~ edges_container(edges).swap(edges);
					//~ empty_container(empty_edges).swap(empty_edges);
					
					//Downsize data
					//std::cout << "Downsize edge data" << std::endl;
					ReallocateData(EDGE);
				}
			}
			else
			{
				edges.clear();
				empty_edges.clear();
			}
			//printf("factor %g edges size %ld capacity %ld\n",factor,edges.size(),edges.capacity());
		}
		if( etype & NODE )
		{
			int cend = nodes.size()-1;
			while( cend >= 0 && !empty_nodes.empty())
			{
				if( nodes[cend] != NULL )
				{
					if( static_cast<int>(empty_nodes.back()) < cend )
					{
						MoveStorage(nodes[cend],empty_nodes.back());
						cend--;
					}
					empty_nodes.pop_back();
				}
				else cend--;
			}
			double factor = static_cast<double>(nodes.capacity())/static_cast<double>(cend);
			if( cend >= 0 )
			{
				//std::cout << " resize nodes " << cells.size() << " to " << cend+1 << std::endl;
				nodes.resize(cend+1);
				empty_nodes.clear();
				if( factor > DOWNSIZE_FACTOR )
				{
					//~ nodes_container(nodes).swap(nodes);
					//~ empty_container(empty_nodes).swap(empty_nodes);
					
					//Downsize data
					//std::cout << "Downsize node data" << std::endl;
					ReallocateData(NODE);
				}
			}
			else
			{
				nodes.clear();
				empty_nodes.clear();
			}
			//printf("factor %g nodes size %ld capacity %ld\n",factor,nodes.size(),nodes.capacity());
		}
		if( etype & ESET )
		{
			int cend = sets.size()-1;
			while( cend >= 0 && !empty_sets.empty())
			{
				if( sets[cend] != NULL )
				{
					if( static_cast<int>(empty_sets.back()) < cend )
					{
						MoveStorage(sets[cend],empty_sets.back());
						cend--;
					}
					empty_sets.pop_back();
				}
				else cend--;
			}
			double factor = static_cast<double>(sets.capacity())/static_cast<double>(cend);
			if( cend >= 0 )
			{
				sets.resize(cend+1);
				empty_sets.clear();
				if( factor > DOWNSIZE_FACTOR )
				{
					//~ sets_container(sets).swap(sets);
					//~ empty_container(empty_sets).swap(empty_sets);
					//Downsize data
					ReallocateData(ESET);
				}
			}
			else
			{
				sets.clear();
				empty_sets.clear();
			}
			//printf("factor %g sets size %ld capacity %ld\n",factor,sets.size(),sets.capacity());
		}
		//renew references and sets
		//~ {
			//~ std::vector< std::pair<Element *,Element *> > mapping;
			//~ for(Mesh::iteratorElement it = BeginElement(etype); it != EndElement(); ++it) 
				//~ mapping.push_back(std::pair<Element *,Element *>(it->ReferenceDF(temp),&*it));
			//~ std::sort(mapping.begin(),mapping.end());
			//~ for(Mesh::iteratorTag it = BeginTag(); it != EndTag(); ++it)
				//~ if( it->GetDataType() == DATA_REFERENCE )
				//~ {
					//~ for(ElementType ietype = NODE; ietype <= CELL; ietype = ietype << 1)
						//~ if( it->isDefined(ietype) )
						//~ {
							//~ if( it->isSparse(ietype) )
							//~ {
								//~ for(Mesh::iteratorElement jt = BeginElement(ietype); jt != EndElement(); ++jt)
									//~ if( jt->HaveData(*it) )
									//~ {
										//~ Storage::reference_array arr = jt->ReferenceArray(*it);
										//~ for(Storage::reference_array::iterator qt = arr.begin(); qt != arr.end(); ++qt)
											//~ if( (*qt)->GetElementType() & etype )
											//~ {
												//~ std::vector< std::pair<Element *,Element *> >::iterator find = std::lower_bound(mapping.begin(),mapping.end(),std::pair<Element *,Element *>(*qt,NULL));
												//~ assert(*qt == find->first);
												//~ *qt = find->second;
											//~ }
									//~ }
							//~ }
							//~ else
							//~ {
								//~ for(Mesh::iteratorElement jt = BeginElement(etype); jt != EndElement(); ++jt)
								//~ {
									//~ Storage::reference_array arr = jt->ReferenceArray(*it);
									//~ for(Storage::reference_array::iterator qt = arr.begin(); qt != arr.end(); ++qt)
										//~ if( (*qt)->GetElementType() & etype )
										//~ {
											//~ std::vector< std::pair<Element *,Element *> >::iterator find = std::lower_bound(mapping.begin(),mapping.end(),std::pair<Element *,Element *>(*qt,NULL));
											//~ assert(*qt == find->first);
											//~ *qt = find->second;
										//~ }
								//~ }
							//~ }
						//~ }
				//~ }
			//~ for(Mesh::iteratorSet it = BeginSet(); it != EndSet(); it++)
			//~ {
				//~ ElementSet add;
				//~ ElementSet::iterator jt = it->begin();
				//~ while(jt != it->end())
					//~ if( jt->GetElementType() & etype )
					//~ {
						//~ std::vector< std::pair<Element *,Element *> >::iterator find = std::lower_bound(mapping.begin(),mapping.end(),std::pair<Element *,Element *>(&*jt,NULL));
						//~ assert(&*jt == find->first);
						//~ if( find->first != find->second )
						//~ {
							//~ it->Erase(jt++);
							//~ add.Insert(find->second);
						//~ }
						//~ else jt++;
					//~ }
				//~ if( !add.empty() ) it->Union(add);
			//~ }
		//~ }
		//~ DeleteTag(temp,etype);
	}

	MarkerType Mesh::CreateMarker()
	{
#if defined(NEW_MARKERS)
		Storage::bulk * marker_space = static_cast<Storage::bulk * >(GetDenseLink(GetMeshLink()->MarkersTag()));
		INMOST_DATA_ENUM_TYPE ret;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < MarkerFields; ++k)
		{
			Storage::bulk mask = ((~marker_space[k]) & (-(~marker_space[k])));
			if( mask )
			{
				ret = (k << MarkerShift) | mask;
				marker_space[k] |= mask;
				return ret;
			}
		}
#else
		MarkerType marker_space = markers;
		MarkerType ret = ((~marker_space) & (-(~marker_space)));
		if( ret )
		{
			SetMarker(ret);
			return ret;
		}
#endif
		throw NoSpaceForMarker;
	}
	void Mesh::ReleaseMarker(MarkerType n)
	{
		RemMarker(n);
	}
	
	INMOST_DATA_ENUM_TYPE Mesh::GetArrayCapacity(ElementType etype)
	{
		assert(OneType(etype));
		switch(etype)
		{
			case NODE: return nodes.capacity(); break;
			case EDGE: return edges.capacity(); break;
			case FACE: return faces.capacity(); break;
			case CELL: return cells.capacity(); break;
			case ESET: return sets.capacity(); break;
		}
		return 0;
	}

	
	bool Mesh::isOriginal(Element * e)
	{
		switch(e->GetElementType())
		{
			case CELL: if( cells[e->LocalID()] == e ) return true; else return false;
			case FACE: if( faces[e->LocalID()] == e ) return true; else return false;
			case EDGE: if( edges[e->LocalID()] == e ) return true; else return false;
			case NODE: if( nodes[e->LocalID()] == e ) return true; else return false;
		}
		return false;
	}
}
#endif
