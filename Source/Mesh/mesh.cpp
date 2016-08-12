#include "inmost.h"
#include <sstream>
#if defined(USE_MESH)
#define WAITNL 	{char c;scanf("%c",&c);}

namespace INMOST
{
  static std::vector<Mesh *> allocated_meshes;


#if defined(USE_PARALLEL_WRITE_TIME)
	void Mesh::AtExit(void)
	{
		while(!allocated_meshes.empty())
		{
			if( allocated_meshes.back() != NULL )
			{
				allocated_meshes.back()->FinalizeFile();
			}
			allocated_meshes.pop_back();
		}
	}
#endif //USE_PARALLEL_WRITE_TIME

  std::string Mesh::GetMeshName()
  {
    return name;
  }
  void Mesh::SetMeshName(std::string new_name)
  {
    name = new_name;
  }
  Mesh * Mesh::GetMesh(std::string name)
  {
    for(int q = 0; q < (int)allocated_meshes.size(); ++q)
      if( allocated_meshes[q]->GetMeshName() == name )
        return allocated_meshes[q];
    return NULL;
  }

	const char * TopologyCheckNotifyString(TopologyCheck c)
	{
		switch(c)
		{
			case THROW_EXCEPTION:        return "exception thrown";
			case PRINT_NOTIFY:           return "print notify";
			case DELETE_ON_ERROR:        return "element should be deleted on error";
			case MARK_ON_ERROR:          return "element is marked on error";
			case DUPLICATE_EDGE:         return "search for duplicate edge";
			case DUPLICATE_FACE:         return "search for duplicate face";
			case DUPLICATE_CELL:         return "search for duplicate cell";
			case DEGENERATE_EDGE:        return "TOPOLOGY ERROR: curvilinear edge found"; 
			case DEGENERATE_FACE:        return "TOPOLOGY ERROR: degenerate face found";
			case DEGENERATE_CELL:        return "TOPOLOGY ERROR: degenerate cell found";
			case FACE_ORIENTATION:       return "TOPOLOGY ERROR: bad face orientation";
			case FACE_PLANARITY:         return "TOPOLOGY ERROR: non-planar face found";
			case INTERLEAVED_FACES:      return "TOPOLOGY ERROR: interleaving faces found";
			case TRIPLE_SHARED_FACE:     return "TOPOLOGY ERROR: face have more then two neighbours"; 
			case FLATTENED_CELL:         return "TOPOLOGY ERROR: flattened cell found"; 
			case ADJACENT_DUPLICATE:     return "TOPOLOGY ERROR: duplicates in adjacent elements";
			case ADJACENT_HIDDEN:        return "TOPOLOGY ERROR: hidden element is used as adjacent"; 
			case ADJACENT_VALID:         return "TOPOLOGY ERROR: invalid handle is used as adjacent"; 
			case ADJACENT_DIMENSION:     return "TOPOLOGY ERROR: wrong dimension of adjacent elements";
			case PROHIBIT_MULTILINE:     return "TOPOLOGY ERROR: multiline is prohibited"; 
			case PROHIBIT_POLYGON:       return "TOPOLOGY ERROR: polygon is prohibited"; 
			case PROHIBIT_MULTIPOLYGON:  return "TOPOLOGY ERROR: multipolygon is prohibited"; 
			case PROHIBIT_POLYHEDRON:    return "TOPOLOGY ERROR: polyhedron is prohibited"; 
			case FACE_EDGES_ORDER:       return "TOPOLOGY ERROR: no order in face edges"; 
			case PROHIBIT_CONCAVE_FACE:  return "TOPOLOGY ERROR: concave faces are prohibited"; 
			case PROHIBIT_CONCAVE_CELL:  return "TOPOLOGY ERROR: concave cells are prohibited"; 
			case PROHIBIT_NONSTAR_FACE:  return "TOPOLOGY ERROR: non star-shaped faces are prohibited"; 
			case PROHIBIT_NONSTAR_CELL:  return "TOPOLOGY ERROR: non star-shpaed cells are prohibited"; 
			case FACE_SELF_INTERSECTION: return "TOPOLOGY ERROR: self intersection of face edges detected"; 
			case CELL_SELF_INTERSECTION: return "TOPOLOGY ERROR: self intersection of cell faces detected"; 
			case DISABLE_2D:             return "TOPOLOGY ERROR: 2d mesh support is disabled"; 
			default: return "unknown";
		}
	}

	const char * ElementTypeName(ElementType t)
	{
		switch(t)
		{
			case NONE: return "NONE";
			case NODE: return "NODE";
			case EDGE: return "EDGE";
			case FACE: return "FACE";
			case CELL: return "CELL";
			case ESET: return "ESET";
			case MESH: return "MESH";
		}
		return "UNKNOWN";
	}

	void Mesh::Init(std::string name)
	{
#if defined(CHECKS_MARKERS)
		check_shared_mrk = check_private_mrk = false;
#endif
		m_link = this;
		integer selfid = 1;
		selfid = TieElement(5);
		assert(selfid == 0);
		dim = 3;
		have_global_id = NONE;
		checkset = DEFAULT_CHECK;
		errorset = 0;
		new_element = hide_element = 0;

		memset(hidden_count,0,sizeof(integer)*6);
		memset(hidden_count_zero,0,sizeof(integer)*6);

		memset(remember,0,sizeof(remember));
		tag_coords        = CreateTag("PROTECTED_COORD",DATA_REAL, NODE,NONE,dim);
		tag_high_conn     = CreateTag("PROTECTED_HIGH_CONN",DATA_REFERENCE,ESET|CELL|FACE|EDGE|NODE,NONE);
		tag_low_conn      = CreateTag("PROTECTED_LOW_CONN",DATA_REFERENCE,ESET|CELL|FACE|EDGE|NODE,NONE);
		tag_markers       = CreateTag("PROTECTED_MARKERS",DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFields);
		tag_geom_type     = CreateTag("PROTECTED_GEOM_TYPE",DATA_BULK,CELL|FACE|EDGE|NODE,NONE,1);
		tag_setname       = CreateTag("PROTECTED_SET_NAME",DATA_BULK,ESET,NONE);
		tag_setcomparator = CreateTag("PROTECTED_SET_COMPARATOR",DATA_BULK,ESET,NONE,1);
		AllocatePrivateMarkers();
		for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
			ReallocateData(ElementNum(etype),GetArrayCapacity(ElementNum(etype)));
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
    allocated_meshes.push_back(this);
  }
		
	Mesh::Mesh()
	:TagManager(), Storage(NULL,ComposeHandle(MESH,0))
	{
		std::stringstream tmp;
		tmp << "Mesh" << allocated_meshes.size();
		name = tmp.str();
		Init(name);
	}

	Mesh::Mesh(std::string name)
	:TagManager(), Storage(NULL,ComposeHandle(MESH,0))
	{
		Init(name);
	}

  bool Mesh::GetPrivateMarker(HandleType h, MarkerType n) const
  {
    assert(isPrivate(n));
    n &= ~MarkerPrivateBit;
    int thread = GetLocalProcessorRank();
    const bulk * mem = static_cast<const bulk *>(MGetDenseLink(h,tag_private_markers[thread]));
    return (mem[n >> MarkerShift] & static_cast<bulk>(n & MarkerMask)) != 0;
  }

  void Mesh::SetPrivateMarker(HandleType h,MarkerType n)  
  {
    assert(isPrivate(n));
    n &= ~MarkerPrivateBit;
    int thread = GetLocalProcessorRank();
    bulk * mem = static_cast<bulk *>(MGetDenseLink(h,tag_private_markers[thread]));
    mem[n >> MarkerShift] |= static_cast<bulk>(n & MarkerMask);
  }

  void Mesh::RemPrivateMarker(HandleType h,MarkerType n) 
  {
    assert(isPrivate(n));
    n &= ~MarkerPrivateBit;
    int thread = GetLocalProcessorRank();
    bulk * mem = static_cast<bulk *>(MGetDenseLink(h,tag_private_markers[thread]));
    mem[n >> MarkerShift] &= ~static_cast<bulk>(n & MarkerMask);
  }

  void Mesh::AllocatePrivateMarkers()
  {
#if defined(USE_OMP)
#pragma omp parallel
    {
#pragma omp single
      {
		  tag_private_markers = new Tag[GetLocalProcessorNumber()];
      }
//#pragma omp ordered
      {
        std::stringstream name;
		name << "PROTECTED_PRIVATE_MARKERS_" << GetLocalProcessorRank();
        tag_private_markers[GetLocalProcessorRank()] = CreateTag(name.str(),DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFieldsPrivate);
      }
    }
#else
    tag_private_markers = new Tag;
    tag_private_markers[0] = CreateTag("PROTECTED_PRIVATE_MARKERS_0",DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFieldsPrivate);
#endif
  }

  int Mesh::GetLocalProcessorNumber() const
  {
#if defined(USE_OMP)
	  return omp_get_num_threads();
#else
	  return 1;
#endif
  }

  int Mesh::GetLocalProcessorRank() const
  {
#if defined(USE_OMP)
	  return omp_get_thread_num();
#else
	  return 0;
#endif
  }

  void Mesh::DeallocatePrivateMarkers()
  {
#if defined(USE_OMP)
#pragma omp parallel
    {
		//retrive tag before it was erased
		Tag del = tag_private_markers[GetLocalProcessorRank()];
#pragma omp barrier
		//delete tag, it will erase the tag
		DeleteTag(del);
#pragma omp barrier
		//deallocate space
#pragma omp single
		{
			delete [] tag_private_markers;
		}
    }
#else
    DeleteTag(tag_private_markers[0]);
    delete tag_private_markers;
#endif
  }

	Storage::enumerator Mesh::MemoryUsage(HandleType h)
	{
		if( isValidHandle(h) )
		{
			integer etypenum = GetHandleElementNum(h);
			enumerator ret = 2*sizeof(integer); //link and address occupied
			for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); ++t)
			{
				if( t->isDefinedByDim(etypenum) )
				{
					bool have_data = true;
					if( t->isSparseByDim(etypenum) )
					{
						have_data = HaveData(h,*t);
						if( have_data ) 
							ret += sizeof(sparse_sub_record); //size occupied for storage of data link
					}
					if( have_data )
					{
						ret += t->GetRecordSize(); //all utilized data for fixed data, size of support structure for variable data
						if( t->GetSize() == ENUMUNDEF )
							ret += GetDataCapacity(h,*t); //actually occupied size for sparse data
					}
				}
				if( !sparse_data[etypenum].empty() ) ret += sizeof(sparse_sub_type); //size needed to support sparse data
			}
			//Any additional impact of supporting huge structures over all elements may be added later
			return ret;
		}
		else return 0;
	}
	
	Mesh::Mesh(const Mesh & other)
	:TagManager(other),Storage(NULL,ComposeHandle(MESH,0))
	{
		{
			std::stringstream tmp;
			tmp << other.name << "_copy";
			name = tmp.str();
		}
#if defined(CHECKS_MARKERS)
		check_shared_mrk = other.check_shared_mrk;
		check_private_mrk = other.check_private_mrk;
#endif
		m_link = this;
		integer selfid = 1;
		selfid = TieElement(5);
		assert(selfid == 0);
		//TagManager constuctor copied only tags
		//copy links:
		for(int i = 0; i < 5; i++)
		{
			links[i] = other.links[i];
			empty_links[i] = other.empty_links[i];
			empty_space[i] = other.empty_space[i];
		}
		//this should alocate space for data and copy it including markers and connections
		for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
		{
			ReallocateData(ElementNum(etype),GetArrayCapacity(ElementNum(etype)));
			for(tag_array_type::size_type i = 0; i < tags.size(); ++i)
			{
				if( tags[i].isDefined(etype) )
				{
					if( tags[i].isSparse(etype) )
					{
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) if( isValidElement(etype,lid) )
						{
							HandleType h = ComposeHandle(etype,lid);
							if( other.HaveData(h,tags[i]) )
								TagManager::CopyData(tags[i],MGetLink(h,tags[i]),other.MGetLink(h,other.tags[i]));
						}
					}
					else
					{
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) if( isValidElement(etype,lid) )
						{
							HandleType h = ComposeHandle(etype,lid);
							TagManager::CopyData(tags[i],MGetLink(h,tags[i]),other.MGetLink(h,other.tags[i]));
						}
					}
				}
			}
		}
		//setup system tags shortcuts
		dim = other.dim;
		tag_coords        = CreateTag("PROTECTED_COORD",DATA_REAL, NODE,NONE,dim);
		tag_high_conn     = CreateTag("PROTECTED_HIGH_CONN",DATA_REFERENCE,ESET|CELL|FACE|EDGE|NODE,NONE);
		tag_low_conn      = CreateTag("PROTECTED_LOW_CONN",DATA_REFERENCE,ESET|CELL|FACE|EDGE|NODE,NONE);
		tag_markers       = CreateTag("PROTECTED_MARKERS",DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFields);
		tag_geom_type     = CreateTag("PROTECTED_GEOM_TYPE",DATA_BULK,CELL|FACE|EDGE|NODE,NONE,1);
		tag_setname       = CreateTag("PROTECTED_SET_NAME",DATA_BULK,ESET,NONE);
		tag_setcomparator = CreateTag("PROTECTED_SET_COMPARATOR",DATA_BULK,ESET,NONE,1);
    AllocatePrivateMarkers();
		//copy supplimentary values
		m_state = other.m_state;
		checkset = other.checkset;
		errorset = other.errorset;
		new_element = other.new_element;
		hide_element = other.hide_element;
		epsilon = other.epsilon;
		have_global_id = other.have_global_id;
		// copy communicator
		if( m_state == Mesh::Parallel ) SetCommunicator(other.comm); else comm = INMOST_MPI_COMM_WORLD;
		// reestablish geometric tags and table
		memcpy(remember,other.remember,sizeof(remember));
		RestoreGeometricTags();
		//this is not needed as it was copied with all the other data
		//recompute global ids
		//AssignGlobalID(other.have_global_id);
    allocated_meshes.push_back(this);
	}
	
	Mesh & Mesh::operator =(Mesh const & other)
	{
		if( this == &other ) return *this; //don't do anything
		{
			std::stringstream tmp;
			tmp << other.name << "_copy";
			name = tmp.str();
		}
#if defined(CHECKS_MARKERS)
		check_shared_mrk = other.check_shared_mrk;
		check_private_mrk = other.check_private_mrk;
#endif
		//first delete everything
		//delete parallel vars
#if defined(USE_MPI)
#if defined(USE_MPI_P2P)
		if( m_state == Mesh::Parallel )
		{
			MPI_Free_mem(shared_space);
			MPI_Win_free(&window);
		}
#endif
#endif
		//clear all data fields
		for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
		{
			for(tag_array_type::size_type i = 0; i < tags.size(); ++i)
			{
				if( tags[i].isDefined(etype) )
				{
					if( tags[i].isSparse(etype) )
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer lid = 0; lid < LastLocalID(etype); ++lid)
							if( isValidElement(etype,lid) )
								DelSparseData(ComposeHandle(etype,lid),tags[i]);
					}
					else if( tags[i].GetSize() == ENUMUNDEF )
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) 
							if( isValidElement(etype,lid) )
								DelDenseData(ComposeHandle(etype,lid),tags[i]);
					}
				}
			}
		}
		//clear links
		for(int i = 0; i < 5; i++)
		{
			links[i].clear();
			empty_links[i].clear();
			empty_space[i].clear();
		}
		DeallocatePrivateMarkers();
		//this should copy tags, clear sparse data, set up dense links
		TagManager::operator =(other);
		//set up new links
		for(int i = 0; i < 5; i++)
		{
			links[i] = other.links[i];
			empty_links[i] = other.empty_links[i];
			empty_space[i] = other.empty_space[i];
		}
		//this should alocate space for data and copy it including markers and connections
		for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
		{
			ReallocateData(ElementNum(etype),GetArrayCapacity(ElementNum(etype)));
			for(tag_array_type::size_type i = 0; i < tags.size(); ++i)
			{
				if( tags[i].isDefined(etype) )
				{
					if( tags[i].isSparse(etype) )
					{
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) if( isValidElement(etype,lid) )
						{
							HandleType h = ComposeHandle(etype,lid);
							if( other.HaveData(h,tags[i]) )
								TagManager::CopyData(tags[i],MGetLink(h,tags[i]),other.MGetLink(h,other.tags[i]));
						}
					}
					else
					{
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) if( isValidElement(etype,lid) )
						{
							HandleType h = ComposeHandle(etype,lid);
							TagManager::CopyData(tags[i],MGetLink(h,tags[i]),other.MGetLink(h,other.tags[i]));
						}
					}
				}
			}
		}
		//setup system tags shortcuts
		dim = other.dim;
		tag_coords        = CreateTag("PROTECTED_COORD",DATA_REAL, NODE,NONE,dim);
		tag_high_conn     = CreateTag("PROTECTED_HIGH_CONN",DATA_REFERENCE,ESET|CELL|FACE|EDGE|NODE,NONE);
		tag_low_conn      = CreateTag("PROTECTED_LOW_CONN",DATA_REFERENCE,ESET|CELL|FACE|EDGE|NODE,NONE);
		tag_markers       = CreateTag("PROTECTED_MARKERS",DATA_BULK,CELL|FACE|EDGE|NODE|ESET|MESH,NONE,MarkerFields);
		tag_geom_type     = CreateTag("PROTECTED_GEOM_TYPE",DATA_BULK,CELL|FACE|EDGE|NODE,NONE,1);
		tag_setname       = CreateTag("PROTECTED_SET_NAME",DATA_BULK,ESET,NONE);
		tag_setcomparator = CreateTag("PROTECTED_SET_COMPARATOR",DATA_BULK,ESET,NONE,1);
    AllocatePrivateMarkers();
		//copy supplimentary values
		m_state = other.m_state;
		checkset = other.checkset;
		errorset = other.errorset;
		new_element = other.new_element;
		hide_element = other.hide_element;
		epsilon = other.epsilon;
		have_global_id = other.have_global_id;
		// copy communicator
		if( m_state == Mesh::Parallel ) SetCommunicator(other.comm); else comm = INMOST_MPI_COMM_WORLD;
		// reestablish geometric tags and table
		memcpy(remember,other.remember,sizeof(remember));
		RestoreGeometricTags();
		//this is not needed as it was copied with all the other data
		//recompute global ids
		//AssignGlobalID(other.have_global_id);
		return *this;
	}

	void Mesh::Clear()
	{
		for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
		{
			for(tag_array_type::size_type i = 0; i < tags.size(); ++i)
			{
				if( tags[i].isDefined(etype) )
				{
					if( tags[i].isSparse(etype) )
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) 
							if( isValidElement(etype,lid) )
								DelSparseData(ComposeHandle(etype,lid),tags[i]);
					}
					else if( tags[i].GetSize() == ENUMUNDEF )
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) 
							if( isValidElement(etype,lid) )
								DelDenseData(ComposeHandle(etype,lid),tags[i]);
					}
				}
			}
		}
		DeallocatePrivateMarkers();
		memset(remember,0,sizeof(bool)*15);
		tags.clear();
		//clear links
		dense_data.clear();
		for(int i = 0; i < 5; i++)
		{
			links[i].clear();
			empty_links[i].clear();
			empty_space[i].clear();
		}
		for(int i = 0; i < 6; i++)
		{
			sparse_data[i].clear();
			back_links[i].clear();
		}
		RemTopologyCheck(ENUMUNDEF);
		SetTopologyCheck(DEFAULT_CHECK);

	}
	
	Mesh::~Mesh()
	{
		//clear all data fields
		for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
		{
			for(tag_array_type::size_type i = 0; i < tags.size(); ++i)
			{
				if( tags[i].isDefined(etype) )
				{
					if( tags[i].isSparse(etype) )
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) 
							if( isValidElement(etype,lid) )
								DelSparseData(ComposeHandle(etype,lid),tags[i]);
					}
					else if( tags[i].GetSize() == ENUMUNDEF )
					{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
						for(integer lid = 0; lid < LastLocalID(etype); ++lid) 
							if( isValidElement(etype,lid) )
								DelDenseData(ComposeHandle(etype,lid),tags[i]);
					}
				}
			}
		}
		DeallocatePrivateMarkers();
		//clear links
		for(int i = 0; i < 5; i++)
		{
			links[i].clear();
			empty_links[i].clear();
			empty_space[i].clear();
		}
#if defined(USE_MPI)
#if defined(USE_MPI_P2P)
		if( m_state == Mesh::Parallel )
		{
			int test = 0;
			MPI_Finalized(&test);
			if( !test )
			{
				MPI_Free_mem(shared_space);
				MPI_Win_free(&window);
			}
			else
			{
				std::cout << "Cannot release memory and window allocated by MPI" << std::endl;
				std::cout << "since MPI was already finalized. Most likely this" << std::endl;
				std::cout << "happens when you define class Mesh in main() function" << std::endl;
				std::cout << "so that destructor get called after Mesh::Finalize()" << std::endl;
				std::cout << "Please enclose area where you use Mesh class" << std::endl;
				std::cout << "with scopes, so that destructor get called when" << std::endl;
				std::cout << "execution goes out of the scope or dynamically" << std::endl;
				std::cout << "allocate and explicitly delete Mesh class." << std::endl;
				std::cout << "Thank you!" << std::endl;
			}
			//~ MPI_Comm_free(&comm);
		}
#endif //USE_MPI_P2P
#endif //USE_MPI
#if defined(USE_PARALLEL_WRITE_TIME)
		FinalizeFile();
#endif //USE_PARALLEL_WRITE_TIME
		{
			for(size_t q = 0; q < allocated_meshes.size(); ++q)
				if (allocated_meshes[q] == this)
					allocated_meshes[q] = NULL;
			std::sort(allocated_meshes.rbegin(), allocated_meshes.rend());
			while(!allocated_meshes.empty() && allocated_meshes.back() == NULL) allocated_meshes.pop_back();
		}
		//arrays for data are deallocated inside ~TagManager()
	}
	
	
	
	Tag Mesh::CreateTag(std::string name, DataType dtype, ElementType etype,ElementType sparse, INMOST_DATA_ENUM_TYPE size)
	{
		Tag ret = TagManager::CreateTag(this,name,dtype,etype,sparse,size);
		return ret;
	}
	Tag Mesh::DeleteTag(Tag tag, ElementType type_mask)
	{
		//deallocate data on elements
		for(ElementType etype = NODE; etype <= MESH; etype = NextElementType(etype))
		{
			if( (etype & type_mask) && tag.isDefined(etype) )
			{
				if( tag.isSparse(etype) )
				{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for(integer lid = 0; lid < LastLocalID(etype); ++lid) 
						if( isValidElement(etype,lid) )
							DelSparseData(ComposeHandle(etype,lid),tag);
				}
				else if( tag.GetSize() == ENUMUNDEF )
				{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
					for(integer lid = 0; lid < LastLocalID(etype); ++lid) 
						if( isValidElement(etype,lid) )
							DelDenseData(ComposeHandle(etype,lid),tag);
				}
			}
		}
#if defined(USE_OMP)
#pragma omp critical
#endif
		{
			tag = TagManager::DeleteTag(tag,type_mask);
		}
		return tag;
	}
	
	
	
	HandleType Mesh::FindSharedAdjacency(const HandleType * arr, enumerator s) const
	{
		if( s == 0 ) return InvalidHandle();
		if( !HideMarker() )
		{
			{
				enumerator flag0, flag1, i;
				dynarray<Element::adj_type const *, 64> hcarr(s);
				hcarr[0] = &HighConn(arr[0]);
				if( !hcarr[0]->empty() ) 
				{
					for(i = 1; i < s; i++) hcarr[i] = &HighConn(arr[i]);
				}
				else return InvalidHandle();
				enumerator it, iend = static_cast<enumerator>(hcarr[0]->size()), jt, jend;
				for(it = 0; it < iend; ++it)
				{
					flag0 = 0;
					for(i = 1; i < s; ++i)
					{
						jend = static_cast<enumerator>(hcarr[i]->size());
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
					if( flag0 == s-1 ) 
					{
						return hcarr[0]->at(it);
					}
				}
			}
		}
		else
		{
			enumerator ss = Count(arr,s,HideMarker()), nextk;
			enumerator k = ENUMUNDEF, i;
			k = getNext(arr,s,k,HideMarker());
			nextk = getNext(arr,s,k,HideMarker());
			Element::adj_type const & hc0 = HighConn(arr[k]);
			Element::adj_type::size_type it, iend = hc0.size(), jt, jend, flag0, flag1;
			for(it = 0; it < iend; ++it) if( !GetMarker(hc0[it],HideMarker()) )
			{
				flag0 = 0;
				i = nextk;
				while( i < s )
				{
					flag1 = 0;
					Element::adj_type const & ihc = HighConn(arr[i]);
					jend = ihc.size();
					for(jt = 0; jt < jend; ++jt) if( !GetMarker(ihc[jt],HideMarker()) )
					{
						if( hc0[it] == ihc[jt] )
						{
							flag0++;
							flag1 = 1;
							break;
						}
					}
					if( flag1 == 0 ) break;
					i = getNext(arr,s,i,HideMarker());
				}
				if( flag0 == ss-1 ) return hc0[it];
			}
		}
		return InvalidHandle();
	}
	
	
	TopologyCheck Mesh::BeginTopologyCheck(ElementType etype, const HandleType * adj, enumerator s)
	{
		enumerator i,j,d = ENUMUNDEF;
		TopologyCheck chk = 0;
		if( GetTopologyCheck(ADJACENT_VALID) )
		{
			for(i = 0; i < s; i++)
				if( !isValidHandle(adj[i]) )
				{
					chk |= ADJACENT_VALID;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_VALID) << std::endl;
				}
		}
		if( GetTopologyCheck(ADJACENT_HIDDEN) )
		{
			for(i = 0; i < s; i++)
				if( GetMarker(adj[i],HideMarker()) )
				{
					chk |= ADJACENT_HIDDEN;
					if( GetTopologyCheck(PRINT_NOTIFY) ) std::cerr << TopologyCheckNotifyString(ADJACENT_HIDDEN) << std::endl;
				}
		}
		if( GetTopologyCheck(ADJACENT_DUPLICATE) )
		{
			bool have_dup = false;
			MarkerType dup = CreatePrivateMarker();
			for(i = 0; i < s; i++)
			{
				if( GetPrivateMarker(adj[i],dup) )
				{
					have_dup = true; // duplication of element
					break;
				}
				else SetPrivateMarker(adj[i],dup);
			}
			for(i = 0; i < s; i++) RemPrivateMarker(adj[i],dup);
			ReleasePrivateMarker(dup);
			
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
				j = ElementByHandle(adj[i])->GetElementDimension();
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
					if( ElementByHandle(adj[i])->nbAdjElements(CELL) == 2 )
					{
						check = true;
						break;
					}
				}
				if( check )
				{
					bool happen = true;
					if( (GetTopologyCheck(DUPLICATE_CELL)) && (FindSharedAdjacency(adj,s) != InvalidHandle()) ) happen = false;
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
				ElementArray<Node> allnodes(this);
				for(i = 0; i < s; i++) allnodes.Unite(ElementByHandle(adj[i])->getNodes());
				ElementArray<Face> faces = allnodes[0].getFaces();
				for(i = 1; i < static_cast<enumerator>(allnodes.size()) && !faces.empty(); i++) faces.Intersect(allnodes[i].getFaces());
				for(i = 0; i < static_cast<enumerator>(faces.size()); i++)
				{
					if( faces[i]->nbAdjElements(NODE) != allnodes.size() )
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
	
	TopologyCheck Mesh::EndTopologyCheck(HandleType he)
	{
		TopologyCheck chk = 0;
		Element e = ElementByHandle(he);
		switch(e->GetGeometricDimension(e->GetGeometricType()))
		{
		case 3:
			if( GetTopologyCheck(FLATTENED_CELL) )
			{
				ElementArray<Face> faces = e->getFaces();
				enumerator num = e->nbAdjElements(NODE);
				for(enumerator i = 0; i < static_cast<enumerator>(faces.size()); i++)
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
				ElementArray<Face> faces = e->getFaces();
				for(enumerator i = 0; i < static_cast<enumerator>(faces.size()); i++)
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
	
	Node Mesh::CreateNode(const real * coords)
	{
		integer id = TieElement(0);
		HandleType h = ComposeHandleNum(0,id);
		SetGeometricType(h,Element::Vertex);
		real * v = static_cast<Storage::real *>(MGetDenseLink(h,CoordsTag()));
		for(integer i = 0; i < dim; i++) v[i] = coords[i];
		SetMarker(h,NewMarker());
		return Node(this,h);
	}

	void PrintHandle(HandleType h)
	{
		std::cout << ElementTypeName(GetHandleElementType(h)) << " " << GetHandleID(h);
	}

	void PrintAdjElem(Element::adj_type const & arr)
	{
		for(Element::adj_type::size_type it = 0; it < arr.size(); ++it)
		{
			PrintHandle(arr[it]);
			std::cout << " ";
		}
		std::cout << std::endl;
	}

	void PrintHandles(const HandleType * beg, const HandleType * end)
	{
		for(const HandleType * it = beg; it != end; ++it)
		{
			PrintHandle(*it);
			std::cout << " ";
		}
		std::cout << std::endl;
	}

	std::pair<Edge,bool> Mesh::CreateEdge(const ElementArray<Node> & nodes)
	{
		HandleType he = InvalidHandle();
		if( !nodes.empty() )
		{
			TopologyCheck chk = 0;
			chk |= BeginTopologyCheck(EDGE,nodes.data(),static_cast<enumerator>(nodes.size()));
			if( chk != 0 )
			{
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
			if( GetTopologyCheck() & DUPLICATE_EDGE )
			{
				HandleType test = FindSharedAdjacency(nodes.data(),static_cast<enumerator>(nodes.size()));
				if (test != InvalidHandle()) return std::make_pair(Edge(this,test),false);
			}
			integer id = TieElement(1);
			he = ComposeHandleNum(1,id);
			for(ElementArray<Node>::size_type i = 0; i < nodes.size(); i++)
			{
				Element::adj_type & hc = HighConn(nodes.at(i));

				//PrintHandle(nodes.at(i)); std::cout << " current: "; PrintAdjElem(hc); 
#if defined(USE_OMP)
#pragma omp critical (node_high_conn)
#endif
				{
					hc.push_back(he);
				}

				//std::cout << "add "; PrintHandle(he); std::cout << " result: "; PrintAdjElem(hc);
			}
			Element::adj_type & lc = LowConn(he);

			//PrintHandle(he); std::cout << " current: "; PrintAdjElem(lc);
			
			lc.insert(lc.end(),nodes.data(),nodes.data()+nodes.size());

			//std::cout << "add "; PrintHandles(nodes.data(),nodes.data()+nodes.size());
			//std::cout << "result: "; PrintAdjElem(lc);
			//DEBUG
			/*
			bool halt = false;
			if(!Element(this,he)->CheckElementConnectivity()) 
			{
				Element(this,he)->PrintElementConnectivity();
				halt = true;
			}
			for(ElementArray<Node>::size_type i = 0; i < nodes.size(); i++) 
				if(!nodes[i]->CheckElementConnectivity()) 
				{
					nodes[i]->PrintElementConnectivity();
					halt = true;
				}
			assert(!halt);
			*/
			//DEBUG END
			ComputeGeometricType(he);
			SetMarker(he,NewMarker());
			RecomputeGeometricData(he);
			chk |= EndTopologyCheck(he);
			if( chk != 0 )
			{
				if( GetTopologyCheck(MARK_ON_ERROR) ) Integer(he,TopologyErrorTag()) = chk;
				if( GetTopologyCheck(DELETE_ON_ERROR) ) { Destroy(he); he = InvalidHandle();}
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
		}
		return std::make_pair(Edge(this,he),true);
	}
	
	
	std::pair<Face,bool> Mesh::CreateFace(const ElementArray<Node> & f_nodes)
	{
		ElementArray<Edge> f_edges(this,f_nodes.size());
		ElementArray<Node> e_nodes(this,2);
		if( f_nodes.size() == 2 ) //This is an edge for 2d!
		{
			e_nodes.resize(1);
			e_nodes.at(0) = f_nodes.at(0);
			f_edges.at(0) = CreateEdge(e_nodes).first->GetHandle();
			e_nodes.at(0) = f_nodes.at(1);
			f_edges.at(1) = CreateEdge(e_nodes).first->GetHandle();
		}
		else
		{
			for(ElementArray<Node>::size_type i = 0; i < f_nodes.size(); i++)
			{
				e_nodes.at(0) = f_nodes.at(i);
				e_nodes.at(1) = f_nodes.at((i+1)%f_nodes.size());
				f_edges.at(i) = CreateEdge(e_nodes).first->GetHandle();
			}
		}
		return CreateFace(f_edges);
	}
	
	std::pair<Face,bool> Mesh::CreateFace(const ElementArray<Edge> & f_edges)
	{
		HandleType he = InvalidHandle();
		if( !f_edges.empty() )
		{
			TopologyCheck chk = 0;
			chk |= BeginTopologyCheck(FACE,f_edges.data(),static_cast<enumerator>(f_edges.size()));
			if( chk != 0 )
			{
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
			if( GetTopologyCheck(DUPLICATE_FACE) )
			{
				HandleType test = FindSharedAdjacency(f_edges.data(),static_cast<enumerator>(f_edges.size()));
				if (test != InvalidHandle()) return std::make_pair(Face(this,test),false);
			}
			integer id = TieElement(2);
			he = ComposeHandleNum(2,id);
			for(ElementArray<Edge>::size_type i = 0; i < f_edges.size(); i++)
			{
				Element::adj_type & hc = HighConn(f_edges.at(i));

				//PrintHandle(f_edges.at(i)); std::cout << " current: "; PrintAdjElem(hc); 
#if defined(USE_OMP)
#pragma omp critical (edge_high_conn)
#endif
				hc.push_back(he);

				//std::cout << "add "; PrintHandle(he); std::cout << " result: "; PrintAdjElem(hc);
			}
			Element::adj_type & lc = LowConn(he);

			//PrintHandle(he); std::cout << "current: "; PrintAdjElem(lc);

			lc.insert(lc.end(),f_edges.data(),f_edges.data()+f_edges.size());

			//std::cout << "add "; PrintHandles(f_edges.data(),f_edges.data()+f_edges.size());
			//std::cout << "result: "; PrintAdjElem(lc);
			//DEBUG
			/*
			bool halt = false;
			if(!Element(this,he)->CheckElementConnectivity())
			{
				Element(this,he)->PrintElementConnectivity();
				halt = true;
			}
			for(ElementArray<Edge>::size_type i = 0; i < f_edges.size(); i++) 
			{
				if(!f_edges[i]->CheckElementConnectivity())
				{
					f_edges[i]->PrintElementConnectivity();
					halt = true;
				}
			}
			assert(!halt);
			*/
			//DEBUG END
			ComputeGeometricType(he);
			SetMarker(he,NewMarker());
			RecomputeGeometricData(he);
			chk |= EndTopologyCheck(he);
			if( chk != 0 )
			{
				if( GetTopologyCheck(MARK_ON_ERROR)   ) Integer(he,TopologyErrorTag()) = chk;
				if( GetTopologyCheck(DELETE_ON_ERROR) ) { Destroy(he); he = InvalidHandle();}
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
		}
		return std::make_pair(Face(this,he),true);
	}
	
	
	std::pair<Cell,bool> Mesh::CreateCell(const ElementArray<Node> & c_f_nodes, const integer * c_f_sizes, integer s, const ElementArray<Node> & suggest_nodes_order)
	{
		ElementArray<Face> c_faces(this,s);
		ElementArray<Node>::size_type j = 0;
		for(integer i = 0; i < s; i++)
		{
			c_faces.at(i) = CreateFace(ElementArray<Node>(this, c_f_nodes.data()+j, c_f_nodes.data()+j + c_f_sizes[i])).first->GetHandle();
			j += c_f_sizes[i];
		}
		return CreateCell(c_faces,suggest_nodes_order);
	}
	

	std::pair<Cell, bool> Mesh::CreateCell(const ElementArray<Node> & c_f_nodes, const integer * c_f_nodeinds, const integer * c_f_numnodes, integer s, const ElementArray<Node> & suggest_nodes_order)
	{
		integer j = 0;
		ElementArray<Node> temp(this);
		ElementArray<Face> c_faces(this,s);
		for(integer i = 0; i < s; i++)
		{
			temp.resize(c_f_numnodes[i]);
			for(integer k = j; k < j+c_f_numnodes[i]; k++)
				temp.at(k-j) = c_f_nodes.at(c_f_nodeinds[k]);
			c_faces.at(i) = CreateFace(temp).first->GetHandle();
			j += c_f_numnodes[i];
			temp.clear();
		}
		return CreateCell(c_faces,suggest_nodes_order);
	}
	
	
	void Mesh::RestoreCellNodes(HandleType hc, ElementArray<Node> & ret)
	{
		ret.clear();
		Cell c(this,hc);
		switch(GetGeometricType(hc))
		{
			case Element::Vertex:
			case Element::Line:
			{
				assert( false );
			}
			case Element::MultiLine:
			case Element::Tri:
			case Element::Quad:
			case Element::Polygon:
			{
				ElementArray<Edge> edges = c->getEdges();
				ret.reserve(edges.size());
				for(ElementArray<Edge>::iterator it = edges.begin(); it != edges.end(); it++)
				{
					ElementArray<Node> nodes = it->getNodes();
					assert( nodes.size() == 1 );
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
				MarkerType cemrk = CreatePrivateMarker();
				MarkerType femrk = CreatePrivateMarker();
				ElementArray<Face> faces = c->getFaces();
				Face face = faces[0];
				ret.reserve(8);
				ElementArray<Node> verts = face->getNodes();
				if( face->BackCell() == c )
					ret.insert(ret.end(),verts.begin(),verts.end());
				else
					ret.insert(ret.end(),verts.rbegin(),verts.rend());
				ElementArray<Edge> c_edges = c->getEdges();
				ElementArray<Edge> f_edges = face->getEdges();
				c_edges.SetPrivateMarker(cemrk);
				f_edges.SetPrivateMarker(femrk);
				for(unsigned int k = 0; k < 4; k++)
				{
					ElementArray<Edge> v_edges = ret[k]->getEdges(cemrk);
					for(ElementArray<Edge>::iterator it = v_edges.begin(); it != v_edges.end(); it++)
					{
						if( !it->GetPrivateMarker(femrk) )
						{
							if( it->getBeg() == ret[k] )
								ret.push_back(it->getEnd());
							else
								ret.push_back(it->getBeg());
							break;
						}
					}
				}
				c_edges.RemPrivateMarker(cemrk);
				f_edges.RemPrivateMarker(femrk);
				ReleasePrivateMarker(cemrk);
				ReleasePrivateMarker(femrk);
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
				MarkerType cemrk = CreatePrivateMarker();
				MarkerType femrk = CreatePrivateMarker();
				ret.reserve(6);
				Face face;
				ElementArray<Face> faces = c->getFaces();
				for(ElementArray<Face>::size_type i = 0; i < faces.size(); i++) //iterate over faces
					if( faces[i].nbAdjElements(EDGE) == 3 ) //number of edges in i-th face
					{
						face = faces[i];
						break;
					}
				ElementArray<Node> verts = face->getNodes();
				if( face->BackCell() == c )
					ret.insert(ret.end(),verts.begin(),verts.end());
				else
					ret.insert(ret.end(),verts.rbegin(),verts.rend());
				ElementArray<Edge> c_edges = c->getEdges();
				ElementArray<Edge> f_edges = face->getEdges();
				c_edges.SetPrivateMarker(cemrk);
				f_edges.SetPrivateMarker(femrk);
				for(ElementArray<Cell>::size_type k = 0; k < 3; k++)
				{
					ElementArray<Edge> v_edges = ret[k]->getEdges(cemrk);
					for(ElementArray<Edge>::iterator it = v_edges.begin(); it != v_edges.end(); it++)
						if( !it->GetPrivateMarker(femrk) )
						{
							if( it->getBeg() == ret[k] )
								ret.push_back(it->getEnd());
							else
								ret.push_back(it->getBeg());
							break;
						}
				}
				c_edges.RemPrivateMarker(cemrk);
				f_edges.RemPrivateMarker(femrk);
				ReleasePrivateMarker(cemrk);
				ReleasePrivateMarker(femrk);
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
				Face quad, triangle;
				MarkerType mrk = CreatePrivateMarker();
				ElementArray<Face> faces = c->getFaces();
				for(ElementArray<Face>::size_type i = 0; i < faces.size(); i++) //go over faces
				{
					if( faces[i].nbAdjElements(EDGE) == 4 ) //check if number of edges = 4
					{
						quad = faces[i];
						break;
					}
				}
				for(ElementArray<Face>::size_type i = 0; i < faces.size(); i++) //go over faces
				{
					if( faces[i].nbAdjElements(EDGE) == 3 ) //check if number of edges = 3
					{
						triangle = faces[i];
						break;
					}
				}
				ElementArray<Node> base_nodes = quad->getNodes();
				if( quad->BackCell() == c )
					for(ElementArray<Node>::iterator it = base_nodes.begin(); it != base_nodes.end(); it++)
					{
						ret.push_back(*it);
						it->SetPrivateMarker(mrk);
					}
				else
					for(ElementArray<Node>::reverse_iterator it = base_nodes.rbegin(); it != base_nodes.rend(); it++)
					{
						ret.push_back(*it);
						it->SetPrivateMarker(mrk);
					}
				ElementArray<Node> tri_nodes = triangle->getNodes();
				for(ElementArray<Node>::iterator it = tri_nodes.begin(); it != tri_nodes.end(); it++)
				{
					if( !it->GetPrivateMarker(mrk) )
					{
						ret.push_back(*it);
						break;
					}
				}
				for(ElementArray<Node>::iterator it = ret.begin(); it != ret.end(); it++)
					it->RemPrivateMarker(mrk);
				ReleasePrivateMarker(mrk);
				break;
			}
			default: //Tet, MultiPolygon, Polyhedron
			{
				MarkerType mrk = CreatePrivateMarker();
				if( !HideMarker() )
				{
					Element::adj_type & lc = LowConn(hc);
					for(Element::adj_type::size_type i = 0; i < lc.size(); i++) //iterate over faces
					{
						Element::adj_type & ilc = LowConn(lc[i]);
						for(Element::adj_type::size_type j = 0; j < ilc.size(); j++) //iterate over face edges
						{
							if( !GetPrivateMarker(ilc[j],mrk) )
							{
								Element::adj_type & jlc = LowConn(ilc[j]);
								for(Element::adj_type::size_type k = 0; k < jlc.size(); k++) //iterator over edge nodes
								{
									if( !GetPrivateMarker(jlc[k],mrk) )
									{
										SetPrivateMarker(jlc[k],mrk);
										ret.push_back(jlc[k]);
									}
								}
								SetPrivateMarker(ilc[j],mrk);
							}
						}
					}
					for(Element::adj_type::size_type i = 0; i < lc.size(); i++) //iterate over faces
					{
						Element::adj_type & ilc = LowConn(lc[i]);
						for(Element::adj_type::size_type j = 0; j < ilc.size(); j++) //iterate over face edges
							RemPrivateMarker(ilc[j],mrk);
					}
				}
				else
				{
					Element::adj_type & lc = LowConn(hc);
					for(Element::adj_type::size_type i = 0; i < lc.size(); i++) if( !Hidden(lc[i]) )  //iterate over faces
					{
						Element::adj_type & ilc = LowConn(lc[i]);
						for(Element::adj_type::size_type j = 0; j < ilc.size(); j++) if( !Hidden(ilc[j]) ) //iterate over face edges
						{
							if( !GetPrivateMarker(ilc[j],mrk) )
							{
								Element::adj_type & jlc = LowConn(ilc[j]);
								for(Element::adj_type::size_type k = 0; k < jlc.size(); k++) if( !Hidden(jlc[k]) ) //iterator over edge nodes
								{
									if( !GetPrivateMarker(jlc[k], mrk) )
									{
										SetPrivateMarker(jlc[k],mrk);
										ret.push_back(jlc[k]);
									}
								}
								SetPrivateMarker(ilc[j],mrk);
							}
						}
					}
					for(Element::adj_type::size_type i = 0; i < lc.size(); i++) if( !Hidden(lc[i]) )  //iterate over faces
					{
						Element::adj_type & ilc = LowConn(lc[i]);
						for(Element::adj_type::size_type j = 0; j < ilc.size(); j++) if( !Hidden(ilc[j]) ) //iterate over face edges
							RemPrivateMarker(ilc[j],mrk);
					}
				}
				for(ElementArray<Node>::size_type it = 0; it < ret.size(); it++) ret[it]->RemPrivateMarker(mrk);
				ReleasePrivateMarker(mrk);
				break;
			}
		}

	}
	
	std::pair<Cell,bool> Mesh::CreateCell(const ElementArray<Face> & c_faces, const ElementArray<Node> & c_nodes)
	{
		HandleType he = InvalidHandle();
		if( !c_faces.empty() )
		{
			TopologyCheck chk = 0;
			chk |= BeginTopologyCheck(CELL,c_faces.data(),static_cast<enumerator>(c_faces.size()));
			if( chk != 0 )
			{
				if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
			}
			if( GetTopologyCheck(DUPLICATE_CELL) )
			{
				HandleType test = FindSharedAdjacency(c_faces.data(),static_cast<enumerator>(c_faces.size()));
				if (test != InvalidHandle()) return std::make_pair(Cell(this,test),false);
			}
			integer id = TieElement(3);
			he = ComposeHandleNum(3,id);
			for(ElementArray<Face>::size_type i = 0; i < c_faces.size(); i++)
			{
				Element::adj_type & hc = HighConn(c_faces.at(i));
#if defined(USE_OMP)
#pragma omp critical (face_high_conn)
#endif
				hc.push_back(he);
			}
			Element::adj_type & lc = LowConn(he);
			lc.insert(lc.begin(),c_faces.begin(),c_faces.end());
			ComputeGeometricType(he);		
			{
				//bool halt = false; //DEBUG
				if( c_nodes.empty() ) 
				{
					ElementArray<Node> nodes(this);
					RestoreCellNodes(he,nodes);
					for(ElementArray<Node>::size_type k = 0; k < nodes.size(); k++)
					{
						Element::adj_type & lc = LowConn(nodes.at(k));
#if defined(USE_OMP)
#pragma omp critical (node_low_conn)
#endif
						lc.push_back(he);
					}
					Element::adj_type & hc = HighConn(he);
					hc.insert(hc.begin(),nodes.begin(),nodes.end());

					//DEBUG
					/*
					for(ElementArray<Node>::size_type i = 0; i < nodes.size(); i++) 
					{
						if(!nodes[i]->CheckElementConnectivity())
						{
							nodes[i]->PrintElementConnectivity();
							halt = true;
						}
					}
					*/
					//END DEBUG
				}
				else
				{
					for(ElementArray<Node>::size_type k = 0; k < c_nodes.size(); k++)
					{
						Element::adj_type & lc = LowConn(c_nodes.at(k));
#if defined(USE_OMP)
#pragma omp critical (node_low_conn)
#endif
						lc.push_back(he);
					}
					Element::adj_type & hc = HighConn(he);
					hc.insert(hc.begin(),c_nodes.begin(),c_nodes.end());
					//DEBUG
					/*
					for(ElementArray<Face>::size_type i = 0; i < c_nodes.size(); i++) 
					{
						if(!c_nodes[i]->CheckElementConnectivity())
						{
							c_nodes[i]->PrintElementConnectivity();
							halt = true;
						}
					}
					*/
					//END DEBUG
				}
				//DEBUG
				/*
				if(!Element(this,he)->CheckElementConnectivity())
				{
					Element(this,he)->PrintElementConnectivity();
					halt = true;
				}
				for(ElementArray<Face>::size_type i = 0; i < c_faces.size(); i++)
				{
					if( !c_faces[i]->CheckElementConnectivity() )
					{
						c_faces[i]->PrintElementConnectivity();
						halt = true;
					}
				}
				assert(!halt);
				*/
				//END DEBUG
				SetMarker(he,NewMarker());
				RecomputeGeometricData(he);
				chk |= EndTopologyCheck(he);
				if( chk != 0 )
				{
					if( GetTopologyCheck(MARK_ON_ERROR)   ) Integer(he, TopologyErrorTag()) = chk;
					if( GetTopologyCheck(DELETE_ON_ERROR) ) { Destroy(he); he = InvalidHandle();}
					if( GetTopologyCheck(THROW_EXCEPTION) ) throw TopologyCheckError;
				}
			}

			
			
		}
		return std::make_pair(Cell(this,he),true);
	}
	
  std::pair<ElementSet,bool> Mesh::CreateSetUnique(std::string name)
	{
		HandleType he = ComposeHandleNum(4,TieElement(4));
		bulk_array set_name = BulkArrayDV(he,SetNameTag());
		set_name.resize(static_cast<bulk_array::size_type>(name.size()));
		memcpy(set_name.data(),name.c_str(),name.size());
		HighConn(he).resize(ElementSet::high_conn_reserved); //Allocate initial space for parent/sibling/child/unsorted info
		BulkDF(he, SetComparatorTag()) = ElementSet::UNSORTED_COMPARATOR;
		return std::make_pair(ElementSet(this,he),true);
	}
	
	std::pair<ElementSet,bool> Mesh::CreateSet(std::string name)
	{
		for(integer it = 0; it < EsetLastLocalID(); ++it) if( isValidElement(ESET,it) )
		{
			ElementSet e = EsetByLocalID(it);
			if( e->GetName() == name )
				return std::make_pair(e->self(),false);
		}
		return CreateSetUnique(name);
	}
	
	
	
	void Mesh::SetDimensions(integer dims)
	{
		if( dim == dims ) return;
		if( dims < 2 || dims > 3 ) throw DimensionIsNotSupportedByGeometry;

		if( NumberOfNodes() > 0)
		{
			std::vector<Storage::real> temp(dims*NodeLastLocalID());
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(Storage::integer k = 0; k < NodeLastLocalID(); ++k) if( isValidElementNum(0,k) )
				memcpy(temp.data()+k*dims,MGetDenseLink(ComposeHandleNum(0,k),CoordsTag()),sizeof(Storage::real)*dims);
			DeleteTag(tag_coords);
			tag_coords = CreateTag("COORD",DATA_REAL,NODE,NONE,dims);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(Storage::integer k = 0; k < NodeLastLocalID(); ++k) if( isValidElementNum(0,k) )
				memcpy(MGetDenseLink(ComposeHandleNum(0,k),CoordsTag()),temp.data()+k*dims,sizeof(Storage::real)*dims);
			dim = dims;
		}
		else
		{
			dim = dims;
			DeleteTag(tag_coords);
			tag_coords = CreateTag("PROTECTED_COORD",DATA_REAL, NODE,NONE,dim);
		}
	}
	
#if defined(USE_PARALLEL_WRITE_TIME)
#define REPORT_STR(x) {WriteTab(out_time) << "<TEXT><![CDATA[" << x << "]]></TEXT>" << std::endl;}
#define REPORT_VAL(str,x) {WriteTab(out_time) << "<VALUE name=\"" << str << "\"> <CONTENT><![CDATA[" << x << "]]></CONTENT> <CODE><![CDATA[" << #x << "]]></CODE></VALUE>" << std::endl;}
#endif

	void Mesh::UntieElement(integer etypenum, integer ID)
	{
#if defined(USE_OMP)
#pragma omp critical (links_interraction)
#endif
		{
			integer ADDR = links[etypenum][ID];
			links[etypenum][ID] = -1;
			back_links[etypenum][ADDR] = -1;
			empty_space[etypenum].push_back(ADDR);
			empty_links[etypenum].push_back(ID);
			//REPORT_VAL("destroyed",ComposeHandleNum(etypenum,ID) << " " << etypenum << " " << ADDR << " " << ID);
		}
	}
	



	Storage::integer Mesh::TieElement(integer etypenum)
	{
		integer ID = -1, ADDR = -1;
#if defined(USE_OMP)
#pragma omp critical (links_interraction)
#endif
		{
			INMOST_DATA_ENUM_TYPE old_size = GetArrayCapacity(etypenum), new_size;
			//ADDR points to gap in data
			if( !empty_space[etypenum].empty() && !isMeshModified() )
			{
				ADDR = empty_space[etypenum].back();
				empty_space[etypenum].pop_back();
			}
			else  //no gaps in data space, number of data entries is equal to total number of links
			{
				ADDR = static_cast<integer>(links[etypenum].size()-empty_links[etypenum].size());
			}
			// ID points to gap in links
			if( !empty_links[etypenum].empty() )
			{
				ID = empty_links[etypenum].back();
				empty_links[etypenum].pop_back();
				links[etypenum][ID] = ADDR;
			}
			else
			{
				ID = static_cast<integer>(links[etypenum].size());
				links[etypenum].push_back(ADDR);
			}
			new_size = GetArrayCapacity(etypenum);
			if( new_size != old_size ) ReallocateData(etypenum,new_size);
			back_links[etypenum][ADDR] = ID;
			last_created = ComposeHandleNum(etypenum,ID);
			//REPORT_VAL("created",last_created << " " << etypenum << " " << ADDR << " " << ID);
		}
		return ID;
	}
	
	
	void Mesh::MoveStorage(integer etypenum, integer old_addr, integer new_addr)
	{
		assert(old_addr != -1);
		assert(new_addr != -1);
		assert(old_addr < static_cast<integer>(back_links[etypenum].size()));
		assert(new_addr < static_cast<integer>(back_links[etypenum].size()));
		integer ID = back_links[etypenum][old_addr];
		assert(ID != -1);
		back_links[etypenum][old_addr] = -1;
		back_links[etypenum][new_addr] = ID;
		links[etypenum][ID] = new_addr;
#if !defined(LAZY_SPARSE_ALLOCATION)
		if( !sparse_data[etypenum].empty() )
#endif
		sparse_data[etypenum][new_addr].swap(sparse_data[etypenum][old_addr]);
		for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); t++) 
		{
			if( !t->isSparseByDim(etypenum) )
			{
				INMOST_DATA_ENUM_TYPE data_pos = t->GetPositionByDim(etypenum);
				if( data_pos == ENUMUNDEF ) continue;
				TagManager::dense_sub_type & arr = GetDenseData(data_pos);
				INMOST_DATA_ENUM_TYPE record_size = t->GetRecordSize();
				memcpy(&arr[new_addr],&arr[old_addr],record_size);
				memset(&arr[old_addr],0,record_size);
			}
		}
	}
	
	void Mesh::ReorderApply(Tag index, ElementType mask)
	{
		(void) index;
		(void) mask;
		throw NotImplemented;
	}
	
	void Mesh::ReorderEmpty(ElementType etype)
	{
		for(int etypenum = 0; etypenum < ElementNum(MESH); etypenum++) if( ElementTypeFromDim(etypenum) & etype )
		{
			integer cend = static_cast<integer>(back_links[etypenum].size())-1;
			//Very likely we have leading free space
			while( cend >= 0 && back_links[etypenum][cend] == -1 ) 
				cend--;
			while( cend >= 0 && !empty_space[etypenum].empty() )
			{
				integer last = empty_space[etypenum].back();
				if( last >= cend )
					empty_space[etypenum].pop_back();
				else if( back_links[etypenum][cend] != -1 )
				{
					MoveStorage(etypenum,cend,last);
					empty_space[etypenum].pop_back();
				}
				else cend--;
			}
			while( cend >= 0 && back_links[etypenum][cend] == -1 ) 
				cend--;
			//back_links[etypenum].resize(cend); //those should not be needed
			empty_space[etypenum].clear();
			ReallocateData(etypenum,GetArrayCapacity(etypenum));				
		}
	}

	MarkerType Mesh::CreateMarker()
	{
		Storage::bulk * marker_space = static_cast<Storage::bulk * >(MGetDenseLink(GetHandle(),tag_markers));
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
		assert(false); //if you reached here then you either don't release markers (it's your bug) or you should increase MarkerFields const in inmost_mesh.h
		return InvalidMarker();
	}
	MarkerType Mesh::CreatePrivateMarker()
	{
		int thread = GetLocalProcessorRank();
		Storage::bulk * marker_space = static_cast<Storage::bulk * >(MGetDenseLink(GetHandle(),tag_private_markers[thread]));
		INMOST_DATA_ENUM_TYPE ret;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < MarkerFields; ++k)
		{
			Storage::bulk mask = ((~marker_space[k]) & (-(~marker_space[k])));
			if( mask )
			{
				ret = (k << MarkerShift) | mask;
				marker_space[k] |= mask;
				ret |= MarkerPrivateBit;
				return ret;
			}
		}
		assert(false); //if you reached here then you either don't release markers (it's your bug) or you should increase MarkerFields const in inmost_mesh.h
		return InvalidMarker();
	}
  
	void Mesh::ReleaseMarker(MarkerType n, ElementType cleanup)
	{
		assert(!isPrivate(n));
		if( cleanup )
		{
			for(ElementType etype = NODE; etype < MESH; etype = NextElementType(etype)) if( cleanup & etype )
			{
				integer end = LastLocalID(etype);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer id = 0; id < end; ++id) 
					if( isValidElement(etype,id) )
						RemMarker(ComposeHandle(etype,id),n);
			}
		}
#if defined(CHECKS_MARKERS)
#ifndef NDEBUG
		if( check_shared_mrk )
		{
			for(int etypenum = 0; etypenum < ElementNum(MESH); ++etypenum)
			{
				integer end = LastLocalIDNum(etypenum);
				for(integer id = 0; id < end; ++id)
					if( isValidElementNum(etypenum,id) )
						assert((static_cast<const bulk *>(MGetDenseLink(etypenum,id,MarkersTag()))[n >> MarkerShift] & static_cast<bulk>(n & MarkerMask)) == 0 && "marker was not properly cleared from elements");
			}
		}
#endif
#endif
		Storage::RemMarker(n);
	}

	void Mesh::ReleasePrivateMarker(MarkerType n, ElementType cleanup)
	{
		assert(isPrivate(n));
		if( cleanup )
		{
			for(ElementType etype = NODE; etype < MESH; etype = NextElementType(etype)) if( cleanup & etype )
			{
				integer end = LastLocalID(etype);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer id = 0; id < end; ++id)
					if( isValidElement(etype,id) )
						RemPrivateMarker(ComposeHandle(etype,id),n);
			}
		}
#if defined(CHECKS_MARKERS)
#ifndef NDEBUG
		if( check_private_mrk )
		{
			int thread = GetLocalProcessorRank();
			MarkerType pn = n & ~MarkerPrivateBit;
			for(int etypenum = 0; etypenum < ElementNum(MESH); ++etypenum)
			{
				integer end = LastLocalIDNum(etypenum);
				for(integer id = 0; id < end; ++id)
					if( isValidElementNum(etypenum,id) )
						assert((static_cast<const bulk *>(MGetDenseLink(etypenum,id,tag_private_markers[thread]))[pn >> MarkerShift] & static_cast<bulk>(pn & MarkerMask)) == 0 && "marker was not properly cleared from elements");
			}
		}
#endif
#endif
		Storage::RemPrivateMarker(n);
	}
	
	//this follows chunk_array definition
	INMOST_DATA_ENUM_TYPE Mesh::GetArrayCapacity(integer etypenum)
	{
		INMOST_DATA_ENUM_TYPE occupied = static_cast<INMOST_DATA_ENUM_TYPE>(links[etypenum].size() - empty_links[etypenum].size() + empty_space[etypenum].size());
		INMOST_DATA_ENUM_TYPE chunks = (occupied >> chunk_bits_elems) + ((occupied & ((1 << chunk_bits_elems)-1))?1:0);
		return chunks * (1 << chunk_bits_elems);
	}
#if defined(USE_AUTODIFF)
  Storage::var & Mesh::Variable(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_VARIABLE);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() != ENUMUNDEF ) 
			return static_cast<var*>(p)[0]; 
		else 
			return static_cast<inner_variable_array*>(p)->at_safe(0);
	}
  Storage::var_array Mesh::VariableArray(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_VARIABLE);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() == ENUMUNDEF ) 
			return var_array(*static_cast<inner_variable_array*>(p)); 
		else 
			return var_array(static_cast<var*>(p),tag.GetSize());
	}
#endif
	Storage::real & Mesh::Real(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_REAL);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() != ENUMUNDEF ) 
			return static_cast<real*>(p)[0]; 
		else 
			return static_cast<inner_real_array*>(p)->at_safe(0);
	}
	Storage::integer & Mesh::Integer(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_INTEGER);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() != ENUMUNDEF ) 
			return static_cast<integer*>(p)[0]; 
		else 
			return static_cast<inner_integer_array*>(p)->at_safe(0);}
	Storage::bulk & Mesh::Bulk(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_BULK);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() != ENUMUNDEF ) 
			return static_cast<bulk*>(p)[0]; 
		else 
			return static_cast<inner_bulk_array*>(p)->at_safe(0);
	}
	Storage::reference & Mesh::Reference(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_REFERENCE);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() != ENUMUNDEF ) 
			return static_cast<reference*>(p)[0]; 
		else 
			return static_cast<inner_reference_array*>(p)->at_safe(0);
	}
  Storage::remote_reference & Mesh::RemoteReference(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_REMOTE_REFERENCE);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() != ENUMUNDEF ) 
			return static_cast<remote_reference*>(p)[0]; 
		else 
			return static_cast<inner_remote_reference_array*>(p)->at_safe(0);
	}
	Storage::real_array Mesh::RealArray(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_REAL);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() == ENUMUNDEF ) 
			return real_array(*static_cast<inner_real_array*>(p)); 
		else 
			return real_array(static_cast<real*>(p),tag.GetSize());
	}
  
	Storage::integer_array  Mesh::IntegerArray(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_INTEGER);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() == ENUMUNDEF ) 
			return integer_array(*static_cast<inner_integer_array*>(p)); 
		else 
			return integer_array(static_cast<integer*>(p),tag.GetSize());
	}
	Storage::bulk_array Mesh::BulkArray(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_BULK);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() == ENUMUNDEF ) 
			return bulk_array(*static_cast<inner_bulk_array     *>(p)); 
		else 
			return bulk_array(static_cast<bulk*>(p),tag.GetSize());
  }
	Storage::reference_array Mesh::ReferenceArray(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_REFERENCE);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() == ENUMUNDEF ) 
			return reference_array(this,*static_cast<inner_reference_array*>(p)); 
		else 
			return reference_array(this,static_cast<reference *>(p),tag.GetSize());
	}
  Storage::remote_reference_array Mesh::RemoteReferenceArray(HandleType h, const Tag & tag) 
	{
		Asserts(h,tag,DATA_REMOTE_REFERENCE);
		void * p = MGetLink(h,tag); 
		if( tag.GetSize() == ENUMUNDEF ) 
			return remote_reference_array(*static_cast<inner_remote_reference_array*>(p)); 
		else 
			return remote_reference_array(static_cast<remote_reference *>(p),tag.GetSize());
	}
	
	void Mesh::AssertsDV(HandleType h, const Tag & tag, DataType expected) const
	{
		Asserts(h,tag,expected);
		assert(tag.GetSize() == ENUMUNDEF);              //data is of variable size
		assert(!tag.isSparseByDim(GetHandleElementNum(h)));//tag is not sparse
	}

	void Mesh::AssertsDF(HandleType h, const Tag & tag, DataType expected) const
	{
		Asserts(h,tag,expected);
		assert(tag.GetSize() != ENUMUNDEF);              //data is of variable size
		assert(!tag.isSparseByDim(GetHandleElementNum(h)));//tag is not sparse
	}

	void Mesh::Asserts(HandleType h, const Tag & tag, DataType expected) const
	{
		assert(isValidHandleRange(h));                   //handle doesn't point out of range
		assert(isValidHandle(h));                        //handle doesn't point to deleted element
		assert(tag.isValid());                           //tag was allocated
		assert(this == tag.GetMeshLink());               //tag is not mine
		assert(tag.GetDataType() == expected);           //tag data type coinside with expected data type
		assert(tag.isDefinedByDim(GetHandleElementNum(h)));           //tag data type coinside with expected data type
	}
	
	void Mesh::ClearMarkerSpace(HandleType h) 
	{
		Storage::bulk * marker_space = static_cast<Storage::bulk *>(MGetDenseLink(h,MarkersTag()));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < MarkerFields; ++k) marker_space[k] = 0;
	}
	void Mesh::GetMarkerSpace(HandleType h, Storage::bulk copy[MarkerFields]) const 
	{
		const Storage::bulk * marker_space = static_cast<const Storage::bulk *>(MGetDenseLink(h,MarkersTag()));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < MarkerFields; ++k) copy[k] = marker_space[k];
	}
	void Mesh::SetMarkerSpace(HandleType h, Storage::bulk source[MarkerFields]) 
	{
		Storage::bulk * marker_space = static_cast<Storage::bulk *>(MGetDenseLink(h,MarkersTag()));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < MarkerFields; ++k) marker_space[k] = source[k];
	}

	INMOST_DATA_ENUM_TYPE Mesh::GetDataSize(HandleType h,const Tag & tag) const
	{
		assert( tag.GetMeshLink() == this );
		if( tag.GetSize() == ENUMUNDEF )
		{
			const void * adata = MGetLink(h,tag);
			assert( adata != NULL );
			switch(tag.GetDataType())
			{
				case DATA_REAL:     return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_real_array     *>(adata)->size());
				case DATA_INTEGER:  return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_integer_array  *>(adata)->size());
				case DATA_BULK:     return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_bulk_array     *>(adata)->size());
				case DATA_REFERENCE:return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_reference_array*>(adata)->size());
        case DATA_REMOTE_REFERENCE:
                            return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_remote_reference_array*>(adata)->size());
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE: return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_variable_array *>(adata)->size());
#endif
			}
			throw BadTag;
		}
		return tag.GetSize();
	}
  INMOST_DATA_ENUM_TYPE Mesh::GetDataCapacity(const INMOST_DATA_BULK_TYPE * adata, INMOST_DATA_ENUM_TYPE size, const Tag & tag) const
  {
    assert( tag.GetMeshLink() == this );
#if defined(USE_AUTODIFF)
    if( tag.GetDataType() == DATA_VARIABLE )
    {
      INMOST_DATA_ENUM_TYPE ret = 0;
      const Sparse::Row::entry * arr = static_cast<const Sparse::Row::entry *>(static_cast<const void *>(adata));
      for(INMOST_DATA_ENUM_TYPE k = 0; k < size; ++k)
        ret += variable::RetriveSize(arr+ret);
      return ret*sizeof(Sparse::Row::entry);
    }
    else
#endif
      return size*tag.GetBytesSize();
    assert(false);
    return 0;
  }
  INMOST_DATA_ENUM_TYPE Mesh::GetDataCapacity(HandleType h,const Tag & tag) const
	{
		assert( tag.GetMeshLink() == this );
		if( tag.GetSize() == ENUMUNDEF )
		{
			const void * adata = MGetLink(h,tag);
			assert( adata != NULL );
			switch(tag.GetDataType())
			{
        case DATA_REAL:     return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_real_array     *>(adata)->size())*tag.GetBytesSize();
				case DATA_INTEGER:  return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_integer_array  *>(adata)->size())*tag.GetBytesSize();
				case DATA_BULK:     return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_bulk_array     *>(adata)->size())*tag.GetBytesSize();
				case DATA_REFERENCE:return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_reference_array*>(adata)->size())*tag.GetBytesSize();
        case DATA_REMOTE_REFERENCE:
                            return static_cast<INMOST_DATA_ENUM_TYPE>(static_cast<const inner_remote_reference_array*>(adata)->size())*tag.GetBytesSize();
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE: 
          {
            INMOST_DATA_ENUM_TYPE ret = 0;
            const inner_variable_array * arr = static_cast<const inner_variable_array *>(adata);
            for(inner_variable_array::size_type k = 0; k < arr->size(); ++k)
              ret += (*arr)[k].RecordSize();//(*arr)[k].GetRow().Size();
            return ret*sizeof(Sparse::Row::entry_s);
          }
#endif
			}
			throw BadTag;
		}
#if defined(USE_AUTODIFF)
    if( tag.GetDataType() == DATA_VARIABLE )
      return static_cast<const var *>(MGetLink(h,tag))->RecordSize()*sizeof(Sparse::Row::entry_s);
#endif
		return tag.GetSize()*tag.GetBytesSize();
	}
	void Mesh::SetDataSize(HandleType h,const Tag & tag,INMOST_DATA_ENUM_TYPE new_size)
	{
		assert( tag.GetMeshLink() == this );
		void * adata = MGetLink(h,tag);
		assert( adata != NULL );
		if( tag.GetSize() == ENUMUNDEF )
		{
			switch(tag.GetDataType())
			{
				case DATA_REAL:     static_cast<inner_real_array      *>(adata)->resize(new_size); break;
				case DATA_INTEGER:  static_cast<inner_integer_array   *>(adata)->resize(new_size); break;
				case DATA_BULK:     static_cast<inner_bulk_array      *>(adata)->resize(new_size); break;
				case DATA_REFERENCE:static_cast<inner_reference_array *>(adata)->resize(new_size); break;
        case DATA_REMOTE_REFERENCE:
                            static_cast<inner_remote_reference_array *>(adata)->resize(new_size); break;
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE: static_cast<inner_variable_array  *>(adata)->resize(new_size); break;
#endif
			}
			return;
		}
		else if( tag.GetSize() == new_size )
			return;
		assert(false);
	}
	void Mesh::GetData(HandleType h,const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, void * data_out) const
	{
		assert( tag.GetMeshLink() == this );
		const void * adata = MGetLink(h,tag);
		assert( adata != NULL );
		INMOST_DATA_ENUM_TYPE data_size = tag.GetSize();
		INMOST_DATA_ENUM_TYPE bytes = tag.GetBytesSize();
		if( data_size == ENUMUNDEF )
		{
			switch(tag.GetDataType())
			{
				case DATA_REAL:     memcpy(data_out,&(*static_cast<const inner_real_array      *>(adata))[shift],bytes*size); break;
				case DATA_INTEGER:  memcpy(data_out,&(*static_cast<const inner_integer_array   *>(adata))[shift],bytes*size); break;
				case DATA_BULK:     memcpy(data_out,&(*static_cast<const inner_bulk_array      *>(adata))[shift],bytes*size); break;
				case DATA_REFERENCE:memcpy(data_out,&(*static_cast<const inner_reference_array *>(adata))[shift],bytes*size); break;
        case DATA_REMOTE_REFERENCE:
                            memcpy(data_out,&(*static_cast<const inner_remote_reference_array *>(adata))[shift],bytes*size); break;
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE:
          {
            const inner_variable_array * arr = static_cast<const inner_variable_array      *>(adata);
            Sparse::Row::entry_s * data = static_cast<Sparse::Row::entry_s *>(data_out);
            int k = 0;
            for(INMOST_DATA_ENUM_TYPE r = 0; r < size; ++r)
              k += (*arr)[r+shift].Record(data+k);
          }
          break;
#endif
			}
		}
#if defined(USE_AUTODIFF)
    else if( tag.GetDataType() == DATA_VARIABLE )
    {
      Sparse::Row::entry_s * data = static_cast<Sparse::Row::entry_s *>(data_out);
      const var * v = static_cast<const var *>(MGetLink(h,tag));
      int k = 0;
      for(INMOST_DATA_ENUM_TYPE r = 0; r < size; ++r)
        k += v[r+shift].Record(data+k);
    }
#endif
		else memcpy(data_out,static_cast<const INMOST_DATA_BULK_TYPE *>(adata)+shift*bytes,size*bytes);
		return;
	}
	void Mesh::SetData(HandleType h, const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, const void * data_in)
	{
		assert( tag.GetMeshLink() == this );
		void * adata = MGetLink(h,tag);
		assert( adata != NULL );
		INMOST_DATA_ENUM_TYPE data_size = tag.GetSize();
		INMOST_DATA_ENUM_TYPE bytes = tag.GetBytesSize();
		if( data_size == ENUMUNDEF )
		{
			switch(tag.GetDataType())
			{
				case DATA_REAL:     memcpy(&(*static_cast<inner_real_array     *>(adata))[shift],data_in,bytes*size); break;
				case DATA_INTEGER:  memcpy(&(*static_cast<inner_integer_array  *>(adata))[shift],data_in,bytes*size); break;
				case DATA_BULK:     memcpy(&(*static_cast<inner_bulk_array     *>(adata))[shift],data_in,bytes*size); break;
				case DATA_REFERENCE:memcpy(&(*static_cast<inner_reference_array*>(adata))[shift],data_in,bytes*size); break;
        case DATA_REMOTE_REFERENCE:
                            memcpy(&(*static_cast<inner_reference_array*>(adata))[shift],data_in,bytes*size); break;
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE:
          {
            inner_variable_array * arr = static_cast<inner_variable_array*>(adata);
            const Sparse::Row::entry_s * data = static_cast<const Sparse::Row::entry_s *>(data_in);
            int k = 0;
            for(INMOST_DATA_ENUM_TYPE r = 0; r < size; ++r)
              k += (*arr)[r+shift].Retrive(data+k);
          }
          break;
#endif
			}
		}
#if defined(USE_AUTODIFF)
    else if( tag.GetDataType() == DATA_VARIABLE )
    {
      const Sparse::Row::entry_s * data = static_cast<const Sparse::Row::entry_s *>(data_in);
      var * v = static_cast<var *>(MGetLink(h,tag));
      int k = 0;
      for(INMOST_DATA_ENUM_TYPE r = 0; r < size; ++r)
        k += v[r+shift].Retrive(data+k);
    }
#endif
		else memcpy(static_cast<INMOST_DATA_BULK_TYPE *>(adata)+shift*bytes,data_in,size*bytes);
	}


	void Mesh::DelDenseData(HandleType h, const Tag & tag)
	{
		assert( tag.GetMeshLink() == this );
		assert( !tag.isSparseByDim(GetHandleElementNum(h)) );
		void * data = MGetLink(h,tag);
		if( data != NULL )
		{
			if( tag.GetSize() == ENUMUNDEF )
				TagManager::DestroyVariableData(tag,data);
#if defined(USE_AUTODIFF)
      else if( tag.GetDataType() == DATA_VARIABLE ) //Have to deallocate the structure to remove inheritance
      {
        for(INMOST_DATA_ENUM_TYPE k = 0; k < tag.GetSize(); ++k)
          (static_cast<variable *>(data)[k]).~variable();
      }
#endif
			//else if( tag.GetDataType() == DATA_REFERENCE )
			//	 memset(data,0xff,tag.GetRecordSize());
			else memset(data,0,tag.GetRecordSize());
		}
	}

	void Mesh::DelSparseData(HandleType h,const Tag & tag)
	{
		assert( tag.GetMeshLink() == this );
		assert( tag.isSparseByDim(GetHandleElementNum(h)) );
		sparse_type & s = MGetSparseLink(h);
		for(sparse_type::size_type i = 0; i < s.size(); ++i) if( s[i].tag == tag.mem )
		{
			if( tag.GetSize() == ENUMUNDEF ) 
				TagManager::DestroyVariableData(tag,s[i].rec);
#if defined(USE_AUTODIFF)
      else if( tag.GetDataType() == DATA_VARIABLE ) //Have to deallocate the structure to remove inheritance
      {
        for(INMOST_DATA_ENUM_TYPE k = 0; k < tag.GetSize(); ++k)
          (static_cast<variable *>(s[i].rec)[k]).~variable();
      }
#endif
			free(s[i].rec);
			s.erase(s.begin()+i);
			break;
		}
	}
	
	void Mesh::DelData(HandleType h,const Tag & tag)
	{
		assert( tag.GetMeshLink() == this );
		if( tag.isSparse(GetHandleElementType(h)) ) 
			DelSparseData(h,tag);
		else 
			DelDenseData(h,tag);
	}

	bool  Mesh::HaveData(HandleType h,const Tag & tag) const 
	{
		integer n = GetHandleElementNum(h);
		if(tag.isSparseByDim(n)) 
		{ 
			if( MGetSparseLink(h,tag) != NULL ) 
				return true; 
			return false; 
		} 
		else 
		{
			if( tag.GetPositionByDim(n) != ENUMUNDEF ) 
				return true; 
			return false;
		}
	}

	bool Mesh::isValidHandleRange (HandleType h) const
	{
		if( !isValidHandle(h) ) return false;
		if( GetHandleElementNum(h) >= 0 && GetHandleElementNum(h) < 5 )
			return GetHandleID(h) < static_cast<integer>(links[GetHandleElementNum(h)].size());
		else if( GetHandleElementNum(h) == 5 )
			return GetHandleID(h) == 0;
		else
			return false;
	}


	ElementSet Mesh::GetSet(std::string name)
	{
		for(integer it = 0; it < EsetLastLocalID(); ++it) if( isValidElement(ESET,it) )
		{
			ElementSet e = EsetByLocalID(it);
			if( e->GetName() == name )
				return e;
		}
		return ElementSet(this,InvalidHandle());
	}

	ElementArray<ElementSet> Mesh::GetSetsByPrefix(std::string name)
	{
		ElementArray<ElementSet> ret(this);
		for(integer it = 0; it < EsetLastLocalID(); ++it) if( isValidElement(ESET,it) )
		{
			ElementSet e = EsetByLocalID(it);
			if( e->GetName().substr(0,name.size()) == name )
				ret.push_back(e->self());
		}
		return ret;
	}

	void Mesh::SetTopologyCheck   (TopologyCheck mask) 
	{
		checkset = checkset | mask;
		if( mask & MARK_ON_ERROR ) 
			tag_topologyerror = CreateTag("TOPOLOGY_ERROR_TAG",DATA_INTEGER,CELL | FACE | EDGE,CELL | FACE | EDGE,1);
	}
	void Mesh::RemTopologyCheck   (TopologyCheck mask) 
	{
		checkset = checkset & ~mask;
		if( mask & MARK_ON_ERROR ) 
			tag_topologyerror = DeleteTag(tag_topologyerror);
	}

  Mesh::Random::Random(unsigned int seed)
	{
		n = seed;
		a = 1103515245;
		c = 12345;
		m = 1u << (sizeof(unsigned int)*8-1);
	}
	Mesh::Random::Random(const Random & other)
	{
		n = other.n;
		a = other.a;
		c = other.c;
		m = other.m;
	}
	unsigned int Mesh::Random::Number()
	{
		n = (a*n + c)%m;
		return (n << 2) >> 18;
	}
}
#endif
