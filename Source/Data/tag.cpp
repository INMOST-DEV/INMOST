#include "inmost.h"
#if defined(USE_MESH)

namespace INMOST
{
	const char * DataTypeName(DataType t)
	{
		switch(t)
		{
			case DATA_REAL:      return "REAL";
			case DATA_INTEGER:   return "INTEGER";
			case DATA_BULK:      return "BULK";
			case DATA_REFERENCE: return "REFERENCE";
#if defined(USE_AUTODIFF)
      case DATA_VARIABLE:  return "VARIABLE";
#endif
		}
		return "UNKNOWN";
	}


	__INLINE static INMOST_DATA_ENUM_TYPE DataTypeBytesSize(DataType t)
	{
		switch(t)
		{
			case DATA_BULK:      return sizeof(INMOST_DATA_BULK_TYPE);
			case DATA_INTEGER:   return sizeof(INMOST_DATA_INTEGER_TYPE);
			case DATA_REAL:      return sizeof(INMOST_DATA_REAL_TYPE);
			case DATA_REFERENCE: return sizeof(HandleType);
#if defined(USE_AUTODIFF)
      case DATA_VARIABLE:  return sizeof(variable);
#endif
		}
		return 0;
	}


	__INLINE static INMOST_DATA_ENUM_TYPE VariableDataSize(DataType t)
	{
		switch(t)
		{
			case DATA_REAL:      return sizeof(inner_real_array);
			case DATA_INTEGER:   return sizeof(inner_integer_array);
			case DATA_BULK:      return sizeof(inner_bulk_array);
			case DATA_REFERENCE: return sizeof(inner_reference_array);
#if defined(USE_AUTODIFF)
      case DATA_VARIABLE:  return sizeof(inner_variable_array);
#endif
		}
		return 0;
	}


  Tag::~Tag() 
  {
    mem = NULL;
  }

  Tag::Tag() 
  {
    mem = NULL;
  }

  Tag::Tag(const Tag & other) 
  {
    mem = other.mem;
  }

	void TagManager::CopyData(const Tag & t, void * adata, const void * bdata)
	{
		
		INMOST_DATA_ENUM_TYPE data_size = t.GetSize();
		INMOST_DATA_ENUM_TYPE bytes = t.GetBytesSize();
    DataType type = t.GetDataType();
		if( data_size == ENUMUNDEF ) //variable size array
		{
			if( adata != NULL ) TagManager::DestroyVariableData(t,adata);
			if( type == DATA_REAL )           new (adata) inner_real_array     (*static_cast<const inner_real_array      * >(bdata));
			else if( type == DATA_INTEGER )   new (adata) inner_integer_array  (*static_cast<const inner_integer_array   * >(bdata));
			else if( type == DATA_BULK )      new (adata) inner_bulk_array     (*static_cast<const inner_bulk_array      * >(bdata));
			else if( type == DATA_REFERENCE ) new (adata) inner_reference_array(*static_cast<const inner_reference_array * >(bdata));
#if defined(USE_AUTODIFF)
      else if( type == DATA_VARIABLE )  new (adata) inner_variable_array (*static_cast<const inner_variable_array  * >(bdata));
#endif
		}
#if defined(USE_AUTODIFF)
    else if( type == DATA_VARIABLE ) //have to call constructor
    {
      for(INMOST_DATA_ENUM_TYPE k = 0; k < data_size; ++k)
        new (static_cast<variable *>(adata)+k) variable(*(static_cast<const variable *>(bdata)+k));
    }
#endif
		else // fixed size array
			memcpy(adata,bdata,data_size*bytes);
	}
	void TagManager::DestroyVariableData(const Tag & t, void * adata)
	{
		INMOST_DATA_ENUM_TYPE data_size = t.GetSize();
		if( data_size == ENUMUNDEF && adata != NULL ) //variable size array
		{
			DataType type = t.GetDataType();
			if( type == DATA_REAL ) 
			{
				(*static_cast<inner_real_array *> (adata)).~inner_real_array();
				new (adata) inner_real_array(); //reinitialize for reuse
			}
			else if( type == DATA_INTEGER ) 
			{
				(*static_cast<inner_integer_array *> (adata)).~inner_integer_array();
				new (adata) inner_integer_array();
			}
			else if( type == DATA_BULK ) 
			{
				(*static_cast<inner_bulk_array *> (adata)).~inner_bulk_array();
				new (adata) inner_bulk_array();
			}
			else if( type == DATA_REFERENCE ) 
			{
				(*static_cast<inner_reference_array *> (adata)).~inner_reference_array();
				new (adata) inner_reference_array();
			}
#if defined(USE_AUTODIFF)
      else if( type == DATA_VARIABLE ) 
      {
        (*static_cast<inner_variable_array *> (adata)).~inner_variable_array();
        new (adata) inner_variable_array();
      }
#endif
		}
	}
	
	TagMemory::TagMemory(Mesh * m, const TagMemory & other)
	{
		for(int i = 0; i < NUM_ELEMENT_TYPS; i++)
		{
			pos[i]	= other.pos[i];
			sparse[i] = other.sparse[i];
		}
		tagname		   = other.tagname;
		dtype		   = other.dtype;
		bulk_data_type = other.bulk_data_type;
		m_link         = m;
		size           = other.size;
		record_size    = other.record_size;
		bytes_size     = other.bytes_size;
	}
	
	TagMemory & TagMemory::operator =(TagMemory const & other)
	{
		for(int i = 0; i < NUM_ELEMENT_TYPS; i++)
		{
			pos[i]      = other.pos[i];
			sparse[i]   = other.sparse[i];
		}
		tagname         = other.tagname;
		dtype           = other.dtype;
		bulk_data_type  = other.bulk_data_type;
		m_link          = other.m_link;
		size            = other.size;
		record_size     = other.record_size;
		bytes_size      = other.bytes_size;
		return *this;	
	}
	
	TagMemory::TagMemory()
	{
		for(int i = 0; i < NUM_ELEMENT_TYPS; i++)
		{
			pos[i]	= ENUMUNDEF;
			sparse[i] = false;
		}
		tagname = "";
	}
	
	Tag::Tag(Mesh * m, std::string name, DataType _dtype,INMOST_DATA_ENUM_TYPE size)
	{
		mem = new TagMemory();
		for(int i = 0; i < NUM_ELEMENT_TYPS; i++)
		{
			mem->pos[i]	= ENUMUNDEF;
			mem->sparse[i] = false;
		}
		mem->tagname	= name;
		mem->dtype		= _dtype;
		mem->size		= size;
		switch(mem->dtype)
		{
			case DATA_BULK:      mem->bulk_data_type = INMOST_MPI_DATA_BULK_TYPE;    break;
			case DATA_REAL:      mem->bulk_data_type = INMOST_MPI_DATA_REAL_TYPE;    break;
			case DATA_INTEGER:   mem->bulk_data_type = INMOST_MPI_DATA_INTEGER_TYPE; break;
			case DATA_REFERENCE: mem->bulk_data_type = INMOST_MPI_DATA_ENUM_TYPE;    break;
#if defined(USE_AUTODIFF)
      case DATA_VARIABLE:
        if( !Sparse::HaveRowEntryType() ) Sparse::CreateRowEntryType();
        mem->bulk_data_type = Sparse::GetRowEntryType();
        break; 
#endif
		}
		mem->bytes_size = DataTypeBytesSize(mem->dtype);
		if(mem->size == ENUMUNDEF )
			mem->record_size = VariableDataSize(mem->dtype);
		else
			mem->record_size = mem->size * mem->bytes_size;
		mem->m_link = m;
	}
	TagManager::TagManager()
	{
	}
	TagManager::TagManager(const TagManager & other)
	{
		tags.resize(other.tags.size());
		dense_data.resize(other.dense_data.size(),dense_sub_type(0));
		for(tag_array_type::size_type i = 0; i < other.tags.size(); i++)
		{
			tags[i].mem = new TagMemory(dynamic_cast<Mesh *>(this),*other.tags[i].mem);
			for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if( tags[i].isDefined(etype) && !tags[i].isSparse(etype) )
				{
					dense_data[tags[i].GetPosition(etype)] = dense_sub_type(tags[i].GetRecordSize());
					//tags[i].AllocateData(etype);
				}
		}
		
		//The rest of the data should be allocated and copied while copying data of
		//individual elements
	}

	TagManager & TagManager::operator =(TagManager const & other)
	{

		//for(tag_iterator it = tags.begin(); it != tags.end(); it++)
		//{
		//	for(ElementType etype = NODE; etype <= MESH; etype = etype << 1 )
		//		if( it->isDefined(etype) && !it->isSparse(etype) && it->GetSize() == ENUMUNDEF )
		//		{
		//			INMOST_DATA_ENUM_TYPE record_size = it->GetRecordSize();
		//			TagManager::dense_sub_type & arr = dense_data[it->GetPosition(etype)];
		//			TagManager::dense_sub_type::iterator jt = arr.begin();
		//			switch(it->GetDataType())
		//			{
		//				case DATA_REAL:      while( jt != arr.end() ) { void * p = static_cast<void *>(&*jt); if( p != NULL ) (*static_cast<inner_real_array      *>( p )).~inner_real_array();      jt+=record_size; }  break;
		//				case DATA_INTEGER:   while( jt != arr.end() ) { void * p = static_cast<void *>(&*jt); if( p != NULL ) (*static_cast<inner_integer_array   *>( p )).~inner_integer_array();   jt+=record_size; }  break;
		//				case DATA_BULK:      while( jt != arr.end() ) { void * p = static_cast<void *>(&*jt); if( p != NULL ) (*static_cast<inner_bulk_array      *>( p )).~inner_bulk_array();      jt+=record_size; }  break;
		//				case DATA_REFERENCE: while( jt != arr.end() ) { void * p = static_cast<void *>(&*jt); if( p != NULL ) (*static_cast<inner_reference_array *>( p )).~inner_reference_array(); jt+=record_size; }  break;
		//			}
		//		}
		//	delete it->mem;
		//}
		tags.resize(other.tags.size());
		dense_data.clear();
		dense_data.resize(other.dense_data.size(),dense_sub_type(0));
		for(tag_array_type::size_type i = 0; i < other.tags.size(); i++)
		{
			tags[i].mem = new TagMemory(dynamic_cast<Mesh *>(this),*other.tags[i].mem);
			for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if( tags[i].isDefined(etype) && !tags[i].isSparse(etype) )
				{
					dense_data[tags[i].GetPosition(etype)] = dense_sub_type(tags[i].GetRecordSize());
					//tags[i].AllocateData(etype);
				}
		}
		for(int i = 0; i < 6; i++)
			sparse_data[i].clear();
		//~ sparse_data.resize(other.sparse_data.size());
		//The rest of the data should be allocated and copied while copying data of
		//individual elements
		return *this;
	}
	
	TagManager::~TagManager()
	{
		dense_data.clear();
		for(int i = 0; i < 6; i++) sparse_data[i].clear();
		for(tag_iterator it = tags.begin(); it != tags.end(); it++) delete it->mem;
		tags.clear();
	}
	
	
	
	Tag TagManager::CreateTag(Mesh *m, std::string name, DataType dtype, ElementType etype,ElementType sparse, INMOST_DATA_ENUM_TYPE size)
	{
		Tag new_tag;
#if !defined(LAZY_SPARSE_ALLOCATION)
		bool need_sparse[6] = {false,false,false,false,false,false};
#endif
		for(tag_array_type::size_type i = 0; i < tags.size(); i++)
    {
			if( tags[i].GetTagName() == name )
			{
				assert( tags[i].GetDataType() == dtype && (size == ENUMUNDEF || size == tags[i].GetSize()) );
				//if( tags[i].GetDataType() != dtype || (size != ENUMUNDEF && size != tags[i].GetSize()) )
				//{
				//	throw TagExists;
				//}
				new_tag = tags[i];
				break;
			}
    }
		if( !new_tag.isValid() )
		{
			new_tag = Tag(m,name,dtype,size);
#if defined(USE_OMP)
#pragma omp critical
#endif
      {
			  tags.push_back(new_tag);
      }
		}
		for(ElementType mask = NODE; mask <= MESH; mask = mask << 1)
		{
			if( (mask & etype) && new_tag.GetPosition(etype&mask) == ENUMUNDEF ) 
			{
				if( sparse & mask )
				{
					new_tag.SetPosition(ENUMUNDEF-1,mask);
					new_tag.SetSparse(mask);
#if !defined(LAZY_SPARSE_ALLOCATION)
					need_sparse[ElementNum(mask)] = true;
#endif
				}
				else
				{
          INMOST_DATA_ENUM_TYPE new_pos = ENUMUNDEF;
#if defined(USE_OMP)
#pragma omp critical
#endif
          {
					  new_pos = static_cast<INMOST_DATA_ENUM_TYPE>(dense_data.size());
					  if( !empty_dense_data.empty() )
					  {
						  new_pos = empty_dense_data.back();
						  empty_dense_data.pop_back();
						  dense_data[new_pos] = dense_sub_type(new_tag.GetRecordSize());
					  }
					  else dense_data.push_back(dense_sub_type(new_tag.GetRecordSize()));
          }
					new_tag.SetPosition(new_pos,mask);
					INMOST_DATA_ENUM_TYPE new_size = dynamic_cast<Mesh *>(this)->GetArrayCapacity(ElementNum(mask));
					if( new_size < 1024 && mask != MESH ) new_size = 1024;
					if( new_size != 1   && mask == MESH ) new_size = 1;
					ReallocateData(new_tag,ElementNum(mask),new_size);
				}
			}
		}
#if !defined(LAZY_SPARSE_ALLOCATION)
		for(int j = 0; j < 6; j++) 
			if( need_sparse[j] && sparse_data[j].empty() )
			{
				INMOST_DATA_ENUM_TYPE new_size = dynamic_cast<Mesh *>(this)->GetArrayCapacity(j);
				if( new_size < 1024 && j != ElementNum(MESH) ) new_size = 1024;
				if( new_size != 1   && j == ElementNum(MESH) ) new_size = 1;
#if defined(USE_OMP)
#pragma omp critical
#endif
        {
				  sparse_data[j].resize(new_size);
        }
			}
#endif
		return new_tag;
	}
	Tag TagManager::GetTag(std::string name) const
	{
		for(tag_array_type::size_type i = 0; i < tags.size(); i++)
			if( tags[i].GetTagName() == name )
				return tags[i];
		assert(false);
		return Tag();
	}
	bool TagManager::HaveTag(std::string name) const
	{
		for(tag_array_type::size_type i = 0; i < tags.size(); i++)
			if( tags[i].GetTagName() == name )
				return true;
		return false;
	}
	Tag TagManager::DeleteTag(Tag tag, ElementType type_mask)
	{
		bool delete_entirely = true;
#if !defined(LAZY_SPARSE_ALLOCATION)
		bool was_sparse[6] = {false,false,false,false,false,false};
#endif
		INMOST_DATA_ENUM_TYPE tpos;//,ipos;
		for(ElementType mask = NODE; mask <= MESH; mask = mask << 1 )
		{
			tpos = tag.GetPosition(mask);
			if( tpos == ENUMUNDEF ) continue;
			if( mask & type_mask )
			{
				if( !tag.isSparse(mask) ) 
				{
          dense_data[tpos].clear(); //here all data should be deleted
				  empty_dense_data.push_back(tpos);
				}
#if !defined(LAZY_SPARSE_ALLOCATION)
				else was_sparse[ElementNum(mask)] = true;
#endif
				tag.SetPosition(ENUMUNDEF,mask);
			}
			else delete_entirely = false;
		}
		if( delete_entirely )
		{
			bool flag = false;
#if !defined(LAZY_SPARSE_ALLOCATION)
			bool have_sparse[6] = {false,false,false,false,false,false};
      for(int j = 0; j < 6; j++)
      {
        for(tag_array_type::size_type i = 0; i < tags.size() && !have_sparse[j]; i++)
			    if( tags[i] != tag && tags[i].isSparseByDim(j) ) have_sparse[j] = true;
      }
      for(int j = 0; j < 6; j++) 
				if( was_sparse[j] && !have_sparse[j] )
          sparse_data[j].clear();
#endif
			for(tag_array_type::size_type i = 0; i < tags.size(); i++)
			{
				if( tags[i] == tag )
				{
				  tags.erase(tags.begin()+i);
      		flag = true;
					break;
				}
			}
			assert(flag);
			delete tag.mem;
			tag.mem = NULL;
		}
		return tag;
	}
	
	
	void TagManager::ListTagNames(std::vector<std::string> & list) const
	{
		for(tag_const_iterator it = tags.begin(); it != tags.end(); it++)
			list.push_back(it->GetTagName());
	}
		
	
	bool TagManager::ElementDefined(Tag const & tag, ElementType etype) const
	{
		INMOST_DATA_ENUM_TYPE pos = tag.GetPosition(etype);
		if( pos == ENUMUNDEF ) return false;
		return true;
	}

	void TagManager::ReallocateData(INMOST_DATA_INTEGER_TYPE etypenum, INMOST_DATA_ENUM_TYPE new_size)
	{
		if( new_size < 1024 && etypenum != ElementNum(MESH) ) new_size = 1024;
		if( new_size != 1   && etypenum == ElementNum(MESH) ) new_size = 1;
		if( back_links [etypenum].size() != new_size ) 
			back_links [etypenum].resize(new_size,-1);
#if defined(LAZY_SPARSE_ALLOCATION)
		if( sparse_data[etypenum].size() != new_size ) 
			sparse_data[etypenum].resize(new_size);
#else
		bool need_sparse = false;
#endif
		for(iteratorTag t = tags.begin(); t != tags.end(); ++t)
		{
			if( t->isDefinedByDim(etypenum) ) 
			{
				if( !t->isSparseByDim(etypenum) ) 
					ReallocateData(*t,etypenum,new_size);
#if !defined(LAZY_SPARSE_ALLOCATION)
				else 
					need_sparse = true;
#endif
			}
		}
#if !defined(LAZY_SPARSE_ALLOCATION)
		if( need_sparse )
		{
			if( sparse_data[etypenum].size() != new_size ) 
				sparse_data[etypenum].resize(new_size);
		}
		else if( !sparse_data[etypenum].empty() )
		{
			sparse_data[etypenum].clear();
		}
#endif
	}
	
	void TagManager::ReallocateData(const Tag & t, INMOST_DATA_INTEGER_TYPE etypenum, INMOST_DATA_ENUM_TYPE new_size)
	{
		INMOST_DATA_ENUM_TYPE        data_pos    = t.GetPositionByDim(etypenum);
		INMOST_DATA_ENUM_TYPE        data_size   = t.GetSize();
		TagManager::dense_sub_type & arr         = GetDenseData(data_pos);
		INMOST_DATA_ENUM_TYPE        old_size    = static_cast<INMOST_DATA_ENUM_TYPE>(arr.size());
		DataType                     data_type   = t.GetDataType();
		if( data_size == ENUMUNDEF )
		{
			if( new_size < old_size )
			{
				switch(data_type)
				{
					case DATA_REAL:      for(INMOST_DATA_ENUM_TYPE it = new_size; it < old_size; ++it) {void * p = static_cast<void *>(&arr[it]); if( p != NULL ) (*static_cast<inner_real_array      *>( p )).~inner_real_array();     } break;
					case DATA_INTEGER:   for(INMOST_DATA_ENUM_TYPE it = new_size; it < old_size; ++it) {void * p = static_cast<void *>(&arr[it]); if( p != NULL ) (*static_cast<inner_integer_array   *>( p )).~inner_integer_array();  } break;
					case DATA_BULK:      for(INMOST_DATA_ENUM_TYPE it = new_size; it < old_size; ++it) {void * p = static_cast<void *>(&arr[it]); if( p != NULL ) (*static_cast<inner_bulk_array      *>( p )).~inner_bulk_array();     } break;
					case DATA_REFERENCE: for(INMOST_DATA_ENUM_TYPE it = new_size; it < old_size; ++it) {void * p = static_cast<void *>(&arr[it]); if( p != NULL ) (*static_cast<inner_reference_array *>( p )).~inner_reference_array();} break;
#if defined(USE_AUTODIFF)
          case DATA_VARIABLE:  for(INMOST_DATA_ENUM_TYPE it = new_size; it < old_size; ++it) {void * p = static_cast<void *>(&arr[it]); if( p != NULL ) (*static_cast<inner_variable_array  *>( p )).~inner_variable_array(); } break;
#endif
				}
			}
		}
#if defined(USE_OMP)
#pragma omp critical
#endif
    {
		  arr.resize(new_size);
    }
		if(  data_size == ENUMUNDEF ) //Initialize variable-sized data
		{
			switch(data_type)
			{
				case DATA_REAL:      for(INMOST_DATA_ENUM_TYPE it = old_size; it < new_size; ++it) new ( &arr[it] ) inner_real_array();      break;
				case DATA_INTEGER:   for(INMOST_DATA_ENUM_TYPE it = old_size; it < new_size; ++it) new ( &arr[it] ) inner_integer_array();   break;
				case DATA_BULK:      for(INMOST_DATA_ENUM_TYPE it = old_size; it < new_size; ++it) new ( &arr[it] ) inner_bulk_array();      break;
				case DATA_REFERENCE: for(INMOST_DATA_ENUM_TYPE it = old_size; it < new_size; ++it) new ( &arr[it] ) inner_reference_array(); break;
#if defined(USE_AUTODIFF)
        case DATA_VARIABLE:  for(INMOST_DATA_ENUM_TYPE it = old_size; it < new_size; ++it) new ( &arr[it] ) inner_variable_array();  break;
#endif
			}
		}
	}
	
}
#endif
