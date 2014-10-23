#include "inmost.h"
#if defined(USE_MESH)

namespace INMOST
{
	
	
	TagMemory::TagMemory(Mesh * m, const TagMemory & other)
	{
		for(int i = 0; i < __NET; i++)
		{
			pos[i]	= other.pos[i];
			sparse[i] = other.sparse[i];
		}
		tagname		= other.tagname;
		dtype		= other.dtype;
		bulk_data_type = other.bulk_data_type;
		tag_manager = other.tag_manager;
		m_link = m;
		size = other.size;
		record_size = other.record_size;
	}
	
	TagMemory & TagMemory::operator =(TagMemory const & other)
	{
		for(int i = 0; i < __NET; i++)
		{
			pos[i]  = other.pos[i];
			sparse[i] = other.sparse[i];
		}
		tagname         = other.tagname;
		dtype           = other.dtype;
		bulk_data_type = other.bulk_data_type;
		tag_manager = other.tag_manager;
		m_link = other.m_link;
		size = other.size;
		record_size = other.record_size;
		return *this;	
	}
	
	
	TagMemory::~TagMemory()
	{
	}
	
	TagMemory::TagMemory()
	{
		for(int i = 0; i < __NET; i++)
		{
			pos[i]	= ENUMUNDEF;
			sparse[i] = false;
		}
		tagname = "";
	}
	
	Tag::Tag(Mesh * m, std::string name, DataType _dtype,INMOST_DATA_ENUM_TYPE size)
	{
		mem = new TagMemory();
		for(int i = 0; i < __NET; i++)
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
			case DATA_REFERENCE: mem->bulk_data_type = INMOST_MPI_DATA_ENUM_TYPE; break;
		}
		if(mem->size == ENUMUNDEF )
			mem->record_size = VariableDataSize(mem->dtype);
		else
			mem->record_size = mem->size * DataTypeBytesSize(mem->dtype);
		mem->tag_manager = NULL;
		mem->m_link = m;
	}
	TagManager::TagManager()
	{
	}
	TagManager::TagManager(Mesh * m, const TagManager & other)
		:tags(64), empty_dense_data(64), dense_data(64)
	{
		tags.resize(other.tags.size());
		dense_data.resize(other.dense_data.size());
		for(unsigned int i = 0; i < other.tags.size(); i++)
		{
			tags[i].mem = new TagMemory(m,*other.tags[i].mem);
			tags[i].SetTagManager(this);
			for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if( tags[i].isDefined(etype) && !tags[i].isSparse(etype) )
					tags[i].AllocateData(etype);
		}
		
		//The rest of the data should be allocated and copied while copying data of
		//individual elements
	}

	TagManager & TagManager::assign(Mesh * m, TagManager const & other)
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
		dense_data.resize(other.dense_data.size());
		for(unsigned int i = 0; i < other.tags.size(); i++)
		{
			tags[i].mem = new TagMemory(m,*other.tags[i].mem);
			tags[i].SetTagManager(this);
			for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
				if( tags[i].isDefined(etype) && !tags[i].isSparse(etype) )
					tags[i].AllocateData(etype);
		}
		
		//~ sparse_data.resize(other.sparse_data.size());
		//The rest of the data should be allocated and copied while copying data of
		//individual elements
		return *this;
	}
	
	TagManager::~TagManager()
	{
		for(tag_iterator it = tags.begin(); it != tags.end(); it++) delete it->mem;
		tags.clear();
	}
	
	
	
	Tag TagManager::CreateTag(Mesh *m, std::string name, DataType dtype, ElementType etype,ElementType sparse, INMOST_DATA_ENUM_TYPE size)
	{
		Tag new_tag;
		for(INMOST_DATA_ENUM_TYPE i = 0; i < tags.size(); i++)
			if( tags[i].GetTagName() == name )
			{
				if( tags[i].GetDataType() != dtype || (size != ENUMUNDEF && size != tags[i].GetSize()) )
				{
					throw TagExists;
				}
				new_tag = tags[i];
				break;
				//for(ElementType mask = NODE; mask <= MESH; mask = mask << 1)
				//	if( (etype&mask) && tags[i].GetPosition(etype&mask) == ENUMUNDEF )
				//	{
				//		if( sparse & mask ) 
				//		{
				//			//~ tags[i].SetPosition(sparse_data.size(),mask);
				//			//~ sparse_data.resize(sparse_data.size()+1);
				//			tags[i].SetSparse(mask);
				//		}
				//		else 
				//		{
				//			INMOST_DATA_ENUM_TYPE new_pos = dense_data.size();
				//			if( !empty_dense_data.empty() )
				//			{
				//				new_pos = empty_dense_data.back();
				//				empty_dense_data.pop_back();
				//				dense_data[new_pos] = chunk_array<INMOST_DATA_BULK_TYPE>(tags[i].GetRecordSize()*4096);
				//			}
				//			else dense_data.push_back(chunk_array<INMOST_DATA_BULK_TYPE>(tags[i].GetRecordSize()*4096));
				//			tags[i].SetPosition(new_pos,mask);
				//			tags[i].AllocateData(mask);
				//		}
				//	}
				//return tags[i];
			}
		if( !new_tag.isValid() )
		{
			new_tag = Tag(m,name,dtype,size);
			new_tag.SetTagManager(this);
			tags.push_back(new_tag);
		}
		for(ElementType mask = NODE; mask <= MESH; mask = mask << 1)
		{
			if( (mask & etype) && new_tag.GetPosition(etype&mask) == ENUMUNDEF ) 
			{
				if( sparse & mask )
				{
					//~ new_tag.SetPosition(sparse_data.size(),mask);
					//~ sparse_data.resize(sparse_data.size()+1);
					new_tag.SetPosition(ENUMUNDEF-1,mask);
					new_tag.SetSparse(mask);
				}
				else
				{
					INMOST_DATA_ENUM_TYPE new_pos = dense_data.size();
					if( !empty_dense_data.empty() )
					{
						new_pos = empty_dense_data.back();
						empty_dense_data.pop_back();
						dense_data[new_pos] = chunk_array<INMOST_DATA_BULK_TYPE>(new_tag.GetRecordSize()*8192);
					}
					else dense_data.push_back(chunk_array<INMOST_DATA_BULK_TYPE>(new_tag.GetRecordSize()*8192));
					new_tag.SetPosition(new_pos,mask);
					new_tag.AllocateData(mask);
				}
			}
		}
		
		
		return new_tag;
	}
	Tag TagManager::GetTag(std::string name) const
	{
		for(INMOST_DATA_BIG_ENUM_TYPE i = 0; i < tags.size(); i++)
			if( tags[i].GetTagName() == name )
				return tags[i];
		throw TagNotFound;
	}
	bool TagManager::HaveTag(std::string name) const
	{
		for(INMOST_DATA_BIG_ENUM_TYPE i = 0; i < tags.size(); i++)
			if( tags[i].GetTagName() == name )
				return true;
		return false;
	}
	Tag TagManager::DeleteTag(Tag tag, ElementType type_mask)
	{
		bool delete_entirely = true;
		INMOST_DATA_ENUM_TYPE tpos;//,ipos;
		for(ElementType mask = NODE; mask <= MESH; mask = mask << 1 )
		{
			tpos = tag.GetPosition(mask);
			if( tpos == ENUMUNDEF ) continue;
			if( mask & type_mask )
			{
				//~ if( tag.isSparse(mask) ) {}//sparse_data.erase(sparse_data.begin()+tpos);
				//~ else dense_data.erase(dense_data.begin()+tpos);
				//~ for(size_t i = 0; i < tags.size(); i++) 
				//~ {
					//~ for(ElementType imask = NODE; imask <= MESH; imask = imask << 1 )
						//~ if( tag.isSparse(mask) == tags[i].isSparse(imask) )
						//~ {
							//~ ipos = tags[i].GetPosition(imask);
							//~ if( ipos == ENUMUNDEF ) continue;
							//~ if( ipos > tpos ) 
							//~ {
								//~ tags[i].SetPosition(ipos-1,imask);
							//~ }
						//~ }
				//~ }
				if( !tag.isSparse(mask) ) 
				{
					dense_data[tpos].clear(); //here all data should be deleted
					empty_dense_data.push_back(tpos);
				}
				tag.SetPosition(ENUMUNDEF,mask);
			}
			else delete_entirely = false;
		}
		if( delete_entirely )
		{
			bool flag = false;
			for(size_t i = 0; i < tags.size(); i++)
				if( tags[i] == tag )
				{
					tags.erase(tags.begin()+i);
					flag = true;
					break;
				}
			if( !flag ) throw TagNotFound;
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
	
	void Tag::ShrinkData(ElementType type, size_t new_size)
	{
		size_t old_size;
		if( isSparse(type) ) return;
		INMOST_DATA_ENUM_TYPE data_pos = GetPosition(type);
		if( data_pos == ENUMUNDEF ) return;
		INMOST_DATA_ENUM_TYPE data_size = GetSize(), record_size;
		size_t bytes = GetBytesSize();
		record_size = (data_size == ENUMUNDEF ? VariableDataSize(GetDataType()) : bytes * data_size);
		TagManager::dense_sub_type & arr = GetTagManager()->GetDenseData(data_pos);
		old_size = arr.size()/record_size;
		//if( new_size < 1024 ) new_size = 1024;
		if( data_size == ENUMUNDEF )
		{
			if( new_size*record_size < arr.size() )
			{
				TagManager::dense_sub_type::iterator it = arr.begin() + new_size*record_size;
				switch(GetDataType())
				{
					case DATA_REAL:      while( it < arr.end() ) { void * p = static_cast<void *>(&*it); if( p != NULL ) (*static_cast<inner_real_array      *>( p )).~inner_real_array(); it+=record_size; }  break;
					case DATA_INTEGER:   while( it < arr.end() ) { void * p = static_cast<void *>(&*it); if( p != NULL ) (*static_cast<inner_integer_array   *>( p )).~inner_integer_array(); it+=record_size; }  break;
					case DATA_BULK:      while( it < arr.end() ) { void * p = static_cast<void *>(&*it); if( p != NULL ) (*static_cast<inner_bulk_array      *>( p )).~inner_bulk_array(); it+=record_size; }  break;
					case DATA_REFERENCE: while( it < arr.end() ) { void * p = static_cast<void *>(&*it); if( p != NULL ) (*static_cast<inner_reference_array *>( p )).~inner_reference_array(); it+=record_size; }  break;
				}
			}
		}
		arr.resize(new_size*record_size,0);
		if( data_size == ENUMUNDEF )
		{
			if( old_size*record_size < arr.size() )
			{
				TagManager::dense_sub_type::iterator it = arr.begin() + old_size*record_size;
				switch(GetDataType())
				{
					case DATA_REAL:      while( it < arr.end() ) { new (&*it) inner_real_array();      it+=record_size; }  break;
					case DATA_INTEGER:   while( it < arr.end() ) { new (&*it) inner_integer_array();   it+=record_size; }  break;
					case DATA_BULK:      while( it < arr.end() ) { new (&*it) inner_bulk_array();      it+=record_size; }  break;
					case DATA_REFERENCE: while( it < arr.end() ) { new (&*it) inner_reference_array(); it+=record_size; }  break;
				}
			}
		}
		//std::cout << "Shrink " << GetTagName() << " for " << ElementTypeName(type) << " from " << old_size << " to " << new_size << std::endl;
		//~ TagManager::dense_sub_type(arr).swap(arr);
	}
	
	void Tag::AllocateData(ElementType t)
	{
		Mesh * m = GetMeshLink();
		if( m == NULL ) return;
		TagManager::dense_sub_type & arr = GetTagManager()->GetDenseData(GetPosition(t));
		INMOST_DATA_ENUM_TYPE record_size = GetRecordSize();
		INMOST_DATA_ENUM_TYPE old_size = arr.size();
		INMOST_DATA_ENUM_TYPE new_size = m->GetArrayCapacity(t);
		if( new_size < 1024 ) new_size = 1024;
		new_size *= record_size;
		arr.resize(new_size);
		//std::cout << "tag " << GetTagName() << " was " << old_size << " now " << new_size << std::endl;
		if( GetSize() == ENUMUNDEF ) //Initialize variable-sized data
		{
			TagManager::dense_sub_type::iterator it = arr.begin() + old_size;
			switch(GetDataType())
			{
				case DATA_REAL:      while( it < arr.end() ) {new ( &*it ) inner_real_array();      it+=record_size; }  break;
				case DATA_INTEGER:   while( it < arr.end() ) {new ( &*it ) inner_integer_array();   it+=record_size; }  break;
				case DATA_BULK:      while( it < arr.end() ) {new ( &*it ) inner_bulk_array();      it+=record_size; }  break;
				case DATA_REFERENCE: while( it < arr.end() ) {new ( &*it ) inner_reference_array(); it+=record_size; }  break;
			}
		}
		//else if( new_size-old_size > 0 ) memset(&arr[old_size],0,new_size-old_size);
	}
	
	
	
	
}
#endif
