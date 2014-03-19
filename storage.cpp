#include "inmost.h"

#if defined(USE_MESH)
#include <new>

namespace INMOST
{
	void Storage::CopyData(Tag t, void * adata, void * bdata)
	{
		
		INMOST_DATA_ENUM_TYPE data_size = t.GetSize();
		size_t bytes = t.GetBytesSize();
		if( data_size == ENUMUNDEF ) //variable size array
		{
			DataType type = t.GetDataType();
			if( adata != NULL ) DestroyVariableData(t,adata);
			if( type == DATA_REAL )           new (adata) inner_real_array     (*static_cast< inner_real_array      * >(bdata));
			else if( type == DATA_INTEGER )   new (adata) inner_integer_array  (*static_cast< inner_integer_array   * >(bdata));
			else if( type == DATA_BULK )      new (adata) inner_bulk_array     (*static_cast< inner_bulk_array      * >(bdata));
			else if( type == DATA_REFERENCE ) new (adata) inner_reference_array(*static_cast< inner_reference_array * >(bdata));
		}
		else // fixed size array
			memcpy(adata,bdata,data_size*bytes);
	}
	void Storage::DestroyVariableData(Tag t, void * adata)
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
		}
	}
	
	void Storage::MoveData(int new_local_id)
	{
		Mesh * m = GetMeshLink();
		if( m == NULL ) return;
		for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); t++) if( !t->isSparse(GetElementType()) )
		{
			INMOST_DATA_ENUM_TYPE data_pos = t->GetPosition(GetElementType());
			if( data_pos == ENUMUNDEF ) continue;
			TagManager::dense_sub_type & arr = t->GetTagManager()->GetDenseData(data_pos);
			INMOST_DATA_ENUM_TYPE record_size = t->GetRecordSize();
			INMOST_DATA_ENUM_TYPE from = record_size*local_id, to = record_size*new_local_id;
			memcpy(&arr[to],&arr[from],record_size);
			memset(&arr[from],0,record_size);
		}
	}
	void Storage::SwapData(int new_local_id)
	{
		Mesh * m = GetMeshLink();
		if( m == NULL ) return;
		dynarray<bulk,2048> temp;
		for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); t++) if( !t->isSparse(GetElementType()) )
		{

			INMOST_DATA_ENUM_TYPE data_pos = t->GetPosition(GetElementType());
			if( data_pos == ENUMUNDEF ) continue;
			TagManager::dense_sub_type & arr = t->GetTagManager()->GetDenseData(data_pos);
			INMOST_DATA_ENUM_TYPE record_size = t->GetRecordSize();
			INMOST_DATA_ENUM_TYPE from = record_size*local_id, to = record_size*new_local_id;
			
			//std::cout << "swap " << t->GetTagName() << " from " << local_id << " to " << new_local_id << " record size: " << record_size << " bytes: " << bytes << " data size: " << data_size << std::endl;
			
			temp.resize(record_size);
			memcpy(&temp[0],&arr[to],record_size);
			memcpy(&arr[to],&arr[from],record_size);
			memcpy(&arr[from],&temp[0],record_size);
		}
	}
	Storage::Storage(Mesh * m, ElementType _etype)
	{
		if( sizeof(INMOST_DATA_BULK_TYPE) != 1 ) throw BadBulkType;
		INMOST_DATA_ENUM_TYPE i = 0;
#if !defined(NDEBUG)
		for(ElementType e = NODE; e <= MESH; e = e << 1 ) if( _etype & e ) i++;
		assert( i != 0 );
		assert( i == 1 );
#endif
		etypenum = ElementNum(_etype);
		markers = 0;
		m_link = m;
		if( m != NULL ) m->TieElement(this);
	}
	Storage::Storage(Mesh * m, INMOST_DATA_ENUM_TYPE lid, const Storage & other)
	{
		etypenum = other.etypenum;
		markers = other.markers;
		m_link = m;
		local_id = lid;
		//if( m != NULL ) m->TieElement(this);
		Mesh * other_mesh = other.GetMeshLink();
		void * a, * b;
		if( m == NULL || other_mesh == NULL ) return;

		for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); t++) if( t->isDefined(GetElementType()) )
		{
			Mesh::iteratorTag other_t = other_mesh->BeginTag() + (t-m->BeginTag());
			b = other.GetLink(*other_t);
			if( b != NULL )
			{
				a = GetLink(*t);
				CopyData(*t,a,b);
			}
		}
	}

	Storage::Storage(const Storage & other) {throw NotImplemented;}

	
	Storage & Storage::assign(Mesh * m, INMOST_DATA_ENUM_TYPE lid, Storage const & other)
	{
		Mesh * other_mesh = other.GetMeshLink();
		etypenum = other.etypenum;
		markers = other.markers;
		
		m_link = m;
		local_id = lid;
		
		void * a, * b;
		if( m == NULL || other_mesh == NULL ) return *this;


		for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); t++) if( t->isDefined(GetElementType()) )
		{
			Mesh::iteratorTag other_t = other_mesh->BeginTag() + (t-m->BeginTag());
			b = other.GetLink(*other_t);
			if( b != NULL )
			{
				a = GetLink(*t);
				CopyData(*t,a,b);
			}
		}
		
		return *this;
	}
	
	
	size_t Storage::GetDataSize(Tag tag) const
	{
		if( tag.GetSize() == ENUMUNDEF )
		{
			void * adata = GetLink(tag);
			if( adata == NULL ) throw NoTagPosition;
			switch(tag.GetDataType())
			{
				case DATA_REAL:     return static_cast<inner_real_array     *>(adata)->size();
				case DATA_INTEGER:  return static_cast<inner_integer_array  *>(adata)->size();
				case DATA_BULK:     return static_cast<inner_bulk_array     *>(adata)->size();
				case DATA_REFERENCE:return static_cast<inner_reference_array*>(adata)->size();
			}
			throw BadTag;
		}
		return tag.GetSize();
	}
	void Storage::SetDataSize(Tag tag,size_t new_size)
	{
		assert( tag.GetMeshLink() == GetMeshLink() );
		void * adata = GetLink(tag);
		if( adata == NULL ) throw NoTagPosition;
		if( tag.GetSize() == ENUMUNDEF )
		{
			switch(tag.GetDataType())
			{
				case DATA_REAL:     static_cast<inner_real_array      *>(adata)->resize(new_size); break;
				case DATA_INTEGER:  static_cast<inner_integer_array   *>(adata)->resize(new_size); break;
				case DATA_BULK:     static_cast<inner_bulk_array      *>(adata)->resize(new_size); break;
				case DATA_REFERENCE:static_cast<inner_reference_array *>(adata)->resize(new_size); break;
			}
			return;
		}
		else if( tag.GetSize() == new_size )
			return;
		throw BadTag;
	}
	void Storage::GetData(Tag tag,size_t shift, size_t size, void * data_out) const
	{
		assert( tag.GetMeshLink() == GetMeshLink() );
		void * adata = GetLink(tag);
		if( adata == NULL ) throw NoData;
		INMOST_DATA_ENUM_TYPE data_size = tag.GetSize();
		size_t bytes = tag.GetBytesSize();
		if( data_size == ENUMUNDEF )
		{
			switch(tag.GetDataType())
			{
				case DATA_REAL:     memcpy(data_out,&(*static_cast<inner_real_array      *>(adata))[shift],bytes*size); break;
				case DATA_INTEGER:  memcpy(data_out,&(*static_cast<inner_integer_array   *>(adata))[shift],bytes*size); break;
				case DATA_BULK:     memcpy(data_out,&(*static_cast<inner_bulk_array      *>(adata))[shift],bytes*size); break;
				case DATA_REFERENCE:memcpy(data_out,&(*static_cast<inner_reference_array *>(adata))[shift],bytes*size); break;
			}
		}
		else memcpy(data_out,static_cast<INMOST_DATA_BULK_TYPE *>(adata)+shift*bytes,size*bytes);
		return;
	}
	void Storage::SetData(Tag tag,size_t shift, size_t size, void * data_in)
	{
		assert( tag.GetMeshLink() == GetMeshLink() );
		void * adata = GetLink(tag);
		if( adata == NULL ) throw NoData;
		INMOST_DATA_ENUM_TYPE data_size = tag.GetSize();
		size_t bytes = tag.GetBytesSize();
		if( data_size == ENUMUNDEF )
		{
			switch(tag.GetDataType())
			{
				case DATA_REAL:     memcpy(&(*static_cast<inner_real_array     *>(adata))[shift],data_in,bytes*size); break;
				case DATA_INTEGER:  memcpy(&(*static_cast<inner_integer_array  *>(adata))[shift],data_in,bytes*size); break;
				case DATA_BULK:     memcpy(&(*static_cast<inner_bulk_array     *>(adata))[shift],data_in,bytes*size); break;
				case DATA_REFERENCE:memcpy(&(*static_cast<inner_reference_array*>(adata))[shift],data_in,bytes*size); break;
			}
		}
		else memcpy(static_cast<INMOST_DATA_BULK_TYPE *>(adata)+shift*bytes,data_in,size*bytes);
		return;
	}
	Storage::~Storage()
	{
		Mesh * m = GetMeshLink();
		if( m == NULL ) return;
		for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); t++)
			if( t->isDefined(GetElementType()) ) DelData(*t);
		
		if( !(GetElementType() & MESH) )
			m->UntieElement(this);
	}
	void Storage::DelData(Tag tag)
	{
		assert( tag.GetMeshLink() == GetMeshLink() );
		if( tag.isSparse(GetElementType()) )
		{
			/* old sparse data
			TagManager::sparse_sub_type & arr = tag.GetTagManager()->GetSparseData(tag.GetPosition(GetElementType()));
			TagManager::sparse_sub_type::iterator it = arr.find(this);
			if( it != arr.end() )
			{
				if( tag.GetSize() == ENUMUNDEF ) DestroyVariableData(tag,it->second);
				free(it->second);
				arr.erase(it);
			}
			*/
			for(sparse_data_array_type::iterator it = inner_data.begin(); it != inner_data.end(); ++it)
				if( it->first == tag ) 
				{
					if( tag.GetSize() == ENUMUNDEF ) DestroyVariableData(tag,it->second);
					free(it->second);
					it = inner_data.erase(it);
					break;
				}
		}
		else 
		{
			void * data = GetLink(tag);
			if( data != NULL )
			{
				if( tag.GetSize() == ENUMUNDEF )
					DestroyVariableData(tag,data);
				else memset(data,0,tag.GetRecordSize());
			}
		}
	}

	

}

#endif
