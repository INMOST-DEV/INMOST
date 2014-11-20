#include "inmost.h"

#if defined(USE_MESH)
#include <new>

namespace INMOST
{
	
	INMOST_DATA_ENUM_TYPE Storage::GetDataSize(const Tag & tag) const
	{
		return GetMeshLink()->GetDataSize(GetHandle(),tag);
	}

	void Storage::SetDataSize(const Tag & tag,INMOST_DATA_ENUM_TYPE new_size) const
	{
		GetMeshLink()->SetDataSize(GetHandle(),tag,new_size);
	}
	
	void Storage::GetData(const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, void * data_out) const
	{
		GetMeshLink()->GetData(GetHandle(),tag,shift,size,data_out);
	}
	
	void Storage::SetData(const Tag & tag,INMOST_DATA_ENUM_TYPE shift, INMOST_DATA_ENUM_TYPE size, const void * data_in) const
	{
		GetMeshLink()->SetData(GetHandle(),tag,shift,size,data_in);
	}
	
	void Storage::DelData(const Tag & tag) const
	{
		GetMeshLink()->DelData(GetHandle(),tag);
	}

	void Storage::DelDenseData(const Tag & tag) const
	{
		GetMeshLink()->DelDenseData(GetHandle(),tag);
	}

	void Storage::DelSparseData(const Tag & tag) const
	{
		GetMeshLink()->DelSparseData(GetHandle(),tag);
	}


	void Storage::SetMarker(MarkerType n)  const
	{
		assert( isValid() );
		GetMeshLink()->SetMarker(GetHandle(),n);
	}
	bool Storage::GetMarker(MarkerType n) const  
	{
		assert( isValid() );
		return GetMeshLink()->GetMarker(GetHandle(),n);
	}
	void Storage::RemMarker(MarkerType n)  const
	{
		assert( isValid() );
		GetMeshLink()->RemMarker(GetHandle(),n);
	}
	void Storage::ClearMarkerSpace()  const
	{
		GetMeshLink()->ClearMarkerSpace(GetHandle());
	}
	void Storage::GetMarkerSpace(Storage::bulk copy[MarkerFields]) const 
	{
		GetMeshLink()->GetMarkerSpace(GetHandle(),copy);
	}
	void Storage::SetMarkerSpace(Storage::bulk source[MarkerFields])  const
	{
		GetMeshLink()->SetMarkerSpace(GetHandle(),source);
	}

	bool Storage::HaveData(const Tag & tag) const
	{
		assert(isValid());
		return GetMeshLink()->HaveData(GetHandle(),tag);
	}
	Storage::real      &              Storage::Real            (const Tag & tag) const
	{
		return GetMeshLink()->Real(GetHandle(),tag);
	}
	
	Storage::integer   &              Storage::Integer         (const Tag & tag)  const
	{
		
		return GetMeshLink()->Integer(GetHandle(),tag);
	}
	
	Storage::bulk      &              Storage::Bulk            (const Tag & tag)  const
	{
		return GetMeshLink()->Bulk(GetHandle(),tag);
	}
	
	Storage::reference &              Storage::Reference       (const Tag & tag)  const
	{
		return GetMeshLink()->Reference(GetHandle(),tag);
	}
	
	Storage::real_array               Storage::RealArray       (const Tag & tag)  const
	{
		return GetMeshLink()->RealArray(GetHandle(),tag);
	}
	
	Storage::integer_array            Storage::IntegerArray    (const Tag & tag)  const
	{
		return GetMeshLink()->IntegerArray(GetHandle(),tag);
	}
	
	Storage::bulk_array               Storage::BulkArray       (const Tag & tag)  const
	{
		return GetMeshLink()->BulkArray(GetHandle(),tag);
	}
	
	Storage::reference_array          Storage::ReferenceArray  (const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceArray(GetHandle(),tag);
	}
		
	Storage::real_array               Storage::RealArrayDF     (const Tag & tag)  const
	{
		return GetMeshLink()->RealArrayDF(GetHandle(),tag);
	}
	Storage::integer_array            Storage::IntegerArrayDF  (const Tag & tag)  const
	{
		return GetMeshLink()->IntegerArrayDF(GetHandle(),tag);
	}
	Storage::bulk_array               Storage::BulkArrayDF     (const Tag & tag)  const
	{
		return GetMeshLink()->BulkArrayDF(GetHandle(),tag);
	}
	Storage::reference_array          Storage::ReferenceArrayDF(const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceArrayDF(GetHandle(),tag);
	}
	Storage::real      &              Storage::RealDF          (const Tag & tag)  const
	{
		return GetMeshLink()->RealDF(GetHandle(),tag);
	}
	Storage::integer   &              Storage::IntegerDF       (const Tag & tag)  const
	{
		return GetMeshLink()->IntegerDF(GetHandle(),tag);
	}
	Storage::bulk      &              Storage::BulkDF          (const Tag & tag)  const
	{
		return GetMeshLink()->BulkDF(GetHandle(),tag);
	}
	Storage::reference &              Storage::ReferenceDF     (const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceDF(GetHandle(),tag);
	}
		
	Storage::real_array               Storage::RealArrayDV     (const Tag & tag)  const
	{
		return GetMeshLink()->RealArrayDV(GetHandle(),tag);
	}
	Storage::integer_array            Storage::IntegerArrayDV  (const Tag & tag)  const
	{
		return GetMeshLink()->IntegerArrayDV(GetHandle(),tag);	
	}
	Storage::bulk_array               Storage::BulkArrayDV     (const Tag & tag)  const
	{
		return GetMeshLink()->BulkArrayDV(GetHandle(),tag);
	}
	Storage::reference_array          Storage::ReferenceArrayDV(const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceArrayDV(GetHandle(),tag);
	}
	Storage::real      &              Storage::RealDV          (const Tag & tag)  const
	{
		return GetMeshLink()->RealDV(GetHandle(),tag);
	}
	Storage::integer   &              Storage::IntegerDV       (const Tag & tag)  const
	{
		return GetMeshLink()->IntegerDV(GetHandle(),tag);
	}
	Storage::bulk      &              Storage::BulkDV          (const Tag & tag)  const
	{
		return GetMeshLink()->BulkDV(GetHandle(),tag);
	}
	Storage::reference &              Storage::ReferenceDV     (const Tag & tag)  const
	{
		return GetMeshLink()->ReferenceDV(GetHandle(),tag);
	}
		
	bool Storage::isValid() const 
	{
		return handle != InvalidHandle() && GetMeshLink()->isValidElement(handle);
	}

	Element Storage::reference_array::operator [](size_type n) const 
	{
		return Element(m,shell<HandleType>::operator[](n));
	}

	Element Storage::reference_array::iterator::operator->()
	{
		return Element(m,shell<HandleType>::iterator::operator *());
	}

	Element Storage::reference_array::const_iterator::operator->()
	{
		return Element(m,shell<HandleType>::const_iterator::operator *());
	}

	Element Storage::reference_array::reverse_iterator::operator->()
	{
		return Element(m,shell<HandleType>::reverse_iterator::operator *());
	}

	Element Storage::reference_array::const_reverse_iterator::operator->()
	{
		return Element(m,shell<HandleType>::const_reverse_iterator::operator *());
	}
}

#endif
