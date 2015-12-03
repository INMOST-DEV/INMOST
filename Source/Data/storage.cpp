#include "inmost.h"

#if defined(USE_MESH)
#include <new>

namespace INMOST
{
	Storage & Storage::operator =(Storage const & other) 
	{
		handle = other.handle; 
		if( handle_link != NULL ) *handle_link = handle; 
		//else handle_link = other.handle_link; //if other have remote link this will copy this link and current will also be remote
		m_link = other.m_link; 
		return *this;
	}
	
  Element MakeElement(const RemoteHandleType & rh)
  {
    return Element(rh.first,rh.second);
  }

  Element MakeElementRef(RemoteHandleType & rh)
  {
    return Element(rh.first,&rh.second);
  }
	
	Element Storage::reference_array::operator [](size_type n)  
	{
		return Element(m,&shell<HandleType>::operator[](n));
	}

	Element Storage::reference_array::operator [](size_type n) const 
	{
		return Element(m,shell<HandleType>::operator[](n));
	}

	Element Storage::reference_array::iterator::operator->()
	{
		return Element(m,&shell<HandleType>::iterator::operator *());
	}

	Element Storage::reference_array::const_iterator::operator->()
	{
		return Element(m,shell<HandleType>::const_iterator::operator *());
	}
	
	Element Storage::reference_array::reverse_iterator::operator->()
	{
		return Element(m,&shell<HandleType>::reverse_iterator::operator *());
	}

	Element Storage::reference_array::const_reverse_iterator::operator->()
	{
		return Element(m,shell<HandleType>::const_reverse_iterator::operator *());
	}

  void Storage::reference_array::push_back(const Storage & ref) 
  {
    shell<reference>::push_back(ref->GetHandle());
  }

  Element Storage::remote_reference_array::operator [](size_type n)  
	{
    
		return MakeElementRef(shell<remote_reference>::operator[](n));
	}

	Element Storage::remote_reference_array::operator [](size_type n) const 
	{
		return MakeElement(shell<remote_reference>::operator[](n));
	}

	Element Storage::remote_reference_array::iterator::operator->()
	{
		return MakeElementRef(shell<remote_reference>::iterator::operator *());
	}

	Element Storage::remote_reference_array::const_iterator::operator->()
	{
		return MakeElement(shell<remote_reference>::const_iterator::operator *());
	}
	
	Element Storage::remote_reference_array::reverse_iterator::operator->()
	{
		return MakeElementRef(shell<remote_reference>::reverse_iterator::operator *());
	}

	Element Storage::remote_reference_array::const_reverse_iterator::operator->()
	{
		return MakeElement(shell<remote_reference>::const_reverse_iterator::operator *());
	}

  void Storage::remote_reference_array::push_back(const Storage & ref) 
  {
    shell<remote_reference>::push_back(RemoteHandleType(ref->GetMeshLink(),ref->GetHandle()));
  }
}

#endif
