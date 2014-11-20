#include "inmost_mesh.h"
#if defined(USE_MESH)
namespace INMOST
{
	template<typename StorageType>
	void ElementArray<StorageType>::Unite(const HandleType * h, INMOST_DATA_ENUM_TYPE num)
	{
		assert(GetMeshLink() == other.GetMeshLink());
		if( empty() ) container.insert(container.end(),h,h+num);
		else
		{
			Mesh * m = GetMeshLink();
			size_type s = size();
			MarkerType mrk = m->CreateMarker();
			SetMarker(mrk);
			for(INMOST_DATA_ENUM_TYPE it = 0; it < num; it++) 
				if( !m->GetMarker(h[it],mrk) ) 
				{
					container.push_back(h[it]);
					m->SetMarker(h[it],mrk);
				}
			RemMarker(mrk);
			m->ReleaseMarker(mrk);
		}
	}

	template<typename StorageType>
	void ElementArray<StorageType>::Subtract(const HandleType * h, INMOST_DATA_ENUM_TYPE num)
	{
		assert(GetMeshLink() == other.GetMeshLink());
		if( !empty() )
		{
			Mesh * mesh = GetMeshLink();
			MarkerType mrk = mesh->CreateMarker();
			//other.SetMarker(mrk);
			mesh->SetMarkerArray(h,num,mrk);
			{
				size_type m = 0, n = 0;
				while( m < size() ) 
				{
					if( !mesh->GetMarker(container[m],mrk) )
						container[n++] = container[m];
					m++;
				}
				container.resize(n);
			}
			//other.RemMarker(mrk);
			mesh->RemMarkerArray(h,num,mrk);
			mesh->ReleaseMarker(mrk);
		}
	}
	template <typename StorageType>
	void ElementArray<StorageType>::Intersect(const HandleType * h, INMOST_DATA_ENUM_TYPE num)
	{
		if( !empty() )
		{
			Mesh * mesh = GetMeshLink();
			MarkerType mrk = mesh->CreateMarker();
			//other.SetMarker(mrk);
			mesh->SetMarkerArray(h,num,mrk);
			{
				size_type m = 0, n = 0;
				while( m < size() ) 
				{
					if( mesh->GetMarker(container[m],mrk) )
						container[n++] = container[m];
					m++;
				}
				container.resize(n);
			}
			//other.RemMarker(mrk);
			mesh->RemMarkerArray(h,num,mrk);
			mesh->ReleaseMarker(mrk);
		}
	}

	template <typename StorageType>
	void ElementArray<StorageType>::SetMarker(MarkerType m) const
	{
		Mesh * mesh = GetMeshLink();
		for(size_type it = 0; it < size(); it++) mesh->SetMarker(container[it],m);
	}

	template <typename StorageType>
	void ElementArray<StorageType>::RemMarker(MarkerType m) const
	{
		Mesh * mesh = GetMeshLink();
		for(size_type it = 0; it < size(); it++) mesh->RemMarker(container[it],m);
	}

	//all possible templates
	template class ElementArray<Element>;
	template class ElementArray<Node>;
	template class ElementArray<Edge>;
	template class ElementArray<Face>;
	template class ElementArray<Cell>;
	template class ElementArray<Storage>;
}
#endif