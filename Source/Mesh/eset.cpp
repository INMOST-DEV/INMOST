//#define _HAS_ITERATOR_DEBUGGING 0

#include "inmost.h"
#if defined(USE_MESH)

namespace INMOST
{
	//shortcuts
	

	ElementSet::iterator & ElementSet::iterator::operator ++() 
	{
		++pos;
		Element::adj_type::size_type end = ptr->size();
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			while(pos < end && (ptr->at(pos) == InvalidHandle() || m->GetMarker(ptr->at(pos),hm)) ) ++pos; 
		}
		else 
		{
			while(pos < end && ptr->at(pos) == InvalidHandle()) ++pos; 
		}
		return *this;
	}


	Storage::enumerator ElementSet::nbSorted() const
	{
		return hSorted(GetMeshLink()->HighConn(GetHandle()));
	}
	ElementSet ElementSet::GetParent() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & hc = GetMeshLink()->HighConn(GetHandle());
		return ElementSet(GetMeshLink(),hParent(hc));
	}

	ElementSet ElementSet::GetSibling() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & hc = GetMeshLink()->HighConn(GetHandle());
		return ElementSet(GetMeshLink(),hSibling(hc));
	}

	ElementSet ElementSet::GetChild() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & hc = GetMeshLink()->HighConn(GetHandle());
		return ElementSet(GetMeshLink(),hChild(hc));
	}


	bool ElementSet::HaveParent() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & hc = GetMeshLink()->HighConn(GetHandle());
		return hParent(hc) != InvalidHandle();
	}

	bool ElementSet::HaveSibling() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & hc = GetMeshLink()->HighConn(GetHandle());
		return hSibling(hc) != InvalidHandle();
	}

	bool ElementSet::HaveChild() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & hc = GetMeshLink()->HighConn(GetHandle());
		return hChild(hc) != InvalidHandle();
	}

	void ElementSet::AddChild(const ElementSet & child) const
	{
		assert(isValid());
		assert(GetElementType() == ESET);
		assert(child->GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & parent_hc = m->HighConn(GetHandle());
		if(hChild(parent_hc) == InvalidHandle()) 
		{
			Element::adj_type & child_hc = m->HighConn(child->GetHandle());
			assert(hParent(child_hc) == InvalidHandle()); // child is not connected to another parent
			hParent(child_hc) = GetHandle(); //connect child to me
			hChild(parent_hc) = child->GetHandle(); //connect me to child
		}
		else 
			ElementSet(m,hChild(parent_hc))->AddSibling(child); //connect as sibling to my children
	}

	void ElementSet::AddSibling(const ElementSet & sibling) const
	{
		assert(isValid());
		assert(GetElementType() == ESET);
		assert(sibling->GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		//perform checks
		HandleType h = GetHandle();
		bool done = false;
		while(!done)
		{
			Element::adj_type & current_hc = m->HighConn(h);
			assert(hParent(current_hc) != InvalidHandle()); //no parent for current set
			if (hSibling(current_hc) == InvalidHandle()) //no sibling here
			{
				Element::adj_type& sibling_hc = m->HighConn(sibling->GetHandle());
				assert(hParent(sibling_hc) == InvalidHandle());
				hSibling(current_hc) = sibling->GetHandle(); //connect sibling to me
				hParent(sibling_hc) = hParent(current_hc); //connect sibling to my parent
				done = true;
			}
			else h = hSibling(current_hc);
		}
	}

	void ElementSet::RemChild(const ElementSet & child) const
	{
		assert(isValid());
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type& parent_hc = m->HighConn(GetHandle());
		assert(hChild(parent_hc) != InvalidHandle());
		if( hChild(parent_hc) == child->GetHandle() )
		{
			Element::adj_type & child_hc = m->HighConn(child->GetHandle());
			assert(hParent(child_hc) == GetHandle()); //was connected to me
			hParent(child_hc) = InvalidHandle();
			hChild(parent_hc) = hSibling(child_hc); //get next sibling as my child
		}
		else ElementSet(m,hChild(parent_hc))->RemSibling(child);
	}

	void ElementSet::RemSibling(const ElementSet & sibling) const
	{
		assert(isValid());
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		HandleType h = GetHandle();
		bool done = false;
		while (!done)
		{
			Element::adj_type& current_hc = m->HighConn(h);
			assert(hSibling(current_hc) != InvalidHandle());
			if (hSibling(current_hc) == sibling->GetHandle())
			{
				Element::adj_type& sibling_hc = m->HighConn(sibling->GetHandle());
				assert(hParent(sibling_hc) == hParent(current_hc));
				hParent(sibling_hc) = InvalidHandle(); //remove it's parent
				hSibling(current_hc) = hSibling(sibling_hc); //get next sibling as my sibling
				done = true;
			}
			else h = hSibling(current_hc);
		}
	}

	Storage::enumerator ElementSet::CountChildren() const
	{
		assert(GetElementType() == ESET);
		ElementSet child = GetChild();
		return child->isValid() ? child->CountSiblings() : 0;
	}

	Storage::enumerator ElementSet::CountSiblings() const
	{
		assert(GetElementType() == ESET);
		ElementSet next = GetSibling();
		return next->isValid() ? next->CountSiblings()+1 : 1;
	}

	std::string ElementSet::GetName() const
	{
		assert(GetElementType() == ESET);
		std::string ret;
		bulk_array arr = BulkArrayDV(GetMeshLink()->SetNameTag());
		ret.resize(arr.size());
		for(bulk_array::size_type it = 0; it < arr.size(); ++it)
			ret[it] = arr[it];
		return ret;
	}

	HandleType * ElementSet::getHandles() const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		return lc.data();
	}

	Storage::enumerator ElementSet::nbHandles() const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		return static_cast<enumerator>(lc.size());
	}

	void ElementSet::SetMarkerElements(MarkerType mrk, ElementType etype) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)
				if(*it != 0 && (GetHandleElementType(*it) & etype) && !m->GetMarker(*it,hm) ) //check for type filters invalid handles
					m->SetMarker(*it,mrk);
		}
		else
		{
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)//check for type filters invalid handles
				if(*it != 0 && GetHandleElementType(*it) & etype)
					m->SetMarker(*it,mrk);
		}
	}


  void ElementSet::SetPrivateMarkerElements(MarkerType mrk, ElementType etype) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)
				if( *it != InvalidHandle() && (GetHandleElementType(*it) & etype) && !m->GetMarker(*it,hm) ) //check for type filters invalid handles
					m->SetPrivateMarker(*it,mrk);
		}
		else
		{
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)//check for type filters invalid handles
				if( *it != InvalidHandle() && (GetHandleElementType(*it) & etype) )
					m->SetPrivateMarker(*it,mrk);
		}
	}
	void ElementSet::RemMarkerElements(MarkerType mrk, ElementType etype) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)
				if( *it != InvalidHandle() && (GetHandleElementType(*it) & etype) && !m->GetMarker(*it,hm) )//check for type filters invalid handles
					m->RemMarker(*it,mrk);
		}
		else
		{
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)//check for type filters invalid handles
				if( *it != InvalidHandle() && (GetHandleElementType(*it) & etype) )
					m->RemMarker(*it,mrk);
		}
	}
  void ElementSet::RemPrivateMarkerElements(MarkerType mrk, ElementType etype) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)
				if( *it != InvalidHandle() && (GetHandleElementType(*it) & etype) && !m->GetMarker(*it,hm) )//check for type filters invalid handles
					m->RemPrivateMarker(*it,mrk);
		}
		else
		{
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)//check for type filters invalid handles
				if( *it != InvalidHandle() && (GetHandleElementType(*it) & etype) )
					m->RemPrivateMarker(*it,mrk);
		}
	}

	
	ElementArray<Element> ElementSet::getAdjElements(ElementType etype) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) && !m->GetMarker(lc[it],hm) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)//check for type filters invalid handles
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) )
					ret.push_back(lc[it]);
		}
		return ret;
	}

	Storage::enumerator ElementSet::nbAdjElements(ElementType etype) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type const & lc = m->LowConn(GetHandle());
		enumerator ret = 0;
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) && !m->GetMarker(lc[it],hm) )//check for type filters invalid handles
					++ret;
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)//check for type filters invalid handles
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) )
					++ret;
		}
		return ret;
	}

	ElementArray<Element> ElementSet::getAdjElements(ElementType etype, MarkerType select, bool invert) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) && !m->GetMarker(lc[it],hm) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		return ret;
	}

	Storage::enumerator ElementSet::nbAdjElements(ElementType etype, MarkerType select, bool invert) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		enumerator ret = 0;
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) && !m->GetMarker(lc[it],hm) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					++ret;
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & etype) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					++ret;
		}
		return ret;
	}

	ElementArray<Node> ElementSet::getNodes() const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Node> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & NODE) && !m->GetMarker(lc[it],hm) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)//check for type filters invalid handles
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & NODE) )
					ret.push_back(lc[it]);
		}
		return ret;
	}

	ElementArray<Node> ElementSet::getNodes(MarkerType select, bool invert) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Node> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & NODE) && !m->GetMarker(lc[it],hm) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & NODE) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		return ret;
	}

	ElementArray<Edge> ElementSet::getEdges() const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Edge> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & EDGE) && !m->GetMarker(lc[it],hm) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)//check for type filters invalid handles
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & EDGE) )
					ret.push_back(lc[it]);
		}
		return ret;
	}

	ElementArray<Edge> ElementSet::getEdges(MarkerType select, bool invert) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Edge> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & EDGE) && !m->GetMarker(lc[it],hm) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & EDGE) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		return ret;
	}

	ElementArray<Face> ElementSet::getFaces() const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Face> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & FACE) && !m->GetMarker(lc[it],hm) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)//check for type filters invalid handles
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & FACE) )
					ret.push_back(lc[it]);
		}
		return ret;
	}

	ElementArray<Face> ElementSet::getFaces(MarkerType select, bool invert) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Face> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & FACE) && !m->GetMarker(lc[it],hm) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & FACE) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		return ret;
	}

	ElementArray<Cell> ElementSet::getCells() const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Cell> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & CELL) && !m->GetMarker(lc[it],hm) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)//check for type filters invalid handles
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & CELL) )
					ret.push_back(lc[it]);
		}
		return ret;
	}

	ElementArray<Cell> ElementSet::getCells(MarkerType select, bool invert) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Cell> ret(m);
		Element::adj_type const & lc = m->LowConn(GetHandle());
		if( m->isMeshModified() )
		{
			MarkerType hm = m->HideMarker();
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & CELL) && !m->GetMarker(lc[it],hm) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && (GetHandleElementType(lc[it]) & CELL) && (invert ^ m->GetMarker(lc[it],select)) )//check for type filters invalid handles
					ret.push_back(lc[it]);
		}
		return ret;
	}

	void ElementSet::PutElement(HandleType e) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		switch(BulkDF(m->SetComparatorTag()))
		{
		case UNSORTED_COMPARATOR:
			{
				Element::adj_type & hc = m->HighConn(GetHandle());
				if( hc.size() > high_conn_reserved )
				{
					lc[static_cast<enumerator>(hc.back())] = e;
					hc.pop_back();
				}
				else lc.push_back(e);
			}
			break;
		case   HANDLE_COMPARATOR:
		case  HIERARCHY_COMPARATOR:
		case GLOBALID_COMPARATOR:
		case CENTROID_COMPARATOR:
			lc.push_back(e);
			break;
		}
	}

	void ElementSet::PutElements(const HandleType * h, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		switch(BulkDF(m->SetComparatorTag()))
		{
		case UNSORTED_COMPARATOR:
			{
				const HandleType * h2 = h;
				Element::adj_type & hc = m->HighConn(GetHandle());
				while( hc.size() > high_conn_reserved )
				{
					lc[static_cast<enumerator>(hc.back())] = *h2;
					hc.pop_back();
					h2++;
				}
				lc.insert(lc.end(),h2, h+num);
			}
			break;
		case   HANDLE_COMPARATOR:
		case HIERARCHY_COMPARATOR:
		case GLOBALID_COMPARATOR:
		case CENTROID_COMPARATOR:
			lc.insert(lc.end(),h,h+num);
			break;
		}
	}

	void ElementSet::AddElement(HandleType e) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		Element::adj_type::iterator find;
		Element::adj_type::size_type sorted = static_cast<Element::adj_type::size_type>(hSorted(hc));
		switch(BulkDF(m->SetComparatorTag()))
		{
		case UNSORTED_COMPARATOR:
			for(Element::adj_type::iterator it = lc.begin(); it != lc.end(); ++it)
				if( *it == e ) return;
			if( hc.size() > high_conn_reserved )
			{
				lc[static_cast<enumerator>(hc.back())] = e;
				hc.pop_back();
			}
			else lc.push_back(e);
			return;
		case   HANDLE_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e);
			break;
		case  HIERARCHY_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e,Mesh::HierarchyComparator(m));
			break;
		case GLOBALID_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e,Mesh::GlobalIDComparator(m));
			break;
		case CENTROID_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e,Mesh::CentroidComparator(m));
			break;
		}
		if( find != lc.begin()+sorted )
		{
			if( *find == e )
				return;
		}
		lc.insert(find,1,e);
		hSorted(hc)++;
	}

	void ElementSet::AddElements(const HandleType * h, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		switch(BulkDF(m->SetComparatorTag()))
		{
		case UNSORTED_COMPARATOR:
			{
				Element::adj_type & hc = m->HighConn(GetHandle());
				MarkerType mrk = m->CreatePrivateMarker();
				SetPrivateMarkerElements(mrk);
				enumerator it = 0;
				while( it != num )
				{
					if( !m->GetPrivateMarker(h[it],mrk) )
					{
						if( hc.size() > high_conn_reserved )
						{
							lc[static_cast<enumerator>(hc.back())] = h[it];
							hc.pop_back();
						}
						else lc.push_back(h[it]);
					}
					++it;
				}
				RemPrivateMarkerElements(mrk);
				m->ReleasePrivateMarker(mrk);
			}
			break;
		case HANDLE_COMPARATOR:
			{
				// sort 
				array<HandleType> sorted(h,h+num);
				m->SortHandles(sorted.data(),sorted.size());
				//now merge with sorted part
				Element::adj_type merge(lc.size()+sorted.size());
				Element::adj_type::iterator it = lc.begin(), iend = lc.begin()+hSorted(hc);
				Element::adj_type::iterator jt = sorted.begin(), jend = sorted.end();
				Element::adj_type::iterator mit = merge.begin(), mitprev = merge.begin();
				//insert first element
				while(it != iend && jt != jend)
				{
					if( *it == InvalidHandle() ) ++it;
					else if( *jt == InvalidHandle() ) ++jt;
					else if( *it < *jt ) {*mit++ = *it++; break;}
					else {*mit++ = *jt++; break;}
				}
				//now insertion checks for duplicates
				while(it != iend && jt != jend)
				{
					if( *it == InvalidHandle() || *it == *mitprev ) ++it;
					else if( *jt == InvalidHandle() || *jt == *mitprev ) ++jt;
					else if( *it < *jt ) {mitprev = mit; *mit = *it; ++mit; ++it;}
					else {mitprev = mit; *mit = *jt; ++mit; ++jt;}
				}
				//attach rest
				while(it != iend && *it != InvalidHandle() && *it != *mitprev )
				{
					mitprev = mit;
					*mit = *it;
					++mit;
					++it;
				}
				while(jt != jend && *jt != InvalidHandle() && *jt != *mitprev )
				{
					mitprev = mit;
					*mit = *jt;
					++mit;
					++jt;
				}
				//indicate end of sorted part
				hSorted(hc) = static_cast<HandleType>(mit-merge.begin());
				//attach unsorted part with control for duplicates
				MarkerType mrk = m->CreatePrivateMarker();
				m->SetPrivateMarkerArray(merge.data(),hSorted(hc),mrk);
				while(it != lc.end() && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
				{
					m->SetPrivateMarker(*it,mrk);
					*mit++ = *it++;
				}
				merge.resize(static_cast<Element::adj_type::size_type>(mit-merge.begin()));
				m->RemPrivateMarkerArray(merge.data(),static_cast<enumerator>(merge.size()),mrk);
				m->ReleasePrivateMarker(mrk);
				// set merged elemetns as current
				lc.swap(merge);
				//indicate that there are no empty spaces
				hc.resize(high_conn_reserved);
			}
			break;
		case GLOBALID_COMPARATOR:
			{
				// sort 
				array<HandleType> sorted(h,h+num);
				m->SortByGlobalID(sorted.data(),sorted.size());
				//now merge with sorted part
				Element::adj_type merge(lc.size()+sorted.size());
				Element::adj_type::iterator it = lc.begin(), iend = lc.begin()+hSorted(hc);
				Element::adj_type::iterator jt = sorted.begin(), jend = sorted.end();
				Element::adj_type::iterator mit = merge.begin();
				MarkerType mrk = m->CreatePrivateMarker();
        //now insertion checks for duplicates
				while(it != iend && jt != jend)
				{
					if( *it == InvalidHandle() || m->GetPrivateMarker(*it,mrk) ) ++it;
					else if( *jt == InvalidHandle() || m->GetPrivateMarker(*jt,mrk) ) ++jt;
					else if( Mesh::GlobalIDComparator(m)(*it,*jt) ) 
          {
            *mit = *it; 
            m->SetPrivateMarker(*mit,mrk);
            ++mit; 
            ++it;
          }
					else 
          {
            *mit = *jt;
            m->SetPrivateMarker(*mit,mrk);
            ++mit; 
            ++jt;
          }
				}
				//attach rest
				while(it != iend && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
				{
					*mit = *it;
          m->SetPrivateMarker(*mit,mrk);
					++mit;
					++it;
				}
				while(jt != jend && *jt != InvalidHandle() && !m->GetPrivateMarker(*jt,mrk) )
				{
					*mit = *jt;
          m->SetPrivateMarker(*mit,mrk);
					++mit;
					++jt;
				}
				//indicate end of sorted part
				hSorted(hc) = static_cast<HandleType>(mit-merge.begin());
				//attach unsorted part
				while(it != lc.end() && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
				{
					m->SetPrivateMarker(*it,mrk);
					*mit++ = *it++;
				}
				merge.resize(static_cast<Element::adj_type::size_type>(mit-merge.begin()));
				m->RemPrivateMarkerArray(merge.data(),static_cast<enumerator>(merge.size()),mrk);
				m->ReleasePrivateMarker(mrk);
				// set merged elemetns as current
				lc.swap(merge);
				//indicate that there are no empty spaces
				hc.resize(high_conn_reserved);
			}
			break;
		case HIERARCHY_COMPARATOR:
			{
				// sort 
				array<HandleType> sorted(h,h+num);
				std::sort(sorted.begin(),sorted.end(),Mesh::HierarchyComparator(m));
				//now merge with sorted part
				Element::adj_type merge(lc.size()+sorted.size());
				Element::adj_type::iterator it = lc.begin(), iend = lc.begin()+hSorted(hc);
				Element::adj_type::iterator jt = sorted.begin(), jend = sorted.end();
				Element::adj_type::iterator mit = merge.begin(), mitprev = merge.begin();
				MarkerType mrk = m->CreatePrivateMarker();
				//now insertion checks for duplicates
				while(it != iend && jt != jend)
				{
					if( *it == InvalidHandle() || m->GetPrivateMarker(*it,mrk) ) ++it;
					else if( *jt == InvalidHandle() || m->GetPrivateMarker(*jt,mrk) ) ++jt;
					else if( Mesh::HierarchyComparator(m)(*it,*jt) ) 
          {
            *mit = *it; 
            m->SetPrivateMarker(*mit,mrk);
            ++mit; ++it;
          }
					else 
          {
            *mit = *jt; 
            m->SetPrivateMarker(*mit,mrk);
            ++mit; ++jt;
          }
				}
				//attach rest
				while(it != iend && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
				{
					*mit = *it;
          m->SetPrivateMarker(*mit,mrk);
					++mit;
					++it;
				}
				while(jt != jend && *jt != InvalidHandle() && !m->GetPrivateMarker(*jt,mrk) )
				{
					*mit = *jt;
          m->SetPrivateMarker(*mit,mrk);
					++mit;
					++jt;
				}
				//indicate end of sorted part
				hSorted(hc) = static_cast<HandleType>(mit-merge.begin());
				//attach unsorted part
				while(it != lc.end() && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
				{
					m->SetPrivateMarker(*it,mrk);
					*mit++ = *it++;
				}
				merge.resize(static_cast<Element::adj_type::size_type>(mit-merge.begin()));
				m->RemPrivateMarkerArray(merge.data(),static_cast<enumerator>(merge.size()),mrk);
				m->ReleasePrivateMarker(mrk);
				// set merged elemetns as current
				lc.swap(merge);
				//indicate that there are no empty spaces
				hc.resize(high_conn_reserved);
			}
			break;
		case CENTROID_COMPARATOR:
			{
				// sort 
				array<HandleType> sorted(h,h+num);
				std::sort(sorted.begin(),sorted.end(),Mesh::CentroidComparator(m));
				//now merge with sorted part
				Element::adj_type merge(lc.size()+sorted.size());
				Element::adj_type::iterator it = lc.begin(), iend = lc.begin()+hSorted(hc);
				Element::adj_type::iterator jt = sorted.begin(), jend = sorted.end();
				Element::adj_type::iterator mit = merge.begin(), mitprev = merge.begin();
				MarkerType mrk = m->CreatePrivateMarker();
				//now insertion checks for duplicates
				while(it != iend && jt != jend)
				{
					if( *it == InvalidHandle() || m->GetPrivateMarker(*it,mrk) ) ++it;
					else if( *jt == InvalidHandle() || m->GetPrivateMarker(*jt,mrk) ) ++jt;
					else if( Mesh::CentroidComparator(m)(*it,*jt) ) 
          {
            *mit = *it; 
            m->SetPrivateMarker(*mit,mrk);
            ++mit; ++it;
          }
					else 
          {
            *mit = *jt; 
            m->SetPrivateMarker(*mit,mrk);
            ++mit; ++it;
          }
				}
				//attach rest
				while(it != iend && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
				{
					*mit = *it;
          m->SetPrivateMarker(*mit,mrk);
					++mit;
					++it;
				}
				while(jt != jend && *jt != InvalidHandle() && !m->GetPrivateMarker(*jt,mrk) )
				{
					*mit = *jt;
          m->SetPrivateMarker(*mit,mrk);
					++mit;
					++jt;
				}
				//indicate end of sorted part
				hSorted(hc) = static_cast<HandleType>(mit-merge.begin());
				//attach unsorted part
				while(it != lc.end() && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
				{
					m->SetPrivateMarker(*it,mrk);
					*mit++ = *it++;
				}
				merge.resize(static_cast<Element::adj_type::size_type>(mit-merge.begin()));
				m->RemPrivateMarkerArray(merge.data(),static_cast<enumerator>(merge.size()),mrk);
				m->ReleasePrivateMarker(mrk);
				// set merged elemetns as current
				lc.swap(merge);
				//indicate that there are no empty spaces
				hc.resize(high_conn_reserved);
			}
			break;
		}
	}

	void ElementSet::RemoveElement(const Storage & e) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		Element::adj_type::iterator find;
		Element::adj_type::size_type i, sorted = static_cast<Element::adj_type::size_type>(hSorted(hc));
		switch(BulkDF(m->SetComparatorTag()))
		{
		case UNSORTED_COMPARATOR:
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] == e->GetHandle() )
				{
					lc[it] = InvalidHandle();
					hc.push_back(static_cast<HandleType>(it));
					break;
				}
			return;
		case   HANDLE_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e->GetHandle());
			break;
		case  HIERARCHY_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e->GetHandle(),Mesh::HierarchyComparator(m));
			break;
		case GLOBALID_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e->GetHandle(),Mesh::GlobalIDComparator(m));
			break;
		case CENTROID_COMPARATOR:
			find = std::lower_bound(lc.begin(),lc.begin()+sorted,e->GetHandle(),Mesh::CentroidComparator(m));
			break;
		}
		if( find != lc.begin()+sorted )
		{
			if( *find == e->GetHandle() )
			{
				*find = InvalidHandle();
				hc.push_back(static_cast<HandleType>(find-lc.begin()));
			}
		}
		else for(i = sorted; i < lc.size(); ++i) //unsorted part
			if( lc[i] == e->GetHandle() )
			{
				lc[i] = InvalidHandle();
				hc.push_back(static_cast<HandleType>(i));
			}
	}

	void ElementSet::RemoveElements(const HandleType * h, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		Element::adj_type::size_type i;
		MarkerType mrk = m->CreatePrivateMarker();
		m->SetPrivateMarkerArray(h,num,mrk);
		for(i = 0; i < lc.size(); i++) if( lc[i] != InvalidHandle() )
		{
			if( m->GetPrivateMarker(lc[i],mrk) )
			{
				lc[i] = InvalidHandle();
				hc.push_back(i);
			}
		}
		m->RemPrivateMarkerArray(h,num,mrk);
		m->ReleasePrivateMarker(mrk);
	}

	ElementArray<Element> ElementSet::Union(const ElementSet & other) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type & lc = m->LowConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && !m->GetPrivateMarker(lc[i],mrk))
			{
				ret.push_back(lc[i]);
				m->SetPrivateMarker(ret.atback(),mrk);
			}
		}
		Element::adj_type & other_lc = m->LowConn(other->GetHandle());
		for(const HandleType * h = other_lc.data(); h != other_lc.data() + other_lc.size(); ++h)
		{
			if( *h != InvalidHandle() && !m->GetPrivateMarker(*h,mrk) )
			{
				ret.push_back(*h);
				m->SetPrivateMarker(ret.atback(),mrk);
			}
		}
		ret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return ret;
	}

	ElementArray<Element> ElementSet::Union(const HandleType * handles, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type & lc = m->LowConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && !m->GetPrivateMarker(lc[i],mrk))
			{
				ret.push_back(lc[i]);
				m->SetPrivateMarker(ret.atback(),mrk);
			}
		}
		for(const HandleType * h = handles; h != handles+num; ++h)
		{
			if( *h != InvalidHandle() && !m->GetPrivateMarker(*h,mrk) )
			{
				ret.push_back(*h);
				m->SetPrivateMarker(ret.atback(),mrk);
			}
		}
		ret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return ret;
	}

	ElementArray<Element> ElementSet::Difference(const ElementSet & other) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type & lc = m->LowConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		other->SetPrivateMarkerElements(mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && !m->GetPrivateMarker(lc[i],mrk))
			{
				ret.push_back(lc[i]);
				m->SetPrivateMarker(ret.atback(),mrk);
			}
		}
		other->RemPrivateMarkerElements(mrk);
		ret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return ret;
	}

	ElementArray<Element> ElementSet::Difference(const HandleType * handles, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type & lc = m->LowConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		m->SetPrivateMarkerArray(handles,num,mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && !m->GetPrivateMarker(lc[i],mrk))
			{
				ret.push_back(lc[i]);
				m->SetPrivateMarker(ret.atback(),mrk);
			}
		}
		m->RemPrivateMarkerArray(handles,num,mrk);
		ret.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return ret;
	}

	ElementArray<Element> ElementSet::Intersection(const ElementSet & other) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type & lc = m->LowConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker(), mrk2 = m->CreatePrivateMarker();
		other->SetPrivateMarkerElements(mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && m->GetPrivateMarker(lc[i],mrk) && !m->GetPrivateMarker(lc[i],mrk2))
			{
				ret.push_back(lc[i]);
				m->SetPrivateMarker(lc[i],mrk2);
			}
		}
		other->RemPrivateMarkerElements(mrk);
		ret.RemPrivateMarker(mrk2);
		m->ReleasePrivateMarker(mrk2);
		m->ReleasePrivateMarker(mrk);
		return ret;
	}

	ElementArray<Element> ElementSet::Intersection(const HandleType * handles, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		ElementArray<Element> ret(m);
		Element::adj_type & lc = m->LowConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker(), mrk2 = m->CreatePrivateMarker();
		m->SetPrivateMarkerArray(handles,num,mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && m->GetPrivateMarker(lc[i],mrk) && !m->GetPrivateMarker(lc[i],mrk2))
			{
				ret.push_back(lc[i]);
				m->SetPrivateMarker(lc[i],mrk2);
			}
		}
		m->RemPrivateMarkerArray(handles,num,mrk);
		ret.RemPrivateMarker(mrk2);
		m->ReleasePrivateMarker(mrk2);
		m->ReleasePrivateMarker(mrk);
		return ret;
	}


	void ElementSet::Unite(const ElementSet & other) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		bulk cmp = BulkDF(m->SetComparatorTag());
		Element::adj_type & other_lc = m->LowConn(other->GetHandle());
		if( cmp == UNSORTED_COMPARATOR || other->BulkDF(m->SetComparatorTag()) != cmp )
			Unite(other_lc.data(),static_cast<enumerator>(other_lc.size()));
		else //other set sorted by the same comparator, just merge them
		{
			Element::adj_type & other_hc = m->HighConn(other->GetHandle());
			Element::adj_type & hc = m->HighConn(GetHandle());
			Element::adj_type & lc = m->LowConn(GetHandle());
			Element::adj_type merge(lc.size()+other_lc.size());
			Element::adj_type::iterator it = lc.begin(), iend = lc.begin()+hSorted(hc);
			Element::adj_type::iterator jt = other_lc.begin(), jend = other_lc.begin()+hSorted(other_hc);
			Element::adj_type::iterator mit = merge.begin(), mitprev = merge.begin();
			MarkerType mrk = m->CreatePrivateMarker();
			//now insertion checks for duplicates
			while(it != iend && jt != jend)
			{
				if( *it == InvalidHandle() || m->GetPrivateMarker(*it,mrk) ) ++it;
				else if( *jt == InvalidHandle() || m->GetPrivateMarker(*jt,mrk) ) ++jt;
				else 
				{
					bool cmp_res = false;
					switch(cmp)
					{
					case HANDLE_COMPARATOR: cmp_res = *it < *jt; break;
					case HIERARCHY_COMPARATOR: cmp_res = Mesh::HierarchyComparator(m)(*it,*jt); break;
					case CENTROID_COMPARATOR: cmp_res = Mesh::CentroidComparator(m)(*it,*jt); break;
					case GLOBALID_COMPARATOR: cmp_res = Mesh::GlobalIDComparator(m)(*it,*jt); break;
					default: assert(false);
					}
					if( cmp_res ) 
					{
						*mit = *it;
            m->SetPrivateMarker(*mit,mrk);
						++mit; 
            ++it;
					}
					else 
					{
						*mit = *jt; 
            m->SetPrivateMarker(*mit,mrk);
						++mit; 
						++jt;
					}
				}
			}
			//attach rest one of the sets not reached it's end
			while(it != iend && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
			{
				*mit = *it;
        m->SetPrivateMarker(*mit,mrk);
				++mit;
				++it;
			}
			while(jt != jend && *jt != InvalidHandle() && !m->GetPrivateMarker(*jt,mrk) )
			{
				*mit = *jt;
        m->SetPrivateMarker(*mit,mrk);
				++mit;
				++jt;
			}
			//indicate end of sorted part
			hSorted(hc) = static_cast<HandleType>(mit-merge.begin());
			//attach unsorted part
			while(it != lc.end() && *it != InvalidHandle() && !m->GetPrivateMarker(*it,mrk) )
			{
				m->SetPrivateMarker(*it,mrk);
				*mit++ = *it++;
			}
			//attach unsorted part
			while(jt != other_lc.end() && *jt != InvalidHandle() && !m->GetPrivateMarker(*jt,mrk) )
			{
				m->SetPrivateMarker(*jt,mrk);
				*mit++ = *jt++;
			}
			merge.resize(static_cast<Element::adj_type::size_type>(mit-merge.begin()));
			m->RemPrivateMarkerArray(merge.data(),static_cast<enumerator>(merge.size()),mrk);
			m->ReleasePrivateMarker(mrk);
			// set merged elemetns as current
			lc.swap(merge);
			//indicate that there are no empty spaces
			hc.resize(high_conn_reserved);
		}
	}

	void ElementSet::Unite(const HandleType * handles, enumerator num) const
	{
		assert(GetElementType() == ESET);
		AddElements(handles,num); //These algorithms are identical
	}


	void ElementSet::Subtract(const ElementSet & other) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		other->SetPrivateMarkerElements(mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && m->GetPrivateMarker(lc[i],mrk))
			{
				lc[i] = InvalidHandle();
				hc.push_back(i);
			}
		}
		other->RemPrivateMarkerElements(mrk);
		m->ReleasePrivateMarker(mrk);
	}

	void ElementSet::Subtract(const HandleType * handles, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		m->SetPrivateMarkerArray(handles,num,mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && m->GetPrivateMarker(lc[i],mrk))
			{
				lc[i] = InvalidHandle();
				hc.push_back(i);
			}
		}
		m->RemPrivateMarkerArray(handles,num,mrk);
		m->ReleasePrivateMarker(mrk);
	}

	void ElementSet::Intersect(const ElementSet & other) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		other->SetPrivateMarkerElements(mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && !m->GetPrivateMarker(lc[i],mrk) )
			{
				lc[i] = InvalidHandle();
				hc.push_back(i);
			}
		}
		other->RemPrivateMarkerElements(mrk);
		m->ReleasePrivateMarker(mrk);
	}

	void ElementSet::Intersect(const HandleType * handles, enumerator num) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		MarkerType mrk = m->CreatePrivateMarker();
		m->SetPrivateMarkerArray(handles,num,mrk);
		for(Element::adj_type::size_type i = 0; i < lc.size(); ++i)
		{
			if( lc[i] != InvalidHandle() && !m->GetPrivateMarker(lc[i],mrk) )
			{
				lc[i] = InvalidHandle();
				hc.push_back(i);
			}
		}
		m->RemPrivateMarkerArray(handles,num,mrk);
		m->ReleasePrivateMarker(mrk);
	}
	
	//void ElementSet::SetExchange(ExchangeType comp) const
	//{
	//	assert(GetElementType() == ESET);
	//	Bulk(GetMeshLink()->SetExchangeTag()) = comp;
	//}

	void ElementSet::SortSet(ComparatorType comp) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type * lc = NULL;
		Element::adj_type & hc = m->HighConn(GetHandle());
		if( comp != UNSORTED_COMPARATOR )
			lc = &m->LowConn(GetHandle());
		bulk & comparator = BulkDF(m->SetComparatorTag());
		switch(comp)
		{
		case GLOBALID_COMPARATOR:
			{
				Element::adj_type::size_type start = hSorted(hc);
				Element::adj_type::size_type i, j = 0, k = 0, num = lc->size() - start;
				if( hc.size() > high_conn_reserved ) // there are invalid handles
				{
					//move invalid handles to the end
					for(i = start; i < lc->size(); i++)
					{
						if( lc->at(i) != InvalidHandle() )
						{
							lc->at(j++) = lc->at(i);
							lc->at(i) = InvalidHandle();
						}
						else k++;
					}
					num = lc->size() - start - k; //size of array to sort
				}
				m->SortByGlobalID(lc->data()+start,static_cast<enumerator>(num));
				if( start > 0 ) //need to merge two parts
				{
					Element::adj_type temp(lc->size());
					Element::adj_type::iterator it = lc->begin(), jt = lc->begin()+start;
					Element::adj_type::iterator iend = jt, jend = lc->end();
					Element::adj_type::iterator t = temp.begin();
					while(it != iend && jt != jend)
					{
						if( *it == InvalidHandle() ) ++it;
						else if( *jt == InvalidHandle() ) ++jt;
						else if( Mesh::GlobalIDComparator(m)(*it,*jt) ) *t++ = *it++;
						else *t++ = *jt++;
					}
					assert(it == iend || jt == jend);
					while(it != iend && *it != InvalidHandle() ) *t++ = *it++;
					while(jt != jend && *jt != InvalidHandle() ) *t++ = *jt++;
					//number of deleted items and invalid handles should match
					assert( (temp.size() - (t-temp.begin())) == static_cast<ptrdiff_t>(hc.size()-ElementSet::high_conn_reserved) );
					lc->swap(temp);
				}
			}
			break;
		case CENTROID_COMPARATOR:
			{
				Element::adj_type::size_type start = hSorted(hc);
				std::sort(lc->begin()+start,lc->end(),Mesh::CentroidComparator(m));
				if( start > 0 ) //need to merge two parts
				{
					Element::adj_type temp(lc->size());
					Element::adj_type::iterator it = lc->begin(), jt = lc->begin()+start;
					Element::adj_type::iterator iend = jt, jend = lc->end();
					Element::adj_type::iterator t = temp.begin();
					while(it != iend && jt != jend)
					{
						if( *it == InvalidHandle() ) ++it;
						else if( *jt == InvalidHandle() ) ++jt;
						else if( Mesh::CentroidComparator(m)(*it,*jt) ) *t++ = *it++;
						else *t++ = *jt++;
					}
					assert(it == iend || jt == jend);
					while(it != iend && *it != InvalidHandle() ) *t++ = *it++;
					while(jt != jend && *jt != InvalidHandle() ) *t++ = *jt++;
					//number of deleted items and invalid handles should match
					assert( (temp.size() - (t-temp.begin())) == static_cast<ptrdiff_t>(hc.size()-ElementSet::high_conn_reserved) );
					lc->swap(temp);
				}
			}
			break;
		case  HIERARCHY_COMPARATOR:
			{
				Element::adj_type::size_type start = hSorted(hc);
				std::sort(lc->begin()+start,lc->end(),Mesh::HierarchyComparator(m));
				if( start > 0 ) //need to merge two parts
				{
					Element::adj_type temp(lc->size());
					Element::adj_type::iterator it = lc->begin(), jt = lc->begin()+start;
					Element::adj_type::iterator iend = jt, jend = lc->end();
					Element::adj_type::iterator t = temp.begin();
					while(it != iend && jt != jend)
					{
						if( *it == InvalidHandle() ) ++it;
						else if( *jt == InvalidHandle() ) ++jt;
						else if( Mesh::HierarchyComparator(m)(*it,*jt) ) *t++ = *it++;
						else *t++ = *jt++;
					}
					assert(it == iend || jt == jend);
					while(it != iend && *it != InvalidHandle() ) *t++ = *it++;
					while(jt != jend && *jt != InvalidHandle() ) *t++ = *jt++;
					//number of deleted items and invalid handles should match
					assert( (temp.size() - (t-temp.begin())) == static_cast<ptrdiff_t>(hc.size()-ElementSet::high_conn_reserved) );
					lc->swap(temp);
				}
			}
			break;
		case HANDLE_COMPARATOR:
			{
				Element::adj_type::size_type start = hSorted(hc);
				m->SortHandles(lc->data()+start,static_cast<enumerator>(lc->size()-start));
				if( start > 0 ) //need to merge two parts
				{
					Element::adj_type temp(lc->size());
					Element::adj_type::iterator it = lc->begin(), jt = lc->begin()+start;
					Element::adj_type::iterator iend = jt, jend = lc->end();
					Element::adj_type::iterator t = temp.begin();
					while(it != iend && jt != jend)
					{
						if( *it == InvalidHandle() ) ++it;
						else if( *jt == InvalidHandle() ) ++jt;
						else if( *it < *jt ) *t++ = *it++;
						else *t++ = *jt++;
					}
					assert(it == iend || jt == jend);
					while(it != iend && *it != InvalidHandle() ) *t++ = *it++;
					while(jt != jend && *jt != InvalidHandle() ) *t++ = *jt++;
					//number of deleted items and invalid handles should match
					assert( (temp.size() - (t-temp.begin())) == static_cast<ptrdiff_t>(hc.size()-ElementSet::high_conn_reserved) );
					lc->swap(temp);
				}
			}
			break;
		case UNSORTED_COMPARATOR:
			break;
		}
		if( comp != UNSORTED_COMPARATOR )
		{
			///all invalid handles moved to the end, cut them
			lc->resize(lc->size() - (hc.size()-ElementSet::high_conn_reserved) );
			hc.resize(ElementSet::high_conn_reserved);
			hSorted(hc) = lc->size();
		}
		else hSorted(hc) = 0;
		comparator = comp;
	}
	ElementSet::ComparatorType ElementSet::GetComparator() const
	{
		return BulkDF(GetMeshLink()->SetComparatorTag());
	}
	//ElementSet::ExchangeType ElementSet::GetExchange() const
	//{
	//	return BulkDF(GetMeshLink()->SetExchangeTag());
	//}
	Element ElementSet::FindElementByGlobalID(integer global_id) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		ComparatorType cmp = BulkDF(m->SetComparatorTag());
		if( cmp == GLOBALID_COMPARATOR )
		{
			Element::adj_type::size_type sorted = static_cast<Element::adj_type::size_type>(hSorted(m->HighConn(GetHandle())));
			Element::adj_type::iterator find = std::lower_bound(lc.begin(),lc.begin()+sorted,global_id,Mesh::GlobalIDComparator(m));
			if( find != lc.begin()+sorted ) 
			{
				if(*find != InvalidHandle() && m->GlobalID(*find) == global_id)
					return Element(m,*find);
			}
			else for(Element::adj_type::size_type k = sorted; k < lc.size(); k++)
				if( lc[k] != InvalidHandle() && m->GlobalID(lc[k]) == global_id )
					return Element(m,lc[k]);
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() && m->GlobalID(lc[it]) == global_id )
					return Element(m,lc[it]);
		}
		return Element(m,InvalidHandle());
	}

	Element ElementSet::FindElementByCentroid(real * centroid) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		ComparatorType cmp = BulkDF(m->SetComparatorTag());
		if( cmp == CENTROID_COMPARATOR )
		{
			Element::adj_type::size_type sorted = static_cast<Element::adj_type::size_type>(hSorted(m->HighConn(GetHandle())));
			Element::adj_type::iterator find = std::lower_bound(lc.begin(),lc.begin()+sorted,centroid,Mesh::CentroidComparator(m));
			if( find != lc.begin()+sorted )
			{
				Storage::real cnt[3];
				m->GetGeometricData(*find,CENTROID,cnt);
				if( *find != InvalidHandle() &&  Mesh::CentroidComparator(m).Compare(cnt,centroid) == 0 )
					return Element(m,*find);
			}
			else for(Element::adj_type::size_type k = sorted; k < lc.size(); k++)
				if( lc[k] != InvalidHandle() )
				{
					Storage::real cnt[3];
					m->GetGeometricData(lc[k],CENTROID,cnt);
					if( Mesh::CentroidComparator(m).Compare(cnt,centroid) == 0 )
						return Element(m,lc[k]);
				}
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] != InvalidHandle() )
				{
					Storage::real cnt[3];
					m->GetGeometricData(lc[it],CENTROID,cnt);
					if( Mesh::CentroidComparator(m).Compare(cnt,centroid) == 0 )
						return Element(m,lc[it]);
				}
		}
		return Element(m,InvalidHandle());
	}

	bool ElementSet::FindHandle(HandleType h, bool use_cmp) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		ComparatorType cmp = BulkDF(m->SetComparatorTag());
		if( cmp == HANDLE_COMPARATOR  && use_cmp)
		{
			Element::adj_type::size_type sorted = static_cast<Element::adj_type::size_type>(hSorted(m->HighConn(GetHandle())));
			Element::adj_type::iterator find = std::lower_bound(lc.begin(),lc.begin()+sorted,h);
			if( find != lc.begin()+sorted )
			{
				if( *find == h )
					return true;
			}
			else for(Element::adj_type::size_type i = sorted; i < lc.size(); i++)
				if( lc[i] == h ) return true;
		}
		else if( cmp == GLOBALID_COMPARATOR && use_cmp)
		{
			integer gid = m->GlobalID(h);
			Element::adj_type::size_type sorted = static_cast<Element::adj_type::size_type>(hSorted(m->HighConn(GetHandle())));
			Element::adj_type::iterator find = std::lower_bound(lc.begin(),lc.begin()+sorted,gid,Mesh::GlobalIDComparator(m));
			if( find != lc.begin()+sorted )
			{
				if( *find == h )
					return true;
			}
			else for(Element::adj_type::size_type i = sorted; i < lc.size(); i++)
				if( lc[i] == h ) return true;
		}
		else if( cmp == CENTROID_COMPARATOR && use_cmp )
		{
			Storage::real cnt[3];
			m->GetGeometricData(h,cmp,cnt);
			Element::adj_type::size_type sorted = static_cast<Element::adj_type::size_type>(hSorted(m->HighConn(GetHandle())));
			Element::adj_type::iterator find = std::lower_bound(lc.begin(),lc.begin()+sorted,cnt,Mesh::CentroidComparator(m));
			if( find != lc.begin()+sorted )
			{
				if( *find == h )
					return true;
			}
			else for(Element::adj_type::size_type i = sorted; i < lc.size(); i++)
				if( lc[i] == h ) return true;
		}
		else if( cmp == HIERARCHY_COMPARATOR && use_cmp )
		{
			Element::adj_type::size_type sorted = static_cast<Element::adj_type::size_type>(hSorted(m->HighConn(GetHandle())));
			Element::adj_type::iterator find = std::lower_bound(lc.begin(),lc.end(),h,Mesh::HierarchyComparator(m));
			if( find != lc.end() )
			{
				if( *find == h )
					return true;
			}
			else for(Element::adj_type::size_type i = sorted; i < lc.size(); i++)
				if( lc[i] == h ) return true;
		}
		else
		{
			for(Element::adj_type::size_type it = 0; it < lc.size(); ++it)
				if( lc[it] == h )
					return true;
		}
		return false;
	}

	ElementSet::iterator ElementSet::Begin() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & lc = GetMeshLink()->LowConn(GetHandle());
		return ++iterator(GetMeshLink(),&lc,static_cast<Element::adj_type::size_type>(-1));
	}

	ElementSet::iterator ElementSet::End() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & lc = GetMeshLink()->LowConn(GetHandle());
		return iterator(GetMeshLink(),&lc,lc.size());
	}

	ElementSet::iterator ElementSet::EndSorted() const
	{
		assert(GetElementType() == ESET);
		Element::adj_type & hc = GetMeshLink()->HighConn(GetHandle());
		Element::adj_type & lc = GetMeshLink()->LowConn(GetHandle());
		return ++iterator(GetMeshLink(),&lc,hSorted(hc)-1);
	}

	ElementSet::iterator ElementSet::Erase(iterator pos) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & hc = m->HighConn(GetHandle());
		Element::adj_type & lc = m->LowConn(GetHandle());
		assert(pos.pos < lc.size()); //cannot point over the end
		assert(&lc == pos.ptr);
		assert(lc[pos.pos] != InvalidHandle()); //cannot point to invalid handles
		//assert( !m->isMeshModified() || !m->GetMarker(lc[pos.pos],m->HideMarker())); //cannot point to hidden markers
		lc[pos.pos] = InvalidHandle();
		hc.push_back(pos.pos);
		return ++pos;
	}

	void ElementSet::Erase(iterator beg, iterator end) const
	{
		assert(GetElementType() == ESET);
		Mesh * m = GetMeshLink();
		Element::adj_type & hc = m->HighConn(GetHandle());
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type::size_type ibeg = beg.pos, iend = end.pos, i;
		assert(&lc == beg.ptr);
		assert(&lc == end.ptr);
		assert(ibeg <= lc.size()); //begin cannot point over the end
		assert(iend <= lc.size()); //end cannot point over the end
		for(i = ibeg; i < iend; ++i)
		{
			if( lc[i] != InvalidHandle() && (!m->isMeshModified() || !m->GetMarker(lc[i],m->HideMarker())) ) //delete non-invalid unhidden elements in the range
			{
				lc[i] = InvalidHandle();
				hc.push_back(i);
			}
		}
	}

	void ElementSet::ReorderEmpty() const
	{
		Mesh * m = GetMeshLink();
		Element::adj_type & hc = m->HighConn(GetHandle());
		Element::adj_type & lc = m->LowConn(GetHandle());
		if( BulkDF(m->SetComparatorTag()) == UNSORTED_COMPARATOR  ) //just compact the set
		{
			// This should be O(m) where m is the number of deleted items
			integer cend = static_cast<integer>(lc.size()-1);
			while( cend >= 0 && lc[cend] == InvalidHandle() ) --cend; //find first good handle
			while( cend >= 0 && hc.size() > high_conn_reserved )
			{
				if( static_cast<integer>(hc.back()) >= cend )
					hc.pop_back();
				else if( lc[cend] != InvalidHandle() )
				{
					lc[hc.back()] = lc[cend];
					lc[cend] = InvalidHandle();
				}
				else cend--;
			}
			while( cend >= 0 && lc[cend] == InvalidHandle() ) --cend; //skip bad handles
			lc.resize(cend);
			hc.resize(high_conn_reserved);
		}
		else
		{
			// This should be O(n+m) where n is the number of stored items and m is the number of deleted items
			// but preserves the initial order of elements
			Element::adj_type::size_type sorted = hSorted(hc);
			Element::adj_type::size_type i, j = 0;
			if( hc.size() > high_conn_reserved ) // there are invalid handles
			{
				//move invalid handles to the end
				for(i = 0; i < sorted; ++i)
				{
					if( lc[i] != InvalidHandle() )
					{
						lc[j++] = lc[i];
						lc[i] = InvalidHandle();
					}
				}
				hSorted(hc) = j; //update information about end of the part
				for(i = sorted; i < lc.size(); ++i) //now add the rest
				{
					if( lc[i] != InvalidHandle() )
					{
						lc[j++] = lc[i];
						lc[i] = InvalidHandle();
					}
				}
				lc.resize(j);
				hc.resize(high_conn_reserved);
			}

		}
	}

	bool ElementSet::Empty() const
	{
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		return (lc.size() - (hc.size() - high_conn_reserved)) == 0;
	}

	Storage::enumerator ElementSet::Size() const
	{
		Mesh * m = GetMeshLink();
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		return (lc.size() - (hc.size() - high_conn_reserved));
	}

	void ElementSet::Clear()
	{
		Mesh * m = GetMeshLink();
		if( HaveParent() ) GetParent()->RemChild(self());
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		hc.resize(ElementSet::high_conn_reserved);
		lc.clear();
		BulkDF(m->SetComparatorTag()) = UNSORTED_COMPARATOR;
	}

	bool ElementSet::DeleteSet()
	{
		Clear();
		return Delete();
	}

	bool ElementSet::DeleteSetTree()
	{
		Mesh * m = GetMeshLink();
		if( HaveParent() ) GetParent()->RemChild(self());
		Element::adj_type & lc = m->LowConn(GetHandle());
		Element::adj_type & hc = m->HighConn(GetHandle());
		HandleType hchild = hChild(hc), nchild;
		hChild(hc) = InvalidHandle();
		while(hchild != InvalidHandle())
		{
			Element::adj_type & chc = m->HighConn(hchild);
			hParent(chc) = InvalidHandle();
			nchild = hSibling(chc);
			hSibling(chc) = InvalidHandle();
			ElementSet(m,hchild).DeleteSetTree();
			hchild = nchild;
		}
		hc.resize(ElementSet::high_conn_reserved);
		lc.clear();
		BulkDF(m->SetComparatorTag()) = UNSORTED_COMPARATOR;
		return Delete();
	}
	
	void ElementSet::SynchronizeSetElements()
	{
		Mesh * m = GetMeshLink();
		if( m->GetMeshState() == Mesh::Serial ) return;
		
		//~ if( GetStatus() != Element::Owned )
		{
			Storage::integer_array set_procs = IntegerArray(m->ProcessorsTag());
			//Storage::integer_array elem_procs;
			//std::set<Storage::integer> send_set(set_procs.begin(),set_procs.end());
			//for(iterator it = Begin(); it != End(); ++it)
			//{
			//	elem_procs = m->IntegerArray(*it,m->ProcessorsTag());
			//	send_set.insert(elem_procs.begin(),elem_procs.end());
			//}
			//~ std::cout << "On " << m->GetProcessorRank() << " send set " << GetName() << " to ";
			//~ for(Storage::integer_array::iterator it = set_procs.begin(); it != set_procs.end(); ++it)
				//~ std::cout << *it << " ";
			//SendTo(set_procs);
			//~ std::cout << " with elements";
			//force send the set along with the elements
			Storage::integer_array send_to = IntegerArray(m->SendtoTag());
			std::set<Storage::integer> send_set(send_to.begin(),send_to.end());
			send_set.insert(set_procs.begin(),set_procs.end());
			send_set.erase(m->GetProcessorRank());
			send_to.replace(send_to.begin(),send_to.end(),send_set.begin(),send_set.end());
			for(iterator it = Begin(); it != End(); ++it)
			{
				//~ std::cout << " " << ElementTypeName(it->GetElementType()) << ":" << it->GlobalID();
				//it->SendTo(send_set);
				it->SendTo(set_procs);
			}
			//~ std::cout << std::endl;
		}
	}
	void ElementSet::SynchronizeSharedSetElements()
	{
		Mesh* m = GetMeshLink();
		if (m->GetMeshState() == Mesh::Serial) return;

		{
			Storage::integer_array set_procs = IntegerArray(m->ProcessorsTag());
			//force send the set along with the elements
			Storage::integer_array send_to = IntegerArray(m->SendtoTag());
			std::set<Storage::integer> send_set(send_to.begin(), send_to.end()), send_elem;
			for (iterator it = Begin(); it != End(); ++it) if(it->GetStatus() == Element::Shared)
			{
				Storage::integer_array elem_procs = it->IntegerArray(m->ProcessorsTag());
				for(Storage::integer_array::iterator jt = elem_procs.begin(); jt != elem_procs.end(); ++jt)
					if (*jt != m->GetProcessorRank())
					{
						send_elem.insert(*jt);
						send_set.insert(*jt);
					}
				it->SendTo(send_elem);
				send_elem.clear();
			}
			send_to.replace(send_to.begin(), send_to.end(), send_set.begin(), send_set.end());
			//~ std::cout << std::endl;
		}
	}
	void ElementSet::SynchronizeSetElementsWithOwner()
	{
		Mesh * m = GetMeshLink();
		if( m->GetMeshState() == Mesh::Serial ) return;
		
		if( GetStatus() != Element::Owned )
		{
			//Storage::integer_array set_procs = IntegerArray(m->ProcessorsTag());
			Storage::integer set_owner = Integer(m->OwnerTag());
			//Storage::integer_array elem_procs;
			//std::set<Storage::integer> send_set;//(set_procs.begin(),set_procs.end());
			//send_set.insert(set_owner);
			//for(iterator it = Begin(); it != End(); ++it)
			//{
			//	elem_procs = m->IntegerArray(*it,m->ProcessorsTag());
			//	send_set.insert(elem_procs.begin(),elem_procs.end());
			//}
			//force send the set along with the elements
			Storage::integer_array send_to = IntegerArray(m->SendtoTag());
			std::set<Storage::integer> send_set(send_to.begin(),send_to.end());
			send_set.insert(set_owner);
			send_to.replace(send_to.begin(),send_to.end(),send_set.begin(),send_set.end());
			//SendTo(send_set);
			for(iterator it = Begin(); it != End(); ++it)
				it->SendTo(send_set);
		}
	}
	void ElementSet::SetSendTo(std::set<Storage::integer> & procs, char dir) const
	{
		SendTo(procs);
		if( (dir & 1) && HaveChild() )
		{
			for(ElementSet it = GetChild(); it != InvalidElementSet(); it = it.GetSibling() )
				it->SetSendTo(procs, dir & 1); // don't let children to go upwards
		}
		if( (dir & 2) && HaveParent() )
		{
			GetParent()->SetSendTo(procs,dir & 2); // don't let parent to go downwards
		}
	}
	void ElementSet::SetSendTo(std::vector<Storage::integer> & procs, char dir) const
	{
		SendTo(procs);
		if( (dir & 1) && HaveChild() )
		{
			for(ElementSet it = GetChild(); it != InvalidElementSet(); it = it.GetSibling() )
				it->SetSendTo(procs, dir & 1); // don't let children to go upwards
		}
		if( (dir & 2) && HaveParent() )
		{
			GetParent()->SetSendTo(procs,dir & 2); // don't let parent to go downwards
		}
	}
	void ElementSet::SetSendTo(Storage::integer_array procs, char dir) const
	{
		SendTo(procs);
		if( (dir & 1) && HaveChild() )
		{
			for(ElementSet it = GetChild(); it != InvalidElementSet(); it = it.GetSibling() )
				it->SetSendTo(procs, dir & 1); // don't let children to go upwards
		}
		if( (dir & 2) && HaveParent() )
		{
			GetParent()->SetSendTo(procs,dir & 2); // don't let parent to go downwards
		}
	}
	void ElementSet::CollectProcessors(std::set<Storage::integer> & procs, char dir) const
	{
		{
			Storage::integer_array set_procs = IntegerArray(GetMeshLink()->ProcessorsTag());
			procs.insert(set_procs.begin(),set_procs.end());
		}
		if( (dir & 1) && HaveChild() )
		{
			for(ElementSet it = GetChild(); it != InvalidElementSet(); it = it.GetSibling() )
				it.CollectProcessors(procs,dir & 1); //don't let children to go upwards
		}
		if( (dir & 2) && HaveParent() )
		{
			GetParent().CollectProcessors(procs,dir & 2); //don't let parent to go downwards
		}
	}
	void ElementSet::SynchronizeSetChildren()
	{
		if( GetMeshLink()->GetMeshState() == Mesh::Serial ) return;
		
		if( GetStatus() != Element::Owned )
		{
			std::set<Storage::integer> send_set;
			CollectProcessors(send_set,1); //collect procs from children
			if( !send_set.empty() ) SetSendTo(send_set,1);
		}
	}

	void ElementSet::SynchronizeSetParents()
	{
		if( GetMeshLink()->GetMeshState() == Mesh::Serial ) return;
		if( GetStatus() != Element::Owned )
		{
			std::set<Storage::integer> send_set;
			CollectProcessors(send_set,0); //don't collect procs from parents
			//std::cout << GetMeshLink()->GetProcessorRank() << " SynchronizeSetParents " << GetName() << " procs ";
			//for(std::set<Storage::integer>::iterator it = send_set.begin(); it != send_set.end(); ++it) std::cout << *it << " ";
			if( !send_set.empty() ) SetSendTo(send_set,2);
			//std::cout << std::endl;
		}
	}
	
}

#endif
