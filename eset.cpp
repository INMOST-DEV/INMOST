#include "inmost.h"
#if defined(USE_MESH)
namespace INMOST
{
	
	ElementSet::ElementSet(bool ordered)
	:Storage(NULL,ESET),ordered(ordered),eset(ordered?LessElementsUnique:LessElementsPointer)
	{
	}
	
	ElementSet::ElementSet(Mesh * m, bool ordered)
	:Storage(m,ESET),ordered(ordered),eset(ordered?LessElementsUnique:LessElementsPointer)
	{
	}
	
	
	ElementSet::ElementSet(const ElementSet & other)
		:Storage(other.m_link,other.LocalID(),other),ordered(other.ordered),eset(other.eset)
	{
	}

	ElementSet::ElementSet(Mesh * m, INMOST_DATA_ENUM_TYPE lid,const ElementSet & other)
	:Storage(m,lid,other),ordered(other.ordered)
	{
		//in this case mesh fills eset array
	}
	
	ElementSet & ElementSet::operator =(ElementSet const & other)
	{
		Storage::operator =(other);
		eset = other.eset;
		return *this;
	}


	ElementSet::~ElementSet()
	{
		
#if defined(USE_SET_DELETE_TEMPORARY)
		for(element_set_type::iterator it = eset.begin(); it != eset.end(); it++) 
			if( isTemporalElement((*it)) ) delete (*it);
#endif
		//Mesh::FindMesh(GetMeshID())->UntieElement(this);
		eset.clear();
	}
	
	std::pair< ElementSet::iterator, bool > ElementSet::Insert(const Element * e)
	{
		if( e == NULL ) throw NullInElementSet;
		std::pair<element_set_type::iterator,bool> ret = eset.insert(const_cast<Element *>(e));
		return std::pair<ElementSet::iterator, bool>(ElementSet::iterator(ret.first),ret.second);
	}
	
	void ElementSet::Insert(ElementSet e)
	{
		eset.insert(e.eset.begin(),e.eset.end());
	}
	
	void ElementSet::Insert(std::vector<Element *> ve)
	{
		std::sort(ve.begin(),ve.end(),eset.key_comp());
		eset.insert(ve.begin(),ve.end());
	}
	
	
	
	void ElementSet::Insert(array<Element *> ve)
	{
		std::sort(ve.begin(),ve.end(),eset.key_comp());
		eset.insert(ve.begin(),ve.end());
	}
	
	bool ElementSet::Erase(Element * e)
	{
		element_set_type::iterator it = eset.find(e);
		if( it != eset.end() )
		{
#if defined(USE_SET_DELETE_TEMPORARY)
			if( isTmporalElement(*it) ) delete (*it);
#endif
			eset.erase(it);
			return true;
		}
		return false;
	}
	
	void ElementSet::Erase(iterator e)
	{
		eset.erase(e);
	}
	void ElementSet::Intersection(ElementSet other)
	{
		element_set_type::value_compare comparator = eset.value_comp();
		element_set_type::iterator it = eset.begin(), jt = other.eset.begin();
		while( (it != eset.end()) && (jt != other.eset.end()) )
		{
			if( comparator((*jt),(*it)) ) jt++;
			else if( comparator((*it),(*jt)) ) 
			{
#if defined(USE_SET_DELETE_TEMPORARY)
				if( isTmporalElement(*it) ) delete (*it);
#endif
				eset.erase(it++);
			}
			else{it++; jt++;}
		}
		eset.erase(it,eset.end());
	}
	
	
	void ElementSet::Union(ElementSet other)
	{
		eset.insert(other.eset.begin(),other.eset.end());
		//iterator it = eset.begin();
		//for(iterator jt = other.eset.begin(); jt != other.eset.end(); jt++)
		//	it = eset.insert(it,(*jt));
	}
	
	void ElementSet::Difference(ElementSet other)
	{
		element_set_type::value_compare comparator = eset.value_comp();
		element_set_type::iterator it = eset.begin(), jt = other.eset.begin();
		while( (it != eset.end()) && (jt != other.eset.end()) )
		{
			if( comparator((*jt),(*it)) ) jt++;
			else if( comparator((*it),(*jt)) ) it++;
			else 
			{
#if defined(USE_SET_DELETE_TEMPORARY)
				if( isTmporalElement(*it) ) delete (*it);
#endif
				eset.erase(it++);
			}
		}
	}
	
	ElementSet::iterator ElementSet::begin()
	{
		return ElementSet::iterator(eset.begin());
	}
	
	ElementSet::iterator ElementSet::end()
	{
		return ElementSet::iterator(eset.end());
	}
	
	ElementSet::reverse_iterator ElementSet::rbegin()
	{
		return ElementSet::reverse_iterator(eset.rbegin());
	}
	
	ElementSet::reverse_iterator ElementSet::rend()
	{
		return ElementSet::reverse_iterator(eset.rend());
	}
	
	
	size_t ElementSet::size() const
	{
		return eset.size();
	}
	bool ElementSet::empty() const
	{
		return eset.empty();
	}
	
	void ElementSet::clear()
	{
#if defined(USE_SET_DELETE_TEMPORARY)
		for(iteartor it = eset.begin(); it != eset.end(); it++)
			if( isTmporalElement(*it) ) delete (*it);
#endif
		eset.clear();
	}
	
	ElementSet::iterator ElementSet::find(Element * e)
	{
		return ElementSet::iterator(eset.find(e));
	}
	
	
	void ElementSet::SetElementsMarker(MIDType marker)
	{
		for(ElementSet::iterator it = begin(); it != end(); it++)
			it->SetMarker(marker);
	}
	void ElementSet::RemElementsMarker(MIDType marker)
	{
		for(ElementSet::iterator it = begin(); it != end(); it++)
			it->RemMarker(marker);
	}
}
#endif
