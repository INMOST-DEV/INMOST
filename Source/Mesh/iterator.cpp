#include "inmost.h"
#include <iostream>
#if defined(USE_MESH)





namespace INMOST
{

	ElementType GetNextType(ElementType current, ElementType types)
	{
		ElementType ret = MESH << 1;
		for(ElementType i = current << 1; i <= MESH; i = i << 1)
			if( types & i )
			{
				ret = i;
				break;
			}
		if( ret > MESH ) ret = NONE;
		return ret;
	}

	ElementType GetPrevType(ElementType current, ElementType types)
	{
		ElementType ret = MESH << 1;
		for(ElementType i = current >> 1; i > NONE; i = i >> 1)
			if( types & i )
			{
				ret = i;
				break;
			}
		if( ret > MESH ) ret = NONE;
		return ret;
	}
	
	HandleType  Mesh::NextHandle(HandleType h) const 
	{
		integer num = GetHandleElementNum(h), id = GetHandleID(h);
		++id;
		while( num < 5 )
		{
			while(id < static_cast<integer>(links[num].size()) && (links[num][id] == -1)) ++id;
			if( id == static_cast<integer>(links[num].size()) ) 
			{
				id = 0;
				++num;
			}
			else break;
		}
		if( num == 5 && id > 0 ) id = 1;
		return ComposeHandleNum(num,id);
	}

	HandleType  Mesh::PrevHandle(HandleType h) const 
	{
		integer num = GetHandleElementNum(h), id = GetHandleID(h);
		--id;
		if( num == 5 ) 
		{
			if( id < 0 )
				num = 4;
			else return ComposeHandleNum(ElementNum(MESH),0);
		}
		while( num >= 0 )
		{
			while(id >= 0 && (links[num][id] == -1)) --id;
			if( id == -1 ) 
			{
				--num;
				if( num > 0 ) id = static_cast<integer>(links[num].size())-1;
			}
			else break;
		}
		if( num < 0 ) return InvalidHandle();
		return ComposeHandleNum(num,id);
	}

	HandleType  Mesh::NextHandle(HandleType h, ElementType etype) const 
	{
		integer num = GetHandleElementNum(h), id = GetHandleID(h);
		++id;
		while( num < 5 )
		{
			while(id < static_cast<integer>(links[num].size()) && (links[num][id] == -1)) ++id;
			if( id == static_cast<integer>(links[num].size()) ) 
			{
				bool stop = true;
				for(integer q = num+1; q < 5; q++)
				{
					if( etype & (1 << q) )
					{
						id = 0;
						num = q;
						stop = false;
						break;
					}
				}
				if( stop ) break;
			}
			else break;
		}
		if( num == 5 && id > 0 ) id = 1;
		return ComposeHandleNum(num,id);
	}

	HandleType  Mesh::PrevHandle(HandleType h, ElementType etype) const 
	{
		integer num = GetHandleElementNum(h), id = GetHandleID(h);
		--id;
		if( num == 5 ) 
		{
			if( id < 0 )
			{
				bool stop = true;
				for(integer q = 4; q >= 0; q--)
				{
					if( etype & (1 << q) )
					{
						
						num = q;
						id = static_cast<integer>(links[num].size())-1;
						stop = false;
						break;
					}
				}
				if( stop ) return InvalidHandle();
			}
			else return ComposeHandleNum(ElementNum(MESH),0);
		}
		while( num >= 0 )
		{
			while(id >= 0 && (links[num][id] == -1)) --id;
			if( id == -1 ) 
			{
				bool stop = true;
				for(integer q = num-1; q >= 0; q--)
				{
					if( etype & (1 << q) )
					{
						num = q;
						id = static_cast<integer>(links[num].size())-1;
						stop = false;
						break;
					}
				}
				if( stop ) return InvalidHandle();
			}
			else break;
		}
		if( num < 0 ) return InvalidHandle();
		return ComposeHandleNum(num,id);
	}

	Storage::integer  Mesh::FirstLocalID(ElementType etype) const 
	{
		assert(OneType(etype)); 
		integer ret = 0, n = ElementNum(etype); 
		if(n == 5) return 0; 
		while(ret < static_cast<integer>(links[n].size()) && (links[n][ret] == -1)) ++ret; 
		return ret;
	}

	Storage::integer Mesh::NumberOf(ElementType t) const
	{
		integer ret = 0;
		for(int m = 0; m < 5; m++) if( (1 << m) & t )
			ret += static_cast<integer>(links[m].size() - empty_links[m].size()) - hidden_count[m];
		if( t & MESH ) ret++;
		return ret;
	}
	
	template<typename EType>
	Mesh::base_iterator<EType> & Mesh::base_iterator<EType>::operator ++()
	{
		while( etype != NONE ) 
		{
			lid = m->NextLocalIDIter(etype,lid);
			if( lid == m->LastLocalID(etype) )
			{
				etype = GetNextType(etype,types);
				lid = -1;
				if( !etype ) break;
			}
			else break;
		}
		return *this;
	}

	template<typename EType>
	Mesh::base_iterator<EType> & Mesh::base_iterator<EType>::operator --()
	{
		while( etype != NONE ) 
		{
			lid = m->PrevLocalIDIter(etype,lid);
			if( lid == -1 )
			{
				etype = GetPrevType(etype,types);
				if( etype ) 
					lid = m->LastLocalIDIter(etype);
				else
				{
					lid = -1;
					break;
				}
			}
			else break;
		}
		return *this;
	}
	
	template<typename EType>
	Mesh::base_iterator<EType>::base_iterator(ElementType T, Mesh * _m, bool last)
	{
		m = _m;
		types = T;
		etype = NONE;
		lid = -1;
		if( last )
		{
			for(ElementType i = MESH; i >= NODE; i = i >> 1) if( types & i ) {etype = i; break;}
			lid = m->LastLocalID(etype);
			operator--();
		}
		else
		{
			for(ElementType i = NODE; i <= MESH; i = i << 1) if( types & i ) {etype = i; break;}
			lid = -1;
			operator++();
		}
	}
	template<typename EType>
	void Mesh::base_iterator<EType>::Print()
	{
		printf("Number: %10d CurrentType %x types %x\n",lid,etype,types);
	}

	Storage::integer Mesh::NextLocalIDIter(ElementType etype, integer lid) const 
	{
		integer q = ElementNum(etype); 
		++lid; 
		while(lid < static_cast<integer>(links[q].size()) && (links[q][lid] == -1|| Hidden(ComposeHandleNum(q,lid)))) ++lid; 
		return lid;
	}

	Storage::integer Mesh::PrevLocalIDIter(ElementType etype, integer lid) const 
	{
		integer q = ElementNum(etype); 
		--lid; 
		while(lid > 0 && (links[q][lid] == -1 || Hidden(ComposeHandleNum(q,lid)))) --lid; 
		return lid;
	}
	Storage::integer Mesh::FirstLocalIDIter(ElementType etype) const
	{
		assert(OneType(etype)); 
		integer ret = 0, n = ElementNum(etype); 
		if(n == 5) return 0; 
		while(ret < static_cast<integer>(links[n].size()) && (links[n][ret] == -1|| Hidden(ComposeHandleNum(n,ret)))) ++ret; 
		return ret;
	}
	Storage::integer Mesh::LastLocalIDIter(ElementType etype) const 
	{
		assert(OneType(etype)); 
		int ret = LastLocalIDNum(ElementNum(etype));
		return PrevLocalIDIter(etype,ret);
	}

	Mesh::iteratorStorage Mesh::Begin(ElementType Types)        {return base_iterator<Storage>(Types,this,false);}
	Mesh::iteratorStorage Mesh::End()                           {return base_iterator<Storage>(this);}
	Mesh::iteratorElement Mesh::BeginElement(ElementType Types) {return base_iterator<Element>(Types & (NODE | EDGE | FACE | CELL | ESET),this,false);}
	Mesh::iteratorElement Mesh::EndElement()                    {return base_iterator<Element>(this);}
	Mesh::iteratorSet     Mesh::BeginSet()                      {return base_iterator<ElementSet>(ESET,this,false);}
	Mesh::iteratorSet     Mesh::EndSet()                        {return base_iterator<ElementSet>(this);}
	Mesh::iteratorCell    Mesh::BeginCell()                     {return base_iterator<Cell>(CELL,this,false);}
	Mesh::iteratorCell    Mesh::EndCell()                       {return base_iterator<Cell>(this);}
	Mesh::iteratorFace    Mesh::BeginFace()                     {return base_iterator<Face>(FACE,this,false);}
	Mesh::iteratorFace    Mesh::EndFace()                       {return base_iterator<Face>(this);}
	Mesh::iteratorEdge    Mesh::BeginEdge()                     {return base_iterator<Edge>(EDGE,this,false);}
	Mesh::iteratorEdge    Mesh::EndEdge()                       {return base_iterator<Edge>(this);}
	Mesh::iteratorNode    Mesh::BeginNode()                     {return base_iterator<Node>(NODE,this,false);}
	Mesh::iteratorNode    Mesh::EndNode()                       {return base_iterator<Node>(this);}

	//all possible templates
	template class Mesh::base_iterator<Element>;
	template class Mesh::base_iterator<ElementSet>;
	template class Mesh::base_iterator<Node>;
	template class Mesh::base_iterator<Edge>;
	template class Mesh::base_iterator<Face>;
	template class Mesh::base_iterator<Cell>;
	template class Mesh::base_iterator<Storage>;
}
#endif
