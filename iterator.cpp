#include "inmost.h"
#include <iostream>
#if defined(USE_MESH)


INMOST::ElementType GetNextType(INMOST::ElementType current, INMOST::ElementType types)
{
	INMOST::ElementType ret = INMOST::MESH << 1;
	for(INMOST::ElementType i = current << 1; i <= INMOST::MESH; i = i << 1)
		if( types & i )
		{
			ret = i;
			break;
		}
	if( ret > INMOST::MESH ) ret = INMOST::NONE;
	return ret;
}

INMOST::ElementType GetPrevType(INMOST::ElementType current, INMOST::ElementType types)
{
	INMOST::ElementType ret = INMOST::MESH << 1;
	for(INMOST::ElementType i = current >> 1; i > INMOST::NONE; i = i >> 1)
		if( types & i )
		{
			ret = i;
			break;
		}
	if( ret > INMOST::MESH ) ret = INMOST::NONE;
	return ret;
}


namespace INMOST
{
	int Mesh::GetTypeEnd(ElementType t)
	{
		switch(t)
		{
			case NODE: return nodes.size()-1;
			case EDGE: return edges.size()-1;
			case FACE: return faces.size()-1;
			case CELL: return cells.size()-1;
			case ESET: return sets.size()-1;
			case MESH: return 0;
		}
		return -1;
	}	
	INMOST_DATA_ENUM_TYPE Mesh::NumberOf(ElementType t)
	{
		INMOST_DATA_ENUM_TYPE ret = 0;
		for(ElementType m = NODE; m <= MESH; m = m << 1)
			if( m & t )
			{
				switch(m)
				{
					case NODE: ret += NumberOfNodes(); break;
					case EDGE: ret += NumberOfEdges(); break;
					case FACE: ret += NumberOfFaces(); break;
					case CELL: ret += NumberOfCells(); break;
					case ESET: ret += NumberOfSets(); break;
					case MESH: ret += 1; break;
				}
			}
		return ret;
	}
	
	Mesh::base_iterator & Mesh::base_iterator::operator ++()
	{
		if( CurrentType != NONE ) Number++;
		switch(CurrentType)
		{
			case NODE: 
			{
				while(Number < static_cast<itenum>(m->nodes.size()) && (m->nodes[Number] == NULL || m->nodes[Number]->GetMarker(m->HideMarker()))) Number++;
				if( Number >= static_cast<itenum>(m->nodes.size()) )
				{
					Number = 0;
					CurrentType = GetNextType(CurrentType,types);
					if( CurrentType == NONE )
					{
						Number = -1;
						break;
					}
				}
				else break;
			}
			case EDGE: 
			{
				if( CurrentType == EDGE )
				{
					while(Number < static_cast<itenum>(m->edges.size()) && (m->edges[Number] == NULL || m->edges[Number]->GetMarker(m->HideMarker()))) Number++;
					if( Number >= static_cast<itenum>(m->edges.size()) )
					{
						Number = 0;
						CurrentType = GetNextType(CurrentType,types);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
			case FACE: 
			{
				if( CurrentType == FACE )
				{
					while(Number < static_cast<itenum>(m->faces.size()) && (m->faces[Number] == NULL || m->faces[Number]->GetMarker(m->HideMarker()))) Number++;
					if( Number >= static_cast<itenum>(m->faces.size()) )
					{
						Number = 0;
						CurrentType = GetNextType(CurrentType,types);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
			case CELL:
			{
				if( CurrentType == CELL )
				{
					while(Number < static_cast<itenum>(m->cells.size()) && (m->cells[Number] == NULL || m->cells[Number]->GetMarker(m->HideMarker()))) Number++;
					if( Number == static_cast<itenum>(m->cells.size()) )
					{
						Number = 0;
						CurrentType = GetNextType(CurrentType,types);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
			case ESET:
			{
				if( CurrentType == ESET )
				{
					while(Number < static_cast<itenum>(m->sets.size()) && (m->sets[Number] == NULL || m->sets[Number]->GetMarker(m->HideMarker()))) Number++;
					if( Number == static_cast<itenum>(m->sets.size()) )
					{
						Number = 0;
						CurrentType = GetNextType(CurrentType,types);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
			case MESH:
			{
				if( CurrentType == MESH )
				{
					if( Number == 1 )
					{
						Number = 0;
						CurrentType = GetNextType(CurrentType,types);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
		}
		return *this;
	}
	Mesh::base_iterator & Mesh::base_iterator::operator --()
	{
		if( CurrentType != NONE ) Number--;
		switch(CurrentType)
		{
			case MESH:
			{
				if( Number < 0 )
				{
					CurrentType = GetPrevType(CurrentType,types);
					Number = m->GetTypeEnd(CurrentType);
					if( CurrentType == NONE )
					{
						Number = -1;
						break;
					}
				}
				else break;
			}
			case ESET:
			{
				if( CurrentType == ESET )
				{
					while(Number >= 0 && (m->sets[Number] == NULL || m->sets[Number]->GetMarker(m->HideMarker())))	Number--;
					if( Number < 0 )
					{
						CurrentType = GetPrevType(CurrentType,types);
						Number = m->GetTypeEnd(CurrentType);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;	
				}
			}
			case CELL:
			{
				if( CurrentType == CELL )
				{
					while(Number >= 0 && (m->cells[Number] == NULL || m->cells[Number]->GetMarker(m->HideMarker()))) Number--;
					if( Number < 0 )
					{
						CurrentType = GetPrevType(CurrentType,types);
						Number = m->GetTypeEnd(CurrentType);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
			case FACE: 
			{
				if( CurrentType == FACE )
				{
					while(Number >= 0 && (m->faces[Number] == NULL || m->faces[Number]->GetMarker(m->HideMarker()))) Number--;
					if( Number < 0 )
					{
						CurrentType = GetPrevType(CurrentType,types);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
			case EDGE: 
			{
				if( CurrentType == EDGE )
				{
					while(Number >= 0 && (m->edges[Number] == NULL || m->edges[Number]->GetMarker(m->HideMarker()))) Number--;
					if( Number < 0 )
					{
						CurrentType = GetPrevType(CurrentType,types);
						Number = m->GetTypeEnd(CurrentType);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
			case NODE: 
			{
				if( CurrentType == NODE )
				{
					while(Number >= 0 && (m->nodes[Number] == NULL || m->nodes[Number]->GetMarker(m->HideMarker()))) Number--;
					if( Number < 0 )
					{
						CurrentType = GetPrevType(CurrentType,types);
						Number = m->GetTypeEnd(CurrentType);
						if( CurrentType == NONE )
						{
							Number = -1;
							break;
						}
					}
					else break;
				}
			}
		}
		return *this;
	}
	
	
	Mesh::base_iterator::base_iterator(ElementType T, Mesh * _m, bool last)
	{
		m = _m;
		types = T;
		CurrentType = NONE;
		if( last )
		{
			for(ElementType i = MESH; i >= NODE; i = i >> 1) if( types & i ) {CurrentType = i; break;}
			Number = m->GetTypeEnd(CurrentType)+1;
			operator--();
		}
		else
		{
			for(ElementType i = NODE; i <= MESH; i = i << 1) if( types & i ) {CurrentType = i; break;}
			Number = -1;
			operator++();
		}
	}
	void Mesh::base_iterator::Print()
	{
		printf("Number: %10d CurrentType %x types %x\n",Number,CurrentType,types);
	}

	Mesh::base_iterator   Mesh::Begin(ElementType Types)        {return base_iterator(Types,this,false);}
	Mesh::base_iterator   Mesh::End()                           {return base_iterator(this);}
	Mesh::iteratorElement Mesh::BeginElement(ElementType Types) {return static_cast<iteratorElement>(Begin(Types & (NODE | EDGE | FACE | CELL)));}
	Mesh::iteratorElement Mesh::EndElement()                    {return static_cast<iteratorElement>(End());}
	Mesh::iteratorSet     Mesh::BeginSet()                      {return static_cast<iteratorSet>(Begin(ESET));}
	Mesh::iteratorSet     Mesh::EndSet()                        {return static_cast<iteratorSet>(End());}
	Mesh::iteratorCell    Mesh::BeginCell()                     {return static_cast<iteratorCell>(Begin(CELL));}
	Mesh::iteratorCell    Mesh::EndCell()                       {return static_cast<iteratorCell>(End());}
	Mesh::iteratorFace    Mesh::BeginFace()                     {return static_cast<iteratorFace>(Begin(FACE));}
	Mesh::iteratorFace    Mesh::EndFace()                       {return static_cast<iteratorFace>(End());}
	Mesh::iteratorEdge    Mesh::BeginEdge()                     {return static_cast<iteratorEdge>(Begin(EDGE));}
	Mesh::iteratorEdge    Mesh::EndEdge()                       {return static_cast<iteratorEdge>(End());}
	Mesh::iteratorNode    Mesh::BeginNode()                     {return static_cast<iteratorNode>(Begin(NODE));}
	Mesh::iteratorNode    Mesh::EndNode()                       {return static_cast<iteratorNode>(End());}	
}
#endif
