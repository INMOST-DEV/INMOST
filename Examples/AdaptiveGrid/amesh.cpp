#include "agrid.h"
#include <math.h>


Storage::integer min_nonzero(Storage::integer a, Storage::integer b)
{
	Storage::integer get_max = std::max(a,b);
	Storage::integer get_min = std::min(a,b);
	if( get_min == 0 ) return get_max; else return get_min;
}

void AdaptiveGrid

void AdaptiveGrid::Refine(const TagInteger & indicator)
{
	int schedule_counter = 1;
	bool scheduled = true;
	//0. Extend indicator for edges and faces
	indicator = CreateTag(indicator.GetTagName(),DATA_INTEGER,FACE|EDGE,NONE,1);
	while(scheduled)
	{
		//1.Communicate indicator - it may be not synced
		ExchangeTag(indicator,0,CELL);
		//2.Propogate indicator down to the faces,edges
		//  select latest schedule for them
		//
		//possible problem - a cell needs all it's elements to be splitted in the same schedule
		//so the cell must mark all it's elements without hanging nodes into the same schedule,
		//unless there is an even earlier cell
		//
		//should we select earliest possible for elements without hanging nodes and
		//latest possible for elements with hanging nodes?
		//
		//with loop over cells we would mark elements without hanging nodes of the cells
		//in current schedule
		for(ElementType etype = FACE; etype >= EDGE; etype = PrevElementType(etype))
		{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
			for(Storage::integer it = 0; it < LastLocalID(etype); ++it) if( isValidElement(etype,it) )
			{
				Element e = ElementByLocalID(etype,it);
				ElementArray<Element> adj = f.getAdjElements(NextElementType(etype));
				//here latest schedule is selected
				if( e->nbAdjElements(NODE,hanging_node,0) > 0 ) //latest possible schedule
				{
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
						if( indicator[adj[kt]] )
							indicator[f] = indicator[f] ? std::min(indicator[f],indicator[adj[kt]]) : indicator[adj[kt]];
				}
				else //earliest possible schedule
				{
					for(ElementArray<Element>::size_type kt = 0; kt < adj.size(); ++kt)
						indicator[f] = std::max(indicator[f],indicator[adj[kt]]);
				}
				
			}
		}
		//3.Communicate indicator on faces and edges
		ExchangeTag(indicator,0,FACE|EDGE);
		//4.Check for each cell without indicator if there is
		//  any hanging node with adjacent in a need to refine,
		//  schedule for refinement earlier.
		scheduled = false;
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(Storage::integer it = 0; it < CellLastLocalID(); ++it) if( isValidCell(it) )
		{
			Cell c = CellByLocalID(it);
			if( indicator[c] == 0 )
			{
				bool scheduled_c = false;
				//retrive only hanging nodes, this may be empty
				ElementArray<Node> hanging = c->getNodes(hanging_node,0);
				for(ElementArray<Node>::size_type kt = 0; kt < hanging.size() && !scheduled_c; ++kt)
				{
					//adjacent edges may be scheduled for refinement
					ElementArray<Edge> adj = hanging[kt].getEdges();
					for(ElementArray<Edge>::size_type lt = 0; lt < adj.size() && !scheduled_c; ++lt)
						if( indicator[adj[lt]] != 0 )
						{
							indicator[c] = schedulde_counter+1;
							scheduled = scheduled_c = true;
						}
				}
			}
		}
		//5.Go back to 1 until no new elements scheduled
		if( scheduled ) schedulde_counter++;
	}
	//5.Refine
	BeginModification();
	while(scheduled_counter)
	{
		scheduled_counter--;
	}
	ResolveModification();
	//Let the user update their data
	ApplyModification();
	EndModification();
}

void AdaptiveGrid::Coarse(const TagInteger & indicator)
{
}
