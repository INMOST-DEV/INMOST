#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"
#include "../Mesh/incident_matrix.hpp"
#if defined(USE_MESH)

// coords/zcorn algorithm
// 1) put all block nodes into pillars, sort each pillar nodes by actual depth
// 2) create edges down along pillars and mark them according to blocks
// 3) consider all pairs of pillars in Oxz and Oyz planes and create uncut block edges, mark them in block numbers
// 4) use line-sweep algorithm to intersect and cut all the edges between each pair of pillars, mark them in block numbers
// 5) mark all the nodes with block numbers that they belong by considering union of block numbers on adjacent edges
// 5) use incident_matrix algorithm to create all the faces between pair of pillars from pillar edges and inter-pillar edges
// 6) from intersection of block numbers on nodes figure out which blocks the face belongs
// 7) add top and bottom interface

//todo
// 1. (ok) remove edges on pillars that do not get any associated block
// 2. (ok) do not create cells with actnum = 0
// 2. fix fast intersect() algorithm, check against intersect_naive()
// 3. CheckEdgeOrder reports bad order of edges for certain interfaces - investigate
// 4. when encounter cells with one face with actnum = 1 should disconnect adjacent cells
// 5. (ok) populate keywords data, poro, perm, actnum, satnum, etc...
// 5.1 add saturation, pressure
// 6. local grid refinement
// 6.1 radial local grid refinement
// 6.2 nested local grid refinement
// 7. (ok) omp-parallel loading
// 8. mpi-parallel loading
// 9. (ok) do not convert arrays xyz and zcorn
//10. read wells into sets
//11. optimize:
//11.1 detect that pillars have no faults, skip intersection in this case
//11.2 skip incident_matrix algorithm when no inner edges found

//eclipse states
#define ECLSTRCMP(x,y) strncmp(x,y,8)

#define ECL_NEVER -1
#define ECL_NONE 0
#define ECL_SKIP_SECTION 1
#define ECL_INCLUDE 2
#define ECL_DIMENS 3
#define ECL_DX 4
#define ECL_DY 5
#define ECL_DZ 6
#define ECL_TOPS 7
#define ECL_PERMX 8
#define ECL_PERMY 9
#define ECL_PERMZ 10
#define ECL_PORO 11
#define ECL_MAPAXIS 12
#define ECL_INRAD 13
#define ECL_COORDS 14
#define ECL_ZCORN 15
#define ECL_ACTNUM 16
#define ECL_SATNUM 17

#define ECL_GTYPE_NONE 0
#define ECL_GTYPE_TOPS 1
#define ECL_GTYPE_ZCORN 2
#define ECL_GTYPE_RADIAL 3
#define ECL_GTYPE_CARTESIAN 4


#define ECL_VAR_NONE 0
#define ECL_VAR_REAL 1
#define ECL_VAR_INT 2

#define ECL_IJK_DATA(i,j,k) (i + ((j)+(k)*dims[1])*dims[0])
#define ECL_IJK_COORDS(i,j,l,coord) ((((i) + (j)*(dims[0]+1))*2+(l))*3+(coord))
#define ECL_IJK_ZCORN(i,j,k,l) ((k)*dims[1]*dims[0]*8 + (1-(l)/4)*dims[1]*dims[0]*4 + (j)*dims[0]*4 + ((l)/2%2)*dims[0]*2 + (i)*2 + (l)%2)

//line sweep events
#define SEG_START 0
#define SEG_END 1



namespace INMOST
{
	static std::string GetFolder(std::string file)
	{
		size_t found = file.find_last_of("/\\");
		if( found == std::string::npos )
			return "";
		else return file.substr(0,found);
	}
	//2d point with comparison operator for line sweep algorithm
	struct Point
	{
		const Storage::real eps;
		Storage::real x, y;
		Point & operator = (Point const & b) {  x = b.x; y = b.y; return *this; }
		Point(const Point & b) : eps(1.0e-7), x(b.x), y(b.y) {}
		Point(Storage::real _x, Storage::real _y) : eps(1.0e-7), x(_x), y(_y) {}
		Point(Storage::real v[2]) : eps(1.0e-7), x(v[0]), y(v[1]) {}
		bool operator <(const Point & b) const
		{
			if (y < b.y - eps) return true;
			else if (y > b.y + eps) return false;
			else if (x < b.x - eps) return true;
			else return false;
		}
		bool operator ==(const Point & b) const
		{
			return fabs(y - b.y) < eps && fabs(x - b.x) < eps;
		}
		bool operator !=(const Point & b) const
		{
			return fabs(y - b.y) > eps || fabs(x - b.x) > eps;
		}
	};
	//class that returns a Point into 3d space based on coords of two pillars
	class Unproject
	{
		Storage::real a0[3], a1[3], b0[3], b1[3];
	public:
		Unproject(Storage::real _a0[3], Storage::real _a1[3], Storage::real _b0[3], Storage::real _b1[3])
		{
			memcpy(a0,_a0,sizeof(Storage::real)*3);
			memcpy(a1,_a1,sizeof(Storage::real)*3);
			memcpy(b0,_b0,sizeof(Storage::real)*3);
			memcpy(b1,_b1,sizeof(Storage::real)*3);
		}
		Unproject(const Unproject & b)
		{
			memcpy(a0,b.a0,sizeof(Storage::real)*3);
			memcpy(a1,b.a1,sizeof(Storage::real)*3);
			memcpy(b0,b.b0,sizeof(Storage::real)*3);
			memcpy(b1,b.b1,sizeof(Storage::real)*3);
		}
		Unproject & operator =(Unproject const & b)
		{
			memmove(a0,b.a0,sizeof(Storage::real)*3);
			memmove(a1,b.a1,sizeof(Storage::real)*3);
			memmove(b0,b.b0,sizeof(Storage::real)*3);
			memmove(b1,b.b1,sizeof(Storage::real)*3);
			return *this;
		}
		void Act(const Point & p, Storage::real v[3]) const
		{
			Storage::real a, b;
			for(int k = 0; k < 3; ++k)
			{
				a = (p.x)*a0[k] + (1-p.x)*a1[k];
				b = (p.x)*b0[k] + (1-p.x)*b1[k];
				v[k] = (1-p.y)*a+p.y*b;
			}
		}
	};
	//Comparator for events in line sweep algorithm
	class event_less
	{
	public:
		bool operator()(const std::pair<Storage::real, int> & a, const std::pair<Storage::real, int> & b) const
		{
			const Storage::real eps = 1.0e-7;
			if (a.first < b.first - eps)
				return true;
			else if (a.first > b.first + eps)
				return false;
			else if (a.second < b.second)
				return true;
			return false;
		}
	};
	//Comparator for depth of nodes in pillar
	class pillar_less
	{
	public:
		bool operator()(Storage::real a, Storage::real b) const
		{
			const Storage::real eps = 1.0e-7;
			if( a < b - eps)
				return true;
			return false;
		}
	};
	class index_comparator
	{
		Storage::integer_array & data;
	public:
		index_comparator(Storage::integer_array & data) : data(data) {}
		index_comparator(const index_comparator & b) : data(b.data) {}
		index_comparator & operator =(index_comparator const & b) {data = b.data; return *this;}
		bool operator ()(int a, int b) const {return data[a] < data[b];}
	};
	template<typename T>
	int count_duplicates(ElementArray<T> & array)
	{
		Mesh * m = array.GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		int dups = 0;
		for(int k = 0; k < array.size(); ++k)
		{
			if( array[k].GetPrivateMarker(mrk) )
				dups++;
			array[k].SetPrivateMarker(mrk);
		}
		array.RemPrivateMarker(mrk);
		m->ReleasePrivateMarker(mrk);
		return dups;
	}
	
	template<typename T>
	void make_unique(ElementArray<T> & array)
	{
		Mesh * m = array.GetMeshLink();
		MarkerType mrk = m->CreatePrivateMarker();
		int i = 0, j = 0;
		while(j < array.size())
		{
			if( array[j].GetPrivateMarker(mrk) )
				++j;
			else
			{
				array[j].SetPrivateMarker(mrk);
				array[i++] = array[j++];
			}
		}
		array.RemPrivateMarker(mrk);
		array.resize(i);
		m->ReleasePrivateMarker(mrk);
	}
	
	
	
	std::pair<bool,Node> intersect_segments(Mesh * m, const Edge & a, const Edge & b, std::map<Point,Node> & intersections, Tag pnt, const Unproject & unp, bool print)
	{
		const Storage::real eps = 1.0e-9;
		if( a->getBeg() == b->getBeg() || a->getBeg() == b->getEnd() || a->getEnd() == b->getEnd() || a->getEnd() == b->getBeg() )
			return std::make_pair(false,InvalidNode());
		Storage::real_array abeg = a->getBeg()->Coords();
		Storage::real_array aend = a->getEnd()->Coords();
		Storage::real_array bbeg = b->getBeg()->Coords();
		Storage::real_array bend = b->getEnd()->Coords();
		Storage::real_array _pabeg = a->getBeg()->RealArray(pnt);
		Storage::real_array _paend = a->getEnd()->RealArray(pnt);
		Storage::real_array _pbbeg = b->getBeg()->RealArray(pnt);
		Storage::real_array _pbend = b->getEnd()->RealArray(pnt);
		Point pabeg(_pabeg[0],_pabeg[1]);
		Point paend(_paend[0],_paend[1]);
		Point pbbeg(_pbbeg[0],_pbbeg[1]);
		Point pbend(_pbend[0],_pbend[1]);
		Point pfind(0,0);
		Storage::real find[3];
		Storage::real div = (pabeg.x - paend.x)*(pbbeg.y - pbend.y) - (pabeg.y - paend.y)*(pbbeg.x - pbend.x), t1,t2;
		if (fabs(div) < 1.0e-13)
		{
			if (print) std::cout << "divisor is zero" << std::endl;
			return std::make_pair(false,InvalidNode());
		}
		pfind.x = ((pabeg.x*paend.y - pabeg.y*paend.x)*(pbbeg.x - pbend.x) - (pabeg.x - paend.x)*(pbbeg.x*pbend.y - pbbeg.y*pbend.x)) / div;
		pfind.y = ((pabeg.x*paend.y - pabeg.y*paend.x)*(pbbeg.y - pbend.y) - (pabeg.y - paend.y)*(pbbeg.x*pbend.y - pbbeg.y*pbend.x)) / div;
		//optimization uses information that we stay in unit cube
		if( pfind.x < 0 || pfind.x > 1 || pfind.y < 0 || pfind.y > 1 )
			return std::make_pair(false,InvalidNode());
		if (print) std::cout << "found ("<< pfind.x << ", " << pfind.y << ")" << std::endl;
		//probably some of these tests are redundant
		if (fabs(paend.x - pabeg.x) > eps)
		{
			t1 = (pfind.x - pabeg.x) / (paend.x - pabeg.x);
			if (t1 < eps || t1 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t1 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		if (fabs(paend.y - pabeg.y) > eps)
		{
			t1 = (pfind.y - pabeg.y) / (paend.y - pabeg.y);
			if (t1 < eps || t1 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t1 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		if (fabs(pbend.x - pbbeg.x) > eps)
		{
			t2 = (pfind.x - pbbeg.x) / (pbend.x - pbbeg.x);
			if (t2 < eps || t2 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t2 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		if (fabs(pbend.y - pbbeg.y) > eps)
		{
			t2 = (pfind.y - pbbeg.y) / (pbend.y - pbbeg.y);
			if (t2 < eps || t2 > 1.0 - eps)  { if (print) std::cout << "out of bound: " << t2 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		if (print) std::cout << "intersection accepted (" << pfind.x << "," << pfind.y << ") t1 " << t1 << " t2 " << t2 << std::endl;
		Node I;
		std::map<Point,Node>::iterator search = intersections.find(pfind);
		//check whether intersection already exists
		if( search != intersections.end() ) //there is a node!
			I = search->second;
		else //no node, create one
		{
			//restore coordinate
			//for(int k = 0; k < 3; ++k)
			//	find[k] = 0.5*((1-t1)*abeg[k]+t1*aend[k] + (1-t2)*bbeg[k]+t2*bend[k]);
			unp.Act(pfind,find);
			I = m->CreateNode(find);
			Storage::real_array _pfind = I->RealArray(pnt);
			_pfind[0] = pfind.x;
			_pfind[1] = pfind.y;
			std::pair<std::map<Point,Node>::iterator,bool> ins = intersections.insert(std::make_pair(pfind,I));
			if( !ins.second ) std::cout << "intersection node " << I->GetHandle() << " was not inserted at " << pfind.x << "," << pfind.y << " since there " << ins.first->second->GetHandle() << std::endl;
		}
		if( I == a->getBeg() || I == b->getBeg() || I == a->getEnd() || I == b->getEnd())
			return std::make_pair(false,I);
		return std::make_pair(true,I);
	}

	 
	void check_multimap(std::multimap<std::pair<Storage::real,int>, int,event_less> & events, const std::multimap<std::pair<Storage::real,int>, int,event_less>::iterator & it, bool print)
	{
		if( print )
		{
			{
				std::multimap<std::pair<Storage::real,int>, int,event_less>::iterator jt = it;
				++jt;
				if( jt != events.end() )
				{
					std::cout << "inserted " << std::scientific << it->first.first << " " << it->first.second << " segment " << it->second << std::endl;
					std::cout << "next     " << std::scientific << jt->first.first << " " << jt->first.second << " segment " << jt->second << std::endl;
					std::cout << "compare " << event_less()(it->first,jt->first) << std::endl;
				}
			}
			if( it != events.begin() )
			{
				std::multimap<std::pair<Storage::real,int>, int,event_less>::iterator  jt = it;
				--jt;
				std::cout << "prev     " << std::scientific << jt->first.first << " " << jt->first.second << " segment " << jt->second << std::endl;
				std::cout << "inserted " << std::scientific << it->first.first << " " << it->first.second << " segment " << it->second << std::endl;
				std::cout << "compare " << event_less()(jt->first,it->first) << std::endl;
			}
		}
	}

	void split_edges(Mesh * m, Node I, Edge a, Edge b, ElementArray<Edge> & splitted_a, ElementArray<Edge> & splitted_b,std::vector<Tag> & transfer)
	{
		//storage for data
		std::vector< std::vector<char> > copy(transfer.size()*2);
		//memorize data
		{
			const Edge * s[2] = {&a,&b};
			for(int q = 0; q < 2; ++q) //segments
				for(int k = 0; k < transfer.size(); ++k) //tags
				{
					int size = s[q]->GetDataSize(transfer[k]);
					copy[k + q*transfer.size()].resize(transfer[k].GetBytesSize()*size);
					if( !copy.empty() ) s[q]->GetData(transfer[k],0,size,&copy[k + q*transfer.size()][0]);
				}
		}
		splitted_a = Edge::SplitEdge(a,ElementArray<Node>(m,1,I->GetHandle()),0);
		splitted_b = Edge::SplitEdge(b,ElementArray<Node>(m,1,I->GetHandle()),0);
		//duplicate data
		{
			const Edge splitted[2][2] =
			{
				{splitted_a[0],splitted_a[1]},
				{splitted_b[0],splitted_b[1]}
			};
			for(int q = 0; q < 2;++q) //segments
				for(int k = 0; k < transfer.size();++k)
				{
					int size = (int)copy[k + q*transfer.size()].size() / transfer[k].GetBytesSize();
					if( size ) for(int l = 0; l < 2; ++l) //two parts
					{
						splitted[q][l].SetDataSize(transfer[k],size);
						splitted[q][l].SetData(transfer[k],0,size,&copy[k + q*transfer.size()][0]);
					}
				}
		}
	}
	


	void intersect_naive(Mesh * m, ElementArray<Edge> & segments, ElementArray<Node> & nodes, std::vector<Tag> & transfer, Tag pnt, const Unproject & unp, bool print)
	{
		//Tag pnt = m->CreateTag("PROJ_PNT"+m->GetLocalProcessorRank(),DATA_REAL,NODE,NODE,2);
		std::map<Point,Node> intersections;
		std::vector<HandleType> initials(segments.size()*2);
		MarkerType initial = m->CreatePrivateMarker();
		for (int k = 0; k < (int)segments.size(); ++k)
		{
			initials[k*2+0] = segments[k]->getBeg()->GetHandle();
			initials[k*2+1] = segments[k]->getEnd()->GetHandle();
			segments[k]->getBeg()->SetPrivateMarker(initial);
			segments[k]->getEnd()->SetPrivateMarker(initial);
			Point pbeg(segments[k]->getBeg()->RealArray(pnt).data());
			Point pend(segments[k]->getEnd()->RealArray(pnt).data());
			intersections.insert(std::make_pair(pbeg,segments[k]->getBeg()));
			intersections.insert(std::make_pair(pend,segments[k]->getEnd()));
		}
		for(int i = 0; i < (int)segments.size(); ++i)
		{
			for(int j = i+1; j < (int)segments.size(); ++j)
			{
				std::pair<bool,Node> I = intersect_segments(m,segments[i],segments[j],intersections,pnt,unp,print);
				if( I.first )
				{
					ElementArray<Edge> splitted_a, splitted_b;
					split_edges(m,I.second,segments[i],segments[j],splitted_a,splitted_b,transfer);
					segments[i] = splitted_a[0];
					segments[j] = splitted_b[0];
					segments.push_back(splitted_a[1]);
					segments.push_back(splitted_b[1]);
				}
			}
		}
		nodes.clear();
		for(std::map<Point,Node>::iterator it = intersections.begin(); it != intersections.end(); ++it)
		{
			if( !it->second->GetPrivateMarker(initial) )
				nodes.push_back(it->second);
		}
		for(int k = 0; k < (int)initials.size(); ++k)
			m->RemPrivateMarker(initials[k],initial);
		m->ReleasePrivateMarker(initial);
	}

	//account for intersection event
	void intersect_event(Mesh * m, int a, int b, Node I, ElementArray<Edge> & segments, std::multimap<Point, int> & sweep, std::multimap<std::pair<Storage::real,int>, int,event_less> & events, std::vector<Tag> & transfer, Tag pnt, bool print)
	{
		const bool checkmm = false;
		//remove event of ending of old segment
		{
			int rem_end_events[2];
			rem_end_events[0] = a;
			rem_end_events[1] = b;
			for (int k = 0; k < 2; ++k)
			{
				std::pair< std::multimap<std::pair<Storage::real,int>, int,event_less>::iterator, std::multimap<std::pair<Storage::real,int>,int,event_less>::iterator > del = events.equal_range(std::make_pair(segments[rem_end_events[k]]->getEnd()->RealArray(pnt)[0],SEG_END)); //get all events at position of the end
				bool flag = false;
				for (std::multimap<std::pair<Storage::real,int>, int,event_less>::iterator it = del.first; it != del.second; ++it) //search over all events
				{
					if (it->first.second == SEG_END && it->second == rem_end_events[k]) //event is end of segment and segment matches current
					{
						events.erase(it); //remove that segment
						flag = true;
						break; //do not expect any more
					}
				}
				if (!flag)
				{
					std::cout << "Cannot find proper ending event for segment " << rem_end_events[k] << std::endl;
					std::cout << "Found events:" << std::endl;
					for (std::multimap<std::pair<Storage::real,int>, int,event_less>::iterator it = del.first; it != del.second; ++it)
						std::cout << "x: " << it->first.first << (it->first.second == SEG_END ? "end ":"start ") << " segment " << it->second << std::endl;
					std::cout << "Was looking at " << segments[rem_end_events[k]]->getEnd()->Coords()[2] << " " << segments[rem_end_events[k]]->getEnd()->GetHandle() << std::endl;
					std::cout << "Other side at " << segments[rem_end_events[k]]->getBeg()->Coords()[2] << " " << segments[rem_end_events[k]]->getBeg()->GetHandle() << std::endl;
					std::cout << "All events: " << events.size() << std::endl;
					for (std::multimap<std::pair<Storage::real, int>, int,event_less>::iterator it = events.begin(); it != events.end(); ++it)
						std::cout << "x: " << it->first.first << " type " << (it->first.second == SEG_START ? "start" : "end") << " segment " << it->second << " " << segments[it->second]->GetHandle() << " " << segments[it->second]->getBeg()->GetHandle() << " " << segments[it->second]->getEnd()->GetHandle() <<std::endl;
					
				}
				assert(flag);
			}
		}
		ElementArray<Edge> splitted_a, splitted_b;
		split_edges(m,I,segments[a],segments[b],splitted_a,splitted_b,transfer);
		//replace segment a by new one
		segments[a] = splitted_a[0];
		//add event of ending of old segment
		if(print) std::cout << "1: Add segment " << a << " " << segments[a]->GetHandle() << " end" << " at " << I->RealArray(pnt)[0] << " " << I->GetHandle() << std::endl;
		check_multimap(events,events.insert(std::make_pair(std::make_pair(I->RealArray(pnt)[0],SEG_END), a)),checkmm);
		//put other side of segment
		segments.push_back(splitted_a[1]);
		//detect proper starting event
		if( segments.back()->getBeg()->RealArray(pnt)[0] > segments.back()->getEnd()->RealArray(pnt)[0] )
			segments.back()->SwapEnds();
		//add event of starting of new segment
		if(print) std::cout << "2: Add segment " << segments.size() - 1 << " " << segments.back()->GetHandle() << " start at " << segments.back()->getBeg()->RealArray(pnt)[0] << " " << segments.back()->getBeg()->GetHandle() << std::endl;
		check_multimap(events,events.insert(std::make_pair(std::make_pair(segments.back()->getBeg()->RealArray(pnt)[0],SEG_START), (int)segments.size() - 1)),checkmm);
		//add event of ending of new segment
		if(print) std::cout << "3: Add segment " << segments.size() - 1 << " " << segments.back()->GetHandle() << " end at " << segments.back()->getEnd()->RealArray(pnt)[0] << " " << segments.back()->getEnd()->GetHandle() << std::endl;
		check_multimap(events,events.insert(std::make_pair(std::make_pair(segments.back()->getEnd()->RealArray(pnt)[0],SEG_END),(int)segments.size() - 1)),checkmm);
		//replace segment b by new one
		segments[b] = splitted_b[0];
		//add event of ending of old segment
		if(print) std::cout << "4: Add segment " << b << " " << segments[b]->GetHandle() << " end at " << I->RealArray(pnt)[0] << " " << I->GetHandle() << std::endl;
		check_multimap(events,events.insert(std::make_pair(std::make_pair(I->RealArray(pnt)[0],SEG_END), b)),checkmm);
		//put other side of segment
		segments.push_back(splitted_b[1]);
		//detect proper starting event
		if( segments.back()->getBeg()->RealArray(pnt)[0] > segments.back()->getEnd()->RealArray(pnt)[0] )
			segments.back()->SwapEnds();
		//add event of starting of new segment
		if(print) std::cout << "5: Add segment " << segments.size() - 1 << " " << segments.back()->GetHandle() << " start at " << segments.back()->getBeg()->RealArray(pnt)[0]<< " " << segments.back()->getBeg()->GetHandle() << std::endl;
		check_multimap(events,events.insert(std::make_pair(std::make_pair(segments.back()->getBeg()->RealArray(pnt)[0],SEG_START), (int)segments.size() - 1)),checkmm);
		//add event of ending of new segment
		if(print) std::cout << "6: Add segment " << segments.size() - 1 << " " << segments.back()->GetHandle() << " end at " << segments.back()->getEnd()->RealArray(pnt)[0] << " " << segments.back()->getEnd()->GetHandle() << std::endl;
		check_multimap(events,events.insert(std::make_pair(std::make_pair(segments.back()->getEnd()->RealArray(pnt)[0],SEG_END), (int)segments.size() - 1)),checkmm);
		if (print)
		{
			std::cout << "Number of events: " << events.size() << std::endl;
			for (std::multimap<std::pair<Storage::real, int>, int,event_less>::iterator it = events.begin(); it != events.end(); ++it)
				std::cout << "x: " << it->first.first << " type " << (it->first.second == SEG_START ? "start" : "end") << " segment " << it->second << " " << segments[it->second]->GetHandle() << " " << segments[it->second]->getBeg()->GetHandle() << " " << segments[it->second]->getEnd()->GetHandle() << std::endl;
			//assert(count_duplicates(segments) == 0);
		}
	}

	//intersect many segments
	void intersect(Mesh * m, ElementArray<Edge> & segments, ElementArray<Node> & nodes, std::vector<Tag> & transfer, Tag pnt, const Unproject & unp, bool print)
	{
		std::map<Point,Node> intersections;
		std::multimap<std::pair<Storage::real,int>,int,event_less> events;
		std::multimap<Point,int> sweep;
		
		MarkerType initial = m->CreatePrivateMarker();
	
		if( print )
		{
			std::cout << "Input segments[" << segments.size() << "]: " << std::endl;
			for (ElementArray<Edge>::iterator it = segments.begin(); it != segments.end(); ++it)
				std::cout << "[ (" <<  it->getBeg()->Coords()[0] << "," << it->getBeg()->Coords()[1] << "," << it->getBeg()->Coords()[2] << "), ("<<  it->getEnd()->Coords()[0]   << "," << it->getEnd()->Coords()[1] << "," << it->getEnd()->Coords()[2] << ") ] " << std::endl;
			std::cout << "Create events based on segments." << std::endl;
		}

		for (int k = 0; k < (int)segments.size(); ++k)
		{
			if (segments[k]->getBeg()->RealArray(pnt)[0] > segments[k]->getEnd()->RealArray(pnt)[0])
				segments[k].SwapEnds();
			events.insert(std::make_pair(std::make_pair(segments[k]->getBeg()->RealArray(pnt)[0],SEG_START),k));
			events.insert(std::make_pair(std::make_pair(segments[k]->getEnd()->RealArray(pnt)[0],SEG_END), k));
			segments[k]->getBeg()->SetPrivateMarker(initial);
			segments[k]->getEnd()->SetPrivateMarker(initial);
			intersections.insert(std::make_pair(Point(segments[k]->getBeg()->RealArray(pnt).data()),segments[k]->getBeg()));
			intersections.insert(std::make_pair(Point(segments[k]->getEnd()->RealArray(pnt).data()),segments[k]->getEnd()));
		}


		if (print)
		{
			std::cout << "Number of events: " << events.size() << std::endl;
			for (std::multimap<std::pair<Storage::real, int>, int,event_less>::iterator it = events.begin(); it != events.end(); ++it)
				std::cout << "x: " << it->first.first << " type " << (it->first.second == SEG_START ? "start" : "end") << " segment " << it->second << " " << segments[it->second]->GetHandle() << " " << segments[it->second]->getBeg()->GetHandle() << " " << segments[it->second]->getEnd()->GetHandle() << std::endl;
		
			std::cout << " Start parsing events" << std::endl;
		}
	
		while (!events.empty())
		{
			std::multimap<std::pair<Storage::real,int>,int,event_less>::iterator first = events.begin();
			int t = first->first.second;
			int s = first->second;
			events.erase(first);
			if (t == SEG_START)
			{
				if( print ) std::cout << "Segment " << s << " start" << std::endl;
				//check if there is a line with same position
				Point p(segments[s]->getBeg()->RealArray(pnt).data());
				std::multimap<Point, int>::iterator ins = sweep.insert(std::make_pair(p, s));
				if (print)
				{
					std::cout << "Inserted into sweep" << std::endl;
					for (std::multimap<Point, int>::iterator it = sweep.begin(); it != sweep.end(); ++it)
						std::cout << "(" << it->first.x << "," << it->first.y << ")" << " segment " << it->second << std::endl;
				}
				//check line (or lines above current)
				for (int dir = 0; dir <= 1; ++dir) // look up or down
				{
					if( print ) std::cout << "Looking " << (dir ? "up" : "down") << std::endl;
					std::multimap<Point, int>::iterator iter = ins;
					while ((dir ? ++iter != sweep.end() : iter != sweep.begin())) //y is greater for next
					{
						if( !dir ) --iter;
						if (print) std::cout << "test " << s << " with " << iter->second << std::endl;
						if (segments[s]->getBeg() != segments[iter->second]->getBeg()) //ignore same starting position
						{
							if (print) std::cout << "checking intersection" << std::endl;
							std::pair<bool,Node> I = intersect_segments(m,segments[s], segments[iter->second],intersections,pnt,unp,print);
							if (I.first)
							{
								if( print ) std::cout << "Intersection of " << s << " " <<segments[s]->GetHandle() << " " << segments[s]->getBeg()->GetHandle() << " " << segments[s]->getEnd()->GetHandle() << " and " << iter->second << " " << segments[iter->second]->GetHandle() << " " << segments[iter->second]->getBeg()->GetHandle() << " " << segments[iter->second]->getEnd()->GetHandle() << " at (" << I.second.Coords()[0] << "," << I.second.Coords()[1] << "," << I.second.Coords()[2] << ") " << std::endl;
								intersect_event(m,s, iter->second, I.second, segments, sweep, events,transfer, pnt, print);
								//break;
							}
						}
						else if (print) std::cout << "skipping segments with same starting point" << std::endl;
						if ((2*dir-1)*(iter->first.y - ins->first.y) > 0) //visited line is above (below) current
							break; //stop search
					}
				}
			}
			else if (t == SEG_END)
			{
				if( print ) std::cout << "Segment " << s << " end" << std::endl;
				//remove segment from sweep
				Point p(segments[s]->getBeg()->RealArray(pnt).data());
				std::pair< std::multimap<Point, int>::iterator, std::multimap<Point, int>::iterator > range = sweep.equal_range(p);
				if( print ) std::cout << "Range distance " << std::distance(range.first,range.second) << " sweep size " << sweep.size() << std::endl;
				std::multimap<Point, int>::iterator above = range.second, below = range.first;
				bool flag = false, test = true;
				if( below == sweep.begin() ) test = false;
				else --below;
				if( above == sweep.end() ) test = false;
				if( test && print ) std::cout << "Test will be performed" << std::endl;
				for (std::multimap<Point, int>::iterator it = range.first; it != range.second; ++it) //search over all events
				{
					if( it->second == s) //found necessery segment
					{
						if (print)
						{
							std::cout << "Erase segment " << s << " from sweep: " << std::endl;
							for (std::multimap<Point, int>::iterator it = sweep.begin(); it != sweep.end(); ++it)
								std::cout << "(" << it->first.x << "," << it->first.y << ")" << " segment " << it->second << std::endl;
						}
						sweep.erase(it);
						flag = true;
						break; //do not expect any more
					}
				}
				if (!flag) std::cout << __FILE__ << ":" << __LINE__ <<  " Error: cannot find segment " << s << " in sweep" << std::endl;
				assert(flag);
				if (test)
				{
					if (print) std::cout << "test " << below->second << " with " << above->second << std::endl;
					if (segments[above->second]->getBeg() != segments[below->second]->getBeg())
					{
						if (print) std::cout << "checking intersection" << std::endl;
						std::pair<bool,Node> I = intersect_segments(m, segments[below->second], segments[above->second],intersections,pnt,unp,print);
						if (I.first)
						{
							if( print ) std::cout << "Intersection of " << below->second << " " << segments[below->second]->GetHandle() << " " << segments[below->second]->getBeg()->GetHandle() << " " << segments[below->second]->getEnd()->GetHandle() << " and " << above->second << " " << segments[above->second]->GetHandle() << " " << segments[above->second]->getBeg()->GetHandle() << " " << segments[above->second]->getEnd()->GetHandle() << " at (" << I.second.Coords()[0] << "," << I.second.Coords()[1] << "," << I.second.Coords()[2] << ") " << std::endl;
							intersect_event(m,below->second, above->second, I.second, segments, sweep, events,transfer,pnt, print);
						}
					}
					else if (print) std::cout << "skipping segments with same starting point" << std::endl;
				}
			}
		}
		//copy intersections
		nodes.clear();
		for(std::map<Point,Node>::iterator it = intersections.begin(); it != intersections.end(); ++it)
		{
			if( !it->second->GetPrivateMarker(initial) )
				nodes.push_back(it->second);
			else it->second->RemPrivateMarker(initial);
		}
		m->ReleasePrivateMarker(initial);
	}
	
	void block_number_intersection(ElementArray<Element> & adj, Tag block, std::vector<int> & out)
	{
		out.clear();
		if( !adj.empty() )
		{
			Storage::integer_array be = adj[0]->IntegerArray(block);
			std::vector<int> inter(be.begin(),be.end()), tmp(inter.size());
			for(ElementArray<Edge>::size_type k = 1; k < adj.size(); ++k)
			{
				be = adj[k]->IntegerArray(block);
				tmp.resize(std::set_intersection(inter.begin(),inter.end(),be.begin(),be.end(),tmp.begin())-tmp.begin());
				inter.swap(tmp);
			}
			out.insert(out.end(),inter.begin(),inter.end());
		}
		//Storage::integer_array bn = n->IntegerArray(block);
		//bn.replace(bn.begin(),bn.end(),inter.begin(),inter.end());
	}

	void block_number_union(Element n, ElementArray<Element> & adj, Tag block)
	{
		if( !adj.empty() )
		{
			Storage::integer_array be = adj[0]->IntegerArray(block);
			std::vector<int> uni(be.begin(),be.end()), tmp(uni.size());
			for(ElementArray<Edge>::size_type k = 1; k < adj.size(); ++k)
			{
				be = adj[k]->IntegerArray(block);
				tmp.resize(uni.size()+be.size());
				tmp.resize(std::set_union(uni.begin(),uni.end(),be.begin(),be.end(),tmp.begin())-tmp.begin());
				uni.swap(tmp);
			}
			Storage::integer_array bn = n->IntegerArray(block);
			bn.replace(bn.begin(),bn.end(),uni.begin(),uni.end());
		}
	}
	

	void Mesh::LoadECL(std::string File)
	{
		std::cout << std::scientific;
		FILE * f = fopen(File.c_str(),"r");
		if( f == NULL )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << File << std::endl;
			throw BadFileName;
		}
		std::vector<HandleType> old_nodes(NumberOfNodes());
		{
			unsigned qq = 0;
			for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it)
				old_nodes[qq++] = *it;
		}
		if( !old_nodes.empty() ) 
			std::sort(old_nodes.begin(),old_nodes.end(),CentroidComparator(this));
		std::vector< std::pair< std::pair<FILE *,std::string>, int> > fs(1,std::make_pair(std::make_pair(f,File),0));
		char readline[2048], *p, *pend, rec[2048];
		int text_end, text_start, state = ECL_NONE, nchars;
		int waitlines = 0;
		int have_dimens = 0, totread, downread, numrecs, offset;
		int gtype = ECL_GTYPE_NONE;
		int argtype = ECL_VAR_NONE;
		int radial = ECL_GTYPE_NONE;
		Storage::real * read_arrayf = NULL;
		Storage::integer * read_arrayi = NULL;
		Storage::integer dims[3], mapaxis[6] = {0,1,0,0,1,0};
		Storage::real inrad = 0;
		std::vector<Storage::real> xyz,perm,poro, tops,zcorn;
		std::vector<Storage::integer> actnum, satnum;
		while(!fs.empty())
		{
			while(fgets(readline,2048,fs.back().first.first) != NULL)
			{
				fs.back().second++; //line number
				{
					if( readline[strlen(readline)-1] == '\n' ) readline[strlen(readline)-1] = '\0';
					text_end = static_cast<int>(strlen(readline));
					for(text_start = 0; isspace(readline[text_start]) && text_start < text_end; text_start++);
					if( text_start == text_end ) continue;
					for(text_end = text_end-1; isspace(readline[text_end]) && text_end > text_start; text_end--);
					readline[text_end+1] = '\0';
					p = readline + text_start;
					pend = readline + text_end + 1;
					for(char * q = p; q < pend; q++) *q = toupper(*q);
				}
				if( p[0] == '-' && p[1] == '-' ) continue; //skip comment
				if(waitlines) {waitlines--; continue;} //skip meaningful lines
				switch(state)
				{
				case ECL_NONE:
					if( !ECLSTRCMP(p,"END") ) //end of data - don't read beyond
					{
						goto ecl_exit_loop;
					}
					else if( !ECLSTRCMP(p,"INCLUDE") ) state = ECL_INCLUDE;
					else if( !ECLSTRCMP(p,"DIMENS") || !ECLSTRCMP(p,"SPECGRID") )
					{
						read_arrayi = dims;
						numrecs = 1;
						downread = totread = 3;
						argtype = ECL_VAR_INT;
						offset = state = ECL_DIMENS;
						have_dimens = 1;
					}
					else if( !ECLSTRCMP(p,"MAPAXIS") )
					{
						read_arrayi = mapaxis;
						numrecs = 1;
						downread = totread = 6;
						argtype = ECL_VAR_INT;
						offset = state = ECL_MAPAXIS;
					}
					else if( !ECLSTRCMP(p,"DX") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = xyz.empty()? NULL : &xyz[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_DX;
					}
					else if( !ECLSTRCMP(p,"DY") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_DX;
						state = ECL_DY;
					}
					else if( !ECLSTRCMP(p,"DZ") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_DX;
						state = ECL_DZ;
					}
					else if( !ECLSTRCMP(p,"COORD") )
					{
						assert(have_dimens);
						if( xyz.empty() ) xyz.resize(3*2*(dims[0]+1)*(dims[1]+1));
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
						numrecs = 1;
						downread = totread = 3*2*(dims[0]+1)*(dims[1]+1);
						argtype = ECL_VAR_REAL;
						offset = state = ECL_COORDS;
						gtype = ECL_GTYPE_ZCORN;
					}
					else if( !ECLSTRCMP(p,"ZCORN") )
					{
						assert(have_dimens);
						if( zcorn.empty() ) zcorn.resize(dims[0]*dims[1]*dims[2]*8);
						read_arrayf = zcorn.empty() ? NULL : &zcorn[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2]*8;
						argtype = ECL_VAR_REAL;
						state = offset = ECL_ZCORN;
					}
					else if( !ECLSTRCMP(p,"TOPS") )
					{
						assert(have_dimens);
						tops.resize(dims[0]*dims[1]*dims[2]);
						read_arrayf = tops.empty() ? NULL : &tops[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_TOPS;
						gtype = ECL_GTYPE_TOPS;
					}
					else if( !ECLSTRCMP(p,"PERMX") )
					{
						assert(have_dimens);
						if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = perm.empty() ? NULL : &perm[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_PERMX;
					}
					else if( !ECLSTRCMP(p,"PERMY") )
					{
						assert(have_dimens);
						if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = perm.empty() ? NULL : &perm[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_PERMX;
						state = ECL_PERMY;
					}
					else if( !ECLSTRCMP(p,"PERMZ") )
					{
						assert(have_dimens);
						if( perm.empty() ) perm.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = perm.empty() ? NULL : &perm[0];
						numrecs = 3;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = ECL_PERMX;
						state = ECL_PERMZ;
					}
					else if( !ECLSTRCMP(p,"PORO") )
					{
						assert(have_dimens);
						poro.resize(3*dims[0]*dims[1]*dims[2]);
						read_arrayf = poro.empty() ? NULL : &poro[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_REAL;
						offset = state = ECL_PORO;
					}
					else if( !ECLSTRCMP(p,"ACTNUM") )
					{
						assert(have_dimens);
						actnum.resize(dims[0]*dims[1]*dims[2]);
						read_arrayi = actnum.empty() ? NULL : &actnum[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_INT;
						offset = state = ECL_ACTNUM;
					}
					else if( !ECLSTRCMP(p,"SATNUM") )
					{
						assert(have_dimens);
						satnum.resize(dims[0]*dims[1]*dims[2]);
						read_arrayi = satnum.empty() ? NULL : &satnum[0];
						numrecs = 1;
						downread = totread = dims[0]*dims[1]*dims[2];
						argtype = ECL_VAR_INT;
						offset = state = ECL_SATNUM;
					}
					else if( !ECLSTRCMP(p,"RADIAL") )
					{
						radial = ECL_GTYPE_RADIAL;
					}
					else if( !ECLSTRCMP(p,"CART") )
					{
						radial = ECL_GTYPE_CARTESIAN;
					}
					else if( !ECLSTRCMP(p,"INRAD") )
					{
						if( radial != ECL_GTYPE_RADIAL ) 
						{
							std::cout << __FILE__ << ":" << __LINE__ << " inner radius specified for cartesian grid ";
							std::cout << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
						}
						if( radial == ECL_GTYPE_NONE ) radial = ECL_GTYPE_RADIAL;
						state = ECL_INRAD;
					}
					else
					{
						//std::cout << __FILE__ << ":" << __LINE__ << " skipped " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
					}
					break;
				case ECL_SKIP_SECTION:
					if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE;
					break;
				case ECL_INCLUDE:
					if( 1 == sscanf(p," %s",rec) )
					{
						int shift_one = 0;
						if( (rec[0] == '\'' || rec[0] == '"') && rec[0] == rec[strlen(rec)-1] ) //remove quotes
						{
							rec[strlen(rec)-1] = '\0';
							shift_one = 1;
						}
						f = fopen((GetFolder(fs.back().first.second) + "/" + std::string(rec+shift_one)).c_str(),"r");
						if( f == NULL )
						{
							std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << (GetFolder(fs.back().first.second) + "/" + std::string(rec+shift_one)) << " included from "<< fs.back().first.second << " line " << fs.back().second << std::endl;
							throw BadFileName;
						}
						fs.push_back(std::make_pair(std::make_pair(f,GetFolder(fs.back().first.second) + "/"+std::string(rec+shift_one)),0));
						if( *(pend-1) == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read file name, string " << p << " in " << fs.back().first.second << " line " << fs.back().second << std::endl;
						throw BadFile;
					}
					break;
				case ECL_ZCORN:
				case ECL_COORDS:
				case ECL_MAPAXIS:
				case ECL_DIMENS:
				case ECL_DX:
				case ECL_DY:
				case ECL_DZ:
				case ECL_PERMX:
				case ECL_PERMY:
				case ECL_PERMZ:
				case ECL_PORO:
				case ECL_ACTNUM:
				case ECL_SATNUM:
				case ECL_TOPS:
					while( downread > 0 && p < pend )
					{
						if( 1 == sscanf(p,"%s%n",rec,&nchars) )
						{
							p += nchars;
							while(isspace(*p) && p < pend) ++p;
							int count = 1;
							int begval = 0;
							for(int q = 0; q < static_cast<int>(strlen(rec)); ++q)
								if( rec[q] == '*' )
								{
									begval = q+1;
									rec[q] = '\0';
									break;
								}
							if( begval > 0 ) count = atoi(rec);
							if( argtype == ECL_VAR_REAL )
							{
								Storage::real val = atof(rec+begval);
								while(count)
								{
									read_arrayf[numrecs*(totread-(downread--))+(state-offset)] = val;
									count--;
								}
							}
							else if( argtype == ECL_VAR_INT )
							{
								Storage::integer val = atoi(rec+begval);
								while(count)
								{
									read_arrayi[numrecs*(totread-(downread--))+(state-offset)] = val;
									count--;
								}
							}
							else
							{
								std::cout << __FILE__ << ":" << __LINE__ << " probably forgot to set up argument type to read " << std::endl;
								throw Impossible;
							}
						}
						else
						{
							std::cout << __FILE__ << ":" << __LINE__ << " cannot read data " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
							throw BadFile;
						}
						if( *p == '/' ) 
						{
							if( downread > 0 )
							{
								std::cout << __FILE__ << ":" << __LINE__ << " early data termination, read " << totread-downread << " of " << totread << " records in " << fs.back().first.second << ":" << fs.back().second << std::endl;
								throw BadFile;
							}
						}
					}
					if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE; else if( downread == 0 ) state = ECL_SKIP_SECTION;
					break;
				case ECL_INRAD:
					if( 1 == sscanf(p,"%lf%n",&inrad,&nchars) )
					{
						p += nchars;
						while(isspace(*p) && p < pend) ++p;
					}
					else
					{
						std::cout << __FILE__ << ":" << __LINE__ << " cannot read data " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
						throw BadFile;
					}
					if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
					break;
				}
			}
ecl_exit_loop:
			fclose(fs.back().first.first);
			fs.pop_back();
		}
		if( radial == ECL_GTYPE_RADIAL )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " radial grids not supported yet " << std::endl;
		}
		if( gtype == ECL_GTYPE_TOPS )
		{
			
			std::vector<HandleType> newnodes((dims[0]+1)*(dims[1]+1)*(dims[2]+1));
			Storage::real x, y, z, node_xyz[3];
			x = 0.0;
			int numnode = 0;
			for(int i = 0; i < dims[0]+1; i++)
			{
				Storage::integer pif = std::min(dims[0]-1,i), pib = std::max(i-1,0);
				y = 0.0;
				for(int j = 0; j < dims[1]+1; j++)
				{
					Storage::integer pjf = std::min(dims[1]-1,j), pjb = std::max(j-1,0);
					z = (
							tops[ECL_IJK_DATA(pib,pjb,0)]+
							tops[ECL_IJK_DATA(pib,pjf,0)]+
							tops[ECL_IJK_DATA(pif,pjb,0)]+
							tops[ECL_IJK_DATA(pif,pjf,0)]
						  )*0.25;
					z -= (
							xyz[3*ECL_IJK_DATA(pib,pjb,0)+2]+
							xyz[3*ECL_IJK_DATA(pib,pjf,0)+2]+
							xyz[3*ECL_IJK_DATA(pif,pjb,0)+2]+
							xyz[3*ECL_IJK_DATA(pif,pjf,0)+2]
							)*0.25;
					for(int k = 0; k < dims[2]+1; k++)
					{
						Storage::integer pkf = std::min(dims[2]-1,k), pkb = std::max(k-1,0);
						node_xyz[0] = x;
						node_xyz[1] = y;
						node_xyz[2] = z;
						int find = -1;
						if( !old_nodes.empty() )
						{
							std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),node_xyz,CentroidComparator(this));
							if( it != old_nodes.end() ) 
							{
								Storage::real_array c = RealArrayDF(*it,CoordsTag());
								if( CentroidComparator(this).Compare(node_xyz,c.data()) == 0 )
									find = static_cast<int>(it - old_nodes.begin());
							}
						}
						if( find == -1 ) newnodes[numnode++] = CreateNode(node_xyz)->GetHandle();
						else newnodes[numnode++] = old_nodes[find];
						//std::cout << i << " " << j << " " << k << " ( " << x << " , " << y << " , " << z << ") " << newnodes.back()->LocalID() << std::endl; 
						x += (
								(
								  xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+0]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+0]
								)
							-
								(
								  xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+0]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+0]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+0]
								)
								)*0.25;
						y += (
								(
							      xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+1]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+1]
								)
							-
								(
								  xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+1]+
								  xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+1]+
								  xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+1]
								)
								)*0.25;
						z += (
								xyz[3*ECL_IJK_DATA(pib,pjb,pkb)+2]+
								xyz[3*ECL_IJK_DATA(pib,pjf,pkb)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjb,pkb)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjf,pkb)+2]+
							    xyz[3*ECL_IJK_DATA(pib,pjb,pkf)+2]+
								xyz[3*ECL_IJK_DATA(pib,pjf,pkf)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjb,pkf)+2]+
								xyz[3*ECL_IJK_DATA(pif,pjf,pkf)+2]
								)*0.125;
					}
					y += (
							xyz[3*ECL_IJK_DATA(pib,pjb,0)+1]+
							xyz[3*ECL_IJK_DATA(pif,pjb,0)+1]+
							xyz[3*ECL_IJK_DATA(pib,pjf,0)+1]+
							xyz[3*ECL_IJK_DATA(pif,pjf,0)+1]
							)*0.25; 
				}
				x += (
						xyz[3*ECL_IJK_DATA(pib,0,0)+0]+
						xyz[3*ECL_IJK_DATA(pif,0,0)+0]
						)*0.5; 
			}
			Tag tagporo,tagperm, tagsatnum;
			if( !poro.empty() ) tagporo = CreateTag("PORO",DATA_REAL,CELL,NONE,1);
			if( !perm.empty() ) tagperm = CreateTag("PERM",DATA_REAL,CELL,NONE,3);
			if( !satnum.empty() ) tagsatnum = CreateTag("SATNUM",DATA_INTEGER,CELL,NONE,1);
			const Storage::integer nvf[24] = { 2, 3, 1, 0, 4, 5, 7, 6, 0, 1, 5, 4, 3, 2, 6, 7, 2, 0, 4, 6, 1, 3, 7, 5 };
			const Storage::integer numnodes[6] = { 4, 4, 4, 4, 4, 4 };
			for(int i = 0; i < dims[0]; i++)
				for(int j = 0; j < dims[1]; j++)
					for(int k = 0; k < dims[2]; k++)
					{
						HandleType verts[8];
						verts[0] = newnodes[((i+0)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+0)];
						verts[1] = newnodes[((i+1)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+0)];
						verts[2] = newnodes[((i+0)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+0)];
						verts[3] = newnodes[((i+1)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+0)];
						verts[4] = newnodes[((i+0)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+1)];
						verts[5] = newnodes[((i+1)*(dims[1]+1)+(j+0))*(dims[2]+1)+(k+1)];
						verts[6] = newnodes[((i+0)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+1)];
						verts[7] = newnodes[((i+1)*(dims[1]+1)+(j+1))*(dims[2]+1)+(k+1)];
						//for(int q = 0; q < 8; q++)
						//	std::cout << verts[q]->Coords()[0] << " " << verts[q]->Coords()[1] << " " << verts[q]->Coords()[2] << " " << verts[q]->LocalID() << std::endl;
						Cell c = CreateCell(ElementArray<Node>(this,verts,verts+8),nvf,numnodes,6).first;
						if( !poro.empty() ) c->RealDF(tagporo) = poro[ECL_IJK_DATA(i,j,k)];
						if( !satnum.empty() ) c->IntegerDF(tagsatnum) = satnum[ECL_IJK_DATA(i,j,k)];
						if( !perm.empty() )
						{
							Storage::real_array arr_perm = c->RealArrayDF(tagperm);
							arr_perm[0] = perm[3*(ECL_IJK_DATA(i,j,k))+0];
							arr_perm[1] = perm[3*(ECL_IJK_DATA(i,j,k))+1];
							arr_perm[2] = perm[3*(ECL_IJK_DATA(i,j,k))+2];
						}
					}
		}
		else if( gtype == ECL_GTYPE_ZCORN )
		{
			//SetTopologyCheck(PRINT_NOTIFY | NEED_TEST_CLOSURE | PROHIBIT_MULTIPOLYGON | PROHIBIT_MULTILINE | MARK_ON_ERROR);
			//SetTopologyCheck(DEGENERATE_EDGE | DEGENERATE_FACE | DEGENERATE_CELL);
			//SetTopologyCheck(TRIPLE_SHARED_FACE | FLATTENED_CELL | INTERLEAVED_FACES);
			//SetTopologyCheck(DUPLICATE_EDGE | DUPLICATE_FACE | DUPLICATE_CELL);
			//SetTopologyCheck(ADJACENT_DUPLICATE | ADJACENT_DIMENSION);
			//RemTopologyCheck(THROW_EXCEPTION);
			//actnum.clear();
			if( zcorn.empty() )
			{
				std::cout << "ZCORN was not provided, cannot construct grid" << std::endl;
				throw BadFile;
			}
			if( xyz.empty() )
			{
				std::cout << "COORD was not provided, cannot construct grid" << std::endl;
				throw BadFile;
			}
			//make arrays more human-readable
			//std::vector< std::vector< std::vector< std::vector< Storage::real > > > > zcorn_array;//, coords_array;
			/*
			{
				zcorn_array.resize(dims[0]);
				for(int i = 0; i < dims[0]; i++)
				{
					zcorn_array[i].resize(dims[1]);
					for(int j = 0; j < dims[1]; j++)
					{
						zcorn_array[i][j].resize(dims[2]);
						for(int k = 0; k < dims[2]; k++)
							zcorn_array[i][j][k].resize(8);
					}
				}
				int pos = 0;
				for(int k = 0; k < dims[2]; k++)
				{
					for(int q = 0; q < 2; q++) //top-bottom
						for(int j = 0; j < dims[1]; j++)
							for(int m = 0; m < 2; m++) //near-far
								for(int i = 0; i < dims[0]; i++)
									for(int l = 0; l < 2; l++) //left-right
										zcorn_array[i][j][k][l+m*2+(1-q)*4] = zcorn[pos++];
				}
			}
			 */
			/*
			{
				coords_array.resize(dims[0]+1);
				for(int i = 0; i < dims[0]+1; i++)
				{
					coords_array[i].resize(dims[1]+1);
					for(int j = 0; j < dims[1]+1; j++)
					{
						coords_array[i][j].resize(2);
						for(int l = 0; l < 2; l++)
							coords_array[i][j][l].resize(3);
					}
				}
				int pos = 0;
				for(int j = 0; j < dims[1]+1; j++)
					for(int i = 0; i < dims[0]+1; i++)
						for(int l = 0; l < 2; l++)
							for(int k = 0; k < 3; k++)
								coords_array[i][j][l][k] = xyz[pos++];
			}
			 */
			//assemble pillars
			{
				//Tag node_number = CreateTag("NODE_NUMBER",DATA_INTEGER,NODE,NONE);
				Tag cell_number = CreateTag("CELL_NUMBER",DATA_INTEGER,CELL,NONE,1);
				Tag edge_number = CreateTag("EDGE_NUMBER",DATA_INTEGER,EDGE,NONE);
				Tag block_number = CreateTag("BLOCK_NUMBER",DATA_INTEGER,EDGE|NODE,NONE);
				typedef std::map<Storage::real,Node,pillar_less> pillar;
				std::vector< pillar > pillars((dims[0]+1)*(dims[1]+1));
				//this variant goes over pillars
				printf("create nodes on pillars\n");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(int i = 0; i < dims[0]+1; i++)
				{
					for(int j = 0; j < dims[1]+1; j++)
					{
						pillar & p = pillars[i*(dims[1]+1)+j];
						//loop over nodes on pillar
						for(int k = 0; k < dims[2]+1; k++)
						{
							Storage::real zpos[8], zalpha[8], z;
							int nzpos = 0;
							//loop over 8 blocks around node
							for(int l = 0; l < 8; ++l)
							{
								//block number
								int bi = i + (l%2) - 1;
								int bj = j + ((l/2)%2) - 1;
								int bk = k + (l/4) - 1;
								if( bi >= 0 && bj >= 0 && bk >= 0 && bi < dims[0] && bj < dims[1] && bk < dims[2] && (actnum.empty() || actnum[ECL_IJK_DATA(bi,bj,bk)]) )
								{
									bool found = false;
									z = zcorn[ECL_IJK_ZCORN(bi,bj,bk,7-l)];
									for(int q = 0; q < nzpos; ++q)
									{
										if( fabs(z-zpos[q]) < 1.0e-7 )
											found = true;
									}
									if( !found )
									{
										zpos[nzpos] = z;
										//zalpha[nzpos] = (z-coords_array[i][j][1][2])/(coords_array[i][j][0][2]-coords_array[i][j][1][2]);
										zalpha[nzpos] = (z-xyz[ECL_IJK_COORDS(i,j,1,2)])/(xyz[ECL_IJK_COORDS(i,j,0,2)]-xyz[ECL_IJK_COORDS(i,j,1,2)]);
										nzpos++;
									}
								}
							}
							for(int l = 0; l < nzpos; ++l)
							{
								pillar::iterator search = p.find(zpos[l]);
								Node n; //node added to pillar
								if( search == p.end() )
								{
									Storage::real node_xyz[3];
									Storage::real alpha = zalpha[l];
									/*
									node_xyz[0] = coords_array[i][j][1][0] + alpha * (coords_array[i][j][0][0] - coords_array[i][j][1][0]);
									node_xyz[1] = coords_array[i][j][1][1] + alpha * (coords_array[i][j][0][1] - coords_array[i][j][1][1]);
									 */
									node_xyz[0] = xyz[ECL_IJK_COORDS(i,j,1,0)] + alpha * (xyz[ECL_IJK_COORDS(i,j,0,0)] - xyz[ECL_IJK_COORDS(i,j,1,0)]);
									node_xyz[1] = xyz[ECL_IJK_COORDS(i,j,1,1)] + alpha * (xyz[ECL_IJK_COORDS(i,j,0,1)] - xyz[ECL_IJK_COORDS(i,j,1,1)]);
									node_xyz[2] = zpos[l];//zcorn_array[bi][bj][bk][7-l];
									int find = -1;
									if( !old_nodes.empty() )
									{
										std::vector<HandleType>::iterator it = std::lower_bound(old_nodes.begin(),old_nodes.end(),node_xyz,CentroidComparator(this));
										if( it != old_nodes.end() )
										{
											Storage::real_array c = RealArrayDF(*it,CoordsTag());
											if( CentroidComparator(this).Compare(node_xyz,c.data()) == 0 )
												find = static_cast<int>(it - old_nodes.begin());
										}
									}
									n = find == -1 ? CreateNode(node_xyz) : Node(this,old_nodes[find]);
									p.insert(std::make_pair(zpos[l],n));
								} else n = search->second; //if
							}
						}
					}
				}
				/*      (6)*<<------[3]--------*(7)
						  /|                  /|
						[6]                 [7]|
						/  |                /  |
					(4)*--------[2]------>>*(5)|
					   |  [10]             |  [11]
					   |                   |   |
					   |   |               |   |
					   |                   |   |
					  [8]  |              [9]  |
					   |(2)*- - - - [1]- - |->>*(3)
					   |  /                |  /
					   |[4]                |[5]
					   |/                  |/
					(0)*<<------[0]--------*(1)
				*/
				//all edges along pillars
				std::vector< ElementArray<Edge> > pillar_edges((dims[0]+1)*(dims[1]+1),ElementArray<Edge>(this));
				//all edges of each block
				//std::vector< std::vector<Edge> > block_edges(dims[0]*dims[1]*dims[2]);
				//create edges along pillars
				printf("create edges along pillars\n");
#if defined(USE_OMP)
#pragma omp parallel
#endif
				{
					//structure to create an edge from pair of nodes
					ElementArray<Node> edge_nodes(this,2);
#if defined(USE_OMP)
#pragma omp for
#endif
					for(int i = 0; i < dims[0]+1; ++i)
					{
						for(int j = 0; j < dims[1]+1; ++j)
						{
							pillar & p = pillars[i*(dims[1]+1)+j];
							if( p.size() > 1 )
							{
								ElementArray<Edge> & p_edges = pillar_edges[i*(dims[1]+1)+j];
								pillar::iterator it = p.begin();
								pillar::iterator pre_it = it++;
								while(it != p.end())
								{
									edge_nodes[0] = it->second;
									edge_nodes[1] = pre_it->second;
									p_edges.push_back(CreateEdge(edge_nodes).first);
									pre_it = it++;
								}
							}
						} //j
					} //i
				}
				//mark edges along pillars (vertical edges, 8, 9, 10, 11)
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(int i = 0; i < dims[0]; ++i)
				{
					for(int j = 0; j < dims[1]; ++j)
					{
						for(int k = 0; k < dims[2]; ++k) if( actnum.empty() || actnum[ECL_IJK_DATA(i,j,k)] )
						{
							const int edge_nodes[4][2] = 
							{
								{0,4}, //edge 8
								{1,5}, //edge 9
								{2,6}, //edge 10
								{3,7}  //edge 11
							}; //nodes of edges 8,9,10,11, bottom->top
							for(int l = 0; l < 4; ++l) //loop edges
							{
								HandleType find[2];
								int i2 = i + l%2, j2 = j + l/2;
								pillar & p = pillars[i2*(dims[1]+1) + j2];
								pillar::iterator b,t, it, jt; 
								//find nodes in pillar that belong to current block
								Storage::real pa = zcorn[ECL_IJK_ZCORN(i,j,k,edge_nodes[l][0])];
								Storage::real pb = zcorn[ECL_IJK_ZCORN(i,j,k,edge_nodes[l][1])];
								if(pa > pb) //otherwise run iterator in different direction
								{
									b = p.find(pb); //bottom
									t = p.find(pa); //top
								}
								else
								{
									b = p.find(pa); //bottom
									t = p.find(pb); //top
								}
								//go over edges
								it = b;
								while(it != t)
								{
									jt = it++; // jt<->it forms a segment
									find[0] = jt->second->GetHandle();
									find[1] = it->second->GetHandle();
									Edge e(this,FindSharedAdjacency(find,2)); //find edge
#if defined(USE_OMP)
#pragma omp critical
#endif
									{
										e->IntegerArray(block_number).push_back((i + (j+k*dims[1])*dims[0]));
										e->IntegerArray(edge_number).push_back(8+l);
									}
									//block_edges[(i + (j+k*dims[1])*dims[0])].push_back(e);
								}
							}
						}
					}
				}
				printf("erase unused edges on pillars\n");
				//erase edges on pillars that are not used in any block
#if defined(USE_OMP)
#pragma omp parallel for
#endif

				for(int i = 0; i < dims[0]+1; ++i)
				{
					for(int j = 0; j < dims[1]+1; ++j)
					{
						ElementArray<Edge> & p = pillar_edges[i*(dims[1]+1) + j];
						ElementArray<Edge>::iterator it = p.begin();
						while(it != p.end())
						{
							if( it->IntegerArray(block_number).empty() )
							{
								it->Delete();
								it = p.erase(it);
							}
							else ++it;
						}
					}
				}
				
				//Tag pillar_mark = CreateTag("PILLAR_MARK",DATA_INTEGER,EDGE,NONE,3);
				//tags to be transfered
				std::vector<Tag> transfer(2);
				transfer[0] = edge_number;
				transfer[1] = block_number;
				//transfer[2] = pillar_mark;
				//store faces for each block to assemble block cells later
				std::vector< ElementArray<Face> > block_faces(dims[0]*dims[1]*dims[2],ElementArray<Face>(this));
				//store top-bottom edges to assemble top and bottom faces later
				std::vector< ElementArray<Edge> > block_edges(dims[0]*dims[1]*dims[2],ElementArray<Edge>(this));
				//print out intersection algorithm
				bool print_inter = false;
				//print out block number information
				bool print_bn = false;
				//print out block edges information
				bool print_bedges = false;
				//some intermediate info
				bool print_info = false;
				
				//check_shared_mrk = false;
				//check_private_mrk = false;

				//go over nx pairs of pillars, then ny pairs of pillars
				for(int q = 0; q < 2; ++q)
				{
					printf("started creating faces for pairs of pillar along %s\n",q ? "ny":"nx");
#if defined(USE_OMP)
#pragma omp parallel
#endif
					{
						//structure to create an edge from pair of nodes
						ElementArray<Node> edge_nodes(this,2);
						//set of lines to be intersected
						ElementArray<Edge> edges(this);
						ElementArray<Node> intersections(this);
						//for sorting a pair of data
						std::vector<int> indices_sort, temporary;
						//array to obtain block number
						std::vector<int> block_number_inter;
						//store edges along pillar to assemble faces along pillars, on back side of the pillar
						std::vector< ElementArray<Edge> > pillar_block_edges_back(dims[2],ElementArray<Edge>(this));
						//on front side of the pillar
						std::vector< ElementArray<Edge> > pillar_block_edges_front(dims[2],ElementArray<Edge>(this));
						//mark edges on pillar for adjacency retrival
						MarkerType mrk = CreatePrivateMarker();
						//mark original edges of each block face, so that we know outer boundary on constructed interface
						MarkerType outer = CreatePrivateMarker();
						//remember visited edges
						MarkerType visited = CreatePrivateMarker();
						//create tag that will be used for intersection of segments in projected space of [0,1]^2 quad
						std::stringstream tag_name;
						tag_name << "PROJ_PNT_PROC_" << GetLocalProcessorRank();
						Tag pnt = CreateTag(tag_name.str(),DATA_REAL,NODE,NONE,2);
#if defined(USE_OMP)
#pragma omp for
#endif
						for(int i = 0; i < dims[0]+q; i++)
						{
							//printf("%s %6.2f%%\r",q ? "ny":"nx", ((Storage::real)i)/((Storage::real)dims[0]+q-1)*100);
							//fflush(stdout);
							for(int j = 0; j < dims[1]+!q; ++j)
							{
								if( print_info )
								{
									std::cout << "working on " << (q ? "ny" : "nx") << " pair of pillars: " << i << "," << j << " and " << i+!q << "," << j+q << std::endl;
								}
								//p0 is at i,j
								//p1 is at i+1,j, when q = 0 and at i,j+1, when q = 1
								pillar & p0 = pillars[(i+ 0)*(dims[1]+1)+j+0];
								pillar & p1 = pillars[(i+!q)*(dims[1]+1)+j+q];
								//add projected point information
								for(pillar::iterator it = p0.begin(); it != p0.end(); ++it)
								{
									Storage::real_array pos = it->second->RealArray(pnt);
									pos[0] = (it->first-xyz[ECL_IJK_COORDS(i,j,1,2)])/(xyz[ECL_IJK_COORDS(i,j,0,2)]-xyz[ECL_IJK_COORDS(i,j,1,2)]);
									pos[1] = 0;
								}
								for(pillar::iterator it = p1.begin(); it != p1.end(); ++it)
								{
									Storage::real_array pos = it->second->RealArray(pnt);
									pos[0] = (it->first-xyz[ECL_IJK_COORDS(i+!q,j+q,1,2)])/(xyz[ECL_IJK_COORDS(i+!q,j+q,0,2)]-xyz[ECL_IJK_COORDS(i+!q,j+q,1,2)]);
									pos[1] = 1;
								}
								//preallocate array for edges
								edges.reserve(2*std::max(p0.size(),p1.size()));
								//add far faces of j-1 block
								if( (1-q)*j + q*i > 0 ) //test j when q = 0 and i when q = 1
								{
									if( print_info ) std::cout << "back cell: " << i-q << "," << j-!q << std::endl;
									for(int k = 0; k < dims[2]; ++k) if( actnum.empty() || actnum[ECL_IJK_DATA(i-q,j-!q,k)] )
									{
										for(int l = 0; l < 2; ++l) //top-bottom
										{
											//q = 0 - test 2,3 + 4*l (edge 1 and 3)
											//q = 1 - test 1,3 + 4*l (edge 5 and 7)
											edge_nodes[0] = p0[zcorn[ECL_IJK_ZCORN(i-q,j-!q,k,2 - q + 4*l)]];
											edge_nodes[1] = p1[zcorn[ECL_IJK_ZCORN(i-q,j-!q,k,3 - 0 + 4*l)]];
											if( edge_nodes[0] != edge_nodes[1] )
											{
												//edge_nodes.SetMarker(used);
												Edge e = CreateEdge(edge_nodes).first;
												if( !e->GetPrivateMarker(visited) )
												{
													edges.push_back(e);
													e->SetPrivateMarker(visited);
												}
												e->IntegerArray(block_number).push_back((i-q + (j-!q+k*dims[1])*dims[0]));
												e->IntegerArray(edge_number).push_back(1 + 2*l + 4*q);
												//Storage::integer_array mark = e->IntegerArray(pillar_mark);
												//back block indices
												//mark[0] = i;
												//mark[1] = j;
												//nx or ny
												//mark[2] = q;
											}
										}
									}
								}
								//add near faces of j block
								if( (1-q)*j + q*i  < dims[!q] ) //test j when q = 0 and i when q = 1
								{
									if( print_info ) std::cout << "front cell: " << i << "," << j << std::endl;
									for(int k = 0; k < dims[2]; ++k)  if( actnum.empty() || actnum[ECL_IJK_DATA(i,j,k)] )
									{
										for(int l = 0; l < 2; ++l) //top-bottom
										{
											//q = 0 - test 0,1 + 4*l (edge 0 and 2)
											//q = 1 - test 0,2 + 4*l (edge 4 and 6)
											edge_nodes[0] = p0[zcorn[ECL_IJK_ZCORN(i,j,k,0 + 0 + 4*l)]];
											edge_nodes[1] = p1[zcorn[ECL_IJK_ZCORN(i,j,k,1 + q + 4*l)]];
											if( edge_nodes[0] != edge_nodes[1] )
											{
												//edge_nodes.SetMarker(used);
												Edge e = CreateEdge(edge_nodes).first;
												if( !e->GetPrivateMarker(visited) )
												{
													edges.push_back(e);
													e->SetPrivateMarker(visited);
												}
												e->IntegerArray(block_number).push_back((i + (j+k*dims[1])*dims[0]));
												e->IntegerArray(edge_number).push_back(0 + 2*l + 4*q);
												//Storage::integer_array mark = e->IntegerArray(pillar_mark);
												//front block indices
												//mark[0] = i;
												//mark[1] = j;
												//nx or ny
												//mark[2] = q;
											}
										}
									}
								}
								edges.RemPrivateMarker(visited);
								//produce intersected edges
								if( print_inter )
								{
									std::cout << "input edges: " << edges.size() << std::endl;
									if( false ) for(int k = 0; k < edges.size(); ++k)
									{
										Storage::integer_array bn = edges[k]->IntegerArray(block_number);
										Storage::integer_array en = edges[k]->IntegerArray(edge_number);
										std::cout << "edge " << k << " " << edges[k]->GetHandle() << " " << edges[k]->getBeg()->GetHandle() << "<->" << edges[k]->getEnd()->GetHandle() << " blocks: ";
										for(int l = 0; l < (int)bn.size(); ++l)
											std::cout << bn[l] << "(" << bn[l]%dims[0] << "," << bn[l]/dims[0]%dims[1] << "," << bn[l]/dims[0]/dims[1] << "):" << en[l] << " ";
										std::cout << std::endl;
									}

									for(int k = 0; k < edges.size(); ++k)
									{
										std::cout << "(" << edges[k]->getBeg()->Coords()[0] << "," << edges[k]->getBeg()->Coords()[1] << "," << edges[k]->getBeg()->Coords()[2] << ") <-> (" << edges[k]->getEnd()->Coords()[0] << "," << edges[k]->getEnd()->Coords()[1] << "," << edges[k]->getEnd()->Coords()[2] << ")" << std::endl;
									}
								}
								assert(count_duplicates(edges) == 0);
								//i << "," << j << " and " << i+!q << "," << j+q <<
								Unproject unp(&xyz[ECL_IJK_COORDS(i,j,0,0)],&xyz[ECL_IJK_COORDS(i,j,1,0)],&xyz[ECL_IJK_COORDS(i+!q,j+q,0,0)],&xyz[ECL_IJK_COORDS(i+!q,j+q,1,0)]);
								//intersect(this,edges,intersections,transfer,pnt,unp,print_inter);
								intersect_naive(this,edges,intersections,transfer,pnt,unp,print_inter);
								assert(count_duplicates(edges) == 0);
								if( print_inter ) std::cout << "intersections: " << intersections.size() << std::endl;
								if( print_inter )
								{

									std::cout << "output edges: " << edges.size() << std::endl;
									for(int k = 0; k < edges.size(); ++k)
									{
										//std::cout << "edge " << k << " " << edges[k]->getBeg()->GetHandle() << "<->" << edges[k]->getEnd()->GetHandle() << std::endl;
										std::cout << "(" << edges[k]->getBeg()->Coords()[0] << "," << edges[k]->getBeg()->Coords()[1] << "," << edges[k]->getBeg()->Coords()[2] << ") <-> (" << edges[k]->getEnd()->Coords()[0] << "," << edges[k]->getEnd()->Coords()[1] << "," << edges[k]->getEnd()->Coords()[2] << ")" << std::endl;

									}
								}

								//distribute all the edges among blocks
								//add intersected edges into blocks, so that we can reconstruct top and bottom faces of each block
								for(int k = 0; k < (int)edges.size(); ++k)
								{
									Storage::integer_array b = edges[k].IntegerArray(block_number);
									for(int r = 0; r < (int)b.size(); ++r)
										block_edges[b[r]].push_back(edges[k]);
								}
								//add vertical edges along pillars
								ElementArray<Edge> & p0_edges = pillar_edges[(i+ 0)*(dims[1]+1)+j+0];
								ElementArray<Edge> & p1_edges = pillar_edges[(i+!q)*(dims[1]+1)+j+q];
								edges.reserve(edges.size()+p0_edges.size()+p1_edges.size());
								for(int k = 0; k < (int)p0_edges.size(); ++k)
								{
									//if( p0_edges[k]->nbAdjElements(NODE,used) == 2 ) //this edge is used
									edges.push_back(p0_edges[k]);
								}
								for(int k = 0; k < (int)p1_edges.size(); ++k)
								{
									//if( p1_edges[k]->nbAdjElements(NODE,used) == 2 ) //this edge is used
									edges.push_back(p1_edges[k]);
								}
								assert(count_duplicates(edges) == 0);
								//sort block numbers on edges
								for(int k = 0; k < (int)edges.size(); ++k)
								{
									Storage::integer_array b = edges[k].IntegerArray(block_number);
									Storage::integer_array e = edges[k].IntegerArray(edge_number);
									assert(e.size() == b.size());
									//sort indices according to b
									indices_sort.resize(b.size());
									for(int l = 0; l < indices_sort.size(); ++l) indices_sort[l] = l;
									std::sort(indices_sort.begin(),indices_sort.end(),index_comparator(b));
									//arrange data in b and e arrays according to indices_sort
									temporary.resize(b.size());
									//first b array
									for(int l = 0; l < (int)b.size(); ++l) temporary[l] = b[l];
									for(int l = 0; l < (int)b.size(); ++l) b[l] = temporary[indices_sort[l]];
									//then e array
									for(int l = 0; l < (int)e.size(); ++l) temporary[l] = e[l];
									for(int l = 0; l < (int)e.size(); ++l) e[l] = temporary[indices_sort[l]];
								}
								//put block numbers to all nodes involved in pillars, so that we can figure out block numbers for constructed faces
								edges.SetPrivateMarker(mrk);
								//report block numbers on edges
								if(print_bn)
								{
									for(int k = 0; k < (int)edges.size(); ++k)
									{
										Storage::integer_array b = edges[k]->IntegerArray(block_number);
										Storage::integer_array e = edges[k]->IntegerArray(edge_number);
										std::cout << "edge " << k << " " << edges[k]->GetHandle() << " blocks [" << b.size() << "]: ";
										for(int l = 0; l < (int)b.size()-1; ++l) std::cout << b[l] << "(" << b[l]%dims[0] << "," << b[l]/dims[0]%dims[1] << "," << b[l]/dims[0]/dims[1] << "):" << e[l] << " ";
										std::cout << b.back() << "(" << b.back()%dims[0] << "," << b.back()/dims[0]%dims[1] << "," << b.back()/dims[0]/dims[1] << "):" << e.back() << std::endl;
									}
								}
								//unite block numbers on nodes
								if( print_bn ) std::cout << "pillar 0 nodes: " << std::endl;
								for(pillar::iterator it = p0.begin(); it != p0.end(); ++it) //if( it->second.GetMarker(used) )
								{
									ElementArray<Element> nedges = it->second->getAdjElements(EDGE,mrk);
									block_number_union(it->second.getAsElement(),nedges,block_number);
									if( print_bn )
									{
										Storage::integer_array b = it->second->IntegerArray(block_number);
										std::cout << "node " << it->second->GetHandle() << " blocks [" << b.size() << "]: ";
										for(int l = 0; l < (int)b.size()-1; ++l) std::cout << b[l] << "(" << b[l]%dims[0] << "," << b[l]/dims[0]%dims[1] << "," << b[l]/dims[0]/dims[1] << "), ";
										std::cout << b.back() << "(" << b.back()%dims[0] << "," << b.back()/dims[0]%dims[1] << "," << b.back()/dims[0]/dims[1] << ")" << std::endl;
									}
								}
								if( print_bn ) std::cout << "pillar 1 nodes: " << std::endl;
								for(pillar::iterator it = p1.begin(); it != p1.end(); ++it) //if( it->second.GetMarker(used) )
								{
									ElementArray<Element> nedges = it->second->getAdjElements(EDGE,mrk);
									block_number_union(it->second.getAsElement(),nedges,block_number);
									if( print_bn )
									{
										Storage::integer_array b = it->second->IntegerArray(block_number);
										std::cout << "node " << it->second->GetHandle() << " blocks [" << b.size() << "]: ";
										for(int l = 0; l < (int)b.size()-1; ++l) std::cout << b[l] << "(" << b[l]%dims[0] << "," << b[l]/dims[0]%dims[1] << "," << b[l]/dims[0]/dims[1] << "), ";
										std::cout << b.back() << "(" << b.back()%dims[0] << "," << b.back()/dims[0]%dims[1] << "," << b.back()/dims[0]/dims[1] << ")" << std::endl;
									}
								}
								if( print_bn ) std::cout << "intersection nodes: " << std::endl;
								for(int k = 0; k < (int)intersections.size(); ++k)
								{
									ElementArray<Element> nedges = intersections[k]->getAdjElements(EDGE,mrk);
									block_number_union(intersections[k].getAsElement(),nedges,block_number);
									if( print_bn )
									{
										Storage::integer_array b = intersections[k]->IntegerArray(block_number);
										std::cout << "node " << k << " " << intersections[k]->GetHandle() << " blocks [" << b.size() << "]: ";
										for(int l = 0; l < (int)b.size()-1; ++l) std::cout << b[l] << "(" << b[l]%dims[0] << "," << b[l]/dims[0]%dims[1] << "," << b[l]/dims[0]/dims[1] << "), ";
										std::cout << b.back() << "(" << b.back()%dims[0] << "," << b.back()/dims[0]%dims[1] << "," << b.back()/dims[0]/dims[1] << ")" << std::endl;
									}
								}
								//unmark nodes
								//for(pillar::iterator it = p0.begin(); it != p0.end(); ++it)  it->second.RemMarker(used);
								//for(pillar::iterator it = p1.begin(); it != p1.end(); ++it)  it->second.RemMarker(used);
								//unmark edges
								edges.RemPrivateMarker(mrk);
								//distribute edges to front and back blocks, so that we can assemble faces
								for(int k = 0; k < (int)edges.size(); ++k)
								{
									//intersect block numbers on nodes to detect to which blocks this edge belongs
									ElementArray<Element> nodes = edges[k]->getAdjElements(NODE);
									block_number_intersection(nodes,block_number,block_number_inter);
									for(int m = 0; m < (int)block_number_inter.size(); ++m)
									{
										int bi = block_number_inter[m] % dims[0];
										int bj = block_number_inter[m] / dims[0] % dims[1];
										int bk = block_number_inter[m] / dims[0] / dims[1];
										assert(bi >= 0 && bi < dims[0]);
										assert(bj >= 0 && bj < dims[1]);
										if( bi == i-q && bj == j-!q ) //back blocks
										{
											pillar_block_edges_back[bk].push_back(edges[k]);
											if( print_bn ) std::cout << "add edge " << k << " " << edges[k]->GetHandle() << " to back  block " << block_number_inter[m] << " (" << bi << "," << bj << "," << bk << ")" << std::endl;
										}
										else if(bi == i && bj == j) //front blocks
										{
											pillar_block_edges_front[bk].push_back(edges[k]);
											if( print_bn ) std::cout << "add edge " << k << " " << edges[k]->GetHandle() << " to front block " << block_number_inter[m] << " (" << bi << "," << bj << "," << bk << ")" << std::endl;
										}
										//skip non-current blocks
									}
								}
								std::vector< ElementArray<Edge> > * pillar_block_edges[2] = {&pillar_block_edges_back,&pillar_block_edges_front};
								//(i-q + (j-!q+k*dims[1])*dims[0])
								int blocki[2] = {i-q,i};
								int blockj[2] = {j-!q,j};
								if( print_info )
								{
									std::cout << "back  block " << blocki[0] << "," << blockj[0] << std::endl;
									std::cout << "front block " << blocki[1] << "," << blockj[1] << std::endl;
									std::cout << "dims " << dims[0] << "," << dims[1] << std::endl;
								}
								//construct interfaces (some computations will be redundant)
								for(int m = 0; m < 2; ++m) //back and front
								{
									if( blocki[m] == -1 || blockj[m] == -1 ) {if(print_info) std::cout << "skip " << (m?"front ":"back ") << std::endl; continue;}
									if( blocki[m] == dims[0] || blockj[m] == dims[1] ) {if(print_info) std::cout << "skip " << (m?"front ":"back ") << std::endl; continue;}
									if( print_bedges ) std::cout << (m?"front ":"back ") << " column of blocks: " << blocki[m] << "," << blockj[m] << std::endl;
									for(int k = 0; k < dims[2]; ++k)  if( actnum.empty() || actnum[ECL_IJK_DATA(blocki[m],blockj[m],k)] ) //go down the piller
									{
										//retrive edge for the side
										ElementArray<Edge> & bedges = pillar_block_edges[m]->at(k);
										if( bedges.empty() ) continue;
										//remove any duplicates
										make_unique(bedges);
										int num_outer = 0;
										std::set<int> outer_edge_number;
										for(int l = 0; l < bedges.size(); ++l) //loop through edges of the block
										{
											//retrive block numbers of edges
											Storage::integer_array bn = bedges[l]->IntegerArray(block_number);
											Storage::integer_array en = bedges[l]->IntegerArray(edge_number);
											for(int r = 0; r < (int)bn.size(); ++r)
											{
												if( bn[r] == blocki[m] + (blockj[m]+k*dims[1])*dims[0] ) //this edge originally created by the block
												{
													bedges[l]->SetPrivateMarker(outer); //mark edge
													num_outer++;
													outer_edge_number.insert(en[r]);
												}
											}
										}
										if( outer_edge_number.size() > 2 )
										{
											//move marked edges to the end
											std::sort(bedges.begin(),bedges.end(),Mesh::PrivateMarkerComparator(this,outer));
											if( print_bedges )
											{
												std::cout << (m?"front ":"back ") << "depth " << k << " block " <<blocki[m] + (blockj[m]+k*dims[1])*dims[0] << " edges [" << bedges.size() << "]:" << std::endl;
												for(int l = 0; l < bedges.size(); ++l)
												{
													Storage::integer_array bn = bedges[l]->IntegerArray(block_number);
													Storage::integer_array en = bedges[l]->IntegerArray(edge_number);
													std::cout << "edge " << l << " " << bedges[l]->GetHandle() << " " << (bedges[l]->GetPrivateMarker(outer) ? "outer":"inner");
													std::cout << " blocks ";
													for(int r = 0; r < (int)bn.size(); ++r)
														std::cout << bn[r] << ":" << en[r] << " ";
													std::cout << std::endl;

												}
												assert(count_duplicates(bedges) == 0);
											}
											//there is enough edges to handle triangle
											if( bedges.size() > 2 )
											{
												//form faces out of edges
												incident_matrix<Edge> matrix(this,bedges.data(),bedges.data()+bedges.size(),bedges.size()-num_outer);
												//collect all faces
												ElementArray<Edge> loop(this);
												while(matrix.find_shortest_loop(loop))
												{
													if( loop.size() > 2 ) //at least triangle
													{
														if( print_bedges )
														{
															std::cout << "Found loop of " << loop.size() << " edges:" <<std::endl;
															for(int g = 0; g < loop.size(); ++g)
															{
																std::cout << "edge " << g << " " << loop[g]->GetHandle() << " ";
																Storage::integer_array bn = loop[g]->IntegerArray(block_number);
																Storage::integer_array en = loop[g]->IntegerArray(edge_number);
																std::cout << " blocks ";
																for(int r = 0; r < (int)bn.size(); ++r)
																	std::cout << bn[r] << ":" << en[r] << " ";
																std::cout << std::endl;
															}

														}
														//make face
														Face f = CreateFace(loop).first;
														if(!f->CheckEdgeOrder()) std::cout << __FILE__ << ":" << __LINE__ << " bad edge order, edges " << loop.size() << std::endl;
														//detect block number
														ElementArray<Element> nodes = f->getAdjElements(NODE);
														if( print_bedges )
														{
															std::cout << "nodes of the face " << f->GetHandle() << "(" << f->LocalID() << "): " << std::endl;
															for(int g = 0; g < nodes.size(); ++g)
															{
																std::cout << "node " << g << " " << nodes[g]->GetHandle() << " ";
																Storage::integer_array bn = nodes[g]->IntegerArray(block_number);
																std::cout << " blocks ";
																for(int r = 0; r < (int)bn.size(); ++r)
																	std::cout << bn[r] << " ";
																std::cout << std::endl;

															}
															std::cout << "intersection: ";
														}
														block_number_intersection(nodes,block_number,block_number_inter);
														for(int r = 0; r < (int)block_number_inter.size(); ++r)
														{
															if( print_bedges ) std::cout << block_number_inter[r] << " (" << block_number_inter[r]%dims[0] << "," << block_number_inter[r]/dims[0] % dims[1] << "," << block_number_inter[r]/dims[0]/dims[1] <<") ";
															block_faces[block_number_inter[r]].push_back(f);
														}
														if( print_bedges ) std::cout << std::endl;
													}
													else if( !loop.empty() ) std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate face " << loop.size() << std::endl;
												}
												if( !matrix.all_visited() )
												{
													std::cout << "Not all edges were visited, matrix: " << std::endl;
													matrix.print_matrix();
													std::cout << "Current block " << blocki[m] + (blockj[m]+k*dims[1])*dims[0] << std::endl;
													for(int l = 0; l < bedges.size(); ++l)
													{
														Storage::integer_array bn = bedges[l]->IntegerArray(block_number);
														Storage::integer_array en = bedges[l]->IntegerArray(edge_number);
														std::cout << "edge " << l << " " << bedges[l]->GetHandle() << " " << bedges[l]->getBeg()->GetHandle() << " <-> " << bedges[l]->getEnd()->GetHandle() << " " << (bedges[l]->GetPrivateMarker(outer) ? "outer":"inner");
														std::cout << " blocks ";
														for(int r = 0; r < (int)bn.size(); ++r)
															std::cout << bn[r] << ":" << en[r] << " ";
														std::cout << std::endl;
													}

													std::cout << "Outer edges [" << outer_edge_number.size() << "]:" << std::endl;
													for(std::set<int>::iterator it = outer_edge_number.begin(); it != outer_edge_number.end(); ++it)
														std::cout << *it << " ";
													std::cout << std::endl;

													//form faces out of edgesm
													incident_matrix<Edge> matrix(this,bedges.data(),bedges.data()+bedges.size(),bedges.size()-num_outer,true);
													//collect all faces
													ElementArray<Edge> loop2(this);
													while(matrix.find_shortest_loop(loop2));
													
													std::cout << "All edges: " << std::endl;
													for(int k = 0; k < (int)edges.size(); ++k)
													{
														std::cout << "(" << edges[k]->getBeg()->Coords()[0] << "," << edges[k]->getBeg()->Coords()[1] << "," << edges[k]->getBeg()->Coords()[2] << ") <-> (" << edges[k]->getEnd()->Coords()[0] << "," << edges[k]->getEnd()->Coords()[1] << "," << edges[k]->getEnd()->Coords()[2] << ")" << std::endl;
													}
												}
											}
											//remove marker
											bedges.RemPrivateMarker(outer);
										}
										//cleanup structure for reuse
										bedges.clear();
									}
								}
								//clean-up structures
								edges.clear();
								intersections.clear();
							} //j
						} //i
						DeleteTag(pnt);
						ReleasePrivateMarker(visited);
						ReleasePrivateMarker(mrk);
						ReleasePrivateMarker(outer);
					} //omp parallel
					//printf("\n");
					printf("finished creating faces for pairs of pillar along %s\n",q ? "ny":"nx");
				} //q
				//do not need this tag on nodes
				DeleteTag(block_number, NODE);
				//now construct top and bottom interfaces
				printf("started tops/bottoms/cells\n");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(int i = 0; i < dims[0]; ++i)
				{
					//printf("top/bottom/cells %6.2f%%\r", ((Storage::real)i)/((Storage::real)dims[0]-1)*100);
					//fflush(stdout);
					for(int j = 0; j < dims[1]; ++j)
					{
						for(int k = 0; k < dims[2]; ++k) if( actnum.empty() || actnum[ECL_IJK_DATA(i,j,k)] )
						{
							//current block number
							int cur = (i + (j+k*dims[1])*dims[0]);
							//bottom face -> 0,4,1,5
							//top face -> 2,7,3,6
							const int map[8][2] =
							{
								{0,0}, //edge 0 -> first position, bottom face
								{2,0}, //edge 1 -> third position, bottom face
								{0,1}, //edge 2 -> first position, top face
								{2,1}, //edge 3 -> third position, top face
								{1,0}, //edge 4 -> second position, bottom face
								{3,0}, //edge 5 -> fourth position, bottom face
								{3,1}, //edge 6 -> fourth position, top face
								{1,1} //edge 7 -> second position, top face
							};
							//which edges should be considered in reverse order
							const bool rev[2][4] = 
							{
								{true,false,false,true},
								{false,false,true,true}
							};
							//array of ordered edges, first dimension - top or bottom, second dimension - edge number
							std::vector<HandleType> edges[2][4];
							//retrive set of edges of the block
							ElementArray<Edge> & cbe = block_edges[cur];
							for(int q = 0; q < cbe.size(); ++q)
							{
								Storage::integer_array bn = cbe[q].IntegerArray(block_number);
								Storage::integer_array en = cbe[q].IntegerArray(edge_number);
								assert(bn.size() == en.size());
								for(int r = 0; r < (int)bn.size(); ++r) if( bn[r] == cur ) //how the edge is represented on current block
									edges[map[en[r]][1]][map[en[r]][0]].push_back(cbe[q].GetHandle());
							}
							ElementArray<Edge> face_edges(this);
							//collect edges and create face
							for(int q = 0; q < 2; ++q) //top and bottom
							{
								for(int r = 0; r < 4; ++r) //all sides
								{
									if( rev[q][r] )
										face_edges.insert(face_edges.end(),edges[q][r].rbegin(),edges[q][r].rend());
									else
										face_edges.insert(face_edges.end(),edges[q][r].begin(),edges[q][r].end());
								}
								make_unique(face_edges); //somehow duplicate appears there
								if( face_edges.size() > 2 )
								{
									Face f = CreateFace(face_edges).first;
									
									if(!f->FixEdgeOrder())
									//if(!f->CheckEdgeOrder())
									{
										std::cout << __FILE__ << ":" << __LINE__ << " bad edge order, edges " << face_edges.size() << std::endl;
										std::cout << "block: " << cur << " (" << i << "," << j << "," << k << ") " << (q ? "top" : "bottom") << std::endl;
										for(int l = 0; l < face_edges.size(); ++l)
										{
											Storage::integer_array bn = face_edges[l]->IntegerArray(block_number);
											Storage::integer_array en = face_edges[l]->IntegerArray(edge_number);
											std::cout << "edge " << l << " " << face_edges[l]->GetHandle() << " ";
											std::cout << face_edges[l]->getBeg()->GetHandle() << "<->" << face_edges[l]->getEnd()->GetHandle() << " ";
											std::cout << "blocks: ";
											for(int p = 0; p < (int)bn.size(); ++p)
												std::cout << bn[p] << "(" << bn[p]%dims[0] << "," << bn[p]/dims[0]%dims[1] << "," << bn[p]/dims[0]/dims[1] << "):" << en[p] << " ";
											std::cout << std::endl;
										}
									}
									block_faces[cur].push_back(f);
								}
								else if( !face_edges.empty() ) std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate face " << face_edges.size() << " " << (q?"top":"bottom") << " of block " << cur << " (" << i << "," << j << "," << k << ") actnum " << actnum[ECL_IJK_DATA(i,j,k)] << std::endl;
								face_edges.clear();
							}
							make_unique(block_faces[cur]); //some faces may be added twice?
							if( block_faces[cur].size() > 3 )
							{
								Cell c =CreateCell(block_faces[cur]).first;
								c->Integer(cell_number) = cur+1;
								//c->IntegerArray(cell_number)[1] = i;
								//c->IntegerArray(cell_number)[2] = j;
								//c->IntegerArray(cell_number)[3] = k;
							}
							else if( !block_faces[cur].empty() )
							{
								//depending on actnum, mark faces that are not connected
								//std::cout << __FILE__ << ":" << __LINE__ << " skip degenerate cell " << block_faces[cur].size() << " block " << cur << " (" << i << "," << j << "," << k << ") actnum " << actnum[ECL_IJK_DATA(i,j,k)] << std::endl;
							}
						} //k
					} //j
				} //i
				printf("finished tops/bottoms/cells\n");
				//printf("\n");
				//cleanup data
				DeleteTag(edge_number);
				DeleteTag(block_number);
				//crack up the mesh along degenerate active cells
				//populate properties to blocks
				Tag tagporo, tagsatnum, tagperm;
				if( !poro.empty() ) tagporo = CreateTag("PORO",DATA_REAL,CELL,NONE,1);
				if( !perm.empty() ) tagperm = CreateTag("PERM",DATA_REAL,CELL,NONE,3);
				if( !satnum.empty() ) tagsatnum = CreateTag("SATNUM",DATA_INTEGER,CELL,NONE,1);

#if defined(USE_OMP)
#pragma omp parallel for
#endif
				for(integer it = 0; it < CellLastLocalID(); ++it) if( isValidCell(it) )
				{
					Cell c = CellByLocalID(it);
					if( c->Integer(cell_number) > 0 ) //maybe this cell existed before
					{
						integer q = c->Integer(cell_number)-1;
						if( !poro.empty() ) c->Real(tagporo) = poro[q];
						if( !satnum.empty() ) c->Integer(tagsatnum) = satnum[q];
						if( !perm.empty() )
						{
							Storage::real_array K = c->RealArray(tagperm);
							K[0] = perm[q*3+0];
							K[1] = perm[q*3+1];
							K[2] = perm[q*3+2];
						}
					}
				}
			}
		}
	}
}

#endif