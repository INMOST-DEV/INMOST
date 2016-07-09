#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"

#if defined(USE_MESH)

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

#define ECL_GTYPE_NONE 0
#define ECL_GTYPE_TOPS 1
#define ECL_GTYPE_ZCORN 2
#define ECL_GTYPE_RADIAL 3
#define ECL_GTYPE_CARTESIAN 4


#define ECL_VAR_NONE 0
#define ECL_VAR_REAL 1
#define ECL_VAR_INT 2

#define ECL_IJK_DATA(i,j,k) (i + ((j)+(k)*dims[1])*dims[0])

//line sweep events
#define SEG_START 0
#define SEG_END 1



namespace INMOST
{
	//special structure for array of 3 reals for std::map
	class position
	{
	public:
		Storage::real xyz[3];
		position & operator =(position const & b) {memcpy(xyz,b.xyz,sizeof(Storage::real)*3); return *this;}
		position() {memset(xyz,0,sizeof(Storage::real)*3);}
		position(Storage::real _xyz[3])  {memcpy(xyz,_xyz,sizeof(Storage::real)*3);}
		position(const position & b) {memcpy(xyz,b.xyz,sizeof(Storage::real)*3);}
		operator Storage::real *() {return xyz;}
		Storage::real & operator [](int k){return xyz[k];}
		Storage::real operator [](int k) const{return xyz[k];}
	};
	//2d point with comparison operator for line sweep algorithm
	struct Point
	{
		double x, y;
		Point & operator = (Point const & b) {  x = b.x; y = b.y; return *this; }
		Point(const Point & b) : x(b.x), y(b.y) {}
		Point(double _x, double _y) : x(_x), y(_y) {}
		bool operator <(const Point & b) const
		{
			if (y < b.y - 1.0e-9) return true;
			else if (y > b.y + 1.0e-9) return false;
			else if (x < b.x - 1.0e-9) return true;
			else return false;
		}
		bool operator ==(const Point & b) const
		{
			return fabs(y - b.y) < 1.0e-9 && fabs(x - b.x) < 1.0e-9;
		}
		bool operator !=(const Point & b) const
		{
			return fabs(y - b.y) > 1.0e-9 || fabs(x - b.x) > 1.0e-9;
		}
	};
	//Comparator for map of nodes
	class position_less
	{
	public:
		bool operator()(const position & a, const position & b)
		{
			for(int k = 0; k < 3; ++k)
			{
				if( a[k] < b[k] - 1.0e-9)
					return true;
				else if( a[k] > b[k] + 1.0e-9 )
					return false;
			}
			return false;
		}
	};
	//Comparator for events in line sweep algorithm
	class event_less
	{
	public:
		bool operator()(const std::pair<double, int> & a, const std::pair<double, int> & b)
		{
			if (a.first < b.first - 1.0e-9)
				return true;
			else if (a.first > b.first + 1.0e-9)
				return false;
			else if (a.second < b.second)
				return true;
			return false;
		}
	};

	//intersect a pair of segments
	std::pair<bool,Node> intersect_segments(Mesh * m, const Edge & a, const Edge & b, std::map<position,Node,position_less> & intersections, int comp, bool print)
	{
		Storage::real_array abeg = a->getBeg()->Coords();
		Storage::real_array aend = a->getEnd()->Coords();
		Storage::real_array bbeg = b->getBeg()->Coords();
		Storage::real_array bend = a->getEnd()->Coords();
		position find;
		Storage::real div = (abeg[2] - aend[2])*(bbeg[comp] - bend[comp]) - (abeg[comp] - aend[comp])*(bbeg[2] - bend[2]), t1,t2;
		if (fabs(div) < 1.0e-13)
		{
			if (print) std::cout << "divisor is zero" << std::endl;
			return std::make_pair(false,InvalidNode());
		}
		find[2]    = ((abeg[2]*aend[comp] - abeg[comp]*aend[2])*(bbeg[2]    - bend[2]   ) - (abeg[2]    - aend[2]   )*(bbeg[2]*bend[comp] - bbeg[comp]*bend[2])) / div;
		find[comp] = ((abeg[2]*aend[comp] - abeg[comp]*aend[2])*(bbeg[comp] - bend[comp]) - (abeg[comp] - aend[comp])*(bbeg[2]*bend[comp] - bbeg[comp]*bend[2])) / div;
		if (print) std::cout << "found ("<< (comp ? "y ":"x ") << find[comp] << ",z " << find[2] << ")" << std::endl;
		//probably some of these tests are redundant
		if (fabs(aend[2] - abeg[2]) > 1.0e-9)
		{
			t1 = (find[2] - abeg[2]) / (aend[2] - abeg[2]);
			if (t1 < 1.0e-9 || t1 > 1.0 - 1.0e-9)  { if (print) std::cout << "out of bound: " << t1 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		if (fabs(aend[comp] - abeg[comp]) > 1.0e-9)
		{
			t1 = (find[comp] - abeg[comp]) / (aend[comp] - abeg[comp]);
			if (t1 < 1.0e-9 || t1 > 1.0 - 1.0e-9)  { if (print) std::cout << "out of bound: " << t1 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		if (fabs(bend[2] - bbeg[2]) > 1.0e-9)
		{
			t2 = (find[2] - bbeg[2]) / (bend[2] - bbeg[2]);
			if (t2 < 1.0e-9 || t2 > 1.0 - 1.0e-9)  { if (print) std::cout << "out of bound: " << t2 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		if (fabs(bend[comp] - bbeg[comp]) > 1.0e-9)
		{
			t2 = (find[comp] - bbeg[comp]) / (bend[comp] - bbeg[comp]);
			if (t2 < 1.0e-9 || t2 > 1.0 - 1.0e-9)  { if (print) std::cout << "out of bound: " << t2 << std::endl; return std::make_pair(false,InvalidNode()); }
		}
		//restore third coordinate
		find[!comp] = 0.5*((1-t1)*abeg[!comp]+t1*aend[!comp] + (1-t2)*bbeg[!comp]+t2*bbeg[!comp]);
		if (print) std::cout << "intersection accepted (" << find[0] << "," << find[1] << "," << find[2] << ")" << std::endl;
		Node I;
		std::map<position,Node,position_less>::iterator search = intersections.find(find);
		//check whether intersection already exists
		if( search != intersections.end() ) //there is a node!
			I = search->second;
		else //no node, create one
		{
			I = m->CreateNode(find);
			intersections.insert(std::make_pair(find,I));
		}
		return std::make_pair(true,I);
	}

	//account for intersection event
	void intersect_event(Mesh * m, int a, int b, Node I, std::vector<Edge> & segments, std::multimap<Point, int> & sweep, std::multimap<std::pair<double,int>, int,event_less> & events, std::vector<Tag> & transfer, bool print)
	{
		//remove event of ending of old segment
		{
			int rem_end_events[2];
			rem_end_events[0] = a;
			rem_end_events[1] = b;
			for (int k = 0; k < 2; ++k)
			{
				std::pair< std::multimap<std::pair<double,int>, int,event_less>::iterator, std::multimap<std::pair<double,int>,int,event_less>::iterator > del = events.equal_range(std::make_pair(segments[rem_end_events[k]]->getEnd()->Coords()[2],SEG_END)); //get all events at position of the end
				bool flag = false;
				for (std::multimap<std::pair<double,int>, int,event_less>::iterator it = del.first; it != del.second; ++it) //search over all events
				{
					if (it->first.second == SEG_END && it->second == rem_end_events[k]) //event is end of segment and segment matches current
					{
						events.erase(it); //remove that segment
						flag = true;
						break; //do not expect any more
					}
				}
				if (!flag) std::cout << "Cannot find proper ending event for segment" << std::endl;
			}
		}
		//storage for data
		std::vector< std::vector<char> > copy(transfer.size()*2);
		//memorize data
		{
			const int s[2] = {a,b};
			for(int q = 0; q < 2; ++q) //segments
				for(int k = 0; k < transfer.size(); ++k) //tags
				{
					int size = segments[s[q]].GetDataSize(transfer[k]);
					copy[k + q*transfer.size()].resize(transfer[k].GetBytesSize()*size);
					if( !copy.empty() ) segments[s[q]].GetData(transfer[k],0,size,&copy[k + q*transfer.size()][0]);
				}
		}

		ElementArray<Edge> splitted_a = Edge::SplitEdge(segments[a],ElementArray<Node>(m,1,I->GetHandle()),0);
		ElementArray<Edge> splitted_b = Edge::SplitEdge(segments[b],ElementArray<Node>(m,1,I->GetHandle()),0);
		//duplicate data
		{
			const Edge splitted[2][2] =
			{
				{splitted_a[0],splitted_a[1]},
				{splitted_b[0],splitted_b[1]},
			};
			for(int q = 0; q < 2;++q) //segments
				for(int k = 0; k < transfer.size();++k)
				{
					int size = copy[k + q*transfer.size()].size() / transfer[k].GetBytesSize();
					if( size ) for(int l = 0; l < 2; ++l) //two parts
					{
						splitted[q][l].SetDataSize(transfer[k],size);
						splitted[q][l].SetData(transfer[k],0,size,&copy[k + q*transfer.size()][0]);
					}
				}
		}
		//replace segment a by new one
		segments[a] = splitted_a[0];
		//add event of ending of old segment
		events.insert(std::make_pair(std::make_pair(I->Coords()[2],SEG_END), a));
		//put other side of segment
		segments.push_back(splitted_a[1]);
		//add event of starting of new segment
		events.insert(std::make_pair(std::make_pair(I->Coords()[2],SEG_START), (int)segments.size() - 1));
		//add event of ending of new segment
		events.insert(std::make_pair(std::make_pair(segments.back()->getEnd()->Coords()[2],SEG_END),(int)segments.size() - 1));
		//replace segment b by new one
		segments[b] = splitted_b[0];
		//add event of ending of old segment
		events.insert(std::make_pair(std::make_pair(I->Coords()[2],SEG_END), b));
		//add event of starting of new segment
		events.insert(std::make_pair(std::make_pair(I->Coords()[2],SEG_START), (int)segments.size() - 1));
		//add event of ending of new segment
		events.insert(std::make_pair(std::make_pair(segments.back()->getEnd()->Coords()[2],SEG_END), (int)segments.size() - 1));
		if (print)
		{
			std::cout << "Number of events: " << events.size() << std::endl;
			for (std::multimap<std::pair<double, int>, int,event_less>::iterator it = events.begin(); it != events.end(); ++it)
				std::cout << "x: " << it->first.first << " type " << (it->first.second == SEG_START ? "start" : "end") << " segment " << it->second << std::endl;
		}
	}

	Point make_point(Edge e, int comp)
	{
		return Point(e->getBeg()->Coords()[2],e->getBeg()->Coords()[comp]);
	}

	//intersect many segments
	void intersect(Mesh * m, std::vector<Edge> & segments, int comp, std::vector<Tag> & transfer, bool print)
	{
		std::map<position,Node,position_less> intersections;
		std::multimap<std::pair<double,int>,int,event_less> events;
		std::multimap<Point,int> sweep;
	
		if( print )
		{
			std::cout << "Input segments[" << segments.size() << "]: " << std::endl;
			for (std::vector<Edge>::iterator it = segments.begin(); it != segments.end(); ++it)
				std::cout << "[ (" << (comp ? "y " : "x ") <<  it->getBeg()->Coords()[comp] << ",z" << it->getBeg()->Coords()[2] << "), ("<< (comp ? "y " : "x ") << it->getEnd()->Coords()[comp] << ",z" << it->getEnd()->Coords()[2] << ") ] " << std::endl;
			std::cout << "Create events based on segments." << std::endl;
		}

		for (int k = 0; k < (int)segments.size(); ++k)
		{
			if (segments[k]->getBeg()->Coords()[2] > segments[k]->getEnd()->Coords()[2])
				segments[k].SwapEnds();
			events.insert(std::make_pair(std::make_pair(segments[k]->getBeg()->Coords()[2],SEG_START),k));
			events.insert(std::make_pair(std::make_pair(segments[k]->getEnd()->Coords()[2],SEG_END), k));
		}


		if (print)
		{
			std::cout << "Number of events: " << events.size() << std::endl;
			for (std::multimap<std::pair<double, int>, int,event_less>::iterator it = events.begin(); it != events.end(); ++it)
				std::cout << "x: " << it->first.first << " type " << (it->first.second == SEG_START ? "start" : "end") << " segment " << it->second << std::endl;
		
			std::cout << " Start parsing events" << std::endl;
		}
	
		while (!events.empty())
		{
			std::multimap<std::pair<double,int>,int,event_less>::iterator first = events.begin();
			int t = first->first.second;
			int s = first->second;
			events.erase(first);
			if (t == SEG_START)
			{
				if( print ) std::cout << "Segment " << s << " start" << std::endl;
				//check if there is a line with same position
				std::multimap<Point, int>::iterator ins = sweep.insert(std::make_pair(make_point(segments[s],comp), s));
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
					while ((dir ? ++iter : iter--) != (dir ? sweep.end() : sweep.begin())) //y is greater for next
					{
						if (print) std::cout << "test " << s << " with " << iter->second << std::endl;
						if (segments[s]->getBeg() != segments[iter->second]->getBeg()) //ignore same starting position
						{
							if (print) std::cout << "checking intersection" << std::endl;
							std::pair<bool,Node> I = intersect_segments(m,segments[s], segments[iter->second],intersections,comp,print);
							if (I.first)
							{
								if( print ) std::cout << "Intersection of " << s << " and " << iter->second << " at (" << (comp?"y ":"x ")  << I.second.Coords()[comp] << ",z" << I.second.Coords()[2] << ")" << std::endl;
								intersect_event(m,s, iter->second, I.second, segments, sweep, events,transfer,print);
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
				std::pair< std::multimap<Point, int>::iterator, std::multimap<Point, int>::iterator > range = sweep.equal_range(make_point(segments[s],comp));
				if( print ) std::cout << "Range distance " << std::distance(range.first,range.second) << " sweep size " << sweep.size() << std::endl;
				std::multimap<Point, int>::iterator above = range.second, below = range.first;
				bool flag = false, test = true;
				if( below-- == sweep.begin() ) test = false;
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
				if (test)
				{
					if (print) std::cout << "test " << below->second << " with " << above->second << std::endl;
					if (segments[above->second]->getBeg() != segments[below->second]->getBeg())
					{
						if (print) std::cout << "checking intersection" << std::endl;
						std::pair<bool,Node> I = intersect_segments(m, segments[below->second], segments[above->second],intersections,comp,print);
						if (I.first)
						{
							if( print ) std::cout << "Intersection of " << below->second << " and " << above->second << " at (" << (comp?"y ":"x ")  << I.second.Coords()[comp] << ",z" << I.second.Coords()[2] << ")" << std::endl;
							intersect_event(m,below->second, above->second, I.second, segments, sweep, events,transfer,print);
						}
					}
					else if (print) std::cout << "skipping segments with same starting point" << std::endl;
				}
			}
		}
	
	}


	void Mesh::LoadECL(std::string File)
	{
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
		std::vector<Storage::integer> actnum;
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
						if( xyz.empty() ) xyz.resize(2*(dims[0]+1)*(dims[1]+1));
						read_arrayf = xyz.empty() ? NULL : &xyz[0];
						numrecs = 1;
						downread = totread = 2*(dims[0]+1)*(dims[1]+1);
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
						std::cout << __FILE__ << ":" << __LINE__ << " skipped " << p << " in " << fs.back().first.second << ":" << fs.back().second << std::endl;
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
						f = fopen(rec+shift_one,"r");
						if( f == NULL )
						{
							std::cout << __FILE__ << ":" << __LINE__ << " cannot open file " << rec+shift_one << " included from "<< fs.back().first.second << " line " << fs.back().second << std::endl;
							throw BadFileName;
						}
						fs.push_back(std::make_pair(std::make_pair(f,std::string(rec+shift_one)),0));
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
					if( *(pend-1) == '/' || *p == '/' ) state = ECL_NONE; else state = ECL_SKIP_SECTION;
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
			Tag tagporo,tagperm;
			if( !poro.empty() ) tagporo = CreateTag("PORO",DATA_REAL,CELL,NONE,1);
			if( !perm.empty() ) tagperm = CreateTag("PERM",DATA_REAL,CELL,NONE,3);

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
						if( !poro.empty() ) c->RealDF(tagporo) = poro[(i*dims[1]+j)*dims[2]+k];
						if( !perm.empty() )
						{
							Storage::real_array arr_perm = c->RealArrayDF(tagperm);
							arr_perm[0] = perm[3*((i*dims[1]+j)*dims[2]+k)+0];
							arr_perm[1] = perm[3*((i*dims[1]+j)*dims[2]+k)+1];
							arr_perm[2] = perm[3*((i*dims[1]+j)*dims[2]+k)+2];
						}
					}
		}
		else if( gtype == ECL_GTYPE_ZCORN )
		{
			//make arrays more human-readable
			std::vector< std::vector< std::vector< std::vector< Storage::real > > > > zcorn_array, coords_array;
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
					for(int q = 0; q < 2; q++)
						for(int j = 0; j < dims[1]; j++)
							for(int m = 0; m < 2; m++)
								for(int i = 0; i < dims[0]; i++)
									for(int l = 0; l < 2; l++)
										zcorn_array[i][j][k][l+m*2+(1-q)*4] = zcorn[pos++];
				}
			}
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
			//assemble pillars
			{
				Tag node_number = CreateTag("NODE_NUMBER",DATA_INTEGER,NODE,NONE);
				Tag edge_number = CreateTag("EDGE_NUMBER",DATA_INTEGER,NODE,NONE);
				Tag block_number = CreateTag("BLOCK_NUMBER",DATA_INTEGER,EDGE|NODE,NONE);
				typedef std::map<Storage::real,Node> pillar;
				std::vector< pillar > pillars((dims[0]+1)*(dims[1]+1));
				for(int i = 0; i < dims[0]; i++)
				{
					for(int j = 0; j < dims[1]; j++)
					{
						for(int k = 0; k < dims[2]; k++)
						{
							for(int l = 0; l < 8; ++l)
							{
								//pillar indexes
								int i2 = i + l%2, j2 = j + (l/2)%2;
								pillar & p = pillars[i2*(dims[1]+1) + j2];
								pillar::iterator search = p.find(zcorn_array[i][j][k][l]);
								Node n; //node added to pillar
								if( search == p.end() )
								{
									Storage::real node_xyz[3];
									Storage::real alpha = (zcorn_array[i][j][k][l]-coords_array[i2][j2][1][2])/(coords_array[i2][j2][0][2]-coords_array[i2][j2][1][2]);
									node_xyz[0] = coords_array[i2][j2][1][0] + alpha * (coords_array[i2][j2][0][0] - coords_array[i2][j2][1][0]);
									node_xyz[1] = coords_array[i2][j2][1][1] + alpha * (coords_array[i2][j2][0][1] - coords_array[i2][j2][1][1]);
									node_xyz[2] = zcorn_array[i][j][k][l];
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
									p.insert(std::make_pair(zcorn_array[i][j][k][l],n));
								} else n = search->second; //if
								n->IntegerArray(node_number).push_back(l);
								n->IntegerArray(block_number).push_back((i + (j+k*dims[1])*dims[0]));
							} //l
						} //k
					} //j
				} //i
				/*      (6)*--------[3]--------*(7)
						  /|                  /|
						[6]                 [7]|
						/  |                /  |
					(4)*--------[2]--------*(5)|
					   |  [10]             |  [11]
					   |                   |   |
					   |   |               |   |
					   |                   |   |
					  [8]  |              [9]  |
					   |(2)*- - - - [1]- - |- -*(3)
					   |  /                |  /
					   |[4]                |[5]
					   |/                  |/
					(0)*--------[0]--------*(1)
				*/
				//structure to create an edge from pair of nodes
				ElementArray<Node> edge_nodes(this,2);
				//all edges along pillars
				//std::vector< std::vector<Edge> > pillar_edges((dims[0]+1)*(dims[1]+1));
				//all edges of each block
				std::vector< std::vector<Edge> > block_edges(dims[0]*dims[1]*dims[2]);
				//create edges along pillars
				for(int i = 0; i < dims[0]+1; ++i)
				{
					for(int j = 0; j < dims[1]+1; ++j)
					{
						pillar & p = pillars[i*(dims[1]+1)+j];
						//std::vector<Edge> & p_edges = pillar_edges[i*(dims[1]+1)+j];
						pillar::iterator it = p.begin();
						pillar::iterator pre_it = it++;
						while(it != p.end())
						{
							edge_nodes[0] = it->second;
							edge_nodes[1] = pre_it->second;
							//p_edges.push_back(CreateEdge(edge_nodes).first);
							pre_it = it++;
						}
					} //j
				} //i
				//mark edges along pillars (vertical edges, 8, 9, 10, 11)
				for(int i = 0; i < dims[0]; ++i)
				{
					for(int j = 0; j < dims[1]; ++j)
					{
						for(int k = 0; k < dims[2]; ++k)
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
								assert(zcorn_array[i][j][k][edge_nodes[l][0]] < zcorn_array[i][j][k][edge_nodes[l][1]]); //otherwise run iterator in different direction
								b = p.find(zcorn_array[i][j][k][edge_nodes[l][0]]); //bottom
								t = p.find(zcorn_array[i][j][k][edge_nodes[l][1]]); //top
								//go over edges
								it = b;
								while(it != t)
								{
									jt = it++; // jt<->it forms a segment
									find[0] = jt->second->GetHandle();
									find[1] = it->second->GetHandle();
									Edge e(this,FindSharedAdjacency(find,2)); //find edge
									e->IntegerArray(block_number).push_back((i + (j+k*dims[1])*dims[0]));
									e->IntegerArray(edge_number).push_back(8+l);
									block_edges[(i + (j+k*dims[1])*dims[0])].push_back(e);
								}
							}
						}
					}
				}
				//set of lines to be intersected
				std::vector<Edge> edges;
				//tags to be transfered
				std::vector<Tag> transfer(2);
				transfer[0] = edge_number;
				transfer[1] = block_number;
				//go over nx pairs of pillars, then ny pairs of pillars
				for(int q = 0; q < 2; ++q)
				{
					for(int i = 0; i < dims[0]+q; i++)
					{
						for(int j = 0; j < dims[1]+!q; ++j)
						{
							//p0 is at i,j
							//p1 is at i+1,j, when q = 0 and at i,j+1, when q = 1
							pillar & p0 = pillars[(i+ 0)*(dims[1]+1)+j+0];
							pillar & p1 = pillars[(i+!q)*(dims[1]+1)+j+q];
							//preallocate array for edges
							edges.reserve(2*std::max(p0.size(),p1.size()));
							//add far faces of j-1 block
							if( (1-q)*j + q*i > 0 ) //test j when q = 0 and i when q = 1
							{
								for(int k = 0; k < dims[2]; ++k)
								{
									for(int l = 0; l < 2; ++l) //top-bottom
									{
										//q = 0 - test 2,3 + 4*l (edge 1 and 3)
										//q = 1 - test 1,3 + 4*l (edge 5 and 7)
										edge_nodes[0] = p0[zcorn_array[i-q][j-!q][k][2 - q + 4*l]];
										edge_nodes[1] = p1[zcorn_array[i-q][j-!q][k][3 - 0 + 4*l]];
										Edge e = CreateEdge(edge_nodes).first;
										edges.push_back(e);
										e->IntegerArray(block_number).push_back((i + (j+k*dims[1])*dims[0]));
										e->IntegerArray(edge_number).push_back(1 + 2*l + 4*q);
									}
								}
							}
							//add near faces of j block
							if( (1-q)*j + q*i  < dims[!q] ) //test j when q = 0 and i when q = 1
							{
								for(int k = 0; k < dims[2]; ++k)
								{
									for(int l = 0; l < 2; ++l) //top-bottom
									{
										//q = 0 - test 0,1 + 4*l (edge 0 and 2)
										//q = 1 - test 0,2 + 4*l (edge 4 and 6)
										edge_nodes[0] = p0[zcorn_array[i][j][k][0 + 0 + 4*l]];
										edge_nodes[1] = p1[zcorn_array[i][j][k][1 + q + 4*l]];
										Edge e = CreateEdge(edge_nodes).first;
										edges.push_back(e);
										e->IntegerArray(block_number).push_back((i + (j+k*dims[1])*dims[0]));
										e->IntegerArray(edge_number).push_back(0 + 2*l + 4*q);
									}
								}
							}
							//produce intersected edges
							intersect(this,edges,!q,transfer,false);
							//add edges to block
							for(int r = 0; r < (int)edges.size(); ++r)
							{
								Storage::integer_array blocks = edges[r].IntegerArray(block_number);
								for(int t = 0; t < blocks.size(); ++t)
									block_edges[blocks[t]].push_back(edges[r]);
							} //q
						} //j
					} //i
				} //q



			}
		}
	}
}

#endif