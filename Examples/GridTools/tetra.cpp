#include "inmost.h"
#include <stdio.h>
#include <deque>
//#include "tetgen/tetgen.h"
using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

const bool divide_cells = true;

real dot(real a[3], real b[3])
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void cross(real a[3], real b[3], real out[3])
{
	out[0] = a[1]*b[2] - a[2]*b[1];
	out[1] = a[2]*b[0] - a[0]*b[2];
	out[2] = a[0]*b[1] - a[1]*b[0];		
}

struct node
{
	real coo[2];
	node(const real _coo[2])
	{
		coo[0] = _coo[0];
		coo[1] = _coo[1];
	}
	node(const node & b)
	{
		coo[0] = b.coo[0];
		coo[1] = b.coo[1];
	}
	node & operator = (node const & b)
	{
		coo[0] = b.coo[0];
		coo[1] = b.coo[1];
		return *this;
	}
};

struct segment
{
	real beg[2];
	real end[2];
	segment(const real _beg[2], const real _end[2])
	{
		for(int k = 0; k < 2; ++k)
		{
			beg[k] = _beg[k];
			end[k] = _end[k];
		}
	}
	segment(const segment & b)
	{
		for(int k = 0; k < 2; ++k)
		{
			beg[k] = b.beg[k];
			end[k] = b.end[k];
		}
	}
	segment & operator = (segment const & b)
	{
		for(int k = 0; k < 2; ++k)
		{
			beg[k] = b.beg[k];
			end[k] = b.end[k];
		}
		return *this;
	}
};


//check whether two segments intersect
bool intersect_segments(const segment & a, const segment & b, bool print = false)
{
	const real eps = 1.0e-9;
	real div1 = (a.beg[0] - a.end[0])*(b.beg[1] - b.end[1]);
	real div2 = (a.beg[1] - a.end[1])*(b.beg[0] - b.end[0]);
	real div = div1 - div2;
	real t, find[2];
	if( print )
	{
		std::cout << "(" << a.beg[0] << "," << a.beg[1] << ",0) <-> (" << a.end[0] << "," << a.end[1] << ",0)" << std::endl;
		std::cout << "(" << b.beg[0] << "," << b.beg[1] << ",0) <-> (" << b.end[0] << "," << b.end[1] << ",0)" << std::endl;
	}
	if( print ) std::cout << "div = " << div << " div1 " << div1 << " div2 " << div2 << " against " << eps*(std::abs(div1) + std::abs(div2)) << std::endl;
	if (std::abs(div) <= eps*(std::abs(div1) + std::abs(div2))) //segments are parallel
	{
		//return true;
		if( print ) std::cout << "parallel segments, div = " << div << std::endl;
		if( std::abs(a.beg[0]-a.end[0]) < eps*(std::abs(a.beg[0])+std::abs(a.end[0])) ) //parallel to x
		{
			real a0 = a.beg[1], a1 = a.end[1];
			real b0 = b.beg[1], b1 = b.end[1];
			if( print ) std::cout << "parallel to x, projection on Oy: a [" << a0 << "," << a1 << "], b [" << b0 << "," << b1 << "]" << std::endl;
			if( a0 > a1 ) std::swap(a0,a1);
			if( b0 > b1 ) std::swap(b0,b1);
			if( a1 > b1 ) //change a and b
			{
				std::swap(a0,b0);
				std::swap(a1,b1);
			}
			if( a1 > b0+eps )
			{
				if( print ) std::cout << "intersection" << std::endl;
				return true;
			}
			if( print ) std::cout << "no intersection" << std::endl;
			return false;
		}
		//reconstruct k in y = kx+b
		real ak = (a.end[1]-a.beg[1])/(a.end[0]-a.beg[0]);
		real bk = (b.end[1]-b.beg[1])/(b.end[0]-b.beg[0]);
		//reconstruct bs in y = kx+b
		real ab = (a.end[0]*a.beg[1] - a.beg[0]*a.end[1])/(a.end[0]-a.beg[0]);
		real bb = (b.end[0]*b.beg[1] - b.beg[0]*b.end[1])/(b.end[0]-b.beg[0]);
		if( print )
		{
			std::cout << "parallel lines:" << std::endl;
			std::cout << "a: y = " << ak << "*x + " << ab << std::endl;
			std::cout << "b: y = " << bk << "*x + " << bb << std::endl;
			std::cout << "kdiff: " << ak - bk << std::endl;
			std::cout << "bdiff: " << ab - bb << " against " << eps*(std::abs(ab)+std::abs(bb)) <<  std::endl;
		}
		if( std::abs(ab-bb) <= eps*(std::abs(ab)+std::abs(bb)) )
		{
			//test segments
			int ind = ak > 1.0 ? 1 : 0; //choose projection of segments
			real a0 = a.beg[ind], a1 = a.end[ind];
			real b0 = b.beg[ind], b1 = b.end[ind];
			if( print ) std::cout << "projection on O" << (ind?"y":"x") << ": a [" << a0 << "," << a1 << "], b [" << b0 << "," << b1 << "]" << std::endl;
			if( a0 > a1 ) std::swap(a0,a1);
			if( b0 > b1 ) std::swap(b0,b1);
			if( a1 > b1 ) //change a and b
			{
				std::swap(a0,b0);
				std::swap(a1,b1);
			}
			if( a1 > b0+eps )
			{
				if( print ) std::cout << "intersection" << std::endl;
				return true;
			}
			if( print ) std::cout << "no intersection" << std::endl;
			return false;
		}
		//printf("divisor is zero\n");
		return false;
	}
	find[0] = ((a.beg[0]*a.end[1] - a.beg[1]*a.end[0])*(b.beg[0] - b.end[0]) - (a.beg[0] - a.end[0])*(b.beg[0]*b.end[1] - b.beg[1]*b.end[0])) / div;
	find[1] = ((a.beg[0]*a.end[1] - a.beg[1]*a.end[0])*(b.beg[1] - b.end[1]) - (a.beg[1] - a.end[1])*(b.beg[0]*b.end[1] - b.beg[1]*b.end[0])) / div;
	if( print ) std::cout << "intersection: (" << find[0] << "," << find[1] << ")" << std::endl;
	//probably some of these tests are redundant
	bool inter = true;
	if (std::abs(a.end[0] - a.beg[0]) > eps)
	{
		t = (find[0] - a.beg[0]) / (a.end[0] - a.beg[0]);
		if( print ) std::cout << "coef a on Ox " << t << std::endl;
		if (t < eps || t > 1.0 - eps)  inter = false;
	}
	if (std::abs(a.end[1] - a.beg[1]) > eps)
	{
		t = (find[1] - a.beg[1]) / (a.end[1] - a.beg[1]);
		if( print ) std::cout << "coef a on Oy " << t << std::endl;
		if (t < eps || t > 1.0 - eps)  inter = false;
	}
	if (std::abs(b.end[0] - b.beg[0]) > eps)
	{
		t = (find[0] - b.beg[0]) / (b.end[0] - b.beg[0]);
		if( print ) std::cout << "coef b on Ox " << t << std::endl;
		if (t < eps || t > 1.0 - eps)  inter = false;
	}
	if (std::abs(b.end[1] - b.beg[1]) > eps)
	{
		t = (find[1] - b.beg[1]) / (b.end[1] - b.beg[1]);
		if( print ) std::cout << "coef b on Oy " << t << std::endl;
		if (t < eps || t > 1.0 - eps)  inter = false;
	}
	if( inter )
	{
		if( print ) std::cout << "intersection" << std::endl;
		return true;
	}
	if( print ) std::cout << "no intersection" << std::endl;
	return false;
}

//check whether any of the segments from setb intersect any segments from seta
bool check_intersect(const std::vector<segment> & segmentsa, const std::vector<segment> & segmentsb, bool print = false)
{
	if( print )
	{
		std::cout << "segments [" << segmentsa.size() << "]:" << std::endl;
		for(size_t k = 0; k < segmentsa.size(); ++k)
			std::cout << "(" << segmentsa[k].beg[0] << "," << segmentsa[k].beg[1] << ",0) <-> (" << segmentsa[k].end[0] << "," << segmentsa[k].end[1] << ",0)" << std::endl;
		std::cout << "joints [" << segmentsb.size() << "]:" << std::endl;
		for(size_t k = 0; k < segmentsb.size(); ++k)
			std::cout << "(" << segmentsb[k].beg[0] << "," << segmentsb[k].beg[1] << ",0) <-> (" << segmentsb[k].end[0] << "," << segmentsb[k].end[1] << ",0)" << std::endl;
	}
	for(size_t i = 0; i < segmentsb.size(); ++i)
	{
		for(size_t j = 0; j < segmentsa.size(); ++j)
		{
			if( print ) std::cout << "check segment " << j << " and joint " << i << std::endl;
			if( intersect_segments(segmentsa[j],segmentsb[i],print) )
			{
				if( print ) std::cout << "there is intersection of segment " << j << " and joint " << i << std::endl;
				return true;
			}
		}
	}
	if( print ) std::cout << "there is no intersection" << std::endl;
	return false;	
}

//check whether any of the centers of setb drop outside of loop represented by seta
bool check_dropout(const std::vector<node> & loop, const std::vector<segment> & segmentsb)
{
	for(int i = 0; i < (int)segmentsb.size(); ++i)
	{
		real cnt[2];
		cnt[0] = (segmentsb[i].beg[0]+segmentsb[i].end[0])*0.5;
		cnt[1] = (segmentsb[i].beg[1]+segmentsb[i].end[1])*0.5;

		int    cn = 0;    // the  crossing number counter
		// loop through all edges of the polygon
		for (size_t j = 0; j < loop.size()-1; j++)  // edge from V[i]  to V[i+1]
		{
			real px = cnt[0];
			real py = cnt[1];
			real v0x = loop[j].coo[0];
			real v0y = loop[j].coo[1];
			real v1x = loop[j+1].coo[0];
			real v1y = loop[j+1].coo[1];
			if (   ((v0y <= py) && (v1y >  py)) // an upward crossing
				|| ((v0y  > py) && (v1y <= py)))// a downward crossing
			{
					// compute  the actual edge-ray intersect x-coordinate
					real vt = (py  - v0y) / (v1y - v0y);
					if (px <  v0x + vt * (v1x - v0x)) // P.x < intersect
						++cn;   // a valid crossing of y=P.y right of P.x
			}
		}
		if( !(cn&1) ) //out of polygon
			return true;
	}
	return false;
}


bool Planarity(Face f, double tol)
{
	real nrm[3], cnt[3], area = f->Area();
	f->Centroid(cnt);
	f->UnitNormal(nrm);
	ElementArray<Node> nodes = f->getNodes();
	for(ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		real_array c = it->Coords();
		if( fabs((c[0]-cnt[0])*nrm[0] + (c[1]-cnt[1])*nrm[1] + (c[2]-cnt[2])*nrm[2]) > tol*area )
			return false;
	}
	return true;
}


int main(int argc, char ** argv)
{
	if( argc < 2 ) 
	{
		printf("Usage: %s input_mesh [output_mesh]\n",argv[0]);
		return -1;
	}

	Mesh m;
	m.SetFileOption("ECL_CURVILINEAR","FALSE");
	m.SetFileOption("ECL_SPLIT_GLUED","TRUE");
	m.Load(argv[1]);
	
	Mesh::GeomParam table;
	table[ORIENTATION] = FACE;
	m.PrepareGeometricData(table);
	//reorder nodes according to RCM
	Tag order = m.CreateTag("ORDER",DATA_INTEGER,NODE,NONE,1);
	Tag idnum = m.CreateTag("RCM_ID",DATA_INTEGER,NODE,NONE,1);
	for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
	{
		it->Integer(order) = (int)it->BridgeAdjacencies2Node(EDGE).size();
		it->Integer(idnum) = -1;
	}
	
	int id = 0;
	while( id < m.NumberOfNodes() )
	{
		Node select = InvalidNode();
		Mesh::iteratorNode kt;
		for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
		{
			if( it->Integer(idnum) == -1 )
			{
				select = it->self();
				kt = it;
				break;
			}
		}
		for(Mesh::iteratorNode it = kt; it != m.EndNode(); ++it)
		{
			if( it->Integer(idnum) == -1 && it->Integer(order) < select->Integer(order) )
				select = it->self();
		}
		std::deque<Node> q;
		q.push_back(select);
		select->Integer(idnum) = id++;
		while(!q.empty())
		{
			select = q.front();
			q.pop_front();
			ElementArray<Node> adjacent = select->BridgeAdjacencies2Node(EDGE);
			std::sort(adjacent.begin(),adjacent.end(),Mesh::IntegerComparator(&m,order));
			for(ElementArray<Node>::iterator it = adjacent.begin(); it != adjacent.end(); ++it)
			{
				if( it->Integer(idnum) == -1 )
				{
					it->Integer(idnum) = id++;
					q.push_back(it->self());
				}
			}
		}
	}
	ElementArray<Node> nodes_rcm(&m,m.NumberOfNodes());
	for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
		nodes_rcm[m.NumberOfNodes() - it->Integer(idnum) - 1] = it->self();
	m.DeleteTag(order);
	m.DeleteTag(idnum);

	//choose nodes to divide each cell
	Tag divnode = m.CreateTag("DIVIDING_NODE",DATA_REFERENCE,CELL|FACE,NONE,1);
	Tag proj_coords = m.CreateTag("PROJECTED_COORDS",DATA_REAL,NODE,NONE,2);
	Tag mat = m.CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);
	MarkerType sharenode = m.CreateMarker();
	MarkerType centernode = m.CreateMarker();
	std::vector<segment> segments, joints;
	std::vector<node> loop;
	std::map<HandleType,real> hits;
	real nrm[3], cnt[3], ncnt[3], scnt[3], pcnt[3], fcnt[3], ray[3], orthx[3], orthy[3], d, nd;
	//for(Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
	for(ElementArray<Node>::iterator it = nodes_rcm.begin(); it != nodes_rcm.end(); ++it)
	{
		ElementArray<Face> node_faces = it->getFaces();
		ElementArray<Edge> node_edges = it->getEdges();
		ElementArray<Cell> cells = it->getCells();
		it->SetMarker(sharenode);
		node_faces.SetMarker(sharenode);
		for(ElementArray<Cell>::iterator jt = cells.begin(); jt != cells.end(); ++jt) if( jt->GetGeometricType() != Element::Tet ) //skip tets
		{
			//retrive shared faces
			ElementArray<Face> faces = jt->getFaces();
			ElementArray<Edge> edges = jt->getEdges();
			//check that cell or some face is not already divided
			bool fail = false;
			if( jt->Reference(divnode) != InvalidHandle() ) fail = true;
			for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end() && !fail; ++kt)
				if( kt->GetMarker(sharenode) && kt->GetGeometricType() != Element::Tri && ((kt->Reference(divnode) != InvalidHandle() && kt->Reference(divnode) != *it) /*|| !Planarity(kt->self(),0.0001)*/) ) fail = true;
			if( !fail )
			{
				//resolve a wired degenerate case
				//mark edges of faces not shared by current node
				node_edges.SetMarker(sharenode);
				for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end(); ++kt) if( !kt->GetMarker(sharenode) )
					kt->getEdges().SetMarker(sharenode);
				//count that all edges are marked
				int count = 0;
				for(ElementArray<Edge>::iterator kt = edges.begin(); kt != edges.end(); ++kt) if( kt->GetMarker(sharenode) ) count++;
				edges.RemMarker(sharenode);
				node_edges.RemMarker(sharenode);
				//not all edges were marked - degenerate case
				if( count != (int)edges.size() ) fail = true;
				//check that all the edges are visible from the node
				for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end() && !fail; ++kt) if( kt->GetMarker(sharenode) && kt->GetGeometricType() != Element::Tri ) //skip tris
				{
					if( kt->Reference(divnode) == *it ) continue; //already chosen
					kt->UnitNormal(nrm);
					kt->Centroid(cnt);
					d = dot(cnt,nrm);
					ElementArray<Node> face_nodes = kt->getNodes();
					//project coordinates
					for(ElementArray<Node>::iterator mt = face_nodes.begin(); mt != face_nodes.end(); ++mt)
					{
						mt->Centroid(ncnt);
						nd = dot(ncnt,nrm)-d;
						//project onto plane
						pcnt[0] = ncnt[0] - nrm[0]*nd;
						pcnt[1] = ncnt[1] - nrm[1]*nd;
						pcnt[2] = ncnt[2] - nrm[2]*nd;
						//relative to center
						pcnt[0] -= cnt[0];
						pcnt[1] -= cnt[1];
						pcnt[2] -= cnt[2];
						//choose the basis
						if( mt == face_nodes.begin() )
						{
							orthx[0] = pcnt[0];
							orthx[1] = pcnt[1];
							orthx[2] = pcnt[2];
							cross(orthx,nrm,orthy);
						}
						//compute and record the coords
						real_array c = mt->RealArray(proj_coords);
						c[0] = dot(orthx,pcnt);
						c[1] = dot(orthy,pcnt);
					}
					//gather all the segments in projected coordinates
					ElementArray<Edge> face_edges = kt->getEdges();
					for(ElementArray<Edge>::iterator mt = face_edges.begin(); mt != face_edges.end(); ++mt)
						segments.push_back(segment(mt->getBeg().RealArray(proj_coords).data(),mt->getEnd().RealArray(proj_coords).data()));
					//check that the node sees centers of all segments
					for(ElementArray<Edge>::iterator mt = face_edges.begin(); mt != face_edges.end(); ++mt)
					{
						if( mt->getBeg() != it->self() && mt->getEnd() != it->self() ) 	//skip adjacent to the node
						{
							scnt[0] = (mt->getBeg()->RealArray(proj_coords)[0]+mt->getEnd()->RealArray(proj_coords)[0])*0.5;
							scnt[1] = (mt->getBeg()->RealArray(proj_coords)[1]+mt->getEnd()->RealArray(proj_coords)[1])*0.5;
							joints.push_back(segment(it->RealArray(proj_coords).data(),scnt));
						}
					}
					for(ElementArray<Node>::iterator mt = face_nodes.begin(); mt != face_nodes.end(); ++mt)
						loop.push_back(node(mt->RealArray(proj_coords).data()));
					loop.push_back(loop.front());
					if( check_intersect(segments,joints) || check_dropout(loop,joints) )//,kt->LocalID() == 206357) )
						fail = true;
					segments.clear();
					joints.clear();
					loop.clear();
				}
				if( !fail )
				{
					//check that all the faces, other then shared with the node, are visible from the node
					it->Centroid(cnt);
					for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end() && !fail; ++kt) if( !kt->GetMarker(sharenode) )
					{
						kt->Centroid(fcnt);
						ray[0] = fcnt[0]-cnt[0];
						ray[1] = fcnt[1]-cnt[1];
						ray[2] = fcnt[2]-cnt[2];
						it->CastRay(cnt,ray,hits);
						int counter = 0;
						for(std::map<HandleType,real>::iterator qt = hits.begin(); qt != hits.end(); ++qt)
							if( qt->second > 1.0e-5 ) counter++;
						if( counter > 1 ) //expect only one hit
							fail = true;
						hits.clear();
					}
					if( !fail )
					{
						//select node as the divisor for faces and cell
						for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end(); ++kt) if( kt->GetMarker(sharenode) && kt->GetGeometricType() != Element::Tri )
							kt->Reference(divnode) = *it;
						jt->Reference(divnode) = *it;
					}
				}
			}
		}
		node_faces.RemMarker(sharenode);
		it->RemMarker(sharenode);
	}
	//all faces that do not get divisor node should be attempted to be divided from their center
	int nfcenter = 0;
	int nffail = 0;
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) if( it->GetGeometricType() != Element::Tet && it->Reference(divnode) == InvalidHandle() )
	{
		bool success = false;
		ElementArray<Node> face_nodes = it->getNodes();
		it->UnitNormal(nrm);
		it->Centroid(cnt);
		d = dot(cnt,nrm);
		//project coordinates
		for(ElementArray<Node>::iterator mt = face_nodes.begin(); mt != face_nodes.end(); ++mt)
		{
			mt->Centroid(ncnt);
			nd = dot(ncnt,nrm)-d;
			//project onto plane
			pcnt[0] = ncnt[0] - nrm[0]*nd;
			pcnt[1] = ncnt[1] - nrm[1]*nd;
			pcnt[2] = ncnt[2] - nrm[2]*nd;
			//relative to center
			pcnt[0] -= cnt[0];
			pcnt[1] -= cnt[1];
			pcnt[2] -= cnt[2];
			//choose the basis
			if( mt == face_nodes.begin() )
			{
				orthx[0] = pcnt[0];
				orthx[1] = pcnt[1];
				orthx[2] = pcnt[2];
				cross(orthx,nrm,orthy);
			}
			//compute and record the coords
			real_array c = mt->RealArray(proj_coords);
			c[0] = dot(orthx,pcnt);
			c[1] = dot(orthy,pcnt);
		}
		//gather all the segments in projected coordinates
		ElementArray<Edge> face_edges = it->getEdges();
		for(ElementArray<Edge>::iterator mt = face_edges.begin(); mt != face_edges.end(); ++mt)
			segments.push_back(segment(mt->getBeg().RealArray(proj_coords).data(),mt->getEnd().RealArray(proj_coords).data()));
		//if( Planarity(it->self(),0.0001) )
		{
			//form a loop out of face nodes
			for(ElementArray<Node>::iterator mt = face_nodes.begin(); mt != face_nodes.end() && !success; ++mt)
				loop.push_back(node(mt->RealArray(proj_coords).data()));
			loop.push_back(loop.front());
			
			//check that some face node is suitable
			for(ElementArray<Node>::iterator mt = face_nodes.begin(); mt != face_nodes.end() && !success; ++mt)
			{
				//check that the node sees centers of all segments
				for(ElementArray<Edge>::iterator qt = face_edges.begin(); qt != face_edges.end(); ++qt)
				{
					if( qt->getBeg() != mt->self() && qt->getEnd() != mt->self() )
					{
						scnt[0] = (qt->getBeg()->RealArray(proj_coords)[0]+qt->getEnd()->RealArray(proj_coords)[0])*0.5;
						scnt[1] = (qt->getBeg()->RealArray(proj_coords)[1]+qt->getEnd()->RealArray(proj_coords)[1])*0.5;
						joints.push_back(segment(mt->RealArray(proj_coords).data(),scnt));
					}
				}
				if( !check_intersect(segments,joints) && !check_dropout(loop,joints) )//,it->LocalID() == 206357) )
				{
					it->Reference(divnode) = mt->GetHandle();
					success = true;
				}
				joints.clear();
			}
		}
		//try face center
		if( !success )
		{
			//check that the node sees centers of all segments
			for(int k = 0; k < (int)segments.size(); ++k)
			{
				const real cnt0[2] = {0.0,0.0};
				scnt[0] = (segments[k].beg[0]+segments[k].end[0])*0.5;
				scnt[1] = (segments[k].beg[1]+segments[k].end[1])*0.5;
				joints.push_back(segment(cnt0,scnt));
			}
			if( check_intersect(segments,joints) ) 
			{
				std::cout << "Center of face " << it->LocalID() << " is not suitable to divide it" << std::endl;
				nffail++;
			}
			else
			{
				it->Reference(divnode) = m.CreateNode(cnt).GetHandle();
				nfcenter++;
			}
		}
		segments.clear();
		joints.clear();
		loop.clear();
	}
	m.DeleteTag(proj_coords);
	if( nfcenter ) std::cout << "Dividing " << nfcenter << " faces using it's center." << std::endl;
	if( nffail ) std::cout << "Failed to choose any node to divide " << nffail << " faces" << std::endl;
	//select divisor for cells
	int nccenter = 0;
	int ncfail = 0;
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( it->GetGeometricType() != Element::Tet && it->Reference(divnode) == InvalidHandle() )
	{
		it->Centroid(cnt);
		bool fail = false;
		ElementArray<Face> faces = it->getFaces();
		for(ElementArray<Face>::iterator kt = faces.begin(); kt != faces.end() && !fail; ++kt)
		{
			kt->Centroid(fcnt);
			ray[0] = fcnt[0]-cnt[0];
			ray[1] = fcnt[1]-cnt[1];
			ray[2] = fcnt[2]-cnt[2];
			it->CastRay(cnt,ray,hits);
			int counter = 0;
			for(std::map<HandleType,real>::iterator qt = hits.begin(); qt != hits.end(); ++qt)
				if( qt->second > 1.0e-5 ) counter++;
			if( counter > 1 ) //expect only one hit
				fail = true;
			hits.clear();
		}
		if( !fail )
		{
			//select node as the divisor for cell
			it->Reference(divnode) = m.CreateNode(cnt).GetHandle();
			m.SetMarker(it->Reference(divnode),centernode);
			nccenter++;
		}
		else 
		{
			std::cout << "Center of cell " << it->LocalID() << " is not suitable to divide it" << std::endl;
			ncfail++;
		}
	}
	if( nccenter ) std::cout << "Dividing " << nccenter << " cells using it's center." << std::endl;
	if( ncfail ) std::cout << "Failed to choose any node to divide " << ncfail << " cells" << std::endl;
	//modify all faces according to selected divisors
	m.BeginModification();
	Tag wasid = m.CreateTag("WAS_ID",DATA_INTEGER,CELL|FACE,NONE,1);
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) it->Integer(wasid) = it->LocalID();
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) it->Integer(wasid) = it->LocalID();
	ElementArray<Edge> new_edges(&m);
	ElementArray<Node> edge_nodes(&m,2);
	int total = m.NumberOfFaces(), k = 0;
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) 
	{
		int id = it->LocalID();
		if( it->GetGeometricType() != Element::Tri && !it->New() && it->Reference(divnode) != InvalidHandle() )
		{
			edge_nodes.at(0) = it->Reference(divnode);
			ElementArray<Node> nodes = it->getNodes();
			for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); ++kt) if( kt->self() != edge_nodes[0] )
			{
				edge_nodes[1] = kt->self();
				std::pair<Edge,bool> edge = m.CreateEdge(edge_nodes);
				if( edge.second ) //this is actually a new edge
					new_edges.push_back(edge.first);
			}
			ElementArray<Face> split_faces = Face::SplitFace(it->self(),new_edges,0);
			for(ElementArray<Face>::iterator kt = split_faces.begin(); kt != split_faces.end(); ++kt)
			{
				kt->Integer(wasid) = id;
				/*
				if( m.TopologyErrorTag().isValid() && kt->HaveData(m.TopologyErrorTag() ) )
					std::cout << "face " << kt->LocalID() << " got topology error, was id " << id << std::endl;
				if( kt->GetGeometricType() != Element::Tri )
				{
					ElementArray<Edge> old_edges = it->getEdges();
					std::cout << "face " << kt->LocalID() << " is " << Element::GeometricTypeName(kt->GetGeometricType()) << ", was id " << id << std::endl;
					std::cout << "old edges[" << old_edges.size() << "]:" << std::endl;
					for(ElementArray<Edge>::iterator mt = old_edges.begin(); mt != old_edges.end(); ++mt)
						std::cout << "(" << mt->getBeg()->Coords()[0] << "," << mt->getBeg()->Coords()[1] << "," << mt->getBeg()->Coords()[2] << ") <-> (" << mt->getEnd()->Coords()[0] << "," << mt->getEnd()->Coords()[1] << "," << mt->getEnd()->Coords()[2] << ")" << std::endl;
					std::cout << "new edges[" << new_edges.size() << "]:" << std::endl;
					for(ElementArray<Edge>::iterator mt = new_edges.begin(); mt != new_edges.end(); ++mt)
						std::cout << "(" << mt->getBeg()->Coords()[0] << "," << mt->getBeg()->Coords()[1] << "," << mt->getBeg()->Coords()[2] << ") <-> (" << mt->getEnd()->Coords()[0] << "," << mt->getEnd()->Coords()[1] << "," << mt->getEnd()->Coords()[2] << ")" << std::endl;

				}
				*/
			}
			new_edges.clear();
			//transfer data
		}
		k++;
		if( k % 25 == 0 )
		{
			printf("%06.2f%% %10d / %10d\r",(double)(k+1)*100.0/(double)total,k,total);
			fflush(stdout);
		}
		if( k >= total ) break; //others are new ones
	}
	m.ApplyModification();
	m.EndModification();
	std::cout << std::endl << "Done with faces" << std::endl;
	//modify all cells according to selected divisors
	//m.BeginModification();
	if( divide_cells )
	{
		ElementArray<Face> new_faces(&m);
		ElementArray<Edge> new_internal_edges(&m);
		total = m.NumberOfCells(), k = 0;
		Tag tetgenid = m.CreateTag("TETGEN_ID",DATA_INTEGER,NODE,NONE); //used to provide local node enumeration for tetgen
		Tag bytetgen = m.CreateTag("BY_TETGEN",DATA_INTEGER,CELL,NONE);
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			int id = it->LocalID();
			if( it->GetGeometricType() != Element::Tet && !it->New() && it->Reference(divnode) != InvalidHandle() )
			{
				bool use_tetgen = false;
				Node div(&m,it->Reference(divnode));
				//check that all the faces of the cell are still visible
				ElementArray<Face> node_faces = div->getFaces();
				node_faces.SetMarker(sharenode);
				ElementArray<Face> cell_faces = it->getFaces();
				div->Centroid(cnt);
				
				for(ElementArray<Face>::iterator jt = cell_faces.begin(); jt != cell_faces.end(); ++jt) if( !jt->GetMarker(sharenode) )
				{
					jt->Centroid(fcnt);
					ray[0] = fcnt[0]-cnt[0];
					ray[1] = fcnt[1]-cnt[1];
					ray[2] = fcnt[2]-cnt[2];
					it->CastRay(cnt,ray,hits);
					int counter = 0;
					for(std::map<HandleType,real>::iterator qt = hits.begin(); qt != hits.end(); ++qt)
						if( qt->second > 1.0e-5 ) counter++;
					if( counter > 1 ) //expect only one hit
					{
						std::cout << "FACE:" << jt->LocalID() << " of CELL:" << it->LocalID() << " hits: " << counter << std::endl;
						for(std::map<HandleType,real>::iterator qt = hits.begin(); qt != hits.end(); ++qt) //if( qt->second > 1.0e-5 )
							std::cout << "intersect " << ElementTypeName(GetHandleElementType(qt->first)) << ":" << GetHandleID(qt->first) << " at " << qt->second << std::endl;

						//use_tetgen = true;
					}
					hits.clear();
				}
				
				node_faces.RemMarker(sharenode);
				if( !use_tetgen )
				{
				//~ backup:
					//mark all edges that emerge from the divisor node
					ElementArray<Edge> node_edges = div.getEdges();
					node_edges.SetMarker(sharenode);
					//gather all the edges of the cell that do not contain divisor node
					ElementArray<Edge> cell_edges = it->getEdges(sharenode,1);
					node_edges.RemMarker(sharenode);
					edge_nodes[0] = div;
					for(ElementArray<Edge>::iterator jt = cell_edges.begin(); jt != cell_edges.end(); ++jt)
					{
						new_edges.push_back(jt->self());
						ElementArray<Node> nodes = jt->getNodes();
						for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); ++kt)
						{
							edge_nodes[1] = kt->self();
							std::pair<Edge,bool> e = m.CreateEdge(edge_nodes);
							new_edges.push_back(e.first);
							if( e.second) new_internal_edges.push_back(e.first);
						}
						std::pair<Face,bool> face = m.CreateFace(new_edges);
						if( face.second ) //this is actually a new face
						{
							face.first.FixEdgeOrder();
							new_faces.push_back(face.first);
						}
						new_edges.clear();
					}
					ElementArray<Cell> split_cells = Cell::SplitCell(it->self(),new_faces,0);
					if( div.GetMarker(centernode) )
					{
						for(ElementArray<Cell>::iterator it = split_cells.begin(); it != split_cells.end(); ++it)
							it->Integer(mat) = 2;
					}
					else
					{
						for(ElementArray<Cell>::iterator it = split_cells.begin(); it != split_cells.end(); ++it)
							it->Integer(mat) = 1;
					}
					for(ElementArray<Cell>::iterator kt = split_cells.begin(); kt != split_cells.end(); ++kt)
					{
						kt->Integer(bytetgen) = 0;
						kt->Integer(wasid) = id;
						
						// if( m.TopologyErrorTag().isValid() && kt->HaveData(m.TopologyErrorTag() ) )
						// std::cout << "cell " << kt->LocalID() << " got topology error, was id " << id << std::endl;
						 if( kt->GetGeometricType() != Element::Tet )
						 {
							 ElementArray<Edge> old_edges = it->getEdges();
							 std::cout << "cell " << kt->LocalID() << " is " << Element::GeometricTypeName(kt->GetGeometricType()) << ", was id " << id << " splitted into " << split_cells.size() << " cells, new faces " << new_faces.size() << std::endl;
							 
						 std::cout << "cells: "; for(ElementArray<Cell>::iterator mt = split_cells.begin(); mt != split_cells.end(); ++mt) std::cout << Element::GeometricTypeName(mt->GetGeometricType()) << " "; std::cout << std::endl;
						 std::cout << "faces: "; for(ElementArray<Face>::iterator mt = new_faces.begin(); mt != new_faces.end(); ++mt) std::cout << Element::GeometricTypeName(mt->GetGeometricType()) << "(" << mt->nbAdjElements(NODE) << ") "; std::cout << std::endl;
						 std::cout << "old edges[" << old_edges.size() << "]:" << std::endl;
						 for(ElementArray<Edge>::iterator mt = old_edges.begin(); mt != old_edges.end(); ++mt)
							std::cout << "(" << mt->getBeg()->Coords()[0] << "," << mt->getBeg()->Coords()[1] << "," << mt->getBeg()->Coords()[2] << ") <-> (" << mt->getEnd()->Coords()[0] << "," << mt->getEnd()->Coords()[1] << "," << mt->getEnd()->Coords()[2] << ")" << std::endl;
						 std::cout << "new internal edges[" << new_internal_edges.size() << "]:" << std::endl;
						 for(ElementArray<Edge>::iterator mt = new_internal_edges.begin(); mt != new_internal_edges.end(); ++mt)
							std::cout << "(" << mt->getBeg()->Coords()[0] << "," << mt->getBeg()->Coords()[1] << "," << mt->getBeg()->Coords()[2] << ") <-> (" << mt->getEnd()->Coords()[0] << "," << mt->getEnd()->Coords()[1] << "," << mt->getEnd()->Coords()[2] << ")" << std::endl;
						 ElementArray<Edge> cedges = kt->getEdges();
						 std::cout << "element edges[" << cedges.size() << "]:" << std::endl;
						 for(ElementArray<Edge>::iterator mt = cedges.begin(); mt != cedges.end(); ++mt)
							std::cout << "(" << mt->getBeg()->Coords()[0] << "," << mt->getBeg()->Coords()[1] << "," << mt->getBeg()->Coords()[2] << ") <-> (" << mt->getEnd()->Coords()[0]<< "," << mt->getEnd()->Coords()[1] << "," << mt->getEnd()->Coords()[2] << ")" << std::endl;
							 
						 }
							  
						
					}
					new_faces.clear();
					new_internal_edges.clear();
				}
				/*
				else //use tetgen
				{
					std::cout << "Call tetgen" << std::endl;
					tetgenio input, output;
					ElementArray<Node> cell_nodes = it->getNodes();
					input.firstnumber = 0;
					input.numberofpoints = cell_nodes.size();
					input.pointlist = new REAL[cell_nodes.size()*3];
					for(int k = 0; k < (int)cell_nodes.size(); ++k)
					{
						Storage::real_array c = cell_nodes[k]->Coords();
						input.pointlist[k*3+0] = c[0];
						input.pointlist[k*3+1] = c[1];
						input.pointlist[k*3+2] = c[2];
						cell_nodes[k]->Integer(tetgenid) = k;
					}
					input.numberoffacets = (int)cell_faces.size();
					input.facetlist = new tetgenio::facet[input.numberoffacets];
					for(int k = 0; k < (int)cell_faces.size(); ++k)
					{
						ElementArray<Node> face_nodes = cell_faces[k]->getNodes();
						tetgenio::facet & f = input.facetlist[k];
						f.numberofpolygons = 1;
						f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
						f.numberofholes = 0;
						f.holelist = NULL;
						tetgenio::polygon &p = f.polygonlist[0];
						p.numberofvertices = (int)face_nodes.size();
						p.vertexlist = new int[p.numberofvertices];
						for(int l = 0; l < (int)face_nodes.size(); ++l)
							p.vertexlist[l] = face_nodes[l]->Integer(tetgenid);
					}
					try
					{
						tetrahedralize("zpQ", &input, &output);
					}
					catch(...)
					{
						std::cout << "tetgen fail!" << std::endl;
						goto backup;
					}
					//printf("numberofpoints was %d now %d\n",input.numberofpoints,output.numberofpoints);
					//printf("numberoftrifaces %d",output.numberoftrifaces);

					//assume initial points stay the same
					for(int k = input.numberofpoints; k < output.numberofpoints; ++k)
						cell_nodes.push_back(m.CreateNode(&output.pointlist[k*3]));

					ElementArray<Node> tri_nodes(&m,3);
					for(int k = 0; k < output.numberoftrifaces; ++k)
					{
						for(int l = 0; l < 3; ++l)
							tri_nodes[l] = cell_nodes[output.trifacelist[k*3+l]];
						std::pair<Face,bool> f = m.CreateFace(tri_nodes);
						if( f.second )
							new_faces.push_back(f.first);
					}
					ElementArray<Cell> split_cells = Cell::SplitCell(it->self(),new_faces,0);
					for(ElementArray<Cell>::iterator kt = split_cells.begin(); kt != split_cells.end(); ++kt)
					{
						kt->Integer(bytetgen) = 1;
						kt->Integer(wasid) = id;
					}
					new_faces.clear();
				}
				 */
			}
			k++;
			if( k % 25 == 0 )
			{
				printf("%06.2f%%  %10d / %10d\r",(double)(k+1)*100.0/(double)total,k,total);
				fflush(stdout);
			}
			if( k >= total ) break; //others are new ones
		}
		m.DeleteTag(tetgenid);
		//
		int orphan = 0;
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
			if( it->nbAdjElements(CELL) == 0 ) orphan++;
		if( orphan ) std::cout << "orphan faces: " << orphan << std::endl;
		//check that there are no unsplit cells (may happen due to non-convexity), split them with central node
		/*
		std::cout << std::endl << "Using tetgen on unfinished" << std::endl;
		total = m.NumberOfCells(), k = 0;
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( it->GetGeometricType() != Element::Tet )
		{
			int id = it->LocalID();
			std::cout << "CELL:"<<id<< " is " << Element::GeometricTypeName(it->GetGeometricType()) << " breaking down " << std::endl;
			//make central node
			real cnt[3];
			it->Centroid(cnt);
			edge_nodes[0] = m.CreateNode(cnt);
			//check that all the faces are triangles?
			bool notri = false;
			ElementArray<Face> cell_faces = it->getFaces();
			std::cout << "Faces: ";
			for(ElementArray<Face>::iterator jt = cell_faces.begin(); jt != cell_faces.end(); ++jt)
			{
				std::cout << Element::GeometricTypeName(jt->GetGeometricType()) << " ";
				if( jt->GetGeometricType() != Element::Tri ) notri = true;
			}
			std::cout << std::endl;
			if( notri ) 
			{
				std::cout << "non-triangular face: skip!" << std::endl;
				continue;
			}
			//connect with edges to make internal faces
			ElementArray<Edge> cell_edges = it->getEdges();
			for(ElementArray<Edge>::iterator jt = cell_edges.begin(); jt != cell_edges.end(); ++jt)
			{
				new_edges.push_back(jt->self());
				ElementArray<Node> nodes = jt->getNodes();
				for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); ++kt)
				{
					edge_nodes[1] = kt->self();
					std::pair<Edge,bool> e = m.CreateEdge(edge_nodes);
					new_edges.push_back(e.first);
					//if( e.second) new_internal_edges.push_back(e.first);
				}
				std::pair<Face,bool> face = m.CreateFace(new_edges);
				if( face.second ) //this is actually a new face
				{
					face.first.FixEdgeOrder();
					new_faces.push_back(face.first);
				}
				new_edges.clear();
			}
			ElementArray<Cell> split_cells = Cell::SplitCell(it->self(),new_faces,0);
			std::cout << "Result: ";
			for(ElementArray<Cell>::iterator jt = split_cells.begin(); jt != split_cells.end(); ++jt)
			{
				std::cout << Element::GeometricTypeName(jt->GetGeometricType()) << " ";
				jt->Integer(mat) = 2;
				jt->Integer(wasid) = id;
			}
			std::cout << std::endl;
			new_faces.clear();
			k++;
			if( k % 25 == 0 )
			{
				printf("%06.2f%%  %10d / %10d\r",(double)(k+1)*100.0/(double)total,k,total);
				fflush(stdout);
			}
			if( k >= total ) break; //others are new ones
		}
		 */
		std::cout << std::endl << "Done with cells" << std::endl;
	}
	std::map<Element::GeometricType,int> types;
	for(Mesh::iteratorElement it = m.BeginElement(CELL|FACE); it != m.EndElement(); ++it)
		types[it->GetGeometricType()]++;
	for(std::map<Element::GeometricType,int>::iterator it = types.begin(); it != types.end(); ++it)
		std::cout << Element::GeometricTypeName(it->first) << ": " <<it->second << std::endl;
	
	//m.ApplyModification();
	//m.EndModification();
	m.ReorderEmpty(CELL|FACE|EDGE|NODE);
	
	Tag lid = m.CreateTag("LID",DATA_INTEGER,CELL|FACE,NONE,1);
	//for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
	int q = 0;
	for(int k = 0; k < m.CellLastLocalID(); ++k )if( m.isValidCell(k) )
	{
		
		Cell it = m.CellByLocalID(k);
		it->Integer(lid) = it->LocalID();
		if( it->Volume() < 1.0e-9 )
			std::cout << "Volume for CELL:" << it->LocalID() << "(expected CELL:" << q << ") is " << it->Volume() << std::endl;
		q++;
	}
	q = 0;
	for(int k = 0; k < m.FaceLastLocalID(); ++k )if( m.isValidFace(k) )
	{
		
		Face it = m.FaceByLocalID(k);
		it->Integer(lid) = it->LocalID();
		if( it->Area() < 1.0e-9 )
			std::cout << "Area for FACE:" << it->LocalID() << "(expected FACE:" << q << ") is " << it->Area() << std::endl;
		q++;
	}
	
	m.ReleaseMarker(sharenode);
	m.ReleaseMarker(centernode,NODE);
	//m.DeleteTag(divnode);
	

	m.Save("tet.pmf");


	if( argc > 2 )
		m.Save(argv[2]);
	else
		m.Save("out.vtk");

	return 0;
}
