#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
typedef Storage::real real;
typedef Storage::real_array real_array;

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

static void normalize(real v[3])
{
	real d = sqrt(dot(v, v));
	if (d) for (int k = 0; k < 3; ++k) v[k] /= d;
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


int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		printf("Usage: %s input_mesh [output_mesh]\n",argv[0]);
		return -1;
	}
	
	Mesh A("");
	A.Load(argv[1]);
	
	Mesh::GeomParam table;
	table[ORIENTATION] = FACE;
	table[BARYCENTER] = FACE | CELL | EDGE;
	A.PrepareGeometricData(table);
	
	Tag proj_coords = A.CreateTag("PROJECTED_COORDS",DATA_REAL,NODE,NONE,2);
	MarkerType myedges = A.CreateMarker();
	real nrm[3], cnt[3], ncnt[3], scnt[3], pcnt[3], orthx[3] = {0,0,0}, orthy[3] = {0,0,0}, d, nd;
	(void)nd;
	std::cout << "Start splitting faces" << std::endl;
	int nsplit = 0, had_faces = A.NumberOfFaces();
	for(Mesh::iteratorFace it = A.BeginFace(); it != A.EndFace(); ++it)
	{
		if( !it->Planarity() )
		{
			std::vector<segment> segments;
			it->UnitNormal(nrm);
			it->Centroid(cnt);
			d = dot(cnt,nrm);
			ElementArray<Node> face_nodes = it->getNodes();
			//project coordinates
			for(ElementArray<Node>::iterator mt = face_nodes.begin(); mt != face_nodes.end(); ++mt)
			{
				mt->Centroid(ncnt);
				nd = dot(ncnt,nrm)-d;
				//project onto plane
				pcnt[0] = ncnt[0] - cnt[0];
				pcnt[1] = ncnt[1] - cnt[1];
				pcnt[2] = ncnt[2] - cnt[2];
				//choose the basis
				if( mt == face_nodes.begin() )
				{
					orthx[0] = pcnt[0];
					orthx[1] = pcnt[1];
					orthx[2] = pcnt[2];
					normalize(orthx);
					cross(orthx,nrm,orthy);
				}
				//compute and record the coords
				real_array c = mt->RealArray(proj_coords);
				c[0] = dot(orthx,pcnt);
				c[1] = dot(orthy,pcnt);
			}
			
			//gather all the segments in projected coordinates
			ElementArray<Edge> face_edges = it->getEdges();
			face_edges.SetMarker(myedges);
			for(ElementArray<Edge>::iterator mt = face_edges.begin(); mt != face_edges.end(); ++mt)
				segments.push_back(segment(mt->getBeg().RealArray(proj_coords).data(),mt->getEnd().RealArray(proj_coords).data()));
			//form a loop out of face nodes
			std::vector<node> loop;
			for(ElementArray<Node>::iterator mt = face_nodes.begin(); mt != face_nodes.end(); ++mt)
				loop.push_back(node(mt->RealArray(proj_coords).data()));
			loop.push_back(loop.front());
			//try each adjacent node
			bool success = false;
			for(ElementArray<Node>::iterator jt = face_nodes.begin(); jt != face_nodes.end() && !success; ++jt)
			{
				std::vector<segment> joints;
				//check that the node sees centers of all segments
				for(ElementArray<Edge>::iterator mt = face_edges.begin(); mt != face_edges.end(); ++mt)
				{
					if( mt->getBeg() != jt->self() && mt->getEnd() != jt->self() ) 	//skip adjacent to the node
					{
						scnt[0] = (mt->getBeg()->RealArray(proj_coords)[0]+mt->getEnd()->RealArray(proj_coords)[0])*0.5;
						scnt[1] = (mt->getBeg()->RealArray(proj_coords)[1]+mt->getEnd()->RealArray(proj_coords)[1])*0.5;
						joints.push_back(segment(jt->RealArray(proj_coords).data(),scnt));
					}
				}
				
				if( !check_intersect(segments,joints) && !check_dropout(loop,joints) )
				{
					//found a suitable node
					ElementArray<Edge> new_edges(&A);
					ElementArray<Node> edge_nodes(&A,2);
					edge_nodes[0] = jt->self();
					for(ElementArray<Node>::iterator kt = face_nodes.begin(); kt != face_nodes.end(); ++kt)
					{
						if( kt->self() != jt->self() )
						{
							edge_nodes[1] = kt->self();
							std::pair<Edge,bool> edge = A.CreateEdge(edge_nodes);
							//if( edge.second ) //this is actually a new edge
							if( !edge.first.GetMarker(myedges) )
								new_edges.push_back(edge.first);
						}
					}
					int had_id = it->LocalID();
					int had_nodes = it->nbAdjElements(NODE);
					ElementArray<Face> split_faces = Face::SplitFace(it->self(),new_edges,0);
					for(INMOST_DATA_ENUM_TYPE q = 0; q < split_faces.size(); ++q)
						if( !split_faces[q].Planarity() )
						{
							std::cout << __FILE__ << ":" << __LINE__ << "Face " << split_faces[q].LocalID() << " non-planar after split, nodes " << split_faces[q].nbAdjElements(NODE) << " original face " << had_id << " with " << had_nodes << " nodes, split by " << new_edges.size() << " edges, resulted in " << split_faces.size() << " new faces" << std::endl;
							std::cout << "original nodes:";
							for(ElementArray<Node>::iterator kt = face_nodes.begin(); kt != face_nodes.end(); ++kt)
								std::cout << " " << kt->LocalID();
							std::cout << std::endl;
							std::cout << "split node " << jt->LocalID() << std::endl;
							std::cout << "split edges:";
							for(ElementArray<Edge>::iterator kt = new_edges.begin(); kt != new_edges.end(); ++kt)
								std::cout << " " << kt->LocalID() << " (" << kt->getBeg().LocalID() << "," << kt->getEnd().LocalID() << ")";
							std::cout << std::endl;
						}
					nsplit++;
					success = true;
				}
			}
			face_edges.RemMarker(myedges);
			if( !success )
			{
				std::cout << "Face " << it->LocalID() << ", nodes " << it->nbAdjElements(NODE) << ", testing center" << std::endl;
				std::vector<segment> joints;
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
				}
				else
				{
					//found a suitable node
					ElementArray<Edge> new_edges(&A);
					ElementArray<Node> edge_nodes(&A,2);
					edge_nodes[0] = A.CreateNode(cnt);
					for(ElementArray<Node>::iterator kt = face_nodes.begin(); kt != face_nodes.end(); ++kt)
					{
						edge_nodes[1] = kt->self();
						std::pair<Edge,bool> edge = A.CreateEdge(edge_nodes);
						new_edges.push_back(edge.first);
					}
					ElementArray<Face> split_faces = Face::SplitFace(it->self(),new_edges,0);
					int had_id = it->LocalID();
					int had_nodes = it->nbAdjElements(NODE);
					for(INMOST_DATA_ENUM_TYPE q = 0; q < split_faces.size(); ++q)
						if( !split_faces[q].Planarity() )
							std::cout << __FILE__ << ":" << __LINE__ << "Face " << split_faces[q].LocalID() << " non-planar after split, nodes " << split_faces[q].nbAdjElements(NODE) << " original face " << had_id << " with " << had_nodes << " nodes, split by " << new_edges.size() << " edges, resulted in " << split_faces.size() << " new faces" << std::endl;
					nsplit++;
					success = true;
				}
			}
			if( !success ) std::cout << "Cannot split face " << it->LocalID() << std::endl;
		}
	}
	
	std::cout <<"total split: " << nsplit << " / " << had_faces << std::endl;
	A.ReleaseMarker(myedges);
	for(Mesh::iteratorFace it = A.BeginFace(); it != A.EndFace(); ++it)
		if( !it->Planarity() )
			std::cout << it->LocalID() << " is still non-planar, nodes " << it->nbAdjElements(NODE) << std::endl;
	
	
	if( A.HaveTag("GRIDNAME") )
	{
		Storage::bulk_array nameA = A.self().BulkArray(A.GetTag("GRIDNAME"));
		std::string ins = "_split_nonplanar";
		nameA.insert(nameA.end(),ins.c_str(),ins.c_str()+ins.size());
	}
	
	if( argc > 2 )
	{
		std::cout << "Save to " << argv[2] << std::endl;
		A.Save(argv[2]);
	}
	else
	{
		std::cout << "Save to out.vtk" << std::endl;
		A.Save("out.vtk");
	}

	return 0;
}
