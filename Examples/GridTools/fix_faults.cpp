#include "inmost.h"

using namespace INMOST;

//todo: non-flat faces

static void GetBox(Element e, float minmax[6])
{
	minmax[0] = minmax[2] = minmax[4] = 1.0e20f; //min
	minmax[1] = minmax[3] = minmax[5] = -1.0e20f; //max
	ElementArray<Node> nodes = e->getNodes();
	for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		Storage::real_array c = it->Coords();
		for (int i = 0; i < (int)c.size(); i++)
		{
			float v = static_cast<float>(c[i]);
			if (minmax[1 + 2 * i] < v) minmax[1 + 2 * i] = v; //max
			if (minmax[0 + 2 * i] > v) minmax[0 + 2 * i] = v; //min
		}
	}
	for (int i = 0; i < 3; ++i)
	{
		if (minmax[1 + 2 * i] < minmax[0 + 2 * i])
		{
			minmax[1 + 2 * i] = 1.0e-5f;
			minmax[0 + 2 * i] = -1.0e-5f;
		}
		else if (minmax[1 + 2 * i] == minmax[0 + 2 * i])
		{
			minmax[1 + 2 * i] += 1.0e-5f;
			minmax[0 + 2 * i] += -1.0e-5f;
		}
	}
}

static void MergeBox(const float bboxa[6], const float bboxb[6], float bbox[6])
{
	for (int k = 0; k < 3; k++)
	{
		bbox[0 + 2 * k] = std::min(bboxa[0 + 2 * k], bboxb[0 + 2 * k]);
		bbox[1 + 2 * k] = std::max(bboxa[1 + 2 * k], bboxb[1 + 2 * k]);
	}
}

static bool IntersectBox(const float bboxa[6], const float bboxb[6])
{
	for (int k = 0; k < 3; ++k)
	{
		float la = bboxa[0 + 2 * k], ra = bboxa[1 + 2 * k];
		float lb = bboxa[0 + 2 * k], rb = bboxb[1 + 2 * k];
		if (ra < lb || la > rb) return false;
	}
	return true;
}

static bool MatchPoints(const double v1[3], const double v2[3])
{
	double l = 0;
	for (int k = 0; k < 3; ++k)
		l += (v1[k] - v2[k])*(v1[k] - v2[k]);
	return sqrt(l) < 1.0e-7;
}

//returns distance between line if shortest line is within segments, writes position of distance into last argument
static double SegmentDistance(const double v1[3], const double v2[3], const double v3[3], const double v4[3], double vout[3])
{
	double v13[3], v21[3], v43[3];
	for (int k = 0; k < 3; ++k)
	{
		v13[k] = v1[k] - v3[k];
		v21[k] = v2[k] - v1[k];
		v43[k] = v4[k] - v3[k];
	}
	double d1321 = 0, d2121 = 0, d4321 = 0, d1343 = 0, d4343 = 0;
	for (int k = 0; k < 3; ++k)
	{
		d1321 += v13[k] * v21[k];
		d2121 += v21[k] * v21[k];
		d4321 += v43[k] * v21[k];
		d1343 += v13[k] * v43[k];
		d4343 += v43[k] * v43[k];
	}
	double mu1 = 0, mu2 = 0;
	/*
	double parallel = 0;
	for (int k = 0; k < 3; ++k)
	parallel += v21[k] * v43[k];
	parallel = fabs(parallel);
	parallel /= sqrt(d2121)*sqrt(d4343);
	if (fabs(parallel - 1.0) > 1.0e-6)
	*/
	if (fabs(d2121*d4343 - d4321*d4321) > 1.0e-6)
	{
		mu1 = (d1343*d4321 - d1321*d4343) / (d2121*d4343 - d4321*d4321);
		mu2 = (d1343 + mu1*d4321) / d4343;
		if (mu1 > 1) mu1 = 1;
		if (mu1 < 0) mu1 = 0;
		if (mu2 > 1) mu2 = 1;
		if (mu2 < 0) mu2 = 0;
	}
	else //parallel lines
	{
		//1d problem
		double la = 0, ra = 0;
		double lb = 0, rb = 0;
		for (int k = 0; k < 3; ++k)
		{
			la += v1[k] * v21[k];
			ra += v2[k] * v21[k];
			lb += v3[k] * v21[k];
			rb += v4[k] * v21[k];
		}
		bool invb = false;
		if (lb > rb)
		{
			std::swap(lb, rb);
			invb = true;
		}

		if (ra < lb) // (la,ra) is to the left of (lb,rb), no intersection
		{
			mu1 = 1;
			mu2 = 0;
		}
		else if (la > rb) // (la,ra) is to the right of (lb,rb), no intersection
		{
			mu1 = 0;
			mu2 = 1;
		}
		else if (lb > la) // (lb,rb) intersects to the right of or contains (la,ra)
		{
			mu2 = 0;
			mu1 = (lb - la) / (ra - la);
		}
		else if (rb < ra) // (lb,rb) intersects to the left of or contains (la,ra)
		{
			mu2 = 1;
			mu1 = (rb - la) / (ra - la);
		}
		else // (lb,rb) contains (la,ra)
		{
			mu1 = 0;
			mu2 = (la - lb) / (rb - lb);
		}
		if (invb) mu2 = 1 - mu2;
	}
	double vs1[3], vs2[3], h = 0;
	for (int k = 0; k < 3; ++k)
	{
		vs1[k] = v1[k] + mu1*(v2[k] - v1[k]);
		vs2[k] = v3[k] + mu2*(v4[k] - v3[k]);
		vout[k] = (vs1[k] + vs2[k])*0.5;
		h += (vs2[k] - vs1[k])*(vs2[k] - vs1[k]);
	}
	return sqrt(h);
}


class kdtree
{
	struct entry
	{
		HandleType e;
		float xyz[3];
		struct entry & operator =(const struct entry & other)
		{
			e = other.e;
			xyz[0] = other.xyz[0];
			xyz[1] = other.xyz[1];
			xyz[2] = other.xyz[2];
			return *this;
		}
	} *set;
	Mesh * m;
	INMOST_DATA_ENUM_TYPE size;
	float bbox[6];
	kdtree * children;
	inline static unsigned int flip(const unsigned int * fp)
	{
		unsigned int mask = -((int)(*fp >> 31)) | 0x80000000;
		return *fp ^ mask;
	}
	inline static unsigned int _0(unsigned int x)	{ return x & 0x7FF; }
	inline static unsigned int _1(unsigned int x)	{ return x >> 11 & 0x7FF; }
	inline static unsigned int _2(unsigned int x)   { return x >> 22; }
	template<int pos>
	inline static int cmpElements(const void * a, const void * b)
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[pos];
		float bd = eb->xyz[pos];
		return (ad > bd) - (ad < bd);
	}
	void radix_sort(int dim, struct entry * temp)
	{
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0, 0, sizeof(unsigned int)*kHist * 3);
		for (i = 0; i < size; i++)
		{
			unsigned int fi = flip((unsigned int *)&set[i].xyz[dim]);
			++b0[_0(fi)]; ++b1[_1(fi)]; ++b2[_2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++)
			{
				b0[kHist - 1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist - 1];
				b1[kHist - 1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist - 1];
				b2[kHist - 1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist - 1];
			}
		}
		for (i = 0; i < size; i++) temp[++b0[_0(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[++b1[_1(flip((unsigned int *)&temp[i].xyz[dim]))]] = temp[i];
		for (i = 0; i < size; i++) temp[++b2[_2(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[i] = temp[i];
	}
	void kdtree_build(int dim, int & done, int total, struct entry * temp)
	{
		if (size > 1)
		{
			if (size > 128) radix_sort(dim, temp); else
				switch (dim)
			{
				case 0: qsort(set, size, sizeof(entry), cmpElements<0>); break;
				case 1: qsort(set, size, sizeof(entry), cmpElements<1>); break;
				case 2: qsort(set, size, sizeof(entry), cmpElements<2>); break;
			}
			children = static_cast<kdtree *>(malloc(sizeof(kdtree)* 2));//new kdtree[2];
			assert(children != NULL);
			INMOST_DATA_ENUM_TYPE shift[2];
			shift[0] = 0;
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size / 2;
			children[0].m = m;
			shift[1] = size / 2;
			children[1].children = NULL;
			children[1].set = set + size / 2;
			children[1].size = size - size / 2;
			children[1].m = m;
			//procs may have collision on array temp
			//#if defined(USE_OMP)
			//#pragma omp parallel for
			//#endif
			for (int k = 0; k < 2; ++k)
				children[k].kdtree_build((dim + 1) % 3, done, total, temp + shift[k]);
			MergeBox(children[0].bbox, children[1].bbox, bbox);
		}
		else
		{
			assert(size == 1);
			GetBox(Element(m, set[0].e), bbox);
		}
	}
	void intersect_nonadj_face_sub(Face f, float fbbox[6], MarkerType mrk, ElementArray<Face> & faces_out)
	{
		if (size > 1)
		{
			if (IntersectBox(fbbox, bbox))
			{
				children[0].intersect_nonadj_face_sub(f, fbbox, mrk, faces_out);
				children[1].intersect_nonadj_face_sub(f, fbbox, mrk, faces_out);
			}
		}
		else if (IntersectBox(fbbox, bbox))
		{
			Face ff(m, set[0].e);
			//this is not directly adjacent face
			if (!ff.GetPrivateMarker(mrk))
			{
				//check two faces are almost on the same plane
				Storage::real fnrm[3], ffnrm[3];
				f.UnitNormal(fnrm);
				ff.UnitNormal(ffnrm);
				if (fabs(fnrm[0] * ffnrm[0] + fnrm[1] * ffnrm[1] + fnrm[2] * ffnrm[2]) >= 1.0 - 1.0e-3)
				{
					//todo: check that faces intersect
					int fid = f.LocalID(), ffid = ff.LocalID();
					ElementArray<Edge> fedges = f.getEdges();
					ElementArray<Edge> ffedges = ff.getEdges();
					bool intersect = false;
					for (ElementArray<Edge>::iterator it = fedges.begin(); it != fedges.end() && !intersect; ++it)
					{
						for (ElementArray<Edge>::iterator jt = ffedges.begin(); jt != ffedges.end() && !intersect; ++jt)
						{
							double vs[3];
							double * v1 = it->getBeg().Coords().data();
							double * v2 = it->getEnd().Coords().data();
							double * v3 = jt->getBeg().Coords().data();
							double * v4 = jt->getEnd().Coords().data();
							double l1 = 0, l2 = 0;
							double h = SegmentDistance(v1, v2, v3, v4, vs);
							for (int k = 0; k < 3; ++k)
							{
								l1 += (v2[k] - v1[k])*(v2[k] - v1[k]);
								l2 += (v4[k] - v3[k])*(v4[k] - v3[k]);
							}
							l1 = sqrt(l1);
							l2 = sqrt(l2);
							if (h < (l1 + l2)*1.0e-13)
								intersect = true;
						}
					}
					if (intersect)
						faces_out.push_back(ff);
				}
			}
		}
	}
	kdtree() : set(NULL), size(0), children(NULL) {}
	void clear_children() { if (children) { children[0].clear_children(); children[1].clear_children(); free(children); } }
public:
	kdtree(Mesh * m) :  m(m), children(NULL)
	{
		double tt;
		size = m->NumberOfFaces();
		assert(size > 1);
		set = new entry[size];
		INMOST_DATA_ENUM_TYPE k = 0;
		tt = Timer();
		printf("Prepearing face set.\n");
		m->ReorderEmpty(FACE);
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for (int i = 0; i < m->FaceLastLocalID(); ++i) if (m->isValidFace(i))
		{
			Face it = m->FaceByLocalID(i);
			Storage::real cnt[3];
			it->Centroid(cnt);
			int k = it->DataLocalID();
			set[k].e = it->GetHandle();
			set[k].xyz[0] = static_cast<float>(cnt[0]);
			set[k].xyz[1] = static_cast<float>(cnt[1]);
			set[k].xyz[2] = static_cast<float>(cnt[2]);
		}
		printf("Done. Time %lg\n", Timer() - tt);
		int done = 0, total = size;
		printf("Building KD-tree.\n");
		tt = Timer();
		struct entry *  temp = new entry[size];
		kdtree_build(0, done, total, temp);
		delete[] temp;
		printf("Done. Time %lg\n", Timer() - tt);
	}
	kdtree(Mesh * m, HandleType * eset, INMOST_DATA_ENUM_TYPE size) :  m(m), size(size), children(NULL)
	{
		double tt;
		assert(size > 1);
		set = new entry[size];
		tt = Timer();
		printf("Prepearing elements set.\n");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for (int k = 0; k < (int)size; k++)
		{
			set[k].e = eset[k];
			Storage::real cnt[3];
			m->GetGeometricData(set[k].e, CENTROID, cnt);
			set[k].xyz[0] = static_cast<float>(cnt[0]);
			set[k].xyz[1] = static_cast<float>(cnt[1]);
			set[k].xyz[2] = static_cast<float>(cnt[2]);
			if (k % 150 == 0)
			{
				printf("%3.1f%%\r", (k*100.0) / (size*1.0));
				fflush(stdout);
			}
		}
		printf("Done. Time %lg\n", Timer() - tt);
		int done = 0, total = size;
		printf("Building KD-tree.\n");
		tt = Timer();
		struct entry *  temp = new entry[size];
		kdtree_build(0, done, total, temp);
		delete[] temp;
		printf("Done. Time %lg\n", Timer() - tt);
	}
	void intersect_nonadj_face(Face f, ElementArray<Face> & faces_out)
	{
		faces_out.clear();
		MarkerType adj = m->CreatePrivateMarker();
		ElementArray<Face> local = f.BridgeAdjacencies2Face(NODE);
		local.push_back(f);
		local.SetPrivateMarker(adj);
		float fbbox[6];
		GetBox(f, fbbox);
		intersect_nonadj_face_sub(f, fbbox, adj, faces_out);
		local.RemPrivateMarker(adj);
		m->ReleasePrivateMarker(adj);
	}
	~kdtree()
	{
		delete[] set;
		clear_children();
	}
};

void cross(const double a[3], const double b[3], double out[3])
{
	out[0] = a[1] * b[2] - a[2] * b[1];
	out[1] = a[2] * b[0] - a[0] * b[2];
	out[2] = a[0] * b[1] - a[1] * b[0];
}

double dot(const double a[3], const double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void normalize(double v[3])
{
	double d = sqrt(dot(v, v));
	if (d) for (int k = 0; k < 3; ++k) v[k] /= d;
}


struct coords
{
	double v[3];
	int compare(const double bv[3]) const
	{
		const double eps = 1.0e-11;//pow(10,-14);
		for (int i = 0; i < 3; i++)
		{
			if (fabs(v[i] - bv[i]) > eps)
			{
				if (v[i] > bv[i])
					return 1;
				else
					return -1;
			}
		}
		return 0;
	}
	bool operator <(const coords & b) const { return compare(b.v) < 0; }
	coords(double _v[3]) { memcpy(v, _v, sizeof(double)* 3); }
	coords(const coords &b) { memcpy(v, b.v, sizeof(double)* 3); }
	coords & operator =(coords const & b) { memmove(v, b.v, sizeof(double)* 3); return *this; }
};




int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [mesh_out=grid.vtk]" << std::endl;
		return -1;
	}

	std::string grid_out = "grid.vtk";
	if (argc > 2) grid_out = std::string(argv[2]);


	Mesh m;
	m.Load(argv[1]);

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	double tt = Timer();
	//There seems to be a problem with VTU reader
	int fixed = 0;
	for (Mesh::iteratorFace f = m.BeginFace(); f != m.EndFace(); ++f)
	{
		if (f->FixNormalOrientation())
			fixed++;
	}
	std::cout << "Time to fix normals: " << Timer() - tt << std::endl;
	std::cout << "Total face normals fixed: " << fixed << std::endl;
	

	kdtree t(&m);

	TagReferenceArray new_points = m.CreateTag("EDGEPOINTS", DATA_REFERENCE, EDGE, NONE);
	TagReferenceArray new_edges = m.CreateTag("FACEEDGES", DATA_REFERENCE, FACE, NONE);
	//TagReferenceArray face2face = m.CreateTag("FACE2FACE", DATA_REFERENCE, FACE, NONE);


	std::map<coords, Node> global_new_nodes; //have to track new nodes globally :(
	int marked = 0;
	tt = Timer();
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:marked)
#endif
	for (int it = 0; it < m.FaceLastLocalID(); ++it) if (m.isValidFace(it))
	{
		Face f = m.FaceByLocalID(it);
		int fid = f.LocalID();
		ElementArray<Face> ifaces(&m);
		t.intersect_nonadj_face(f, ifaces);
		if (!ifaces.empty())
		{
			marked++;

			ElementArray<Edge> fedges = f->getEdges();
			ElementArray<Node> fnodes = f->getNodes(); //to find out which node of segment is inside of face
			//std::cout << "Face " << f->LocalID() << " intersects with " << ifaces.size() << " faces: ";
			//record how do we want to split current face based on all intersections
			for (ElementArray<Face>::iterator ff = ifaces.begin(); ff != ifaces.end(); ++ff)
			{
				int ffid = ff->LocalID();
				//face2face[f].push_back(ff->self());
				//std::cout << iface->LocalID() << " ";
				ElementArray<Edge> ffedges = ff->getEdges();
				//find out each intersection with current edges
				//we are not interested of intersections on edges of other faces
				for (ElementArray<Edge>::iterator jt = ffedges.begin(); jt != ffedges.end(); ++jt)
				{
					int jtid = jt->LocalID();
					double * v3 = jt->getBeg().Coords().data();
					double * v4 = jt->getEnd().Coords().data();
					double l2 = 0;
					for (int k = 0; k < 3; ++k)
						l2 += (v4[k] - v3[k])*(v4[k] - v3[k]);
					l2 = sqrt(l2);
					bool _intersect = false; //edge goes through some edge of current face introducing new node
					bool _parallel = false; //edge goes parallel to some edge of current face
					bool _endpoint = false; //intersection is on endpoint of jt segment
					ElementArray<Node> internodes(&m);
					for (ElementArray<Edge>::iterator it = fedges.begin(); it != fedges.end(); ++it)
					{
						int itid = it->LocalID();
						//can put this algorithm to SegmentDistance function in future
						double * v1 = it->getBeg().Coords().data();
						double * v2 = it->getEnd().Coords().data();
						double l1 = 0;
						for (int k = 0; k < 3; ++k)
							l1 += (v2[k] - v1[k])*(v2[k] - v1[k]);
						l1 = sqrt(l1);
						double v13[3], v21[3], v43[3];
						for (int k = 0; k < 3; ++k)
						{
							v13[k] = v1[k] - v3[k];
							v21[k] = v2[k] - v1[k];
							v43[k] = v4[k] - v3[k];
						}
						double d1321 = 0, d2121 = 0, d4321 = 0, d1343 = 0, d4343 = 0;
						for (int k = 0; k < 3; ++k)
						{
							d1321 += v13[k] * v21[k];
							d2121 += v21[k] * v21[k];
							d4321 += v43[k] * v21[k];
							d1343 += v13[k] * v43[k];
							d4343 += v43[k] * v43[k];
						}
						double mu1 = 0, mu2 = 0;
						/*
						double parallel = 0;
						for (int k = 0; k < 3; ++k)
						parallel += v21[k] * v43[k];
						parallel = fabs(parallel);
						parallel /= sqrt(d2121)*sqrt(d4343);
						if( fabs(parallel - 1.0) > 1.0e-6 )
						*/
						if (fabs(d2121*d4343 - d4321*d4321) > 1.0e-6)
						{
							mu1 = (d1343*d4321 - d1321*d4343) / (d2121*d4343 - d4321*d4321);
							mu2 = (d1343 + mu1*d4321) / d4343;
							//here only one or zero intersections, but we have to introduce an edge
							if (mu1 >= 0 - 1.0e-13 && mu1 <= 1 + 1.0e-13 && mu2 >= 0 - 1.0e-13 && mu2 <= 1 + 1.0e-13)
							{
								//find distance and point
								double vs1[3], vs2[3], h = 0, vv[3];
								for (int k = 0; k < 3; ++k)
								{
									vs1[k] = v1[k] + mu1*(v2[k] - v1[k]);
									vs2[k] = v3[k] + mu2*(v4[k] - v3[k]);
									vv[k] = (vs1[k] + vs2[k])*0.5;
									h += (vs2[k] - vs1[k])*(vs2[k] - vs1[k]);
								}
								if (h < (l1 + l2)*1.0e-13)
								{
									Node matchnode = InvalidNode();
									if (MatchPoints(vv, v1))
										matchnode = it->getBeg();
									else if (MatchPoints(vv, v2))
										matchnode = it->getEnd();
									if (matchnode.isValid())
									{
										bool found = false;
										for (int k = 0; k < (int)internodes.size() && !found; ++k)
										{
											if (internodes[k] == matchnode)
												found = true;
										}
										if (!found)
											internodes.push_back(matchnode);
									}
									else
									{
										_intersect = true;
										bool checkinternodes = false;
										//nodes array may be access by another processor
#if defined(USE_OMP)
#pragma omp critical
#endif
										{
											//the node may be already found on another edge
											Storage::reference_array nodes = new_points[it->self()];
											for (int k = 0; k < (int)nodes.size() && !matchnode.isValid(); ++k)
											{
												if (MatchPoints(vv, nodes[k].getAsNode().Coords().data()))
													matchnode = nodes[k].getAsNode();
											}
											if (!matchnode.isValid())
											{
												std::map<coords, Node>::iterator find = global_new_nodes.find(vv);
												//create new node, it may not be found in internodes array
												if (find != global_new_nodes.end())//&& MatchPoints(vv,find->first.v))
													matchnode = find->second;
												else if (MatchPoints(vv, v3))
													matchnode = jt->getBeg();
												else if (MatchPoints(vv, v4))
													matchnode = jt->getEnd();
												else
												{
													matchnode = m.CreateNode(vv);
													global_new_nodes[vv] = matchnode;
												}
												nodes.push_back(matchnode);
											}
											else checkinternodes = true; //the node existed before and may be found in internodes array
										}
										if (checkinternodes)
										{
											bool found = false;
											for (int k = 0; k < (int)internodes.size() && !found; ++k)
											{
												if (internodes[k] == matchnode)
													found = true;
											}
											if (!found) checkinternodes = false;
										}
										//add the point if it is not yet present in the array
										if (!checkinternodes)
										{
											internodes.push_back(matchnode);
											if (MatchPoints(vv, v3) || MatchPoints(vv, v4)) _endpoint = true;
										}
									}
								} // if h
							} //if mu
						}
						else //parallel lines
						{
							//1d problem
							double la = 0, ra = 0;
							double lb = 0, rb = 0;
							double t[2]; //relative positions of intersection points on first segment
							int nt = 0;  //number of intersection points
							for (int k = 0; k < 3; ++k)
							{
								la += v1[k] * v21[k];
								ra += v2[k] * v21[k];
								lb += v3[k] * v21[k];
								rb += v4[k] * v21[k];
							}
							bool invb = false;
							if (lb > rb)
							{
								std::swap(lb, rb);
								invb = true;
							}
							if (ra < lb) nt = 0; // (la,ra) is to the left of (lb,rb), zero intersection points
							else if (la > rb) nt = 0; // (la,ra) is to the right of (lb,rb), zero intersection points
							else if (lb > la) // (lb,rb) intersects to the right of or contained in (la,ra), one or two intersection points
							{
								mu2 = 0;
								mu1 = (lb - la) / (ra - la);
								if (rb > ra) // (lb,rb) intersects to the right of (la,ra), one intersection point
								{
									nt = 1;
									t[0] = mu1;
								}
								else // (lb,rb) is contained in (la,ra), two intersection points
								{
									nt = 2;
									t[0] = mu1;
									t[1] = (rb - la) / (ra - la);
								}
							}
							else if (rb < ra) // (lb,rb) intersects to the left of or contained in (la,ra), one or two intersection points
							{
								mu2 = 1;
								mu1 = (rb - la) / (ra - la);
								if (lb < la) // (lb,rb)  intersection to the left of (la,ra), one intersection point
								{
									nt = 1;
									t[0] = mu1;
								}
								else // (lb,rb) is contained in (la,ra), two intersection points
								{
									nt = 2;
									t[0] = mu1;
									t[1] = (lb - la) / (ra - la);
								}
							}
							else nt = 0; // (la,ra) is contained in (lb,rb), no intersection points
							if (nt)
							{
								if (invb) mu2 = 1 - mu2;
								double h = 0;
								for (int k = 0; k < 3; ++k)
								{
									double vs1 = v1[k] + mu1*(v2[k] - v1[k]);
									double vs2 = v3[k] + mu2*(v4[k] - v3[k]);
									h += (vs2 - vs1)*(vs2 - vs1);
								}
								if (h < (l1 + l2)*1.0e-13)
								{
									Storage::reference_array nodes = new_points[it->self()];
									for (int k = 0; k < nt; ++k)
									{
										//calculate point position
										double vv[3];
										for (int q = 0; q < 3; ++q)
											vv[q] = v1[q] + t[k] * (v2[q] - v1[q]);
										//if we do not coinside with the ends of current edge
										if (!MatchPoints(vv, v1) && !MatchPoints(vv, v2))
										{
											//check this node was not yet introduced
											//nodes array may be access by another processor
#if defined(USE_OMP)
#pragma omp critical
#endif
											{
												_parallel = true;
												bool found = false;
												for (int l = 0; l < (int)nodes.size() && !found; ++l)
												{
													if (MatchPoints(vv, nodes[l].getAsNode().Coords().data()))
														found = true;
												}
												if (!found)
												{
													std::map<coords, Node>::iterator find = global_new_nodes.find(vv);
													Node n;
													if (find != global_new_nodes.end())// && MatchPoints(vv,find->first.v) )
														n = find->second;
													else if (MatchPoints(vv, v3))
														n = jt->getBeg();
													else if (MatchPoints(vv, v4))
														n = jt->getEnd();
													else
													{
														n = m.CreateNode(vv);
														global_new_nodes[vv] = n;
													}
													nodes.push_back(n);
												}
											}
										}
									} // for k
								} // if h
							} //if nt
						} //if not parallel lines
					} //loop on fedges
					if (!internodes.empty())
					{
						Edge e = InvalidEdge();
						//the segment goes through current face
						if (internodes.size() % 2 == 0)
						{
							if (internodes.size() > 2)
								std::cout << __FILE__ << ":" << __LINE__ << " Seems that the face " << f.LocalID() << " is not convex" << std::endl;
							else if (internodes[0] != internodes[1])//we have and edge crossing convex polygon
								e = m.CreateEdge(internodes).first;
						}
						else
						{
							if (internodes.size() > 1)
								std::cout << __FILE__ << ":" << __LINE__ << " Seems that the face " << f.LocalID() << " is not convex" << std::endl;
							else if (!_parallel) //check v3 or v4 is inside of the face, pick the one to construct the edge
							{

								double nrm[3], cnt[3], pcnt[3], d, nd;
								double orthx[3], orthy[3];
								f.UnitNormal(nrm);
								f.Centroid(cnt);
								//compute axises for projection
								fnodes[0].Centroid(pcnt);
								d = dot(cnt, nrm);
								nd = dot(pcnt, nrm) - d;
								for (int k = 0; k < 3; ++k)
									orthx[k] = pcnt[k] - cnt[k];
								normalize(orthx);
								cross(orthx, nrm, orthy);
								//loop over two ends of the edge
								double * v34[2] = { v3, v4 };
								bool inside[2] = { true, true };
								int cn[2] = { 0, 0 };
								for (int q = 0; q < 2; ++q)
								{
									double * v = v34[q];
									nd = dot(v, nrm) - d;
									for (int k = 0; k < 3; ++k)
										pcnt[k] = v[k] - cnt[k];
									//project current point
									double px = dot(pcnt, orthx), py = dot(pcnt, orthy);
									double v0x, v0y, v1x, v1y;
									// loop through all edges of the polygon
									for (int j = 0; j < fnodes.size(); j++)
									{
										v = fnodes[j].Coords().data();
										nd = dot(v, nrm) - d;
										for (int k = 0; k < 3; ++k)
											pcnt[k] = v[k] - cnt[k];
										v0x = dot(orthx, pcnt);
										v0y = dot(orthy, pcnt);
										v = fnodes[(j + 1) % fnodes.size()].Coords().data();
										nd = dot(v, nrm) - d;
										for (int k = 0; k < 3; ++k)
											pcnt[k] = v[k] - cnt[k];
										v1x = dot(orthx, pcnt);
										v1y = dot(orthy, pcnt);
										if ((v0y <= py && v1y >= py) || (v0y >= py && v1y <= py))
										{
											double vt = (py - v0y) / (v1y - v0y);
											double ptx = v0x + vt * (v1x - v0x);
											if (px <= ptx) ++cn[q];
										}
									}
									if (!(cn[q] & 1)) inside[q] = false;
								}
								bool report = false;
								if (inside[0] && inside[1])
								{
									report = true;
									std::cout << __FILE__ << ":" << __LINE__ << " Both points cannot be inside of polygon" << std::endl;
								}
								else if (!inside[0] && !inside[1] && !_endpoint)
								{
									report = true;
									std::cout << __FILE__ << ":" << __LINE__ << " Both points cannot be outside of polygon" << std::endl;
								}
								else if (inside[0]) // v3 is inside
									internodes.push_back(jt->getBeg());
								else if (inside[1]) // v4 is inside
									internodes.push_back(jt->getEnd());

								if (internodes.size() == 2 && internodes[0] != internodes[1])
									e = m.CreateEdge(internodes).first;
								else if (report)
									std::cout << __FILE__ << ":" << __LINE__ << " Cannot create edge due to some problem, internodes " << internodes.size() << " face " << f.LocalID() << " with " << ff->LocalID() << std::endl;
							} // one node in internodes array
						} // uneven size of internodes array
						if (e.isValid())
						{
							Storage::reference_array edges = new_edges[f];
							bool found = false;
							for (int k = 0; k < (int)edges.size() && !found; ++k)
							{
								if (edges[k] == e)
									found = true;
							}
							if (!found) new_edges[f].push_back(e);
						}
					} //there are points in internodes array
				} //loop on ffedges
			} //loop on ifaces
		} // if not empty
	} //loop mesh faces
	std::cout << "time to mark: " << Timer() - tt << std::endl;

	std::cout << "global new nodes: " << global_new_nodes.size() << std::endl;
	for (std::map<coords, Node>::iterator it = global_new_nodes.begin(); it != global_new_nodes.end(); ++it)
	{
		Storage::real_array c = it->second->Coords();
		std::cout << it->second->LocalID() << " " << it->first.v[0] << "," << it->first.v[1] << "," << it->first.v[2] << " " << c[0] << "," << c[1] << "," << c[2] << std::endl;
	}

	tt = Timer();
	TagInteger mark = m.CreateTag("SPLITTED", DATA_INTEGER, FACE | EDGE, NONE, 1);
	//start splitting
	m.BeginModification();
	//split all edges
	for (Mesh::iteratorEdge it = m.BeginEdge(); it != m.EndEdge(); ++it) if (!it->New())
	{
		Storage::reference_array nodes = new_points[it->self()];
		if (!nodes.empty())
		{
			//sort nodes according to distance from getBeg to getEnd
			std::vector< std::pair<double, HandleType> > tnodes;
			double *v1 = it->getBeg().Coords().data();
			double *v2 = it->getEnd().Coords().data();
			double v21[3];
			for (int k = 0; k < 3; ++k) v21[k] = v2[k] - v1[k];
			for (int k = 0; k < (int)nodes.size(); ++k) tnodes.push_back(std::make_pair(dot(v21, nodes[k].getAsNode().Coords().data()), nodes[k].GetHandle()));
			std::sort(tnodes.begin(), tnodes.end());
			ElementArray<Node> split_nodes(&m);
			for (int k = 0; k < (int)tnodes.size(); ++k)
				split_nodes.push_back(tnodes[k].second);
			ElementArray<Edge> new_edges = Edge::SplitEdge(it->self(), split_nodes, 0);
			for (int k = 0; k < (int)new_edges.size(); ++k)
				mark[new_edges[k]] = 1;

		}
	}
	//split all faces
	for (Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it) if (!it->New())
	{
		Storage::reference_array edges = new_edges[it->self()];
		if (!edges.empty())
		{
			ElementArray<Face> new_faces = Face::SplitFace(it->self(), ElementArray<Edge>(&m, edges.begin(), edges.end()), 0);
			for (int k = 0; k < (int)new_faces.size(); ++k)
				mark[new_faces[k]] = 1;
		}
	}
	m.EndModification();


	std::cout << "time to split: " << Timer() - tt << std::endl;

	std::cout << "Fix faults:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	m.Save(grid_out);
	return 0;
}