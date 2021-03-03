#include "clipper.h"
#include "inc_glut.h"
#include "color_bar.h"


#define CLIP_NONE 0
#define CLIP_NODE 1
#define CLIP_FULL 2
#define CLIP_ENDP 3

#define CLIP_FACE_NONE      0
#define CLIP_FACE_INSIDE    1
#define CLIP_FACE_OUTSIDE   2
#define CLIP_FACE_INTERSECT 3

namespace INMOST
{
	kdtree::entry::entry(const entry & other)
	{
		e = other.e;
		xyz[0] = other.xyz[0];
		xyz[1] = other.xyz[1];
		xyz[2] = other.xyz[2];
	}
	kdtree::entry & kdtree::entry::operator =(const entry & other)
	{
		e = other.e;
		xyz[0] = other.xyz[0];
		xyz[1] = other.xyz[1];
		xyz[2] = other.xyz[2];
		return *this;
	}
	static int cmpElements0(const void * a, const void * b)
	{
		const kdtree::entry * ea = ((const kdtree::entry *)a);
		const kdtree::entry * eb = ((const kdtree::entry *)b);
		float ad = ea->xyz[0];
		float bd = eb->xyz[0];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements1(const void * a, const void * b)
	{
		const kdtree::entry * ea = ((const kdtree::entry *)a);
		const kdtree::entry * eb = ((const kdtree::entry *)b);
		float ad = ea->xyz[1];
		float bd = eb->xyz[1];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements2(const void * a, const void * b)
	{
		const kdtree::entry * ea = ((const kdtree::entry *)a);
		const kdtree::entry * eb = ((const kdtree::entry *)b);
		float ad = ea->xyz[2];
		float bd = eb->xyz[2];
		return (ad > bd) - (ad < bd);
	}

	inline static unsigned int flip(const unsigned int * fp) { unsigned int mask = -((int)(*fp >> 31)) | 0x80000000; return *fp ^ mask; }
	inline static unsigned int _0(unsigned int x)	{ return x & 0x7FF; }
	inline static unsigned int _1(unsigned int x)	{ return x >> 11 & 0x7FF; }
	inline static unsigned int _2(unsigned int x)   { return x >> 22; }




	Storage::integer clip_plane_edge(Storage::real sp0[3], Storage::real sp1[3], Storage::real p[3], Storage::real n[3], Storage::real node[3])
	{
		Storage::real u[3], w[3], D, N, sI;
		u[0] = sp1[0] - sp0[0]; u[1] = sp1[1] - sp0[1]; u[2] = sp1[2] - sp0[2];
		w[0] = sp0[0] - p[0];   w[1] = sp0[1] - p[1];   w[2] = sp0[2] - p[2];
		D = (n[0] * u[0] + n[1] * u[1] + n[2] * u[2]);
		N = -(n[0] * w[0] + n[1] * w[1] + n[2] * w[2]);
		if (fabs(D) < 1.0e-9)
		{
			if (fabs(N) < 1.0e-9) return CLIP_FULL;
			else return CLIP_NONE;
		}
		else
		{
			sI = N / D;
			if (sI < 0 - 1.0e-9 || sI > 1 + 1.0e-9) return CLIP_NONE;
			else
			{
				node[0] = sp0[0] + sI * u[0];
				node[1] = sp0[1] + sI * u[1];
				node[2] = sp0[2] + sI * u[2];
				return (sI > 1.0e-9 && sI < 1.0 - 1.0e-9) ? CLIP_NODE : CLIP_ENDP;
			}
		}
	}


	void kdtree::radix_sort(int dim, struct entry * temp)
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

	void kdtree::kdtree_build(int dim, int & done, int total, struct entry * temp)
	{
		if (size > 1)
		{
			if (size > 128) radix_sort(dim, temp); else
				switch (dim)
			{
				case 0: qsort(set, size, sizeof(entry), cmpElements0); break;
				case 1: qsort(set, size, sizeof(entry), cmpElements1); break;
				case 2: qsort(set, size, sizeof(entry), cmpElements2); break;
			}
			children = static_cast<kdtree *>(malloc(sizeof(kdtree)* 2));//new kdtree[2];
			assert(children != NULL);
			children[0].marked = 0;
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size / 2;
			children[0].m = m;
			children[1].marked = 0;
			children[1].children = NULL;
			children[1].set = set + size / 2;
			children[1].size = size - size / 2;
			children[1].m = m;
			children[0].kdtree_build((dim + 1) % 3, done, total, temp);
			children[1].kdtree_build((dim + 1) % 3, done, total, temp);
			for (int k = 0; k < 3; k++)
			{
				bbox[0 + 2 * k] = std::min(children[0].bbox[0 + 2 * k], children[1].bbox[0 + 2 * k]);
				bbox[1 + 2 * k] = std::max(children[0].bbox[1 + 2 * k], children[1].bbox[1 + 2 * k]);
			}
		}
		else
		{
			assert(size == 1);
			if (GetHandleElementType(set[0].e) == EDGE)
			{
				Storage::real_array n1 = Edge(m, set[0].e)->getBeg()->Coords();
				Storage::real_array n2 = Edge(m, set[0].e)->getEnd()->Coords();
				for (int k = 0; k < 3; k++)
				{
					bbox[0 + 2 * k] = std::min(n1[k], n2[k]);
					bbox[1 + 2 * k] = std::max(n1[k], n2[k]);
				}
				done++;
				if (done % 150 == 0)
				{
					printf("%3.1f%%\r", (done*100.0) / (total*1.0));
					fflush(stdout);
				}
			}
			else
			{
				ElementArray<Node> nodes = Element(m, set[0].e)->getNodes();
				bbox[0] = bbox[2] = bbox[4] = 1.0e20;
				bbox[1] = bbox[3] = bbox[5] = -1.0e20;
				for (INMOST_DATA_ENUM_TYPE k = 0; k < nodes.size(); ++k)
				{
					Storage::real_array coords = nodes[k].Coords();
					for (INMOST_DATA_ENUM_TYPE q = 0; q < 3; q++)
					{
						bbox[q * 2 + 0] = std::min<float>(bbox[q * 2 + 0], coords[q]);
						bbox[q * 2 + 1] = std::max<float>(bbox[q * 2 + 1], coords[q]);
					}
				}
			}
		}
	}

	kdtree::kdtree() : set(NULL), marked(0), size(0), children(NULL) {}

	inline int kdtree::plane_bbox(Storage::real p[3], Storage::real n[3]) const
	{
		Storage::real pv[3], nv[3];
		for (int k = 0; k < 3; ++k)
		{
			if (n[k] >= 0)
			{
				pv[k] = bbox[1 + 2 * k]; //max
				nv[k] = bbox[0 + 2 * k]; //min
			}
			else
			{
				pv[k] = bbox[0 + 2 * k]; //min
				nv[k] = bbox[1 + 2 * k]; //max
			}
		}
		Storage::real pvD, nvD;
		pvD = n[0] * (pv[0] - p[0]) + n[1] * (pv[1] - p[1]) + n[2] * (pv[2] - p[2]);
		nvD = n[0] * (nv[0] - p[0]) + n[1] * (nv[1] - p[1]) + n[2] * (nv[2] - p[2]);
		if (nvD*pvD <= 0.0)
			return 2;
		else if (nvD < 0.0)
			return 1;
		else return 0;
	}

	bool kdtree::sub_intersect_plane_edge(Tag clip_point, Tag clip_state, ElementArray<Cell> & cells, MarkerType mrk, Storage::real p[3], Storage::real n[3])
	{
		if (size == 1)
		{
			assert(GetHandleElementType(set[0].e) == EDGE);
			Edge ee = Edge(m, set[0].e);
			Storage::real_array sp0 = ee->getBeg()->Coords();
			Storage::real_array sp1 = ee->getEnd()->Coords();
			Storage::integer & clip = m->IntegerDF(set[0].e, clip_state);
			clip = clip_plane_edge(&sp0[0], &sp1[0], p, n, &ee->RealArrayDF(clip_point)[0]);
			if (clip)
			{
				ElementArray<Cell> ecells = ee->getCells();
				for (INMOST_DATA_ENUM_TYPE k = 0; k < ecells.size(); ++k) if (!ecells[k].GetMarker(mrk))
				{
					ecells[k].SetMarker(mrk);
					cells.push_back(ecells[k]);
				}
				marked = 1;
			}
		}
		else if (plane_bbox(p, n) == 2)
		{
			bool test1 = children[0].sub_intersect_plane_edge(clip_point, clip_state, cells, mrk, p, n);
			bool test2 = children[1].sub_intersect_plane_edge(clip_point, clip_state, cells, mrk, p, n);
			if (test1 || test2) marked = 1;
		}
		return marked != 0;
	}

	void kdtree::sub_intersect_plane_faces(Tag clip_state, Storage::real p[3], Storage::real n[3])
	{
		if (size == 1)
		{
			Storage::integer state;
			Element ee(m, set[0].e);
			assert(ee->GetElementDimension() == 2);
			ElementArray<Node> nodes = ee->getNodes();
			Storage::real_array coords = nodes[0].Coords();
			Storage::real dot0 = n[0] * (coords[0] - p[0]) + n[1] * (coords[1] - p[1]) + n[2] * (coords[2] - p[2]);
			if (dot0 <= 0.0) state = CLIP_FACE_INSIDE; else state = CLIP_FACE_OUTSIDE;
			for (INMOST_DATA_ENUM_TYPE k = 1; k < nodes.size(); k++)
			{
				coords = nodes[k].Coords();
				Storage::real dot = n[0] * (coords[0] - p[0]) + n[1] * (coords[1] - p[1]) + n[2] * (coords[2] - p[2]);
				if (dot*dot0 <= 0.0)
				{
					state = CLIP_FACE_INTERSECT;
					break;
				}
			}
			m->IntegerDF(set[0].e, clip_state) = state;
		}
		else
		{
			marked = plane_bbox(p, n);
			if (marked == 0)
			{
				for (INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) m->IntegerDF(set[k].e, clip_state) = CLIP_FACE_OUTSIDE;
			}
			else if (marked == 1)
			{
				for (INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) m->IntegerDF(set[k].e, clip_state) = CLIP_FACE_INSIDE;
			}
			else
			{
				children[0].sub_intersect_plane_faces(clip_state, p, n);
				children[1].sub_intersect_plane_faces(clip_state, p, n);
			}
		}
	}

	void kdtree::unmark_old_edges(Tag clip_state)
	{
		if (size == 1)
		{
			assert(GetHandleElementType(set[0].e) == EDGE);
			marked = 0;
			if (GetHandleElementType(set[0].e) == EDGE)
				m->IntegerDF(set[0].e, clip_state) = CLIP_NONE;
			else if (GetHandleElementType(set[0].e) == FACE)
				m->IntegerDF(set[0].e, clip_state) = CLIP_FACE_NONE;
		}
		else if (children)
		{
			if (children[0].marked) { children[0].unmark_old_edges(clip_state); marked = 0; }
			if (children[1].marked) { children[1].unmark_old_edges(clip_state); marked = 0; }
		}
	}

	void kdtree::clear_children() { if (children) { children[0].clear_children(); children[1].clear_children(); free(children); } }

	kdtree::kdtree(Mesh * m) : marked(0), m(m), children(NULL)
	{
		double tt;
		size = m->NumberOfEdges();
		assert(size > 1);
		set = new kdtree::entry[size];
		INMOST_DATA_ENUM_TYPE k = 0;
		tt = Timer();
		printf("Prepearing edge set.\n");
		for (Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it)
		{
			set[k].e = *it;
			set[k].xyz[0] = (it->getBeg()->Coords()[0] + it->getEnd()->Coords()[0])*0.5;
			set[k].xyz[1] = (it->getBeg()->Coords()[1] + it->getEnd()->Coords()[1])*0.5;
			set[k].xyz[2] = (it->getBeg()->Coords()[2] + it->getEnd()->Coords()[2])*0.5;
			k++;
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
		for (int k = 0; k < 3; k++)
		{
			bbox[0 + 2 * k] = std::min(children[0].bbox[0 + 2 * k], children[1].bbox[0 + 2 * k]);
			bbox[1 + 2 * k] = std::max(children[0].bbox[1 + 2 * k], children[1].bbox[1 + 2 * k]);
		}
	}

	kdtree::kdtree(Mesh * m, HandleType * eset, INMOST_DATA_ENUM_TYPE size) : marked(0), m(m), size(size), children(NULL)
	{
		double tt;
		assert(size > 1);
		set = new entry[size];
		tt = Timer();
		printf("Prepearing elements set.\n");
		for (INMOST_DATA_ENUM_TYPE k = 0; k < size; k++)
		{
			set[k].e = eset[k];
			Storage::real cnt[3];
			m->GetGeometricData(set[k].e, CENTROID, cnt);
			set[k].xyz[0] = cnt[0];
			set[k].xyz[1] = cnt[1];
			set[k].xyz[2] = cnt[2];
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
		for (int k = 0; k < 3; k++)
		{
			bbox[0 + 2 * k] = std::min(children[0].bbox[0 + 2 * k], children[1].bbox[0 + 2 * k]);
			bbox[1 + 2 * k] = std::max(children[0].bbox[1 + 2 * k], children[1].bbox[1 + 2 * k]);
		}
	}

	void kdtree::intersect_plane_edge(Tag clip_point, Tag clip_state, ElementArray<Cell> & cells, MarkerType mark_cells, Storage::real p[3], Storage::real n[3])
	{
		if (marked)
		{
			unmark_old_edges(clip_state);
			cells.clear();
		}
		sub_intersect_plane_edge(clip_point, clip_state, cells, mark_cells, p, n);
	}

	void kdtree::intersect_plane_face(Tag clip_state, Storage::real p[3], Storage::real n[3])
	{
		sub_intersect_plane_faces(clip_state, p, n);
	}

	kdtree::~kdtree()
	{
		delete[] set;
		clear_children();
	}

	clipper::edge_point::edge_point(Storage::real _xyz[3], Storage::integer n, float v)
	{
		xyz[0] = _xyz[0];
		xyz[1] = _xyz[1];
		xyz[2] = _xyz[2];
		edge = n;
		val = v;
	}
	clipper::edge_point::edge_point(){}

	bool clipper::edge_point::operator ==(const edge_point& b) const
	{
		Storage::real temp = 0.0;
		for (int k = 0; k < 3; k++) temp += (xyz[k] - b.xyz[k])*(xyz[k] - b.xyz[k]);
		if (temp < 1.0e-8) return true; else return false;
	}

	bool clipper::edge_point::operator !=(const edge_point& b) const { return !(operator ==(b)); }

	void clipper::edge_point::print() { printf("%g %g %g e %d\n", xyz[0], xyz[1], xyz[2], (int)edge); }

	clipper::~clipper()
	{
		delete tree;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k++) cells[k]->RemMarker(marker);
		mm->ReleaseMarker(marker);
		mm->DeleteTag(clips);
		mm->DeleteTag(clipsv);
		mm->DeleteTag(clipsn);
		mm->DeleteTag(clip_point);
		mm->DeleteTag(clip_state);
	}

	clipper::clipper(Mesh * m)
	{
		mm = m;
		cells.SetMeshLink(mm);
		tree = new kdtree(m);
		marker = m->CreateMarker();
		clips = m->CreateTag("CLIPS", DATA_REAL, CELL, CELL);
		clipsv = m->CreateTag("CLIPS_VAL", DATA_REAL, CELL, CELL);
		clipsn = m->CreateTag("CLIPS_NUM", DATA_INTEGER, CELL, CELL);
		clip_point = m->CreateTag("CLIP_POINT", DATA_REAL, EDGE, NONE, 3);
		clip_state = m->CreateTag("CLIP_STATE", DATA_INTEGER, EDGE, NONE, 1);
	}

	Storage::real clipper::compute_value(Edge e, Storage::real * pnt)
	{
		if (isColorBarEnabled())
		{
			Storage::real_array c1 = e->getBeg()->Coords();
			Storage::real_array c2 = e->getEnd()->Coords();
			Storage::real d1, d2, t;
			d1 = sqrt((pnt[0] - c1[0])*(pnt[0] - c1[0]) + (pnt[1] - c1[1])*(pnt[1] - c1[1]) + (pnt[2] - c1[2])*(pnt[2] - c1[2]));
			d2 = sqrt((c2[0] - c1[0])*(c2[0] - c1[0]) + (c2[1] - c1[1])*(c2[1] - c1[1]) + (c2[2] - c1[2])*(c2[2] - c1[2]));
			t = d1 / d2; //(pnt == c2, t = 1 : pnt == c1, t = 0)
			return e->getBeg()->RealDF(GetVisualizationTag())*(1 - t) + e->getEnd()->RealDF(GetVisualizationTag())*t;
		}
		else return 0.f;
	}

	void clipper::clip_plane(Storage::real p[3], Storage::real n[3])
	{
		const bool print = false;

		for (INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); ++k)
		if (cells[k]->GetMarker(marker))
		{
			cells[k]->RealArray(clips).clear();
			cells[k]->RealArray(clipsv).clear();
			cells[k]->IntegerArray(clipsn).clear();
			cells[k]->RemMarker(marker);
		}
		tree->intersect_plane_edge(clip_point, clip_state, cells, marker, p, n);
		std::vector<edge_point> clipcoords, loopcoords;
		std::vector<bool> closed;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); ++k)
		{
			//assuming faces are convex we will have at most one clipping edge per polygon
			//otherwise every pair of clipping nodes forming edge should appear consequently
			//as long as we go through face's edges in ordered way
			clipcoords.clear();
			//we will gather all the pairs of nodes, then form closed loop
			ElementArray<Face> faces = cells[k]->getFaces();
			Face full_face;
			int ntotpoints = 0, ntotedges = 0;
			for (INMOST_DATA_ENUM_TYPE q = 0; q < faces.size(); ++q)
			{
				int last_edge_type = CLIP_NONE;
				int nfulledges = 0, npoints = 0, nstartedge = ntotedges;
				ElementArray<Edge> edges = faces[q].getEdges();
				for (INMOST_DATA_ENUM_TYPE r = 0; r < edges.size(); ++r)
				{
					Storage::integer state = edges[r].IntegerDF(clip_state);
					if (state == CLIP_FULL)
					{
						nfulledges++;
						edge_point n1 = edge_point(&edges[r].getBeg()->Coords()[0], ntotedges, isColorBarEnabled() ? edges[r].getBeg()->RealDF(GetVisualizationTag()) : 0.f);
						edge_point n2 = edge_point(&edges[r].getEnd()->Coords()[0], ntotedges, isColorBarEnabled() ? edges[r].getEnd()->RealDF(GetVisualizationTag()) : 0.f);
						if (npoints % 2 == 0) //all privious edges are closed, just add this one
						{
							clipcoords.push_back(n1);
							clipcoords.push_back(n2);
							npoints += 2;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						}
						else if (n1 == clipcoords.back()) //this may be prolongation of one point that hit one edge
						{
							clipcoords.push_back(n2);
							npoints++;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						}
						else if (n2 == clipcoords.back())
						{
							clipcoords.push_back(n1);
							npoints++;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						}
						else printf("%s:%d strange orphan node before me\n", __FILE__, __LINE__);
					}
					else if (state == CLIP_ENDP)
					{
						edge_point n = edge_point(&edges[r].RealArrayDF(clip_point)[0], ntotedges, compute_value(edges[r], &edges[r].RealArrayDF(clip_point)[0]));
						bool add = true;
						if (last_edge_type == CLIP_ENDP)
						{
							if (n == clipcoords.back())
								add = false;
						}
						else if (last_edge_type == CLIP_FULL)
						{
							if (n == clipcoords.back() || n == clipcoords[clipcoords.size() - 2])
								add = false;
						}
						if (add) //this one node should be prolongation of privious edge
						{
							if (print)
							{
								printf("added: ");
								n.print();
							}
							clipcoords.push_back(n);
							npoints++;
							if (npoints % 2 == 0)
							{
								if (print) printf("edge %d accepted\n", ntotedges);
								ntotedges++;
							}
							last_edge_type = CLIP_ENDP;
						}
						else if (print)
						{
							printf("ignored: ");
							n.print();
						}
					}
					else if (state == CLIP_NODE)
					{
						edge_point n = edge_point(&edges[r].RealArrayDF(clip_point)[0], ntotedges, compute_value(edges[r], &edges[r].RealArrayDF(clip_point)[0]));
						if (print)
						{
							printf("added: ");
							n.print();
						}
						clipcoords.push_back(n);
						npoints++;
						if (npoints % 2 == 0)
						{
							if (print) printf("edge %d accepted\n", ntotedges);
							ntotedges++;
						}
						last_edge_type = CLIP_NODE;
					}
				}
				if (npoints % 2 != 0)
				{
					if (print) printf("edge %d not closed - remove\n", ntotedges);
					clipcoords.pop_back();
					npoints--;
					//printf("%s:%d this should not happen!\n",__FILE__,__LINE__);
				}

				if (nfulledges == static_cast<int>(edges.size()))
				{
					full_face = faces[q];
					break;
				}
				if (print)
				{
					printf("nodes on face %d\n", (int)faces[q].LocalID());
					for (int m = nstartedge * 2; m < static_cast<int>(clipcoords.size()); m++) clipcoords[m].print();
				}
				ntotpoints += npoints;
			}
			if (full_face.isValid())
			{
				ElementArray<Node> nodes = full_face->getNodes();
				Storage::real_array cl = cells[k]->RealArray(clips);
				Storage::real_array clv = cells[k]->RealArray(clipsv);
				Storage::integer_array cln = cells[k]->IntegerArray(clipsn);
				cln.resize(1, (int)nodes.size());
				cl.resize(static_cast<Storage::real_array::size_type>(3 * nodes.size()));
				clv.resize(static_cast<Storage::real_array::size_type>(nodes.size()));
				for (INMOST_DATA_ENUM_TYPE r = 0; r < nodes.size(); r++)
				{
					Storage::real_array p = nodes[r].Coords();
					cl[0 + 3 * r] = p[0];
					cl[1 + 3 * r] = p[1];
					cl[2 + 3 * r] = p[2];
					clv[r] = isColorBarEnabled() ? nodes[r].RealDF(GetVisualizationTag()) : 0.0;
				}
				cells[k]->SetMarker(marker);
			}
			else if (ntotedges > 2)
			{
				if (print)
				{
					printf("coords on cell %d\n", (int)cells[k]->LocalID());
					for (int m = 0; m < static_cast<int>(clipcoords.size()); m++) clipcoords[m].print();
				}
				//Can make this faster using hash
				closed.resize(ntotedges);
				std::fill(closed.begin(), closed.end(), false);
				Storage::real_array cl = cells[k]->RealArray(clips);
				Storage::real_array clv = cells[k]->RealArray(clipsv);
				Storage::integer_array cln = cells[k]->IntegerArray(clipsn);
				bool havestart = true;
				cln.clear();
				clv.clear();
				cl.clear();
				while (havestart)
				{
					havestart = false;
					for (int r = 0; r < ntotedges && !havestart; ++r) if (!closed[r])
					{
						loopcoords.push_back(clipcoords[r * 2 + 0]); //this is starting point
						loopcoords.push_back(clipcoords[r * 2 + 1]); //this is next
						closed[r] = true;
						havestart = true;
					}
					if (!havestart) break;
					bool hit = true;
					while (hit)
					{
						hit = false;
						for (int q = 0; q < ntotedges; ++q) if (!closed[q])
						{
							//some end of q-th edge connects to current end point - connect it
							if (clipcoords[q * 2 + 0] == loopcoords.back())
							{
								loopcoords.push_back(clipcoords[q * 2 + 1]);
								closed[q] = true;
								hit = true;
								break;
							}
							else if (clipcoords[q * 2 + 1] == loopcoords.back())
							{
								loopcoords.push_back(clipcoords[q * 2 + 0]);
								closed[q] = true;
								hit = true;
								break;
							}
						}
					}
					//loopcoords.pop_back();
					if (loopcoords.size() > 2)
					{
						cln.push_back((int)loopcoords.size());
						int offset = clv.size();
						cl.resize(static_cast<Storage::real_array::size_type>(offset * 3 + 3 * loopcoords.size()));
						clv.resize(static_cast<Storage::real_array::size_type>(offset + loopcoords.size()));
						for (INMOST_DATA_ENUM_TYPE r = 0; r < loopcoords.size(); ++r)
						{
							cl[(offset + r) * 3 + 0] = loopcoords[r].xyz[0];
							cl[(offset + r) * 3 + 1] = loopcoords[r].xyz[1];
							cl[(offset + r) * 3 + 2] = loopcoords[r].xyz[2];
							clv[offset + r] = loopcoords[r].val;
						}
					}
					loopcoords.clear();
				}
				clipcoords.clear();

				cells[k]->SetMarker(marker);
			}
		}
	}

	void clipper::gen_clip(std::vector<face2gl> & out, Storage::real n[3], bool elevation)
	{
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1, std::min<INMOST_DATA_ENUM_TYPE>(15, size() / 100));
		for (INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k++) if (cells[k]->GetMarker(marker))
		{
			Storage::real_array cl = cells[k]->RealArray(clips);
			Storage::real_array clv = cells[k]->RealArray(clipsv);
			Storage::integer_array cln = cells[k]->IntegerArray(clipsn);
			INMOST_DATA_ENUM_TYPE offset = 0;
			for (int r = 0; r < (int)cln.size(); ++r)
			{
				face2gl f;
				f.set_color(0.6, 0.6, 0.6, 1);
				if (elevation && isColorBarEnabled())
				{
					Storage::real pos[3], t;
					for (INMOST_DATA_ENUM_TYPE q = offset; q < (INMOST_DATA_ENUM_TYPE)offset + cln[r]; q++)
					{
						t = (clv[q] - GetColorBar()->get_min()) / (GetColorBar()->get_max() - GetColorBar()->get_min());
						pos[0] = cl[q * 3 + 0] + t*n[0];
						pos[1] = cl[q * 3 + 1] + t*n[1];
						pos[2] = cl[q * 3 + 2] + t*n[2];
						f.add_vert(pos);
					}
				}
				else
				for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++) f.add_vert(&cl[q * 3]);
				if (isColorBarEnabled())
				{
					for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++)
					{
						//f.add_color(GetColorBar()->pick_color(clv[q]));
						if (GetVisualizationType() == CELL && !isVisualizationSmooth())
							f.add_texcoord(GetColorBar()->pick_texture(cells[k].RealDF(GetVisualizationTag())));
						else
							f.add_texcoord(GetColorBar()->pick_texture(clv[q]));
					}
				}
				f.compute_center();
				f.set_elem(cells[k]->GetElementType(), cells[k]->LocalID());
				if (k%pace == 0) f.set_flag(true);
				out.push_back(f);
				offset += cln[r];
			}
		}
	}


	void clipper::draw_clip(INMOST_DATA_ENUM_TYPE pace, Storage::real n[3], bool elevation)
	{
		if (isColorBarEnabled()) GetColorBar()->BindTexture();
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glBegin(GL_TRIANGLES);
		for (INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k += pace) if (cells[k]->GetMarker(marker))
		{
			Storage::real_array cl = cells[k]->RealArray(clips);
			Storage::real_array clv = cells[k]->RealArray(clipsv);
			Storage::integer_array cln = cells[k]->IntegerArray(clipsn);
			INMOST_DATA_ENUM_TYPE offset = 0;
			for (int r = 0; r < (int)cln.size(); ++r)
			{
				Storage::real cnt[3] = { 0, 0, 0 };
				Storage::real cntv = 0.0;
				for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++)
				{
					cnt[0] += cl[q * 3 + 0];
					cnt[1] += cl[q * 3 + 1];
					cnt[2] += cl[q * 3 + 2];
					cntv += clv[q];
				}
				cnt[0] /= static_cast<Storage::real>(cln[r]);
				cnt[1] /= static_cast<Storage::real>(cln[r]);
				cnt[2] /= static_cast<Storage::real>(cln[r]);
				cntv /= static_cast<Storage::real>(cln[r]);
				if (!isColorBarEnabled())
				{
					for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++)
					{
						glVertexNdv(cnt);
						glVertexNdv(&cl[q * 3]);
						glVertexNdv(&cl[(((q + 1 - offset) % cln[r]) + offset) * 3]);
					}
				}
				else if (elevation)
				{
					Storage::real pos[3], t;
					for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++)
					{
						if (GetVisualizationType() == CELL && !isVisualizationSmooth())
							glTexCoord1d(GetColorBar()->pick_texture(cells[k].RealDF(GetVisualizationTag())));
						else
							glTexCoord1d(GetColorBar()->pick_texture(cntv));
						t = (cntv - GetColorBar()->get_min()) / (GetColorBar()->get_max() - GetColorBar()->get_min());
						pos[0] = cnt[0] + n[0] * t;
						pos[1] = cnt[1] + n[1] * t;
						pos[2] = cnt[2] + n[2] * t;
						glVertexNdv(pos);
						if (!(GetVisualizationType() == CELL && !isVisualizationSmooth()))
							glTexCoord1d(GetColorBar()->pick_texture(clv[q]));
						t = (clv[q] - GetColorBar()->get_min()) / (GetColorBar()->get_max() - GetColorBar()->get_min());
						pos[0] = cl[q * 3 + 0] + n[0] * t;
						pos[1] = cl[q * 3 + 1] + n[1] * t;
						pos[2] = cl[q * 3 + 2] + n[2] * t;
						glVertexNdv(pos);
						if (!(GetVisualizationType() == CELL && !isVisualizationSmooth()))
							glTexCoord1d(GetColorBar()->pick_texture(clv[(q + 1 - offset) % cln[r] + offset]));
						t = (clv[(q + 1 - offset) % cln[r] + offset] - GetColorBar()->get_min()) / (GetColorBar()->get_max() - GetColorBar()->get_min());
						pos[0] = cl[((q + 1 - offset) % cln[r] + offset) * 3 + 0] + n[0] * t;
						pos[1] = cl[((q + 1 - offset) % cln[r] + offset) * 3 + 1] + n[1] * t;
						pos[2] = cl[((q + 1 - offset) % cln[r] + offset) * 3 + 2] + n[2] * t;
						glVertexNdv(pos);
					}
				}
				else
				{
					for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++)
					{
						if (GetVisualizationType() == CELL && !isVisualizationSmooth())
							glTexCoord1d(GetColorBar()->pick_texture(cells[k].RealDF(GetVisualizationTag())));
						else
							glTexCoord1d(GetColorBar()->pick_texture(cntv));
						glVertexNdv(cnt);
						if (!(GetVisualizationType() == CELL && !isVisualizationSmooth()))
							glTexCoord1d(GetColorBar()->pick_texture(clv[q]));
						glVertexNdv(&cl[q * 3]);
						if (!(GetVisualizationType() == CELL && !isVisualizationSmooth()))
							glTexCoord1d(GetColorBar()->pick_texture(clv[(q + 1 - offset) % cln[r] + offset]));
						glVertexNdv(&cl[((q + 1 - offset) % cln[r] + offset) * 3]);
					}
				}
				offset += cln[r];
			}
		}
		glEnd();

		glDisable(GL_POLYGON_OFFSET_FILL);
		if (isColorBarEnabled()) GetColorBar()->UnbindTexture();
	}


	void clipper::draw_clip_edges(INMOST_DATA_ENUM_TYPE pace, Storage::real n[3], bool elevation)
	{
		
		glBegin(GL_LINES);
		if (isColorBarEnabled() && elevation)
		{
			for (INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k += pace) if (cells[k]->GetMarker(marker))
			{
				Storage::real_array cl = cells[k]->RealArray(clips);
				Storage::real_array clv = cells[k]->RealArray(clipsv);
				Storage::integer_array cln = cells[k]->IntegerArray(clipsn);
				INMOST_DATA_ENUM_TYPE offset = 0;
				for (size_t r = 0; r < cln.size(); ++r)
				{
					Storage::real pos[3], t;
					for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++)
					{
						t = (clv[q] - GetColorBar()->get_min()) / (GetColorBar()->get_max() - GetColorBar()->get_min());
						pos[0] = cl[q * 3 + 0] + t*n[0];
						pos[1] = cl[q * 3 + 1] + t*n[1];
						pos[2] = cl[q * 3 + 2] + t*n[2];
						glVertexNdv(pos);
						t = (clv[(q + 1 - offset) % cln[r] + offset] - GetColorBar()->get_min()) / (GetColorBar()->get_max() - GetColorBar()->get_min());
						pos[0] = cl[((q + 1 - offset) % cln[r] + offset) * 3 + 0] + t*n[0];
						pos[1] = cl[((q + 1 - offset) % cln[r] + offset) * 3 + 1] + t*n[1];
						pos[2] = cl[((q + 1 - offset) % cln[r] + offset) * 3 + 2] + t*n[2];
						glVertexNdv(pos);
					}
					offset += cln[r];
				}
			}
		}
		else
		{
			for (INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k += pace) if (cells[k]->GetMarker(marker))
			{
				Storage::real_array cl = cells[k]->RealArray(clips);
				Storage::integer_array cln = cells[k]->IntegerArray(clipsn);
				INMOST_DATA_ENUM_TYPE offset = 0;
				for (size_t r = 0; r < cln.size(); ++r)
				{
					for (INMOST_DATA_ENUM_TYPE q = offset; q < offset + cln[r]; q++)
					{
						glVertexNdv(&cl[q * 3]);
						glVertexNdv(&cl[((q + 1 - offset) % cln[r] + offset) * 3]);
					}
					offset += cln[r];
				}
			}
		}
		glEnd();

	}

	INMOST_DATA_ENUM_TYPE clipper::size() { return (INMOST_DATA_ENUM_TYPE)cells.size(); }


	bnd_clipper::~bnd_clipper()
	{
		mm->DeleteTag(clip_state);
		if( tree ) delete tree;
		delete[]faces;
	}

	bnd_clipper::bnd_clipper(Mesh * m, HandleType * _faces, INMOST_DATA_ENUM_TYPE size)
	{
		mm = m;
		clip_state = m->CreateTag("CLIP_FACE_STATE", DATA_INTEGER, (mm->GetDimensions() == 2 ? CELL : NONE) |FACE, NONE, 1);
		nfaces = size;
		faces = new HandleType[nfaces];
		for (INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++) faces[k] = _faces[k];
		if( mm->GetDimensions() == 2 )
		{
			for (INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++)
				mm->IntegerDF(faces[k],clip_state) = CLIP_FACE_INSIDE;
			tree = NULL;
		}
		else tree = new kdtree(mm, faces, nfaces);
	}

	Storage::real bnd_clipper::compute_value(Node n1, Node n2, Storage::real * c1, Storage::real * c2, Storage::real * pnt)
	{
		if (isColorBarEnabled())
		{
			Storage::real d1, d2, t;
			d1 = sqrt((pnt[0] - c1[0])*(pnt[0] - c1[0]) + (pnt[1] - c1[1])*(pnt[1] - c1[1]) + (pnt[2] - c1[2])*(pnt[2] - c1[2]));
			d2 = sqrt((c2[0] - c1[0])*(c2[0] - c1[0]) + (c2[1] - c1[1])*(c2[1] - c1[1]) + (c2[2] - c1[2])*(c2[2] - c1[2]));
			t = d1 / d2; //(pnt == c2, t = 1 : pnt == c1, t = 0)
			return n1->RealDF(GetVisualizationTag())*(1 - t) + n2->RealDF(GetVisualizationTag())*t;
		}
		else return 0.f;
	}

	void bnd_clipper::clip_plane(Storage::real p[3], Storage::real n[3])
	{
		if( mm->GetDimensions() != 2 )
			tree->intersect_plane_face(clip_state, p, n);
	}

	void bnd_clipper::gen_clip(std::vector<face2gl> & out, Storage::real p[3], Storage::real n[3], bool elevation)
	{
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1, std::min<INMOST_DATA_ENUM_TYPE>(15, nfaces / 100));
		for (INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++)
		{
			int state = mm->IntegerDF(faces[k], clip_state);
			if (state == CLIP_FACE_INSIDE)
			{
				ElementArray<Node> nodes = Element(mm, faces[k])->getNodes();
				face2gl f;
				f.set_color(0.6, 0.6, 0.6, 1);
				for (INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
				{
					f.add_vert(&nodes[q].Coords()[0], mm->GetDimensions());
					if (isColorBarEnabled())
					{
						//f.add_color(GetColorBar()->pick_color(nodes[q].RealDF(GetVisualizationTag())));
						if (GetVisualizationType() == CELL && !isVisualizationSmooth())
						{
							if( GetHandleElementType(faces[k]) == CELL )
								f.add_texcoord(GetColorBar()->pick_texture(Cell(mm, faces[k]).RealDF(GetVisualizationTag())));
							else
								f.add_texcoord(GetColorBar()->pick_texture(Face(mm, faces[k]).BackCell().RealDF(GetVisualizationTag())));
						}
						else
							f.add_texcoord(GetColorBar()->pick_texture(nodes[q].RealDF(GetVisualizationTag())));
					}
				}
				f.compute_center();
				f.set_elem(GetHandleElementType(faces[k]), GetHandleID(faces[k]));
				if (k%pace == 0) f.set_flag(true);
				out.push_back(f);
			}
			else if (state == CLIP_FACE_INTERSECT)
			{
				face2gl f;
				f.set_color(0.6, 0.6, 0.6, 1);
				ElementArray<Node> nodes = Element(mm, faces[k])->getNodes();
				std::vector<bool> nodepos(nodes.size());
				std::vector<Storage::real> faceverts;
				Storage::real_array coords = nodes[0].Coords();
				for (INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
				{
					coords = nodes[q].Coords();
					Storage::real dot = n[0] * (coords[0] - p[0]) + n[1] * (coords[1] - p[1]) + n[2] * (coords[2] - p[2]);
					nodepos[q] = dot < 1.0e-10;
				}
				for (INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
				{
					if (nodepos[q])
					{
						coords = nodes[q].Coords();
						f.add_vert(&coords[0]);
						if (isColorBarEnabled())
						{
							//f.add_color(GetColorBar()->pick_color(nodes[q].RealDF(GetVisualizationTag())));
							if (GetVisualizationType() == CELL && !isVisualizationSmooth())
								f.add_texcoord(GetColorBar()->pick_texture(Face(mm, faces[k]).BackCell().RealDF(GetVisualizationTag())));
							else
								f.add_texcoord(GetColorBar()->pick_texture(nodes[q].RealDF(GetVisualizationTag())));
						}
					}
					if (nodepos[q] != nodepos[(q + 1) % nodes.size()])
					{
						Storage::real_array sp0 = nodes[q].Coords();
						Storage::real_array sp1 = nodes[(q + 1) % nodes.size()].Coords();
						Storage::real node[3], t, c = 0;
						if (clip_plane_edge(&sp0[0], &sp1[0], p, n, node) > CLIP_NONE)
						{
							if (isColorBarEnabled())
							{
								//f.add_color(GetColorBar()->pick_color(compute_value(nodes[q],nodes[(q+1)%nodes.size()],&sp0[0],&sp1[0],node)));
								if (GetVisualizationType() == CELL && !isVisualizationSmooth())
									c = Face(mm, faces[k]).BackCell().RealDF(GetVisualizationTag());
								else
									c = compute_value(nodes[q], nodes[(q + 1) % nodes.size()], &sp0[0], &sp1[0], node);
								//f.add_texcoord(GetColorBar()->pick_texture(Face(mm,faces[k]).BackCell().RealDF(GetVisualizationTag())));
								f.add_texcoord(GetColorBar()->pick_texture(c));
								if (elevation)
								{
									c = compute_value(nodes[q], nodes[(q + 1) % nodes.size()], &sp0[0], &sp1[0], node);
									t = (c - GetColorBar()->get_min()) / (GetColorBar()->get_max() - GetColorBar()->get_min());
									node[0] += t*n[0];
									node[1] += t*n[1];
									node[2] += t*n[2];
								}
							}
							f.add_vert(node);
						}
					}
				}
				f.compute_center();
				f.set_elem(GetHandleElementType(faces[k]), GetHandleID(faces[k]));
				if (k%pace == 0) f.set_flag(true);
				out.push_back(f);
			}
		}
	}

	void bnd_clipper::draw_clip(INMOST_DATA_ENUM_TYPE pace, Storage::real p[3], Storage::real n[3])
	{
		for (INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k += pace)
		{
			int state = CLIP_FACE_NONE;
			ElementArray<Node> nodes = Face(mm,faces[k])->getNodes();
			Storage::real_array coords = nodes[0].Coords();
			Storage::real dot0 = n[0] * (coords[0] - p[0]) + n[1] * (coords[1] - p[1]) + n[2] * (coords[2] - p[2]);
			if (dot0 <= 0.0) state = CLIP_FACE_INSIDE; else state = CLIP_FACE_OUTSIDE;
			if (state == CLIP_FACE_INSIDE)
			{
				for (INMOST_DATA_ENUM_TYPE q = 1; q < nodes.size(); ++q)
				{
					coords = nodes[q].Coords();
					Storage::real dot = n[0] * (coords[0] - p[0]) + n[1] * (coords[1] - p[1]) + n[2] * (coords[2] - p[2]);
					if (dot*dot0 <= 0.0)
					{
						state = CLIP_FACE_INTERSECT;
						break;
					}
				}
				if (state == CLIP_FACE_INSIDE)
				{
					ElementArray<Node> nodes = Face(mm, faces[k])->getNodes();

					glEnable(GL_POLYGON_OFFSET_FILL);
					glPolygonOffset(1.0, 1.0);


					glBegin(GL_POLYGON);
					if (isColorBarEnabled())
					{
						for (INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
						{
							GetColorBar()->pick_color(nodes[q].RealDF(GetVisualizationTag())).set_color();
							glTexCoord1f(GetColorBar()->pick_texture(nodes[q].RealDF(GetVisualizationTag())));
							glVertexNdv(&nodes[q].Coords()[0],mm->GetDimensions());
						}
					}
					else
					{
						for (INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
							glVertexNdv(&nodes[q].Coords()[0],mm->GetDimensions());
					}
					glEnd();

					glDisable(GL_POLYGON_OFFSET_FILL);
				}
			}
			mm->IntegerDF(faces[k], clip_state) = state;
		}
	}

	void bnd_clipper::draw_clip_edges(INMOST_DATA_ENUM_TYPE pace)
	{
		for (INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k += pace)
		{
			if (mm->IntegerDF(faces[k], clip_state) == CLIP_FACE_INSIDE)
			{
				ElementArray<Node> nodes = Element(mm, faces[k])->getNodes();
				glBegin(GL_LINE_LOOP);
				for (INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++) glVertexNdv(&nodes[q].Coords()[0],mm->GetDimensions());
				glEnd();
			}
		}
	}

	INMOST_DATA_ENUM_TYPE bnd_clipper::size() { return nfaces; }
}
