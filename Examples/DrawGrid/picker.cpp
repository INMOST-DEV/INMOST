#include "picker.h"

namespace INMOST
{
	static int cmpElements0(const void * a, const void * b)
	{
		const kdtree_picker::entry * ea = ((const kdtree_picker::entry *)a);
		const kdtree_picker::entry * eb = ((const kdtree_picker::entry *)b);
		float ad = ea->xyz[0];
		float bd = eb->xyz[0];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements1(const void * a, const void * b)
	{
		const kdtree_picker::entry * ea = ((const kdtree_picker::entry *)a);
		const kdtree_picker::entry * eb = ((const kdtree_picker::entry *)b);
		float ad = ea->xyz[1];
		float bd = eb->xyz[1];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements2(const void * a, const void * b)
	{
		const kdtree_picker::entry * ea = ((const kdtree_picker::entry *)a);
		const kdtree_picker::entry * eb = ((const kdtree_picker::entry *)b);
		float ad = ea->xyz[2];
		float bd = eb->xyz[2];
		return (ad > bd) - (ad < bd);
	}

	inline static unsigned int flip(const unsigned int * fp) { unsigned int mask = -((int)(*fp >> 31)) | 0x80000000; return *fp ^ mask; }
	inline static unsigned int _0(unsigned int x)	{ return x & 0x7FF; }
	inline static unsigned int _1(unsigned int x)	{ return x >> 11 & 0x7FF; }
	inline static unsigned int _2(unsigned int x)   { return x >> 22; }



	void kdtree_picker::radix_sort(int dim, struct entry * temp)
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

	void kdtree_picker::kdtree_build(int dim, std::vector<face2gl> & in, struct entry * temp)
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
			children = static_cast<kdtree_picker *>(malloc(sizeof(kdtree_picker)* 2));//new kdtree_picker[2];
			assert(children != NULL);
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size / 2;
			children[1].children = NULL;
			children[1].set = set + size / 2;
			children[1].size = size - size / 2;
			children[0].kdtree_build((dim + 1) % 3, in, temp);
			children[1].kdtree_build((dim + 1) % 3, in, temp);
			for (int k = 0; k < 3; k++)
			{
				bbox[0 + 2 * k] = std::min(children[0].bbox[0 + 2 * k], children[1].bbox[0 + 2 * k]);
				bbox[1 + 2 * k] = std::max(children[0].bbox[1 + 2 * k], children[1].bbox[1 + 2 * k]);
			}
		}
		else
		{
			assert(size == 1);
			bbox[0] = bbox[2] = bbox[4] = 1.0e20;
			bbox[1] = bbox[3] = bbox[5] = -1.0e20;
			for (unsigned k = 0; k < in[set[0].index].size(); ++k)
			{
				double * coords = in[set[0].index].get_vert(k);
				for (unsigned q = 0; q < 3; q++)
				{
					bbox[q * 2 + 0] = std::min(bbox[q * 2 + 0], coords[q]);
					bbox[q * 2 + 1] = std::max(bbox[q * 2 + 1], coords[q]);
				}
			}
		}
	}

	void kdtree_picker::clear_children()
	{
		if (children)
		{
			children[0].clear_children();
			children[1].clear_children();
			free(children);
		}
	}

	int kdtree_picker::raybox(double pos[3], double ray[3], double closest) const
	{
		double tnear = -1.0e20, tfar = 1.0e20, t1, t2, c;
		for (int i = 0; i < 3; i++)
		{
			if (fabs(ray[i]) < 1.0e-15)
			{
				if (pos[i] < bbox[i * 2] || pos[i] > bbox[i * 2 + 1])
					return 0;
			}
			else
			{
				t1 = (bbox[i * 2 + 0] - pos[i]) / ray[i];
				t2 = (bbox[i * 2 + 1] - pos[i]) / ray[i];
				if (t1 > t2) { c = t1; t1 = t2; t2 = c; }
				if (t1 > tnear) tnear = t1;
				if (t2 < tfar) tfar = t2;
				if (tnear > closest) return 0;
				if (tnear > tfar) return 0;
				if (tfar < 0) return 0;
			}
		}
		return 1;
	}

	void kdtree_picker::sub_intersect_ray_faces(std::vector<face2gl> & in, double p[3], double dir[3], std::pair<double, int> & closest) const
	{
		if (size == 1)
		{
			face2gl & f = in[set[0].index];
			double * tri[3], btri[3][3], dot[3], prod[3][3], norm[3], d, proj[3];
			tri[0] = f.get_center();
			for (unsigned i = 0; i < f.size(); i++)
			{
				unsigned j = (i + 1) % f.size();
				tri[1] = f.get_vert(i);
				tri[2] = f.get_vert(j);
				norm[0] = (tri[2][1] - tri[0][1])*(tri[1][2] - tri[0][2]) - (tri[2][2] - tri[0][2])*(tri[1][1] - tri[0][1]);
				norm[1] = (tri[2][2] - tri[0][2])*(tri[1][0] - tri[0][0]) - (tri[2][0] - tri[0][0])*(tri[1][2] - tri[0][2]);
				norm[2] = (tri[2][0] - tri[0][0])*(tri[1][1] - tri[0][1]) - (tri[2][1] - tri[0][1])*(tri[1][0] - tri[0][0]);
				d = norm[0] * (tri[0][0] - p[0]) + norm[1] * (tri[0][1] - p[1]) + norm[2] * (tri[0][2] - p[2]);
				d /= norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2];
				proj[0] = p[0] + d*dir[0];
				proj[1] = p[1] + d*dir[1];
				proj[2] = p[2] + d*dir[2];

				if (d > closest.first) break;

				for (int k = 0; k < 3; k++)
				{
					btri[k][0] = tri[k][0] - proj[0];
					btri[k][1] = tri[k][1] - proj[1];
					btri[k][2] = tri[k][2] - proj[2];
				}
				for (int k = 0; k < 3; k++)
				{
					int l = (k + 1) % 3;
					prod[k][0] = btri[k][1] * btri[l][2] - btri[k][2] * btri[l][1];
					prod[k][1] = btri[k][2] * btri[l][0] - btri[k][0] * btri[l][2];
					prod[k][2] = btri[k][0] * btri[l][1] - btri[k][1] * btri[l][0];
				}

				for (int k = 0; k < 3; k++)
				{
					int l = (k + 1) % 3;
					dot[k] = prod[k][0] * prod[l][0] + prod[k][1] * prod[l][1] + prod[k][2] * prod[l][2];
				}
				if (dot[0] >= 0 && dot[1] >= 0 && dot[2] >= 0)
				{
					closest.first = d;
					closest.second = set[0].index;
					break; //don't expect anything better here
				}
			}
		}
		else
		{
			if (raybox(p, dir, closest.first))
			{
				children[0].sub_intersect_ray_faces(in, p, dir, closest);
				children[1].sub_intersect_ray_faces(in, p, dir, closest);
			}
		}
	}

	kdtree_picker::kdtree_picker() : set(NULL), size(0), children(NULL) {}

	kdtree_picker::~kdtree_picker() { clear_children(); delete[] set; }

	kdtree_picker::kdtree_picker(std::vector<face2gl> & in)
	{
		size = (unsigned)in.size();
		set = new struct entry[size];
		for (unsigned k = 0; k < size; ++k)
		{
			set[k].index = k;
			in[k].get_center(set[k].xyz);
		}
		struct entry * temp = new struct entry[size];
		kdtree_build(0, in, temp);
		delete[] temp;
	}

	int kdtree_picker::ray_faces(std::vector<face2gl> & in, double p[3], double dir[3]) const
	{
		std::pair<double, int> closest(1.0e20, -1);
		sub_intersect_ray_faces(in, p, dir, closest);
		return closest.second;
	}

	picker::~picker() { delete tree; }

	picker::picker(std::vector<face2gl> & _faces)
	{
		faces = &_faces;
		tree = new kdtree_picker(*faces);
	}

	int picker::select(double p[3], double ray[3]) const
	{
		return tree->ray_faces(*faces, p, ray);
	}
}
