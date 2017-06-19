#ifndef _PICKER_H
#define _PICKER_H

#include "face2gl.h"

namespace INMOST
{
	class kdtree_picker
	{
	public:
		struct entry
		{
			int index;
			float xyz[3];
		} *set;
	private:
		unsigned size;
		double bbox[6];
		kdtree_picker * children;
		
		void radix_sort(int dim, struct entry * temp);
		void kdtree_build(int dim, std::vector<face2gl> & in, struct entry * temp);
		void clear_children();
		int raybox(double pos[3], double ray[3], double closest) const;
		void sub_intersect_ray_faces(std::vector<face2gl> & in, double p[3], double dir[3], std::pair<double, int> & closest) const;
	public:
		kdtree_picker();
		~kdtree_picker();
		kdtree_picker(std::vector<face2gl> & in);
		int ray_faces(std::vector<face2gl> & in, double p[3], double dir[3]) const;
	};

	class picker
	{
		kdtree_picker * tree;
		std::vector<face2gl> * faces;
	public:
		~picker();
		picker(std::vector<face2gl> & _faces);
		int select(double p[3], double ray[3]) const;
	};
}
#endif