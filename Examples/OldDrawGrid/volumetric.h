#ifndef _VOLUMETRIC_H
#define _VOLUMETRIC_H

#include "inmost.h"
#include "face2gl.h"

namespace INMOST
{

	typedef struct point
	{
		float coords[3];
		float diam;
		float dist;
		int id;
		point() {}
		point & operator = (point const & b);
		point(const point & b);
		bool operator <(const point & b) { return dist < b.dist; }
	} point_t;


	class volumetric
	{
		std::vector<point_t> points;
		Mesh * m;

		void radix_sort_dist(std::vector<point_t> & set);
	public:

		volumetric(Mesh * _m);
		void camera(double pos[3], int interactive);
		void draw(int interactive) const;
	};

	class volumetric2
	{
		std::vector<face2gl> faces;
		Mesh * m;
	public:
		volumetric2(Mesh * _m);
		void camera(double pos[3], int interactive);
		void draw(int interactive);
	};
}

#endif