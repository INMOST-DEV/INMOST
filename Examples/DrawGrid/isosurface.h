#ifndef _ISOSURF_H
#define _ISOSURF_H

#include "inmost.h"
#include "coord.h"
#include "octree.h"


namespace INMOST
{
	
	struct triangle
	{
		int a, b, c;
	};

	class Isosurface
	{
	private:
		std::vector<coord> points;
		std::vector<triangle> tris;
		std::vector<double> vals;
	public:
		Isosurface() {}
		Isosurface();
		Isosurface(const Isosurface & other) { points = other.points; tris = other.tris; vals = other.vals; }
		Isosurface & operator =(Streamline const & other) { points = other.points; tris = other.tris; vals = other.vals; return *this; }
		~Isosurface() { points.clear(); tris.clear(); vals.clear() }
		void Draw(int reduced);
		void SVGDraw(std::ostream & file, double modelview[16], double projection[16], int viewport[4]);
	};

	void BuildIsosurfaces(Mesh *m, Tag values, double iso, std::vector<Isosurface> & output);
}
#endif
