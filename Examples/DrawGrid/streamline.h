#ifndef _STREAMLINE_H
#define _STREAMLINE_H

#include "inmost.h"
#include "coord.h"
#include "octree.h"


namespace INMOST
{
	

	class Streamline
	{
	private:
		std::vector<coord> points;
		std::vector<double> velarr;
	public:
		Streamline() {}
		Streamline(const Octree & octsearch, coord pos, Tag velocity_tag, ElementType velocity_defined, Tag cell_size, Storage::real velocity_min, Storage::real velocity_max, Storage::real sign, MarkerType visited);
		Streamline(const Streamline & other) { points = other.points; velarr = other.velarr; }
		Streamline & operator =(Streamline const & other) { points = other.points; velarr = other.velarr; return *this; }
		~Streamline() { points.clear(); velarr.clear(); }
		void Draw(int reduced);
		void SVGDraw(std::ostream & file, double modelview[16], double projection[16], int viewport[4]);
	};

	void BuildStreamlines(Mesh *m, Tag vel, ElementType vel_def, std::vector<Streamline> & output);
}
#endif
