#ifndef _FACE2GL_H
#define _FACE2GL_H

#include "inmost.h"
#include "inc_glut.h"
#include "color.h"

namespace INMOST
{

	class face2gl
	{
		float dist;
		double c[4];
		double cnt[3];
		bool flag;
		ElementType etype;
		Storage::integer id;
		std::vector<double> verts;
		std::vector<double> texcoords;
		std::vector<color_t> colors;
		color_t cntcolor;
		double cnttexcoord;
	public:
		void shift(double x, double y, double z);
		static void radix_sort_dist(std::vector<face2gl> & set);
		face2gl();
		face2gl(const face2gl & other);
		face2gl & operator =(face2gl const & other);
		~face2gl();
		void draw_colour() const;
		void svg_draw_colour(std::ostream & file, bool drawedges, double modelview[16], double projection[16], int viewport[4]) const;
		void draw_colour_alpha(double alpha) const;
		void draw() const;
		void svg_draw(std::ostream & file, bool drawedges, double modelview[16], double projection[16], int viewport[4]) const;
		void wgl_draw_bin(std::ostream & file) const;
		void wgl_draw_edges_bin(std::ostream & file) const;
		void wgl_draw(std::ostream & file) const;
		void wgl_draw_edges(std::ostream & file) const;
		void drawedges() const;
		void svg_drawedges(std::ostream & file, double modelview[16], double projection[16], int viewport[4]) const;
		bool operator <(const face2gl & other) const { return dist < other.dist; }
		void set_color(double r, double g, double b, double a) { c[0] = r; c[1] = g; c[2] = b; c[3] = a; }
		void add_vert(double x, double y, double z) { unsigned s = (unsigned)verts.size(); verts.resize(s + 3); verts[s] = x; verts[s + 1] = y; verts[s + 2] = z; }
		void add_vert(double v[3]) { verts.insert(verts.end(), v, v + 3); }
		void add_vert(float v[3]) { verts.insert(verts.end(), v, v + 3); }
		void add_vert(double * v, int N) { double vv[3] = {0,0,0}; for(int k = 0; k < N; ++k) vv[k] = v[k]; verts.insert(verts.end(), vv, vv + 3); }
		void add_vert(float * v, int N) { double vv[3] = {0,0,0}; for(int k = 0; k < N; ++k) vv[k] = v[k]; verts.insert(verts.end(), vv, vv + 3); }
		void add_color(color_t c) { colors.push_back(c); }
		void add_texcoord(double val) { texcoords.push_back(val); }
		double * get_vert(int k) { return &verts[k * 3]; }
		unsigned size() { return (unsigned)verts.size() / 3; }
		void set_center(double _cnt[3], color_t c = color_t(0, 0, 0, 0)) {cnt[0] = _cnt[0]; cnt[1] = _cnt[1]; cnt[2] = _cnt[2]; cntcolor = c;}
		void get_center(float _cnt[3]) {_cnt[0] = cnt[0]; _cnt[1] = cnt[1];	_cnt[2] = cnt[2];}
		void get_center(double _cnt[3]) {_cnt[0] = cnt[0]; _cnt[1] = cnt[1]; _cnt[2] = cnt[2];}
		double * get_center() { return cnt; }
		void compute_center();
		void compute_center_color();
		void compute_center_texcoord();
		void compute_dist(double cam[3]) {dist = sqrt((cnt[0] - cam[0])*(cnt[0] - cam[0]) + (cnt[1] - cam[1])*(cnt[1] - cam[1]) + (cnt[2] - cam[2])*(cnt[2] - cam[2]));}
		void set_flag(bool set) { flag = set; }
		bool get_flag() { return flag; }
		void set_elem(ElementType _etype, Storage::integer _id) { etype = _etype; id = _id; }
		Element get_elem(Mesh * m) { return m->ElementByLocalID(etype, id); }
	};

	face2gl DrawFace(Element f);

	void svg_draw_faces_nc(std::ostream & file, std::vector<face2gl> & set, bool drawedges, double modelview[16], double projection[16], int viewport[4], int highlight = -1);
	void svg_draw_faces(std::ostream & file, std::vector<face2gl> & set, bool drawedges, double modelview[16], double projection[16], int viewport[4], int highlight = -1);
	void svg_draw_edges(std::ostream & file, std::vector<face2gl> & set, double modelview[16], double projection[16], int viewport[4], int highlight = -1);
	void wgl_draw_faces_bin(std::ostream & file, const std::vector<face2gl> & set);
	void wgl_draw_edges_bin(std::ostream & file, const std::vector<face2gl> & set);
	void wgl_draw_faces(std::ostream & file, const std::vector<face2gl> & set);
	void wgl_draw_edges(std::ostream & file, const std::vector<face2gl> & set);

	void draw_faces_nc(std::vector<face2gl> & set, int highlight = -1);
	void draw_faces(std::vector<face2gl> & set, int highlight = -1);
	void draw_faces_alpha(std::vector<face2gl> & set, double alpha);
	void draw_edges(std::vector<face2gl> & set, int highlight = -1);
	void draw_faces_interactive_nc(std::vector<face2gl> & set);
	void draw_faces_interactive(std::vector<face2gl> & set);
	void draw_faces_interactive_alpha(std::vector<face2gl> & set, double alpha);
	void draw_edges_interactive(std::vector<face2gl> & set);
}
#endif
