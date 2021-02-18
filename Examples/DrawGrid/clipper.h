#ifndef _CLIPPER_H
#define _CLIPPER_H

#include "inmost.h"
#include "face2gl.h"

namespace INMOST
{

	class kdtree
	{
	public:
		struct entry
		{
			HandleType e;
			float xyz[3];
			entry(){}
			entry(const entry & other);
			entry & operator =(const entry & other);
		} *set;
	private:
		int marked;
		Mesh * m;
		INMOST_DATA_ENUM_TYPE size;
		float bbox[6];
		kdtree * children;
		void radix_sort(int dim, struct entry * temp);
		void kdtree_build(int dim, int & done, int total, struct entry * temp);
		kdtree();
		inline int plane_bbox(Storage::real p[3], Storage::real n[3]) const;
		bool sub_intersect_plane_edge(Tag clip_point, Tag clip_state, ElementArray<Cell> & cells, MarkerType mrk, Storage::real p[3], Storage::real n[3]);
		void sub_intersect_plane_faces(Tag clip_state, Storage::real p[3], Storage::real n[3]);
		void unmark_old_edges(Tag clip_state);
		void clear_children();
	public:
		kdtree(Mesh * m);
		kdtree(Mesh * m, HandleType * eset, INMOST_DATA_ENUM_TYPE size);
		void intersect_plane_edge(Tag clip_point, Tag clip_state, ElementArray<Cell> & cells, MarkerType mark_cells, Storage::real p[3], Storage::real n[3]);
		void intersect_plane_face(Tag clip_state, Storage::real p[3], Storage::real n[3]);
		~kdtree();
	};

	class clipper
	{
	public:
		struct edge_point
		{
			Storage::real val;
			Storage::real xyz[3];
			Storage::integer edge;
			edge_point();
			edge_point(Storage::real _xyz[3], Storage::integer n, float v);
			bool operator ==(const edge_point& b) const;
			bool operator !=(const edge_point& b) const;
			void print();
		};
	private:
		Tag clip_point, clip_state;
		kdtree * tree;
		Tag clips, clipsv, clipsn;
		MarkerType marker;
		ElementArray<Cell> cells;
		Mesh * mm;
	public:
		~clipper();
		clipper(Mesh * m);
		Storage::real compute_value(Edge e, Storage::real * pnt);
		void clip_plane(Storage::real p[3], Storage::real n[3]);
		void gen_clip(std::vector<face2gl> & out, Storage::real n[3],bool elevation);
		void draw_clip(INMOST_DATA_ENUM_TYPE pace, Storage::real n[3], bool elevation);
		void draw_clip_edges(INMOST_DATA_ENUM_TYPE pace, Storage::real n[3], bool elevation);
		INMOST_DATA_ENUM_TYPE size();
	};


	class bnd_clipper
	{
		Tag clip_state;
		kdtree * tree;	Mesh * mm;
		HandleType * faces;
		INMOST_DATA_ENUM_TYPE nfaces;
	public:
		~bnd_clipper();
		bnd_clipper(Mesh * m, HandleType * _faces, INMOST_DATA_ENUM_TYPE size);
		Storage::real compute_value(Node n1, Node n2, Storage::real * c1, Storage::real * c2, Storage::real * pnt);
		void clip_plane(Storage::real p[3], Storage::real n[3]);
		void gen_clip(std::vector<face2gl> & out, Storage::real p[3], Storage::real n[3], bool elevation);
		void draw_clip(INMOST_DATA_ENUM_TYPE pace, Storage::real p[3], Storage::real n[3]);
		void draw_clip_edges(INMOST_DATA_ENUM_TYPE pace);
		INMOST_DATA_ENUM_TYPE size();
	};
}

#endif
