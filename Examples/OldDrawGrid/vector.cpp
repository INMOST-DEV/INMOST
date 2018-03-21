#include "vector.h"
#include "inc_glut.h"
#include "color_bar.h"
#include "svg_line.h"

namespace INMOST
{
	
	

	static GLUquadric * cylqs = NULL;
	static void drawcylinder(coord a, coord b, double width)
	{
		double matrix[16];
		if (cylqs == NULL)
		{
			cylqs = gluNewQuadric();
			gluQuadricNormals(cylqs, GLU_SMOOTH);
			gluQuadricOrientation(cylqs, GLU_OUTSIDE);
			gluQuadricDrawStyle(cylqs, GLU_FILL);//GLU_SILHOUETTE
		}
		glPushMatrix();
		glTranslated(a[0], a[1], a[2]);
		get_matrix(a, b, matrix);
		glMultMatrixd(matrix);
		gluCylinder(cylqs, width, width, sqrt((b - a) ^ (b - a)), 4, 2);
		glPopMatrix();
	}

	Vectors::Vectors(Mesh * m, TagRealArray t, ElementType etype) : m(m), etype(etype)
	{
		double cmin[3] = { 1.0e20, 1.0e20, 1.0e20};
		double cmax[3] = {-1.0e20,-1.0e20,-1.0e20};
		for(int k = 0; k < m->NodeLastLocalID(); ++k) if( m->isValidNode(k) )
		{
			Node n = m->NodeByLocalID(k);
			for(int l = 0; l < 3; ++l)
			{
				if( cmin[l] > n->Coords()[l] ) cmin[l] = n->Coords()[l];
				if( cmax[l] < n->Coords()[l] ) cmax[l] = n->Coords()[l];
			}
		}
		double diag = 0;
		for(int l = 0; l < 3; ++l)
			diag += (cmax[l]-cmin[l])*(cmax[l]-cmin[l]);
		diag = sqrt(diag); //estimate mesh size
		max_length = 0;
		for(int k = 0; k < m->LastLocalID(etype); ++k) if( m->isValidElement(etype,k) )
		{
			Element e = m->ElementByLocalID(etype,k);
			if( e.HaveData(t) )
			{
				vec_t add;
				e->Centroid(add.cnt.data());
				add.eid = k;
				for(int l = 0; l < std::min(t[e].size(),3u); ++l)
					add.dir[l] = t[e][l];
				double l = add.dir.length();
				if( l > max_length ) max_length = l;
				add.dir /= l;
				add.length = l;
				vecs.push_back(add);
			}
		}
		scale = diag/100.0;
		std::cout << "scale: " << scale << std::endl;
	}
	
	void Vectors::Draw(int reduced)
	{
		glBegin(GL_LINES);
		int pace = 1;
		if( reduced )
			pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,(unsigned)vecs.size()/100));
			for (unsigned int i = 0; i < vecs.size(); i+=pace)
			{
				if (isColorBarEnabled())
				{
					color_t c = GetColorBar()->pick_color(m->ElementByLocalID(etype,vecs[i].eid)->RealDF(GetVisualizationTag()));
					c.set_color();
				}
				glVertex3dv(vecs[i].cnt.data());
				glVertex3dv((vecs[i].cnt+vecs[i].dir*scale*vecs[i].length/max_length).data());
			}
		glEnd();
	}


	void Vectors::SVGDraw(std::ostream & file, double modelview[16], double projection[16], int viewport[4])
	{
		for (unsigned int i = 0; i < vecs.size(); i++)
		{
			coord p2 = (vecs[i].cnt+vecs[i].dir*scale*vecs[i].length/max_length);
			double * v0 = vecs[i].cnt.data();
			double * v1 = p2.data();
			
			if (isColorBarEnabled())
			{
				color_t c = GetColorBar()->pick_color(m->ElementByLocalID(etype,vecs[i].eid)->RealDF(GetVisualizationTag()));
				file << "<g stroke=\"" << c.svg_rgb() << "\">" << std::endl;
			}
			else file << "<g>" << std::endl;
			svg_line(file, v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], modelview, projection, viewport);
			file << "</g>" << std::endl;
		}
	}
}
