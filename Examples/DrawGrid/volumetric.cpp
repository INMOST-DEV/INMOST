#include "volumetric.h"
#include "inc_glut.h"
#include "color_bar.h"

namespace INMOST
{
	point & point::operator = (point const & b)
	{
		coords[0] = b.coords[0];
		coords[1] = b.coords[1];
		coords[2] = b.coords[2];
		diam = b.diam;
		dist = b.dist;
		id = b.id;
		return *this;
	}

	point::point(const point & b)
	{
		coords[0] = b.coords[0];
		coords[1] = b.coords[1];
		coords[2] = b.coords[2];
		diam = b.diam;
		dist = b.dist;
		id = b.id;
	}

	inline static unsigned int flip(const unsigned int * fp) { unsigned int mask = -((int)(*fp >> 31)) | 0x80000000; return *fp ^ mask; }
	inline static unsigned int _0(unsigned int x)	{ return x & 0x7FF; }
	inline static unsigned int _1(unsigned int x)	{ return x >> 11 & 0x7FF; }
	inline static unsigned int _2(unsigned int x)   { return x >> 22; }


	void volumetric::radix_sort_dist(std::vector<point_t> & set)
	{
		static std::vector<point_t> tmp;
		tmp.resize(set.size());
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0, 0, sizeof(unsigned int)*kHist * 3);
		for (i = 0; i < set.size(); i++)
		{
			unsigned int fi = flip((unsigned int *)&set[i].dist);
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
		for (i = 0; i < set.size(); i++) tmp[++b0[_0(flip((unsigned int *)&set[i].dist))]] = set[i];
		for (i = 0; i < set.size(); i++) set[++b1[_1(flip((unsigned int *)&tmp[i].dist))]] = tmp[i];
		for (i = 0; i < set.size(); i++) tmp[++b2[_2(flip((unsigned int *)&set[i].dist))]] = set[i];
		for (i = 0; i < set.size(); i++) set[i] = tmp[set.size() - 1 - i];
	}


	volumetric::volumetric(Mesh * _m)
	{
		m = _m;
		const double opt_points = 1000000;
		double density = opt_points / (double)m->NumberOfCells();
		std::cout << "point density " << density << std::endl;
		if( density > 1 )
			points.reserve(ceil(density*m->NumberOfCells()));
		else
			points.reserve(m->NumberOfCells());
		//int q = 0;
		double div = pow(density,1.0/4.0);
		for (Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			point_t p;
			Storage::real cnt[3], cntf[3];
			it->Centroid(cnt);
			p.coords[0] = cnt[0];
			p.coords[1] = cnt[1];
			p.coords[2] = cnt[2];
			p.id = it->LocalID();
			double diam = 0;
			ElementArray<Face> faces = it->getFaces();
			for (ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); ++f)
			{
				f->Centroid(cntf);
				Storage::real d = sqrt((cnt[0] - cntf[0])*(cnt[0] - cntf[0]) + (cnt[1] - cntf[1])*(cnt[1] - cntf[1]) + (cnt[2] - cntf[2])*(cnt[2] - cntf[2]));
				//diam += d;
				if (diam < d) diam = d;
			}
			
			p.diam = diam / div;// / density;
			
			points.push_back(p);
			
			int kmax = (int)(rand()%2) ? ceil(density) : floor(density);
			
			for(int k = 1; k < kmax; ++k)
			{
				double r = diam*(rand()/(double)RAND_MAX*0.4+0.75);
				double phi = rand()/(double)RAND_MAX*2*3.141592;
				double u = rand()/(double)RAND_MAX*2-1;
				p.coords[0] = cnt[0] + sqrt(1-u*u)*cos(phi)*r;
				p.coords[1] = cnt[1] + sqrt(1-u*u)*sin(phi)*r;
				p.coords[2] = cnt[2] + u*r;
				points.push_back(p);
			}
			
		}
		/*
		points.reserve(m->NumberOfNodes()*20);
		int q = 0;
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
		{
		point_t pnt;
		Storage::real_array cnt = it->Coords();
		pnt.coords[0] = cnt[0];
		pnt.coords[1] = cnt[1];
		pnt.coords[2] = cnt[2];
		pnt.id = it->LocalID();
		points.push_back(pnt);

		ElementArray<Element> adj = it->getAdjElements(CELL);
		for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); ++jt)
		{
		Storage::real cnt2[3];
		jt->Centroid(cnt2);
		const int ncoefs = 1;
		const Storage::real coefs[ncoefs] = {0.82};
		for(int j = 0; j < ncoefs; ++j)
		{
		pnt.coords[0] = cnt[0]*(1-coefs[j]) + cnt2[0]*coefs[j];
		pnt.coords[1] = cnt[1]*(1-coefs[j]) + cnt2[1]*coefs[j];
		pnt.coords[2] = cnt[2]*(1-coefs[j]) + cnt2[2]*coefs[j];
		pnt.id = it->LocalID();
		points.push_back(pnt);
		}
		}

		}
		*/
		std::cout << "number of points " << points.size() << "\n";
	}

	void volumetric::camera(double pos[3], int interactive)
	{
		if (interactive) return;
		float posf[3];
		posf[0] = pos[0];
		posf[1] = pos[1];
		posf[2] = pos[2];
		for (size_t k = 0; k < points.size(); ++k)
		{
			points[k].dist = sqrtf(
				(posf[0] - points[k].coords[0])*(posf[0] - points[k].coords[0]) +
				(posf[1] - points[k].coords[1])*(posf[1] - points[k].coords[1]) +
				(posf[2] - points[k].coords[2])*(posf[2] - points[k].coords[2]));
		}
		double t = Timer();
		radix_sort_dist(points);
		//std::sort(points.rbegin(),points.rend());
		std::cout << "Time to sort " << Timer() - t << "\n";
	}

	void volumetric::draw(int interactive) const
	{
		double origin[3], right[3], up[3];
		GLdouble modelview[16], projection[16];
		GLint viewport[4] = { 0, 0, 1, 1 };
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		GLdouble outx, outy, outz;  // Var's to save the answer in
		gluUnProject(0.5, 0.5, 0.,
			modelview, projection, viewport,
			&outx, &outy, &outz);
		origin[0] = outx;
		origin[1] = outy;
		origin[2] = outz;
		gluUnProject(1.0, 0.5, 0.,
			modelview, projection, viewport,
			&outx, &outy, &outz);
		right[0] = outx;
		right[1] = outy;
		right[2] = outz;
		gluUnProject(0.5, 1.0, 0.,
			modelview, projection, viewport,
			&outx, &outy, &outz);
		up[0] = outx;
		up[1] = outy;
		up[2] = outz;
		right[0] -= origin[0];
		right[1] -= origin[1];
		right[2] -= origin[2];
		up[0] -= origin[0];
		up[1] -= origin[1];
		up[2] -= origin[2];

		double l = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
		if (l)
		{
			right[0] /= l;
			right[1] /= l;
			right[2] /= l;
		}
		l = sqrt(up[0] * up[0] + up[1] * up[1] + up[2] * up[2]);
		if (l)
		{
			up[0] /= l;
			up[1] /= l;
			up[2] /= l;
		}

		const float alpha = 0.008f;
		const float mult = 1.0f;
		const float rmult = 0.7f;
		//glPointSize(5.0);
		glColor4f(0.5f, 0.5f, 0.5f, alpha);
		glEnable(GL_BLEND);
		//glBegin(GL_TRIANGLES);
		//glBegin(GL_QUADS);
		if (interactive)
		{
			for (size_t k = 0; k < points.size(); ++k) if (k % 100 == 0)
			{
				if (isColorBarEnabled())
				{
					color_t c = GetColorBar()->pick_color(m->CellByLocalID(points[k].id)->RealDF(GetVisualizationTag()));
					c.a() = alpha;
					c.set_color();
				}

				glBegin(GL_TRIANGLE_FAN);
				glVertex3f(points[k].coords[0], points[k].coords[1], points[k].coords[2]);
				glVertex3f(points[k].coords[0] + (right[0] * rmult + up[0])*points[k].diam*mult, points[k].coords[1] + (right[1] * rmult + up[1])*points[k].diam*mult, points[k].coords[2] + (right[2] * rmult + up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0] + (right[0] * rmult * 2)*points[k].diam*mult, points[k].coords[1] + (right[1] * rmult * 2)*points[k].diam*mult, points[k].coords[2] + (right[2] * rmult * 2)*points[k].diam*mult);
				glVertex3f(points[k].coords[0] + (right[0] * rmult - up[0])*points[k].diam*mult, points[k].coords[1] + (right[1] * rmult - up[1])*points[k].diam*mult, points[k].coords[2] + (right[2] * rmult - up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0] + (-right[0] * rmult - up[0])*points[k].diam*mult, points[k].coords[1] + (-right[1] * rmult - up[1])*points[k].diam*mult, points[k].coords[2] + (-right[2] * rmult - up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0] + (-right[0] * rmult * 2)*points[k].diam*mult, points[k].coords[1] + (-right[1] * rmult * 2)*points[k].diam*mult, points[k].coords[2] + (-right[2] * rmult * 2)*points[k].diam*mult);
				glVertex3f(points[k].coords[0] + (-right[0] * rmult + up[0])*points[k].diam*mult, points[k].coords[1] + (-right[1] * rmult + up[1])*points[k].diam*mult, points[k].coords[2] + (-right[2] * rmult + up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0] + (right[0] * rmult + up[0])*points[k].diam*mult, points[k].coords[1] + (right[1] * rmult + up[1])*points[k].diam*mult, points[k].coords[2] + (right[2] * rmult + up[2])*points[k].diam*mult);
				glEnd();

				/*
				glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-up[0]*2)*points[k].diam*mult,points[k].coords[1]+(-up[1]*2)*points[k].diam*mult,points[k].coords[2]+(-up[2]*2)*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
				*/
				/*
				glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]-up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]-up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
				*/
				//glVertex3f(points[k].coords[0]+right[0]*points[k].diam*mult,points[k].coords[1]+right[1]*points[k].diam*mult,points[k].coords[2]+right[2]*points[k].diam*mult);
				//glVertex3f(points[k].coords[0]+up[0]*points[k].diam*mult,points[k].coords[1]+up[1]*points[k].diam*mult,points[k].coords[2]+up[2]*points[k].diam*mult);
				//glVertex3f(points[k].coords[0]-right[0]*points[k].diam*mult,points[k].coords[1]-right[1]*points[k].diam*mult,points[k].coords[2]-right[2]*points[k].diam*mult);
				//glVertex3f(points[k].coords[0]-up[0]*points[k].diam*mult,points[k].coords[1]-up[1]*points[k].diam*mult,points[k].coords[2]-up[2]*points[k].diam*mult);
			}
		}
		else
		for (size_t k = 0; k < points.size(); ++k)
		{
			if (isColorBarEnabled())
			{
				color_t c = GetColorBar()->pick_color(m->CellByLocalID(points[k].id)->RealDF(GetVisualizationTag()));
				c.a() = alpha;
				c.set_color();
			}

			glBegin(GL_TRIANGLE_FAN);
			glVertex3f(points[k].coords[0], points[k].coords[1], points[k].coords[2]);
			glVertex3f(points[k].coords[0] + (right[0] + up[0])*points[k].diam*mult, points[k].coords[1] + (right[1] + up[1])*points[k].diam*mult, points[k].coords[2] + (right[2] + up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0] + (right[0] * 2)*points[k].diam*mult, points[k].coords[1] + (right[1] * 2)*points[k].diam*mult, points[k].coords[2] + (right[2] * 2)*points[k].diam*mult);
			glVertex3f(points[k].coords[0] + (right[0] - up[0])*points[k].diam*mult, points[k].coords[1] + (right[1] - up[1])*points[k].diam*mult, points[k].coords[2] + (right[2] - up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0] + (-right[0] - up[0])*points[k].diam*mult, points[k].coords[1] + (-right[1] - up[1])*points[k].diam*mult, points[k].coords[2] + (-right[2] - up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0] + (-right[0] * 2)*points[k].diam*mult, points[k].coords[1] + (-right[1] * 2)*points[k].diam*mult, points[k].coords[2] + (-right[2] * 2)*points[k].diam*mult);
			glVertex3f(points[k].coords[0] + (-right[0] + up[0])*points[k].diam*mult, points[k].coords[1] + (-right[1] + up[1])*points[k].diam*mult, points[k].coords[2] + (-right[2] + up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0] + (right[0] + up[0])*points[k].diam*mult, points[k].coords[1] + (right[1] + up[1])*points[k].diam*mult, points[k].coords[2] + (right[2] + up[2])*points[k].diam*mult);
			glEnd();

			/*
			glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-up[0]*2)*points[k].diam*mult,points[k].coords[1]+(-up[1]*2)*points[k].diam*mult,points[k].coords[2]+(-up[2]*2)*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
			*/
			/*
			glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]-up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]-up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
			*/
			//glVertex3f(points[k].coords[0]+right[0]*points[k].diam*mult,points[k].coords[1]+right[1]*points[k].diam*mult,points[k].coords[2]+right[2]*points[k].diam*mult);
			//glVertex3f(points[k].coords[0]+up[0]*points[k].diam*mult,points[k].coords[1]+up[1]*points[k].diam*mult,points[k].coords[2]+up[2]*points[k].diam*mult);
			//glVertex3f(points[k].coords[0]-right[0]*points[k].diam*mult,points[k].coords[1]-right[1]*points[k].diam*mult,points[k].coords[2]-right[2]*points[k].diam*mult);
			//glVertex3f(points[k].coords[0]-up[0]*points[k].diam*mult,points[k].coords[1]-up[1]*points[k].diam*mult,points[k].coords[2]-up[2]*points[k].diam*mult);
		}
		//glEnd();
		glDisable(GL_BLEND);
		//glPointSize(1.0);
	}



	volumetric2::volumetric2(Mesh * _m)
	{
		m = _m;

		faces.reserve(m->NumberOfFaces());
		int q = 0;
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1, std::min<INMOST_DATA_ENUM_TYPE>(15, (unsigned)m->NumberOfFaces() / 100));
		for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
		{
			faces.push_back(DrawFace(it->self()));
			if (q%pace == 0) faces.back().set_flag(true);
			q++;
		}
		std::cout << "number of faces " << faces.size() << "\n";
	}

	void volumetric2::camera(double pos[3], int interactive)
	{
		if (interactive) return;
		for (size_t k = 0; k < faces.size(); ++k)
			faces[k].compute_dist(pos);
		//face2gl::radix_sort_dist(faces);
		//std::sort(faces.rbegin(), faces.rend());
		std::sort(faces.begin(), faces.end());
	}

	void volumetric2::draw(int interactive)
	{
		if (interactive)
			draw_faces_interactive_alpha(faces, 0.015);
		else
			draw_faces_alpha(faces, 0.015);
	}
}
