#include "face2gl.h"
#include "color_bar.h"

namespace INMOST
{


	face2gl DrawFace(Element f)
	{
		face2gl ret;
		ElementArray<Node> nodes = f->getNodes();

		if (f->nbAdjElements(CELL) == 0) ret.set_color(1, 0, 0, 0.1);
		else ret.set_color(0, 1, 0, 0.1);

		
		for (ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
			ret.add_vert(&kt->Coords()[0],kt->GetMeshLink()->GetDimensions());
		ret.set_elem(f->GetElementType(), f->LocalID());
		ret.compute_center();
		return ret;
	}

	inline static unsigned int flip(const unsigned int * fp) { unsigned int mask = -((int)(*fp >> 31)) | 0x80000000; return *fp ^ mask; }
	inline static unsigned int _0(unsigned int x)	{ return x & 0x7FF; }
	inline static unsigned int _1(unsigned int x)	{ return x >> 11 & 0x7FF; }
	inline static unsigned int _2(unsigned int x)   { return x >> 22; }


	void face2gl::radix_sort_dist(std::vector<face2gl> & set)
	{
		static std::vector<face2gl> tmp;
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
	void face2gl::shift(double x, double y, double z)
	{
		cnt[0] += x;
		cnt[1] += y;
		cnt[2] += z;
		for (size_t k = 0; k < verts.size(); k += 3)
		{
			verts[k + 0] += x;
			verts[k + 1] += y;
			verts[k + 2] += z;
		}
	}

	face2gl::face2gl() 
	{ 
		etype = NONE; 
		id = 0; 
		dist = 0; 
		flag = false; 
		memset(c, 0, sizeof(double)* 4); 
	}

	face2gl::face2gl(const face2gl & other) :verts(other.verts)
	{
		etype = other.etype;
		id = other.id;
		dist = other.dist;
		cnt[0] = other.cnt[0];
		cnt[1] = other.cnt[1];
		cnt[2] = other.cnt[2];
		c[0] = other.c[0];
		c[1] = other.c[1];
		c[2] = other.c[2];
		c[3] = other.c[3];
		flag = other.flag;
		colors = other.colors;
		texcoords = other.texcoords;
		cntcolor = other.cntcolor;
		cnttexcoord = other.cnttexcoord;
	}

	face2gl & face2gl::operator =(face2gl const & other)
	{
		etype = other.etype;
		id = other.id;
		verts = other.verts;
		dist = other.dist;
		cnt[0] = other.cnt[0];
		cnt[1] = other.cnt[1];
		cnt[2] = other.cnt[2];
		c[0] = other.c[0];
		c[1] = other.c[1];
		c[2] = other.c[2];
		c[3] = other.c[3];
		flag = other.flag;
		colors = other.colors;
		cntcolor = other.cntcolor;
		texcoords = other.texcoords;
		cnttexcoord = other.cnttexcoord;
		return *this;
	}

	face2gl::~face2gl() {}

	void face2gl::draw_colour() const
	{
		if (texcoords.empty())
		{
			if (colors.empty())
			{
				glColor4dv(c);
				for (unsigned k = 0; k < verts.size(); k += 3)
				{
					glVertex3dv(cnt);
					glVertex3dv(&verts[k]);
					glVertex3dv(&verts[(k + 3) % verts.size()]);
				}
			}
			else
			{
				for (unsigned k = 0; k < verts.size(); k += 3)
				{
					cntcolor.set_color();
					glVertex3dv(cnt);
					colors[k / 3].set_color();
					glVertex3dv(&verts[k]);
					colors[(k / 3 + 1) % colors.size()].set_color();
					glVertex3dv(&verts[(k + 3) % verts.size()]);
				}
			}
		}
		else
		{
			for (unsigned k = 0; k < verts.size(); k += 3)
			{
				glTexCoord1d(cnttexcoord);
				glVertex3dv(cnt);
				glTexCoord1d(texcoords[k / 3]);
				glVertex3dv(&verts[k]);
				glTexCoord1d(texcoords[(k / 3 + 1) % texcoords.size()]);
				glVertex3dv(&verts[(k + 3) % verts.size()]);
			}
		}
	}

	void face2gl::svg_draw_colour(std::ostream & file, bool drawedges, double modelview[16], double projection[16], int viewport[4]) const
	{
		int height = glutGet(GLUT_WINDOW_HEIGHT);
		//~ static int shape_id = 0;
		if (colors.empty() && texcoords.empty())
		{
			double pcntx, pcnty, z;
			double pverts1x, pverts1y;
			double pverts2x, pverts2y;
			gluProject(cnt[0], cnt[1], cnt[2], modelview, projection, viewport, &pcntx, &pcnty, &z); pcnty = height - pcnty;
			//file << "<g stroke=\"none\" fill=\"rgb("<<floor(c[0]*255)<<","<<floor(c[1]*255)<<","<<floor(c[2]*255)<<")\" fill-opacity=\"" << c[3] << "\">" << std::endl;
			file << "<g stroke=\"none\" fill=\"white\" fill-opacity=\"" << c[3] << "\">" << std::endl;
			for (unsigned k = 0; k < verts.size(); k += 3)
			{
				gluProject(verts[k + 0], verts[k + 1], verts[k + 2], modelview, projection, viewport, &pverts1x, &pverts1y, &z); pverts1y = height - pverts1y;
				gluProject(verts[(k + 3) % verts.size() + 0], verts[(k + 3) % verts.size() + 1], verts[(k + 3) % verts.size() + 2], modelview, projection, viewport, &pverts2x, &pverts2y, &z); pverts2y = height - pverts2y;
				file << "<polygon points=\"" << pcntx << "," << pcnty << " " << pverts1x << "," << pverts1y << " " << pverts2x << "," << pverts2y << "\"/>" << std::endl;
			}
			file << "</g>" << std::endl;
		}

		else if (!colors.empty())
		{
			double pcntx, pcnty, z;
			double pverts1x, pverts1y;
			double pverts2x, pverts2y;
			gluProject(cnt[0], cnt[1], cnt[2], modelview, projection, viewport, &pcntx, &pcnty, &z); pcnty = height - pcnty;
			file << "<g stroke=\"none\" fill=\"" << cntcolor.svg_rgb() << "\">" << std::endl;
			for (unsigned k = 0; k < verts.size(); k += 3)
			{
				gluProject(verts[k + 0], verts[k + 1], verts[k + 2], modelview, projection, viewport, &pverts1x, &pverts1y, &z); pverts1y = height - pverts1y;
				gluProject(verts[(k + 3) % verts.size() + 0], verts[(k + 3) % verts.size() + 1], verts[(k + 3) % verts.size() + 2], modelview, projection, viewport, &pverts2x, &pverts2y, &z); pverts2y = height - pverts2y;
				file << "<polygon points=\"" << pcntx << "," << pcnty << " " << pverts1x << "," << pverts1y << " " << pverts2x << "," << pverts2y << "\"/>" << std::endl;
			}
			file << "</g>" << std::endl;
		}
		else if (!texcoords.empty())
		{
			double pcntx, pcnty, z;
			double pverts1x, pverts1y;
			double pverts2x, pverts2y;
			color_t c0 = GetColorBar()->pick_color(cnttexcoord*(GetColorBar()->get_max() - GetColorBar()->get_min()) + GetColorBar()->get_min());
			gluProject(cnt[0], cnt[1], cnt[2], modelview, projection, viewport, &pcntx, &pcnty, &z); pcnty = height - pcnty;
			//file << "<!-- " << cnttexcoord << " color " << c0.r() << " " << c0.g() << " " << c0.b() << " " << c0.a() << " -->" << std::endl;
			//file << "<g stroke=\"none\" "<< c0.svg_rgba_fill() <<"\">" << std::endl;
			file << "<g stroke=\"none\" fill=\"" << c0.svg_rgb() << "\">" << std::endl;
			for (unsigned k = 0; k < verts.size(); k += 3)
			{
				gluProject(verts[k + 0], verts[k + 1], verts[k + 2], modelview, projection, viewport, &pverts1x, &pverts1y, &z); pverts1y = height - pverts1y;
				gluProject(verts[(k + 3) % verts.size() + 0], verts[(k + 3) % verts.size() + 1], verts[(k + 3) % verts.size() + 2], modelview, projection, viewport, &pverts2x, &pverts2y, &z); pverts2y = height - pverts2y;
				file << "<polygon points=\"" << pcntx << "," << pcnty << " " << pverts1x << "," << pverts1y << " " << pverts2x << "," << pverts2y << "\"/>" << std::endl;
			}
			file << "</g>" << std::endl;
		}
		/*

		else if( !texcoords.empty() )
		{
		double pverts0x,pverts0y,z;
		double pverts1x,pverts1y;
		double pverts2x,pverts2y;
		color_t c0 = GetColorBar()->pick_color(cnttexcoord*(GetColorBar()->get_max()-GetColorBar()->get_min())+GetColorBar()->get_min()),c1,c2;
		gluProject(cnt[0],cnt[1],cnt[2],modelview,projection,viewport,&pverts0x,&pverts0y,&z); pverts0y = height-pverts0y;
		for(unsigned k = 0; k < verts.size(); k+=3)
		{
		c1 = GetColorBar()->pick_color(texcoords[k/3]*(GetColorBar()->get_max()-GetColorBar()->get_min())+GetColorBar()->get_min());
		c2 = GetColorBar()->pick_color(texcoords[(k/3+1)%texcoords.size()]*(GetColorBar()->get_max()-GetColorBar()->get_min())+GetColorBar()->get_min());
		gluProject(verts[k+0],verts[k+1],verts[k+2],modelview,projection,viewport,&pverts1x,&pverts1y,&z); pverts1y = height-pverts1y;
		gluProject(verts[(k+3)%verts.size()+0],verts[(k+3)%verts.size()+1],verts[(k+3)%verts.size()+2],modelview,projection,viewport,&pverts2x,&pverts2y,&z); pverts2y = height-pverts2y;
		file << "<defs>" << std::endl;
		file << "<linearGradient id=\"fadeA"<<shape_id<<"\" gradientUnits=\"userSpaceOnUse\" x1=\"" << pverts0x << "\" y1=\"" << pverts0y << "\" x2=\"" << (pverts1x+pverts2x)*0.5 << "\" y2=\"" << (pverts1y+pverts2y)*0.5 << "\">" << std::endl;
		file << "<stop offset=\"0%\" stop-color=\"" << c0.svg_rgb() << "\"/>" <<std::endl;
		file << "<stop offset=\"100%\" stop-color=\"" << ((c1+c2)*0.5).svg_rgb() << "\"/>" <<std::endl;
		file << "</linearGradient>"<<std::endl;
		file << "<linearGradient id=\"fadeB"<<shape_id<<"\" gradientUnits=\"userSpaceOnUse\" x1=\"" << pverts1x << "\" y1=\"" << pverts1y << "\" x2=\"" << (pverts0x+pverts2x)*0.5 << "\" y2=\"" << (pverts0y+pverts2y)*0.5 << "\">" << std::endl;
		file << "<stop offset=\"0%\" stop-color=\""<< c1.svg_rgb() <<"\"/>" <<std::endl;
		file << "<stop offset=\"100%\" stop-color=\"" << ((c0+c2)*0.5).svg_rgb() << "\"/>" <<std::endl;
		file << "</linearGradient>"<<std::endl;
		file << "<linearGradient id=\"fadeC"<<shape_id<<"\" gradientUnits=\"userSpaceOnUse\" x1=\"" << pverts2x << "\" y1=\"" << pverts2y << "\" x2=\"" << (pverts0x+pverts1x)*0.5 << "\" y2=\"" << (pverts0y+pverts1y)*0.5 << "\">" << std::endl;
		file << "<stop offset=\"0%\" stop-color=\"" << c2.svg_rgb() << "\"/>" <<std::endl;
		file << "<stop offset=\"100%\" stop-color=\"" << ((c0+c1)*0.5).svg_rgb() << "\"/>" <<std::endl;
		file << "</linearGradient>"<<std::endl;
		file << "<path id=\"pathA"<<shape_id<<"\" d=\"M "<<pverts0x<<","<<pverts0y<<" L "<<pverts1x<<","<<pverts1y<<" "<<pverts2x<<","<<pverts2y<<" Z\" fill=\"url(#fadeA"<<shape_id<<")\"/>" <<std::endl;
		file << "<path id=\"pathB"<<shape_id<<"\" d=\"M "<<pverts0x<<","<<pverts0y<<" L "<<pverts1x<<","<<pverts1y<<" "<<pverts2x<<","<<pverts2y<<" Z\" fill=\"url(#fadeB"<<shape_id<<")\"/>" <<std::endl;
		file << "<filter id=\"filter"<<shape_id<<"\">"<<std::endl;
		file << "<feImage xlink:href=\"#pathA"<<shape_id<<"\" result=\"layerA"<<shape_id<<"\" x=\"0\" y=\"0\"/>" <<std::endl;
		file << "<feImage xlink:href=\"#pathB"<<shape_id<<"\" result=\"layerB"<<shape_id<<"\" x=\"0\" y=\"0\"/>" <<std::endl;
		file << "<feComposite in=\"layerA"<<shape_id<<"\" in2=\"layerB"<<shape_id<<"\" operator=\"arithmetic\" k1=\"0\" k2=\"1\" k3=\"1\" k4=\"0\" result=\"temp"<<shape_id<<"\"/>"<<std::endl;
		file << "<feComposite in=\"temp"<<shape_id<<"\" in2=\"graphic"<<shape_id<<"\" operator=\"arithmetic\" k1=\"0\" k2=\"1\" k3=\"1\" k4=\"0\"/>"<<std::endl;
		file << "</filter>"<<std::endl;
		file << "</defs>" << std::endl;
		file << "<g stroke=\"none\">" << std::endl;
		file << "<path d=\"M "<<pverts0x<<","<<pverts0y<<" L "<<pverts1x<<","<<pverts1y<<" "<<pverts2x<<","<<pverts2y<<" Z\" fill=\"url(#fadeC"<<shape_id<<")\" filter=\"url(#filter"<<shape_id<<")\"/>" <<std::endl;
		file << "</g>" << std::endl;
		shape_id++;
		}
		}
		else if( !colors.empty() )
		{
		double pverts0x,pverts0y,z;
		double pverts1x,pverts1y;
		double pverts2x,pverts2y;
		gluProject(cnt[0],cnt[1],cnt[2],modelview,projection,viewport,&pverts0x,&pverts0y,&z); pverts0y = height-pverts0y;
		for(unsigned k = 0; k < verts.size(); k+=3)
		{
		gluProject(verts[k+0],verts[k+1],verts[k+2],modelview,projection,viewport,&pverts1x,&pverts1y,&z); pverts1y = height-pverts1y;
		gluProject(verts[(k+3)%verts.size()+0],verts[(k+3)%verts.size()+1],verts[(k+3)%verts.size()+2],modelview,projection,viewport,&pverts2x,&pverts2y,&z); pverts2y = height-pverts2y;
		file << "<defs>" << std::endl;
		file << "<linearGradient id=\"fadeA"<<shape_id<<"\" gradientUnits=\"userSpaceOnUse\" x1=\"" << pverts0x << "\" y1=\"" << pverts0y << "\" x2=\"" << (pverts1x+pverts2x)*0.5 << "\" y2=\"" << (pverts1y+pverts2y)*0.5 << "\">" << std::endl;
		file << "<stop offset=\"0%\" stop-color=\"" << cntcolor.svg_rgb() << "\"/>" <<std::endl;
		file << "<stop offset=\"100%\" stop-color=\"#000000\"/>" <<std::endl;
		file << "</linearGradient>"<<std::endl;
		file << "<linearGradient id=\"fadeB"<<shape_id<<"\" gradientUnits=\"userSpaceOnUse\" x1=\"" << pverts1x << "\" y1=\"" << pverts1y << "\" x2=\"" << (pverts0x+pverts2x)*0.5 << "\" y2=\"" << (pverts0y+pverts2y)*0.5 << "\">" << std::endl;
		file << "<stop offset=\"0%\" stop-color=\"" << colors[k/3].svg_rgb() << "\"/>" <<std::endl;
		file << "<stop offset=\"100%\" stop-color=\"#000000\"/>" <<std::endl;
		file << "</linearGradient>"<<std::endl;
		file << "<linearGradient id=\"fadeC"<<shape_id<<"\" gradientUnits=\"userSpaceOnUse\" x1=\"" << pverts2x << "\" y1=\"" << pverts2y << "\" x2=\"" << (pverts0x+pverts1x)*0.5 << "\" y2=\"" << (pverts0y+pverts1y)*0.5 << "\">" << std::endl;
		file << "<stop offset=\"0%\" stop-color=\"" << colors[(k/3+1)%colors.size()].svg_rgb() << "\"/>" <<std::endl;
		file << "<stop offset=\"100%\" stop-color=\"#000000\"/>" <<std::endl;
		file << "</linearGradient>"<<std::endl;
		file << "<path id=\"pathA"<<shape_id<<"\" d=\"M "<<pverts0x<<","<<pverts0y<<" L "<<pverts1x<<","<<pverts1y<<" "<<pverts2x<<","<<pverts2y<<" Z\" fill=\"url(#fadeA"<<shape_id<<")\"/>" <<std::endl;
		file << "<path id=\"pathB"<<shape_id<<"\" d=\"M "<<pverts0x<<","<<pverts0y<<" L "<<pverts1x<<","<<pverts1y<<" "<<pverts2x<<","<<pverts2y<<" Z\" fill=\"url(#fadeB"<<shape_id<<")\"/>" <<std::endl;
		file << "<filter id=\"filter"<<shape_id<<"\">"<<std::endl;
		file << "<feImage xlink:href=\"#pathA"<<shape_id<<"\" result=\"layerA"<<shape_id<<"\" x=\"0\" y=\"0\"/>" <<std::endl;
		file << "<feImage xlink:href=\"#pathB"<<shape_id<<"\" result=\"layerB"<<shape_id<<"\" x=\"0\" y=\"0\"/>" <<std::endl;
		file << "<feComposite in=\"layerA"<<shape_id<<"\" in2=\"layerB"<<shape_id<<"\" operator=\"arithmetic\" k1=\"0\" k2=\"1.0\" k4=\"0\" result=\"temp"<<shape_id<<"\"/>"<<std::endl;
		file << "<feComposite in=\"temp"<<shape_id<<"\" in2=\"graphic"<<shape_id<<"\" operator=\"arithmetic\" k1=\"0\" k2=\"1.0\" k4=\"0\"/>"<<std::endl;
		file << "</filter>"<<std::endl;
		file << "</defs>" << std::endl;
		file << "<g stroke=\"none\" stroke-width=\"0\" shape-rendering=\"crispEdges\">" << std::endl;
		file << "<path d=\"M "<<pverts0x<<","<<pverts0y<<" L "<<pverts1x<<","<<pverts1y<<" "<<pverts2x<<","<<pverts2y<<" Z\" fill=\"url(#fadeC"<<shape_id<<")\" filter=\"url(#filter"<<shape_id<<")\"/>" <<std::endl;
		file << "</g>" << std::endl;
		shape_id++;
		}
		}
		*/

		if (drawedges)
		{
			double px, py, z;
			file << "<g stroke=\"black\">" << std::endl;
			file << "<polyline fill=\"none\" points=\"";
			for (unsigned k = 0; k < verts.size() + 3; k += 3)
			{
				gluProject(verts[k%verts.size() + 0], verts[k%verts.size() + 1], verts[k%verts.size() + 2], modelview, projection, viewport, &px, &py, &z); py = height - py;
				file << px << "," << py << " ";
			}
			file << "\" />" << std::endl;
			file << "</g>" << std::endl;
		}

	}

	void face2gl::draw_colour_alpha(double alpha) const
	{
		if (texcoords.empty())
		{
			if (colors.empty())
			{
				//~ double cc[4] = {c[0],c[1],c[2],alpha};
				//glColor4dv(c);
				for (unsigned k = 0; k < verts.size(); k += 3)
				{
					glVertex3dv(cnt);
					glVertex3dv(&verts[k]);
					glVertex3dv(&verts[(k + 3) % verts.size()]);
				}
			}
			else
			{
				for (unsigned k = 0; k < verts.size(); k += 3)
				{
					color_t t = cntcolor;
					t.a() = alpha;
					t.set_color();
					glVertex3dv(cnt);
					t = colors[k / 3];
					t.a() = alpha;
					t.set_color();
					glVertex3dv(&verts[k]);
					t = colors[(k / 3 + 1) % colors.size()];
					t.a() = alpha;
					t.set_color();
					glVertex3dv(&verts[(k + 3) % verts.size()]);
				}
			}
		}
		else
		{
			for (unsigned k = 0; k < verts.size(); k += 3)
			{
				glTexCoord1d(cnttexcoord);
				glVertex3dv(cnt);
				glTexCoord1d(texcoords[k / 3]);
				glVertex3dv(&verts[k]);
				glTexCoord1d(texcoords[(k / 3 + 1) % texcoords.size()]);
				glVertex3dv(&verts[(k + 3) % verts.size()]);
			}
		}
	}

	void face2gl::draw() const
	{
		for (unsigned k = 0; k < verts.size(); k += 3)
		{
			glVertex3dv(cnt);
			glVertex3dv(&verts[k]);
			glVertex3dv(&verts[(k + 3) % verts.size()]);
		}
	}
	
	__INLINE void vfill4(const double vi[3], double tc, float v[4])
	{
		v[0] = vi[0], v[1] = vi[1], v[2] = vi[2], v[3] = tc;
	}
	__INLINE void vfill3(const double vi[3], float v[3])
	{
		v[0] = vi[0], v[1] = vi[1], v[2] = vi[2];
	}
	
	void face2gl::wgl_draw(std::ostream & file) const
	{
		size_t k;
		for (k = 0; k < verts.size()-3; k += 3)
		{
			file << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << (texcoords.empty() ? 1.0 : cnttexcoord) << std::endl;
			file << verts[k+0] << " " << verts[k+1] << " " << verts[k+2] << " " << (texcoords.empty() ? 1.0 : texcoords[k/3+0]) << std::endl;
			file << verts[k+3] << " " << verts[k+4] << " " << verts[k+5] << " " << (texcoords.empty() ? 1.0 : texcoords[k/3+1]) << std::endl;
		}
		file << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << (texcoords.empty() ? 1.0 : cnttexcoord) << std::endl;
		k = verts.size()-3;
		file << verts[k+0] << " " << verts[k+1] << " " << verts[k+2] << " " << (texcoords.empty() ? 1.0 : texcoords[k/3+0]) << std::endl;
		k = 0;
		file << verts[k+0] << " " << verts[k+1] << " " << verts[k+2] << " " << (texcoords.empty() ? 1.0 : texcoords[k/3+0]) << std::endl;
	}
	
	void face2gl::wgl_draw_bin(std::ostream & file) const
	{
		size_t k;
		float v[4];
		for (k = 0; k < verts.size()-3; k += 3)
		{
			vfill4(cnt        ,texcoords.empty() ? 1.0 : cnttexcoord,v);
			file.write(reinterpret_cast<char *>(v),sizeof(float)*4);
			vfill4(&verts[k+0],texcoords.empty() ? 1.0 : texcoords[k/3+0],v);
			file.write(reinterpret_cast<char *>(v),sizeof(float)*4);
			vfill4(&verts[k+3],texcoords.empty() ? 1.0 : texcoords[k/3+1],v);
			file.write(reinterpret_cast<char *>(v),sizeof(float)*4);
		}
		vfill4(cnt        ,texcoords.empty() ? 1.0 : cnttexcoord,v);
		file.write(reinterpret_cast<char *>(v),sizeof(float)*4);
		k = verts.size()-3;
		vfill4(&verts[k+0],texcoords.empty() ? 1.0 : texcoords[k/3+0],v);
		file.write(reinterpret_cast<char *>(v),sizeof(float)*4);
		k = 0;
		vfill4(&verts[k+0],texcoords.empty() ? 1.0 : texcoords[k/3+0],v);
		file.write(reinterpret_cast<char *>(v),sizeof(float)*4);
	}
	
	void face2gl::wgl_draw_edges(std::ostream & file) const
	{
		size_t k;
		for (k = 0; k < verts.size()-3; k += 3)
		{
			file << verts[k+0] << " " << verts[k+1] << " " << verts[k+2] << std::endl;
			file << verts[k+3] << " " << verts[k+4] << " " << verts[k+5] << std::endl;
		}
		k = verts.size()-3;
		file << verts[k+0] << " " << verts[k+1] << " " << verts[k+2] << std::endl;
		k = 0;
		file << verts[k+0] << " " << verts[k+1] << " " << verts[k+2] << std::endl;
	}
	
	void face2gl::wgl_draw_edges_bin(std::ostream & file) const
	{
		size_t k;
		float v[3];
		for (k = 0; k < verts.size()-3; k += 3)
		{
			vfill3(&verts[k+0],v);
			file.write(reinterpret_cast<char *>(v),sizeof(float)*3);
			vfill3(&verts[k+3],v);
			file.write(reinterpret_cast<char *>(v),sizeof(float)*3);
		}
		k = verts.size()-3;
		vfill3(&verts[k+0],v);
		file.write(reinterpret_cast<char *>(v),sizeof(float)*3);
		k = 0;
		vfill3(&verts[k+0],v);
		file.write(reinterpret_cast<char *>(v),sizeof(float)*3);
	}
	
	

	void face2gl::svg_draw(std::ostream & file, bool drawedges, double modelview[16], double projection[16], int viewport[4]) const
	{
		double pcntx, pcnty, z;
		double pverts1x, pverts1y;
		double pverts2x, pverts2y;
		int height = glutGet(GLUT_WINDOW_HEIGHT);
		gluProject(cnt[0], cnt[1], cnt[2], modelview, projection, viewport, &pcntx, &pcnty, &z); pcnty = height - pcnty;
		file << "<g stroke=\"none\">" << std::endl;
		for (unsigned k = 0; k < verts.size(); k += 3)
		{
			gluProject(verts[k + 0], verts[k + 1], verts[k + 2], modelview, projection, viewport, &pverts1x, &pverts1y, &z); pverts1y = height - pverts1y;
			gluProject(verts[(k + 3) % verts.size() + 0], verts[(k + 3) % verts.size() + 1], verts[(k + 3) % verts.size() + 2], modelview, projection, viewport, &pverts2x, &pverts2y, &z); pverts2y = height - pverts2y;
			file << "<polygon points=\"" << pcntx << "," << pcnty << " " << pverts1x << "," << pverts1y << " " << pverts2x << "," << pverts2y << "\"/>" << std::endl;
		}
		file << "</g>" << std::endl;
		if (drawedges)
		{
			double px, py, z;
			file << "<g stroke=\"black\">" << std::endl;
			file << "<polyline fill=\"none\" points=\"";
			for (unsigned k = 0; k < verts.size(); k += 3)
			{
				gluProject(verts[k%verts.size() + 0], verts[k%verts.size() + 1], verts[k%verts.size() + 2], modelview, projection, viewport, &px, &py, &z); py = height - py;
				file << px << "," << py << " ";
			}
			file << "\" />" << std::endl;
			file << "</g>" << std::endl;
		}
	}

	void face2gl::drawedges() const
	{
		for (unsigned k = 0; k < verts.size(); k += 3)
		{
			glVertex3dv(&verts[k]);
			glVertex3dv(&verts[(k + 3) % verts.size()]);
		}
	}

	void face2gl::svg_drawedges(std::ostream & file, double modelview[16], double projection[16], int viewport[4]) const
	{
		double px, py, z;
		int height = glutGet(GLUT_WINDOW_HEIGHT);
		file << "<polyline fill=\"none\" points=\"";
		for (unsigned k = 0; k < verts.size(); k += 3)
		{
			gluProject(verts[k%verts.size() + 0], verts[k%verts.size() + 1], verts[k%verts.size() + 2], modelview, projection, viewport, &px, &py, &z); py = height - py;
			file << px << "," << py << " ";
		}
		file << "\" />" << std::endl;
	}

	void face2gl::compute_center()
	{
		cnt[0] = cnt[1] = cnt[2] = 0;
		for (INMOST_DATA_ENUM_TYPE k = 0; k < verts.size(); k += 3)
		{
			cnt[0] += verts[k + 0];
			cnt[1] += verts[k + 1];
			cnt[2] += verts[k + 2];
		}
		cnt[0] /= (verts.size() / 3)*1.0;
		cnt[1] /= (verts.size() / 3)*1.0;
		cnt[2] /= (verts.size() / 3)*1.0;
		compute_center_color();
		compute_center_texcoord();
	}

	void face2gl::compute_center_color()
	{
		if (!colors.empty())
		{
			cntcolor.r() = 0;
			cntcolor.g() = 0;
			cntcolor.b() = 0;
			cntcolor.a() = 0;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < colors.size(); k++)
				cntcolor = cntcolor + colors[k];
			cntcolor = cntcolor*(1.0f / static_cast<float>(colors.size()));
		}
	}

	void face2gl::compute_center_texcoord()
	{
		if (!texcoords.empty())
		{
			cnttexcoord = 0.0;
			for (INMOST_DATA_ENUM_TYPE k = 0; k < texcoords.size(); k++)
				cnttexcoord += texcoords[k];
			cnttexcoord /= static_cast<double>(texcoords.size());
		}
	}


	void draw_faces_nc(std::vector<face2gl> & set, int highlight)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glColor4f(0, 1, 0, 0.1);
		glBegin(GL_TRIANGLES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) set[q].draw();
		glEnd();
		if (highlight != -1)
		{
			glColor4f(1, 0, 0, 1);
			glBegin(GL_TRIANGLES);
			set[highlight].draw();
			glEnd();
		}
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	void svg_draw_faces_nc(std::ostream & file, std::vector<face2gl> & set, bool drawedges, double modelview[16], double projection[16], int viewport[4], int highlight)
	{
		
		file << "<g stroke=\"none\" fill=\"green\" fill-opacity=\"0.1\" stroke-opacity=\"0.1\">" << std::endl;
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) set[q].svg_draw(file, drawedges, modelview, projection, viewport);
		file << "</g>" << std::endl;
		if (highlight != -1)
		{
			file << "<g fill=\"red\" fill-opacity=\"1\">" << std::endl;
			set[highlight].svg_draw(file, drawedges, modelview, projection, viewport);
			file << "</g>" << std::endl;
		}
	}

	void draw_faces(std::vector<face2gl> & set, int highlight)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPrintError();
		glPolygonOffset(1.0, 1.0);
		glPrintError();
		if (isColorBarEnabled()) GetColorBar()->BindTexture();
		glPrintError();
		glBegin(GL_TRIANGLES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) 
			set[q].draw_colour();
		glEnd();
		glPrintError();
		if (isColorBarEnabled()) GetColorBar()->UnbindTexture();
		glPrintError();
		if (highlight != -1)
		{
			glColor4f(1, 0, 0, 1);
			glPrintError();
			glBegin(GL_TRIANGLES);
			set[highlight].draw();
			glEnd();
			glPrintError();
		}

		glDisable(GL_POLYGON_OFFSET_FILL);
		glPrintError();

	}
	
	void wgl_draw_faces(std::ostream & file, const std::vector<face2gl> & set)
	{
		for (size_t q = 0; q < set.size(); q++)
			set[q].wgl_draw(file);
	}
	
	void wgl_draw_edges(std::ostream & file, const std::vector<face2gl> & set)
	{
		for (size_t q = 0; q < set.size(); q++)
			set[q].wgl_draw_edges(file);
	}
	
	void wgl_draw_faces_bin(std::ostream & file, const std::vector<face2gl> & set)
	{
		for (size_t q = 0; q < set.size(); q++)
			set[q].wgl_draw_bin(file);
	}
	
	void wgl_draw_edges_bin(std::ostream & file, const std::vector<face2gl> & set)
	{
		for (size_t q = 0; q < set.size(); q++)
			set[q].wgl_draw_edges_bin(file);
	}


	void svg_draw_faces(std::ostream & file, std::vector<face2gl> & set, bool drawedges, double modelview[16], double projection[16], int viewport[4], int highlight)
	{
		

		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) set[q].svg_draw_colour(file, drawedges, modelview, projection, viewport);


		if (highlight != -1)
		{
			file << "<g color=\"red\" fill-opacity=\"1\">" << std::endl;
			set[highlight].svg_draw(file, drawedges, modelview, projection, viewport);
			file << "</g>" << std::endl;
		}



	}

	void draw_faces_alpha(std::vector<face2gl> & set, double alpha)
	{
		if (isColorBarEnabled()) GetColorBar()->BindTexture();

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glBegin(GL_TRIANGLES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) set[q].draw_colour_alpha(alpha);
		glEnd();
		if (isColorBarEnabled()) GetColorBar()->UnbindTexture();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	void draw_edges(std::vector<face2gl> & set, int highlight)
	{
		glBegin(GL_LINES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) set[q].drawedges();
		glEnd();
		if (highlight != -1)
		{
			glColor4f(0, 1, 0, 1);
			glBegin(GL_LINES);
			set[highlight].drawedges();
			glEnd();
		}
	}

	void svg_draw_edges(std::ostream & file, std::vector<face2gl> & set, double modelview[16], double projection[16], int viewport[4], int highlight)
	{
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) set[q].svg_drawedges(file, modelview, projection, viewport);
		if (highlight != -1)
		{
			file << "<g fill=\"green\">" << std::endl;
			set[highlight].svg_drawedges(file, modelview, projection, viewport);
			file << "</g>" << std::endl;
		}
	}

	void draw_faces_interactive_nc(std::vector<face2gl> & set)
	{
		glColor4f(0, 1, 0, 0.1);

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);


		glBegin(GL_TRIANGLES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) if (set[q].get_flag()) set[q].draw();
		glEnd();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	void draw_faces_interactive(std::vector<face2gl> & set)
	{
		if (isColorBarEnabled()) GetColorBar()->BindTexture();

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);


		glBegin(GL_TRIANGLES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) if (set[q].get_flag()) set[q].draw_colour();
		glEnd();
		if (isColorBarEnabled()) GetColorBar()->UnbindTexture();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	void draw_faces_interactive_alpha(std::vector<face2gl> & set, double alpha)
	{
		if (isColorBarEnabled()) GetColorBar()->BindTexture();

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);


		glBegin(GL_TRIANGLES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) if (set[q].get_flag()) set[q].draw_colour_alpha(alpha);
		glEnd();
		if (isColorBarEnabled()) GetColorBar()->UnbindTexture();

		glDisable(GL_POLYGON_OFFSET_FILL);

	}

	void draw_edges_interactive(std::vector<face2gl> & set)
	{
		glBegin(GL_LINES);
		for (INMOST_DATA_ENUM_TYPE q = 0; q < set.size(); q++) if (set[q].get_flag()) set[q].drawedges();
		glEnd();
	}
}
