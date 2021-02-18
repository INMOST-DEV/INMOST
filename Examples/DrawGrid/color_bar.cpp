#include "color_bar.h"

#include "inc_glut.h"
#include "svg_line.h"
#include <algorithm>

namespace INMOST
{
	color_bar * color_bar::CommonColorBar = NULL;
	Tag color_bar::vtag = Tag();
	ElementType color_bar::vtype = NONE;
	bool color_bar::smooth = false;

	ElementType GetVisualizationType()
	{
		return color_bar::GetVisualizationType();
	}

	bool isVisualizationSmooth()
	{
		return color_bar::isVisualizationSmooth();
	}

	void color_bar::InitColorBar()
	{
		if (CommonColorBar)
			delete CommonColorBar;
		CommonColorBar = new color_bar;
	}
	void color_bar::DestroyColorBar()
	{
		delete CommonColorBar;
	}

	color_bar * color_bar::GetColorBar()
	{
		return CommonColorBar;
	}

	color_bar::~color_bar()
	{
		if (glIsTexture(texture))
			glDeleteTextures(1, &texture);
	}

	color_bar::color_bar()
	{
		min = 0;
		max = 1;
		comment = "";

		/*
		ticks.push_back(0.f);
		ticks.push_back(0.2f);
		ticks.push_back(0.4f);
		ticks.push_back(0.6f);
		ticks.push_back(0.8f);
		ticks.push_back(1.f);


		//colors.push_back(color_t(1,0,0));
		//colors.push_back(color_t(1,1,0));
		//colors.push_back(color_t(0,1,0));
		//colors.push_back(color_t(0,1,1));
		//colors.push_back(color_t(0,0,1));
		//colors.push_back(color_t(1,0,1));

		colors.push_back(color_t(1,0,1));
		colors.push_back(color_t(0,0,1));
		colors.push_back(color_t(0,1,1));
		colors.push_back(color_t(0,1,0));
		colors.push_back(color_t(1,1,0));
		colors.push_back(color_t(1,0,0));
		*/

		//inversed gnuplot color scheme
		ticks.push_back(0.f);
		ticks.push_back(0.05f);
		ticks.push_back(0.5f);
		ticks.push_back(0.75f);
		ticks.push_back(0.95f);
		ticks.push_back(1.f);

		colors.push_back(color_t(1, 1, 1));
		colors.push_back(color_t(1, 1, 0));
		colors.push_back(color_t(0.85, 0, 0));
		colors.push_back(color_t(0.65, 0.25, 0.85));
		colors.push_back(color_t(0.45, 0, 0.55));
		colors.push_back(color_t(0, 0, 0));

		
		//colors.push_back(color_t(1,0,0));
	}


	color_t color_bar::pick_color(float value) const
	{
		if (value < min)
			return color_t(0.4, 1.0, 0.4);
		if (value > max)
			return color_t(0, 0.6, 0);
		float t = (value - min) / (max - min);
		std::vector<float>::const_iterator it = std::lower_bound(ticks.begin(), ticks.end(), t);
		size_t pos = it - ticks.begin();
		if (it == ticks.end() || pos >= ticks.size())
		{
			return colors.back();
		}
		if (pos == 0)
		{
			return colors[0];
		}
		float interp = (t - ticks[pos - 1]) / (ticks[pos] - ticks[pos - 1]);
		return (colors[pos] * interp + colors[pos - 1] * (1 - interp));
	}
	
	
	void color_bar::InitTexture()
	{
		samples = 512	;

		float * pixel_array = new float[(samples + 2) * 4];



		for (int q = 0; q < samples + 2; ++q)
		{
			float t = 1.0f*q / static_cast<float>(samples + 1);
			color_t c = pick_color(t*(max-min)+min);
			//countour lines
			//if( ((q+1) % 128 == 0 || (q+1) % 128 == 127) && (q+1) < samples )
			//	c = pick_color(1-t) + color_t(0,2*t*(1-t),0);

			pixel_array[(q)* 4 + 0] = c.r();
			pixel_array[(q)* 4 + 1] = c.g();
			pixel_array[(q)* 4 + 2] = c.b();
			pixel_array[(q)* 4 + 3] = c.a();
		}

		pixel_array[0] = 0;
		pixel_array[1] = 1;
		pixel_array[2] = 0;
		pixel_array[3] = 1;

		pixel_array[(samples + 1) * 4 + 0] = 0;
		pixel_array[(samples + 1) * 4 + 1] = 1;
		pixel_array[(samples + 1) * 4 + 2] = 0;
		pixel_array[(samples + 1) * 4 + 3] = 1;

		glPrintError();


		//glEnable(GL_TEXTURE);
		//glPrintError();
		glEnable(GL_TEXTURE_1D);
		glPrintError();
		glGenTextures(1, &texture);
		glPrintError();
		glBindTexture(GL_TEXTURE_1D, texture);
		glPrintError();
		//glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

		glTexImage1D(GL_TEXTURE_1D, 0, 4, samples + 2, 1, GL_RGBA, GL_FLOAT, pixel_array);
		glPrintError();

		std::cout << "Created texture " << texture << std::endl;

		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glPrintError();
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glPrintError();
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glPrintError();
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

		glPrintError();


		UnbindTexture();

		delete[] pixel_array;
	}

	void color_bar::BindTexture()
	{
		//glDisable(GL_TEXTURE_GEN_S ); 
		//glDisable(GL_TEXTURE_2D);
		glEnable(GL_TEXTURE_1D);
		glPrintError();
		glBindTexture(GL_TEXTURE_1D, texture);
		glPrintError();

	}

	void color_bar::UnbindTexture()
	{
		//glDisable(GL_TEXTURE_1D);
		glDisable(GL_TEXTURE_1D);
		glPrintError();

	}

	double color_bar::pick_texture(double value) const
	{
		double eps = 1.0 / static_cast<double>(samples);
		return (value - min) / (max - min)*(1 - 2 * eps) + eps;
		//return std::max(std::min((value-min)/(max-min),0.99),0.01);
	}

	void color_bar::Draw()
	{
		float text_pos = -0.89;
		float left = -0.95;
		float right = -0.9;
		float bottom = -0.75;
		float top = 0.75;
		BindTexture();
		glBegin(GL_QUADS);

		glTexCoord1d(1.0 / 1024.0);
		glVertex2f(left, bottom);
		glVertex2f(right, bottom);
		glTexCoord1d(1.0);
		glVertex2f(right, top);
		glVertex2f(left, top);
		/*
		for(int i = 0; i < ticks.size()-1; ++i)
		{
		colors[i].set_color();
		glVertex2f(left,bottom+ticks[i]*(top-bottom));
		glVertex2f(right,bottom+ticks[i]*(top-bottom));
		colors[i+1].set_color();
		glVertex2f(right,bottom+ticks[(i+1)]*(top-bottom));
		glVertex2f(left,bottom+ticks[(i+1)]*(top-bottom));
		}
		*/
		glEnd();
		UnbindTexture();

		int tickmarks = 11;

		glColor4f(0, 0, 0, 1);
		for (int i = 0; i < tickmarks; ++i)
		{
			float t = 1.0f*i / static_cast<float>(tickmarks - 1);
			glRasterPos2f(text_pos, bottom + t*(top - bottom));
			printtext("%g", min + t*(max - min));
		}
		if (comment != "")
		{
			glRasterPos2f(left, bottom - 0.04);
			printtext("%s", comment.c_str());
		}

		glBegin(GL_LINE_LOOP);
		glVertex2f(left, bottom);
		glVertex2f(right, bottom);
		glVertex2f(right, top);
		glVertex2f(left, top);
		glEnd();

		glBegin(GL_LINES);
		for (int i = 0; i < tickmarks; ++i)
		{
			float t = 1.0f*i / static_cast<float>(tickmarks - 1);
			float pos = bottom + t*(top - bottom);
			glVertex2f(left, pos);
			glVertex2f(left + (right - left)*0.2, pos);

			glVertex2f(right + (left - right)*0.25, pos);
			glVertex2f(right, pos);
		}
		glEnd();
	}

	void color_bar::DrawSVG(std::ostream & file, double modelview[16], double projection[16], int viewport[4])
	{
		float text_pos = -0.89;
		float left = -0.95;
		float right = -0.9;
		float bottom = -0.75;
		float top = 0.75;
		double px, py;
		double px1, py1, z;
		double px2, py2;
		int height = glutGet(GLUT_WINDOW_HEIGHT);
		
		for (size_t i = 0; i < ticks.size() - 1; ++i)
		{
			file << "<g>" << std::endl;
			colors[i].set_color();
			gluProject(left, bottom + ticks[i] * (top - bottom), 0, modelview, projection, viewport, &px1, &py1, &z); py1 = height - py1;
			gluProject(right, bottom + ticks[i + 1] * (top - bottom), 0, modelview, projection, viewport, &px2, &py2, &z); py2 = height - py2;
			file << "<defs>" << std::endl;
			file << "<linearGradient id=\"grad" << i << "\" x1=\"0%\" y1=\"100%\" x2=\"0%\" y2=\"0%\">" << std::endl;
			file << "<stop offset=\"0%\" style=\"stop-color:" << colors[i].svg_rgb() << ";stop-opacity:1\"/>" << std::endl;
			file << "<stop offset=\"100%\" style=\"stop-color:" << colors[i + 1].svg_rgb() << ";stop-opacity:1\"/>" << std::endl;
			file << "</linearGradient>" << std::endl;
			file << "</defs>" << std::endl;
			file << "<rect stroke=\"none\" x=\"" << px1 << "\" y=\"" << py2 << "\" width=\"" << fabs(px2 - px1) << "\" height=\"" << fabs(py1 - py2) << "\" fill=\"url(#grad" << i << ")\"/>" << std::endl;
			file << "</g>" << std::endl;
		}

		int tickmarks = 11;

		file << "<g>" << std::endl;
		for (int i = 0; i < tickmarks; ++i)
		{
			float t = 1.0f*i / static_cast<float>(tickmarks - 1);
			gluProject(text_pos, bottom + t*(top - bottom), 0, modelview, projection, viewport, &px, &py, &z); py = height - py;
			file << "<text x=\"" << px << "\" y=\"" << py << "\">" << min + t*(max - min) << "</text>" << std::endl;
		}
		if (comment != "")
		{
			gluProject(left, bottom - 0.04, 0, modelview, projection, viewport, &px, &py, &z); py = height - py;
			file << "<text x=\"" << px << "\" y=\"" << py << "\">" << comment.c_str() << "</text>" << std::endl;
		}

		gluProject(left, bottom, 0, modelview, projection, viewport, &px1, &py1, &z); py1 = height - py1;
		gluProject(right, top, 0, modelview, projection, viewport, &px2, &py2, &z); py2 = height - py2;
		file << "<rect stroke=\"black\" fill=\"none\" x=\"" << px1 << "\" y=\"" << py2 << "\" width=\"" << px2 - px1 << "\" height=\"" << py1 - py2 << "\"/>" << std::endl;
		file << "</g>" << std::endl;

		file << "<g stroke=\"black\">" << std::endl;
		for (int i = 0; i < tickmarks; ++i)
		{
			float t = 1.0f*i / static_cast<float>(tickmarks - 1);
			float pos = bottom + t*(top - bottom);
			svg_line(file, left, pos, 0, left + (right - left)*0.2, pos, 0, modelview, projection, viewport);
			svg_line(file, right + (left - right)*0.25, pos, 0, right, pos, 0, modelview, projection, viewport);
		}
		file << "</g>" << std::endl;
		//file << "</g>" << std::endl;
	}


	color_bar * GetColorBar() { return color_bar::GetColorBar(); }
	bool isColorBarEnabled() { return color_bar::isColorBarEnabled(); }
	Tag GetVisualizationTag() { return color_bar::GetVisualizationTag(); }
}
