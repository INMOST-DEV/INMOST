#ifndef _COLOR_H
#define _COLOR_H

#include "inc_glut.h"
#include <iostream>
#include <sstream>

struct color_t
{
	float c[4];
	color_t() { memset(c, 0, sizeof(float)* 4); }
	color_t(float r, float g, float b)
	{
		c[0] = r;
		c[1] = g;
		c[2] = b;
		c[3] = 1.0;
	}
	color_t(float r, float g, float b, float a)
	{
		c[0] = r;
		c[1] = g;
		c[2] = b;
		c[3] = a;
	}
	color_t(const color_t & other) { memcpy(c, other.c, sizeof(float)* 4); }
	color_t & operator =(color_t const & other)
	{
		memmove(c, other.c, sizeof(float)* 4);
		return *this;
	}
	void set_color() const { glColor4fv(c); }
	float & r() { return c[0]; }
	float & g() { return c[1]; }
	float & b() { return c[2]; }
	float & a() { return c[3]; }
	float r() const { return c[0]; }
	float g() const { return c[1]; }
	float b() const { return c[2]; }
	float a() const { return c[3]; }
	color_t operator *(float mult)const { return color_t(c[0] * mult, c[1] * mult, c[2] * mult, c[3] * mult); }
	color_t operator +(color_t other)const{ return color_t(c[0] + other.c[0], c[1] + other.c[1], c[2] + other.c[2], c[3] + other.c[3]); }
	color_t operator -(color_t other)const { return color_t(c[0] - other.c[0], c[1] - other.c[1], c[2] - other.c[2], other.c[3]); }
	std::string svg_rgb() const
	{
		std::stringstream out;
		out << "rgb(" << floor(r() * 255) << ", " << floor(g() * 255) << ", " << floor(b() * 255) << ")";
		return out.str();
	}

	std::string svg_rgba_fill() const
	{
		std::stringstream out;
		out << "fill=\"" << svg_rgb() << "\" fill-opacity=\"" << a() << "\"";
		return out.str();
	}

	std::string svg_rgba_stroke() const
	{
		std::stringstream out;
		out << "stroke=\"" << svg_rgb() << "\" stroke-opacity=\"" << a() << "\"";
		return out.str();
	}
};

#endif