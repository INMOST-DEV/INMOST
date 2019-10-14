#ifndef _COLOR_BAR_H
#define _COLOR_BAR_H

#include "inmost.h"
#include "color.h"
#include <vector>

namespace INMOST
{
	class color_bar
	{
		static color_bar * CommonColorBar;
		static Tag vtag;
		static ElementType vtype;
		static bool smooth;
		float min, max;
		std::vector<float> ticks; //ticks from 0 to 1 for each color
		std::vector<color_t> colors; //4 floats for each tick
		std::string comment;
		unsigned texture;
		int samples;
		void InitTexture();
	public:
		color_bar();
		~color_bar();
		void set_comment(std::string text) { comment = text; }
		void set_min(float newmin) { min = newmin; }
		void set_max(float newmax) { max = newmax; }
		float get_min() { return min; }
		float get_max() { return max; }
		color_t pick_color(float value) const;
		void BindTexture();
		void UnbindTexture();
		double pick_texture(double value) const;
		void Draw();
		void DrawSVG(std::ostream & file, double modelview[16], double projection[16], int viewport[4]);
		static void InitColorBar();
		static void DestroyColorBar();
		static color_bar * GetColorBar();
		static bool isColorBarEnabled() { return vtag.isValid(); }
		static Tag GetVisualizationTag() { return vtag; }
		static ElementType GetVisualizationType() { return vtype; }
		static bool isVisualizationSmooth() { return smooth; }
		static void SetVisualizationTag(Tag t, ElementType et, bool st) { vtag = t; vtype = et; smooth = st; }
		static void UnsetVisualizationTag() { vtag = Tag(); vtype = NONE;}
		static void InitColorBarTexture() { CommonColorBar->InitTexture(); }
	};

	Tag GetVisualizationTag();
	ElementType GetVisualizationType();
	bool isVisualizationSmooth();
	color_bar * GetColorBar();
	bool isColorBarEnabled();
}

#endif
