#include "svg_line.h"
#include "inc_glut.h"


void svg_line(std::ostream & file, double x1, double y1, double z1, double x2, double y2, double z2, double modelview[16], double projection[16], int viewport[4])
{
	double px1, py1, z;
	double px2, py2;
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	gluProject(x1, y1, z1, modelview, projection, viewport, &px1, &py1, &z); py1 = height - py1;
	gluProject(x2, y2, z2, modelview, projection, viewport, &px2, &py2, &z); py2 = height - py2;
	file << "<line x1=\"" << px1 << "\" y1=\"" << py1 << "\" x2=\"" << px2 << "\" y2=\"" << py2 << "\"/>" << std::endl;
}
