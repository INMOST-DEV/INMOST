#ifndef _SVG_H
#define _SVG_H


#include <iostream>
void svg_line(std::ostream & file, double x1, double y1, double z1, double x2, double y2, double z2, double modelview[16], double projection[16], int viewport[4]);



#endif