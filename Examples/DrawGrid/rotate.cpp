/*********************************
 Implementation of quaternion-based
 object rotation
 
 Functions
 clickmotion - call when holded mouse moved
 click       - call when user clicks
 motion      - call when mouse just moves
 quatinit    - flush rotation value
 rotate      - multiply GL matrix
 
 Dependency: rotate.h,
 Standard: math.h 
 Specific: glut.h
 **********************************/

#include "my_glut.h"
#include "rotate.h"
#include "math.h"

struct quaternion
{
	double x,y,z,w;
};

struct vector
{
	double x,y,z;
};

// Rotation
struct quaternion q;
struct vector drag, onclick;
double mx,my;
extern int width, height;
extern int interactive;
//


void clickmotion(int nmx, int nmy) // Mouse
{
	struct vector n;
	double norm,length,t;
	mx = 2.*(nmx/(double)width - 0.5);
	my = 2.*(0.5 - nmy/(double)height);
	norm = mx*mx + my*my;
	if( norm > 1.0 )
	{
		length = sqrt(norm);
		drag.x = mx/length;
		drag.y = my/length;
		drag.z = 0.0;
	}
	else
	{
		drag.x = mx;
		drag.y = my;
		drag.z = sqrt(1.0-norm);
	}
	n.x = drag.y*onclick.z - drag.z*onclick.y;
	n.y = drag.z*onclick.x - drag.x*onclick.z;
	n.z = drag.x*onclick.y - drag.y*onclick.x;
	if ( n.x*n.x + n.y*n.y + n.z*n.z > 10e-7 )
	{
		t = drag.x*onclick.x + drag.y*onclick.y + drag.z*onclick.z;
		q.x = + q.x*t + q.y*n.z - q.z*n.y + q.w*n.x;
		q.y = - q.x*n.z + q.y*t + q.z*n.x + q.w*n.y;
		q.z = + q.x*n.y - q.y*n.x + q.z*t + q.w*n.z;
		q.w = - q.x*n.x - q.y*n.y - q.z*n.z + q.w*t;
		onclick.x = drag.x;
		onclick.y = drag.y;
		onclick.z = drag.z;
	}
	glutPostRedisplay();
}
void motion(int nmx, int nmy) // Mouse
{
	mx = 2.*(nmx/(double)width - 0.5);
	my = 2.*(0.5 - nmy/(double)height);
}
void click(int b, int s, int nmx, int nmy) // Mouse
{
	double norm,length;
	switch(b)
	{
		case GLUT_LEFT_BUTTON:
			if( s == GLUT_DOWN ) interactive = 1;
			else interactive = 0;
			mx = 2.*(nmx/(double)width - 0.5);
			my = 2.*(0.5 - nmy/(double)height);
			norm = mx*mx + my*my;
			if( norm > 1.0 )
			{
				length = sqrt(norm);
				drag.x = mx/length;
				drag.y = my/length;
				drag.z = 0.0;
			}
			else	
			{
				drag.x = mx;
				drag.y = my;
				drag.z = sqrt(1.0-norm);
			}
			onclick.x = drag.x;
			onclick.y = drag.y;
			onclick.z = drag.z;
			break;
	}
	glutPostRedisplay();
}


void quatinit()
{
	q.x = 0.0;
	q.y = 0.0;
	q.z = 0.0;
	q.w = 1.0;
}

void rotatevector(double * vec)
{
	int i;
	double rot[16];
	double temp[4] = {vec[0],vec[1],vec[2],1.0};
	double ret[4];
	rot[ 0] = (q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z);
	rot[ 1] = 2.*(q.x*q.y - q.w*q.z);
	rot[ 2] = 2.*(q.x*q.z + q.w*q.y);
	rot[ 3] = 0.0;
	rot[ 4] = 2.*(q.x*q.y + q.w*q.z);
	rot[ 5] = (q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z);
	rot[ 6] = 2.*(q.y*q.z - q.w*q.x);
	rot[ 7] = 0.0;
	rot[ 8] = 2.*(q.x*q.z - q.w*q.y);
	rot[ 9] = 2.*(q.y*q.z + q.w*q.x);
	rot[10] = (q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z);
	rot[11] = 0.0;
	rot[12] = 0.0;
	rot[13] = 0.0;
	rot[14] = 0.0;
	rot[15] = (q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
	for(i=0; i < 4; i++)
	{
		ret[i]  = temp[0] * rot[i*4];
		ret[i] += temp[1] * rot[i*4+1];
		ret[i] += temp[2] * rot[i*4+2];
		ret[i] += temp[3] * rot[i*4+3];
	}
	vec[0] = ret[0]/ret[3];
	vec[1] = ret[1]/ret[3];
	vec[2] = ret[2]/ret[3];
}

void rotate()
{
	double rot[16];
	rot[ 0] = (q.w*q.w + q.x*q.x - q.y*q.y - q.z*q.z);
	rot[ 1] = 2.*(q.x*q.y - q.w*q.z);
	rot[ 2] = 2.*(q.x*q.z + q.w*q.y);
	rot[ 3] = 0.0;
	rot[ 4] = 2.*(q.x*q.y + q.w*q.z);
	rot[ 5] = (q.w*q.w - q.x*q.x + q.y*q.y - q.z*q.z);
	rot[ 6] = 2.*(q.y*q.z - q.w*q.x);
	rot[ 7] = 0.0;
	rot[ 8] = 2.*(q.x*q.z - q.w*q.y);
	rot[ 9] = 2.*(q.y*q.z + q.w*q.x);
	rot[10] = (q.w*q.w - q.x*q.x - q.y*q.y + q.z*q.z);
	rot[11] = 0.0;
	rot[12] = 0.0;
	rot[13] = 0.0;
	rot[14] = 0.0;
	rot[15] = (q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
	glMultMatrixd(rot);
}

