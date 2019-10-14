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

#include "inc_glut.h"
#include "rotate.h"
#include "math.h"
#include <vector>

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
extern int interactive;
//
std::vector<quaternion> storage;

void clickmotion(int nmx, int nmy) // Mouse
{
	struct vector n;
	double norm,length,t;
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
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
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	mx = 2.*(nmx/(double)width - 0.5);
	my = 2.*(0.5 - nmy/(double)height);
}
void click(int b, int s, int nmx, int nmy) // Mouse
{
	double norm,length;
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
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

void quatpush()
{
	storage.push_back(q);
}

void quatpop()
{
	if( !storage.empty() )
	{
		q = storage.back();
		storage.pop_back();
	}
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

void reverse_rotatevector(double * vec)
{
	int i;
	double rot[16];
	double temp[4] = {vec[0],vec[1],vec[2],1.0};
	double ret[4];
	q.x = -q.x;
	q.y = -q.y;
	q.z = -q.z;
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
	q.x = -q.x;
	q.y = -q.y;
	q.z = -q.z;
}

void revrotatevector(double * vec)
{
	int i;
	double rot[16];
	double temp[4] = {vec[0],vec[1],vec[2],1.0};
	double ret[4];
	q.w = -q.w;
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
	q.w = -q.w;
}


void rotatevector_from_stack(double * vec)
{
	int i;
	struct quaternion q = storage.back();
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

void quatget(double *vec)
{
	vec[0] = q.x / q.w;
	vec[1] = q.y / q.w;
	vec[2] = q.z / q.w;
}

void quatget4(double gq[4])
{
	gq[0] = q.x;
	gq[1] = q.y;
	gq[2] = q.z;
	gq[3] = q.w;
}

void quatget4_from_stack(double gq[4])
{
	gq[0] = storage.back().x;
	gq[1] = storage.back().y;
	gq[2] = storage.back().z;
	gq[3] = storage.back().w;
}

void quatset4(double gq[4])
{
	q.x = gq[0];
	q.y = gq[1];
	q.z = gq[2];
	q.w = gq[3];
}


void reverse_rotatevector_from_stack(double * vec)
{
	int i;
	struct quaternion q = storage.back();
	q.x = -q.x;
	q.y = -q.y;
	q.z = -q.z;
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


void rotate_from_stack()
{
	double rot[16];
	struct quaternion q = storage.back();
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

