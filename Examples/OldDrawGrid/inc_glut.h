#ifndef _INC_GLUT_H
#define _INC_GLUT_H

#if defined (__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#endif
#if defined(_WIN32)
//#include <GL/glut.h>
// В windows следует скачать тот самый glut по ссылке
// http://www.opengl.org/resources/libraries/glut/
// и положить в папку с компилирующимся файлом
// так же следует установить пару библиотек
#define NOMINMAX
#include <windows.h>
#include <GL/glut.h>
//#include "glut.h"
#pragma comment(lib,"glut32.lib")
#endif
#if defined(__linux__)
#include <GL/glut.h>
#endif

void printtext(const char * fmt, ...); //in printtext.cpp

static void glVertexNdv(double * v, int N)
{
	if( N == 2 ) glVertex2dv(v);
	else glVertex3dv(v);
}

#endif
