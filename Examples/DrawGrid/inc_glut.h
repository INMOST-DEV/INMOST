#ifndef _INC_GLUT_H
#define _INC_GLUT_H

#define MAC_WORKAROUND

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
//#pragma comment(lib,"glut32.lib")
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
#include <iostream>
static void glPrintError()
{
	unsigned err = GL_NO_ERROR;
	do
	{
		err = glGetError();
		switch (err)
		{
		case GL_INVALID_ENUM: std::cout << "invalid enum" << std::endl; break;
		case GL_INVALID_VALUE: std::cout << "invalid value" << std::endl; break;
		case GL_INVALID_OPERATION: std::cout << "invalid operation" << std::endl; break;
		//case GL_INVALID_FRAMEBUFFER_OPERATION:
		case GL_OUT_OF_MEMORY: std::cout << "out of memory" << std::endl; break;
		case GL_STACK_UNDERFLOW: std::cout << "stack underflow" << std::endl; break;
		case GL_STACK_OVERFLOW: std::cout << "stack overflow" << std::endl; break;
		case GL_NO_ERROR: break;
		default: std::cout << "unknown error " << err << std::endl; break;
		}
		if (err != GL_NO_ERROR) throw err;
	} while (err != GL_NO_ERROR);
}

#endif
