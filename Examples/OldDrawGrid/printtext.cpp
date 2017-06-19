#include "inc_glut.h"
#include <stdarg.h>
#include <stdio.h>
void printtext(const char * fmt, ...)
{

	unsigned int i;
	char stext[131072];
	va_list ap;
	if (fmt == NULL) return;
	va_start(ap, fmt);
	vsprintf(stext, fmt, ap);
	va_end(ap);
	for (i = 0; i<strlen(stext); i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10,
			stext[i]);
	}
}