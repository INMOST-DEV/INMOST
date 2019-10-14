#include "inc_glut.h"
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
int font_size = 10;
void printtext(const char * fmt, ...)
{
	void * font;
	if( font_size <= 10 )
		font = GLUT_BITMAP_TIMES_ROMAN_10;
	else if ( font_size <= 12 )
		font = GLUT_BITMAP_HELVETICA_12;
	else if ( font_size <= 18 )
		font = GLUT_BITMAP_HELVETICA_18;
	else
		font = GLUT_BITMAP_TIMES_ROMAN_24;
	unsigned int i;
	char stext[131072];
	va_list ap;
	if (fmt == NULL) return;
	va_start(ap, fmt);
	vsprintf(stext, fmt, ap);
	va_end(ap);
	for (i = 0; i<strlen(stext); i++)
	{
		glutBitmapCharacter(font,
			stext[i]);
	}
}
