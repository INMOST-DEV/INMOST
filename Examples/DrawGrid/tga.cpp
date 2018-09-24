#include "tga.h"
#include <stdio.h>

bool write_tga(const char *filename, int w, int h, char *buffer)
{
	FILE *f = fopen(filename, "wb");
	if (!f)
		return false;
	putc(0, f);
	putc(0, f);
	putc(2, f);
	putc(0, f); putc(0, f);
	putc(0, f); putc(0, f);
	putc(0, f);
	putc(0, f); putc(0, f);
	putc(0, f); putc(0, f);
	putc((w & 0x00ff), f);
	putc((w & 0xff00) / 256, f);
	putc((h & 0x00ff), f);
	putc((h & 0xff00) / 256, f);
	putc(24, f);
	putc(0, f);
	size_t buflen = w * h * 3;
	fwrite(buffer, 1, buflen, f);
	fclose(f);
	return true;
}