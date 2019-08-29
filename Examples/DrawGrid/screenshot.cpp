#include "inc_glut.h"
#include "tga.h"
#include <stdio.h>

extern void draw_screen(); //global openl drawing routine

void screenshot(int tiles)
{
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int oldwidth = width;
	int oldheight = height;
	width *= tiles;
	height *= tiles;

	char * pixelbuffer = new char[width*height * 3 + oldwidth*oldheight * 3];
	char * tempbuffer = pixelbuffer + width*height * 3;
	//int WH[2];
	//glGetIntegerv(GL_VIEWPORT,WH);
	//printf("W %d H %d\n",WH[0],WH[1]);


	for (int i = 0; i < tiles; ++i)
	{
		for (int j = 0; j < tiles; ++j)
		{
			glViewport(-oldwidth*i, -oldheight*j, width, height);
			draw_screen();
			glReadBuffer(GL_BACK);
			glReadPixels(0, 0, oldwidth, oldheight, GL_BGR_EXT, GL_UNSIGNED_BYTE, tempbuffer);

			int koff = oldwidth*(i);
			int loff = oldheight*(j);

			for (int l = 0; l < oldheight; ++l)
			for (int k = 0; k < oldwidth; ++k)
			for (int m = 0; m < 3; ++m)
				pixelbuffer[((koff + k) + (loff + l)*width) * 3 + m] = tempbuffer[(k + l*oldwidth) * 3 + m];

			//filename[0] += i;
			//filename[1] += j;
			//write_tga(filename,width,height,pixelbuffer);
			//filename[0] = filename[1] = '0';
		}
	}

	/*
	glViewport(-oldwidth*3,-oldheight*3,width,height);
	draw_screen();
	glReadBuffer(GL_BACK);
	glReadPixels(0,0,width,height,GL_BGR_EXT,GL_UNSIGNED_BYTE,pixelbuffer);
	*/





	write_tga("screenshot.tga", width, height, pixelbuffer);
	delete[] pixelbuffer;
	width = oldwidth;
	height = oldheight;
//#if !(defined(MAC_WORKAROUND) && (defined (__APPLE__) || defined(MACOSX)))
	glViewport(0, 0, width, height);
//#endf //MAC_WORKAROUND
}
