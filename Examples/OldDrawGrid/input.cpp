#include "input.h"
#include "inc_glut.h"
#include <stdlib.h>

Input::Input(int * val, std::string comment) : comment(comment) 
{ 
	input_link = val; 
	type = Integer; 
	canceled = false; 
	done = false; 
	str = ""; 
}
Input::Input(double * val, std::string comment) : comment(comment) 
{ 
	input_link = val; 
	type = Double; 
	canceled = false; 
	done = false; 
	str = ""; 
}
Input::Input(char * val, std::string comment) : comment(comment) 
{ 
	input_link = val; 
	type = String; 
	canceled = false; 
	done = false; 
	str = ""; 
}
Input::Input(void * link, InputType type, std::string comment) : comment(comment), input_link(link), type(type) 
{ 
	canceled = false; 
	done = false; 
	str = ""; 
}
Input::Input(const Input & other) :str(other.str), comment(other.comment), input_link(other.input_link), type(other.type), done(other.done), canceled(other.canceled) 
{
}
Input & Input::operator = (Input const & other) 
{ 
	comment = other.comment; 
	input_link = other.input_link; 
	str = other.str; 
	type = other.type; 
	canceled = other.canceled; 
	done = other.done; 
	return *this; 
}


void Input::KeyPress(char c)
{
	if (c == 13)
	{

		done = true;
		if (type == Double) *((double *)input_link) = atof(str.c_str());
		else if (type == Integer) *((int *)input_link) = atoi(str.c_str());
		else if (type == String) strcpy((char *)input_link, str.c_str());
		glutPostRedisplay();
	}
	else if (c == 8)
	{
		if (!str.empty()) str.erase(str.size() - 1);
		glutPostRedisplay();
	}
	else if (c == 27)
	{
		canceled = true;
		done = true;
		glutPostRedisplay();
	}
	else if (type == String || ((c >= '0' && c <= '9') || ((str.empty() || tolower(*str.rbegin()) == 'e') && (c == '+' || c == '-')) || (type == Double && (c == '.' || c == 'e' || c == 'E'))))
	{
		str += c;
		glutPostRedisplay();
	}
}

void Input::Draw()
{
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	float h = 26.0f / (float)height;

	glColor3f(1, 1, 1);
	glBegin(GL_QUADS);
	glVertex3f(-0.99, -0.99, 1);
	glVertex3f(-0.99, -0.99 + h, 1);
	glVertex3f(0.99, -0.99 + h, 1);
	glVertex3f(0.99, -0.99, 1);
	glEnd();
	glColor3f(0, 0, 0);
	glBegin(GL_LINE_LOOP);
	glVertex3f(-0.99, -0.99, 1);
	glVertex3f(-0.99, -0.99 + h, 1);
	glVertex3f(0.99, -0.99 + h, 1);
	glVertex3f(0.99, -0.99, 1);
	glEnd();

	glColor4f(0, 0, 0, 1);
	glRasterPos2d(-0.985, -0.98);
	//printtext(str.c_str());
	char oldval[4096];
	if (type == Double) sprintf(oldval, "%g", *(double*)input_link);
	else if (type == Integer) sprintf(oldval, "%d", *(int*)input_link);
	else if (type == String) sprintf(oldval, "%s", (char *)input_link);
	printtext("input number (%s[%s]:%s): %s", comment.c_str(), oldval, type == Integer ? "integer" : (type == Double ? "double" : "string"), str.c_str());
}