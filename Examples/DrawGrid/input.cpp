#include "input.h"
#include "inc_glut.h"
#include <stdlib.h>
#include "clipboard.h"
#include <iostream>
/*
#if defined(__APPLE__) || defined(MACOSX)
#include <Carbon/Carbon.h>

KeyMap keyStates;
bool IS_KEYDOWN( uint16_t vKey )
{
	uint8_t index = vKey / 32 ;
	uint8_t shift = vKey % 32 ;
	return keyStates[index].bigEndianValue & (1 << shift) ;
}

bool getGetCommandModifier()
{
	// This grabs all key states, then checks if you were holding down command or not
	GetKeys(keyStates) ;
	if( IS_KEYDOWN( kVK_Command ) )
		return true;
	return false;
}
#endif
*/
Input::Input(int * val, std::string comment) : comment(comment) 
{ 
	input_link = val; 
	type = Integer; 
	canceled = false; 
	done = false;
	noctrl = false;
	str = ""; 
}
Input::Input(double * val, std::string comment) : comment(comment) 
{ 
	input_link = val; 
	type = Double; 
	canceled = false; 
	done = false;
 noctrl = false;
	str = ""; 
}
Input::Input(char * val, std::string comment) : comment(comment) 
{ 
	input_link = val; 
	type = String; 
	canceled = false; 
	done = false;
	noctrl = false;
	str = ""; 
}
Input::Input(void * link, InputType type, std::string comment) : comment(comment), input_link(link), type(type) 
{ 
	canceled = false; 
	done = false;
	noctrl = false;
	str = ""; 
}
Input::Input(const Input & other) :str(other.str), comment(other.comment), input_link(other.input_link), type(other.type), done(other.done), canceled(other.canceled), noctrl(other.noctrl)
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
	noctrl = other.noctrl;
	return *this; 
}


void Input::KeyPress(char c)
{
	bool ctrl = false;
//#if defined(__APPLE__) || defined(MACOSX)
//	ctrl = getGetCommandModifier();
//#else
	ctrl = glutGetModifiers() & (GLUT_ACTIVE_CTRL);
//#endif
	if( (c == 'v' || c == 'V' || c == 22) && ctrl && !noctrl ) //paste
	{
		std::string paste = getTextFromPasteboard();
		std::cout << "paste: " << paste << std::endl;
		if( !paste.empty() )
		{
			noctrl = true;
			for(size_t k = 0; k < paste.length(); ++k) KeyPress(paste[k]);
			noctrl = false;
		}
	}
	else if( (c == 'c' || c == 'C' || c == 3) && ctrl ) //copy
	{
		std::string copy = GetString();
		std::cout << "copy: " << copy << std::endl;
		setTextToPasteboard(copy);
	}
	else if (c == 13)
	{

		done = true;
		if (type == Double) *((double *)input_link) = atof(str.c_str());
		else if (type == Integer) *((int *)input_link) = atoi(str.c_str());
		else if (type == String) strcpy((char *)input_link, str.c_str());
		glutPostRedisplay();
	}
#if defined(__APPLE__) || defined(MACOSX)
	else if (c == 127 || c == 8)
#else
	else if (c == 8)
#endif
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
