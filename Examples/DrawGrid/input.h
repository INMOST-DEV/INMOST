#ifndef _INPUT_H
#define _INPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

class Input
{
public:
	enum InputType { Double, Integer, String };
private:
	std::string str;
	std::string comment;
	void * input_link;
	InputType type;
	bool done;
	bool canceled;
	bool noctrl;
public:
	Input(int * val, std::string comment);
	Input(double * val, std::string comment);
	Input(char * val, std::string comment);
	Input(void * link, InputType type, std::string comment);
	Input(const Input & other);
	Input & operator =(Input const & other);
	~Input() {}
	void KeyPress(char c);
	bool Done() { return done; }
	bool Canceled() { return canceled; }
	void Draw();
	std::string GetString() { return str; }
};

#endif
