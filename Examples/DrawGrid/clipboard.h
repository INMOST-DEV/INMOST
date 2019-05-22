#ifndef _CLIPBOARD_H
#define _CLIPBOARD_H

#include <string>

bool setTextToPasteboard(std::string str);
std::string getTextFromPasteboard();

#endif
