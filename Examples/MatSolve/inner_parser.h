#ifndef INMOST_INNER_PARSER_H_H
#define INMOST_INNER_PARSER_H_H

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cstdio>

enum OptionType {
    REAL,
    ENUM
};

class InnerOption {
public:

    InnerOption(const std::string &name, const std::string &value, OptionType type) : name(name), value(value),
                                                                                      type(type) { }
    std::string name;
    std::string value;
    OptionType type;
};

class InnerOptions {
public:
    std::vector<InnerOption *> options;
};

char *findInnerOptions(const char *databaseFilePath);
InnerOptions *parseInnerDatabaseOptions(char *optionsFile);

#endif //INMOST_INNER_PARSER_H_H
