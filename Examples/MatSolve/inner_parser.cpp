#include "inner_parser.h"

void removeSpaces(char* source) {
    char* i = source;
    char* j = source;
    while(*j != 0)
    {
        *i = *j++;
        if(*i != ' ')
            i++;
    }
    *i = 0;
}

InnerOptions *parseInnerDatabaseOptions(char *optionsFile) {
    FILE *databaseFile = fopen(optionsFile, "r");
    if (!databaseFile) {
        std::cout << "Inner options file not found" << std::endl;
        return nullptr;
    }
    InnerOptions *options = new InnerOptions();
    char *tmp = (char *) calloc(256, sizeof(char));
    char *parameterName = (char *) calloc(128, sizeof(char));
    char *parameterValue = (char *) calloc(128, sizeof(char));
    char *type = (char *) calloc(36, sizeof(char));
    while (!feof(databaseFile) && fgets(tmp, 256, databaseFile)) {
        //removeSpaces(tmp);
        char *line = tmp;
        //Comment str
        if (line[0] == '#') continue;
        //First 4 chars is 'real' or 'enum'
        bool isReal = false, isEnum = false;
        sscanf(line, "%s %s %s", type, parameterName, parameterValue);
        if (strncmp(type, "real", 4) == 0) {
            //fprintf(stdout, "Real parameter:");
            isReal = true;
        } else if (strncmp(type, "enum", 4) == 0) {
            //fprintf(stdout, "Enum parameter:");
            isEnum = true;
        } else {
            fprintf(stderr, "Skipping bad line: %s", line);
            continue;
        }
        options->options.push_back(new InnerOption(std::string(parameterName), std::string(parameterValue), isReal ? REAL : ENUM));

    }
    free(parameterValue);
    free(parameterName);
    free(tmp);
    free(type);
    return options;
}

char *findInnerOptions(const char *databaseFilePath) {
    FILE *databaseFile = fopen(databaseFilePath, "r");
    if (databaseFile == NULL) {
        fprintf(stderr, "Database file not found\n");
        #if defined(USE_MPI)
        MPI_Finalize();
        #endif
        exit(1);
    }
    char *tmp = (char *) calloc(256, sizeof(char));
    char *parameterName = (char *) calloc(64, sizeof(char));
    char *parameterValue = (char *) calloc(64, sizeof(char));
    char *fileName = NULL;
    while (!feof(databaseFile) && fgets(tmp, 256, databaseFile)) {
        removeSpaces(tmp);
        char *line = tmp;
        if (strncmp(line, "inner", 5) != 0) continue;
        line += 6;
        size_t fileNameLength = 0;
        while (!isspace(*(line + fileNameLength))) fileNameLength += 1;
        fileName = (char *) calloc(fileNameLength + 1, sizeof(char));
        memcpy(fileName, line, fileNameLength);
        fileName[fileNameLength] = 0;
    }
    free(tmp);
    return fileName;
}

