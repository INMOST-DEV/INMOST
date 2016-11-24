#include "Utils.h"

namespace INMOST {

    template<>
    int from_string<int>(std::string str) {
        return atoi(str.c_str());
    }

    template<>
    double from_string<double>(std::string str) {
        return atof(str.c_str());
    }

    template<>
    std::string from_string<std::string>(std::string str) {
        return str;
    }

    std::string string_to_lower(const std::string &str) {
        std::string lower = std::string(str);
        for (int i = 0; i < lower.length(); i++) {
            lower[i] = (char) tolower(lower[i]);
        }
        return lower;
    }

}