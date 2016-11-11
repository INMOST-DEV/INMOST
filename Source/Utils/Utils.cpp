#include "inmost_utils.h"

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

}