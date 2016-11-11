#ifndef INMOST_UTILS_H
#define INMOST_UTILS_H

#include <sstream>
#include <string>

namespace INMOST {

    template<typename T>
    T from_string(std::string str) {
        T v;
        std::istringstream ss(str);
        ss >> v;
        return v;
    }

}

#endif //INMOST_UTILS_H
