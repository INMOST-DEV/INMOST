#include "utils.h"
#include <cstdlib>

// temporary fix for GeRa

bool IsMaster() {
    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank == 0;
}

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
        for (int    i     = 0; i < lower.length(); i++)
            lower[i] = (char) tolower(lower[i]);
        return lower;
    }

    void MPIBarrier() {
#if defined(USE_MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    int MPIGetRank() {
        int rank = 0;
#if defined(USE_MPI)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        return rank;
    }

}
