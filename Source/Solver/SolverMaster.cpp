#include "inmost.h"


namespace INMOST {

    std::map<std::string, SolverBaseFactory *> SolverMaster::solvers = std::map<std::string, SolverBaseFactory *>();

    SolverInterface *SolverMaster::getSolver(std::string name) {
        solvers_map_iterator_t iterator = SolverMaster::solvers.find(name);
        if (iterator != SolverMaster::solvers.end()) {
            return iterator->second->create();
        } else {
            throw INMOST::SolverNotFound;
        }
    }

    SolverInterface *SolverMaster::copySolver(const SolverInterface *other) {
        solvers_map_iterator_t iterator = SolverMaster::solvers.find(other->SolverName());
        if (iterator != SolverMaster::solvers.end()) {
            return iterator->second->copy(other);
        } else {
            throw INMOST::SolverNotFound;
        }
    }

    std::vector<std::string> SolverMaster::getAvailableSolvers() {
        std::vector<std::string> s;
        for (solvers_map_iterator_t bi = SolverMaster::solvers.begin(); bi != SolverMaster::solvers.end(); bi++) {
            s.push_back(bi->first);
        }
        return s;
    }

    bool SolverMaster::isSolverAvailable(std::string name) {
        return SolverMaster::solvers.find(name) != SolverMaster::solvers.end();
    }


}

