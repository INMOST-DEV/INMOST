//
// Created by Dmitri Bagaev on 28.09.16.
//

#include "inmost.h"


namespace INMOST {

    std::map<std::string, SolverBaseFactory *> SolverFactory::solvers = std::map<std::string, SolverBaseFactory *>();

    SolverInterface *SolverFactory::getSolver(std::string name) {
        auto iterator = SolverFactory::solvers.find(name);
        if (iterator != SolverFactory::solvers.end()) {
            return iterator->second->create();
        } else {
            throw INMOST::SolverNotFound;
        }
    }

    SolverInterface *SolverFactory::copySolver(const SolverInterface *other) {
        auto iterator = SolverFactory::solvers.find(other->SolverName());
        if (iterator != SolverFactory::solvers.end()) {
            return iterator->second->copy(other);
        } else {
            throw INMOST::SolverNotFound;
        }
    }

    std::vector<std::string> SolverFactory::getAvailableSolvers() {
        std::vector<std::string> s;
        for (auto bi = SolverFactory::solvers.begin(); bi != SolverFactory::solvers.end(); bi++) {
            s.push_back(bi->first);
        }
        return s;
    }

    bool SolverFactory::isSolverAvailable(std::string name) {
        return SolverFactory::solvers.find(name) != SolverFactory::solvers.end();
    }


}

