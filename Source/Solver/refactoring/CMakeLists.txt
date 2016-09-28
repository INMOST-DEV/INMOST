if(USE_SOLVER_PETSC)
    add_subdirectory(solver_petsc)
endif()


set(HEADER
        ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverInterface.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverFactory.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Solver2.h)

set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/Solver2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverFactory.cpp)

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)