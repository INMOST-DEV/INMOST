set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/solver.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/OrderInfo.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverMaster.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverParameters.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sparse.cpp
        )

set(HEADER
        ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverMaster.h
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverInterface.h)

add_subdirectory(solver_inner)
add_subdirectory(solver_fcbiilu2)
add_subdirectory(solver_k3biilu2)

if (USE_SOLVER_PETSC)
    add_subdirectory(solver_petsc)
endif()

if (USE_SOLVER_TRILINOS AND USE_MPI)
    add_subdirectory(solver_trilinos)
endif()

if (USE_SOLVER_ANI AND TEST_FORTRAN_ANI3D)
    add_subdirectory(solver_ani)
endif()

if (USE_SOLVER_SUPERLU)
    add_subdirectory(solver_superlu)
endif()


set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)
set(SOLVER_DEFINITIONS ${SOLVER_DEFINITIONS} PARENT_SCOPE)