set(SOURCE ${SOURCE} ${CMAKE_CURRENT_SOURCE_DIR}/SolverTrilinos.cpp)
set(HEADER ${HEADER} ${CMAKE_CURRENT_SOURCE_DIR}/SolverTrilinos.h)

add_subdirectory(solver_aztec)
add_subdirectory(solver_belos)
add_subdirectory(solver_ml)
add_subdirectory(solver_ifpack)

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)