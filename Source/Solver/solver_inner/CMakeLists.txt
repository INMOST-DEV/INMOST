set(HEADER ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/solver_prototypes.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/solver_bcgsl.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverInner.h)

set(SOURCE ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverInner.cpp)

add_subdirectory(solver_ilu2)
add_subdirectory(solver_ddpqiluc2)
add_subdirectory(solver_mptiluc)
add_subdirectory(solver_mptilu2)

set(SOURCE ${SOURCE} PARENT_SCOPE)
set(HEADER ${HEADER} PARENT_SCOPE)