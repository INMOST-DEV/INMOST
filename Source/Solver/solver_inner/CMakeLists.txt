set(HEADER ${HEADER} ${CMAKE_CURRENT_SOURCE_DIR}/solver_prototypes.hpp
					 ${CMAKE_CURRENT_SOURCE_DIR}/solver_bcgsl.hpp)

add_subdirectory(solver_ilu2)

set(SOURCE ${SOURCE} PARENT_SCOPE)
set(HEADER ${HEADER} PARENT_SCOPE)