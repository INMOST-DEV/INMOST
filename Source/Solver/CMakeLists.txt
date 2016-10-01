set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/Solver.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/OrderInfo.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverMaster.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SolverParameters.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/sparse.cpp
        )

add_subdirectory(solver_inner)

if (USE_SOLVER_PETSC)
    add_subdirectory(solver_petsc)
endif()

if (USE_SOLVER_TRILINOS)
    add_subdirectory(solver_trilinos)
endif()

if (USE_SOLVER_ANI)
    add_subdirectory(solver_ani)
endif()

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/solver_fcbiilu2.h" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/solver_fcbiilu2.cpp" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fcbiilu2.cpp")
    set(SOLVER_DEFINITIONS ${SOLVER_DEFINITIONS} -DHAVE_SOLVER_FCBIILU2)
    set(HAVE_SOLVER_FCBIILU2 TRUE PARENT_SCOPE)
    list(APPEND HEADER ${CMAKE_CURRENT_SOURCE_DIR}/solver_fcbiilu2.h)
    list(APPEND SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/solver_fcbiilu2.cpp ${CMAKE_CURRENT_SOURCE_DIR}/fcbiilu2.cpp)
else ()
    set(HAVE_SOLVER_FCBIILU2 FALSE)
endif ()
if (USE_MPI)
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.h" AND EXISTS
            "${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.cpp" AND EXISTS
            "${CMAKE_CURRENT_SOURCE_DIR}/k3d.h" AND EXISTS
            "${CMAKE_CURRENT_SOURCE_DIR}/k3d.cpp")
        set(SOLVER_DEFINITIONS ${SOLVER_DEFINITIONS} -DHAVE_SOLVER_K3BIILU2)
        set(HAVE_SOLVER_K3BIILU2 TRUE PARENT_SCOPE)
        list(APPEND HEADER ${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.h ${CMAKE_CURRENT_SOURCE_DIR}/k3d.h)
        list(APPEND SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.cpp ${CMAKE_CURRENT_SOURCE_DIR}/k3d.cpp)
    else ()
        set(HAVE_SOLVER_K3BIILU2 FALSE PARENT_SCOPE)
    endif ()
endif ()


set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)
set(SOLVER_DEFINITIONS ${SOLVER_DEFINITIONS} PARENT_SCOPE)