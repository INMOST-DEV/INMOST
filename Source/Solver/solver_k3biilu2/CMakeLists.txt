if (USE_MPI)
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.h" AND EXISTS
            "${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.cpp" AND EXISTS
            "${CMAKE_CURRENT_SOURCE_DIR}/k3d.h" AND EXISTS
            "${CMAKE_CURRENT_SOURCE_DIR}/k3d.cpp")
        set(SOLVER_DEFINITIONS ${SOLVER_DEFINITIONS} -DHAVE_SOLVER_K3BIILU2 PARENT_SCOPE)
        set(HAVE_SOLVER_K3BIILU2 TRUE PARENT_SCOPE)
        set(HEADER ${HEADER} 
            ${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.h 
            ${CMAKE_CURRENT_SOURCE_DIR}/k3d.h
            ${CMAKE_CURRENT_SOURCE_DIR}/SolverK3BIILU2.h PARENT_SCOPE)

        set(SOURCE ${SOURCE} 
            ${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.cpp 
            ${CMAKE_CURRENT_SOURCE_DIR}/k3d.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/SolverK3BIILU2.cpp PARENT_SCOPE)
    else ()
        set(HAVE_SOLVER_K3BIILU2 FALSE PARENT_SCOPE)
    endif ()
endif()