if(USE_SOLVER_K3BIILU2)
if (USE_FP64)
if (USE_MPI)
    if (USE_SOLVER_METIS)
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.h" AND EXISTS
                "${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.cpp" AND EXISTS
                "${CMAKE_CURRENT_SOURCE_DIR}/k3d_base.hxx" AND EXISTS
                "${CMAKE_CURRENT_SOURCE_DIR}/k3d_base.cxx" AND USE_SOLVER_METIS)
            set(SOLVER_DEFINITIONS ${SOLVER_DEFINITIONS} -DHAVE_SOLVER_K3BIILU2 PARENT_SCOPE)
            set(HAVE_SOLVER_K3BIILU2 TRUE PARENT_SCOPE)
            set(HEADER ${HEADER}
                    ${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/k3d_base.hxx
                    ${CMAKE_CURRENT_SOURCE_DIR}/k3d_slv.hxx
                    ${CMAKE_CURRENT_SOURCE_DIR}/SolverK3BIILU2.h PARENT_SCOPE)

            set(SOURCE ${SOURCE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/solver_k3biilu2.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/k3d_base.cxx
                    ${CMAKE_CURRENT_SOURCE_DIR}/k3d_slv.cxx
                    ${CMAKE_CURRENT_SOURCE_DIR}/SolverK3BIILU2.cpp PARENT_SCOPE)
    endif()
    else ()
        set(HAVE_SOLVER_K3BIILU2 FALSE PARENT_SCOPE)
        message ("Flag USE_SOLVER_METIS=ON is required for solver K3BIILU")
    endif ()
else ()
	set(HAVE_SOLVER_K3BIILU2 FALSE PARENT_SCOPE)
	message ("Flag USE_MPI=ON is required for solver K3BIILU")
endif ()
else()
	set(HAVE_SOLVER_K3BIILU2 FALSE PARENT_SCOPE)
	message ("Flag USE_FP64=ON is required for solver K3BIILU")
endif()
endif()
