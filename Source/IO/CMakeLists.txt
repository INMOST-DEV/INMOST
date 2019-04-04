set(SOURCE
    ${SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_pmf_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_xml_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_vtk_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_vtu_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_pvtk_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_gmsh_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_gmv_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_mkf_file.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_ecl_file.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_openfoam_file.cpp
    PARENT_SCOPE
)

set(HEADER
    ${HEADER}
    ${CMAKE_CURRENT_SOURCE_DIR}/io.hpp
    PARENT_SCOPE
)
