set(SOURCE ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/optimization_parameters.cpp)

set(HEADER ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer_interface.h
        ${CMAKE_CURRENT_SOURCE_DIR}/optimization_parameters.h)

add_subdirectory(optimizers)

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)