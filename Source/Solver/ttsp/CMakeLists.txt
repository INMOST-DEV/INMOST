set(SOURCE ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/ttsp.cpp)

set(HEADER ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/ttsp.h)

add_subdirectory(optimizers)

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)