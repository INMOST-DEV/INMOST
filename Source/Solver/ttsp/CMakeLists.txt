set(HEADER ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/ttsp_configuration.h)

set(SOURCE ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/ttsp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ttsp_configuration.cpp)

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)