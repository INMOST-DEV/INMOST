set(SOURCE ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/static_sin.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/static_sin_r2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/dynamic_x2.cpp)

set(HEADER ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/dynamic_function.h
        ${CMAKE_CURRENT_SOURCE_DIR}/static_sin.h
        ${CMAKE_CURRENT_SOURCE_DIR}/static_sin_r2.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dynamic_x2.h)

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)