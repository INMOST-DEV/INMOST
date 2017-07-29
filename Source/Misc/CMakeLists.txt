set(SOURCE
        ${SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/timer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/xml.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/base64.cpp
    PARENT_SCOPE
)

set(HEADER
        ${HEADER}
        ${CMAKE_CURRENT_SOURCE_DIR}/utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/base64.h
        PARENT_SCOPE
)