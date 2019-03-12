if (USE_TTSP)
    set(SOURCE ${SOURCE}
            ${CMAKE_CURRENT_SOURCE_DIR}/ttsp.cpp)

    add_subdirectory(optimizers)
endif ()

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)