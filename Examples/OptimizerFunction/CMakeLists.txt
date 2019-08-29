project(OptimizerFunction)
set(SOURCE main.cpp)

add_subdirectory(function)
include_directories(function)

add_executable(OptimizerFunction ${SOURCE})

target_link_libraries(OptimizerFunction inmost)

if(USE_MPI)
    message("linking OptimizerFunction with MPI")
    target_link_libraries(OptimizerFunction ${MPI_LIBRARIES})
    if(MPI_LINK_FLAGS)
        set_target_properties(OptimizerFunction PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif(USE_MPI)

install(TARGETS OptimizerFunction EXPORT inmost-targets RUNTIME DESTINATION bin)
