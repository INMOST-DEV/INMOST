project(TTSPFunction)
set(SOURCE main.cpp)

add_subdirectory(function)
include_directories(function)

add_executable(TTSPFunction ${SOURCE})

target_link_libraries(TTSPFunction inmost)

if(USE_MPI)
    message("linking TTSPFunction with MPI")
    target_link_libraries(TTSPFunction ${MPI_LIBRARIES})
    if(MPI_LINK_FLAGS)
        set_target_properties(TTSPFunction PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif(USE_MPI)

install(TARGETS TTSPFunction EXPORT inmost-targets RUNTIME DESTINATION bin)
