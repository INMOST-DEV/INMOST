project(xml_reader_test000)
set(SOURCE main.cpp)

add_executable(xml_reader_test000 ${SOURCE})

target_link_libraries(xml_reader_test000 inmost)

if(USE_MPI)
    message("linking xml_reader_test000 with MPI")
    target_link_libraries(xml_reader_test000 ${MPI_LIBRARIES})
    if(MPI_LINK_FLAGS)
        set_target_properties(xml_reader_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif(USE_MPI)

if(USE_ZLIB)
  message("linking xml_reader_test000 with ZLIB")
  target_link_libraries(xml_reader_test000 ${ZLIB_LIBRARIES}) 
endif(USE_ZLIB)


add_test(NAME xml_reader_test000_all COMMAND $<TARGET_FILE:xml_reader_test000>)
