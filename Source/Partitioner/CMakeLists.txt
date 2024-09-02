
if( USE_PARTITIONER_PARMETIS OR USE_PARTITIONER_METIS )
  set(SOURCE
    ${SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/partitioner.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AdaptiveRepartitioner.cpp
    PARENT_SCOPE
  )
else()
  set(SOURCE
    ${SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/partitioner.cpp
    PARENT_SCOPE
  )
endif()

set(HEADER
    ${HEADER}
    PARENT_SCOPE
)
