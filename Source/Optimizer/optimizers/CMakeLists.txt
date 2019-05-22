add_subdirectory(noop)
add_subdirectory(bruteforce)
add_subdirectory(alternating)
add_subdirectory(annealing)

if(USE_OPTIMIZER_BAYESIAN)
    add_subdirectory(bayesian)
endif()

set(HEADER ${HEADER} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)