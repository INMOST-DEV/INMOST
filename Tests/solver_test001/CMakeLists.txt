project(solver_test001)
set(SOURCE main.cpp)

add_executable(solver_test001 ${SOURCE})
target_link_libraries(solver_test001 inmost)

if(USE_MPI)
  message("linking solver_test001 with MPI")
  target_link_libraries(solver_test001 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(solver_test001 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)



if(USE_SOLVER_ANI)
  message("linking solver_test001 with ani3d")
  target_link_libraries(solver_test001 ani3d)
endif()
if(USE_SOLVER_PETSC)
  message("linking solver_test001 with PETSc")
  target_link_libraries(solver_test001 ${PETSC_LIBRARIES})
endif()
if(USE_SOLVER_TRILINOS)
  message("linking solver_test001 with Trilinos")
  target_link_libraries(solver_test001 ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()
if(USE_SOLVER_METIS)
  message("linking solver_test001 with Metis")
  target_link_libraries(solver_test001 ${METIS_LIBRARIES})
endif()
if(USE_SOLVER_MONDRIAAN)
  message("linking solver_test001 with Mondriaan")
  target_link_libraries(solver_test001 ${MONDRIAAN_LIBRARIES})
endif()
if(USE_SOLVER_SUPERLU)
  message("linking solver_test001 with SuperLU")
  target_link_libraries(solver_test001 ${SUPERLU_LIBRARIES})
endif()

set(TOKAMAK_TESTS utm300
                  utm1700a
                  utm1700b
                  utm3060
                  utm5940)

foreach(test ${TOKAMAK_TESTS})
  add_test(NAME solver_test001_sparsekit_tokamak_${test}  COMMAND $<TARGET_FILE:solver_test001> ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/tokamak/${test}.mtx ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/tokamak/${test}_rhs1.mtx)
endforeach()


set(DRIVCAV_OLD_TESTS cavity01
                      cavity02
                      cavity03
                      cavity04
                      cavity05
                      cavity06
                      cavity07
                      cavity08
                      cavity09
                      cavity10
                      cavity11
                      cavity12
                      cavity13
                      cavity14
                      cavity15
                      cavity16
                      cavity17
                      cavity18
                      cavity19
                      cavity20
                      cavity21
                      cavity22
                      cavity23
                      cavity24
                      cavity25
                      cavity26)

foreach(test ${DRIVCAV_OLD_TESTS})
  add_test(NAME solver_test001_sparsekit_drivcav_old_${test}     COMMAND $<TARGET_FILE:solver_test001> ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/drivcav_old/${test}.mtx ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/drivcav_old/${test}_rhs1.mtx)
endforeach()



set(DRIVCAV_TESTS e05r0000
                  e05r0100
                  e05r0200
                  e05r0300
                  e05r0400
                  e05r0500
                  e20r0000
                  e20r0100
                  e20r0500
                  e20r1000
                  e20r2000
                  e20r3000
                  e20r4000
                  e20r5000
                  e30r0000
                  e30r0100
                  e30r0500
                  e30r1000
                  e30r2000
                  e30r3000
                  e30r4000
                  e30r5000
                  e40r0000
                  e40r0100
                  e40r0500
                  e40r1000
                  e40r2000
                  e40r3000
                  e40r4000
                  e40r5000)


foreach(test ${DRIVCAV_TESTS})
  add_test(NAME solver_test001_sparsekit_drivcav_${test}     COMMAND $<TARGET_FILE:solver_test001> ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/drivcav/${test}.mtx ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/drivcav/${test}_rhs1.mtx)
endforeach()


set(FIDAP_TESTS fidap001
                fidap002
                fidap003
                fidap004
                fidap005
                fidap006
                fidap007
                fidap008
                fidap009
                fidap010
                fidap011
                fidap012
                fidap013
                fidap014
                fidap015
                fidap018
                fidap019
                fidap020
                fidap021
                fidap022
                fidap023
                fidap024
                fidap025
                fidap026
                fidap027
                fidap028
                fidap029
                fidap031
                fidap032
                fidap033
                fidap035
                fidap036
                fidap037
                fidapm02
                fidapm03
                fidapm05
                fidapm07
                fidapm08
                fidapm09
                fidapm10
                fidapm11
                fidapm13
                fidapm15
                fidapm29
                fidapm33
                fidapm37)



foreach(test ${FIDAP_TESTS})
  add_test(NAME solver_test001_sparsekit_fidap_${test}     COMMAND $<TARGET_FILE:solver_test001> ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/fidap/${test}.mtx ${CMAKE_CURRENT_SOURCE_DIR}/matrices/sparskit/fidap/${test}_rhs1.mtx)
endforeach()

