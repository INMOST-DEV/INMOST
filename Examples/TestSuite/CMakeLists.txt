project(TestSuite)

add_executable(single_well single_well.cpp)
add_executable(circular_perm circular_perm.cpp)
add_executable(channels_perm channels_perm.cpp)
add_executable(locking_perm locking_perm.cpp)
add_executable(fvca5_test1_1 fvca5_test1_1.cpp)
add_executable(fvca5_test1_2 fvca5_test1_2.cpp)
add_executable(fvca5_test3 fvca5_test3.cpp)
add_executable(fvca5_test4 fvca5_test4.cpp)
add_executable(fvca5_test5 fvca5_test5.cpp)
add_executable(fvca5_test6 fvca5_test6.cpp)
add_executable(fvca5_test7 fvca5_test7.cpp)
add_executable(fvca5_test8 fvca5_test8.cpp)
add_executable(fvca5_test9 fvca5_test9.cpp)
add_executable(wugao_test2 wugao_test2.cpp)
add_executable(wugao_test3 wugao_test3_unsymmetric.cpp)
add_executable(wugao_test4 wugao_test4.cpp)
add_executable(adv_test00 adv00.cpp)
add_executable(adv_test01 adv01.cpp)
add_executable(adv_test02 adv02.cpp)
add_executable(adv_test03 adv03.cpp)
add_executable(fvca5_test7_adv fvca5_test7_with_advection.cpp)
add_executable(wugao_test3_adv wugao_test3_unsymmetric_with_advection.cpp)
add_executable(transform_grid transform_grid.cpp)
add_executable(two_wells_3d two_wells_3d.cpp)
add_executable(edwards_3dtest4 edwards_3dtest4.cpp)
add_executable(discontinuous_4zones discontinuous_4zones.cpp)
add_executable(linear_sol linear_sol.cpp)





target_link_libraries(single_well inmost)
if(USE_MPI)
  message("linking single_well with MPI")
  target_link_libraries(single_well ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(single_well PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS single_well EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(circular_perm inmost)
if(USE_MPI)
  message("linking circular_perm with MPI")
  target_link_libraries(circular_perm ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(circular_perm PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS circular_perm EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(channels_perm inmost)
if(USE_MPI)
  message("linking channels_perm with MPI")
  target_link_libraries(channels_perm ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(channels_perm PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS channels_perm EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(locking_perm inmost)
if(USE_MPI)
  message("linking locking_perm with MPI")
  target_link_libraries(locking_perm ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(locking_perm PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS locking_perm EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test1_1 inmost)
if(USE_MPI)
  message("linking fvca5_test1_1 with MPI")
  target_link_libraries(fvca5_test1_1 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test1_1 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test1_1 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test1_2 inmost)
if(USE_MPI)
  message("linking fvca5_test1_2 with MPI")
  target_link_libraries(fvca5_test1_2 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test1_2 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test1_2 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(fvca5_test3 inmost)
if(USE_MPI)
  message("linking fvca5_test3 with MPI")
  target_link_libraries(fvca5_test3 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test3 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test3 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test4 inmost)
if(USE_MPI)
  message("linking fvca5_test4 with MPI")
  target_link_libraries(fvca5_test4 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test4 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test4 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test5 inmost)
if(USE_MPI)
  message("linking fvca5_test5 with MPI")
  target_link_libraries(fvca5_test5 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test5 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test5 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test6 inmost)
if(USE_MPI)
  message("linking fvca5_test6 with MPI")
  target_link_libraries(fvca5_test6 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test6 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test6 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test7 inmost)
if(USE_MPI)
  message("linking fvca5_test7 with MPI")
  target_link_libraries(fvca5_test7 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test7 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test7 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test7_adv inmost)
if(USE_MPI)
  message("linking fvca5_test7_adv with MPI")
  target_link_libraries(fvca5_test7_adv ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test7_adv PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test7_adv EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test8 inmost)
if(USE_MPI)
  message("linking fvca5_test8 with MPI")
  target_link_libraries(fvca5_test8 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test8 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test8 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(fvca5_test9 inmost)
if(USE_MPI)
  message("linking fvca5_test9 with MPI")
  target_link_libraries(fvca5_test9 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(fvca5_test9 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS fvca5_test9 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(wugao_test2 inmost)
if(USE_MPI)
  message("linking wugao_test2 with MPI")
  target_link_libraries(wugao_test2 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(wugao_test2 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS wugao_test2 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(wugao_test3 inmost)
if(USE_MPI)
  message("linking wugao_test3 with MPI")
  target_link_libraries(wugao_test3 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(wugao_test3 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS wugao_test3 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(wugao_test3_adv inmost)
if(USE_MPI)
  message("linking wugao_test3_adv with MPI")
  target_link_libraries(wugao_test3_adv ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(wugao_test3_adv PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS wugao_test3_adv EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(wugao_test4 inmost)
if(USE_MPI)
  message("linking wugao_test4 with MPI")
  target_link_libraries(wugao_test4 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(wugao_test4 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS wugao_test4 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(adv_test00 inmost)
if(USE_MPI)
  message("linking adv_test00 with MPI")
  target_link_libraries(adv_test00 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(adv_test00 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS adv_test00 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(adv_test01 inmost)
if(USE_MPI)
  message("linking adv_test01 with MPI")
  target_link_libraries(adv_test01 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(adv_test01 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS adv_test01 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(adv_test02 inmost)
if(USE_MPI)
  message("linking adv_test02 with MPI")
  target_link_libraries(adv_test02 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(adv_test02 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS adv_test02 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(adv_test03 inmost)
if(USE_MPI)
  message("linking adv_test03 with MPI")
  target_link_libraries(adv_test03 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(adv_test03 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS adv_test03 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)



target_link_libraries(transform_grid inmost)
if(USE_MPI)
  message("linking transform_grid with MPI")
  target_link_libraries(transform_grid ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(transform_grid PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS transform_grid EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(two_wells_3d inmost)
if(USE_MPI)
  message("linking two_wells_3d with MPI")
  target_link_libraries(two_wells_3d ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(two_wells_3d PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS two_wells_3d EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(edwards_3dtest4 inmost)
if(USE_MPI)
  message("linking edwards_3dtest4 with MPI")
  target_link_libraries(edwards_3dtest4 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(edwards_3dtest4 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS edwards_3dtest4 EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

target_link_libraries(discontinuous_4zones inmost)
if(USE_MPI)
  message("linking discontinuous_4zones with MPI")
  target_link_libraries(discontinuous_4zones ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(discontinuous_4zones PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS discontinuous_4zones EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)


target_link_libraries(linear_sol inmost)
if(USE_MPI)
  message("linking linear_sol with MPI")
  target_link_libraries(linear_sol ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(linear_sol PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
install(TARGETS linear_sol EXPORT projects-targets RUNTIME DESTINATION bin/TestSuite)

install(FILES readme.txt DESTINATION bin/TestSuite)
