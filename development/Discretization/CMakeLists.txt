project(Discretization)
set(SOURCE discr_common.cpp
			discr_tpfa.cpp
			discr_tpfa2.cpp
			discr_ntpfa_a.cpp
			discr_ntpfa_b.cpp
			discr_ntpfa_c.cpp
			discr_mpfa_a.cpp
			discr_mpfa_b.cpp
			discr_nmpfa_a.cpp
			discr_nmpfa_b.cpp
			discr_mpfa_o.cpp
			discr_mpfa_l.cpp
			discr_mpfa_c.cpp
			discr_triplet_basic.cpp
			discr_basic_a.cpp
			discr.h)

add_library(Discretization ${SOURCE})
if(USE_MPI)
  target_link_libraries(Discretization ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(Discretization PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)