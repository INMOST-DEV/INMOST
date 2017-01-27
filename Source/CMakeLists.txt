
add_subdirectory(Data)
add_subdirectory(Mesh)
add_subdirectory(IO)
add_subdirectory(Solver)
add_subdirectory(NonlinearSolver)
add_subdirectory(Autodiff)
add_subdirectory(Partitioner)
add_subdirectory(Misc)
add_subdirectory(Headers)


set(SOLVER_DEFINITIONS ${SOLVER_DEFINITIONS} PARENT_SCOPE)
set(SOURCE ${SOURCE} PARENT_SCOPE)
set(HEADER ${HEADER} PARENT_SCOPE)
