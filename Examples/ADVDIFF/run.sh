GRIDGEN=../INMOST/Examples/GridGen/GridGen
GENERATORS_FOLDER=../DiscretizationToolkit/Grids/generators

# quad grid
#${GRIDGEN} 4 $1 $1 1 grid.pmf

# triangular grid
#${GRIDGEN} 3 $1 $1 1 grid.pmf

# hex grid
${GENERATORS_FOLDER}/hex_grid $1

# acute grid
#${GENERATORS_FOLDER}/acute_grid $1

# nonconvex grid
#${GENERATORS_FOLDER}/nonconvex_grid $1

# discontinuous dirichlet condition problem
#${GENERATORS_FOLDER}/adv_test00 grid.pmf grid_out.pmf

# zalesak disc rotation
${GENERATORS_FOLDER}/adv_test01 grid.pmf grid_out.pmf 4

./NFVADV grid_out.pmf $2

