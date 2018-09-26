SET GRIDGEN=..\build64_2010_INMOST\Examples\GridGen\Release\GridGen.exe
SET GENERATORS_FOLDER= ..\build64_2010_DiscretizationToolkit\Grids\generators\Release

:: quad grid
 %GRIDGEN% 4 %1 %1 1 grid.pmf
:: %GRIDGEN% 4 80 81 1 grid.pmf

:: triangular grid
:: %GRIDGEN% 3 %1 %1 1 grid.pmf

:: hex grid
:: %GENERATORS_FOLDER%\hex_grid.exe %1

:: acute grid
:: %GENERATORS_FOLDER%\acute_grid.exe %1

:: nonconvex grid
:: %GENERATORS_FOLDER%\nonconvex_grid.exe %1

:: (advection only) discontinuous dirichlet condition problem
:: %GENERATORS_FOLDER%\adv_test00.exe grid.pmf grid_out.pmf

:: (advection only) zalesak disc rotation
:: %GENERATORS_FOLDER%\adv_test01.exe grid.pmf grid_out.pmf 4

:: (advection only) enright disc
:: %GENERATORS_FOLDER%\adv_test02.exe grid.pmf grid_out.pmf 4

:: (diffusion only) oblique barrier problem, linear solution
:: %GENERATORS_FOLDER%\fvca5_test7.exe grid.pmf grid_out.pmf

:: (advection-diffusion) oblique barrier problem, linear solution
:: %GENERATORS_FOLDER%\fvca5_test7_adv.exe grid.pmf grid_out.pmf

:: (diffusion only) mild anisotropy test2
:: %GENERATORS_FOLDER%\fvca5_test1_2.exe grid.pmf grid_out.pmf

:: (advection-reaction-diffusion) discontinuous
%GENERATORS_FOLDER%\adv_test03.exe grid.pmf grid_out.pmf

:: (diffusion with unsymmetric tensor)
:: %GENERATORS_FOLDER%\wugao_test3.exe grid.pmf grid_out.pmf

:: (advection-diffusion with unsymmetric tensor)
:: %GENERATORS_FOLDER%\wugao_test3_adv.exe grid.pmf grid_out.pmf

:: (dmp test)
:: %GENERATORS_FOLDER%\dmp_grid.exe .\grid.pmf 0.6 0 1 1 1 1000

::mpiexec -np 4 .\Release\NFVADV.exe .\grid_out.pmf %2 
.\Release\NFVADV.exe .\grid_out.pmf %2
