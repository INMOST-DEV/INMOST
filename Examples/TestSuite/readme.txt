The executables here generate data for diffusion and advection tests.

circular_perm  - Populates a permeability tensor that runs the flow in a circular fashion.
                 The new tensor is named PERM and has 6 entries per element.
                 If there was any data named PERM on the mesh it will be deleted.


channels_perm     - Populates a permeability tensor in pattern forming narrow channels. 5 channels.

locking_perm      - Permeability from parametric locking test problem. Populates permeability, 
                    reference solution and reference flux.
                    Run with Models\onephase.
                    With third parameter 0 reproduces problem as considered in:
                    "Mesh locking effects in the finite volume solution of 2-D anisotropic 
                    diffusion equations" Gianmarco Manzini, Mario Putti
                    With third parameter 1 reproduces problems:
		    test5 in "A LINEARITY PRESERVING CELL-CENTERED SCHEME FOR THE HETEROGENEOUS 
                    AND ANISOTROPIC DIFFUSION EQUATIONS ON GENERAL MESHES" Gao Zhiming, Wu Jiming
                    test2 in "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert

discontinuous_4zones - Problem with discontinuous full anisotropic tensor in 4 zones and continuous 
                       linear solution.


fvca5_test1_1     - test1.1 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert
                    test1 in
		    "A LINEARITY PRESERVING CELL-CENTERED SCHEME FOR THE HETEROGENEOUS 
                    AND ANISOTROPIC DIFFUSION EQUATIONS ON GENERAL MESHES" Gao Zhiming, Wu Jiming


fvca5_test1_2     - test1.2 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert

fvca5_test3       - Monotonicity test. Solution should be in [0,1].
                    test3 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert

fvca5_test4       - Monotonicity test with discontinuous permeability. Solution should be in [0,1].
                    test4 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert

fvca5_test5       - Locking and monotonicity test with full-tensor anisotropic permeability. Solution should be in [0,1].
                    test5 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert

fvca5_test6       - Analytic test.
                    test6 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert

fvca5_test7       - Oblique barrier test with linear solution.
                    test7 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert

fvca5_test8       - Perturbed parallelogram test
                    test8 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert
                    also in
                    "COMPACT MULTIPOINT FLUX APPROXIMATION METHOD" paper 
                    by Aavatsmark, Eigestad, Mallison, Nordbotten.

fvca5_test9       - Two wells test
                    test9 in
                    "Benchmark on Discretization Schemes for Anisotropic Diffusion Problems 
                    on General Grids" Raphaele Herbin and Florence Hubert
                    also in
                    "COMPACT MULTIPOINT FLUX APPROXIMATION METHOD" paper 
                    by Aavatsmark, Eigestad, Mallison, Nordbotten.

wugao_test2       - test2 in
		    "A LINEARITY PRESERVING CELL-CENTERED SCHEME FOR THE HETEROGENEOUS 
                    AND ANISOTROPIC DIFFUSION EQUATIONS ON GENERAL MESHES" Gao Zhiming, Wu Jiming

wugao_test3       - Test with unsymmetric permeability tensor.
                    test3 in
		    "A LINEARITY PRESERVING CELL-CENTERED SCHEME FOR THE HETEROGENEOUS 
                    AND ANISOTROPIC DIFFUSION EQUATIONS ON GENERAL MESHES" Gao Zhiming, Wu Jiming

wugao_test4       - Good for numerical locking.
                    test4 in
		    "A LINEARITY PRESERVING CELL-CENTERED SCHEME FOR THE HETEROGENEOUS 
                    AND ANISOTROPIC DIFFUSION EQUATIONS ON GENERAL MESHES" Gao Zhiming, Wu Jiming

adv_test00        - Steady-state test for transverse diffusion from 
                    "On spurious oscillations at layers diminishing (SOLD) methods
                    for convection-diffusion equations: Part II - Analysis for P1 and 
                    Q1 finite elements" by Volker John and Petr Knobloch
                    better description in
                    

adv_test01        - Test for shape preservation (Zalesak disc) from
		    "Fully Multidimensional Flux-Corrected Transport Algorithms for Fluids"
                    by Steven T. Zalesak

adv_test02        - Test for shape preservation (Enright's sphere) from
                    "Use of the particle level set method for enhanced resolution of 
                    free surface flows" by Douglas Patrick Enright

adv_test03        - Test for advection with discontinuous locally degenerate diffusion from
                    "Discontinuous Galerkin Methods for Anisotropic Semidefinite Diffusion with Advection"
                    by Di Pietro, Ern, Guermond, section 6.1

fvca5_test7_adv   - Oblique barrier test with advection, solution is linear. Similar to fvca5_test7.

wugao_test3_adv   - Test with unsymmetric diffusion tensor and advection. Similar to wugao_test3.

two_wells_3d      - similar to fvca5_test9 but in 3d with vertical wells and full tensor anisotropy.

single_well       - in 3d with single vertical well and full tensor anisotropy.

edwards_3dtest4   - case 4 in 3d numerical tests from “q-Families of 
                    CVD(MPFA) Schemes on General Elements:Numerical 
                    Convergence and the Maximum Principle” 
                    by Mayur Pal and Michael Edwards
