fix_faults - Intersects faces of the mesh along unresolved faults. Introduce a mesh with a set of blocks, 
             shifted against each other and this tool will introduce edges and matching faces along
             faces of shifted blocks. Only for convex mesh elements.

mesh_input - general polyhedral grid with unresolved faults.
mesh_output - output general polyhedral grid with resolved faults, default grid.vtk

fix_tiny - Collapses small cells, faces and edges in the mesh. 
           Cells may collapse into face, edge or node, face into edge or node and edge into node.
           Only for convex mesh elements.
           Warning: cells are not implemented

mesh_input - general polyhedral grid
mesh_output - output general polyhedral grid with no tiny elements, default grid.vtk

unite_faces - Unite faces of the cell that share the same adjacent cell, some edges and are approximately
              coplanar.

mesh_input - general polyhedral grid
mesh_output - output general polyhedral grid with united faces, default grid.vtk


difference_same - Computes difference between the real data on all elements of two similar grids
                  and writes the grid with the difference. It writes absolute value of difference.
                  Outputs L_inf, L_1, L_2 norms of the difference.
                  Similar grids means their grid blocks order should match between each other

Usage:

difference_same mesh_file1 mesh_file2

mesh_file1 - first mesh file
mesh_file2 - second mesh file


difference_map - Computes difference between the real data on cells of two different grids
                 and writes both grids with the difference. It writes absolute value of difference.
                 Outputs L_inf, L_1, L_2 norms of the difference.
                 It will map cells between two grids by checking which cell centers on one grid
                 fit into which cells on the other and vice verse, resulting in bijective mapping.
                 Then for each cell on one grid difference between it's value and mean value
                 of it's mapping is calculated.

Usage:

difference_map mesh_input1 mesh_input2 [mesh_output1] [mesh_output2]

mesh_input1  - first mesh file
mesh_input2  - second mesh file
mesh_output1 - where to write first mesh with the difference, if not specified, writes to mesh1diff.vtk
mesh_output2 - where to write second mesh with the difference, if not specified, writes to mesh2diff.vtk