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
