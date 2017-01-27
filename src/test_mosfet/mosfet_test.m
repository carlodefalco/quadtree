addpath(canonicalize_file_name("../"));

x = linspace (0, 3e-6, 20);
y = linspace (-8e-6, 2e-6, 20);
msh = msh2m_quadtree (x, y, 1, 1:4);

for i = 1 : 5
    # Semiconductor nodes.
    sc_nodes = semic_nodes(msh);
    
    # Insulator elements.
    ins_elems = insulator_elems(msh);
    
    [Nd, Na] = mosfet_doping (msh);
    Nd(!sc_nodes) = Na(!sc_nodes) = 0;
    
    marker = doping_driven_refinement (msh, Nd - Na, 0.1);
    marker(ins_elems) = 0;
    
    fclose all;
    filename = sprintf("sol_mosfet_%d", i);
    delete([filename ".vtu"]);
    fpl_vtk_write_field_quadmesh (filename, msh,
                                  {Na, "Na"; Nd, "Nd"; (Nd - Na), "D"}, {}, 1);

    msh = msh2m_quadtree_refine(msh, find(marker));
    
    fprintf("%d\n", sum(marker));
end
