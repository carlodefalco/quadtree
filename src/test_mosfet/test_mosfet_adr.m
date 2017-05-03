clear all;
close all;
clc;

addpath(canonicalize_file_name("../"));

l = 0;
r = 3e-6;
L =  r - l;

x = unique([linspace(l, l + L/4, 3), ...
            linspace(l + L/4, r - L/4, 4), ...
            linspace(r - L/4, r, 3)]);

y = union(linspace(-8e-6, 0, 10),
          linspace(0, 2e-6, 6));

msh = msh2m_quadtree(x, y, 1, 1:4);
msh.l = l;
msh.r = r;
msh.L = L;

msh = mark_regions(msh);
msh = msh2m_quadtree_submesh(msh, 1:2);

for i = 1 : 4
    msh.l = l;
    msh.r = r;
    msh.L = L;
    
    # Semiconductor nodes.
    sc_nodes = semic_nodes(msh);
    
    # Insulator elements.
    ins_elems = insulator_elems(msh);
    
    [Nd, Na] = mosfet_doping (msh);
    Nd(!sc_nodes) = Na(!sc_nodes) = 0;
    
    marker = doping_driven_refinement (msh, Nd - Na, 0.1);
    marker(ins_elems) = 0;
    
    filename = sprintf ("sol_mosfet_%d", i);
    delete ([filename ".vtu"]);
    fpl_vtk_write_field_quadmesh (filename, msh,
                                  {Na, "Na"; Nd, "Nd"; (Nd - Na), "D"}, {}, 1);

    msh = msh2m_quadtree_refine(msh, find(marker));
end

fprintf("Mesh building complete!\n\n");
