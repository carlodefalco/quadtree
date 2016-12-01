clear all;
close all;
clc;

# Mesh definition.
x = linspace(0, 5, 10);
y = linspace(0, 3,  5);

region = 1;
sides = 1:4;

msh = msh2m_quadtree(x, y, region, sides);

# Mesh refinement.
msh = msh2m_refine_quadtree(msh, [4, 8]);
msh = msh2m_refine_quadtree(msh, [40 48 52]);
msh = msh2m_refine_quadtree(msh, [44 12 43 59]);
msh = msh2m_refine_quadtree(msh, [33 74 78]);

# Show refined mesh.
figure; quadmesh(msh, "show_cell_numbers", "show_node_numbers");
