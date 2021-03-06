clear all;
close all;
clc;

n = 10;

# Mesh definition.
x = linspace(0, 1, n);
y = linspace(0, 1, n);

region = 1;
sides = 1:4;

msh = msh2m_quadtree(x, y, region, sides);
#msh = msh2m_quadtree_refine(msh, [22 23 31 32]);

x = msh.p(1, :).';
y = msh.p(2, :).';

u_ex = (x - 1/2).^2 .* (y - 1/2).^2;
du_x_ex = 2 * (x - 1/2) .* (y - 1/2).^2;
du_y_ex = 2 * (x - 1/2).^2 .* (y - 1/2);

du = bim2c_quadtree_pde_edge_gradient(msh, u_ex);
[du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du);

#boundary = msh2m_nodes_on_sides(msh, 1:4);
#du_x(boundary) = du_x_ex(boundary);
#du_y(boundary) = du_y_ex(boundary);

# Save solution to file.
fclose all;
filename = sprintf("./sol_%d", 1);
if (exist([filename ".vtu"], "file"))
    delete([filename ".vtu"]);
endif
fpl_vtk_write_field_quadmesh(filename, msh, {du_x_ex, "du_x_ex";
                                             du_y_ex, "du_y_ex";
                                             du_x, "du_x";
                                             du_y, "du_y"}, {}, 1);
