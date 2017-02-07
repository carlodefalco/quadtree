clear all;
close all;
clc;

addpath(canonicalize_file_name("./test_sol/"));

function [P] = projector(msh)
    Nnodes = columns(msh.p);
    
    hanging = find(any(msh.hanging));
    non_hanging = find(!any(msh.hanging));
    
    II = non_hanging;
    JJ = msh.full_to_reduced(non_hanging);
    VV = ones(size(II));
    
    P = sparse(II, JJ, VV, Nnodes, numel(non_hanging));
    
    II = repmat(hanging, [2, 1])(:);
    JJ = msh.full_to_reduced(msh.hanging(:, hanging))(:);
    VV = 0.5 * ones(size(II));
    
    P += sparse(II, JJ, VV, Nnodes, numel(non_hanging));
endfunction

# Mesh definition.
x = linspace(0, 1, 3);
y = linspace(0, 1, 3);

msh = msh2m_quadtree(x, y);
msh = msh2m_quadtree_refine(msh, 1);

msh.P = projector(msh);

i = 1;

# Build global matrix.
Nnodes = columns(msh.p);
Nelements = columns(msh.t);

x = msh.p(1, :).';
y = msh.p(2, :).';

# Define parameters and exact solution.
epsilon = 1;

u_ex = x .* y;
du_x_ex = y;
du_y_ex = x;

# Assemble system.
alpha = epsilon * ones(Nelements, 1);

A = bim2a_quadtree_laplacian(msh, alpha);

dnodes = msh2m_nodes_on_sides(msh, 1:4);

f = zeros(Nelements, 1);
g = ones(Nnodes, 1);
rhs = bim2a_quadtree_rhs(msh, f, g);

# Compute solution and error.
u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);

[du] = bim2c_quadtree_pde_edge_gradient(msh, u);
[du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du);
[estimator] = bim2c_quadtree_pde_ZZ_estimator_du(msh, u);

err(i) = norm(u - u_ex, inf);

# Save solution to file.
fclose all;
filename = sprintf("sol_%d", i);
if (exist([filename ".vtu"], "file"))
    delete([filename ".vtu"]);
endif
fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex";
                                             u - u_ex, "err";
                                             du_x, "du_x"; du_y, "du_y";
                                             du_x_ex, "du_x_ex"; du_y_ex, "du_y_ex";
                                             du_x - du_x_ex, "err_du_x";
                                             du_y - du_y_ex, "err_du_y"}, {estimator.', "estimator"}, 1);

function [phi] = plot_phi5(msh, x, y)
    phi5 = bim2c_quadtree_eval_fun(msh, 2, x, y)(:, 2) ...
           + bim2c_quadtree_eval_fun(msh, 3, x, y)(:, 4) ...
           + bim2c_quadtree_eval_fun(msh, 4, x, y)(:, 1) ...
           + bim2c_quadtree_eval_fun(msh, 7, x, y)(:, 3);
    
    phi11 = bim2c_quadtree_eval_fun(msh, 6, x, y)(:, 3) ...
            + bim2c_quadtree_eval_fun(msh, 7, x, y)(:, 2);
    
    phi12 = bim2c_quadtree_eval_fun(msh, 7, x, y)(:, 4) ...
            + bim2c_quadtree_eval_fun(msh, 8, x, y)(:, 3);
    
    phi = phi5 + 0.5 * phi11 + 0.5 * phi12;
    
    phi(x == 0.5 & y != 0.25) /= 2;
    phi(y == 0.5 & x != 0.25) /= 2;
    
    phi(x == 0.5 & y == 0.25) /= 3;
    phi(y == 0.5 & x == 0.25) /= 3;
    
    phi(x == 0.25 & y < 0.5) /= 2;
    phi(y == 0.25 & x < 0.5) /= 2;
endfunction

x = linspace(0, 1, 101);
y = x;

[X, Y] = meshgrid(x, y);

phi5 = reshape(plot_phi5(msh, X(:), Y(:)), size(X));
surf(X, Y, phi5);
