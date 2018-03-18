clear all;
close all;
clc;

addpath(canonicalize_file_name("../"));

n = 3;

# Mesh definition.
x = linspace(0, 1, n);
y = linspace(0, 1, n);

msh = msh2m_quadtree(x, y);

tol_max = 1e-3;
Nelems_max = 10000;

for i = 1 : 10
    fprintf("i = %d\n", i);
    [msh, edge_space, node_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    
    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);

    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    # Define parameters and exact solution.
    lambda = 10;
    
    u_ex = (exp(lambda * x) - 1) / (exp(lambda) - 1) .* (exp(lambda * y) - 1) / (exp(lambda) - 1);
    du_x_ex = (lambda * exp(lambda * x)) / (exp(lambda) - 1) .* (exp(lambda * y) - 1) / (exp(lambda) - 1);
    du_y_ex = (exp(lambda * x) - 1) / (exp(lambda) - 1) .* (lambda * exp(lambda * y)) / (exp(lambda) - 1);
    
    # Assemble system.
    alpha = @(msh) ones(columns(msh.t), 1);
    psi = @(msh) lambda * (msh.p(1, :) + msh.p(2, :)).';
    
    A = bim2a_quadtree_advection_diffusion(msh, alpha(msh), psi(msh));

    dnodes = msh2m_nodes_on_sides(msh, 1:4);

    f = @(msh) zeros(columns(msh.t), 1);
    g = @(msh) zeros(columns(msh.p), 1);
    rhs = bim2a_quadtree_rhs(msh, f(msh), g(msh));

    u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);
    
    refineable_elements = find(!any(msh.children));
    
    # Compute numerical gradient on edge and node space.
    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u);
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);

    err_edge = bim2c_quadtree_pde_error_semiH1_edge(msh, du_edge, du_x_ex, du_y_ex);
    err_edge_norm(i) = norm(err_edge, 2);
    
    err_node = bim2c_quadtree_pde_error_semiH1_node(msh, du_x, du_y, du_x_ex, du_y_ex);
    err_node_norm(i) = norm(err_node, 2);
    
    # Determine elements to be refined.
    to_refine = false(1, Nelems);
    
    estimator = bim2c_quadtree_pde_ZZ_estimator_du(msh, u);
    tol = mean(estimator);
    
    to_refine(refineable_elements) = true;
    
    # Save solution to file.
    fclose all;
    basename = "./sol_advection_uniform/sol";
    filename = sprintf([basename "_%d"], i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex"; du_x, "du_x"; du_x_ex, "du_x_ex"},
                                 {err_edge.', "err_edge"; err_node.', "err_node"; estimator.', "estimator"}, 1);
    
    n_dofs(i) = sum(!any(msh.hanging));
    n_elems(i) = numel(refineable_elements);
    n_to_refine(i) = sum(to_refine);
    global_estimator(i) = norm(estimator, 2);
    global_error(i) = norm(u - u_ex, inf);
    
    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", ...
         "err_edge_norm", "err_node_norm", ...
         "global_estimator", "global_error");
    
    # Do refinement.
    if (n_elems(i) >= Nelems_max || (tol <= tol_max / n_elems(i)))
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
