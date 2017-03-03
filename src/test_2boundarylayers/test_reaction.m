clear all;
close all;
clc;

addpath(canonicalize_file_name("../"));

n = 6;

# Mesh definition.
x = linspace(0, 1, n);
y = linspace(0, 1, n);

msh = msh2m_quadtree(x, y);

tol_max = 1e-3;
Nelems_max = 10000;

for i = 1 : 10
    fprintf("i = %d\n", i);
    msh = bim2c_quadtree_mesh_properties(msh);

    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);

    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    # Define parameters and exact solution.
    epsilon = 1e-2;
    
    u_ex = (1 - sinh(x / sqrt(epsilon)) / sinh(1 / sqrt(epsilon))) .* (1 - sinh(y / sqrt(epsilon)) / sinh(1 / sqrt(epsilon)));
    du_x_ex = (cosh(x / sqrt(epsilon)) / sinh(1 / sqrt(epsilon)) .* (-1 + sinh(y / sqrt(epsilon)) / sinh(1 / sqrt(epsilon)))) / sqrt(epsilon);
    du_y_ex = (cosh(y / sqrt(epsilon)) / sinh(1 / sqrt(epsilon)) .* (-1 + sinh(x / sqrt(epsilon)) / sinh(1 / sqrt(epsilon)))) / sqrt(epsilon);
    
    # Assemble system.
    alpha = @(msh) epsilon * ones(columns(msh.t), 1);
    delta = @(msh) ones(columns(msh.t), 1);
    zeta  = @(msh) ones(columns(msh.p), 1);
    
    A = bim2a_quadtree_laplacian(msh, alpha(msh)) + bim2a_quadtree_reaction(msh, delta(msh), zeta(msh));

    dnodes = msh2m_nodes_on_sides(msh, 1:4);

    f = @(msh) ones(columns(msh.t), 1);
    g = @(msh) 1 - sinh(msh.p(1, :).' / sqrt(epsilon)) .* sinh(msh.p(2, :).' / sqrt(epsilon)) / sinh(1 / sqrt(epsilon))^2;
    rhs = bim2a_quadtree_rhs(msh, f(msh), g(msh));

    u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);
    
    err_L2 = bim2c_quadtree_pde_error_L2_node(msh, u, u_ex);
    
    refineable_elements = find(!any(msh.children));
    
    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u);
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);

    err_edge = bim2c_quadtree_pde_error_semiH1_edge(msh, du_edge, du_x_ex, du_y_ex);
    err_edge_norm(i) = norm(err_edge, 2);
    
    err_node = bim2c_quadtree_pde_error_semiH1_node(msh, du_x, du_y, du_x_ex, du_y_ex);
    err_node_norm(i) = norm(err_node, 2);
    
    # Determine elements to be refined.
    estimator = bim2c_quadtree_pde_ZZ_estimator_du(msh, u);
    to_refine = bim2c_quadtree_pde_ZZ_to_refine(msh, estimator, mean(estimator), 1);
    
    # Save solution to file.
    fclose all;
    basename = "./sol_reaction/sol";
    filename = sprintf([basename "_%d"], i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex"},
                                 {err_edge.', "err_edge"; err_node.', "err_node"; estimator.', "estimator"}, 1);
    
    n_dofs(i) = sum(!any(msh.hanging));
    n_elems(i) = numel(refineable_elements);
    n_to_refine(i) = sum(to_refine);
    global_estimator(i) = norm(estimator, 2);
    global_error(i) = norm(u - u_ex, inf);
    
    save("-binary", sprintf([basename "_%d_data.mat"], i), "*");
    
    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", ...
         "err_edge_norm", "err_node_norm", ...
         "global_estimator", "global_error");
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    # Do refinement.
    if (n_elems(i) >= Nelems_max)
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
