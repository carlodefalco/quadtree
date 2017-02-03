clear all;
close all;
clc;

addpath(canonicalize_file_name("../"));

tol_max = 1e-3;
Nelems_max = 10000;

for i = 1 : 10
    fprintf("i = %d\n", i);
    
    n(i) = 10 * 2^(i-1) + 1;
    
    x = linspace(0, 1, n(i));
    y = linspace(0, 1, n(i));

    msh = msh2m_quadtree(x, y);
    
    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);
    
    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    # Define parameters.
    data.eps1 = 0.5;
    data.eps2 = 1;
    
    omega1 = omega1_nodes(msh);
    omega2 = setdiff(1:Nnodes, omega1);
    
    omega2_el = omega2_elems(msh);
    epsilon = diffusion_coefficient(msh, data, omega2_el);
    
    c = -0.4375 * data.eps2 / (0.5 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + data.eps2 * sinh(0.5 / sqrt(data.eps1)));
    d = 1.75 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) / (sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + 2 * data.eps2 * sinh(0.5 / sqrt(data.eps1)));
    
    u_ex = 1 + 2 * c * sinh(x / sqrt(data.eps1));
    u_ex(omega2) = -0.5 * (x(omega2) - 1) .* (x(omega2) + 2 * d);
    
    du_x_ex1 = 2 * c * cosh(x / sqrt(data.eps1)) / sqrt(data.eps1);
    du_x_ex2 = -x - d + 0.5;
    
    du_x_ex = zeros(size(u_ex));
    du_x_ex(omega1) = du_x_ex1(omega1);
    du_x_ex(omega2) = du_x_ex2(omega2);
    
    du_y_ex = zeros(size(du_x_ex));
    
    # Assemble system.
    A = bim2a_quadtree_laplacian(msh, epsilon) + ...
        bim2a_quadtree_reaction(msh, !omega2_el, ones(Nnodes, 1));
    
    f = ones(Nelems, 1);
    f(omega2_el) = data.eps2;
    g = ones(Nnodes, 1);
    rhs = bim2a_quadtree_rhs (msh, f, g);
    
    # Compute solution and error.
    dnodes = msh2m_nodes_on_sides(msh, [2 4]);
    u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);
    
    refineable_elements = find(!any(msh.children));
    
    msh_new = mark_regions(msh);
    
    [msh1, omega1, omega1_el] = msh2m_quadtree_submesh(msh_new, 1);
    [msh2, omega2, omega2_el] = msh2m_quadtree_submesh(msh_new, 2);
    
    err_L2 = bim2c_quadtree_pde_error_L2_node(msh, u, u_ex);
    err_L2_1 = bim2c_quadtree_pde_error_L2_node(msh1, u(omega1), u_ex(omega1));
    err_L2_2 = bim2c_quadtree_pde_error_L2_node(msh2, u(omega2), u_ex(omega1));
    
    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u);
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);
    
    du_edge1 = bim2c_quadtree_pde_edge_gradient(msh1, u(omega1));
    [du_x1, du_y1] = bim2c_quadtree_pde_reconstructed_gradient(msh1, du_edge1);
    
    du_edge2 = bim2c_quadtree_pde_edge_gradient(msh2, u(omega2));
    [du_x2, du_y2] = bim2c_quadtree_pde_reconstructed_gradient(msh2, du_edge2);
    
    err_edge = bim2c_quadtree_pde_error_semiH1_edge(msh, du_edge, du_x_ex, du_y_ex);
    err_edge_norm(i) = norm(err_edge, 2);
    err_node = bim2c_quadtree_pde_error_semiH1_node(msh, du_x, du_y, du_x_ex, du_y_ex);
    err_node_norm(i) = norm(err_node, 2);
    
    err_edge1 = bim2c_quadtree_pde_error_semiH1_edge(msh1, du_edge1, du_x_ex1(omega1), du_y_ex(omega1));
    err_edge_norm1(i) = norm(err_edge1, 2);
    err_node1 = bim2c_quadtree_pde_error_semiH1_node(msh1, du_x1, du_y1, du_x_ex1(omega1), du_y_ex(omega1));
    err_node_norm1(i) = norm(err_node1, 2);
    
    err_edge2 = bim2c_quadtree_pde_error_semiH1_edge(msh2, du_edge2, du_x_ex2(omega2), du_y_ex(omega2));
    err_edge_norm2(i) = norm(err_edge2, 2);
    err_node2 = bim2c_quadtree_pde_error_semiH1_node(msh2, du_x2, du_y2, du_x_ex2(omega2), du_y_ex(omega2));
    err_node_norm2(i) = norm(err_node2, 2);
    
    # Determine elements to be refined.
    to_refine = false(1, Nelems);
    
    estimator = zeros(1, Nelems);
    estimator(omega1_el) = bim2c_quadtree_pde_ZZ_estimator_du(msh1, u(omega1));
    estimator(omega2_el) = bim2c_quadtree_pde_ZZ_estimator_du(msh2, u(omega2));
    tol = mean(estimator);
    
    to_refine(refineable_elements) = (estimator > tol);
    
    # Save solution to file.
    fclose all;
    basename = "./sol_discontinuous_uniform_piecew_estimator/sol";
    filename = sprintf([basename "_%d_omega"], i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh,
                                 {u, "u"; u_ex, "u_ex";
                                  du_x_ex, "du_x_ex"; du_y_ex, "du_y_ex";
                                  du_x, "du_x_node"; du_y, "du_y_node"},
                                 {err_edge.', "err_edge"; err_node.', "err_node"; estimator.', "estimator"}, 1);
    
    filename1 = sprintf([basename "_%d_omega1"], i);
    if (exist([filename1 ".vtu"], "file"))
        delete([filename1 ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename1, msh1,
                                 {u(omega1), "u"; u_ex(omega1), "u_ex";
                                  du_x_ex1(omega1), "du_x_ex"; du_y_ex(omega1), "du_y_ex";
                                  du_x1, "du_x_node"; du_y1, "du_y_node"},
                                 {err_edge1(:), "err_edge"; err_node1(:), "err_node"; estimator(omega1_el)(:), "estimator"}, 1);
    
    filename2 = sprintf([basename "_%d_omega2"], i);
    if (exist([filename2 ".vtu"], "file"))
        delete([filename2 ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename2, msh2,
                                 {u(omega2), "u"; u_ex(omega2), "u_ex";
                                  du_x_ex2(omega2), "du_x_ex"; du_y_ex(omega2), "du_y_ex";
                                  du_x2, "du_x_node"; du_y2, "du_y_node"},
                                 {err_edge2(:), "err_edge"; err_node2(:), "err_node"; estimator(omega2_el)(:), "estimator"}, 1);
    
    filename3 = sprintf([basename "_%d_edge"], i);
    if (exist([filename3 ".vtu"], "file"))
        delete([filename3 ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh_nedelec(filename3, msh, {du_edge, "du_edge"}, 1);
    
    save("-binary", sprintf([basename "_%d_data.mat"], i), "*");
    
    n_dofs(i) = sum(!any(msh.hanging));
    n_elems(i) = numel(refineable_elements);
    n_to_refine(i) = sum(to_refine);
    global_estimator(i) = norm(estimator, 2);
    global_error(i) = norm(err_L2, 2);
    
    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", ...
         "err_node_norm1", "err_node_norm2", ...
         "err_edge_norm1", "err_edge_norm2", ...
         "global_estimator", "global_error");
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
endfor
