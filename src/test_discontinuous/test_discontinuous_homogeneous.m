clear all;
close all;
clc;

addpath(canonicalize_file_name("../"));

x = linspace(0, 1, 11);
y = linspace(0, 1, 11);

msh = msh2m_quadtree(x, y);

tol_max = 1e-3;
Nelems_max = 10000;

for i = 1 : 10
    fprintf("i = %d\n", i);
    
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
    
    c = -0.5 * data.eps2 / (0.5 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + data.eps2 * sinh(0.5 / sqrt(data.eps1)));
    d = 2 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) / (sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + 2 * data.eps2 * sinh(0.5 / sqrt(data.eps1)));

    u_ex = 1 + 2 * c * sinh(x / sqrt(data.eps1));
    u_ex(omega2) = -d * (x(omega2) - 1);
    
    du_x_ex = 2 * c * cosh(x / sqrt(data.eps1)) / sqrt(data.eps1);
    du_x_ex(omega2) = -d;
    
    du_y_ex = zeros(size(du_x_ex));
    
    # Assemble system.
    A = bim2a_quadtree_laplacian(msh, epsilon) + ...
        bim2a_quadtree_reaction(msh, !omega2_el, ones(Nnodes, 1));
    
    f = ones(Nelems, 1);
    f(omega2_el) = 0;
    g = ones(Nnodes, 1);
    rhs = bim2a_quadtree_rhs (msh, f, g);
    
    # Compute solution and error.
    dnodes = msh2m_nodes_on_sides(msh, [2 4]);
    u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);
    
    refineable_elements = find(!any(msh.children));
    
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
    
    to_refine(refineable_elements) = (estimator > tol);
    
    # Save solution to file.
    fclose all;
    basename = "./sol_discontinuous_homogeneous/sol";
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
    global_error(i) = norm(bim2c_quadtree_pde_error_L2_node(msh, u, u_ex), 2);
    capacitance(i) = compute_capacitance(msh, A, u);
    
    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", ...
         "err_edge_norm", "err_node_norm", ...
         "global_estimator", "global_error", "capacitance");
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    # Do refinement.
    if (n_elems(i) >= Nelems_max || (tol <= tol_max / n_elems(i)))
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
