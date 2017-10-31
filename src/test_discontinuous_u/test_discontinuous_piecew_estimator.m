clear all;
close all;
clc;

addpath(canonicalize_file_name("../"));

x = linspace(0, 1, 11);
y = linspace(0, 1, 11);

msh = msh2m_quadtree(x, y);

Nelems_max = 10000;

for i = 1 : 10
    fprintf("i = %d\n", i);
    msh = bim2c_quadtree_mesh_properties(msh);
    
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
    
    [msh1, omega1, omega1_el_all] = msh2m_quadtree_submesh(msh_new, 1);
    [msh2, omega2, omega2_el_all] = msh2m_quadtree_submesh(msh_new, 2);
    
    [~, omega1_el] = intersect(refineable_elements, omega1_el_all);
    [~, omega2_el] = intersect(refineable_elements, omega2_el_all);
    
    err_L2 = bim2c_quadtree_pde_error_L2_node(msh, u, u_ex);
    err_L2_1 = bim2c_quadtree_pde_error_L2_node(msh1, u(omega1), u_ex(omega1));
    err_L2_2 = bim2c_quadtree_pde_error_L2_node(msh2, u(omega2), u_ex(omega2));
    
    # Determine elements to be refined.
    estimator = zeros(1, numel(refineable_elements));
    estimator(omega1_el) = bim2c_quadtree_pde_ZZ_estimator_u(msh1, u(omega1));
    estimator(omega2_el) = bim2c_quadtree_pde_ZZ_estimator_u(msh2, u(omega2));
    
    crit = 1;
    to_refine = bim2c_quadtree_pde_ZZ_to_refine(msh, estimator, 1e-5 / sqrt(Nnodes), crit);
    
    # Save solution to file.
    fclose all;
    basename = sprintf("./sol_discontinuous_piecew_estimator_crit%d/sol", crit);
    filename = sprintf([basename "_%d_omega"], i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh,
                                 {u, "u"; u_ex, "u_ex"},
                                 {estimator.', "estimator"}, 1);
    
    filename1 = sprintf([basename "_%d_omega1"], i);
    if (exist([filename1 ".vtu"], "file"))
        delete([filename1 ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename1, msh1,
                                 {u(omega1), "u"; u_ex(omega1), "u_ex"},
                                 {estimator(omega1_el)(:), "estimator"}, 1);
    
    filename2 = sprintf([basename "_%d_omega2"], i);
    if (exist([filename2 ".vtu"], "file"))
        delete([filename2 ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename2, msh2,
                                 {u(omega2), "u"; u_ex(omega2), "u_ex"},
                                 {estimator(omega2_el)(:), "estimator"}, 1);
    
    n_dofs(i) = sum(!any(msh.hanging));
    n_elems(i) = numel(refineable_elements);
    n_to_refine(i) = sum(to_refine);
    global_estimator(i) = norm(estimator, 2);
    global_error(i) = norm(err_L2, 2);
    capacitance(i) = compute_capacitance(msh, A, u);
    
    save("-binary", sprintf([basename "_%d_data.mat"], i), "*");
    
    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", ...
         "global_estimator", "global_error", "capacitance");
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    # Do refinement.
    if (n_elems(i) >= Nelems_max || all(!to_refine))
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
