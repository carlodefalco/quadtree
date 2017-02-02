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
    data.eps1 = 1;
    data.eps2 = 0.5;
    
    omega1 = omega1_nodes(msh);
    omega2 = setdiff(1:Nnodes, omega1);
    
    omega2_el = omega2_elems(msh);
    epsilon = diffusion_coefficient(msh, data, omega2_el);
    
    c = -0.4375 * data.eps2 / (0.5 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + data.eps2 * sinh(0.5 / sqrt(data.eps1)));
    d = 1.75 * sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) / (sqrt(data.eps1) * cosh(0.5 / sqrt(data.eps1)) + 2 * data.eps2 * sinh(0.5 / sqrt(data.eps1)));
    
    u_ex = 1 + 2 * c * sinh(x / sqrt(data.eps1));
    u_ex(omega2) = -0.5 * (x(omega2) - 1) .* (x(omega2) + 2 * d);
    
    du_x_ex = 2 * c * cosh(x / sqrt(data.eps1)) / sqrt(data.eps1);
    du_x_ex(omega2) = x(omega2) - d - 0.5;
    
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
    
    [msh, edge_space, node_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    edge_space.shp = edge_space.shp(:, :, :, refineable_elements);
    edge_space.connectivity = edge_space.connectivity(:, refineable_elements);
    
    node_space.shp = node_space.shp(:, :, refineable_elements);
    node_space.connectivity = node_space.connectivity(:, refineable_elements);
    
    # Compute numerical gradient on edge and node space.
    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u);
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);

    Nquad = size(edge_space.shp, 2); # No. of quadrature nodes.
    
    # Evaluate exact gradient on quadrature nodes.
    Nconn = size(node_space.connectivity);
    
    du_x_ex = reshape(du_x_ex(node_space.connectivity), [1, Nconn]);
    du_x_ex = repmat(du_x_ex, [Nquad, 1, 1]);
    
    du_y_ex = reshape(du_y_ex(node_space.connectivity), [1, Nconn]);
    du_y_ex = repmat(du_y_ex, [Nquad, 1, 1]);
    
    du_x_ex = squeeze(sum(node_space.shp .* du_x_ex, 2));
    du_y_ex = squeeze(sum(node_space.shp .* du_y_ex, 2));
    
    # Evaluate edge gradient on quadrature nodes.
    Nconn = size(edge_space.connectivity);
    du_edge = reshape(du_edge(edge_space.connectivity), [1, Nconn]);
    du_edge = repmat(du_edge, [Nquad, 1, 1]);
    
    du_x_edge = squeeze(sum(squeeze(edge_space.shp(1, :, :, :)) .* du_edge, 2));
    du_y_edge = squeeze(sum(squeeze(edge_space.shp(2, :, :, :)) .* du_edge, 2));
    
    # Evaluate node gradient on quadrature nodes.
    Nconn = size(node_space.connectivity);
    du_x = reshape(du_x(node_space.connectivity), [1, Nconn]);
    du_x = repmat(du_x, [Nquad, 1, 1]);
    
    du_y = reshape(du_y(node_space.connectivity), [1, Nconn]);
    du_y = repmat(du_y, [Nquad, 1, 1]);
    
    du_x_node = squeeze(sum(node_space.shp .* du_x, 2));
    du_y_node = squeeze(sum(node_space.shp .* du_y, 2));
    
    # Compute errors and estimator.
    err = (du_x_edge - du_x_ex).^2 + (du_y_edge - du_y_ex).^2;
    err_edge = sqrt(sum(err .* msh.wjacdet(:, refineable_elements), 1));
    err_edge_norm(i) = norm(err_edge, 2);
    
    err = (du_x_node - du_x_ex).^2 + (du_y_node - du_y_ex).^2;
    err_node = sqrt(sum(err .* msh.wjacdet(:, refineable_elements), 1));
    err_node_norm(i) = norm(err_node, 2);
    
    # Determine elements to be refined.
    to_refine = false(1, Nelems);
    
    estimator = bim2c_quadtree_pde_ZZ_estimator_du(msh, u);
    tol = mean(estimator);
    
    to_refine(refineable_elements) = (estimator > tol);
    
    # Save solution to file.
    fclose all;
    basename = "./sol_discontinuous/sol";
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
    
    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", ...
         "err_edge_norm", "err_node_norm", ...
         "global_estimator", "global_error");
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    # Do refinement.
    if (n_elems(i) >= Nelems_max || (tol <= tol_max / n_elems(i)))
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
