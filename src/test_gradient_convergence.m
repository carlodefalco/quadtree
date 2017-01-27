clear all;
close all;
clc;

n = 10 * 2.^(0:3);

for i = 1 : numel(n)
    fprintf("i = %d, n = %d\n", i, n(i));
    
    # Mesh definition.
    x = linspace(0, 2 * pi, n(i));
    y = linspace(0, 2 * pi, n(i));

    msh = msh2m_quadtree(x, y);

    [msh, edge_space, node_space] = bim2c_quadtree_mesh_properties(msh);
    Nquad = size(edge_space.shp, 2); # No. of quadrature nodes.

    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    u_ex = sin(x) .* cos(2 * y);
    du_x_ex = cos(x) .* cos(2 * y);
    du_y_ex = -2 * sin(x) .* sin(2 * y);
    
    # Define parameters and exact solution.
    epsilon = 1;
    
    # Assemble system.
    alpha = @(msh) epsilon * ones(columns(msh.t), 1);
    psi = @(x, y) x / epsilon;
    
    A = bim2a_quadtree_advection_diffusion(msh, alpha(msh), psi(x, y));

    dnodes = msh2m_nodes_on_sides(msh, 1:4);

    f = @(msh) ones(columns(msh.t), 1);
    g = @(x, y) cos(2 * y) .* (cos(x) + 5 * epsilon * sin(x));
    rhs = bim2a_quadtree_rhs(msh, f(msh), g(x, y));

    u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);
    
    # Compute numerical gradient on edge and node space.
    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u);
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);
    
    # Fix boundary values.
    boundary = msh2m_nodes_on_sides(msh, 1:4);
    du_x(boundary) = du_x_ex(boundary);
    du_y(boundary) = du_y_ex(boundary);

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
    err_edge = sqrt(sum(err .* msh.wjacdet, 1));
    err_edge_norm(i) = sqrt(sum(err_edge.^2));
    
    err = (du_x_node - du_x_ex).^2 + (du_y_node - du_y_ex).^2;
    err_node = sqrt(sum(err .* msh.wjacdet, 1));
    err_node_norm(i) = sqrt(sum(err_node.^2));
    
    estimator = bim2c_quadtree_pde_ZZ_estimator_du(msh, u);
    
    # Save solution to file.
    fclose all;
    filename = sprintf("sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex"},
                                 {err_edge.', "err_edge"; err_node.', "err_node"; estimator.', "estimator"}, 1);
endfor

h = 1 ./ n;
figure;
loglog(h, err_node_norm, h, h.^2);
xlabel("h");
ylabel("err");
legend("Error on node space", "Order 2");

figure;
loglog(h, err_edge_norm, h, h);
xlabel("h");
ylabel("err");
legend("Error on edge space", "Order 1");
