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

    x = msh.p(1, :).';
    y = msh.p(2, :).';

    u_ex = sin(x) .* cos(2 * y);
    du_x_ex = cos(x) .* cos(2 * y);
    du_y_ex = -2 * sin(x) .* sin(2 * y);

    du_edge_ex = zeros(columns(msh.sides), 1);
    orien_x = (msh.orien == 1);
    orien_y = (msh.orien == 0);
    du_edge_ex(orien_x) = mean(du_x_ex(msh.sides(:, orien_x)));
    du_edge_ex(orien_y) = mean(du_y_ex(msh.sides(:, orien_y)));
    #du_edge_ex = bim2c_quadtree_pde_edge_gradient(msh, u_ex);
    #[du_x_ex, du_y_ex] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge_ex);
    
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

    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u);
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);
    
    ## Compute errors.
    err_node(i) = norm([du_x; du_y] - [du_x_ex; du_y_ex], inf);
    err_edge(i) = norm(du_edge - du_edge_ex, inf);
endfor

h = 1 ./ n;
figure;
loglog(h, err_node, h, h);
xlabel("h");
ylabel("err");
legend("Error on node space", "Order 1");

figure;
loglog(h, err_edge, h, h.^2);
xlabel("h");
ylabel("err");
legend("Error on edge space", "Order 2");
