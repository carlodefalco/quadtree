clear all;
close all;
clc;

n = 10 * 2.^(0:4);

for i = 1 : numel(n)
    fprintf("i = %d, n = %d\n", i, n(i));
    
    # Mesh definition.
    x = linspace(0, 1, n(i));
    y = linspace(0, 1, n(i));

    msh = msh2m_quadtree(x, y);

    x = msh.p(1, :).';
    y = msh.p(2, :).';

    u_ex = sin(2 * pi * x) .* sin(2 * pi * y);
    du_x_ex = 2 * pi * cos(2 * pi * x) .* sin(2 * pi * y);
    du_y_ex = 2 * pi * sin(2 * pi * x) .* cos(2 * pi * y);

    du_edge_ex = zeros(columns(msh.sides), 1);
    orien_x = (msh.orien == 1);
    orien_y = (msh.orien == 0);
    du_edge_ex(orien_x) = mean(du_x_ex(msh.sides(:, orien_x)));
    du_edge_ex(orien_y) = mean(du_y_ex(msh.sides(:, orien_y)));

    ## Compute gradient on edge space.
    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u_ex);

    ## Compute reconstructed gradient on node space.
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);
    
    ## Compute errors.
    err_edge(i) = norm(du_edge - du_edge_ex, inf);
    err_node(i) = norm([du_x; du_y] - [du_x_ex; du_y_ex], inf);
endfor

h = 1 ./ n;
figure;
loglog(h, err_node, h, err_edge, h, h.^2);
legend("Error on node space", "Error on edge space", "Order 2");
