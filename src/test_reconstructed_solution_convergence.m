clear all;
close all;
clc;

n = 10 * 2.^(0:3);

for i = 1 : numel(n)
    fprintf("i = %d, n = %d\n", i, n(i));
    
    # Mesh definition.
    x = linspace(0, 1, n(i));
    y = linspace(0, 1, n(i));
    
    msh = msh2m_quadtree(x, y);
    [msh, ~, node_space, q2_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    
    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    # Define parameters and exact solution.
    u_ex = sin(x) .* cos(2 * y);
    
    # Assemble system.
    alpha = @(msh) ones(columns(msh.t), 1);
    psi = @(x, y) x;
    
    A = bim2a_quadtree_advection_diffusion(msh, alpha(msh), psi(x, y));

    dnodes = msh2m_nodes_on_sides(msh, 1:4);

    f = @(msh) ones(columns(msh.t), 1);
    g = @(x, y) cos(2 * y) .* (cos(x) + 5 * sin(x));
    rhs = bim2a_quadtree_rhs(msh, f(msh), g(x, y));

    u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);
    
    # Evaluate u on quadrature nodes for each element.
    Nquad = size(node_space.shp, 1); # No. of quadrature nodes.
    Nconn = size(node_space.connectivity);
    u_ex_ = reshape(u_ex(node_space.connectivity), [1, Nconn]);
    u_ex_ = repmat(u_ex_, [Nquad, 1, 1]);
    
    u_ex_node = squeeze(sum(node_space.shp .* u_ex_, 2));
    
    # Evaluate u_star on quadrature nodes for each element.
    [u_star_node, u_star_edge, u_star_center] = bim2c_quadtree_pde_reconstructed_solution(msh, u);
    
    Nconn = size(q2_space.connectivity);
    u_star = reshape([u_star_node(q2_space.connectivity(1:4, :));
                      u_star_edge(q2_space.connectivity(5:8, :));
                      u_star_center(q2_space.connectivity(9, :)).'], [1, Nconn]);
    u_star = repmat(u_star, [Nquad, 1, 1]);
    
    u_q2 = squeeze(sum(q2_space.shp .* u_star, 2));
    
    err = (u_ex_node - u_q2).^2;
    err = sqrt(sum(err .* msh.wjacdet, 1));
    err_norm(i) = norm(err, 2);
    
    estimator = bim2c_quadtree_pde_ZZ_estimator_u(msh, u);
    
    # Save solution to file.
    fclose all;
    filename = sprintf("sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, ...
                                 {u, "u";
                                  u_ex, "u_ex";
                                  abs(u - u_ex), "|u - u_ex|"},
                                 {err(:), "err";
                                  estimator(:), "estimator";
                                  u_star_center(q2_space.connectivity(9, :)), "uu"}, 1);
endfor

h = 1 ./ n;
figure;
loglog(h, err_norm, h, h.^2);
xlabel("h");
ylabel("err");
legend("Error", "Order 2");
