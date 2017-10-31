clear all;
close all;
clc;

n = 10 * 2.^(0:3);

for i = 1 : numel(n)
    fprintf("i = %d, n = %d\n", i, n(i));
    
    # Mesh definition.
    x = linspace(0, 2*pi, n(i));
    y = linspace(0, 2*pi, n(i));
    
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
    
    err = bim2c_quadtree_pde_error_L2_node_q2(msh, u, u_ex);
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
                                  estimator(:), "estimator"}, 1);
endfor

h = 1 ./ n;
figure;
loglog(h, err_norm, h, h.^2);
xlabel("h");
ylabel("err");
legend("Error", "Order 2");
