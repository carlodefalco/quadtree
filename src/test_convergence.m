clear all;
close all;
clc;

n = 10 * 2.^(0:3);

for i = 1 : length(n)
    fprintf("i = %d\n", i);
    
    # Mesh definition.
    x = linspace(0, 2*pi, n(i));
    y = linspace(0, 2*pi, n(i));

    region = 1;
    sides = 1:4;

    msh = msh2m_quadtree(x, y, region, sides);
    
    # Build global matrix.
    Nnodes = columns(msh.p);
    Nelements = columns(msh.t);

    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    # Define parameters and exact solution.
    epsilon = 1;

    u_ex = sin(x) .* cos(2*y);
    
    # Assemble system.
    alpha = @(msh) epsilon * ones(columns(msh.t), 1);
    psi = @(x, y) x / epsilon;
    
    A = bim2a_quadtree_advection_diffusion(msh, alpha(msh), psi(x, y));

    dnodes = msh2m_nodes_on_sides(msh, 1:4);

    f = @(msh) ones(columns(msh.t), 1);
    g = @(x, y) cos(2*y) .* (cos(x) + 5 * epsilon * sin(x));
    rhs = bim2a_quadtree_rhs(msh, f(msh), g(x, y));
    
    # Compute solution and error.
    u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);

    err(i) = norm(u - u_ex, inf);
    
    # Save solution to file.
    fclose all;
    filename = sprintf("sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex"}, {}, 1);
endfor
