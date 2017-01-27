clear all;
close all;
clc;

n = 10;

# Mesh definition.
x = linspace(0, 2 * pi, n);
y = linspace(0, 2 * pi, n);

region = 1;
sides = 1:4;

msh = msh2m_quadtree(x, y, region, sides);

for i = 1:10
    fprintf("i = %d\n", i);
    
    # Build global matrix.
    Nnodes = columns(msh.p);
    Nelements = columns(msh.t);

    # Define parameters and exact solution.
    epsilon = 1;

    u_ex = (sin(msh.p(1, :)) .* cos(2*msh.p(2, :))).';

    # Assemble system.
    alpha = @(msh) epsilon * ones(columns(msh.t), 1);
    psi = @(msh) msh.p(1, :) / epsilon;

    A = @(msh) bim2a_quadtree_advection_diffusion(msh, alpha(msh), psi(msh));
    
    f = @(msh) ones(columns(msh.t), 1);
    g = @(msh) cos(2*msh.p(2, :)) .* (cos(msh.p(1, :)) + 5 * epsilon * sin(msh.p(1, :)));
    
    rhs = @(msh) bim2a_quadtree_rhs(msh, f(msh), g(msh));

    # Compute solution and error.
    dnodes = msh2m_nodes_on_sides(msh, 1:4); # Dirichlet nodes.
    
    u = bim2a_quadtree_solve(msh, A(msh), rhs(msh), u_ex, dnodes);

    # Save solution to file.
    fclose all;
    filename = sprintf("sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex"; u - u_ex, "err"}, {}, 1);

    # Determine elements to be refined.
    tol = 1e-2;
    refineable_elements = find(!any(msh.children));
    
    to_refine = false(1, Nelements);
    to_refine(refineable_elements) = ...
        parcellfun(4, @(iel) msh2m_to_refine(msh, A, rhs, u, iel, tol),
                   num2cell(refineable_elements));
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    if (!any(to_refine))
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
