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

    x = msh.p(1, :).';
    y = msh.p(2, :).';

    # Define parameters and exact solution.
    epsilon = 1;
    lambda = 1 / epsilon;

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

    # Save solution to file.
    fclose all;
    filename = sprintf("sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex"}, {}, 1);

    # Determine elements to be refined.
    tol = 1e-3;
    to_refine = false(1, Nelements);
    
    refineable_elements = find(!any(msh.children));
    
    for iel = refineable_elements
        nodes = msh.t(1:4, iel);

        x_iel = [min(msh.p(1, nodes)), max(msh.p(1, nodes))];
        y_iel = [min(msh.p(2, nodes)), max(msh.p(2, nodes))];

        msh_iel = msh2m_quadtree(x_iel, y_iel);
        msh_iel = msh2m_refine_quadtree(msh_iel, 1);

        # Build local matrix.
        x_iel = msh_iel.p(1, :).';
        y_iel = msh_iel.p(2, :).';
        
        # Assemble system.
        A_iel = bim2a_quadtree_advection_diffusion(msh_iel, alpha(msh_iel), psi(x_iel, y_iel));

        dnodes_iel = msh2m_nodes_on_sides(msh_iel, 1:4);

        rhs_iel = bim2a_quadtree_rhs(msh_iel, f(msh_iel), g(x_iel, y_iel));

        # Compute boundary conditions.
        # Original vertices.
        u_dnodes = zeros(2 * numel(nodes) + 1, 1);
        u_dnodes(msh_iel.t(1:4, 1)) = u(nodes);

        # Added vertices (side midpoints).
        hanging_sides = [1, 3, 4, 2;
                         3, 4, 2, 1];
        u_dnodes((numel(nodes)+1) : end - 1) = mean(u_dnodes(hanging_sides));
        
        # Vertex in the middle of the element.
        u_dnodes(end) = 0;

        # Compute solution and evaluate error.
        u_iel = bim2a_quadtree_solve(msh_iel, A_iel, rhs_iel, u_dnodes, dnodes_iel);
        
        if (abs(u_iel(end) - mean(u(nodes))) > tol)
            to_refine(iel) = true;
        endif
    endfor
    
    fprintf("Elements to refine = %d\n\n", sum(to_refine));
    
    if (!any(to_refine))
        break;
    else
        msh = msh2m_refine_quadtree_recursive(msh, find(to_refine));
    endif
endfor

# Plot refined mesh.
quadmesh(msh);
title("New mesh");
xlabel("x");
ylabel("y");
