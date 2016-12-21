function to_refine = msh2m_to_refine(msh, A_fun, rhs_fun, u, iel, tol)
    nodes = msh.t(1:4, iel);

    x_iel = [min(msh.p(1, nodes)), max(msh.p(1, nodes))];
    y_iel = [min(msh.p(2, nodes)), max(msh.p(2, nodes))];

    hx = diff(x_iel);
    hy = diff(y_iel);
    
    msh_iel = msh2m_quadtree(x_iel, y_iel);
    msh_iel = msh2m_refine_quadtree(msh_iel, 1);

    # Build local matrix.
    x_iel = msh_iel.p(1, :).';
    y_iel = msh_iel.p(2, :).';
    
    # Assemble system.
    A_iel = A_fun(msh_iel);


    rhs_iel = rhs_fun(msh_iel);

    # Boundary conditions.
    dnodes_iel = msh2m_nodes_on_sides(msh_iel, 1:4); # Dirichlet nodes.
    
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
    
    to_refine = (abs(u_iel(end) - mean(u(nodes))) > tol * sqrt(hx * hy));
endfunction
