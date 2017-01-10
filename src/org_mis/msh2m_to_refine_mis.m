function [to_refine] = msh2m_to_refine_mis(msh, material, constants,
                                           A_fun, M_fun, phi,
                                           charge_n, iel, tol)
    nodes = msh.t(1:4, iel);

    x_iel = [min(msh.p(1, nodes)), max(msh.p(1, nodes))];
    y_iel = [min(msh.p(2, nodes)), max(msh.p(2, nodes))];

    hx = diff(x_iel);
    hy = diff(y_iel);
    
    msh_iel = msh2m_quadtree(x_iel, y_iel);
    msh_iel = msh2m_refine_quadtree(msh_iel, 1);
    
    msh_iel.x_min    = msh.x_min;
    msh_iel.x_sc_min = msh.x_sc_min;
    msh_iel.x_sc_max = msh.x_sc_max;
    msh_iel.x_max    = msh.x_max;
    msh_iel.y_sc     = msh.y_sc;
    msh_iel.y_ins    = msh.y_ins;

    # Build local matrix.
    x_iel = msh_iel.p(1, :).';
    y_iel = msh_iel.p(2, :).';
    
    # Assemble system.
    A_iel = A_fun(msh_iel);
    M_iel = M_fun(msh_iel);

    # Boundary conditions.
    dnodes_iel = msh2m_nodes_on_sides(msh_iel, 1:4); # Dirichlet nodes.
    
    # Original vertices.
    phi_dnodes = zeros(2 * numel(nodes) + 1, 1);
    phi_dnodes(msh_iel.t(1:4, 1)) = phi(nodes);

    # Added vertices (side midpoints).
    hanging_sides = [1, 3, 4, 2;
                     3, 4, 2, 1];
    phi_dnodes((numel(nodes)+1) : end - 1) = mean(phi_dnodes(hanging_sides));
    
    # Compute solution and evaluate error.
    phi_iel = nlpoisson(msh_iel, phi_dnodes, A_iel, M_iel, dnodes_iel, charge_n);
    
    to_refine = (abs(phi_iel(end) - mean(phi(nodes))) > tol * sqrt(hx * hy));
endfunction
