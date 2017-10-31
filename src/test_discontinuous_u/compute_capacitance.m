function [C] = compute_capacitance(msh, mat, u)
    non_hanging = find(msh.full_to_reduced);
    
    # Boundary conditions.
    left = msh2m_nodes_on_sides(msh, 4);
    right = msh2m_nodes_on_sides(msh, 2);
    dnodes = union(left, right);
    
    x = msh.p(1, :).';
    u0 = x;
    
    # Compute solution.
    rhs = zeros(size(u));
    delta_u = bim2a_quadtree_solve(msh, mat, rhs, u0, dnodes);
    
    # Compute capacitance.
    contact = msh.full_to_reduced(right);
    C = sum(mat(contact, :) * delta_u(non_hanging));
endfunction
