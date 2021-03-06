function u = bim2a_quadtree_solve(msh, A, rhs, u, dnodes)
    # Solve for non-hanging and non-Dirichlet nodes.
    non_hanging = msh.full_to_reduced(msh.full_to_reduced != 0);
    reduced_dnodes = msh.full_to_reduced(dnodes);
    
    intnodes = setdiff(non_hanging, reduced_dnodes);
    
    u(msh.reduced_to_full(intnodes)) = ...
        A(intnodes, intnodes) \ ...
        (rhs(intnodes) - A(intnodes, reduced_dnodes) * u(dnodes));
    
    # Interpolate solution at the hanging nodes.
    hanging = find(!msh.full_to_reduced);
    u(hanging) = mean(u(msh.hanging(:, hanging)));
endfunction
