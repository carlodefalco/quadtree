function u = bim2a_quadtree_solve(msh, A, f, u, dnodes)
    # Solve for non-hanging and non-Dirichlet nodes.
    non_hanging = msh.full_to_reduced(msh.full_to_reduced != 0);
    reduced_dnodes = msh.full_to_reduced(dnodes);
    
    intnodes = setdiff(non_hanging, reduced_dnodes);
    
    u(intnodes) = A(intnodes, intnodes) \ ...
                  (f(intnodes) - A(intnodes, reduced_dnodes) * u(reduced_dnodes));
    
    # Interpolate solution at the hanging nodes.
    hanging = find(!msh.full_to_reduced);
    
    for i = 1 : length(hanging)
        idx = hanging(i);
        
        u(idx) = 0.5 * (u(msh.hanging(1, idx)) + u(msh.hanging(2, idx)));
    endfor
endfunction
