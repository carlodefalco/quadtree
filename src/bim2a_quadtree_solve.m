function u = bim2a_quadtree_solve(msh, A, f, u, dnodes)
    # Solve for non-hanging and non-Dirichlet nodes.
    non_hanging = find(msh.full_to_reduced);
    intnodes = setdiff(non_hanging, dnodes);
    
    u(intnodes) = A(intnodes, intnodes) \ ...
                  (f(intnodes) - A(intnodes, dnodes) * u(dnodes));
    
    # Interpolate solution at the hanging nodes.
    hanging = find(!msh.full_to_reduced);
    
    for i = 1 : length(hanging)
        idx = hanging(i);
        
        u(idx) = 0.5 * (u(msh.hanging(1, idx)) + u(msh.hanging(2, idx)));
    endfor
endfunction
