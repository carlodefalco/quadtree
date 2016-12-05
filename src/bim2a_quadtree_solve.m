function u = bim2a_quadtree_solve(A, f, msh)
    u = zeros(columns(msh.p), 1);
    
    # Solve for non-hanging nodes.
    non_hanging = find(msh.full_to_reduced);
    
    u(non_hanging) = A \ f;
    
    # Interpolate solution at the hanging nodes.
    hanging = find(!msh.full_to_reduced);
    
    for i = 1 : length(hanging)
        idx = hanging(i);
        
        u(idx) = 0.5 * (u(msh.hanging(1, idx)) + u(msh.hanging(2, idx)));
    endfor
endfunction
