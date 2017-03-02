function A = test_sol2_projection(msh, D)
    refineable_elements = find(!any(msh.children));
    
    persistent A_loc = [1 -0.5 0 -0.5;
                        -0.5 1 -0.5 0;
                        0 -0.5 1 -0.5;
                        -0.5 0 -0.5 1];
     
     Nnodes = columns(msh.p);
 
     for inode = 1 : 4
         for jnode = 1 : 4
             II(inode, jnode, :) = msh.t(inode, refineable_elements);
             JJ(inode, jnode, :) = msh.t(jnode, refineable_elements);
             VV(inode, jnode, :) = D(refineable_elements) .* A_loc(inode, jnode);
         endfor
     endfor
     
     A = sparse(II, JJ, VV, Nnodes, Nnodes);
     A = msh.P.' * A * msh.P;
endfunction
