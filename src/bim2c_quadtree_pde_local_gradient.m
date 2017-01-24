function [du_x, du_y] = bim2c_quadtree_pde_local_gradient(msh, u, iel)
    if (!isscalar(iel))
        error ("bim2c_quadtree_pde_local_gradient: iel must be a scalar.");
    endif
    
    nodes = msh.t(1:4, iel)(:);
    
    hx = diff(msh.p(1, nodes([1, 2])));
    hy = diff(msh.p(2, nodes([1, 3])));
    
    u_loc = u(nodes);
    
    # grad(u) = sum(u_i * J^(-T) * grad(hat(phi)_i)).
    du_x = u_loc .* msh.shg(1, :, iel)(:) / hx;
    du_y = u_loc .* msh.shg(2, :, iel)(:) / hy;
endfunction
