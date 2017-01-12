function [du_x, du_y] = bim2c_quadtree_pde_local_gradient(msh, u, iel)
    nodes = msh.t(1:4, iel)(:);
    
    u_loc = u(nodes);
    
    du_x = u_loc .* msh.shg(1, :, iel)(:);
    du_y = u_loc .* msh.shg(2, :, iel)(:);
endfunction
