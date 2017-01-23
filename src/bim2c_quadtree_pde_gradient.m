function [du] = bim2c_quadtree_pde_gradient(msh, u)
    Nsides = columns(msh.sides);
    du_x = du_y = zeros(Nsides, 1);
    
    for ii = 1 : Nsides
        side = msh.sides(:, ii);
        
        if (msh.orien(ii))
            h = diff(msh.p(1, side)); # hx.
        else
            h = diff(msh.p(2, side)); # hy.
        endif
        
        du(ii) = diff(u(side)) / h;
    endfor
endfunction
