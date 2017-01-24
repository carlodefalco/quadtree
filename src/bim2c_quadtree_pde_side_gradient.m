function [du] = bim2c_quadtree_pde_side_gradient(msh, u)
    Nsides = columns(msh.sides);
    du = zeros(Nsides, 1);
    
    for ii = 1 : Nsides
        side = msh.sides(:, ii);
        
        # h = hx if msh.orien(ii) == true;
        # h = hy if msh.orien(ii) == false.
        orien = double(msh.orien(ii));
        h = diff(msh.p(2 - orien, side));
        
        ## Compute tangential component of gradient at side ii.
        du(ii) = diff(u(side)) / h;
    endfor
endfunction
