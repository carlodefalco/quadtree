function [phi] = bim2c_quadtree_eval_fun(msh, iel, x, y)
    x = x(:).';
    y = y(:).';
    
    x1 = msh.p(1, msh.t(1, iel));
    y1 = msh.p(2, msh.t(1, iel));
    
    x2 = msh.p(1, msh.t(2, iel));
    y2 = msh.p(2, msh.t(3, iel));
    
    hx = x2 - x1;
    hy = y2 - y1;
    
    phi = [ (x - x2) .* (y - y2);
           -(x - x1) .* (y - y2);
            (x - x1) .* (y - y1);
           -(x - x2) .* (y - y1)].' / (hx * hy);
endfunction
