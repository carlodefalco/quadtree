function [dphi_x, dphi_y] = bim2c_quadtree_eval_dfun_nedelec(msh, iel, x, y)
    x = x(:).';
    y = y(:).';
    
    x1 = msh.p(1, msh.t(1, iel));
    y1 = msh.p(2, msh.t(1, iel));
    
    x2 = msh.p(1, msh.t(2, iel));
    y2 = msh.p(2, msh.t(3, iel));
    
    hx = x2 - x1;
    hy = y2 - y1;
    
    dphi_x = [     -(y - y2);
              zeros(size(y));
                    (y - y1);
              zeros(size(y))].' / hy;
    
    dphi_y = [zeros(size(x));
                    (x - x1);
              zeros(size(x));
                   -(x - x2)].' / hx;
endfunction
