function [phi] = bim2c_quadtree_eval_fun_q2(msh, iel, x, y)
    x = x(:).';
    y = y(:).';
    
    x1 = msh.p(1, msh.t(1, iel));
    y1 = msh.p(2, msh.t(1, iel));    
    
    x3 = msh.p(1, msh.t(2, iel));
    y3 = msh.p(2, msh.t(3, iel));

    x2 = (x1 + x3) / 2;
    y2 = (y1 + y3) / 2;
    
    hx = x3 - x1;
    hy = y3 - y1;
    
    phi = [ 4 * (x - x2) .* (x - x3) .* (y - y2) .* (y - y3);
            4 * (x - x1) .* (x - x2) .* (y - y2) .* (y - y3);
            4 * (x - x1) .* (x - x2) .* (y - y1) .* (y - y2);
            4 * (x - x2) .* (x - x3) .* (y - y1) .* (y - y2);
           -8 * (x - x1) .* (x - x3) .* (y - y2) .* (y - y3);
           -8 * (x - x1) .* (x - x2) .* (y - y1) .* (y - y3);
           -8 * (x - x1) .* (x - x3) .* (y - y1) .* (y - y2);
           -8 * (x - x2) .* (x - x3) .* (y - y1) .* (y - y3);
           16 * (x - x1) .* (x - x3) .* (y - y1) .* (y - y3)].' / (hx^2 * hy^2);
    
    # Fix zero values outside iel.
    phi(x < x1 | x > x3, :) = 0;
    phi(y < y1 | y > y3, :) = 0;
endfunction
