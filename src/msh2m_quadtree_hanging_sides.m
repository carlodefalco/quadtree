function [hanging_sides] = msh2m_quadtree_hanging_sides(msh)
    hanging_sides = zeros(1, columns(msh.sides));
    
    for ii = 1 : columns(msh.p)
        idx = find(all(msh.sides == msh.hanging(:, ii)
                       | flip(msh.sides) == msh.hanging(:, ii)));
        
        if (!isempty(idx))
            sidelist = find(msh.sides(1, :) == ii
                            | msh.sides(2, :) == ii);
            
            sidelist = sidelist(msh.orien(sidelist) == msh.orien(idx));
            hanging_sides(sidelist) = idx;
        endif
    endfor
endfunction
