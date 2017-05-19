function [msh] = mark_contacts(msh)
    # Mark edge labels.
    for iedge = 1 : columns(msh.e)
        coords = msh.p(:, msh.e(1:2, iedge));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        
        if (msh.e(5, iedge) == 1)
            if (x_max <= msh.dim.x_source)
                msh.e(5, iedge) = 5;
            elseif (x_min >= msh.dim.x_drain)
                msh.e(5, iedge) = 6;
            end
        endif
    endfor
    msh.onboundary(msh2m_nodes_on_sides(msh, 5)) = 5;
    msh.onboundary(msh2m_nodes_on_sides(msh, 6)) = 6;
endfunction
