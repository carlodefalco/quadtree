function [msh] = mark_regions(msh)
    region = 1;
    
    # Mark elements to discard.
    for iel = 1 : columns(msh.t)
        coords = msh.p(:, msh.t(1:4, iel));
        
        if (min(coords(1, :)) >= 0.5)
            msh.t(5, iel) = region + 1;
        endif
    endfor

    # Mark edges to discard.
    for iedge = 1 : columns(msh.e)
        x = msh.p(1, msh.e(1:2, iedge));
        
        if (min(x) >= 0.5)
            msh.e(7, iedge) = region + 1;
        endif
    endfor

    # Mark internal edges.
    l1 = find(msh.p(1, :) == 0.5);

    e1 = l1(1:end-1);
    e2 = l1(2:end);

    ne = numel(e1);
    newside = 5;

    newedges = [e1;
                e2;
                zeros(2, ne);
                newside * ones(1, ne);
                region(ones(1, ne));
                region(ones(1, ne)) + 1];

    msh.e = [msh.e newedges];

    msh.onboundary(l1) = newside;
endfunction
