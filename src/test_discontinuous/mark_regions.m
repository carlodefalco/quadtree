function [msh] = mark_regions(msh)
    region = 1;
    
    # Mark elements in (region + 1).
    for iel = 1 : columns(msh.t)
        x = msh.p(1, msh.t(1:4, iel));
        
        if (min(x) >= 0.5)
            msh.t(5, iel) = region + 1;
        endif
    endfor

    # Mark edges in (region + 1).
    for iedge = 1 : columns(msh.e)
        x = msh.p(1, msh.e(1:2, iedge));
        
        if (min(x) >= 0.5)
            msh.e(7, iedge) = region + 1;
        endif
    endfor

    # Mark internal edges.
    l1 = find(msh.p(1, :) == 0.5);
    [~, idx] = sort(msh.p(1, l1));
    l1 = l1(idx);

    e1 = l1(1:end-1);
    e2 = l1(2:end);

    ne = numel(e1);
    newside = 5;

    newedges = [e1;
                e2;
                zeros(2, ne);
                newside * ones(1, ne);
                region * ones(1, ne);
                region * ones(1, ne) + 1];
    
    # Unmark nodes hanging from one side but not from the other.
    hanging_nodes = find(any(msh.hanging(:, l1)));
    for i = 1 : numel(hanging_nodes)
        idx = hanging_nodes(i);
        
        elements = find(any(msh.t(1:4, :) == l1(idx)));
        regions = unique(msh.t(5, elements));
        
        if (numel(regions) == 1)
            row = 6 * (regions == region + 1) + 7 * (regions == region);
            cols = find(any(newedges(1:2, :) == l1(idx)));
            
            newedges(row, cols) = 0;
        endif
    endfor
    
    msh.e = [msh.e newedges];

    msh.onboundary(l1) = newside;
endfunction
