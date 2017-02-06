function [msh] = mark_regions_2D(msh)
    region = 1;
    
    # Mark elements to discard.
    for iel = 1 : columns(msh.t)
        coords = msh.p(:, msh.t(1:4, iel));
        
        if (max(coords(1, :)) <= msh.dim.x_sc_min ||
            min(coords(1, :)) >= msh.dim.x_sc_max ||
            min(coords(2, :)) >= 0)
            
            msh.t(5, iel) = region + 1;
        endif
    endfor

    # Mark edges to discard.
    for iedge = 1 : columns(msh.e)
        coords = msh.p(:, msh.e(1:2, iedge));
        
        if (max(coords(1, :)) <= msh.dim.x_sc_min ||
            min(coords(1, :)) >= msh.dim.x_sc_max ||
            min(coords(2, :)) >= 0)
            
            msh.e(7, iedge) = region + 1;
        endif
    endfor

    # Mark internal edges.
    l1 = find(msh.p(1, :) >= msh.dim.x_sc_min &
              msh.p(1, :) <= msh.dim.x_sc_max & 
              msh.p(2, :) == 0);
    [~, idx] = sort(msh.p(1, l1));
    l1 = l1(idx);
    
    l2 = find(msh.p(1, :) == msh.dim.x_sc_max & msh.p(2, :) <= 0);
    [~, idx] = sort(msh.p(2, l2));
    l2 = l2(idx);

    e1 = [l1(1:end-1) l2(1:end-1)];
    e2 = [l1(2:end)   l2(2:end)  ];

    ne = numel(e1);
    newside = 5;

    newedges = [e1;
                e2;
                zeros(2, ne);
                newside * ones(1, ne);
                region * ones(1, ne);
                region * ones(1, ne) + 1];
    
    # Unmark nodes hanging from one side but not from the other.
    l = union(l1, l2);
    hanging_nodes = find(any(msh.hanging(:, l)));
    for i = 1 : numel(hanging_nodes)
        idx = hanging_nodes(i);
        
        elements = find(any(msh.t(1:4, :) == l(idx)));
        regions = unique(msh.t(5, elements));
        
        if (isscalar(regions))
            row = 6 * (regions == region + 1) + 7 * (regions == region);
            cols = find(any(newedges(1:2, :) == l(idx)));
            
            newedges(row, cols) = 0;
        endif
    endfor
    
    msh.e = [msh.e newedges];

    msh.onboundary(l) = newside;
endfunction
