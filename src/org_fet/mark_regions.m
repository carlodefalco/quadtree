function [msh] = mark_regions(msh)
    region = 1;
    
    # Mark elements in (region + 1).
    for iel = 1 : columns(msh.t)
        coords = msh.p(:, msh.t(1:4, iel));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        y_min = min(coords(2, :));
        
        if (x_min >= msh.dim.x_sc_max ||
            x_max <= msh.dim.x_sc_min ||
            y_min >= 0)
            
            msh.t(5, iel) = region + 1;
        endif
    endfor

    # Mark edges in (region + 1).
    for iedge = 1 : columns(msh.e)
        coords = msh.p(:, msh.e(1:2, iedge));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        y_min = min(coords(2, :));
        
        if (x_min >= msh.dim.x_sc_max ||
            x_max <= msh.dim.x_sc_min ||
            y_min >= 0)
            
            msh.e(7, iedge) = region + 1;
        endif
    endfor

    # Mark internal edges.
    # Semiconductor left.
    l1 = find(msh.p(1, :) == msh.dim.x_sc_min & msh.p(2, :) <= 0);
    [~, idx] = sort(msh.p(2, l1));
    l1 = l1(idx);
    
    # Semiconductor right.
    l2 = find(msh.p(1, :) == msh.dim.x_sc_max & msh.p(2, :) <= 0);
    [~, idx] = sort(msh.p(2, l2));
    l2 = l2(idx);

    # Semiconductor top.
    l3 = find(msh.p(1, :) >= msh.dim.x_sc_min &
              msh.p(1, :) <= msh.dim.x_sc_max & 
              msh.p(2, :) == 0);
    [~, idx] = sort(msh.p(1, l3));
    l3 = l3(idx);
    
    e1 = [l1(1:end-1) l2(1:end-1) l3(1:end-1)];
    e2 = [l1(2:end)   l2(2:end)   l3(2:end)  ];

    ne = numel(e1);
    newside = 7;

    newedges = [e1;
                e2;
                zeros(2, ne);
                newside * ones(1, ne);
                region * ones(1, ne);
                region * ones(1, ne) + 1];
    
    # Unmark nodes hanging from one side but not from the other.
    l = unique([l1 l2 l3]);
    
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
