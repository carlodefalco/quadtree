function [msh] = mark_contacts(msh)
    region = 2;
    
    # Mark elements in (region + 1).
    for iel = 1 : columns(msh.t)
        coords = msh.p(:, msh.t(1:4, iel));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        y_max = max(coords(2, :));
        
        if (# Source.
            (x_max <= msh.dim.x_source &&
             y_max <= msh.dim.y_sc) ||
            # Drain.
            (x_min >= msh.dim.x_drain &&
             y_max <= msh.dim.y_sc))
            
            msh.t(5, iel) = region + 1;
        endif
    endfor

    # Mark edges in (region + 1).
    for iedge = 1 : columns(msh.e)
        coords = msh.p(:, msh.e(1:2, iedge));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        y_max = max(coords(2, :));
        
        if (# Source.
            (x_max <= msh.dim.x_source &&
             y_max <= msh.dim.y_sc) ||
            # Drain.
            (x_min >= msh.dim.x_drain &&
             y_max <= msh.dim.y_sc))
            
            msh.e(7, iedge) = region + 1;
        endif
    endfor

    # Mark internal edges.
    # Source top.
    l1 = find(msh.p(1, :) <= msh.dim.x_source &
              msh.p(2, :) == msh.dim.y_sc);
    [~, idx] = sort(msh.p(1, l1));
    l1 = l1(idx);
    
    s1 = [l1(1:end-1)];
    s2 = [l1(2:end)  ];

    ne = numel(s1);
    newside_s1 = 5;

    newedges_s1 = [s1;
                   s2;
                   zeros(2, ne);
                   newside_s1 * ones(1, ne);
                   region * ones(1, ne) - 1;
                   region * ones(1, ne) + 1];
    
    # Source right.
    l2 = find(msh.p(1, :) == msh.dim.x_source &
              msh.p(2, :) <= msh.dim.y_sc);
    [~, idx] = sort(msh.p(1, l2));
    l2 = l2(idx);
    
    s1 = [l2(1:end-1)];
    s2 = [l2(2:end)  ];

    ne = numel(s1);
    newside_s2 = 6;

    newedges_s2 = [s1;
                   s2;
                   zeros(2, ne);
                   newside_s2 * ones(1, ne);
                   region * ones(1, ne) - 1;
                   region * ones(1, ne) + 1];
    
    # Drain left.
    l3 = find(msh.p(1, :) == msh.dim.x_drain &
              msh.p(2, :) <= msh.dim.y_sc);
    [~, idx] = sort(msh.p(1, l3));
    l3 = l3(idx);
    
    d1 = [l3(1:end-1)];
    d2 = [l3(2:end)  ];

    ne = numel(d1);
    newside_d1 = 7;

    newedges_d1 = [d1;
                   d2;
                   zeros(2, ne);
                   newside_d1 * ones(1, ne);
                   region * ones(1, ne) - 1;
                   region * ones(1, ne) + 1];
    
    # Drain top.
    l4 = find(msh.p(1, :) >= msh.dim.x_drain &
              msh.p(2, :) == msh.dim.y_sc);
    [~, idx] = sort(msh.p(1, l4));
    l4 = l4(idx);
    
    d1 = [l4(1:end-1)];
    d2 = [l4(2:end)  ];

    ne = numel(d1);
    newside_d2 = 8;

    newedges_d2 = [d1;
                   d2;
                   zeros(2, ne);
                   newside_d2 * ones(1, ne);
                   region * ones(1, ne) - 1;
                   region * ones(1, ne) + 1];
    
    newedges = [newedges_s1 newedges_s2 newedges_d1 newedges_d2];
    
    # Unmark nodes hanging from one side but not from the other.
    l_s = unique([l1 l2]);
    l_d = unique([l3 l4]);
    
    l = union(l_s, l_d);
    
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

    msh.onboundary(l1) = newside_s1;
    msh.onboundary(l2) = newside_s2;
    msh.onboundary(l3) = newside_d1;
    msh.onboundary(l4) = newside_d2;
endfunction
