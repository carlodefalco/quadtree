function [msh] = mark_contacts(msh)
    region = 2;
    
    # Mark elements in (region + 1).
    for iel = 1 : columns(msh.t)
        coords = msh.p(:, msh.t(1:4, iel));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        y_max = max(coords(2, :));
        
        if (# Source.
            (x_min >= msh.dim.x_source_min &&
             x_max <= msh.dim.x_source_max &&
             y_max <= msh.dim.y_contact) ||
            # Drain.
            (x_min >= msh.dim.x_drain_min &&
             x_max <= msh.dim.x_drain_max &&
             y_max <= msh.dim.y_contact))
            
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
            (x_min >= msh.dim.x_source_min &&
             x_max <= msh.dim.x_source_max &&
             y_max <= msh.dim.y_contact) ||
            # Drain.
            (x_min >= msh.dim.x_drain_min &&
             x_max <= msh.dim.x_drain_max &&
             y_max <= msh.dim.y_contact))
            
            msh.e(7, iedge) = region + 1;
        endif
    endfor

    # Mark internal edges.
    # Source left.
    l1 = find(msh.p(1, :) == msh.dim.x_source_min &
              msh.p(2, :) <= msh.dim.y_contact);
    [~, idx] = sort(msh.p(1, l1));
    l1 = l1(idx);
    
    # Source right.
    l2 = find(msh.p(1, :) == msh.dim.x_source_max &
              msh.p(2, :) <= msh.dim.y_contact);
    [~, idx] = sort(msh.p(1, l2));
    l2 = l2(idx);
    
    # Source top.
    l3 = find(msh.p(1, :) >= msh.dim.x_source_min &
              msh.p(1, :) <= msh.dim.x_source_max &
              msh.p(2, :) == msh.dim.y_contact);
    [~, idx] = sort(msh.p(1, l3));
    l3 = l3(idx);
    
    # Drain left.
    l4 = find(msh.p(1, :) == msh.dim.x_drain_min &
              msh.p(2, :) <= msh.dim.y_contact);
    [~, idx] = sort(msh.p(1, l4));
    l4 = l4(idx);
    
    # Drain right.
    l5 = find(msh.p(1, :) == msh.dim.x_drain_max &
              msh.p(2, :) <= msh.dim.y_contact);
    [~, idx] = sort(msh.p(1, l5));
    l5 = l5(idx);
    
    # Drain top.
    l6 = find(msh.p(1, :) >= msh.dim.x_drain_min &
              msh.p(1, :) <= msh.dim.x_drain_max &
              msh.p(2, :) == msh.dim.y_contact);
    [~, idx] = sort(msh.p(1, l6));
    l6 = l6(idx);
    
    # Source.
    s1 = [l1(1:end-1) l2(1:end-1) l3(1:end-1)];
    s2 = [l1(2:end)   l2(2:end)   l3(2:end)  ];

    ne = numel(s1);
    newside_s = 5;

    newedges_s = [s1;
                  s2;
                  zeros(2, ne);
                  newside_s * ones(1, ne);
                  region * ones(1, ne) - 1;
                  region * ones(1, ne) + 1];
    
    # Drain.
    d1 = [l4(1:end-1) l5(1:end-1) l6(1:end-1)];
    d2 = [l4(2:end)   l5(2:end)   l6(2:end)  ];

    ne = numel(d1);
    newside_d = 6;

    newedges_d = [d1;
                  d2;
                  zeros(2, ne);
                  newside_d * ones(1, ne);
                  region * ones(1, ne);
                  region * ones(1, ne) + 1];
    
    newedges = [newedges_s newedges_d];
    
    # Unmark nodes hanging from one side but not from the other.
    l_s = unique([l1 l2 l3]);
    l_d = unique([l4 l5 l6]);
    
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

    msh.onboundary(l_s) = newside_s;
    msh.onboundary(l_d) = newside_d;
endfunction
