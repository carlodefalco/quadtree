function [msh] = mark_regions(msh)
    region = 1;
    
    # Mark elements in (region + 1).
    for iel = 1 : columns(msh.t)
        coords = msh.p(:, msh.t(1:4, iel));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        y_min = min(coords(2, :));
        
        if ((x_max <= msh.l + msh.L/4 ||
             x_min >= msh.r - msh.L/4) &&
            y_min >= 0)
            
            msh.t(5, iel) = region + 2;
        elseif (y_min >= 0)
            msh.t(5, iel) = region + 1;
        endif
    endfor

    # Mark edges in (region + 1).
    for iedge = 1 : columns(msh.e)
        coords = msh.p(:, msh.e(1:2, iedge));
        
        x_min = min(coords(1, :));
        x_max = max(coords(1, :));
        y_min = min(coords(2, :));
        
        if ((x_max <= msh.l + msh.L/4 ||
             x_min >= msh.r - msh.L/4) &&
            y_min >= 0)
            
            msh.e(7, iedge) = region + 2;
        elseif (y_min >= 0)
            msh.e(7, iedge) = region + 1;
        endif
    endfor

    # Mark internal edges.
    l1 = find(msh.p(1, :) <= msh.l + msh.L/4 & 
              msh.p(2, :) == 0);
    [~, idx] = sort(msh.p(1, l1));
    l1 = l1(idx);

    e1 = [l1(1:end-1)];
    e2 = [l1(2:end)  ];

    ne = numel(e1);
    newside = 6;

    newedges1 = [e1;
                 e2;
                 zeros(2, ne);
                 newside * ones(1, ne);
                 region * ones(1, ne);
                 region * ones(1, ne) + 2];

    msh.onboundary(l1) = newside;
    
    
    l2 = find(msh.p(1, :) == msh.l + msh.L/4 & 
              msh.p(2, :) >= 0);
    [~, idx] = sort(msh.p(2, l2));
    l2 = l2(idx);

    e1 = [l2(1:end-1)];
    e2 = [l2(2:end)  ];

    ne = numel(e1);
    newside = 8;

    newedges2 = [e1;
                 e2;
                 zeros(2, ne);
                 newside * ones(1, ne);
                 region * ones(1, ne) + 1;
                 region * ones(1, ne) + 2];

    msh.onboundary(l2) = newside;
    
    
    l3 = find(msh.p(1, :) == msh.r - msh.L/4 & 
              msh.p(2, :) >= 0);
    [~, idx] = sort(msh.p(2, l3));
    l3 = l3(idx);

    e1 = [l3(1:end-1)];
    e2 = [l3(2:end)  ];

    ne = numel(e1);
    newside = 8;

    newedges3 = [e1;
                 e2;
                 zeros(2, ne);
                 newside * ones(1, ne);
                 region * ones(1, ne) + 1;
                 region * ones(1, ne) + 2];

    msh.onboundary(l3) = newside;
    
    
    l4 = find(msh.p(1, :) >= msh.r - msh.L/4 & 
              msh.p(2, :) == 0);
    [~, idx] = sort(msh.p(1, l4));
    l4 = l4(idx);

    e1 = [l4(1:end-1)];
    e2 = [l4(2:end)  ];

    ne = numel(e1);
    newside = 5;

    newedges4 = [e1;
                 e2;
                 zeros(2, ne);
                 newside * ones(1, ne);
                 region * ones(1, ne) + 1;
                 region * ones(1, ne) + 2];

    msh.onboundary(l4) = newside;
    
    
    l5 = find(msh.p(1, :) >= msh.l + msh.L/4 & 
              msh.p(1, :) <= msh.r - msh.L/4 & 
              msh.p(2, :) == 0);
    [~, idx] = sort(msh.p(1, l5));
    l5 = l5(idx);

    e1 = [l5(1:end-1)];
    e2 = [l5(2:end)  ];

    ne = numel(e1);
    newside = 7;

    newedges5 = [e1;
                 e2;
                 zeros(2, ne);
                 newside * ones(1, ne);
                 region * ones(1, ne);
                 region * ones(1, ne) + 1];

    msh.onboundary(l5) = newside;
    

    # Unmark nodes hanging from one side but not from the other.
    newedges = [newedges1, newedges2, newedges3, newedges4, newedges5];
    l = unique([l1, l2, l3, l4, l5]);
    
    hanging_nodes = find(any(msh.hanging(:, l)));
    for i = 1 : numel(hanging_nodes)
        idx = hanging_nodes(i);
        
        elements = find(any(msh.t(1:4, :) == l(idx)));
        regions = unique(msh.t(5, elements));
        
        if (isscalar(regions))
            switch (regions)
                case region
                    row = 7;
                case region + 1
                    row = 6;
                case region + 2
                    row = 6;
            endswitch
            
            cols = find(any(newedges(1:2, :) == l(idx)));
            
            newedges(row, cols) = 0;
        endif
    endfor
    
    msh.e = [msh.e newedges];
endfunction
