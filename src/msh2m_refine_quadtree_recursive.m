function msh = msh2m_refine_quadtree_recursive(msh, refinelist)
    for ii = 1 : numel (refinelist)
        msh = do_refinement_recursive (msh, refinelist(ii));
    endfor
endfunction

function msh = do_refinement_recursive(msh, iel)
    if (! any(msh.children (:, iel)))
        if (! any(msh.hanging (:, msh.t(1:4, iel)(:))))
            msh = msh2m_refine_quadtree(msh, iel);
            return;
        endif
    endif
    
    nodes = msh.t(1:4, iel);
    
    ## Compute neighbors elements.
    # Linear indexing.
    tmp = arrayfun(@(i) find(msh.t(1:4, :) == i),
                   nodes, "UniformOutput", false);
    
    # Convert "tmp" to subscript indexing.
    [~, neighbors] = arrayfun(@(i) ind2sub([4, columns(msh.t)], tmp{i}),
                              1:4, "UniformOutput", false);
    neighbors = unique([neighbors{1}; neighbors{2}; neighbors{3}; neighbors{4}]);
    
    ## Ignore neighbors from an incompatible level.
    curr_level = msh.level(iel);
    neighbors(msh.level(neighbors) != curr_level - 1) = [];
    # Ignore parent and add current element as the last one to be refined.
    neighbors = setdiff(neighbors, [msh.parent(iel), iel]);
    neighbors(end + 1) = iel;
    
    msh = msh2m_refine_quadtree_recursive(msh, neighbors);
endfunction

%!demo
%! n = 10;
%! 
%! # Mesh definition.
%! x = linspace(0, 1, n);
%! y = linspace(0, 2, n);
%! 
%! region = 1;
%! sides = 1:4;
%! 
%! msh = msh2m_quadtree(x, y, region, sides);
%! 
%! # Try recursive refinement.
%! msh = msh2m_refine_quadtree_recursive (msh, 22);
%! msh = msh2m_refine_quadtree_recursive (msh, [84 81 103 106]);
%! quadmesh(msh, "show_cell_numbers", "show_node_numbers");
%! 
