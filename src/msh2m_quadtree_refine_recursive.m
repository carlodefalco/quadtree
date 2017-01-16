function msh = msh2m_quadtree_refine_recursive(msh, refinelist)
    for ii = 1 : numel(refinelist)
        msh = do_refinement_recursive(msh, refinelist(ii));
    endfor
endfunction

function msh = do_refinement_recursive(msh, iel)
    nodes = msh.t(1:4, iel);
    hanging_nodes = msh.hanging(:, nodes);
    
    if (!any(hanging_nodes(:)))
        msh = msh2m_quadtree_refine(msh, iel);
        return;
    endif
    
    ## Compute neighbor elements (those sharing a node with current element).
    # Linear indexing.
    tmp = arrayfun(@(i) find(msh.t(1:4, :) == nodes(i)),
                   1:numel(nodes), "UniformOutput", false);
    
    # Convert "tmp" to subscript indexing.
    [~, neighbors] = ind2sub([4, columns(msh.t)], vertcat(tmp{:}));
    neighbors = unique(neighbors);
    
    ## Ignore neighbors from an incompatible level.
    neighbors(msh.level(neighbors) != msh.level(iel) - 1) = [];
    
    # Ignore current element and its parent.
    neighbors = setdiff(neighbors, [iel, msh.parent(iel)]);
    
    # Ignore neighbors not containing current element hanging nodes.
    is_hanging_neighbor = false(size(neighbors));
    
    for i = 1 : numel(neighbors)
        for j = 1 : columns(hanging_nodes)
            # Mark i-th neighbor if it contains j-th hanging node.
            if (sum(ismember(hanging_nodes(:, j), msh.t(1:4, neighbors(i)))) == 2)
                is_hanging_neighbor(i) = true;
                break;
            endif
        endfor
    endfor
    
    neighbors = neighbors(is_hanging_neighbor);
    
    ## Refine all neighbors.
    msh = msh2m_quadtree_refine_recursive(msh, neighbors);
    
    ## Refine current element.
    msh = msh2m_quadtree_refine_recursive(msh, iel);
endfunction

%!demo
%! # Mesh definition.
%! 
%! n = 10;
%! 
%! x = linspace(0, 1, n);
%! y = linspace(0, 2, n);
%! 
%! region = 1;
%! sides = 1:4;
%! 
%! msh = msh2m_quadtree(x, y, region, sides);
%! 
%! # Try recursive refinement.
%! msh = msh2m_quadtree_refine_recursive (msh, 22);
%! msh = msh2m_quadtree_refine_recursive (msh, [84 81 98 106]);
%! quadmesh(msh, "show_cell_numbers", "show_node_numbers");
%! 
