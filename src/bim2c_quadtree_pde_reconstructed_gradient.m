function [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du)
    Nnodes = columns(msh.p);
    
    du_x = du_y = zeros(Nnodes, 1);
    
    non_hanging = find(msh.full_to_reduced);
    hanging = find(!msh.full_to_reduced);
    
    for ii = non_hanging
        ## Compute adjacent sides (those containing current node).
        [~, sides] = find(msh.sides == ii);
        
        ## Ignore hanging sides.
        sides(find(msh.hanging_sides(sides))) = [];
        
        ## Identify horizontal and vertical sides.
        sides_x = sides( msh.orien(sides));
        sides_y = sides(!msh.orien(sides));
        
        ## Determine the length of each side.
        sides_x_nodes = msh.sides(:, sides_x);
        sides_y_nodes = msh.sides(:, sides_y);
        
        hx = diff(reshape(msh.p(1, sides_x_nodes), size(sides_x_nodes))).';
        hy = diff(reshape(msh.p(2, sides_y_nodes), size(sides_y_nodes))).';
        
        ## Compute the reconstructed gradient along x.
        if (numel(sides_x) == 1)
            neighbor_node = setdiff(sides_x_nodes, ii);
            [~, neighbors] = find(msh.sides == neighbor_node);
            
            # Ignore current side.
            neighbors = setdiff(neighbors, sides_x);
            
            # Ignore vertical sides.
            neighbors = neighbors(msh.orien(neighbors));
            
            # Ignore hanging sides.
            neighbors(find(msh.hanging_sides(neighbors))) = [];
            
            # Compute hx.
            neighbors_nodes = msh.sides(:, neighbors);
            hx_neighbor = diff(reshape(msh.p(1, neighbors_nodes), ...
                                       size(neighbors_nodes))).';
            
            # Add neighbor to sides_x list and compute weights.
            sides_x = [sides_x; neighbors];
            weights_x = [1/hx + 2/hx_neighbor; -1/hx_neighbor];
        else
            weights_x = 1 ./ hx;
        endif
        
        ## Compute the reconstructed gradient along y.
        if (numel(sides_y) == 1) # Boundary nodes.
            neighbor_node = setdiff(sides_y_nodes, ii);
            [~, neighbors] = find(msh.sides == neighbor_node);
            
            # Ignore current side.
            neighbors = setdiff(neighbors, sides_y);
            
            # Ignore horizontal sides.
            neighbors = neighbors(!msh.orien(neighbors));
            
            # Ignore hanging sides.
            neighbors(find(msh.hanging_sides(neighbors))) = [];
            
            # Compute hy.
            neighbor_nodes = msh.sides(:, neighbors);
            hy_neighbor = diff(reshape(msh.p(2, neighbor_nodes), ...
                                       size(neighbor_nodes))).';
            
            # Add neighbor to sides_y list and compute weights.
            sides_y = [sides_y; neighbors];
            weights_y = [1/hy + 2/hy_neighbor; -1/hy_neighbor];
        else
            weights_y = 1 ./ hy;
        endif
        
        du_x(ii) = sum(du(sides_x) .* weights_x) / sum(weights_x);
        du_y(ii) = sum(du(sides_y) .* weights_y) / sum(weights_y);
    endfor
    
    ## Interpolate gradient at the hanging nodes.
    du_x(hanging) = mean(du_x(msh.hanging(:, hanging)));
    du_y(hanging) = mean(du_y(msh.hanging(:, hanging)));
endfunction
