function msh = msh2m_polymesh(msh, fractures)
    assert (numel (fractures) == 2);
    
    # Compute boundaries.
    x = [min(msh.p(1, :)) max(msh.p(1, :))];
    y = [min(msh.p(2, :)) max(msh.p(2, :))];
    
    n_elems = columns (msh.t);
    
    # Initialize elements to remove.
    to_remove = [];
    
    for el = 1 : n_elems
        x1 = msh.p(1, msh.t{el}(1));
        x2 = msh.p(1, msh.t{el}(2));
        
        y1 = msh.p(2, msh.t{el}(1));
        y2 = msh.p(2, msh.t{el}(4));
        
        # Initialize newly added nodes.
        new_nodes = [];
        
        # Find intersection with fractures.
        for s = 1 : numel (fractures)
            # Identify fracture AB.
            A = fractures{s}{1};
            B = fractures{s}{2};
            
            xA = A(1); yA = A(2);
            xB = B(1); yB = B(2);
            
            # Ignore fractures that don't intersect current element.
            if ((xA <= x1 && xB <= x1) ||
                (xA >= x2 && xB >= x2) ||
                (yA <= y1 && yB <= y1) ||
                (yA >= y2 && yB >= y2))
                continue;
            endif
            
            # If AB is not horizontal.
            if (yA != yB)
                # Bottom side.
                lambda = (y1 - yB) / (yA - yB); # Evaluate line at y=y1.
                  
                x_intersect = xA * lambda + xB * (1 - lambda); # And find corresponding x value.
                
                if (x_intersect >= x1 && x_intersect <= x2)
                    new_nodes(:, end + 1) = [x_intersect, y1, s]';
                endif
                
                # Top side.
                lambda = (y2 - yB) / (yA - yB);
                  
                x_intersect = xA * lambda + xB * (1 - lambda);
                  
                if (x_intersect >= x1 && x_intersect <= x2)
                    new_nodes(:, end + 1) = [x_intersect, y2, s]';
                endif
            endif
            
            # If AB is not vertical.
            if (xA != xB)
                # Left side.
                lambda = (x1 - xB) / (xA - xB);
                  
                y_intersect = yA * lambda + yB * (1 - lambda);
                  
                if (y_intersect >= y1 && y_intersect <= y2)
                    new_nodes(:, end + 1) = [x1, y_intersect, s]';
                endif
                
                # Right side.
                lambda = (x2 - xB) / (xA - xB);
                
                y_intersect = yA * lambda + yB * (1 - lambda);
                  
                if (y_intersect >= y1 && y_intersect <= y2)
                    new_nodes(:, end + 1) = [x2, y_intersect, s]';
                endif
                
            endif
        endfor
        
        if (isempty(new_nodes))
            continue;
        endif
        
        new_nodes = unique(new_nodes', "rows")';
        
        # Identify whether an intersection point coincides (up to a tolerance)
        # with an already existing vertex of current element.
        old_nodes = [];
        
        for k = 1 : columns(new_nodes)
            for n = 1 : 4
                if (norm(new_nodes(1:2, k) - msh.p(:, msh.t{el}(n)), inf) < 1e-15)
                    old_nodes(:, end + 1) = [k, n, new_nodes(3, k)]';
                    break;
                endif
            endfor
        endfor
        
        # Remove existing vertices from the list.
        if (! isempty(old_nodes))
            new_nodes(:, old_nodes(1, :)) = [];
        endif
        
        # If all the intersections are existing vertices, skip current element.
        if (isempty(new_nodes))
            continue;
        endif
        
        n_new = columns(new_nodes);
        
        nodes.p = [msh.p(:, msh.t{el}(1:4)), new_nodes(1:2, :)]; # Node coordinates.
        
        nodes.is_vertex = [true(1, 4), false(1, n_new)]; # 1 for vertices, 0 for intersections.
        nodes.is_new    = [false(1, 4), true(1, n_new)]; # 0 for vertices, 1 for intersections.
        
        nodes.fracture_id = [zeros(1, 4), new_nodes(3, :)]; # 0 for vertices, 1 for intersections.
        nodes.id = [msh.t{el}(1:4)', columns(msh.p) + (1:n_new)]; # Global node ids.
        
        # Mark existing nodes also as intersections.
        if (! isempty(old_nodes))
            nodes.is_new     (old_nodes(2, :)) = true;
            nodes.fracture_id(old_nodes(2, :)) = old_nodes(3, :);
        endif
        
        # Sort nodes counterclockwise (starting from west).
        [~, idx] = sort(angle(complex(center(nodes.p(1, :)),
                                      center(nodes.p(2, :)))));
        
        nodes.p           = nodes.p          (:, idx);
        nodes.is_vertex   = nodes.is_vertex  (:, idx);
        nodes.is_new      = nodes.is_new     (:, idx);
        nodes.fracture_id = nodes.fracture_id(:, idx);
        nodes.id          = nodes.id         (:, idx);
        
        # Add new vertices and elements.
        n_new += 4;
        new_id = nodes.id(:, !nodes.is_vertex & nodes.is_new);
        new_p  = nodes.p(:, !nodes.is_vertex & nodes.is_new);
        new_t  = {};
        
        for k = 1 : n_new
            # Add actual node.
            new_t{end + 1} = nodes.id(k);
            
            stop = false;
            
            next = mod(k, n_new) + 1;
            goto_sibling = false;
            
            # While the loop hasn't come back to the actual node.
            while (next != k)
              # Add next node.
              new_t{end}(end + 1, 1) = nodes.id(next);
              
              # Find next node to add.
              if ((nodes.is_vertex(next) && !nodes.is_new(next)) ||
                  goto_sibling)
                  next = mod(next, n_new) + 1;
                  goto_sibling = false;
              else
                  # Find sibling.
                  next = setdiff( find(nodes.fracture_id == nodes.fracture_id(next)), next );
                  goto_sibling = true;
              endif
            endwhile
            
            # Shift nodes so that the first one is the one with the lowest id.
            [~, idx] = min(new_t{end});
            new_t{end} = circshift(new_t{end}, -(idx - 1));
            
            # Add region index.
            new_t{end}(end + 1, 1) = 1;
        endfor
        
        # Extract unique values in new_t.
        new_t_unique = {};
        
        for k = 1 : numel(new_t)
            to_add = true;
            
            for n = 1 : numel(new_t_unique)
                if ((size(new_t{k}) == size(new_t_unique{n})) &&
                    (new_t{k} == new_t_unique{n}))
                    to_add = false;
                    break;
                endif
            endfor
            
            if (to_add)
                new_t_unique{end + 1} = new_t{k};
            endif
        endfor
        
        # Update mesh.
        msh.p(:, new_id) = new_p;
        to_remove(end + 1) = el;
        msh.t(end + 1 : end + numel(new_t_unique)) = new_t_unique;
    endfor
    
    msh.t(to_remove) = [];
    
    n_elems = columns(msh.t);
    
    # Delete duplicate nodes.
    [msh.p, r, c] = unique(msh.p', "rows");
    msh.p = msh.p';
    msh.t = cellfun(@(x) [c(x(1:end-1)); x(end)], msh.t, "UniformOutput", false);
    
    # Compute region indices and find boundary nodes.
    msh.e = [];
    
    for el = 1 : n_elems
        nodes = msh.t{el}([1:end-1 1]);
        
        # Compute region indices.
        A1 = fractures{1}{1};
        B1 = fractures{1}{2};
        
        xA1 = A1(1); yA1 = A1(2);
        xB1 = B1(1); yB1 = B1(2);
        
        A2 = fractures{2}{1};
        B2 = fractures{2}{2};
        
        xA2 = A2(1); yA2 = A2(2);
        xB2 = B2(1); yB2 = B2(2);
        
        frac1 = (yB1 - yA1) / (xB1 - xA1) * (msh.p(1, nodes) - xA1) + yA1;
        frac2 = (yB2 - yA2) / (xB2 - xA2) * (msh.p(1, nodes) - xA2) + yA2;
        
        if (all(msh.p(2, nodes) <= frac1 & msh.p(2, nodes) <= frac2))
            msh.t{el}(end) = 1;
        elseif (all(msh.p(2, nodes) >= frac1 & msh.p(2, nodes) <= frac2) ||
                 all(msh.p(2, nodes) <= frac1 & msh.p(2, nodes) >= frac2))
            msh.t{el}(end) = 2;
        elseif (all(msh.p(2, nodes) >= frac1 & msh.p(2, nodes) >= frac2))
            msh.t{el}(end) = 3;
        endif
        
        # Add boundary nodes.
        for n = 1 : numel(nodes) - 1
            x1 = msh.p(1, nodes(n));
            x2 = msh.p(1, nodes(n + 1));
            
            y1 = msh.p(2, nodes(n));
            y2 = msh.p(2, nodes(n + 1));
            
            if (y1 == y(1) && y2 == y(1))
                msh.e(:, end + 1) = [nodes(n); nodes(n + 1); 1];
            end
            
            if (x1 == x(2) && x2 == x(2))
                msh.e(:, end + 1) = [nodes(n); nodes(n + 1); 2];
            end
            
            if (y1 == y(2) && y2 == y(2))
                msh.e(:, end + 1) = [nodes(n); nodes(n + 1); 3];
            end
            
            if (x1 == x(1) && x2 == x(1))
                msh.e(:, end + 1) = [nodes(n); nodes(n + 1); 4];
            end
        endfor
    endfor
endfunction
