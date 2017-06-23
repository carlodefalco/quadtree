function [estimator] = bim2c_quadtree_pde_Kopteva_estimator(msh, u, epsilon)
    refineable_elements = find(!any(msh.children));
    Nelem = numel(refineable_elements);
    
    # Compute derivatives.
    du = bim2c_quadtree_pde_edge_gradient(msh, u);
    
    [d2u_x, d2u_y] = second_order_derivative(msh, du);
    [d3u_x, d3u_y] = third_order_derivative (msh, d2u_x, d2u_y);
    
    # Get absolute value.
    du = abs(du);
    [d2u_x, d2u_y] = deal(abs(d2u_x), abs(d2u_y));
    [d3u_x, d3u_y] = deal(abs(d3u_x), abs(d3u_y));
    
    # Compute coefficients.
    M0 = ones(1, Nelem);
    M1 = max(du(msh.ts(1:4, refineable_elements)).^2, [], 1);
    M2 = max([d2u_x(msh.t(1:4, refineable_elements));
              d2u_y(msh.t(1:4, refineable_elements))], [], 1);
    M3 = epsilon(refineable_elements) .* ...
         max([d3u_x(msh.ts(1:4, refineable_elements));
              d3u_y(msh.ts(1:4, refineable_elements))], [], 1);
    
    # Compute estimator.
    h = max(msh.size(:, refineable_elements));
    estimator = h .* max([M0; M1; M2; M3], [], 1);
endfunction

function [d2u_x, d2u_y] = second_order_derivative(msh, du)
    Nnodes = columns(msh.p);
    
    d2u_x = d2u_y = zeros(Nnodes, 1);
    
    non_hanging = find(msh.full_to_reduced);
    hanging = find(!msh.full_to_reduced);
    
    for ii = non_hanging
        ## Compute adjacent sides (those containing current node).
        [~, sides] = find(msh.sides == ii);
        
        ## Identify horizontal and vertical sides.
        sides_x = sides( msh.orien(sides));
        sides_y = sides(!msh.orien(sides));
        
        if (numel(sides_x) > 1 # Ignore boundary nodes.
            && !any(msh.hanging_sides(sides_x)))
            
            ## Determine the length of each side.
            sides_x_nodes = msh.sides(:, sides_x);
            
            hx = diff(reshape(msh.p(1, sides_x_nodes), size(sides_x_nodes))).';
            
            d2u_x(ii) = diff(du(sides_x)) / mean(hx);
        endif
        
        if (numel(sides_y) > 1 # Ignore boundary nodes.
            && !any(msh.hanging_sides(sides_y)))
            
            ## Determine the length of each side.
            sides_y_nodes = msh.sides(:, sides_y);
            
            hy = diff(reshape(msh.p(2, sides_y_nodes), size(sides_y_nodes))).';
            
            d2u_y(ii) = diff(du(sides_y)) / mean(hy);
        endif
    endfor
endfunction

function [d3u_x, d3u_y] = third_order_derivative(msh, d2u_x, d2u_y)
    d3u_x = bim2c_quadtree_pde_edge_gradient(msh, d2u_x);
    d3u_y = bim2c_quadtree_pde_edge_gradient(msh, d2u_y);
    
    # Only the tangential component is used.
    d3u_x(!msh.orien) = 0;
    d3u_y( msh.orien) = 0;
    
    # Third order derivative is zero on edges containing hanging nodes.
    hanging = find(!msh.full_to_reduced);
    sides_with_hanging = any(ismember(msh.sides, hanging));
    
    d3u_x(sides_with_hanging) = 0;
    d3u_y(sides_with_hanging) = 0;
endfunction
