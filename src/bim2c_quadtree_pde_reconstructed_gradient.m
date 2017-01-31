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
        
        h_x = abs(diff(reshape(msh.p(1, sides_x_nodes), size(sides_x_nodes)))).';
        h_y = abs(diff(reshape(msh.p(2, sides_y_nodes), size(sides_y_nodes)))).';
        
        ## Compute the reconstructed gradient.
        du_x(ii) = sum(du(sides_x) ./ h_x) / sum(1 ./ h_x);
        du_y(ii) = sum(du(sides_y) ./ h_y) / sum(1 ./ h_y);
    endfor
    
    ## Interpolate gradient at the hanging nodes.
    du_x(hanging) = mean(du_x(msh.hanging(:, hanging)));
    du_y(hanging) = mean(du_y(msh.hanging(:, hanging)));
endfunction
