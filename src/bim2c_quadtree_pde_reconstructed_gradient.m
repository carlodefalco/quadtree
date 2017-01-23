function [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du)
    Nnodes = columns(msh.p);
    
    du_x = du_y = zeros(Nnodes, 1);
    
    non_hanging = find(msh.full_to_reduced);
    hanging = find(!msh.full_to_reduced);
    
    for ii = non_hanging
        ## Compute adjacent sides (those sharing a node with current element).
        [~, sides] = ind2sub(size(msh.sides), find(msh.sides == ii));
        
        ## Ignore hanging edge parent.
        sides(logical(msh.hanging_sides(sides))) = [];
        
        ## TODO.
    endfor
    
    ## Interpolate gradient at the hanging nodes.
    du_x(hanging) = mean(du_x(msh.hanging(:, hanging)));
    du_y(hanging) = mean(du_y(msh.hanging(:, hanging)));
endfunction
