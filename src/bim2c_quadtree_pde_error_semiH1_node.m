function [seminorm_H1] = bim2c_quadtree_pde_error_semiH1_node(msh, du_x, du_y, du_x_ex, du_y_ex)
    refineable_elements = find(!any(msh.children));
    
    [msh, ~, node_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    
    node_space.shp = node_space.shp(:, :, refineable_elements);
    node_space.connectivity = node_space.connectivity(:, refineable_elements);
    
    Nquad = size(node_space.shp, 1); # No. of quadrature nodes.
    
    # Evaluate du_edge on quadrature nodes for each mesh element.
    Nconn = size(node_space.connectivity);
    
    du_x = reshape(du_x(node_space.connectivity), [1, Nconn]);
    du_x = repmat(du_x, [Nquad, 1, 1]);
    
    du_y = reshape(du_y(node_space.connectivity), [1, Nconn]);
    du_y = repmat(du_y, [Nquad, 1, 1]);
    
    du_x_node = squeeze(sum(node_space.shp .* du_x, 2));
    du_y_node = squeeze(sum(node_space.shp .* du_y, 2));
    
    # Evaluate exact gradient on quadrature nodes.
    du_x_ex = reshape(du_x_ex(node_space.connectivity), [1, Nconn]);
    du_x_ex = repmat(du_x_ex, [Nquad, 1, 1]);
    
    du_y_ex = reshape(du_y_ex(node_space.connectivity), [1, Nconn]);
    du_y_ex = repmat(du_y_ex, [Nquad, 1, 1]);
    
    du_x_ex = squeeze(sum(node_space.shp .* du_x_ex, 2));
    du_y_ex = squeeze(sum(node_space.shp .* du_y_ex, 2));
    
    ## Compute norm for each mesh element.
    err = (du_x_node - du_x_ex).^2 + (du_y_node - du_y_ex).^2;
    seminorm_H1 = sqrt(sum(err .* msh.wjacdet(:, refineable_elements), 1));
endfunction
