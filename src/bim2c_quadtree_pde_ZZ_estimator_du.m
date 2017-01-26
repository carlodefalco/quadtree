function [estimator] = bim2c_quadtree_pde_ZZ_estimator_du(msh, u)
    refineable_elements = find(!any(msh.children));
    
    [msh, edge_space, node_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    edge_space.shp = edge_space.shp(:, :, :, refineable_elements);
    edge_space.connectivity = edge_space.connectivity(:, refineable_elements);
    
    node_space.shp = node_space.shp(:, :, refineable_elements);
    node_space.connectivity = node_space.connectivity(:, refineable_elements);
    
    Nquad = size(edge_space.shp, 2); # No. of quadrature nodes.
    
    ## Compute gradient on edge space.
    Nconn = size(edge_space.connectivity);
    
    du_edge = bim2c_quadtree_pde_edge_gradient(msh, u);
    
    ## Compute reconstructed gradient on node space.
    Nconn = size(node_space.connectivity);
    
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du_edge);
    
    
    # Evaluate du_edge on quadrature nodes for each mesh element.
    du_edge = reshape(du_edge(edge_space.connectivity), [1, Nconn]);
    du_edge = repmat(du_edge, [Nquad, 1, 1]);
    
    du_x_edge = squeeze(sum(squeeze(edge_space.shp(1, :, :, :)) .* du_edge, 2));
    du_y_edge = squeeze(sum(squeeze(edge_space.shp(2, :, :, :)) .* du_edge, 2));
    
    # Evaluate du_node on quadrature nodes for each element.
    du_x = reshape(du_x(node_space.connectivity), [1, Nconn]);
    du_x = repmat(du_x, [Nquad, 1, 1]);
    
    du_y = reshape(du_y(node_space.connectivity), [1, Nconn]);
    du_y = repmat(du_y, [Nquad, 1, 1]);
    
    du_x_node = squeeze(sum(node_space.shp .* du_x, 2));
    du_y_node = squeeze(sum(node_space.shp .* du_y, 2));
    
    
    ## Compute error estimator for each mesh element.
    err = (du_x_edge - du_x_node).^2 + (du_y_edge - du_y_node).^2;
    estimator = sqrt(sum(err .* msh.wjacdet(:, refineable_elements), 1));
endfunction
