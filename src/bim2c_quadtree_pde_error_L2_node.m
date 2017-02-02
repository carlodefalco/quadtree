function [error_L2] = bim2c_quadtree_pde_error_L2_node(msh, u, u_ex)
    refineable_elements = find(!any(msh.children));
    
    [msh, ~, node_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    
    node_space.shp = node_space.shp(:, :, refineable_elements);
    node_space.connectivity = node_space.connectivity(:, refineable_elements);
    
    Nquad = size(node_space.shp, 1); # No. of quadrature nodes.
    
    # Evaluate u on quadrature nodes for each element.
    Nconn = size(node_space.connectivity);
    u = reshape(u(node_space.connectivity), [1, Nconn]);
    u = repmat(u, [Nquad, 1, 1]);
    u = squeeze(sum(node_space.shp .* u, 2));
    
    # Evaluate u_ex on quadrature nodes for each element.
    u_ex = reshape(u_ex(node_space.connectivity), [1, Nconn]);
    u_ex = repmat(u_ex, [Nquad, 1, 1]);
    u_ex = squeeze(sum(node_space.shp .* u_ex, 2));
    
    ## Compute norm for each mesh element.
    err = (u - u_ex).^2;
    error_L2 = sqrt(sum(err .* msh.wjacdet(:, refineable_elements), 1));
endfunction
