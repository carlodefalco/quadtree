function [error_L2] = bim2c_quadtree_pde_error_L2_node_q2(msh, u, u_ex)
    refineable_elements = find(!any(msh.children));
    
    [msh, ~, node_space, q2_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    
    node_space.shp = node_space.shp(:, :, refineable_elements);
    node_space.connectivity = node_space.connectivity(:, refineable_elements);
    
    Nquad = size(node_space.shp, 1); # No. of quadrature nodes.
    
    # Evaluate u_star on quadrature nodes for each element.
    Nconn = size(q2_space.connectivity);
    
    [u_star_node, u_star_edge, u_star_center] = bim2c_quadtree_pde_reconstructed_solution(msh, u);
    u_star = reshape([u_star_node(q2_space.connectivity(1:4, :));
                      u_star_edge(q2_space.connectivity(5:8, :));
                      u_star_center(q2_space.connectivity(9, :)).'], [1, Nconn]);
    u_star = repmat(u_star, [Nquad, 1, 1]);
    
    u_q2 = squeeze(sum(q2_space.shp .* u_star, 2));
    
    # Evaluate u_ex on quadrature nodes for each element.
    Nconn = size(node_space.connectivity);
    u_ex = reshape(u_ex(node_space.connectivity), [1, Nconn]);
    u_ex = repmat(u_ex, [Nquad, 1, 1]);
    u_ex = squeeze(sum(node_space.shp .* u_ex, 2));
    
    ## Compute norm for each mesh element.
    err = (u_q2 - u_ex).^2;
    error_L2 = sqrt(sum(err .* msh.wjacdet(:, refineable_elements), 1));
endfunction
