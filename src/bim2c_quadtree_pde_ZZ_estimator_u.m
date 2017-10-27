function [estimator] = bim2c_quadtree_pde_ZZ_estimator_u(msh, u)
    refineable_elements = find(!any(msh.children));
    
    [msh, ~, node_space, q2_space] = bim2c_quadtree_mesh_properties(msh, [], []);
    node_space.shp = node_space.shp(:, :, refineable_elements);
    node_space.connectivity = node_space.connectivity(:, refineable_elements);
    
    q2_space.shp = q2_space.shp(:, :, refineable_elements);
    q2_space.connectivity = q2_space.connectivity(:, refineable_elements);
    
    Nquad = size(node_space.shp, 1); # No. of quadrature nodes.
    
    ## Compute reconstructed solution on q2 space.
    [u_star_node, u_star_edge, u_star_center] = bim2c_quadtree_pde_reconstructed_solution(msh, u);
    
    # Evaluate u on quadrature nodes for each element.
    Nconn = size(node_space.connectivity);
    u = reshape(u(node_space.connectivity), [1, Nconn]);
    u = repmat(u, [Nquad, 1, 1]);
    
    u_node = squeeze(sum(node_space.shp .* u, 2));
    
    # Evaluate u_star on quadrature nodes for each element.
    Nconn = size(q2_space.connectivity);
    u_star = reshape([u_star_node(q2_space.connectivity(1:4, :));
                      u_star_edge(q2_space.connectivity(5:8, :));
                      u_star_center(q2_space.connectivity(9, :)).'], [1, Nconn]);
    u_star = repmat(u_star, [Nquad, 1, 1]);
    
    u_q2 = squeeze(sum(q2_space.shp .* u_star, 2));
    
    
    ## Compute error estimator for each mesh element.
    err = (u_node - u_q2).^2;
    estimator = sqrt(sum(err .* msh.wjacdet(:, refineable_elements), 1));
endfunction
