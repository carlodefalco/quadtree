function [msh, edge_space, node_space, q2_space, rt_space] = ...
         bim2c_quadtree_mesh_properties(msh, nodes=[], weights=[])
         
  ## Check inputs.
  if (nargin > 3)
    error("bim2c_quadtree_mesh_properties: wrong number of input parameters.");
  elseif (!(isstruct(msh) && isfield(msh, "p") && isfield (msh, "t") && isfield(msh, "e")))
    error("bim2c_quadtree_mesh_properties: first input is not a valid mesh structure.");
  endif

  ## Compute mesh properties.
  [msh.t, msh.wjacdet, msh.nodes, msh.area, msh.size, msh.B] = ...
      msh2m_quadtree_geometrical_properties(msh, nodes, weights, ...
                                            "wjacdet", "nodes", "area", ...
                                            "hx_hy", "change_basis");
  
  if (nargout > 1)
    ## Compute shape functions and connectivity for edge space.
    for iel = 1 : columns(msh.t)
        [dphi_x, dphi_y] = bim2c_quadtree_eval_dfun_nedelec ...
                (msh, iel, msh.nodes(1, :, iel), msh.nodes(2, :, iel));
        
        edge_space.shp(1, :, :, iel) = dphi_x;
        edge_space.shp(2, :, :, iel) = dphi_y;
    endfor
    
    edge_space.connectivity = msh.ts;
  endif
  
  if (nargout > 2)
    ## Compute shape functions and connectivity for node space.
    for iel = 1 : columns(msh.t)
        phi = bim2c_quadtree_eval_fun ...
            (msh, iel, msh.nodes(1, :, iel), msh.nodes(2, :, iel));
        
        node_space.shp(:, :, iel) = phi;
    endfor
    
    node_space.connectivity = msh.t(1:4, :);
  endif
  
  if (nargout > 3)
    ## Compute shape functions and connectivity for Q2 space.
    for iel = 1 : columns(msh.t)
        phi = bim2c_quadtree_eval_fun_q2 ...
            (msh, iel, msh.nodes(1, :, iel), msh.nodes(2, :, iel));
        
        q2_space.shp(:, :, iel) = phi;
    endfor
    
    nel = columns(msh.t);
    
    q2_space.connectivity = [msh.t(1:4, :); msh.ts; 1:nel];
  endif
  
  if (nargout > 4)
    ## Compute shape functions and connectivity for Raviart-Thomas edge space.
    for iel = 1 : columns(msh.t)
        [dphi_x, dphi_y] = bim2c_quadtree_eval_dfun_rt ...
                (msh, iel, msh.nodes(1, :, iel), msh.nodes(2, :, iel));
        
        rt_space.shp(1, :, :, iel) = dphi_x;
        rt_space.shp(2, :, :, iel) = dphi_y;
    endfor
    
    rt_space.connectivity = msh.ts;
  endif
endfunction
