function [msh, edge_space, node_space] = ...
         bim2c_quadtree_mesh_properties(msh, nodes=[], weights=[])
         
  ## Check inputs.
  if (nargin > 3)
    error("bim2c_quadtree_mesh_properties: wrong number of input parameters.");
  elseif (!(isstruct(msh) && isfield(msh, "p") && isfield (msh, "t") && isfield(msh, "e")))
    error("bim2c_quadtree_mesh_properties: first input is not a valid mesh structure.");
  endif

  ## Compute mesh properties.
  [msh.t, msh.wjacdet, msh.nodes, msh.area] = ...
      msh2m_quadtree_geometrical_properties(msh, nodes, weights, "wjacdet", "nodes", "area");
  
  ## Compute shape functions and connectivity for edge space.
  for iel = 1 : columns(msh.t)
    [dphi_x, dphi_y] = bim2c_quadtree_eval_dfun_nedelec ...
            (msh, iel, msh.nodes(1, :, iel), msh.nodes(2, :, iel));
    
    edge_space.shp(1, :, :, iel) = dphi_x;
    edge_space.shp(2, :, :, iel) = dphi_y;
  endfor
  
  edge_space.connectivity = msh.ts;
  
  ## Compute shape functions and connectivity for node space.
  for iel = 1 : columns(msh.t)
    phi = bim2c_quadtree_eval_fun ...
          (msh, iel, msh.nodes(1, :, iel), msh.nodes(2, :, iel));
    
    node_space.shp(:, :, iel) = phi;
  endfor
  
  node_space.connectivity = msh.t(1:4, :);
endfunction
