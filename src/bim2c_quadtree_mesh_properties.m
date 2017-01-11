function [msh] = bim2c_quadtree_mesh_properties(msh)
  ## Check input.
  if nargin != 1
    error("bim2c_quadtree_mesh_properties: wrong number of input parameters.");
  elseif !(isstruct(msh) && isfield(msh, "p") && isfield (msh, "t") && isfield(msh, "e"))
    error("bim2c_quadtree_mesh_properties: first input is not a valid mesh structure.");
  endif

  [msh.wjacdet, msh.area, msh.shg] = ...
      msh2m_quatree_geometrical_properties(msh, "wjacdet", "area", "shg");
endfunction
