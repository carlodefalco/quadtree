function A = bim2a_quadtree_reaction (msh, delta, zeta)
  real_elem = find (! any (msh.children));
  
  nprocs = 4;
  
  # Serial process if mesh has few elements.
  if (columns(msh.t) < 100)
    nprocs = 1;
  endif
  
  [II, JJ, VV] = parcellfun(nprocs, @(iel) local_to_global(msh, delta, zeta, iel),
                            num2cell(real_elem), 'UniformOutput', false,
                            'VerboseLevel', 0);
  
  II = vertcat(II{:});
  JJ = vertcat(JJ{:});
  VV = vertcat(VV{:});
  
  A = sparse (II, JJ, VV, numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));
endfunction

function [II, JJ, VV] = local_to_global (msh, delta, zeta, iel)
  coords = msh.p(:, msh.t(1:4, iel));
  
  hx = diff(coords(1, [1, 2]));
  hy = diff(coords(2, [1, 3]));
  
  ## A_loc = local_matrix (hx, hy, delta(iel), zeta(msh.t(1:4, iel)));
  A_loc = __reaction_local_matrix__ (hx, hy, delta(iel), zeta(msh.t(1:4, iel)));
  
  [II, JJ, VV] = bim2a_quadtree_local_to_global(msh, A_loc, iel);
endfunction

function A_loc = local_matrix (hx, hy, delta_loc, zeta_loc)
  A_loc = delta_loc * sparse(diag(zeta_loc)) * hx * hy / 4;
endfunction
