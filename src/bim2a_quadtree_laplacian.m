function A = bim2a_quadtree_laplacian (msh, D)
  real_elem = find (! any (msh.children));
  
  nprocs = 4;
  
  # Serial process if mesh has few elements.
  if (columns(msh.t) < 100)
    nprocs = 1;
  endif
  
  [II, JJ, VV] = parcellfun(nprocs, @(iel) local_to_global(msh, D, iel),
                            num2cell(real_elem), 'UniformOutput', false,
                            'VerboseLevel', 0);
  
  II = vertcat(II{:});
  JJ = vertcat(JJ{:});
  VV = vertcat(VV{:});
  
  A = sparse (II, JJ, VV, numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));
endfunction

function [II, JJ, VV] = local_to_global (msh, D, iel)
  coords = msh.p(:, msh.t(1:4, iel));
  
  hx = diff(coords(1, [1, 2]));
  hy = diff(coords(2, [1, 3]));
  
  ## A_loc = local_matrix (hx, hy, D(iel));
  A_loc = __laplacian_local_matrix__ (hx, hy, D(iel));
  
  [II, JJ, VV] = bim2a_quadtree_local_to_global(msh, A_loc, iel);
endfunction

function A_loc = local_matrix (hx, hy, D_loc)
  # The diagonal elements are accounted for twice, by
  # A_loc += A_loc';
  # so they are divided by 2.
  A_loc = speye(4) * (hx^2 + hy^2) / (4 * hx * hy);
  
  A_loc(1, 2) = -hy / (2 * hx);
  A_loc(1, 4) = -hx / (2 * hy);
  A_loc(2, 3) = -hx / (2 * hy);
  A_loc(3, 4) = -hy / (2 * hx);
  A_loc += A_loc';
  
  A_loc *= D_loc;
endfunction
