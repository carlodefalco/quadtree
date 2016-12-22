function A = bim2a_quadtree_advection_diffusion (msh, alpha, psi)
  real_elem = find (! any (msh.children));
  
  nprocs = 4;
  
  # Serial process if mesh has few elements.
  if (columns(msh.t) < 100)
    nprocs = 1;
  endif
  
  [II, JJ, VV] = parcellfun(nprocs, @(iel) local_to_global(msh, alpha, psi, iel),
                            num2cell(real_elem), 'UniformOutput', false,
                            'VerboseLevel', 0);
  
  II = vertcat(II{:});
  JJ = vertcat(JJ{:});
  VV = vertcat(VV{:});
  
  A = sparse (II, JJ, VV, numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));
endfunction

function [II, JJ, VV] = local_to_global (msh, alpha, psi, iel)
  coords = msh.p(:, msh.t(1:4, iel));
  
  hx = diff(coords(1, [1, 2]));
  hy = diff(coords(2, [1, 3]));
  
  ## A_loc = local_matrix (hx, hy, alpha(iel), psi(msh.t(1:4, iel)));
  A_loc = __advdiff_local_matrix__ (hx, hy, alpha(iel), psi(msh.t(1:4, iel)));
  
  [II, JJ, VV] = bim2a_quadtree_local_to_global(msh, A_loc, iel);
endfunction

function A_loc = local_matrix (hx, hy, alpha_loc, psi_loc)
  psi12 = psi_loc(2) - psi_loc(1);
  psi23 = psi_loc(3) - psi_loc(2);
  psi34 = psi_loc(4) - psi_loc(3);
  psi41 = psi_loc(1) - psi_loc(4);
  
  [bp12, bm12] = bimu_bernoulli(psi12);
  [bp23, bm23] = bimu_bernoulli(psi23);
  [bp34, bm34] = bimu_bernoulli(psi34);
  [bp41, bm41] = bimu_bernoulli(psi41);
  
  bp12 = alpha_loc * bp12 * hy / (2 * hx);
  bm12 = alpha_loc * bm12 * hy / (2 * hx);
  bp23 = alpha_loc * bp23 * hx / (2 * hy);
  bm23 = alpha_loc * bm23 * hx / (2 * hy);
  bp34 = alpha_loc * bp34 * hy / (2 * hx);
  bm34 = alpha_loc * bm34 * hy / (2 * hx);
  bp41 = alpha_loc * bp41 * hx / (2 * hy);
  bm41 = alpha_loc * bm41 * hx / (2 * hy);

  A_loc(1, 1) = bm12 + bp41;
  A_loc(1, 2) = -bp12;
  A_loc(1, 4) = -bm41;
  
  A_loc(2, 1) = -bm12;
  A_loc(2, 2) = bp12 + bm23;
  A_loc(2, 3) = -bp23;
  
  A_loc(3, 2) = -bm23;
  A_loc(3, 3) = bp23 + bm34;
  A_loc(3, 4) = -bp34;
  
  A_loc(4, 1) = -bp41;
  A_loc(4, 3) = -bm34;
  A_loc(4, 4) = bp34 + bm41;
endfunction
