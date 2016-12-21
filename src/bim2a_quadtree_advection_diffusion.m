function A = bim2a_quadtree_advection_diffusion (msh, alpha, psi)

  real_elem = find (! any (msh.children));

  II = JJ = VV = zeros(1, 16 * columns(msh.t));
  idx = 1;
  
  for iel = real_elem
    coords = msh.p(:, msh.t(1:4, iel));
    
    hx = diff(coords(1, [1, 2]));
    hy = diff(coords(2, [1, 3]));
    
    ##    A_loc = local_matrix (hx, hy, alpha(iel), psi(msh.t(1:4, iel)));
    A_loc = __advdiff_local_matrix__ (hx, hy, alpha(iel), psi(msh.t(1:4, iel)));
    
    for inode = 1:4
      for jnode = 1:4
        if (! any (msh.hanging(:, msh.t(inode, iel))))
          loci = msh.full_to_reduced(msh.t(inode, iel));
        else
          loci = msh.full_to_reduced(msh.hanging(:, msh.t(inode, iel)).');
        endif
        
        if (! any (msh.hanging(:, msh.t(jnode, iel))))
          locj = msh.full_to_reduced(msh.t(jnode, iel));
        else
          locj = msh.full_to_reduced(msh.hanging(:, msh.t(jnode, iel)).');
        endif
        
        if (! any (msh.hanging(:, msh.t([inode jnode], iel))))
          locv = 1;
        else
          locv = [1/2 1/2];
        endif
        
        II(idx : (idx + numel(locv) - 1)) = loci;
        JJ(idx : (idx + numel(locv) - 1)) = locj;
        VV(idx : (idx + numel(locv) - 1)) = A_loc(inode, jnode)*locv;
        idx += numel(locv);
      endfor
    endfor
  endfor
  
  where = 1:(idx-1);
  A = sparse (II(where), JJ(where), VV(where), numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));
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
