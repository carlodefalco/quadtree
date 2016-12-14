function A = bim2a_quadtree_advection_diffusion (msh, alpha, psi)

  real_elem = find (! any (msh.children));

  II = JJ = VV = zeros(1, 16 * columns(msh.t));
  idx = 1;
  
  for iel = real_elem
    A_loc = local_matrix (msh, iel, alpha, psi);

    for inode = 1:4
      if (! any (msh.hanging(:, msh.t(inode, iel))))
        loci = msh.full_to_reduced(msh.t(inode, iel));
        for jnode = 1:4
          if (! any (msh.hanging(:, msh.t(jnode, iel))))
            locj = msh.full_to_reduced(msh.t(jnode, iel));
            locv = 1;
          else
            locj = msh.full_to_reduced(msh.hanging(:, msh.t(jnode, iel)).');
            locv = [1/2 1/2];
          end
          
          II(idx : (idx + numel(locj) - 1)) = loci;
          JJ(idx : (idx + numel(locj) - 1)) = locj;
          VV(idx : (idx + numel(locj) - 1)) = A_loc(inode, jnode)*locv;
          idx += numel(locj);
        endfor
      endif
    endfor
  endfor
  
  where = 1:(idx-1);
  A = sparse (II(where), JJ(where), VV(where), numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));
endfunction

function A_loc = local_matrix (msh, iel, alpha, psi)
  x = msh.p(1, msh.t(1:4, iel));
  y = msh.p(2, msh.t(1:4, iel));
  
  hx = diff(x([1, 2]));
  hy = diff(y([1, 3]));
  
  psiloc = psi(msh.t(1:4, iel));
  psi12 = psiloc(2) - psiloc(1);
  psi23 = psiloc(3) - psiloc(2);
  psi34 = psiloc(4) - psiloc(3);
  psi41 = psiloc(1) - psiloc(4);
  
  [bp12, bm12] = bimu_bernoulli(psi12);
  [bp23, bm23] = bimu_bernoulli(psi23);
  [bp34, bm34] = bimu_bernoulli(psi34);
  [bp41, bm41] = bimu_bernoulli(psi41);
  
  bp12 = alpha(iel) * bp12 * hy / (2 * hx);
  bm12 = alpha(iel) * bm12 * hy / (2 * hx);
  bp23 = alpha(iel) * bp23 * hx / (2 * hy);
  bm23 = alpha(iel) * bm23 * hx / (2 * hy);
  bp34 = alpha(iel) * bp34 * hy / (2 * hx);
  bm34 = alpha(iel) * bm34 * hy / (2 * hx);
  bp41 = alpha(iel) * bp41 * hx / (2 * hy);
  bm41 = alpha(iel) * bm41 * hx / (2 * hy);

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
