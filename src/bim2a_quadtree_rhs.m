function rhs = bim2a_quadtree_rhs (msh, f, g)
  real_elem = find (! any (msh.children));

  rhs = sparse (numel (msh.reduced_to_full), 1);

  II = JJ = VV = [];
  for iel = real_elem
    rhs_loc = local_rhs (msh, iel, f, g);

    for inode = 1:4
      if (! any (msh.hanging(:, msh.t(inode, iel))))
        loci = msh.full_to_reduced(msh.t(inode, iel));
        
        II = [II, loci];
        VV = [VV, rhs_loc(inode)];
      endif
    endfor
  endfor

  rhs = sparse (II, 1, VV, numel (msh.reduced_to_full), 1);
endfunction

function rhs_loc = local_rhs (msh, iel, f, g)
  x = msh.p(1, msh.t(1:4, iel));
  y = msh.p(2, msh.t(1:4, iel));
  
  hx = diff(x([1, 2]));
  hy = diff(y([1, 3]));
  
  rhs_loc = f(iel) * g(msh.t(1:4, iel)) * hx * hy / 4;
endfunction
