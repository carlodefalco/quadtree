function rhs = bim2a_quadtree_rhs (msh, f, g)
  real_elem = find (! any (msh.children));

  II = VV = zeros(1, 4 * numel(real_elem));
  idx = 1;
  
  for iel = real_elem
    coords = msh.p(:, msh.t(1:4, iel));
    
    hx = diff(coords(1, [1, 2]));
    hy = diff(coords(2, [1, 3]));
    
    rhs_loc = local_rhs (hx, hy, f(iel), g(msh.t(1:4, iel)));

    for inode = 1:4
      if (! any (msh.hanging(:, msh.t(inode, iel))))
        loci = msh.full_to_reduced(msh.t(inode, iel));
        locv = 1;
      else
        loci = msh.full_to_reduced(msh.hanging(:, msh.t(inode, iel)).');
        locv = [1/2 1/2];
      endif
        
      II(idx : (idx + numel(locv) - 1)) = loci;
      VV(idx : (idx + numel(locv) - 1)) = rhs_loc(inode) * locv;
      idx += numel(loci);
    endfor
  endfor

  where = 1:(idx-1);
  rhs = sparse (II(where), 1, VV(where), numel (msh.reduced_to_full), 1);
endfunction

function rhs_loc = local_rhs (hx, hy, f_loc, g_loc)
  rhs_loc = f_loc * g_loc * hx * hy / 4;
endfunction
