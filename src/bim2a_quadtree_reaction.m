function A = bim2a_quadtree_reaction (msh, delta)

  real_elem = find (! any (msh.children));

  II = JJ = VV = zeros(1, 16 * columns(msh.t));
  idx = 1;
  
  for iel = real_elem
    coords = msh.p(:, msh.t(1:4, iel));
    
    hx = diff(coords(1, [1, 2]));
    hy = diff(coords(2, [1, 3]));
    
    A_loc = local_matrix (hx, hy, delta(iel));

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

function A_loc = local_matrix (hx, hy, delta_loc)
  A_loc = delta_loc * speye(4) * hx * hy / 4;
endfunction
