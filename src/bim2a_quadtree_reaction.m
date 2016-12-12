function A = bim2a_quadtree_reaction (msh, delta)

  real_elem = find (! any (msh.children));

  A = sparse (numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));


  II = JJ = VV = [];
  for iel = real_elem
    
    A_loc = local_matrix (msh, iel, delta);

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
          
          II = [II, loci*ones(1, numel(locj))];
          JJ = [JJ, locj];
          VV = [VV, A_loc(inode, jnode)*locv];
          
        endfor
      endif
    endfor
  endfor

  A = sparse (II, JJ, VV, numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));
endfunction

function A_loc = local_matrix (msh, iel, delta)
  x = msh.p(1, msh.t(1:4, iel));
  y = msh.p(2, msh.t(1:4, iel));
  
  hx = diff(x([1, 2]));
  hy = diff(y([1, 3]));
  
  A_loc = delta(iel) * speye(4) * hx * hy / 4;
endfunction
