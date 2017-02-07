function [II, JJ, VV] = bim2a_quadtree_local_to_global(msh, A_loc, iel)
  II = JJ = VV = zeros(16, 1);
  idx = 1;
  
  for inode = 1:4
    if (! any (msh.hanging(:, msh.t(inode, iel))))
      loci = msh.full_to_reduced(msh.t(inode, iel));
    else
      loci = msh.full_to_reduced(msh.hanging(:, msh.t(inode, iel)).');
    endif
    
    for jnode = 1:4
      if (! any (msh.hanging(:, msh.t(jnode, iel))))
        locj = msh.full_to_reduced(msh.t(jnode, iel));
      else
        locj = msh.full_to_reduced(msh.hanging(:, msh.t(jnode, iel)).');
      endif
      
      hanging_ij = msh.hanging(:, msh.t([inode jnode], iel));
      
      if (!any(hanging_ij))
        # Neither i nor j are hanging.
        ii = loci;
        jj = locj;
        
        vv = 1;
      elseif (!all(any(hanging_ij)))
        # Either i or j is hanging.
        ii = loci;
        jj = locj;
        
        vv = [1/2 1/2];
      else
        # Both i and j are hanging.
        [idx_ii, idx_jj] = meshgrid(loci, locj);
        
        ii = idx_ii(:);
        jj = idx_jj(:);
        
        vv = [1/4 1/4 1/4 1/4];
      endif
      
      II(idx : (idx + numel(vv) - 1)) = ii;
      JJ(idx : (idx + numel(vv) - 1)) = jj;
      VV(idx : (idx + numel(vv) - 1)) = A_loc(inode, jnode) * vv;
      idx += numel(vv);
    endfor
  endfor
  
  II = II(1 : (idx-1));
  JJ = JJ(1 : (idx-1));
  VV = VV(1 : (idx-1));
endfunction
