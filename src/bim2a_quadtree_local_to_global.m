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
      
      if (! any (msh.hanging(:, msh.t([inode jnode], iel))))
        loci_new = loci;
        locj_new = locj;
        
        locv = 1;
      elseif (all(any(msh.hanging(:, msh.t([inode jnode], iel)))))
        [idx1, idx2] = meshgrid(loci, locj);
        
        loci_new = idx1(:);
        locj_new = idx2(:);
        
        locv = [1/4 1/4 1/4 1/4];
      else
        loci_new = loci;
        locj_new = locj;
        
        locv = [1/2 1/2];
      endif
      
      II(idx : (idx + numel(locv) - 1)) = loci_new;
      JJ(idx : (idx + numel(locv) - 1)) = locj_new;
      VV(idx : (idx + numel(locv) - 1)) = A_loc(inode, jnode) * locv;
      idx += numel(locv);
    endfor
  endfor
  
  II = II(1 : (idx-1));
  JJ = JJ(1 : (idx-1));
  VV = VV(1 : (idx-1));
endfunction
