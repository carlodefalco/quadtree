function A = bim2a_quadtree_laplacian (msh, D)

  real_elem = find (! any (msh.children));

  A = sparse (numel (msh.reduced_to_full),
              numel (msh.reduced_to_full));


  II = JJ = VV = [];
  for iel = real_elem
    
    A_loc = D(iel) * local_matrix (msh, iel);

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

function a = local_matrix (msh, iel)

  a = [ 2 -1  0 -1;
       -1  2 -1  0;
        0 -1  2 -1;
       -1  0 -1  2];

endfunction

%!demo
%! msh = msh2m_quadtree (1:3, 1:3, 1, 1:4);
%! msh = msh2m_refine_quadtree(msh, 1);
%! A = bim2a_quadtree_laplacian (msh, ones (columns (msh.t), 1));
%! dnodes = msh2m_nodes_on_sides (msh, 1);
%! A(msh.full_to_reduced(dnodes), :) = [];
%! A(:, msh.full_to_reduced(dnodes)) = [];
%! spy (A);


