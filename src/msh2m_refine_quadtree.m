function msh = msh2m_refine_quadtree (msh, refinelist)

  for ii = 1 : numel (refinelist)
    msh = do_refinement (msh, refinelist(ii));
  endfor
  
endfunction

function msh = do_refinement (msh, iel);
  
  nel = columns (msh.t);
  nn  = columns (msh.p);

  if (! any (msh.children (:, iel)))
    if (! any (msh.hanging (:, msh.t(1:4, iel))(:)))
      p = [];
      
      hanging(:, 1) = [msh.t(1, iel); msh.t(2, iel)];
      hanging(:, 2) = [msh.t(2, iel); msh.t(3, iel)];
      hanging(:, 3) = [msh.t(3, iel); msh.t(4, iel)];
      hanging(:, 4) = [msh.t(4, iel); msh.t(1, iel)];

      padd = 0;
      rot = [1 2 3 4 1];

      for is = 1 : 4        
        hh{is} = find ((hanging(1, is) == msh.hanging(1, :)
                        & hanging(2, is) == msh.hanging(2, :))
                       | (hanging(2, is) == msh.hanging(1, :)
                          & hanging(1, is) == msh.hanging(2, :)));
        if (isempty (hh{is}))
          ++padd;
          nni{is} = nn + padd;
          p(:, padd) =  (msh.p(:, msh.t(rot(is), iel)) +                                       
                         msh.p(:, msh.t(rot(is+1), iel)))/2;
        else
          assert (numel (hh{is}) == 1)
          hanging(:, is) = 0;
          nni{is} = hh{is};
        endif    
      endfor
      
      p(:, ++padd) = sum (msh.p(:, msh.t (1:4, iel)), 2) / 4;
      nni{5} = nn + padd;
      
      t(:, 1) = [msh.t(1, iel); nni{1}; nni{5}; nni{4}; msh.t(5, iel)];
      t(:, 2) = [nni{1}; msh.t(2, iel); nni{2}; nni{5}; msh.t(5, iel)];
      t(:, 3) = [nni{5}; nni{2}; msh.t(3, iel); nni{3}; msh.t(5, iel)];
      t(:, 4) = [nni{4}; nni{5}; nni{3}; msh.t(4, iel); msh.t(5, iel)];
      
      level  = (msh.level(:, iel) + 1) * ones (1, 4);
      parent = iel * ones (1, 4);
      
      msh.p        = cat (2, msh.p, p);
      msh.t        = cat (2, msh.t, t);
      msh.level    = cat (2, msh.level, level);
      msh.parent   = cat (2, msh.parent, parent);
      msh.onboundary =  [msh.onboundary zeros(1, columns(p))];
      
      sides = unique(msh.e(5, :));
      
      for ih = 1:4
        ihh = hanging(:, ih);
        
        if (any (ihh))
          ob  = msh.onboundary(ihh);
          internal_boundary = false;
          
          if (any (ob))
            oc = oncorner (ihh, msh.e);
            if (all (ob))
              if (oc(1) && ! oc(2))
                msh.onboundary (nni{ih}) = ob(2);
              elseif (oc(2) && ! oc(1))
                msh.onboundary (nni{ih}) = ob(1);
              else
                idx = find ((msh.e(1, :) == ihh(1)
                             | (msh.e(2, :) == ihh(1)))
                            & (msh.e(1, :) == ihh(2)
                               | (msh.e(2, :) == ihh(2))));
                if (! isempty (idx))
                  msh.onboundary (nni{ih}) = msh.e (5, idx);
                else
                  msh.onboundary (nni{ih}) = [];
                  internal_boundary = true;
                endif
              endif
              msh.e(1, end+1) = ihh(1);
              msh.e(2, end)   = nni{ih};
              msh.e(5, end)   = msh.onboundary(nni{ih});
              msh.e(1, end+1) = nni{ih};
              msh.e(2, end)   = ihh(2);
              msh.e(5, end)   = msh.onboundary(nni{ih});
            else
              msh.hanging (:, nni{ih}) = 0;
            endif
          else
            msh.onboundary (nni{ih}) = 0;
          endif
          if (all (ob) && ! internal_boundary)
            msh.hanging (:, nni{ih}) = 0;
          else
            msh.hanging (:, nni{ih}) = ihh;
          end
        else
          msh.onboundary (nni{ih}) = 0;
          msh.hanging (:, nni{ih}) = 0;
        endif
      endfor
      
      msh.hanging(:, nni{5})  = [0; 0];
      msh.children = cat (2, msh.children, zeros (4,columns(t)));
      msh.children(:, iel) = nel + [1:4] .';
    endif
  endif


  msh.reduced_to_full = find (! any (msh.hanging));
  tmp = 1 : columns (msh.reduced_to_full);

  msh.full_to_reduced = zeros (1, columns (msh.p));
  msh.full_to_reduced(1, msh.reduced_to_full) = tmp;

endfunction

function res = oncorner(ihh, e)
  res = zeros(size(ihh));
  
  for i = 1 : numel(res)
    edges = find(e(1, :) == ihh(i) | e(2, :) == ihh(i));
    if (numel(unique(e(5, edges))) != 1)
      res(i) = 1;
    endif
  endfor
endfunction

%!test
%! x = y = 1:2;
%! msh = msh2m_refine_quadtree ...
%!   (msh2m_quadtree (x, y, 1, 1:4), 1);
%! assert (size (msh.p), [2 9]);

%!demo
%! # Mesh definition.
%! x = linspace(0, 5, 10);
%! y = linspace(0, 3,  5);
%! 
%! region = 1;
%! sides = 1:4;
%! 
%! msh = msh2m_quadtree(x, y, region, sides);
%! 
%! # Mesh refinement.
%! msh = msh2m_refine_quadtree(msh, [4, 8]);
%! msh = msh2m_refine_quadtree(msh, [40 48 52]);
%! msh = msh2m_refine_quadtree(msh, [44 12 43 59]);
%! msh = msh2m_refine_quadtree(msh, [33 74 78]);
%! 
%! # Show refined mesh.
%! figure; quadmesh(msh, "show_cell_numbers", "show_node_numbers");
