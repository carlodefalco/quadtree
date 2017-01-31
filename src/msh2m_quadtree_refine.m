function msh = msh2m_quadtree_refine(msh, refinelist)
    for ii = 1 : numel(refinelist)
        msh = do_refinement_recursive(msh, refinelist(ii));
    endfor
    
    # Rebuild mesh sides (only for refineable elements).
    refineable_elements = find(!any(msh.children));
    
    s1 = sort(msh.t(1:2,    refineable_elements), 1);
    s2 = sort(msh.t(2:3,    refineable_elements), 1);
    s3 = sort(msh.t(3:4,    refineable_elements), 1);
    s4 = sort(msh.t([4, 1], refineable_elements), 1);
    
    o1 = true (1, columns(s1)); # true : horizontal.
    o2 = false(1, columns(s2)); # false: vertical.
    o3 = true (1, columns(s3));
    o4 = false(1, columns(s4));
    
    allsides = [s1 s2 s3 s4].';
    allorien = [o1 o2 o3 o4];
    
    [sides, ii, jj] = unique(allsides, "rows");
    msh.sides = sides.';
    msh.orien = allorien(ii);
    
    # Sort sides to be oriented along positive x and y directions.
    sides_x = find(msh.orien);
    
    sides_x_nodes = msh.sides(:, sides_x);
    [~, idx] = sort(reshape(msh.p(1, sides_x_nodes), ...
                            size(sides_x_nodes)));
    
    idx1 = sub2ind(size(msh.sides), idx(1, :), sides_x);
    idx2 = sub2ind(size(msh.sides), idx(2, :), sides_x);
    [msh.sides(1, sides_x), msh.sides(2, sides_x)] = ...
        deal(msh.sides(idx1), msh.sides(idx2));
    
    sides_y = find(!msh.orien);
    sides_y_nodes = msh.sides(:, sides_y);
    [~, idx] = sort(reshape(msh.p(2, sides_y_nodes), ...
                            size(sides_y_nodes)));
    
    idx1 = sub2ind(size(msh.sides), idx(1, :), sides_y);
    idx2 = sub2ind(size(msh.sides), idx(2, :), sides_y);
    [msh.sides(1, sides_y), msh.sides(2, sides_y)] = ...
        deal(msh.sides(idx1), msh.sides(idx2));
    
    # Build sides connectivity.
    msh.ts = zeros(4, columns(msh.t));
    msh.ts(:, refineable_elements) = reshape(jj, [], 4).';
    
    ## Compute hanging sides.
    msh.hanging_sides = zeros(1, columns(msh.sides));
    
    for ii = 1 : columns(msh.p)
        idx = find(all(msh.sides == msh.hanging(:, ii)
                       | flip(msh.sides) == msh.hanging(:, ii)));
        
        # Determine index of the parent side.
        if (!isempty(idx))
            sidelist = find(any(ismember(msh.sides, ii)));
            
            sidelist = sidelist(msh.orien(sidelist) == msh.orien(idx));
            msh.hanging_sides(sidelist) = idx;
        endif
    endfor
endfunction

function msh = do_refinement_recursive(msh, iel)
    nodes = msh.t(1:4, iel);
    hanging_nodes = msh.hanging(:, nodes);
    
    if (!any(hanging_nodes(:)))
        msh = do_refinement(msh, iel);
        return;
    endif
    
    ## Compute neighbor elements (those sharing a node with current element).
    [~, neighbors] = arrayfun(@(i) find(msh.t(1:4, :) == nodes(i)),
                              1:numel(nodes), "UniformOutput", false);
    neighbors = unique(vertcat(neighbors{:}));
    
    ## Ignore neighbors from an incompatible level.
    neighbors(msh.level(neighbors) != msh.level(iel) - 1) = [];
    
    # Ignore current element and its parent.
    neighbors = setdiff(neighbors, [iel, msh.parent(iel)]);
    
    # Ignore neighbors not containing current element hanging nodes.
    is_hanging_neighbor = false(size(neighbors));
    
    for i = 1 : numel(neighbors)
        for j = 1 : columns(hanging_nodes)
            # Mark i-th neighbor if it contains j-th hanging node.
            if (sum(ismember(hanging_nodes(:, j), msh.t(1:4, neighbors(i)))) == 2)
                is_hanging_neighbor(i) = true;
                break;
            endif
        endfor
    endfor
    
    neighbors = neighbors(is_hanging_neighbor);
    
    ## Refine all neighbors.
    msh = msh2m_quadtree_refine(msh, neighbors);
    
    ## Refine current element.
    msh = msh2m_quadtree_refine(msh, iel);
endfunction

function msh = do_refinement (msh, iel);
  nn  = columns (msh.p);
  nel = columns (msh.t);
  ns  = columns (msh.sides);

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
      
      msh.p = cat (2, msh.p, p);
      msh.t = cat (2, msh.t, t);
      
      msh.level    = cat (2, msh.level, level);
      msh.parent   = cat (2, msh.parent, parent);
      msh.onboundary =  [msh.onboundary zeros(1, columns(p))];
      
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
%! msh = msh2m_quadtree_refine ...
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
%! msh = msh2m_quadtree_refine(msh, [4, 8]);
%! msh = msh2m_quadtree_refine(msh, [40 48 52]);
%! msh = msh2m_quadtree_refine(msh, [44 12 43 59]);
%! msh = msh2m_quadtree_refine(msh, [33 74 78]);
%! 
%! # Show refined mesh.
%! figure; quadmesh(msh, "show_cell_numbers", "show_node_numbers");


%!demo
%! # Mesh definition.
%! 
%! n = 10;
%! 
%! x = linspace(0, 1, n);
%! y = linspace(0, 2, n);
%! 
%! region = 1;
%! sides = 1:4;
%! 
%! msh = msh2m_quadtree(x, y, region, sides);
%! 
%! # Refinement.
%! msh = msh2m_quadtree_refine (msh, 22);
%! msh = msh2m_quadtree_refine (msh, [84 81 98 106]);
%! quadmesh(msh, "show_cell_numbers", "show_node_numbers");
%! 
