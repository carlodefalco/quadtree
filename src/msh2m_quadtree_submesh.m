function [omesh, nodelist, elementlist, sidelist] = msh2m_quadtree_submesh(imesh, sdl)
  ## Check input.
  if (nargin != 2)
    error("msh2m_submesh: wrong number of input parameters.");
  endif
  if (!isstruct(imesh))
    error("msh2m_submesh: first input is not a valid mesh structure.");
  endif
  if (!isvector(sdl))
    error("msh2m_submesh: second input is not a valid vector.");
  endif
  
  ## Extract sub-mesh.
  nsd = length(sdl); # Number of subdomains.

  ## Set list of output triangles.
  elementlist = [];
  for isd = 1 : nsd
    elementlist = [elementlist find(imesh.t(5, :) == sdl(isd))];
  endfor

  omesh.t = imesh.t(:, elementlist);

  ## Set list of output nodes.
  nodelist = unique(reshape(imesh.t(1:4, elementlist), 1, []));
  omesh.p  = imesh.p(:, nodelist);

  ## Use new node numbering in connectivity matrix.
  indx(nodelist) = 1 : length(nodelist);
  omesh.t(1:4, :) = indx(omesh.t(1:4, :));

  ## Set list of output edges.
  omesh.e = [];
  for isd=1:nsd
    omesh.e = [omesh.e imesh.e(:, imesh.e(7, :) == sdl(isd))];
    omesh.e = [omesh.e imesh.e(:, imesh.e(6, :) == sdl(isd))];
  endfor

  omesh.e = unique(omesh.e', "rows")';

  ## Use new node numbering in boundary segment list.
  omesh.e(1:2, :) = indx(omesh.e(1:2, :));

  ## Update quadtree mesh properties.
  omesh.level  = imesh.level (:, elementlist);
  omesh.parent = imesh.parent(:, elementlist);
  
  omesh.children = imesh.children(:, elementlist);
  
  ## Use new node numbering in children list.
  indx_elem(elementlist) = 1 : length(elementlist);
  
  where = all(omesh.children != 0);
  omesh.children(:, where) = indx_elem(omesh.children(:, where));
  
  omesh.hanging    = imesh.hanging   (:, nodelist);
  omesh.onboundary = imesh.onboundary(:, nodelist);
  
  omesh.reduced_to_full = find(!any(omesh.hanging));
  tmp = 1 : columns(omesh.reduced_to_full);

  omesh.full_to_reduced = zeros (1, columns(omesh.p));
  omesh.full_to_reduced(1, omesh.reduced_to_full) = tmp;
  
  ## Rebuild side structures.
  omesh.ts = imesh.ts(:, elementlist);
  
  sidelist = setdiff(unique(omesh.ts), 0);
  indx_sides(sidelist) = 1 : length(sidelist);
  
  where = all(omesh.ts != 0);
  omesh.ts(:, where) = indx_sides(omesh.ts(:, where));
  
  omesh.sides = indx(imesh.sides(:, sidelist));
  omesh.orien = imesh.orien(sidelist);
  
  omesh.hanging_sides = imesh.hanging_sides(sidelist);
  
  where = (omesh.hanging_sides != 0);
  omesh.hanging_sides(where) = indx_sides(omesh.hanging_sides(where));
endfunction

%!demo
%! # Define mesh.
%! n = 11;
%! 
%! x = linspace(0, 1, n);
%! y = linspace(0, 1, n);
%! 
%! region = 1;
%! sides = 1:4;
%! 
%! msh = msh2m_quadtree(x, y, region, sides);
%! 
%! # Mark elements to discard.
%! for iel = 1 : columns(msh.t)
%!     coords = msh.p(:, msh.t(1:4, iel));
%!     
%!     if (min(coords(1, :)) >= 0.5
%!         && min(coords(2, :)) >= 0.5)
%!         
%!         msh.t(5, iel) = region + 1;
%!     endif
%! endfor
%! 
%! # Mark edges to discard.
%! for iedge = 1 : columns(msh.e)
%!     coords = msh.p(:, msh.e(1:2, iedge));
%!     
%!     if (min(coords(1, :)) >= 0.5
%!         && min(coords(2, :)) >= 0.5)
%!         
%!         msh.e(7, iedge) = region + 1;
%!     endif
%! endfor
%! 
%! # Mark internal edges.
%! l1 = flip(find(msh.p(1, :) >= 0.5 & msh.p(2, :) == 0.5));
%! l2 = find(msh.p(1, :) == 0.5 & msh.p(2, :) >= 0.5);
%! 
%! e1 = [l1(1:end-1) l2(1:end-1)];
%! e2 = [l1(2:end)   l2(2:end)  ];
%! 
%! ne = numel(e1);
%! newside = 5;
%! 
%! newedges = [e1;
%!             e2;
%!             zeros(2, ne);
%!             newside * ones(1, ne);
%!             region(ones(1, ne));
%!             region(ones(1, ne))];
%! 
%! msh.e = [msh.e newedges];
%! 
%! msh.onboundary(unique([l1 l2])) = newside;
%! 
%! # Extract and show submesh.
%! msh = msh2m_quadtree_submesh(msh, region);
%! figure; quadmesh(msh, "show_cell_numbers", "show_node_numbers");
