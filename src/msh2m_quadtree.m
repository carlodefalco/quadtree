function msh = msh2m_quadtree (x, y, region = 1, sides = [1:4]);

  nx = numel (x);
  ny = numel (y);
  
  [XX,YY] = meshgrid (x, y);
  msh.p = [XX(:), YY(:)].';

  iiv(ny,nx) = 0;
  iiv(:) = 1 : nx*ny;
  iiv(end,:) = [];
  iiv(:,end) = [];
  iiv = iiv(:).';

  msh.t = [iiv; iiv+ny; iiv+ny+1; iiv+1];
  msh.t(5,:) = region;

  l1 = 1 + ny * ([1:nx] - 1);
  l4 = 1:ny;
  l2 = ny * (nx-1) + 1:nx*ny;
  l3 = ny + l1 - 1;

  mr1 = [l1([1:end-1]), l2([1:end-1]), l3([1:end-1]), l4([1:end-1])];
  mr2 = [l1([2:end]),   l2([2:end]),   l3([2:end]),   l4([2:end])];

  assert (numel (mr1) == numel (mr2));
  ne = numel (mr1);
  
  msh.e = [mr1;
           mr2;
           zeros(2, ne);
           sides(1*ones(1,numel(l1)-1)), sides(2*ones(1,numel(l2)-1)), sides(3*ones(1,numel(l3)-1)), sides(4*ones(1,numel(l4)-1))
           zeros(1, ne);
           region(ones(1,ne))];
  
  s1 = sort(msh.t(1:2,    :), 1); o1 = ones (1, columns(s1)); # 1: horizontal.
  s2 = sort(msh.t(2:3,    :), 1); o2 = zeros(1, columns(s2)); # 0: vertical.
  s3 = sort(msh.t(3:4,    :), 1); o3 = ones (1, columns(s3));
  s4 = sort(msh.t([4, 1], :), 1); o4 = zeros(1, columns(s4));
  
  allsides = [s1 s2 s3 s4].';
  allorien = [o1 o2 o3 o4];
  
  [sides, ii, jj] = unique(allsides, "rows");
  msh.sides = sides.';
  msh.orien = allorien(ii);
  
  msh.ts = reshape(jj, [], 4).';
  
  msh.level    = ones  (1, columns (msh.t));
  msh.parent   = zeros (1, columns (msh.t));
  msh.children = zeros (4, columns (msh.t));
  msh.hanging  = zeros (2, columns (msh.p));
  
  msh.hanging_sides = zeros (1, columns (msh.sides));

  msh.onboundary = zeros (1, columns (msh.p));
  for ib = 1 : 4
    bnd = msh.e (1:2, msh.e(5, :) == ib);
    ee = [unique(bnd(1,:)), unique(bnd(2,:))];
    msh.onboundary (ee) = ib;
  endfor
  
  msh.reduced_to_full = 1 : columns(msh.hanging);
  tmp = 1 : columns (msh.reduced_to_full);

  msh.full_to_reduced = zeros (1, columns (msh.p));
  msh.full_to_reduced(1, msh.reduced_to_full) = tmp;
endfunction

%!test
%! x = y = 1:2;
%! msh = msh2m_quadtree (x, y, 1, 1:4);
%! assert (size (msh.p), [2 4]);
