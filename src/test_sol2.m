clear all;
close all;
clc;

function [P] = projector(msh)
    Nnodes = columns(msh.p);
    
    hanging = find(any(msh.hanging));
    non_hanging = find(!any(msh.hanging));
    
    II = non_hanging;
    JJ = msh.full_to_reduced(non_hanging);
    VV = ones(size(II));
    
    P = sparse(II, JJ, VV, Nnodes, numel(non_hanging));
    
    II = repmat(hanging, [2, 1])(:);
    JJ = msh.full_to_reduced(msh.hanging(:, hanging))(:);
    VV = 0.5 * ones(size(II));
    
    P += sparse(II, JJ, VV, Nnodes, numel(non_hanging));
endfunction

# Mesh definition.
msh.p = [0 0 1 1 0.5 1 0 0.5; 0 1.5 0 1.5 0 0.5 0.5 0.5];
msh.t = [1 5 8 7; 5 3 6 8; 7 6 4 2]';
msh.e = [1 5; 5 3; 3 6; 6 4; 4 2; 2 7; 7 1]';
msh.e(end+1:end+5, :) = 0;
msh.e(3:7, 1:2) = [0 0 1 0 1; 0 0 1 0 1]';
msh.e(3:7, 3:4) = [0 0 2 0 1; 0 0 2 0 1]';
msh.e(3:7, 5) = [0 0 3 0 1]';
msh.e(3:7, 6:7) = [0 0 4 0 1; 0 0 4 0 1]';

msh.children = zeros(4, 3);
msh.hanging = zeros(2, 8);
msh.hanging(:, 8) = [7 6];

msh.reduced_to_full = find (! any (msh.hanging));
tmp = 1 : columns (msh.reduced_to_full);

msh.full_to_reduced = zeros (1, columns (msh.p));
msh.full_to_reduced(1, msh.reduced_to_full) = tmp;

# Build mesh sides (only for refineable elements).
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

# Sort horizontal sides to be oriented along positive x direction.
sides_x = find(msh.orien);
sides_x_nodes = msh.sides(:, sides_x);
[~, idx] = sort(reshape(msh.p(1, sides_x_nodes), ...
                        size(sides_x_nodes)));

idx1 = sub2ind(size(msh.sides), idx(1, :), sides_x);
idx2 = sub2ind(size(msh.sides), idx(2, :), sides_x);
[msh.sides(1, sides_x), msh.sides(2, sides_x)] = ...
    deal(msh.sides(idx1), msh.sides(idx2));

# Sort vertical sides to be oriented along positive y direction.
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

hanging_nodes = find(any(msh.hanging));

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

# x = [0 1];
# y = [0 1 2];
# msh = msh2m_quadtree(x, y);
# msh = msh2m_quadtree_refine(msh, 1);

# Build global matrix.
msh.P = projector(msh);

Nnodes = columns(msh.p);
Nelements = columns(msh.t);

x = msh.p(1, :).';
y = msh.p(2, :).';

# Define parameters and exact solution.
u_ex = x .* y;

# Assemble system.
D = ones(Nelements, 1);

A = bim2a_quadtree_laplacian(msh, D);
A2 = test_sol2_projection(msh, D);

dnodes = msh2m_nodes_on_sides(msh, 1:4);

f = zeros(Nelements, 1);
g = ones(Nnodes, 1);
rhs = bim2a_quadtree_rhs(msh, f, g);

# Compute solution and error.
u = bim2a_quadtree_solve(msh, A, rhs, u_ex, dnodes);

err = norm(u - u_ex, inf);

# Save solution to file.
fclose all;
filename = "sol_1";
if (exist([filename ".vtu"], "file"))
    delete([filename ".vtu"]);
endif
fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"; u_ex, "u_ex";
                                             u - u_ex, "err"}, {}, 1);
