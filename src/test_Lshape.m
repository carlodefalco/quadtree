clear all;
close all;
clc;

# Mesh definition.
n = 11;

x = linspace(0, 1, n);
y = linspace(0, 1, n);

region = 1;
sides = 1:4;

msh = msh2m_quadtree(x, y, region, sides);

# Mark elements to discard.
for iel = 1 : columns(msh.t)
    coords = msh.p(:, msh.t(1:4, iel));
    
    if (min(coords(1, :)) >= 0.5
        && min(coords(2, :)) >= 0.5)
        
        msh.t(5, iel) = region + 1;
    endif
endfor

# Mark edges to discard.
for iedge = 1 : columns(msh.e)
    coords = msh.p(:, msh.e(1:2, iedge));
    
    if (min(coords(1, :)) >= 0.5
        && min(coords(2, :)) >= 0.5)
        
        msh.e(7, iedge) = region + 1;
    endif
endfor

# Mark internal edges.
l1 = flip(find(msh.p(1, :) >= 0.5 & msh.p(2, :) == 0.5));
l2 = find(msh.p(1, :) == 0.5 & msh.p(2, :) >= 0.5);

e1 = [l1(1:end-1) l2(1:end-1)];
e2 = [l1(2:end)   l2(2:end)  ];

ne = numel(e1);
newside = 5;

newedges = [e1;
            e2;
            zeros(2, ne);
            newside * ones(1, ne);
            region(ones(1, ne));
            region(ones(1, ne))];

msh.e = [msh.e newedges];

msh.onboundary(unique([l1 l2])) = newside;

# Extract submesh.
msh = msh2m_quadtree_submesh(msh, region);

for i = 1:10
    fprintf("i = %d\n", i);
    
    # Build global matrix.
    Nnodes = columns(msh.p);
    Nelements = columns(msh.t);

    # Assemble system.
    D = @(msh) 1e-2 * ones(columns(msh.t), 1);

    A = @(msh) bim2a_quadtree_laplacian(msh, D(msh));
    
    f = @(msh) ones(columns(msh.t), 1);
    g = @(msh) ones(columns(msh.p), 1);
    
    rhs = @(msh) bim2a_quadtree_rhs(msh, f(msh), g(msh));

    # Compute solution and error.
    dnodes = msh2m_nodes_on_sides(msh, 1:5); # Dirichlet nodes.
    
    u0 = zeros(Nnodes, 1);
    u = bim2a_quadtree_solve(msh, A(msh), rhs(msh), u0, dnodes);

    # Save solution to file.
    fclose all;
    filename = sprintf("sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"}, {}, 1);

    # Determine elements to be refined.
    tol = 1e-2;
    refineable_elements = find(!any(msh.children));
    
    to_refine = false(1, Nelements);
    to_refine(refineable_elements) = ...
        parcellfun(4, @(iel) msh2m_to_refine(msh, A, rhs, u, iel, tol),
                   num2cell(refineable_elements));
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    if (!any(to_refine))
        break;
    else
        msh = msh2m_quadtree_refine_recursive(msh, find(to_refine));
    endif
endfor
