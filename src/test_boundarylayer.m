clear all;
close all;
clc;

n = 10;

# Mesh definition.
x = linspace(0, 1, n);
y = linspace(0, 1, n);

region = 1;
sides = 1:4;

msh = msh2m_quadtree(x, y, region, sides);

tol_max = 1e-3;
Nelems_max = 10000;

for i = 1:10
    fprintf("i = %d\n", i);
    
    # Build global matrix.
    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);

    # Define parameters and exact solution.
    epsilon = 1e-2;

    # Assemble system.
    alpha = @(msh) epsilon * ones(columns(msh.t), 1);
    delta = @(msh) ones(columns(msh.t), 1);
    zeta  = @(msh) 1 + msh.p(1, :).^2 .* msh.p(2, :).^2;
    
    A = @(msh) bim2a_quadtree_laplacian(msh, alpha(msh)) + ...
               bim2a_quadtree_reaction (msh, delta(msh), zeta(msh));
    
    f = @(msh) ones(columns(msh.t), 1);
    g = @(msh) 1 + 2 * msh.p(1, :) .* msh.p(2, :);
    
    rhs = @(msh) bim2a_quadtree_rhs(msh, f(msh), g(msh));

    # Compute solution and error.
    u = zeros(columns(msh.p), 1);
    
    # Boundary conditions.
    d14 = msh2m_nodes_on_sides(msh, [1 4]);
    d2 = msh2m_nodes_on_sides(msh, 2);
    d3 = msh2m_nodes_on_sides(msh, 3);
    
    dnodes = unique([d14 d2 d3]);
    
    u(d14) = 1;
    u(d2) = 1 - msh.p(2, d2).^2;
    u(d3) = 1 - msh.p(1, d3).^2;
    
    u = bim2a_quadtree_solve(msh, A(msh), rhs(msh), u, dnodes);

    # Determine elements to be refined.
    to_refine = false(1, Nelems);
    
    estimator = bim2c_quadtree_pde_ZZ_estimator_du(msh, u);
    tol = mean(estimator);
    
    refineable_elements = find(!any(msh.children));
    to_refine(refineable_elements) = (estimator > tol);

    # Save solution to file.
    fclose all;
    basename = "./soluzione_boundary_layer_1e-2_ZZ_du/sol";
    filename = sprintf([basename "_%d"], i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {u, "u"}, {estimator.', "estimator"}, 1);
    
    n_dofs(i) = sum(!any(msh.hanging));
    n_elems(i) = numel(refineable_elements);
    n_to_refine(i) = sum(to_refine);
    global_estimator(i) = norm(estimator, 2);

    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", "global_estimator");
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    if (n_elems(i) >= Nelems_max || (tol <= tol_max / n_elems(i)))
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
