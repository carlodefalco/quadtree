clear all;
close all;
clc;

addpath(canonicalize_file_name("../"));

PhiB = 0.54;
sigman = 2.6;

constants = physical_constants_fun (295);
material  = material_properties_N2200 (constants, PhiB, sigman);

# Quadrature nodes and weights
[gx, gw] = ghrule (101);
quadrature.gx = gx;
quadrature.gw = gw;

# Initialize gaussian charge rules
charge_n = @(phi) gaussian_charge_n (phi, material, constants, quadrature);

# Create mesh.
x_min      = 0.70e-6;
x_sc_max   = 1.15e-6;
x_bulk_max = 1.20e-6;
x_max      = 1.40e-6;

y_sc  = -35e-9;
y_ins = 441e-9;

n = 2.^(0:5);

tol_max = 1e-3;
Nelems_max = 50000;

for i = 1 : numel(n)
    fprintf("i = %d\n", i);
    x = union([linspace(x_min, x_sc_max, 10 * n(i)), ...
           linspace(x_sc_max, x_bulk_max, 2 * n(i))],
           linspace(x_bulk_max, x_max, 5 * n(i)));
    y = union(linspace(y_sc, 0, 2 * n(i)), linspace(0, y_ins, 10 * n(i)));
    
    msh = msh2m_quadtree(x, y);
    
    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);
    
    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    msh.dim.x_min      = x_min;
    msh.dim.x_sc_min   = x_min;
    msh.dim.x_sc_max   = x_sc_max;
    msh.dim.x_bulk_max = x_bulk_max;
    msh.dim.x_max      = x_max;
    msh.dim.y_sc       = y_sc;
    msh.dim.y_ins      = y_ins;
    
    # Define parameters.
    scnodes = semic_nodes(msh);
    insulator = @(msh) insulator_elems(msh);
    epsilon = @(msh) electrical_permittivity(msh, material, insulator(msh));
    
    # Assemble system.
    A = @(msh) bim2a_quadtree_laplacian(msh, epsilon(msh));

    M = @(msh) bim2a_quadtree_reaction(msh, !insulator(msh), ones(columns(msh.p), 1));

    # Initial guess.
    Vg = 10; # [V].
    phi0 = ((y - msh.dim.y_sc) * Vg - (y - msh.dim.y_ins) * material.PhiB) ./ ...
           (msh.dim.y_ins - msh.dim.y_sc);
    
    # Bulk and gate contacts.
    bulk = intersect(msh2m_nodes_on_sides(msh, 1), find(x <= msh.dim.x_bulk_max));
    gate = msh2m_nodes_on_sides(msh, 3);
    dnodes = union(bulk, gate);
    
    # Compute solution and error.
    [phi, res, niter, C] = nlpoisson(msh, phi0, A(msh), M(msh), dnodes, charge_n);
    
    nout = zeros(size(phi));
    nout(scnodes) = -charge_n(phi(scnodes)) / constants.q;

    # Determine elements to be refined.
    refineable_elements = find(!any(msh.children));
    
    msh_new = mark_regions_2D(msh);
    
    [msh1, omega1, omega1_el_all] = msh2m_quadtree_submesh(msh_new, 1);
    [msh2, omega2, omega2_el_all] = msh2m_quadtree_submesh(msh_new, 2);
    
    [~, omega1_el] = intersect(refineable_elements, omega1_el_all);
    [~, omega2_el] = intersect(refineable_elements, omega2_el_all);
    
    estimator = zeros(1, numel(refineable_elements));
    estimator(omega1_el) = bim2c_quadtree_pde_ZZ_estimator_du(msh1, phi(omega1));
    estimator(omega2_el) = bim2c_quadtree_pde_ZZ_estimator_du(msh2, phi(omega2));
    
    crit = 4;
    to_refine = bim2c_quadtree_pde_ZZ_to_refine(msh, estimator, 0.75, crit);
    
    # Save solution to file.
    fclose all;
    basename = "./sol_MIS_2D_uniform/sol";
    filename = sprintf([basename "_%d"], i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {phi, "phi"; nout, "n"}, ...
                                                {estimator.', "estimator"}, 1);
    
    n_dofs(i) = sum(!any(msh.hanging));
    n_elems(i) = numel(refineable_elements);
    n_to_refine(i) = sum(to_refine);
    global_estimator(i) = norm(estimator, 2);
    capacitance(i) = C;
    
    save("-binary", sprintf([basename "_%d_data.mat"], i), "*");
    
    save("-text", [basename "_results.txt"], ...
         "n_dofs", "n_elems", "n_to_refine", "global_estimator", "capacitance");
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    # Do refinement.
    if (n_elems(i) >= Nelems_max)
        break;
    else
        #msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
