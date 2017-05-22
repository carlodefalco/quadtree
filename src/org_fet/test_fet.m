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
L = 10e-6; # Channel length.
contact = 10e-6;

x_min    = 0;
x_source = x_min + contact;
x_drain  = x_source + L;
x_max    = x_drain + contact;

y_sc      = -35e-9;
y_ins     = 441e-9;

x = unique([linspace(x_min, x_source, 4), ...
            linspace(x_source, x_drain, 3), ...
            linspace(x_drain, x_max, 4)]);
y = unique([linspace(y_sc, 0, 3), ...
            linspace(0, y_ins, 10)]);

msh = msh2m_quadtree(x, y);

# Assign mesh geometrical parameters.
msh.dim.x_min    = x_min;
msh.dim.x_source = x_source;
msh.dim.x_drain  = x_drain;
msh.dim.x_max    = x_max;

msh.dim.y_sc      = y_sc;
msh.dim.y_ins     = y_ins;

# Mark source and drain.
[msh] = mark_contacts(msh);

Nelems_max = 15000;

for i = 1 : 10
    fprintf("i = %d\n", i);
    msh = bim2c_quadtree_mesh_properties(msh, [], []);
    
    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);
    
    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    msh.dim.x_min    = x_min;
    msh.dim.x_source = x_source;
    msh.dim.x_drain  = x_drain;
    msh.dim.x_max    = x_max;
    
    msh.dim.y_sc      = y_sc;
    msh.dim.y_ins     = y_ins;
    
    # Define parameters.
    scnodes = semic_nodes(msh);
    insulator = @(msh) insulator_elems(msh);
    epsilon = @(msh) electrical_permittivity(msh, material, insulator(msh));
    
    # Assemble system.
    A = @(msh) bim2a_quadtree_laplacian(msh, epsilon(msh));

    M = @(msh) bim2a_quadtree_reaction(msh, !insulator(msh), ones(columns(msh.p), 1));
    
    # Source, drain and gate contacts.
    gate   = msh2m_nodes_on_sides(msh, 3);
    source = msh2m_nodes_on_sides(msh, 5);
    drain  = msh2m_nodes_on_sides(msh, 6);
    
    # Initial guess.
    Vg = -20; # [V].
    phi0 = ((y - msh.dim.y_sc) * Vg - (y - msh.dim.y_ins) * material.PhiB) ./ ...
            (msh.dim.y_ins - msh.dim.y_sc);
    
    # Compute solution and error.
    [phi, res, niter, C] = nlpoisson_fet(msh, phi0, A(msh), M(msh), ...
                                         gate, source, drain, ...
                                         material, constants, charge_n);
    
    n = zeros(size(phi));
    n(scnodes) = -charge_n(phi(scnodes)) / constants.q;
    
    E_edge = bim2c_quadtree_pde_edge_gradient(msh, phi);
    [E_x, E_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, E_edge);
    
    # Determine elements to be refined.
    refineable_elements = find(!any(msh.children));
    
    msh_new = mark_regions(msh);
    
    [msh1, omega1, omega1_el_all] = msh2m_quadtree_submesh(msh_new, 1);
    [msh2, omega2, omega2_el_all] = msh2m_quadtree_submesh(msh_new, 2);
    
    [~, omega1_el] = intersect(refineable_elements, omega1_el_all);
    [~, omega2_el] = intersect(refineable_elements, omega2_el_all);
    
    estimator = zeros(1, numel(refineable_elements));
    estimator(omega1_el) = bim2c_quadtree_pde_ZZ_estimator_du(msh1, phi(omega1));
    estimator(omega2_el) = bim2c_quadtree_pde_ZZ_estimator_du(msh2, phi(omega2));
    
    to_refine = false(1, columns(msh.t));
    
    if (mod(i, 2) != 0)
        crit = 1;
        to_refine = bim2c_quadtree_pde_ZZ_to_refine(msh, estimator, mean(estimator), crit);
    else
        to_refine(refineable_elements) = true;
    end
    
    # Save solution to file.
    fclose all;
    basename = "./sol_FET/sol";
    filename = sprintf([basename "_%d"], i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {phi, "phi"; n, "n"}, ...
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
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
