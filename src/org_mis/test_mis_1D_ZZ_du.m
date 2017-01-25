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
x_min    =     0;
x_sc_min =     0;
x_sc_max = 35e-9;
x_max    = 35e-9;

y_sc  = -35e-9;
y_ins = 441e-9;

x = linspace(x_sc_min, x_sc_max, 3);
y = union(linspace(y_sc, 0, 10), linspace(0, y_ins, 10));

msh = msh2m_quadtree(x, y);

for i = 1 : 10
    fprintf("i = %d\n", i);
    
    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);
    
    x = msh.p(1, :).';
    y = msh.p(2, :).';
    
    msh.dim.x_min    = x_min;
    msh.dim.x_sc_min = x_sc_min;
    msh.dim.x_sc_max = x_sc_max;
    msh.dim.x_max    = x_max;
    msh.dim.y_sc     = y_sc;
    msh.dim.y_ins    = y_ins;
    
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
    bulk = msh2m_nodes_on_sides(msh, 1);
    gate = msh2m_nodes_on_sides(msh, 3);
    dnodes = union(bulk, gate);
    
    # Compute solution and error.
    [phi, res, niter, C] = nlpoisson(msh, phi0, A(msh), M(msh), dnodes, charge_n);
    
    n = zeros(size(phi));
    n(scnodes) = -charge_n(phi(scnodes)) / constants.q;
    
    # Save solution to file.
    fclose all;
    filename = sprintf("./sol_MIS_1D_ZZ_du/sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {phi, "phi"; n, "n"}, {}, 1);
    save("-text", [filename "_capacitance.txt"], "C");

    # Determine elements to be refined.
    tol = 1e-2;
    refineable_elements = find(!any(msh.children));
    
    estimator = bim2c_quadtree_pde_ZZ_estimator_du(msh, phi);
    
    to_refine = false(1, Nelems);
    to_refine(refineable_elements) = (estimator > tol / (2^msh.level(refineable_elements)));
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    # Do refinement.
    if (!any(to_refine))
        break;
    else
        msh = msh2m_quadtree_refine_recursive(msh, find(to_refine));
    endif
endfor
