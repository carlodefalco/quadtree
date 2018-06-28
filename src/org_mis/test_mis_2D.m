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
x_min      = 0.70e-3;
x_sc_max   = 1.15e-3;
x_bulk_max = 1.20e-3;
x_max      = 1.40e-3;

y_sc  = -35e-9;
y_ins = 441e-9;

x = union([linspace(x_min, x_sc_max, 10), ...
           linspace(x_sc_max, x_bulk_max, 2)],
           linspace(x_bulk_max, x_max, 5));
y = union(linspace(y_sc, 0, 2), linspace(0, y_ins, 10));

msh = msh2m_quadtree(x, y);

# Mark edge labels.
for iedge = 1 : columns(msh.e)
    coords = msh.p(:, msh.e(1:2, iedge));
    
    x_min = min(coords(1, :));
    
    if (msh.e(5, iedge) == 1 &&
        x_min >= x_bulk_max)
        
        msh.e(5, iedge) = 5;
    endif
endfor
msh.onboundary(msh2m_nodes_on_sides(msh, 5)) = 5;

for i = 1 : 15
    fprintf("i = %d\n", i);
    msh = bim2c_quadtree_mesh_properties(msh, [], []);
    
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
    Vg = 20; # [V].
    phi0 = ((y - msh.dim.y_sc) * Vg - (y - msh.dim.y_ins) * material.PhiB) ./ ...
           (msh.dim.y_ins - msh.dim.y_sc);
    
    # Bulk and gate contacts.
    bulk = msh2m_nodes_on_sides(msh, 1);
    gate = msh2m_nodes_on_sides(msh, 3);
    
    # Compute solution and error.
    [phi, res, niter, C] = nlpoisson(msh, phi0, A(msh), M(msh), bulk, gate, constants, material, charge_n);
    
    n = zeros(size(phi));
    n(scnodes) = -charge_n(phi(scnodes)) / constants.q;
    
    # Save solution to file.
    fclose all;
    filename = sprintf("./sol_MIS_2D/sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, msh, {phi, "phi"; n, "n"}, {}, 1);
    save("-text", [filename "_capacitance.txt"], "C");

    # Determine elements to be refined.
    tol = 1e-2;
    refineable_elements = find(!any(msh.children));
    
    to_refine = false(1, Nelems);
    to_refine(refineable_elements) = ...
        parcellfun(4, @(iel) msh2m_to_refine_mis(msh, material, constants,
                                                 A, M, phi, charge_n, iel, tol),
                   num2cell(refineable_elements));
    
    fprintf("Elements to refine = %d / %d\n\n", sum(to_refine), numel(refineable_elements));
    
    # Do refinement.
    if (!any(to_refine))
        break;
    else
        msh = msh2m_quadtree_refine(msh, find(to_refine));
    endif
endfor
