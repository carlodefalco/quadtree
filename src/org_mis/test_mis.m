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
L_sc  =  35e-9;
L_ins = 441e-9;
L_ov  =  10e-6;

x = linspace(0, L_ov, 100);
y = linspace(-L_sc, L_ins, 50);

device.msh = msh2m_quadtree(x, y);

for i = 1 : 1
    nnodes = columns(device.msh.p);
    nelems  = columns(device.msh.t);

    device.scnodes = false(nnodes, 1);
    device.scnodes(device.msh.p(2, :) <= 0) = true;

    device.insulator = false(nelems, 1);
    for iel = 1 : nelems
        coords = device.msh.p(:, device.msh.t(1:4, iel));
        
        if (min(coords(2, :)) >= 0)
            device.insulator(iel) = true;
        endif
    endfor

    # Assemble system matrices.
    epsilon = material.eps_semic * ones(nelems, 1);
    epsilon(device.insulator) = material.eps_ins;

    A = bim2a_quadtree_laplacian(device.msh, epsilon);

    M = bim2a_quadtree_reaction(device.msh, !device.insulator, ones(nnodes, 1));

    ## Initial guess
    phi0 = repmat(linspace(material.PhiB, material.PhiB + 10, numel(y)), size(x)).';

    [phi, n, res, niter] = nlpoisson(device, material, constants, ...
                                     phi0, A, M, charge_n);

    fclose all;
    filename = sprintf("sol_%d", i);
    if (exist([filename ".vtu"], "file"))
        delete([filename ".vtu"]);
    endif
    fpl_vtk_write_field_quadmesh(filename, device.msh, {phi, "phi"; n, "n"}, {}, 1);
endfor
