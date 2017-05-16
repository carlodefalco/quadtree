## Solve the non-linear Poisson problem using Newton's algorithm.
function [phiout, resnrm, iter, C] = ...
         nlpoisson_fet(msh, phi0, A, M, ...
                       gate, source1, source2, drain1, drain2, ...
                       material, constants, charge_n)
    
    dnodes = unique([gate, source1, source2, drain1, drain2]);
    non_hanging = find(msh.full_to_reduced);
    
    # Newton's algorithm.
    maxit = 1000;
    tol = 1e-6;
    
    phi    = phi0;
    phiout = phi0;
    
    # Boundary conditions.
    dphi_bc = zeros(size(phi));
    
    M_source1 = material.eps_semic * bim2a_quadtree_boundary_mass(msh, 5);
    corner1 = find(ismember(source1, source2));
    M_source1(corner1, :) = [];
    M_source1(:, corner1) = [];
    
    M_source2 = material.eps_semic * bim2a_quadtree_boundary_mass(msh, 6);
    
    M_drain1 = material.eps_semic * bim2a_quadtree_boundary_mass(msh, 7);
    
    M_drain2 = material.eps_semic * bim2a_quadtree_boundary_mass(msh, 8);
    corner2 = find(ismember(drain2, drain1));
    M_drain2(corner2, :) = [];
    M_drain2(:, corner2) = [];
    
    source1 = setdiff(source1, source2);
    drain2  = setdiff(drain2 , drain1 );
    
    source1_red = msh.full_to_reduced(source1);
    source2_red = msh.full_to_reduced(source2);
    drain1_red  = msh.full_to_reduced(drain1 );
    drain2_red  = msh.full_to_reduced(drain2 );
    
    for iter = 1 : maxit
        # Assemble system.
        phiout = phi(non_hanging);
        
        [rho, drho] = charge_n(phiout);
        
        res = (A * phiout - M * rho);
        
        jac = A - M .* sparse(diag(drho));
        
        # Source contact top.
        E_source1 = -M_source1 \ res(source1_red); # Outward normal electric field.
        
        [PhiBcorr_source1, dPhiBcorr_source1] = barrier_lowering(E_source1, material, constants);
        
        res(source1_red) = M_source1 * (phi(source1) - (material.PhiB + constants.Vth * PhiBcorr_source1));
        
        for col = 1 : columns(jac)
            jac(source1_red, col) .*= constants.Vth * dPhiBcorr_source1;
        end
        jac(source1_red, source1_red) += M_source1;
        
        # Source contact right.
        E_source2 = -M_source2 \ res(source2_red); # Outward normal electric field.
        
        [PhiBcorr_source2, dPhiBcorr_source2] = barrier_lowering(E_source2, material, constants);
        
        res(source2_red) = M_source2 * (phi(source2) - (material.PhiB + constants.Vth * PhiBcorr_source2));
        
        for col = 1 : columns(jac)
            jac(source2_red, col) .*= constants.Vth * dPhiBcorr_source2;
        end
        jac(source2_red, source2_red) += M_source2;
        
        
        # Drain contact left.
        E_drain1 = -M_drain1 \ res(drain1_red); # Outward normal electric field.
        
        [PhiBcorr_drain1, dPhiBcorr_drain1] = barrier_lowering(E_drain1, material, constants);
        
        res(drain1_red) = M_drain1 * (phi(drain1) - (material.PhiB + constants.Vth * PhiBcorr_drain1));
        
        for col = 1 : columns(jac)
            jac(drain1_red, col) .*= constants.Vth * dPhiBcorr_drain1;
        end
        jac(drain1_red, drain1_red) += M_drain1;
        
        # Drain contact top.
        E_drain2 = -M_drain2 \ res(drain2_red); # Outward normal electric field.
        
        [PhiBcorr_drain2, dPhiBcorr_drain2] = barrier_lowering(E_drain2, material, constants);
        
        res(drain2_red) = M_drain2 * (phi(drain2) - (material.PhiB + constants.Vth * PhiBcorr_drain2));
        
        for col = 1 : columns(jac)
            jac(drain2_red, col) .*= constants.Vth * dPhiBcorr_drain2;
        end
        jac(drain2_red, drain2_red) += M_drain2;
        
        
        # Compute solution.
        dphi = -bim2a_quadtree_solve(msh, jac, res, dphi_bc, gate);
        
        clamp = 0.01;
        dphi(dphi >  clamp) =  clamp;
        dphi(dphi < -clamp) = -clamp;
        
        phi += dphi;
        
        resnrm(iter) = norm(dphi, inf);
        fprintf("%d, %g\n", iter, resnrm(iter));
        
        if resnrm(iter) < tol
            break;
        endif
    endfor
    
    phiout = phi;
    
    if (nargout == 4)
        # Assemble system.
        [~, drho]  = charge_n(phi(non_hanging));
        
        mat = A - M .* sparse(diag(drho));
        rhs = zeros(size(phiout));
        
        # Boundary conditions.
        y = msh.p(2, :).';
        delta_phi0 = (y - msh.dim.y_contact) ./ (msh.dim.y_ins - msh.dim.y_contact);
        delta_phi0(y < msh.dim.y_contact) = 0;
        
        # Compute solution.
        delta_phi = bim2a_quadtree_solve(msh, mat, rhs, delta_phi0, dnodes);
        
        # Compute capacitance.
        delta_Q = sum(mat(msh.full_to_reduced(gate), :) * delta_phi(non_hanging));
        C = 8.1e-7 * delta_Q / (msh.dim.x_max - msh.dim.x_min);
    endif
endfunction
