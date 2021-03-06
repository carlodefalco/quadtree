## Solve the non-linear Poisson problem using Newton's algorithm.
function [phiout, resnrm, iter, C] = ...
         nlpoisson_fet(msh, phi0, A, M, gate, source, drain, material, constants, charge_n)
    
    dnodes = unique([gate, source, drain]);
    non_hanging = find(msh.full_to_reduced);
    
    # Newton's algorithm.
    maxit = 1000;
    tol = 1e-6;
    
    phi    = phi0;
    phiout = phi0;
    
    # Boundary conditions.
    dphi_bc = zeros(size(phi));
    
    M_source = material.eps_semic * bim2a_quadtree_boundary_mass(msh, 5);
    source_red = msh.full_to_reduced(source);
    
    M_drain = material.eps_semic * bim2a_quadtree_boundary_mass(msh, 6);
    drain_red = msh.full_to_reduced(drain);
    
    for iter = 1 : maxit
        # Assemble system.
        phiout = phi(non_hanging);
        
        [rho, drho] = charge_n(phiout);
        
        res = (A * phiout - M * rho);
        
        jac = A - M .* sparse(diag(drho));
        
        # Source contact.
        E = -M_source \ res(source_red); # Outward normal electric field.
        
        [PhiBcorr, dPhiBcorr] = barrier_lowering(E, material, constants);
        
        res(source_red) = M_source * (phi(source) - (material.PhiB + constants.Vth * PhiBcorr));
        
        for col = 1 : columns(jac)
            jac(source_red, col) .*= constants.Vth * dPhiBcorr;
        end
        jac(source_red, source_red) += M_source;
        
        # Drain contact.
        E = -M_drain \ res(drain_red); # Outward normal electric field.
        
        [PhiBcorr, dPhiBcorr] = barrier_lowering(E, material, constants);
        
        res(drain_red) = M_drain * (phi(drain) - (material.PhiB + constants.Vth * PhiBcorr));
        
        for col = 1 : columns(jac)
            jac(drain_red, col) .*= constants.Vth * dPhiBcorr;
        end
        jac(drain_red, drain_red) += M_drain;
        
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
        C = 8.1e-7 * delta_Q / (msh.dim.x_gate_max - msh.dim.x_gate_min);
    endif
endfunction
