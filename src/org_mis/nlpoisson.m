## Solve the non-linear Poisson problem using Newton's algorithm.
function [phiout, resnrm, iter, C] = ...
         nlpoisson(msh, phi0, A, M, bulk, gate, constants, material, charge_n)
    
    non_hanging = find(msh.full_to_reduced);
    
    # Newton's algorithm.
    maxit = 1000;
    tol = 1e-6;
    
    phi    = phi0;
    phiout = phi0;
    
    # Boundary conditions.
    dphi_bc = zeros(size(phi));
    
    M_bulk = material.eps_semic * bim2a_quadtree_boundary_mass(msh, 1);
    
    bulk_red = msh.full_to_reduced(bulk);
    
    for iter = 1 : maxit
        # Assemble system.
        phiout = phi(non_hanging);
        
        [rho, drho] = charge_n(phiout);
        
        res = (A * phiout - M * rho);
        
        jac = A - M .* sparse(diag(drho));
        
        # Bulk contact.
        E = -M_bulk \ res(bulk_red); # Outward normal electric field.
        
        [PhiBcorr, dPhiBcorr] = barrier_lowering(E, material, constants);
        
        res(bulk_red) = M_bulk * (phi(bulk) - (material.PhiB + constants.Vth * PhiBcorr));
        
        for col = 1 : columns(jac)
            jac(bulk_red, col) .*= constants.Vth * dPhiBcorr;
        end
        jac(bulk_red, bulk_red) += M_bulk;
        
        # Compute solution.
        dphi = -bim2a_quadtree_solve(msh, jac, res, dphi_bc, gate);
        
        phi += dphi;
        
        resnrm(iter) = norm(dphi, inf);
        
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
        phi0 = (y - msh.dim.y_sc) ./ (msh.dim.y_ins - msh.dim.y_sc);
        
        # Compute solution.
        delta_phi = bim2a_quadtree_solve(msh, mat, rhs, phi0, union(bulk, gate));
        
        # Compute capacitance.
        x = msh.p(1, :).';
        gate1 = msh.full_to_reduced(
                    intersect(msh2m_nodes_on_sides(msh, 3), find(x <= msh.dim.x_sc_max))
                );
        gate2 = msh.full_to_reduced(
                    intersect(msh2m_nodes_on_sides(msh, 3), find(x >= msh.dim.x_sc_max & x <= msh.dim.x_bulk_max))
                );
        gate3 = msh.full_to_reduced(
                    intersect(msh2m_nodes_on_sides(msh, 3), find(x >= msh.dim.x_bulk_max))
                );
        
        # Integrate along x-direction.
        delta_Q1 = 2 * sum(mat(gate1, :) * delta_phi(non_hanging));
        delta_Q2 = 2 * sum(mat(gate2, :) * delta_phi(non_hanging));
        delta_Q3 = 2 * sum(mat(gate3, :) * delta_phi(non_hanging));
        
        # Integrate along z-direction.
        C1 = delta_Q1 * 2 * (msh.dim.x_sc_max - msh.dim.x_min);
        C2 = delta_Q2 * 2 * (msh.dim.x_bulk_max - msh.dim.x_min);
        C3 = delta_Q3 * 2 * (msh.dim.x_max - msh.dim.x_min);
        
        C = C1 + C2 + C3;
    endif
endfunction
