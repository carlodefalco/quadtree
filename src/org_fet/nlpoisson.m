## Solve the non-linear Poisson problem using Newton's algorithm.
function [phiout, resnrm, iter, C] = ...
         nlpoisson(msh, phi0, A, M, gate, source, drain, charge_n)
    
    dnodes = unique([gate, source, drain]);
    non_hanging = find(msh.full_to_reduced);
    
    # Newton's algorithm.
    maxit = 1000;
    tol = 1e-6;
    
    phi    = phi0;
    phiout = phi0;
    
    # Boundary conditions.
    dphi_bc = zeros(size(phi));
    
    for iter = 1 : maxit
        # Assemble system.
        phiout = phi(non_hanging);
        
        [rho, drho] = charge_n(phiout);
        
        res = (A * phiout - M * rho);
        
        jac = A - M .* sparse(diag(drho));
        
        # Compute solution.
        dphi = -bim2a_quadtree_solve(msh, jac, res, dphi_bc, dnodes);
        
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
        delta_phi0 = (y - msh.dim.y_contact) ./ (msh.dim.y_ins - msh.dim.y_sc);
        delta_phi0(y < msh.dim.y_contact) = 0;
        
        # Compute solution.
        delta_phi = bim2a_quadtree_solve(msh, mat, rhs, delta_phi0, dnodes);
        
        # Compute capacitance.
        delta_Q = sum(mat(msh.full_to_reduced(gate), :) * delta_phi(non_hanging));
        C = delta_Q / (msh.dim.x_max - msh.dim.x_min); # Per unit length.
    endif
endfunction
