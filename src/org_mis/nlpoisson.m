## Solve the non-linear Poisson problem using Newton's algorithm.
function [phiout, nout, resnrm, iter] = ...
         nlpoisson (device, material, constants, phi0, A, M, charge_n)
  
  # Newton's algorithm.
  maxit = 1000;
  tol = 1e-6;
  
  phi    = phi0;
  phiout = phi0;
  
  for iter = 1 : maxit
      phiout = phi;
      
      [rho, drho] = charge_n (phiout);
      
      res = (A * phiout - M * rho);
      
      jac = A - M .* diag (drho);
      
      dphi_bc = zeros(size(phiout));
      
      # Bulk and gate contacts.
      dnodes = msh2m_nodes_on_sides(device.msh, [1 3]);
      dphi_bc(dnodes) = phi(dnodes) - phi0(dnodes);

      dphi = -bim2a_quadtree_solve(device.msh, jac, res, dphi_bc, dnodes);
      
      phi += dphi;
      
      resnrm(iter) = norm(dphi, inf);
      
      if resnrm(iter) < tol
        break;
      endif
  endfor
  
  phiout = phi;
  
  # Post-processing.
  nout = zeros(size(phiout));
  nout(device.scnodes) = -charge_n(phiout(device.scnodes)) / constants.q;
endfunction
