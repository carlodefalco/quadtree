## Solve the non-linear Poisson problem using Newton's algorithm.
function [phiout, resnrm, iter] = ...
         nlpoisson(msh, phi0, A, M, dnodes, charge_n)
  
  non_hanging = find(msh.full_to_reduced);
  
  # Newton's algorithm.
  maxit = 1000;
  tol = 1e-6;
  
  phi    = phi0;
  phiout = phi0;
  
  for iter = 1 : maxit
      phiout = phi(non_hanging);
      
      [rho, drho] = charge_n (phiout);
      
      res = (A * phiout - M * rho);
      
      jac = A - M .* diag(drho);
      
      dphi_bc = zeros(size(phi));

      dphi = -bim2a_quadtree_solve(msh, jac, res, dphi_bc, dnodes);
      
      phi += dphi;
      
      resnrm(iter) = norm(dphi, inf);
      
      if resnrm(iter) < tol
          break;
      endif
  endfor
  
  phiout = phi;
endfunction
