## Gaussian charge rule for electrons.
function [rhon, drhon_dV] = gaussian_charge_n (V, material, constants, quadrature)
  rhon = -constants.q * n_approx (V, material, constants, quadrature);
  
  if nargout > 1
    drhon_dV = -constants.q * dn_dV_approx (V, material, constants, quadrature);
  endif
endfunction

function n = n_approx (V, material, constants, quadrature)
  persistent N0     = material.N0;
  persistent sigman = material.sigman;
  persistent kT     = constants.Kb * constants.T0;
  
  persistent gx = quadrature.gx;
  persistent gw = quadrature.gw;
  
  n = zeros (size (V));
  
  for i = 1:numel (gx)
    coeff = (sqrt(2) * sigman * gx(i) - constants.q * V) / kT;
    
    n += N0 / sqrt(pi) * gw(i) ./ (1 + exp (coeff));
  endfor
endfunction

function dn_dV = dn_dV_approx (V, material, constants, quadrature)
  persistent N0     = material.N0;
  persistent sigman = material.sigman;
  persistent kT     = constants.Kb * constants.T0;
  
  persistent gx = quadrature.gx;
  persistent gw = quadrature.gw;
  
  dn_dV = zeros (size (V));
  
  for i = 1:numel (gx)
    coeff = (sqrt(2) * sigman * gx(i) - constants.q * V) / kT;
    
    dn_dV += -constants.q * N0 / sigman * sqrt(2 / pi) * ...
              gw(i) * gx(i) ./ (1 + exp (coeff));
  endfor
endfunction
