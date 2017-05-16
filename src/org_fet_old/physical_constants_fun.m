function constants = physical_constants_fun (T0)
  constants.Kb   = 1.380648813131e-23;
  constants.q    = 1.602176620898e-19;
  constants.eps0 = 8.854187817620e-12;
  
  if (nargin == 1)
    constants.T0 = T0;
  else
    constants.T0 = 300;
  endif
  
  constants.Vth = constants.Kb * constants.T0 / constants.q;  
endfunction
