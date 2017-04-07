function material = material_properties_N2200 (constants, PhiB, sigman)
  material.eps_semic_r = 2.9; # [~]
  material.eps_ins_r   = 3.6; # [~]
  
  material.eps_semic = constants.eps0 * material.eps_semic_r; # [F / m]
  material.eps_ins   = constants.eps0 * material.eps_ins_r; # [F / m]
  
  material.PhiB = -PhiB; # Metal to semiconductor barrier. [eV]
  
  material.N0 = 1e27; # [m^{-3}]
  
  material.sigman = sigman * constants.Kb * 300; # [J]
  material.sigman_kT = material.sigman / (constants.Kb * constants.T0); # [~]
endfunction
