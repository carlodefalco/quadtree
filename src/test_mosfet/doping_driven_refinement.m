function [marker, delta] = doping_driven_refinement (msh, doping, treshold)

  y     = signedlog (doping);
  delta = max (y(msh.t(1:4, :)), [], 1) - min (y(msh.t(1:4, :)), [], 1);
  marker  = delta > treshold;
  
endfunction


function y = signedlog (x)
  y = asinh (x / 2) / log(10);
endfunction
  
