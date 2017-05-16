## Physical models for Schottky barrier lowering.
function [PhiBcorr, dPhiBcorr] = barrier_lowering(E, material, constants)
    persistent schottky = sqrt(constants.q / (4 * pi * material.eps_semic)) / constants.Vth;
    
    PhiBcorr_inj = @(f) sqrt(f);
    PhiBcorr_ext = @(f) f / 4;
    
    f = schottky^2 * E;
    
    PhiBcorr = dPhiBcorr = zeros(size(f));
    
    if (any(f > 0))
        PhiBcorr (f > 0) = PhiBcorr_inj(f(f > 0));
        dPhiBcorr(f > 0) = schottky^2 / (2 * sqrt(f(f > 0)));
    end
    
    if (any(f < 0))
        PhiBcorr (f < 0) = PhiBcorr_ext(f(f < 0));
        dPhiBcorr(f < 0) = schottky^2 / 4;
    end
endfunction
