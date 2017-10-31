function [epsilon] = diffusion_coefficient(msh, data, omega2)
    Nelems = columns(msh.t);
    
    epsilon = data.eps1 * ones(Nelems, 1);
    epsilon(omega2) = data.eps2;
endfunction
