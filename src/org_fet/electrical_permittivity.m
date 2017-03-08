function [epsilon] = electrical_permittivity(msh, material, insulator)
    Nelems = columns(msh.t);
    
    epsilon = material.eps_semic * ones(Nelems, 1);
    epsilon(insulator) = material.eps_ins;
endfunction
