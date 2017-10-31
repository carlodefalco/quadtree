function [omega2] = omega2_elems(msh)
    Nelems = columns(msh.t);
    
    omega2 = false(1, Nelems);
    for iel = 1 : Nelems
        x = msh.p(1, msh.t(1:4, iel));
        
        if (min(x) >= 0.5)
            omega2(iel) = true;
        endif
    endfor
endfunction
