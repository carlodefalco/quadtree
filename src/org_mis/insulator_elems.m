function [insulator] = insulator_elems(msh)
    Nelems = columns(msh.t);
    
    insulator = false(1, Nelems);
    for iel = 1 : Nelems
        coords = msh.p(:, msh.t(1:4, iel));
        
        if (max(coords(1, :)) <= msh.x_sc_min ||
            min(coords(1, :)) >= msh.x_sc_max ||
            min(coords(2, :)) >= 0)
            
            insulator(iel) = true;
        endif
    endfor
endfunction
