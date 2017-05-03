function [insulator] = insulator_elems(msh)
    Nelems = columns(msh.t);

    insulator = false(1, Nelems);
    for iel = 1 : Nelems
        coords = msh.p(:, msh.t(1:4, iel));
        
        if (
            mean (coords(2, :)) >= 0
            && mean (coords(1, :)) >= msh.l + msh.L/4
            && mean (coords(1, :)) <= msh.r - msh.L/4
           )
            insulator(iel) = true;
        endif
    endfor
endfunction
