function [insulator] = insulator_elems(msh)
    Nelems = columns(msh.t);

    l = min (msh.p (1, :));
    r = max (msh.p (1, :));
    L =  r - l;
    
    insulator = false(1, Nelems);
    for iel = 1 : Nelems
        coords = msh.p(:, msh.t(1:4, iel));
        
        if (
            mean (coords(2, :)) >= 0
            && mean (coords(1, :)) >= l + L/4
            && mean (coords(1, :)) <= l + L - L/4
           )
            insulator(iel) = true;
        endif
    endfor
endfunction
