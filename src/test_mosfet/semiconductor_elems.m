function [semiconductor] = semiconductor_elems(msh)
  Nelems = columns(msh.t);
  
  semiconductor = false (1, Nelems);
  for iel = 1 : Nelems
    coords = msh.p(:, msh.t(1:4, iel));
    
    if (mean (coords(2, :)) <= 0)
      semiconductor (iel) = true;
    endif
  endfor
endfunction
