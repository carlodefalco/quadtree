function [omega1] = omega1_nodes(msh)
    Nnodes = columns(msh.p);

    omega1 = find(msh.p(1, :) <= 0.5);
endfunction
