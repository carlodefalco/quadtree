function [scnodes] = semic_nodes(msh)
    Nnodes = columns(msh.p);

    scnodes = false(1, Nnodes);
    scnodes(msh.p(2, :) <= 0) = true;
endfunction
