function [scnodes] = semic_nodes(msh)
    Nnodes = columns(msh.p);

    scnodes = false(Nnodes, 1);
    semic = (msh.p(1, :) >= msh.x_sc_min) & (msh.p(1, :) <= msh.x_sc_max) & (msh.p(2, :) <= 0);
    scnodes(semic) = true;
endfunction
