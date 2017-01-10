function [scnodes] = semic_nodes(msh)
    Nnodes = columns(msh.p);

    scnodes = false(1, Nnodes);
    semic = (msh.p(1, :) >= msh.x_sc_min) & (msh.p(1, :) <= msh.x_sc_max) & (msh.p(2, :) <= 0);
    scnodes(semic) = true;
endfunction
