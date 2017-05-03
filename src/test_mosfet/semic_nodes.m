function [scnodes] = semic_nodes (msh)
  scnodes_idx = unique (msh.t(1:4, msh.t(5, :) == 1)(:));
  scnodes = false (1, columns (msh.p));
  scnodes(scnodes_idx) = true;
endfunction
