function [M] = bim2a_quadtree_boundary_mass(msh, sidelist)
    ## Check input
    if nargin > 3
        error("bim2a_quadtree_boundary_mass: wrong number of input parameters.");
    elseif !(isstruct(msh)     && isfield(msh,"p") &&
             isfield (msh,"t") && isfield(msh,"e"))
        error("bim2a_quadtree_boundary_mass: first input is not a valid mesh structure.");
    elseif !( isvector(sidelist) && isnumeric(sidelist) )
        error("bim2a_quadtree_boundary_mass: second input is not a valid numeric vector.");
    endif
    
    nodelist = msh2m_nodes_on_sides(msh, sidelist);
    
    edges = [];
    for is = sidelist
        edges = [edges, msh.e([1:2 5], msh.e(5, :) == is)];
    endfor
    
    non_hanging_edges = [];
    for col = 1 : columns(edges)
        if all(!(all(msh.sides == edges(1:2, col)
                     | flip(msh.sides) == edges(1:2, col))))
            non_hanging_edges(end + 1) = col;
        endif
    endfor
    edges(:, non_hanging_edges) = [];
    
    l = sqrt((msh.p(1, edges(1, :)) - msh.p(1, edges(2, :))).^2 +
             (msh.p(2, edges(1, :)) - msh.p(2, edges(2, :))).^2);
    
    dd = zeros(size(nodelist));
    
    for in = 1 : numel(nodelist)
        dd(in) = (sum(l(edges(1, :) == nodelist(in))) +
                  sum(l(edges(2, :) == nodelist(in)))) / 2;
    endfor
    
    M = sparse(diag(dd));
endfunction
