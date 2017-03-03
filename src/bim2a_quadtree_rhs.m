function rhs = bim2a_quadtree_rhs (msh, f, g)
    if (!isfield(msh, "B") || !isfield(msh, "size"))
        error("bim2a_quadtree_rhs: call bim2c_quadtree_mesh_properties on msh first.");
    end

    refineable_elements = find(!any(msh.children));
    elems = msh.t(1:4, refineable_elements);

    Nnodes = columns(msh.p);
    Nelems = numel(refineable_elements);
    
    II = elems;
    
    VV = zeros(4, Nelems);
    
    hx = msh.size(1, refineable_elements);
    hy = msh.size(2, refineable_elements);
    
    g = g(elems);
    
    VV(1, :) = g(1) .* hx .* hy / 4;
    VV(2, :) = g(2) .* hx .* hy / 4;
    VV(3, :) = g(3) .* hx .* hy / 4;
    VV(4, :) = g(4) .* hx .* hy / 4;
    
    VV .*= repmat(reshape(f(refineable_elements), 1, Nelems), 4, 1);
    
    rhs = sparse(II, 1, VV, Nnodes, 1);
    
    # Change basis.
    rhs = msh.B.' * rhs;
endfunction
