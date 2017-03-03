function A = bim2a_quadtree_reaction (msh, delta, zeta)
    if (!isfield(msh, "B") || !isfield(msh, "size"))
        error("bim2a_quadtree_laplacian: call bim2c_quadtree_mesh_properties on msh first.");
    end

    refineable_elements = find(!any(msh.children));
    elems = msh.t(1:4, refineable_elements);

    Nnodes = columns(msh.p);
    Nelems = numel(refineable_elements);
    
    II = reshape(repmat(elems, 4, 1), 4, 4, Nelems);
    JJ = reshape(repmat(elems(:)', 4, 1), 4, 4, Nelems);
    
    VV = zeros(4, 4, Nelems);
    
    hx = msh.size(1, refineable_elements);
    hy = msh.size(2, refineable_elements);
    
    zeta = zeta(elems);
    
    VV(1, 1, :) = zeta(1, :) .* hx .* hy / 4;
    VV(2, 2, :) = zeta(2, :) .* hx .* hy / 4;
    VV(3, 3, :) = zeta(3, :) .* hx .* hy / 4;
    VV(4, 4, :) = zeta(4, :) .* hx .* hy / 4;
    
    VV .*= repmat(reshape(delta(refineable_elements), 1, 1, Nelems), 4, 4, 1);
    
    A = sparse(II, JJ, VV, Nnodes, Nnodes);
    
    # Change basis.
    A = msh.B.' * A * msh.B;
endfunction
