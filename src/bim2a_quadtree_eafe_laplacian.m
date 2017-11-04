function A = bim2a_quadtree_eafe_laplacian (msh, D)
    if (!isfield(msh, "B") || !isfield(msh, "size"))
        error("bim2a_quadtree_laplacian: call bim2c_quadtree_mesh_properties on msh first.");
    end

    refineable_elements = find(!any(msh.children));
    elems = msh.t(1:4, refineable_elements);

    Nnodes = columns(msh.p);
    Nelems = numel(refineable_elements);
    
    II = reshape(repmat(elems, 4, 1), 4, 4, Nelems);
    JJ = permute (II, [2, 1, 3]);    
    
    VV = zeros(4, 4, Nelems);
    
    hx = msh.size(1, refineable_elements);
    hy = msh.size(2, refineable_elements);

    VV(1, 2, :) = VV(2, 1, :) = -hm (D(elems(1,:)), D(elems(2,:))).' .* hy ./ (2 * hx);
    VV(1, 4, :) = VV(4, 1, :) = -hm (D(elems(1,:)), D(elems(4,:))).' .* hx ./ (2 * hy);
    VV(2, 3, :) = VV(3, 2, :) = -hm (D(elems(2,:)), D(elems(3,:))).' .* hx ./ (2 * hy);
    VV(3, 4, :) = VV(4, 3, :) = -hm (D(elems(3,:)), D(elems(4,:))).' .* hy ./ (2 * hx);

    VV(1, 1, :) = -VV(1, 2, :) -VV(1, 4, :);
    VV(2, 2, :) = -VV(2, 1, :) -VV(2, 3, :);
    VV(3, 3, :) = -VV(3, 2, :) -VV(3, 4, :);
    VV(4, 4, :) = -VV(4, 1, :) -VV(4, 3, :);
    
    A = sparse(II, JJ, VV, Nnodes, Nnodes);
    
    # Change basis.
    A = msh.B.' * A * msh.B;
endfunction

function m = hm (a, b)
  m = 2 ./ (1./a + 1./b);
endfunction
