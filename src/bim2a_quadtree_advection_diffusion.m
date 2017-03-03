function A = bim2a_quadtree_advection_diffusion (msh, alpha, psi)
    if (!isfield(msh, "B") || !isfield(msh, "size"))
        error("bim2a_quadtree_advection_diffusion: call bim2c_quadtree_mesh_properties on msh first.");
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
    
    psi = psi(elems);
    
    psi12 = psi(2, :) - psi(1, :);
    psi23 = psi(3, :) - psi(2, :);
    psi34 = psi(4, :) - psi(3, :);
    psi41 = psi(1, :) - psi(4, :);

    [bp12, bm12] = bimu_bernoulli(psi12);
    [bp23, bm23] = bimu_bernoulli(psi23);
    [bp34, bm34] = bimu_bernoulli(psi34);
    [bp41, bm41] = bimu_bernoulli(psi41);

    bp12 = bp12 .* hy ./ (2 * hx);
    bm12 = bm12 .* hy ./ (2 * hx);
    bp23 = bp23 .* hx ./ (2 * hy);
    bm23 = bm23 .* hx ./ (2 * hy);
    bp34 = bp34 .* hy ./ (2 * hx);
    bm34 = bm34 .* hy ./ (2 * hx);
    bp41 = bp41 .* hx ./ (2 * hy);
    bm41 = bm41 .* hx ./ (2 * hy);
    
    VV(1, 1, :) = bm12 + bp41;
    VV(1, 2, :) = -bp12;
    VV(1, 4, :) = -bm41;
    
    VV(2, 1, :) = -bm12;
    VV(2, 2, :) = bp12 + bm23;
    VV(2, 3, :) = -bp23;
    
    VV(3, 2, :) = -bm23;
    VV(3, 3, :) = bp23 + bm34;
    VV(3, 4, :) = -bp34;
    
    VV(4, 1, :) = -bp41;
    VV(4, 3, :) = -bm34;
    VV(4, 4, :) = bp34 + bm41;
    
    VV .*= repmat(reshape(alpha(refineable_elements), 1, 1, Nelems), 4, 4, 1);
    
    A = sparse(II, JJ, VV, Nnodes, Nnodes);
    
    # Change basis.
    A = msh.B.' * A * msh.B;
endfunction
