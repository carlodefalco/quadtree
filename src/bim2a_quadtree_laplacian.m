function A = bim2a_quadtree_laplacian (msh, D)
    A = bim2a_quadtree_advection_diffusion(msh, D, zeros(columns(msh.p), 1));
endfunction
