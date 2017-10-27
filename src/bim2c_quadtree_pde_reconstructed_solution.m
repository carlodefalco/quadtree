function [u_star_node, u_star_edge, u_star_center] = bim2c_quadtree_pde_reconstructed_solution(msh, u)
    du = bim2c_quadtree_pde_edge_gradient(msh, u);
    [du_x, du_y] = bim2c_quadtree_pde_reconstructed_gradient(msh, du);
    
    Nnodes = columns(msh.p);
    Nelems = columns(msh.t);
    refineable_elements = find(!any(msh.children));
    
    u_star_node = u_star_edge = zeros(numel(u), 1);
    u_star_center = zeros(Nelems, 1);
    
    for iel = refineable_elements
        x1 = msh.p(1, msh.t(1, iel));
        y1 = msh.p(2, msh.t(1, iel));    
        
        x2 = msh.p(1, msh.t(2, iel));
        y2 = msh.p(2, msh.t(3, iel));
        
        hx = x2 - x1;
        hy = y2 - y1;
        
        u_star_node(msh.t(1:4, iel)) = u(msh.t(1:4, iel));
        
        u_star_edge(msh.ts(1, iel)) = mean(u(msh.t( 1:2 , iel))) + hx * (du_x(msh.t(1, iel)) - du_x(msh.t(2, iel))) / 8;
        u_star_edge(msh.ts(2, iel)) = mean(u(msh.t( 2:3 , iel))) + hy * (du_y(msh.t(2, iel)) - du_y(msh.t(3, iel))) / 8;
        u_star_edge(msh.ts(3, iel)) = mean(u(msh.t( 3:4 , iel))) + hx * (du_x(msh.t(4, iel)) - du_x(msh.t(3, iel))) / 8;
        u_star_edge(msh.ts(4, iel)) = mean(u(msh.t([1 4], iel))) + hy * (du_y(msh.t(1, iel)) - du_y(msh.t(4, iel))) / 8;
        
        u_star_center(iel) = mean(u_star_edge(msh.ts(1:4, iel))) + ...
                                  hx * (mean(du_x(msh.t([1 4], iel))) - mean(du_x(msh.t(2:3, iel)))) / 16 + ...
                                  hy * (mean(du_y(msh.t( 1:2 , iel))) - mean(du_y(msh.t(3:4, iel)))) / 16;

    endfor
endfunction
