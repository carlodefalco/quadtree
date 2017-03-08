function fpl_vtk_write_field_quadmesh_nedelec (basename, msh, sidedata, endfile)
    refineable_elements = find(!any(msh.children));
        
    for ii = 1 : size(sidedata, 1)
        nodedata_x{ii, 2} = [sidedata{ii, 2} "_x"];
        nodedata_y{ii, 2} = [sidedata{ii, 2} "_y"];
    endfor
    
    for iel = refineable_elements
        nodes = msh.p(:, msh.t(1:4, iel));
        
        new_msh.p(:, end+1 : end+4) = nodes;
        
        [dphi_x, dphi_y] = bim2c_quadtree_eval_dfun_nedelec(msh, iel, nodes(1, :), nodes(2, :));
        
        for ii = 1 : size(sidedata, 1)
            data = sidedata{ii, 1}(msh.ts(:, iel));
            
            nodedata_x{ii, 1}(end+1 : end+4, 1) = dphi_x * data(:);
            nodedata_y{ii, 1}(end+1 : end+4, 1) = dphi_y * data(:);
        endfor
    endfor
    
    new_msh.t = reshape(1:columns(new_msh.p), [4, numel(refineable_elements)]);
    new_msh.children = zeros(4, columns(new_msh.t));
    
    fpl_vtk_write_field_quadmesh(basename, new_msh, [nodedata_x; nodedata_y], {}, endfile);
endfunction
