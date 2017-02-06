function [to_refine] = bim2c_quadtree_pde_ZZ_to_refine(msh, estimator, tol, criterion)
    refineable_elements = find(!any(msh.children));
    
    to_refine = false(1, columns(msh.t));
    
    switch criterion
        case 1
            to_refine(refineable_elements) = (estimator > tol);
        case 2
            if (!isfield(msh, "area"))
                [~, msh.area] = ...
                    msh2m_quadtree_geometrical_properties(msh, [], [], "area");
            endif
            
            to_refine(refineable_elements) = (estimator ./ msh.area(refineable_elements) > tol);
        case 3
            [estimator_sorted, idx_sorted] = sort(estimator, "descend");
            estimator_sum = cumsum(estimator_sorted);
            
            idx = find(estimator_sum >= tol * estimator_sum(end), 1);
            
            to_refine(refineable_elements(idx_sorted(1:idx))) = true;
        case 4
            [estimator_sorted, idx_sorted] = sort(estimator, "descend");
            estimator_sum = cumsum(estimator_sorted);
            
            idx = find(estimator_sum >= tol * estimator_sum(end), 1);
            
            threshold = estimator_sorted(idx) * 0.95;

            to_refine(refineable_elements) = (estimator >= threshold);
        otherwise
            error("bim2c_quadtree_pde_ZZ_to_refine: wrong criterion selected.");
    endswitch
endfunction
