function [t, varargout] = msh2m_quadtree_geometrical_properties (msh, varargin)
  ## Check inputs.
  if (nargin < 2)
    error (["msh2m_quadtree_geometrical_properties: ", ...
            "wrong number of input parameters."]);
  elseif (! (isstruct (msh) && isfield (msh, "p") 
             && isfield (msh, "t") && isfield (msh, "e")))
    error (["msh2m_quadtree_geometrical_properties: ", ...
            "first input is not a valid mesh structure."]);
  elseif (! iscellstr (varargin))
    error (["msh2m_quadtree_geometrical_properties: ", ... 
            "only string value admitted for properties."]);
  endif
  
  p = msh.p;
  e = msh.e;
  t = msh.t;
  
  ## Sort vertices counterclockwise (starting from west).
  for iel = 1 : columns(t)
    nodes = t(1:4, iel);
    
    x = p(1, nodes);
    y = p(2, nodes);
    
    [~, idx] = sort(angle(complex(x - mean(x), y - mean(y))));
    
    t(1:4, iel) = t(idx, iel);
  endfor
  
  ## Compute properties.
  
  for nn = 1 : length(varargin)
    
    request = varargin{nn};
    
    switch request
      # Weighted Jacobian determinant.
      case "wjacdet"
        varargout{nn} = computearea(p, e, t, "wjac");
      
      # Element areas.
      case "area"
        varargout{nn} = computearea(p, e, t, "area");
      
      # Gradient of hat functions.
      case "shg"
        varargout{nn} = shapegrad(p, t);
        
      otherwise
        warning (["msh2m_quadtree_geometrical_properties: ", ...
                  "unexpected value in property string. ", ...
                  "Empty vector will be returned."])
        varargout{nn} = [];
    endswitch

  endfor

endfunction

function [out] = computearea(p, e, t, request)
  weight = 1/4 * ones(1, 4);
  
  hx = p(1, t(1, :)) - p(1, t(2, :));
  hy = p(2, t(1, :)) - p(2, t(3, :));
  
  jacdet = hx .* hy;

  for i = 1 : numel(weight)
    wjacdet(i, :) = jacdet .* weight(i);
  endfor

  if (request == "wjac")
    out = wjacdet;
  elseif (request == "area")
    out = sum(wjacdet, 1)';
  endif
endfunction

function [shg] = shapegrad(p, t)
  x1 = p(1, t(1, :));
  y1 = p(2, t(1, :));
  
  x2 = p(1, t(2, :));
  y2 = p(2, t(2, :));
  
  x3 = p(1, t(3, :));
  y3 = p(2, t(3, :));
  
  x4 = p(1, t(4, :));
  y4 = p(2, t(4, :));

  shg(1, 1, :) =  (y1 - y2);
  shg(2, 1, :) =  (x1 - x2);
  
  shg(1, 2, :) = -(y2 - y2);
  shg(2, 2, :) =  (x2 - x1);
  
  shg(1, 3, :) =  (y3 - y1);
  shg(2, 3, :) = -(x3 - x1);
  
  shg(1, 4, :) =  (y4 - y1);
  shg(2, 4, :) = -(x4 - x2);
endfunction
