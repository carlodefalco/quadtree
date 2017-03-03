function [t, varargout] = ...
  msh2m_quadtree_geometrical_properties (msh, nodes=[], weights=[], varargin)
  
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
    x = p(1, t(1:4, iel));
    y = p(2, t(1:4, iel));
    
    [~, idx] = sort(angle(complex(x - mean(x), y - mean(y))));
    
    t(1:4, iel) = t(idx, iel);
  endfor
  
  ## Define nodes and weights.
  if (isempty(nodes) || isempty(weights))
    [nodes, weights] = grule(3);
    
    nodes = 0.5 + 0.5 * nodes;
    weights /= 2;
  endif
  
  if (!issorted(nodes))
    [nodes, idx] = sort(nodes);
    weights = weights(idx);
  endif
  
  [X, Y] = meshgrid(nodes);
  nodes = [X(:) Y(:)].';
  
  [W1, W2] = meshgrid(weights);
  weights = (W1(:) .* W2(:)).';
  
  ## Compute properties.
  for nn = 1 : length(varargin)
    request = varargin{nn};
    
    switch request
      # Weighted Jacobian determinant.
      case "wjacdet"
        varargout{nn} = computearea(p, e, t, nodes, weights, "wjac");
        
      # Quadrature nodes.
      case "nodes"
        varargout{nn} = computearea(p, e, t, nodes, weights, "quad");
        
      # Element sizes [hx; hy].
      case "hx_hy"
        varargout{nn} = computearea(p, e, t, nodes, weights, "size");
        
      # Element areas.
      case "area"
        varargout{nn} = computearea(p, e, t, nodes, weights, "area");
        
      # Change basis.
      case "change_basis"
        varargout{nn} = changebasis(msh);
      
      otherwise
        warning (["msh2m_quadtree_geometrical_properties: ", ...
                  "unexpected value in property string. ", ...
                  "Empty vector will be returned."])
        varargout{nn} = [];
    endswitch

  endfor
endfunction

function [out] = computearea(p, e, t, x, w, request)
  x1 = p(1, t(1, :));
  x2 = p(1, t(2, :));
  
  y1 = p(2, t(1, :));
  y2 = p(2, t(3, :));
  
  hx = x2 - x1;
  hy = y2 - y1;
  
  jacdet = hx .* hy;

  wjacdet = zeros(numel(w), numel(jacdet));
  
  for i = 1 : numel(w)
    wjacdet(i, :) = jacdet .* w(i);
  endfor

  if (request == "wjac")
    out = wjacdet;
  elseif (request == "quad")
    nodes = zeros(2, numel(w), numel(jacdet));
    
    for i = 1 : numel(w)
      nodes(1, i, :) = x1 + hx .* x(1, i);
      nodes(2, i, :) = y1 + hy .* x(2, i);
    endfor
    
    out = nodes;
  elseif (request == "size")
    out = [hx; hy];
  elseif (request == "area")
    out = sum(wjacdet, 1);
  endif
endfunction

function [B] = changebasis(msh)
    Nnodes = columns(msh.p);
    
    hanging = find(any(msh.hanging));
    non_hanging = find(!any(msh.hanging));
    
    # Identity block.
    II = non_hanging;
    JJ = msh.full_to_reduced(non_hanging);
    VV = ones(size(II));
    
    B = sparse(II, JJ, VV, Nnodes, numel(non_hanging));
    
    # Contributions due to hanging nodes.
    II = repmat(hanging, [2, 1])(:);
    JJ = msh.full_to_reduced(msh.hanging(:, hanging))(:);
    VV = 0.5 * ones(size(II));
    
    B += sparse(II, JJ, VV, Nnodes, numel(non_hanging));
endfunction
