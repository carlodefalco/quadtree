
function quadmesh (varargin)
  p = inputParser ();
  p.addRequired ("msh", @isstruct);
  p.addSwitch ("show_cell_numbers");
  p.addSwitch ("show_node_numbers");
  p.parse (varargin{:});
  
  for ii = 1 : columns (p.Results.msh.t)
    if (! any (p.Results.msh.children(:, ii)))
      plot (p.Results.msh.p(1, p.Results.msh.t([1:4 1], ii)),
            p.Results.msh.p(2, p.Results.msh.t([1:4 1], ii)), 'x-');
      hold on
    endif
  endfor

  if (p.Results.show_cell_numbers)
    for ii = 1 : columns (p.Results.msh.t)
      coords = ...
      sum (p.Results.msh.p(:, p.Results.msh.t(1:4, ii)).* ...
           [1.2 1.2 0.8 0.8], 2) / 4;
      text (coords(1), coords(2), num2str (ii),
            'backgroundcolor', [1 0 0],
            'fontweight', 'bold',
            'color', [0 1 0]);
      hold on
    endfor
  endif

  if (p.Results.show_node_numbers)
    for ii = 1 : columns (p.Results.msh.p)
      coords = p.Results.msh.p(:, ii);
      text (coords(1), coords(2), num2str (ii),
            'backgroundcolor', [0 1 0],
            'fontweight', 'bold',
            'color', [1 0 0]);
      hold on
    endfor
  endif

  hold off
  
endfunction
