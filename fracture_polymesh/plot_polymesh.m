function plot_polymesh(varargin)
    p = inputParser();
    p.addRequired("msh", @isstruct);
    p.addSwitch("show_cell_numbers");
    p.addSwitch("show_node_numbers");
    p.parse(varargin{:});
    
    msh = p.Results.msh;
    
    for ii = 1 : numel(msh.t)
        plot (msh.p(1, msh.t{ii}([1:end-1 1])),
              msh.p(2, msh.t{ii}([1:end-1 1])), 'x-');
        hold on;
    endfor
    
    if (p.Results.show_cell_numbers)
        for ii = 1 : numel(msh.t)
            coords = mean(msh.p(:, msh.t{ii}(1:end-1)), 2);
            
            text(coords(1), coords(2), num2str (ii),
                 'BackgroundColor', [1 0 0],
                 'FontWeight', 'Bold',
                 'Color', [0 1 0]);
            hold on;
        endfor
    endif
    
    if (p.Results.show_node_numbers)
        for ii = 1 : columns(msh.p)
            coords = msh.p(:, ii);
            text(coords(1), coords(2), num2str (ii),
                 'BackgroundColor', [0 1 0],
                 'FontWeight', 'Bold',
                 'Color', [1 0 0]);
            hold on;
        endfor
    endif
    
    hold off;
endfunction
