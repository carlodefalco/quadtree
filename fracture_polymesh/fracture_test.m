clear all;
close all;
clc;

x = linspace(0, 1, 20);
y = linspace(0, 1, 20);

qmsh = msh2m_quadtree(x, y);

msh.p = qmsh.p;
msh.t = mat2cell(qmsh.t, 5, ones(columns(qmsh.t), 1));

clear qmsh;

# List of segments, each one represented
# as a pair of nodes {A, B}.
fractures = {{[0 0.48], [1, 0.15]},
             {[0 0.5 ], [1, 0.8 ]}};

figure;
plot_polymesh(msh, "show_cell_numbers", "show_node_numbers");
title("Original mesh");

msh = msh2m_polymesh(msh, fractures);

figure;
plot_polymesh(msh, "show_cell_numbers", "show_node_numbers");
title("Polygonal mesh");
