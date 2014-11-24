

[x,y] = meshgrid(linspace(0,1,100));
z = peaks(100);
clev = [-5 -3 -2:.5:2 3 5];

h = plotgrid('setup', cell(1,2), [],[],'sp', 0.05, 'mar', 0.05, 'mt', 0.25, 'mb', 0.25);

axes(h.ax(1));
hc = contourfcmap(x,y,z,clev,jet(12), ...
     [.8 .8 .8], [.2 .2 .2], 'eastoutside')
           
axes(h.ax(2));
contourf(x,y,z,clev);
cb = colorbar('eastoutside');
colormap(h.ax(2), jet);

set([h.ax hc.cbax cb], 'fontsize', 8);

set(h.fig, 'color', 'w');

title(h.ax(1), 'contourfcmap');
title(h.ax(2), 'contourf with same contour levels');

export_fig('contourfcmap_example', h.fig, '-png');
