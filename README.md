
# contourfcmap.m: filled contour plot with precise colormap


Author: Kelly Kearney


This repository includes the code for the `contourfcmap.m` Matlab function, along with all dependent functions required to run it.


This function creates a shaded contour map, similar to that created by the contourf function.  However, the relationship between a contourf plot and its colormap (i.e. exactly which color corresponds to each contour interval), can often be confusing and inconsistent.  This function instead allows the user to specify exactly which colors to use in each interval, and also to choose colors for regions that exceed the contour line limits.



## Contents


- Getting started        
- Syntax        
- Examples        
- Contributions

## Getting started


**Prerequisites**


This function requires Matlab R14 or later.


**Downloading and installation**


This code can be downloaded from [Github](https://github.com/kakearney/contourfcmap-pkg/) or the [MatlabCentral File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/29638).  The File Exchange entry is updated daily from the GitHub repository.


**Matlab Search Path**


The following folders need to be added to your Matlab Search path (via `addpath`, `pathtool`, etc.):



```matlab
contourfcmap-pkg/FEX-function_handle
contourfcmap-pkg/contourcs
contourfcmap-pkg/contourfcmap
contourfcmap-pkg/fillnan
contourfcmap-pkg/multiplepolyint
contourfcmap-pkg/parsepv
```



## Syntax



```
h = contourfcmap(x,y,z,clev,cmap,lo,hi,cbarloc)
h = contourfcmap(x,y,z,clev,cmap)
h = contourfcmap(x,y,z,clev,cmap,lo)
h = contourfcmap(x,y,z,clev,cmap,lo,hi)
h = contourfcmap(x,y,z,clev,cmap, param1, val1, ...)
```


Input variables



- `x`:          x-coordinates of grid, following size restrictions of surf               command
- `y`:          y-coordinates of grid, following size restrictions of surf               command
- `z`:          array of data to be contoured, following size restritions               of surf command
- `clev`:       vector of length n, contour levels, must be monotonically               increasing
- `cmap`:       n-1 x 3 colormap array, specifying colors used to shade               each contour interval

Optional input variables (passed as parameter/value pairs)



- `lo`:         1 x 3 colormap array, specifying color used for all data               falling below the first contour level.  If not included or               empty, default will be to white.


- `hi`:         1 x 3 colormap array, specifying color used for all data               falling above the last contour level.  If not included or               empty, default will be to white.


- `cbarloc`:    string specifying colorbar location (see colorbar),  or a               1 x 4 position vector for colorbar.  If not included, no               colorbar will be created.  Note that the colorbar created               is not a colorbar object, but just an axis with plotted               patches; it is for labeling only and is not linked to the               "peer axis" in any way. Default is no colorbar.


- `evencb`:     logical scalar.  If true, intervals on the colorbar will               be evenly spaced, regardless of value.  If false, ticks               will correspond to clev values. If not included or empty,               default is false.


- `method`:     string specifying calculation method.               'recolor' creates a contourf object, and                               then change the color properties of the                               underlying components. Note: In 2014b,                               recolor will not persist when saving to                               file via anything that uses print.m                               (including print, saveas, export_fig, etc).               'calccontour' creates new patch and line objects from                               scratch.  This method requires the Mapping                               Toolbox.  It can be beneficial if you want                               more consistency between output regardless                               of which version of Matlab is being run.                               It also properly fills regions falling                               below the lowest contour level and above                               the highest contour level, which isn't                               always the case with contourf-generated                               contour objects.


- `flag`:       logical flag indicating whether to use the fast version of               multiplepolyint.m (only applicable if method is               'calccontour').  In plots with a large number of contour               lines, this significantly speeds up the calculation.               However, it works by directly calling a private mex               function in the Mapping Toolbox, so it might be a bit less               stable than the slower version, and may break in future               releases. Default is true.

Output variables:



- `h.c`: contour matrix for filled contour plot ('recolor'                       method only)
- `h.h`:      handle to contourgroup or contour object ('recolor'                       method only)
- `h.p`:      handles to patch objects (the filled contour                       regions) ('calccontour' method only)
- `h.l`:      handles to line objects (the contour lines)                       ('calccontour' method only)
- `h.cbax`:   handle to axis of colorbar


## Examples


First we'll plot a contourf plot using the standard Matlab functions. Without labeling the contour lines, it's difficult to tell which values in the colorbar correspond to the color intervals in the contuorf plot.



```matlab
count = 0;

% Data

[x,y] = meshgrid(linspace(0,1,100));
z = peaks(100);
clev = [-5 -3 -2:.5:2 3 5];

% Set up axes

h.fig = figure('color', 'w');
h.ax(1) = axes('position', [0.05  0.25 0.425 0.5]);
h.ax(2) = axes('position', [0.525 0.25 0.425 0.5]);

% Plot

axes(h.ax(1));
contourf(x,y,z,clev);
cb = colorbar('eastoutside');
colormap(h.ax(1), jet);
title(h.ax(1), 'contourf');
```


![](./readmeExtras/README_01.png)

Using contourfcmap, we can set the colors explictly, so it's much easier to tell exactly which color corresponds to which value.



```matlab
axes(h.ax(2));
hc = contourfcmap(x,y,z,clev,jet(12), ...
'lo', [.8 .8 .8], 'hi', [.2 .2 .2], 'cbarloc', 'eastoutside', ...
'method', 'calccontour');
title(h.ax(2), 'contourfcmap');
```


![](./readmeExtras/README_02.png)

If you prefer, you can set the colorbar to show the contour intervals evenly-spaced, even if the values aren't.



```matlab
delete(hc.cbax);
set(h.ax(2), 'position', [0.525 0.25 0.425 0.5]);

hc = contourfcmap(x,y,z,clev,jet(12), ...
'lo', [.8 .8 .8], 'hi', [.2 .2 .2], 'cbarloc', 'eastoutside', ...
'method', 'calccontour', 'evencb', true);
title(h.ax(2), 'contourfcmap');
```


![](./readmeExtras/README_03.png)

**2014b updates**


In the 2014b release, Matlab introduced new handle graphics objects. Contour objects got a complete makeover with this update, and I had to completely rewrite contourfcmap to accomodate these changes.  The new version offers two methods to create the contours.


`recolor`: This is similar to the old pre-2014b version of contourfcmap, and is the default option.  It takes the existing contour objects and simply recolors the children patches (pre-2014b) or TriangleStrips (2014b+).  This method is quick, efficient, and requires no outside toolboxes.  But it has a few drawbacks:



- It inherits some weaknesses from contourf itself in the way it handles regions lower than the lowest specified contour or higher than the highest specified contour.  Sometimes these regions are left empty, and whether or not the regions are shaded vary by Matlab release.
- In 2014b+, the recolor option has a tendency to get undone very easily.  Changing the colormap, adding a colorbar, and most frustratingly, printing figures to file (using anything that relies on the print command, including print, saveas, export_fig, etc) undoes it.  The only way I've found to save these figures is by doing a screenshot (i.e. `im = frame2im(getframe(gcf); imwrite(im, 'myfile.png', 'png')`).
- Like contourf itself, the contour object handle returned by pre-2014b and 2014b+ are different classes (double contourgroup object vs graphics object Contour object) with different properties, and in some applications this inconsistency can make programming across Matlab versions difficult.
- In 2014b+, the recoloring processes applies to undocumented grpahics object properties, so this code could break in future updates.

`calccontour`: With this method, I calculate contour patches and lines from scratch, and plot these lines and patches in place of a contour object.  This method was originally added as a stopgap solution for the 2014b release, before I managed to explore the undocumented properties of the new contour objects.  It's more resource intensive than the original, and requires the Mapping Toolbox (for its many polygon-geometry-related functions), but it solves the problems encountered with the recolor method.


This example compares the two options. The results you see will be slightly different in you're running a R2014a or earlier.  Note that to get the recolor method to show up, I have to create those two plots last (otherwise the colorbars will revert the colors).



```matlab
[x,y] = meshgrid(linspace(0,1,100));
z = peaks(100);
clev = [-5 -3 -2:.5:2 3 5];
cmap = cool(length(clev)-1);

v = ver('matlab');
v = regexprep(v.Release, '[\(\)]', '');
clf;
ax(1) = subplot(2,2,1);
hh = contourfcmap(x,y,z,clev,cmap, ...
'lo', [0 0 1], ...
'hi', [1 1 0], ...
'method', 'calccontour', ...
'cbarloc', 'eastoutside');
title(ax(1), {sprintf('calccontour method (%s)', v), 'Enclosed too-low region'});
set(hh.cbax, 'fontsize', 8);

clev = 1:5;
cmap = cool(length(clev)-1);

ax(3) = subplot(2,2,3);
hh = contourfcmap(x,y,z,clev,cmap, ...
'lo', [0 0 1], ...
'hi', [1 1 0], ...
'method', 'calccontour', ...
'cbarloc', 'eastoutside');
title(ax(3), {sprintf('calccontour method (%s)', v), 'Open too-low region'});
set(hh.cbax, 'fontsize', 8);

ax(4) = subplot(2,2,4);
hh = contourfcmap(x,y,z,clev,cmap, ...
'lo', [0 0 1], ...
'hi', [1 1 0], ...
'method', 'recolor');
title(ax(4), {sprintf('recolor method (%s)', v), 'Open too-low region'});

clev = [-5 -3 -2:.5:2 3 5];
cmap = cool(length(clev)-1);

ax(2) = subplot(2,2,2);
hh = contourfcmap(x,y,z,clev,cmap, ...
'lo', [0 0 1], ...
'hi', [1 1 0], ...
'method', 'recolor');
title(ax(2), {sprintf('recolor method (%s)', v), 'Enclosed too-low region'});

set(ax([2 4]), 'dataaspectratio', ax(1,1).DataAspectRatio)

% Take snapshot of the figure and display as image, since the recolor
% method will reset on printing on my computer

im = getframe(gcf);
clf;
axes('position', [0 0 1 1]);
image(im.cdata);
set(gca, 'visible', 'off');
```


![](./readmeExtras/README_04.png)


## Contributions


Community contributions to this package are welcome!


To report bugs, please submit [an issue](https://github.com/kakearney/example-pkg/issues) on GitHub and include:



- your operating system
- your version of Matlab and all relevant toolboxes (type `ver` at the Matlab command line to get this info)
- code/data to reproduce the error or buggy behavior, and the full text of any error messages received

Please also feel free to submit enhancement requests, or to send pull requests (via GitHub) for bug fixes or new features.


I do monitor the MatlabCentral FileExchange entry for any issues raised in the comments, but would prefer to track issues on GitHub.



<sub>[Published with MATLAB R2016a]("http://www.mathworks.com/products/matlab/")</sub>
