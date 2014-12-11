function varargout = contourfcmap(x,y,z,clev,cmap,lo,hi,cbarloc,evencb)
%CONTOURFCMAP Filled contour plot with specified colors
%
% h = contourfcmap(x,y,z,clev,cmap,lo,hi,cbarloc)
% h = contourfcmap(x,y,z,clev,cmap)
% h = contourfcmap(x,y,z,clev,cmap,lo)
% h = contourfcmap(x,y,z,clev,cmap,lo,hi)
% 
% This function creates a shaded contour map, similar to that created by
% the contourf function.  However, the relationship between a contourf plot
% and its colormap (i.e. exactly which color corresponds to each contour
% interval), can often be confusing and inconsistent, in my opinion.  This
% function instead allows the user to specify exactly which colors to use
% in each interval, and also to choose colors for regions that exceed the
% contour line limits.
%
% Update 7/18/2014: Due to changes in the way contourf objects are created
% in R2014b, I have rewritten this function completely to manually
% calculate child lines and patches that make up a filled contour plot.  To
% do this, I've relied on several of the polygon calculation functions
% found in the Mapping Toolbox.  If you don't have the Mapping Toolbox and
% are running an older version of Matlab, it will fall back on the old
% version of this function.  Unfortunately I haven't had the time write
% out a solution for R2014b+ without the Mapping Toolbox.
%
% Input variables
%
%   x:          x-coordinates of grid, following size restrictions of surf
%               command 
%
%   y:          y-coordinates of grid, following size restrictions of surf
%               command 
%
%   z:          array of data to be contoured, following size restritions
%               of surf command 
%
%   clev:       vector of length n, contour levels, must be monotonically
%               increasing 
%
%   cmap:       n-1 x 3 colormap array, specifying colors used to shade
%               each contour interval 
%
%   lo:         1 x 3 colormap array, specifying color used for all data
%               falling below the first contour level.  If not included or
%               empty, default will be to white. 
%
%   hi:         1 x 3 colormap array, specifying color used for all data
%               falling above the last contour level.  If not included or
%               empty, default will be to white.
%
%   cbarloc:    string specifying colorbar location (see colorbar),  or a 
%               1 x 4 position vector for colorbar.  If not included, no
%               colorbar will be created.  Note that the colorbar created
%               is not a colorbar object, but just an axis with plotted
%               patches; it is for labeling only and is not linked to the
%               "peer axis" in any way.
%
%   evencb:     logiscal scalar.  If true, intervals on the colorbar will
%               be evenly spaced, regardless of value.  If false, ticks
%               will correspond to clev values. If not included or empty,
%               default is false.
%
% Output variables:
%
%   h:          1 x 1 structure:
%
%               p:      handles to patch objects
%               
%               l:      handles to line objects
%
%               cbax:   handle to axis of colorbar
%
% Output variables (if running pre2014b w/o Mapping Toolbox):
%
%   h:          1 x 1 structure
%
%               c:      contour matrix for filled contour plot
%
%               h:      handle to contourgroup
%
%               cbax:   handle to axis of colorbar
%
% Example:
%
% [x,y] = meshgrid(linspace(0,1,100));
% z = peaks(100);
% contourfcmap(x,y,z,[-5 -3 -2:.5:2 3 5],jet(12), ...
%              [.8 .8 .8], [.2 .2 .2], 'eastoutside')
% 

% Copyright 2010-2014 Kelly Kearney

% New version is a little more stable than the old one, but requires the
% Mapping Toolbox for all the polygon operations.  Fall back on the old
% version if Mapping Toolbox not found.

if isempty(ver('map'))
    if verLessThan('matlab', '8.1.0')
        varargout = contourfcmap_pre2014b(x,y,z,clev,cmap,lo,hi,cbarloc);
    else
        error('Sorry, this function currently requires the Mapping Toolbox to run in R2014b or later');
    end
end

%------------------------
% Parse input
%------------------------

narginchk(5,9);

% Check clev

if ~isvector(clev)
    error('clev must be a vector');
end
nlev = length(clev);
clev = reshape(clev, 1, []); 

% Check cmap

if size(cmap,2) ~=3 || size(cmap,1) ~= (nlev-1) || any(cmap(:)<0|cmap(:)>1)
    error('cmap must be nlev-1 x 3 colormap array');
end

% Colors for too-low and too-high patches

if nargin < 6 || isempty(lo)
    lo = [1 1 1];
end

if nargin < 7 || isempty(hi)
    hi = [1 1 1];
end

if ~isequal(size(lo),[1 3])
    error('lo must be 1 x 3 colormap array');
end

if ~isequal(size(hi),[1 3])
    error('hi must be 1 x 3 colormap array');
end
    
% Colorbar location

if nargin >= 8 && ~isempty(cbarloc)
    pos = {'north', 'south', 'east', 'west', 'northoutside', 'southoutside', 'eastoutside', 'westoutside'};
    if ischar(cbarloc)
        if ~any(strcmp(lower(cbarloc), pos))
            error('Unrecognizd colorbar position');
        end
    elseif ~isequal(size(cbarloc), [1 4]) || any(cbarloc > 1 | cbarloc < 0)
        error('cbarloc must be position string or  1 x 4 normalized position');
    end
    showcb = true;
else
    showcb = false;
end

if nargin < 9 || isempty(cbarloc)
    evencb = false;
end

% Axis

ax = gca;

%------------------------
% Contour calculations
%------------------------

% Get contour lines

if isvector(x)
    x = reshape(x, 1, []);
end
if isvector(y)
    y = reshape(y, [], 1);
end

S = contourcs(x(1,:),y(:,1),z,clev);

% For the lines that hit the boundaries, figure out if corners need to be
% included

isclosed = arrayfun(@(A) isequal([A.X(1) A.Y(1)], [A.X(end) A.Y(end)]), S);

xlim = minmax(x);
ylim = minmax(y);

xcorner = xlim([1 2 2 1 1]);
ycorner = ylim([2 2 1 1 2]);

for ii = find(~isclosed)'
    
    % Which wall is each point on?
    
    xtmp = S(ii).X([1 end]);
    ytmp = S(ii).Y([1 end]);

    wall = zeros(1,2);
    wall(xtmp == xlim(1)) = 1; % left
    wall(ytmp == ylim(2)) = 2; % top
    wall(xtmp == xlim(2)) = 3; % right
    wall(ytmp == ylim(1)) = 4; % bottom
    
    if wall(1) == wall(2) % Same wall, just connect
        S(ii).X = [S(ii).X S(ii).X(1)];
        S(ii).Y = [S(ii).Y S(ii).Y(1)];
    else
        tbl = {...
            [1 2]   1
            [2 1]   1
            [1 3]   [2 1]
            [3 1]   [1 2]
            [1 4]   4
            [4 1]   4
            [2 3]   2
            [3 2]   2
            [2 4]   [3 2]
            [4 2]   [2 3]
            [3 4]   3
            [4 3]   3};
        
        
        [tf, loc] = ismember(wall, cat(1, tbl{:,1}), 'rows');
        S(ii).X = [S(ii).X xcorner(tbl{loc,2}) S(ii).X(1)];
        S(ii).Y = [S(ii).Y ycorner(tbl{loc,2}) S(ii).Y(1)];
        
    end
end

% Triangulate patches, eliminating overlap

[xc,yc] = poly2cw({S.X}, {S.Y});

xc = [xc xcorner];
yc = [yc ycorner];

[xnew, ynew] = multiplepolyint(xc,yc);

np = length(xnew);

[f,v] = deal(cell(np,1));
[xsamp, ysamp] = deal(zeros(np,1));
for ip = 1:length(xnew)
    [f{ip},v{ip}] = poly2fv(xnew{ip}, ynew{ip});
    
    vx = v{ip}(:,1);
    vy = v{ip}(:,2);
    xsamp(ip) = mean(vx((f{ip}(1,:))));
    ysamp(ip) = mean(vy((f{ip}(1,:))));
    
end

% There's probably a more elegant way to figure out which color goes where,
% but this works for now

val = interp2(x,y,z,xsamp,ysamp);
[~, lev] = histc(val, clev);

cmap2 = [cmap; lo; hi];
lev(val < clev(1)) = size(cmap2,1)-1;
lev(val > clev(end)) = size(cmap2,1);

% Plot contour lines and patches

hold(ax, 'on');
for ip = 1:np
    hp(ip) = patch('faces', f{ip}, 'vertices', v{ip}, 'facecolor', cmap2(lev(ip),:), 'edgecolor', 'none');
end
for il = 1:length(xnew)
    hl(il) = line(xnew{il}, ynew{il}, 'color', 'k');
end

%------------------------
% Create colorbar
%------------------------

if showcb
    if evencb
        clevcbar = linspace(0,1,length(clev));
    else
        clevcbar = reshape(clev,1,[]);
    end
        
    height = 0.1 * (max(clevcbar) - min(clevcbar));
    y1 = [clevcbar(1)-height; clevcbar(1); clevcbar(1); clevcbar(1)-height; clevcbar(1)-height]; 
    y2 = [clevcbar(end); clevcbar(end)+height; clevcbar(end)+height; clevcbar(end); clevcbar(end)];

    yp = [y1 [clevcbar(1:end-1); clevcbar(2:end); clevcbar(2:end); clevcbar(1:end-1); clevcbar(1:end-1)] y2];
    xp = [0;0;1;1;0] * ones(1,nlev+1);
    cp = [lo; cmap; hi];   
    cp = permute(cp, [3 1 2]);

    if ~ischar(cbarloc)
        cbarcoord = cbarloc;
        hascbcoord = true;
        if cbarcoord(3)-cbarcoord(1) < cbarcoord(4)-cbarcoord(2)
            isvert = true;
            cbarloc = 'east';
        else
            isvert = false;
            cbarloc = 'south';
        end
    else
        hascbcoord = false;
        if any(strcmp(pos([3 4 7 8]), lower(cbarloc)))
            isvert = true;
        else
            isvert = false;
        end
    end

    if ~isvert
        tmp = xp;
        xp = yp;
        yp = tmp;
    end

    
    cbax = colorbar(cbarloc);
    drawnow;
    axpos = get(ax, 'position');
    cbpos = get(cbax, 'position');

    delete(cbax);
    drawnow; 
    set(ax, 'position', axpos);
    cbax = axes('position', cbpos);
    set(cbax, 'userdata', 'contourfcolorbar');

    patch(xp,yp,cp);
    if isvert
        set(cbax, 'ytick', clevcbar, 'ylim', minmax(yp), 'xlim', [0 1], 'xtick', []);
        if evencb
            set(cbax, 'yticklabel', strtrim(cellstr(num2str(clev'))));
        end
    else
        set(cbax, 'xtick', clevcbar, 'xlim', minmax(xp), 'ylim', [0 1], 'ytick', []);
        if evencb
            set(cbax, 'xticklabel', strtrim(cellstr(num2str(clev'))));
        end
    end
    
    switch cbarloc
        case {'east', 'westoutside'}
            set(cbax, 'yaxislocation', 'left');
        case {'eastoutside', 'west'}
            set(cbax, 'yaxislocation', 'right');
        case {'north', 'southoutside'}
            set(cbax, 'xaxislocation', 'bottom');
        case {'northoutside', 'south'}
            set(cbax, 'xaxislocation', 'top');
    end
    
end

if showcb && hascbcoord
    set(cbax, 'position', cbarcoord);
end

%------------------------
% Output
%------------------------

hndl.p = hp;
hndl.l = hl;
if showcb
    % Return focus to axis (needed if colorbar created), but keep colorbar
    % on top
    hndl.cbax = cbax;
    axes(ax);
    uistack(cbax, 'top');
end

if nargout > 0
    varargout{1} = hndl;
end



function a = minmax(b)
a = [min(b(:)) max(b(:))];




