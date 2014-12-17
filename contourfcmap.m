function varargout = contourfcmap(x,y,z,clev,cmap,varargin);% lo,hi,cbarloc,evencb)
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
% Optional input variables (passed as parameter/value pairs)
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
%               "peer axis" in any way. Default is no colorbar.
%
%   evencb:     logiscal scalar.  If true, intervals on the colorbar will
%               be evenly spaced, regardless of value.  If false, ticks
%               will correspond to clev values. If not included or empty,
%               default is false.
%
%   method:     string specifying calculation method:
%
%               'recolor':      the default, create a contourf object, and
%                               then change the color properties of the
%                               underlying components. Note: In 2014b,
%                               recolor will not persist when saving to
%                               file via anything that uses print.m
%                               (including print, saveas, export_fig, etc)
%
%               'calccontour':  Create new patch and line objects from
%                               scratch.  This method requires the Mapping
%                               Toolbox.  It can be beneficial if you want
%                               more consistency between output regardless
%                               of which version of Matlab is being run.
%                               It also properly fills regions falling
%                               below the lowest contour level and above
%                               the highest contour level, which isn't
%                               always the case with contourf-generated
%                               contour objects.
%
% Output variables:
%
%   h:          1 x 1 structure, with the following fields (some may not be
%               present, depending on input options)
%
%               c:      contour matrix for filled contour plot ('recolor'
%                       method only)
%
%               h:      handle to contourgroup or contour object ('recolor'
%                       method only)
%
%               p:      handles to patch objects (the filled contour
%                       regions) ('calccontour' method only)
%               
%               l:      handles to line objects (the contour lines)
%                       ('calccontour' method only) 
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

%------------------------
% Parse input
%------------------------

% New syntax uses parameter/value options, but we'll still accept the old
% syntax 

Opt.lo = [1 1 1];
Opt.hi = [1 1 1];
Opt.cbarloc = [];
Opt.evencb = false;
Opt.method = 'recolor';

isc = cellfun(@ischar, varargin);
if ~mod(length(varargin),2) & all(isc(1:2:end))
    
    % New syntax
    
    Opt = parsepv(Opt, varargin);
    
else
    % Old syntax
    
    narginchk(5,9);
    
    fld = fieldnames(Opt);
    for ii = 1:length(varargin)
        Opt.(fld{ii}) = varargin{ii};
    end
    
end
    
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

if isempty(Opt.lo)
    Opt.lo = [1 1 1];
end

if isempty(Opt.hi)
    Opt.hi = [1 1 1];
end

if ~isequal(size(Opt.lo),[1 3])
    error('Opt.lo must be 1 x 3 colormap array');
end

if ~isequal(size(Opt.hi),[1 3])
    error('Opt.hi must be 1 x 3 colormap array');
end
    
% Colorbar location

if nargin >= 8 && ~isempty(Opt.cbarloc)
    pos = {'north', 'south', 'east', 'west', 'northoutside', 'southoutside', 'eastoutside', 'westoutside'};
    if ischar(Opt.cbarloc)
        if ~any(strcmp(lower(Opt.cbarloc), pos))
            error('Unrecognizd colorbar position');
        end
    elseif ~isequal(size(Opt.cbarloc), [1 4]) || any(Opt.cbarloc > 1 | Opt.cbarloc < 0)
        error('cbarloc must be position string or  1 x 4 normalized position');
    end
    showcb = true;
else
    showcb = false;
end

if nargin < 9 || isempty(Opt.cbarloc)
    Opt.evencb = false;
end

% Axis

ax = gca;


%------------------------
% Create colorbar
%------------------------

% Need to create the colorbar up front, otherwise doing so will mess up
% recoloring

if showcb
    if Opt.evencb
        clevcbar = linspace(0,1,length(clev));
    else
        clevcbar = reshape(clev,1,[]);
    end
        
    height = 0.1 * (max(clevcbar) - min(clevcbar));
    y1 = [clevcbar(1)-height; clevcbar(1); clevcbar(1); clevcbar(1)-height; clevcbar(1)-height]; 
    y2 = [clevcbar(end); clevcbar(end)+height; clevcbar(end)+height; clevcbar(end); clevcbar(end)];

    yp = [y1 [clevcbar(1:end-1); clevcbar(2:end); clevcbar(2:end); clevcbar(1:end-1); clevcbar(1:end-1)] y2];
    xp = [0;0;1;1;0] * ones(1,nlev+1);
    cp = [Opt.lo; cmap; Opt.hi];   
    cp = permute(cp, [3 1 2]);

    if ~ischar(Opt.cbarloc)
        cbarcoord = Opt.cbarloc;
        hascbcoord = true;
        if cbarcoord(3)-cbarcoord(1) < cbarcoord(4)-cbarcoord(2)
            isvert = true;
            Opt.cbarloc = 'east';
        else
            isvert = false;
            Opt.cbarloc = 'south';
        end
    else
        hascbcoord = false;
        if any(strcmp(pos([3 4 7 8]), lower(Opt.cbarloc)))
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

    
    cbax = colorbar(Opt.cbarloc);
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
        if Opt.evencb
            set(cbax, 'yticklabel', strtrim(cellstr(num2str(clev'))));
        end
    else
        set(cbax, 'xtick', clevcbar, 'xlim', minmax(xp), 'ylim', [0 1], 'ytick', []);
        if Opt.evencb
            set(cbax, 'xticklabel', strtrim(cellstr(num2str(clev'))));
        end
    end
    
    switch Opt.cbarloc
        case {'east', 'westoutside'}
            set(cbax, 'yaxislocation', 'left');
        case {'eastoutside', 'west'}
            set(cbax, 'yaxislocation', 'right');
        case {'north', 'southoutside'}
            set(cbax, 'xaxislocation', 'bottom');
        case {'northoutside', 'south'}
            set(cbax, 'xaxislocation', 'top');
    end
    
    hout.cbax = cbax;
    
end
    
%------------------------
% Contour calculations
%------------------------

axes(ax);

switch Opt.method
    case 'recolor'
        if verLessThan('matlab', '8.4.0')
            
            [c,h] = contourf(x,y,z,clev);

            hpatch = get(h, 'children');

            cdata = cell2mat(get(hpatch, 'cdata'));

            % Mark too-high contours

            isabove = cdata == max(clev);

            % Distinguish between too-lo contours and NaN contours

            if ~any(isnan(z(:)))
                isbelow = isnan(cdata);
            else
                warning('NaN in data; haven''t figured NaN vs lo yet');
                isbelow = isnan(cdata);
            end

            level = get(h, 'LevelList');
            if ~isequal(level, clev)
                error('oops, something new with levels, check this');
            end

            [tf, idx] = ismember(cdata, clev);

            isweirdextra = idx == 0 & ~isbelow & ~isabove; % Why do these show up?!

            for ip = 1:length(hpatch)
                if isbelow(ip)
                        set(hpatch(ip), 'facecolor', Opt.lo);
                elseif isabove(ip)
                        set(hpatch(ip), 'facecolor', Opt.hi);
                elseif isweirdextra(ip)
                    dist = cdata(ip) - clev;
                    dist(dist<0) = Inf;
                    [blah, imin] = min(dist);
                    set(hpatch(ip), 'facecolor', cmap(imin,:));
                else
                    set(hpatch(ip), 'facecolor', cmap(idx(ip),:));
                end
            end
            
            hout.h = h;
            hout.c = c;
            
        else
            
            [c,h] = contourf(x,y,z,clev);
            
            drawnow;
            Fp = h.FacePrims;
            np = length(Fp);
            xyz = zeros(3, np);
            for ii = 1:np
                
                % Can't seem to find the property that links each
                % TrangleStrip to a particular contour level, and hence
                % color.  So I'm just going to find in as I do in the
                % calccontour method, by interpolating the value at a point
                % inside the first triangle.
                
                xyz(:,ii) = mean(Fp(ii).VertexData(:,Fp(ii).StripData(1:3)),2);

            end
            val = interp2(x,y,z,xyz(1,:), xyz(2,:));
            [n, idx] = histc(val, clev);
            idx(val < clev(1)) = nlev;
            idx(val > clev(end)) = nlev + 1;
            
            cdata = [cmap; Opt.lo; Opt.hi] * 255;
            
            newcol = cdata(idx,:);
            
            for ii = 1:np
                Fp(ii).ColorData = uint8([newcol(ii,:) 255]');
            end
            
            hout.h = h;
            hout.c = c;
            
        end
        
        
    case 'calccontour'
        
        if isempty(ver('map'))
            error('The calccontour method requires the Mapping Toolbox');
        end
        
        % Get contour lines

        if isvector(x)
            x = reshape(x, 1, []);
        end
        if isvector(y)
            y = reshape(y, [], 1);
        end

        S = contourcs(x(1,:),y(:,1),z,clev);

        % For the lines that hit the boundaries, figure out if corners need
        % to be included 

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

        % There's probably a more elegant way to figure out which color
        % goes where, but this works for now 

        val = interp2(x,y,z,xsamp,ysamp);
        [~, lev] = histc(val, clev);

        cmap2 = [cmap; Opt.lo; Opt.hi];
        lev(val < clev(1)) = size(cmap2,1)-1;
        lev(val > clev(end)) = size(cmap2,1);

        % Plot contour lines and patches

        hout.p = gobjects(np,1);
        hout.l = gobjects(np,1);
        
        hold(ax, 'on');
        for ip = 1:np
            hout.p(ip) = patch('faces', f{ip}, 'vertices', v{ip}, 'facecolor', cmap2(lev(ip),:), 'edgecolor', 'none');
        end
        for il = 1:length(xnew)
            hout.l(il) = line(xnew{il}, ynew{il}, 'color', 'k');
        end    
        
end




if showcb && hascbcoord
    set(cbax, 'position', cbarcoord);
end

%------------------------
% Output
%------------------------

if nargout > 0
    varargout{1} = hout;
end

% hndl.p = hp;
% hndl.l = hl;
% if showcb
%     % Return focus to axis (needed if colorbar created), but keep colorbar
%     % on top
%     hndl.cbax = cbax;
%     axes(ax);
%     uistack(cbax, 'top');
% end
% 
% if nargout > 0
%     varargout{1} = hndl;
% end



function a = minmax(b)
a = [min(b(:)) max(b(:))];




