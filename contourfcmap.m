function varargout = contourfcmap(x,y,z,clev,cmap,varargin);% lo,hi,cbarloc,evencb)
%CONTOURFCMAP Filled contour plot with specified colors
%
% h = contourfcmap(x,y,z,clev,cmap,lo,hi,cbarloc)
% h = contourfcmap(x,y,z,clev,cmap)
% h = contourfcmap(x,y,z,clev,cmap,lo)
% h = contourfcmap(x,y,z,clev,cmap,lo,hi)
% h = contourfcmap(x,y,z,clev,cmap, param1, val1, ...)
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
%   evencb:     logical scalar.  If true, intervals on the colorbar will
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
%   flag:       logical flag indicating whether to use the fast version of
%               multiplepolyint.m (only applicable if method is
%               'calccontour').  In plots with a large number of contour
%               lines, this significantly speeds up the calculation.
%               However, it works by directly calling a private mex
%               function in the Mapping Toolbox, so it might be a bit less
%               stable than the slower version, and may break in future
%               releases. Default is true.
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

% Copyright 2010-2016 Kelly Kearney

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
Opt.flag = true;

isc = cellfun(@ischar, varargin);
if ~mod(length(varargin),2) && all(isc(1:2:end))
    
    % New syntax
    
    Opt = parsepv(Opt, varargin);
    
else
    % Old syntax
    
    if verLessThan('matlab', 'R2011b')
        error(nargchk(5,9,nargin));
    else
        narginchk(5,9);
    end
    
    fld = fieldnames(Opt);
    for ii = 1:length(varargin)
        Opt.(fld{ii}) = varargin{ii};
    end
    
end
    
% Check version and dependencies

if strcmp(Opt.method, 'calccontour')
    if verLessThan('matlab', '9.3.0') % R2017b
        vmap = ver('map');
        if isempty(vmap)
            error('The calccontour method requires either Matlab R2017b or later, or the Mapping Toolbox');
        end
        cmethod = 'map';
    else
        cmethod = 'polyshape';
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
            error('Unrecognized colorbar position');
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
% X/Y pattern
%------------------------

% X and Y can either be vectors, or matrices, and the matrices can be
% irregular.  Here I check whether the input coordinate arrays are vectors,
% meshgrid-style matrices, ndgrid-style matrices, or an irregular (but
% still structured!) grid.

isvec = isvector(x) && ...
        isvector(y) && ...
        ((length(x)==size(z,2) && ...
         length(y)==size(z,1)) || ...
        (length(x)==size(z,1) && ...
         length(y)==size(z,2)));
    
ismgrid = isequal(size(x),size(y),size(z)) && ...
          ~any(reshape(bsxfun(@minus, x, x(1,:)),[],1)) && ...
          ~any(reshape(bsxfun(@minus, y, y(:,1)),[],1));
isngrid = isequal(size(x),size(y),size(z)) && ...
          ~any(reshape(bsxfun(@minus, x, x(:,1)),[],1)) && ...
          ~any(reshape(bsxfun(@minus, y, y(1,:)),[],1));
      
isigrid = isequal(size(x),size(y),size(z)) && ...
          ~ismgrid && ...
          ~isngrid;
      
if ~(isvec || ismgrid || isngrid || isigrid)
    htmp = contourf(x,y,z);
    
    % If this works, it means I've found an input format I hadn't
    % considered in my code.  Might work, might not.  Email me if you hit
    % this.
    
    warning('You''ve found an accepatable x/y/z format that I haven''t tested... may have unexpected results');
    delete(htmp);
end
       
% Convert all rectilinear input format to vector x and y, with z rows
% corresponding to y and columns to x.  Irregular grids will stay as they
% are.

if isigrid
    irrflag = true; % Flag for irregular grid
else
    irrflag = false;
    if isvec
        if size(z,2) ~= length(x)
            z = z';
        end
    elseif ismgrid
        x = x(1,:);
        y = y(:,1);
    elseif isngrid
        x = x(:,1);
        y = y(1,:);
        z = z';
    end
end

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

            if isscalar(hpatch)
                cdata = get(hpatch, 'cdata');
            else
                cdata = cell2mat(get(hpatch, 'cdata'));
            end

            % Mark too-high contours

            isabove = cdata == max(clev);

            % Distinguish between too-lo contours and NaN contours

            if ~any(isnan(z(:)))
                isbelow = isnan(cdata);
                isn = false(size(cdata));
            else
                idxtmp = find(isnan(cdata));
                xy = get(hpatch(idxtmp), {'xdata','ydata'});
                xy = cellfun(@(x) x(1), xy);
                if irrflag
                    F = scatteredInterpolant(x(:), y(:), z(:));
                    ztmp = F(xy(:,1),xy(:,2));
                else
                    ztmp = interp2(x,y,z,xy(:,1),xy(:,2));
                end
                
                isn = false(size(cdata));
                isn(idxtmp(isnan(ztmp))) = true;
                
                isbelow = isnan(cdata) & ~isn;
                
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
                elseif isn(ip)
                    set(hpatch(ip), 'facecolor', [1 1 1]);
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
            
            % I can't seem to find the property that links each
            % TriangleStrip to a particular contour level (and hence the
            % appropriate color).  So I'm going to have to determine that
            % by checking the value in one triangle per strip.
            
            vd = get(Fp, 'VertexData');
            
            % Vertex data is sometimes scaled.  I'm not quite sure what
            % determines this.  My best bet so far is to just check to see
            % if the vd data is in the axis limits.  If not, assume it's
            % scaled.  (I used to also check that potentially-scaled values
            % were between 0 and 1, but I've hit a few examples where
            % scaled values exceed that by a few thousandths... still not
            % very clear on what's going on here).
            
            if ~iscell(vd)
                vd = {vd};
            end
            lims = [min(cat(2, vd{:}), [], 2) max(cat(2, vd{:}), [], 2)];
            axlims = get(ax, {'xlim', 'ylim','zlim'});
            axlims = cat(1, axlims{:});
            
            isin = bsxfun(@ge, lims, axlims(:,1)) & bsxfun(@le, lims, axlims(:,2));
%             is01 = lims >= 0 & lims <= 1;
            
            if ~all(isin(:))
%                 if all(is01(:)) % Scaled
                s = diff(axlims,1,2);
                o = axlims(:,1);
                vd = cellfun(@(v) bsxfun(@plus, bsxfun(@times, double(v), s), o), vd, 'uni', 0); 
            end

            
            sd = get(Fp, 'StripData');
            if ~iscell(sd)
                sd = {sd};
            end
            xyz = cellfun(@(v,s) double(v(:,s(1:end-1))), vd, sd, 'uni', 0);
            idx = zeros(np,1);
            
            if irrflag
%                 error('Still figuring this one out for irregular grids');
                F = scatteredInterpolant(x(:), y(:), z(:));
            else
                F = griddedInterpolant({x,y},z');
            end
               
            for ii = 1:np
                tmp = F(xyz{ii}(1,:), xyz{ii}(2,:));
%                 tmp = interp2(x,y,z,xyz{ii}(1,:), xyz{ii}(2,:));
                [ntmp, bin] = histc(tmp, [-Inf clev Inf]);
                [~,idx(ii)] = max(ntmp);
            end
            
            cdata = [Opt.lo; cmap; Opt.hi] * 255;
            newcol = cdata(idx,:);
            
            for ii = 1:np
                Fp(ii).ColorData = uint8([newcol(ii,:) 255]');
            end
            
            hout.h = h;
            hout.c = c;
            
        end
        
        
    case 'calccontour'
        
        % Get contour lines

        if isvector(x)
            x = reshape(x, 1, []);
        end
        if isvector(y)
            y = reshape(y, [], 1);
        end
        
        % Calculate contour lines
        
        if irrflag
            [c, htmp] = contour(x,y,z,clev);
            S = contourcs(c, 'cmat');
            delete(htmp);
            
            F = scatteredInterpolant(x(:), y(:), z(:));
        else
            S = contourcs(x,y,z,clev);
            F = griddedInterpolant({x,y},z');
        end
        
        % Coordinates for wall
        
        if isvector(x)
            xlim = [min(x) max(x)];
            ylim = [min(y) max(y)];
            xwall = xlim([1 1 2 2 1]);
            ywall = ylim([1 2 2 1 1]);
        else
            xwall = [x(end:-1:1,1)' x(2,:), x(2:end,end)', x(end,end-1:-1:1)];
            ywall = [y(end:-1:1,1)' y(2,:), y(2:end,end)', y(end,end-1:-1:1)];
        end
        
        % If there are NaNs in the dataset, calculate contour lines
        % for the NaN regions and for the inpainted data
        
        nflag = any(isnan(z(:)));
        if nflag
            
            % Inpaint the NaN holes using nearest neighbor values
            
            if irrflag
                [zi, A] = fillnan(z, {x,y});
                F = scatteredInterpolant(x(:), y(:), zi(:));
                Fn = scatteredInterpolant(x(:), y(:), double(isnan(z(:))));
            else
                [zi, A] = fillnan(z', {x, y});
                zi = zi';
                F = griddedInterpolant({x,y},zi');
                Fn = griddedInterpolant({x,y},double(isnan(z')));
            end
            
            % Calculate the polygons that surround all NaNs
            
            [xv,yv] = voronoigrid(x,y,xwall,ywall);
            
            xn = xv(isnan(z));
            yn = yv(isnan(z));

            [xn,yn] = poly2cw(xn,yn);
            [xn, yn] = polyjoin(xn, yn);
            
            % Calculate contours for now-filled dataset
            
            if irrflag
                [c, htmp] = contour(x,y,zi,clev);
                S = contourcs(c, 'cmat');
                delete(htmp);     
            else
                S  = contourcs(x, y, zi, clev);
            end
            
        end
        
       
        % Extend contours to wall

        S = extendtowall(S,xwall,ywall);
        [xc,yc] = poly2cw({S.X}, {S.Y});
        
        % Remove overlap and triangulate
        
        if nflag
            
            [xc, yc] = removeoverlap(xc, yc, xwall, ywall, Opt.flag);
            for ic = 1:length(xc)
                [xc{ic}, yc{ic}] = polybool('-', xc{ic}, yc{ic}, xn, yn);
            end
            isemp = cellfun(@isempty, xc);
            xc = xc(~isemp);
            yc = yc(~isemp);
            
            [f,v,lev] = poly2faces(xc, yc, F, clev);
                
        else
            [xc, yc] = removeoverlap(xc, yc, xwall, ywall, Opt.flag);
            [f,v,lev] = poly2faces(xc, yc, F, clev);
        end
        np = length(f);
        
        
        % Plot contour lines and patches

        cmap2 = [Opt.lo; cmap; Opt.hi; 1 1 1];

        hout.p = gobjects(np,1);
        hout.l = gobjects(np,1);
        
        hold(ax, 'on');
        for ip = 1:np
            hout.p(ip) = patch('faces', f{ip}, 'vertices', v{ip}, 'facecolor', cmap2(lev(ip),:), 'edgecolor', 'none');
        end
        [hout.c, hout.l] = contour(x,y,z,clev,'k'); 
        
        % A few axes properties changes similar to those made by contourf
        
        axis(ax, 'tight');
        set(ax, 'layer', 'top', 'box', 'on');
        
end

if showcb && hascbcoord
    set(cbax, 'position', cbarcoord);
end

if showcb
    uistack(cbax, 'top'); % Has to be done after plotting, but resets color in R2014b-recolor... aaargh!
end

%------------------------
% Output
%------------------------

if nargout > 0
    varargout{1} = hout;
end

%*********** Subfunctions ************************************************

%------------
% Minmax
%------------

% function a = minmax(b)
% a = [min(b(:)) max(b(:))];

%--------------------
% Extend contours to 
% "walls" of the grid
%--------------------

function S = extendtowall(S, xwall, ywall)

isclosed = arrayfun(@(A) isequal([A.X(1) A.Y(1)], [A.X(end) A.Y(end)]), S);

[arclen, seglen] = arclength(xwall(:), ywall(:));
twall = [0; cumsum(seglen)]./arclen;

% plot(xwall, ywall);

for ii = find(~isclosed)'

    % Where along the wall does the open part of the contour hit?
    
    [xy,d,t] = distance2curve([xwall(:) ywall(:)], [S(ii).X([end 1])', S(ii).Y([end 1])']);
    
%     plot(S(ii).X, S(ii).Y);
%     plot(xy(:,1), xy(:,2), '*');
    
    % Figure out which wall points need to be included, if any
    
    if t(1) < t(2)
        isin = twall > t(1) & twall < t(2);
        tnew = [t(1); twall(isin); t(2)];
        pt = interparc(tnew, xwall, ywall, 'linear');
    else
        isin1 = twall > t(1);
        isin2 = twall < t(2);
        tnew1 = [t(1); twall(isin1)];
        tnew2 = [twall(isin2); t(2)];
        pt = [interparc(tnew1, xwall, ywall, 'linear'); interparc(tnew2, xwall, ywall, 'linear')];
        
    end
    
%     plot(pt(:,1), pt(:,2), '--');
    
    S(ii).X = [S(ii).X(1:end-1) pt(:,1)'];
    S(ii).Y = [S(ii).Y(1:end-1) pt(:,2)'];
    
%     plot(S(ii).X, S(ii).Y, '--');
    
end

%--------------------
% Calculate overlap
% and triangulate
%--------------------

function [xnew, ynew] = removeoverlap(xc, yc, xwall, ywall, flag)

% Eliminate any empty contours or duplicates
        
[xwall, ywall] = poly2cw(xwall, ywall);

xc = [xwall; xc(:)];
yc = [ywall; yc(:)];
isemp = cellfun(@isempty, xc);

% for ii = 1:length(xc)
%     for jj = 1:length(xc)
%         iseq(ii,jj) = isequal(xc{ii},xc{jj}) & isequal(yc{ii},yc{jj});
%     end
% end
% [~,iunq] = unique(iseq, 'rows');
% 
% iseq = cellfun(@isequal, xc, xc') & cellfun(@isequal, yc, yc');
% 
% xunq = unique(xc);
% yunq = unique(yc);
% [~,ix] = ismember(xc, xunq);
% [~,iy] = ismember(yc, yunq);


xc = xc(~isemp);
yc = yc(~isemp);

% Calculate overlap

[xnew, ynew] = multiplepolyint(xc,yc);


function [f,v,lev] = poly2faces(xnew, ynew, F, clev)

np = length(xnew);

% Triangulate

[f,v] = deal(cell(np,1));
for ip = 1:np
    [f{ip},v{ip}] = poly2fv(xnew{ip}, ynew{ip});
end
isemp = cellfun('isempty', f);
f = f(~isemp);
v = v(~isemp);
np = length(f);

% There's probably a more elegant way to figure out which color
% goes where, but this works for now 

lev = zeros(np,1);
for ii = 1:np
    vx = v{ii}(:,1);
    vy = v{ii}(:,2);
    fx = mean(vx(f{ii}),2);
    fy = mean(vy(f{ii}),2);

    tmp = F(fx,fy);
%             tmp = interp2(x,y,z,fx,fy);
    [ntmp, bin] = histc(tmp, [-Inf clev Inf]);
    if ~any(ntmp)
        lev(ii) = length(clev)+2;
    else
        [~,lev(ii)] = max(ntmp);
    end
end

%--------------------
% Merge faces that 
% share an edge
%--------------------

function [xout,yout] = mergefaces(x,y)

[x,y] = polyjoin(x,y);
[vxy, ~, vidx] = unique([x y], 'rows');
nvert = sum(~isnan(vxy(:,1)));
isn = all(isnan(vxy),2);
vxy = vxy(1:nvert,:);
vidx(vidx > nvert) = NaN;

segs = [vidx(1:end-1) vidx(2:end)];
isn = any(isnan(segs),2);
segs = segs(~isn,:);

isdup = ismember(segs, fliplr(segs), 'rows');


xseg = [x(1:end-1) x(2:end)];
yseg = [y(1:end-1) y(2:end)];

isn = any(isnan(xseg),2);
xseg = xseg(~isn,:);
yseg = yseg(~isn,:);

xyunq = unique([xseg yseg], 'rows');


isdup = ismember([xseg yseg], [fliplr(xseg) fliplr(yseg)], 'rows');

%--------------------
% Voronoi regions for
% each point
%--------------------

function [xv,yv] = voronoigrid(x,y, xwall, ywall)

if isvector(x) && isvector(y)
    [xg, yg] = meshgrid(x,y);
else
    xg = x;
    yg = y;
end

[nr,nc] = size(xg);


xmid = filter2([1 1 ; 1 1]./4, padarray(xg, [1 1], 'replicate', 'both'), 'valid');
ymid = filter2([1 1 ; 1 1]./4, padarray(yg, [1 1], 'replicate', 'both'), 'valid');

[xv,yv] = deal(cell(size(xg)));
for ir = 1:nr
    for ic = 1:nc
        ridx = [0 1 1 0 0]+ir;
        cidx = [0 0 1 1 0]+ic;
        idx = sub2ind([nr+1 nc+1], ridx, cidx);
        xv{ir,ic} = xmid(idx);
        yv{ir,ic} = ymid(idx);
    end
end



