function [xnew, ynew, indices] = multiplepolyint(x, y)
%MULTIPLEPOLYINT Multiple polygon intersection
%
% [xnew, ynew, indices] = multiplepolyint(x, y)
%
% Determines the regions where polygons overlap, distinguishing between
% each set of overlaps.
%
% Input variables:
%
%   x:          cell array of x vertices for each input polygon
%
%   y:          cell array of y vertices for each input polygon
%
% Output variables:
%
%   xnew:       cell array of x vertices for each resulting polygon
%
%   ynew:       cell array of y vertices for each resulting polygon
%
%   indices:    cell array holding indices of input polygons that
%               correspond to each output polygon.  For example, the output
%               polygon resulting from the overlap of polygons 1 and 2 will
%               has indices [1 2].

% Copyright 2006 Kelly Kearney

%------------------------------
% Check input
%------------------------------

if ~iscell(x) || ~iscell(y) || ~isequal(size(x), size(y))
    error('x and y must be cell arrays with the same dimensions');
end
if ~all(cellfun(@isvector, x)) || ~all(cellfun(@isvector, y))
    error('Contents of x and y must be vectors');
end
    
x = x(:);
y = y(:);

%------------------------------
% Find intersections
%------------------------------

xnew = x(1);
ynew = y(1);
indices = {1};

for ipoly = 2:length(x)
    
    x1 = x{ipoly};
    y1 = y{ipoly};
    
    for icomp = 1:length(xnew)
        
        x2 = xnew{icomp};
        y2 = ynew{icomp};
        
        % Intersecting 
        [xint{icomp}, yint{icomp}] = polybool('&', x1, y1, x2, y2);
        
        noint = isempty(xint{icomp});
        
        % Only in 2
        
        if noint
            xxor{icomp} = x2;
            yxor{icomp} = y2;
        else
            [xxor{icomp}, yxor{icomp}] = polybool('xor', x2, y2, xint{icomp}, yint{icomp});
        end
        
        noxor = isempty(xxor{icomp});
        
        
        % Indices
        
        if noint
            indint{icomp} = [];
        else
            indint{icomp} = [indices{icomp} ipoly];
        end
        
        if noxor
            indxor{icomp} = [];
        else
            indxor{icomp} = indices{icomp};
        end
        
    end
    
    % Only in 1
    
    [xallint, yallint] = polyjoin(xint, yint);
    if isempty(xallint)
        xout = x1;
        yout = y1;
    else
        [xout, yout] = polybool('xor', x1, y1, xallint, yallint);
    end
    
    if isempty(xout)
        indout = [];
    else
        indout = ipoly;
    end
    
    xtemp = [xint xxor {xout}];
    ytemp = [yint yxor {yout}];
    indtemp = [indint indxor {indout}];
    
    isbad = cellfun(@(a,b,c) isempty(a) & isempty(b) & isempty(c), xtemp, ytemp, indtemp);
    xnew = xtemp(~isbad);
    ynew = ytemp(~isbad);
    indices = indtemp(~isbad);
end
      
   