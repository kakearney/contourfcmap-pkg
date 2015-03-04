function [xnew, ynew, indices] = multiplepolyint(x, y, flag)
%MULTIPLEPOLYINT Multiple polygon intersection
%
% [xnew, ynew, indices] = multiplepolyint(x, y)
% [xnew, ynew, indices] = multiplepolyint(x, y, flag)
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
%   flag:       logical scalar, indicating whether to use fast method
%               (true) or not (false).  Default is false.  The fast method
%               basically skips over polybool and uses gpcmex directly;
%               because it accesses private functions in the Mapping
%               Toolbox, it may be more fragile than the default method.
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
%               have indices [1 2].

% Copyright 2006-2015 Kelly Kearney

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

if nargin < 3
    flag = false;
end
if ~(isscalar(flag) && islogical(flag))
    error('flag should be a logical value');
end

% If the fast flag is on, we need access to two functions hidden in a
% private directory of the mapping toolbox.

if flag
    v2gpcpath = fullfile(matlabroot, 'toolbox', 'map', 'map', 'private','vectorsToGPC.m');
    vfgpcpath = fullfile(matlabroot, 'toolbox', 'map', 'map', 'private','vectorsFromGPC.m');
    gpcmexpath = fullfile(matlabroot, 'toolbox', 'map', 'map', 'private','gpcmex.mexmaci64');
    if ~exist(vfgpcpath, 'file') || ~exist(gpcmexpath, 'file')
        error('multiplepolyint:privatepath', ...
            ['Please modify the paths in multiplepolyint.m (above this) to point to\n', ...
             'your copies of vectorsToGPC.m and the mex function gpcmex.  These can\n', ...
             'be found in the toolbox/map/map/private folder of the Mapping Toolbox']);
    end
    vectorsToGPC = function_handle(v2gpcpath);
    vectorsFromGPC = function_handle(vfgpcpath);
    gpcmex = function_handle(gpcmexpath);
end

%------------------------------
% Find intersections
%------------------------------

if ~flag % The original way, using polybool

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

else % The riskier way, going straight to gpcmex
    
    % Convert all polygons to GPC-style structures
    
    p = cell(length(x),1);
    for ii = 1:length(x)
        [xtmp, ytmp] = polysplit(x{ii},y{ii});
        p{ii} = struct('x', xtmp, 'y', ytmp, 'ishole', num2cell(~ispolycw(xtmp, ytmp)));
    end
    
    % Start with first polygon
    
    pnew = p(1);
    indices = {1};
    
    isemp = @(p) isempty([p.x]);
    
    for ipoly = 2:length(x)
        
        p1 = p{ipoly};
        for icomp = 1:length(pnew)
            
            p2 = pnew{icomp};
            
            % Intersecting
            
            pint{icomp} = gpcmex('int', p1, p2);
            noint = isemp(pint{icomp});
            
            % Only in 2
            
            if noint
                pxor{icomp} = p2;
            else
                pxor{icomp} = gpcmex('xor', p2, pint{icomp});
            end

            noxor = isemp(pxor{icomp});
            
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
        
        pallint = cat(1, pint{:});
        if isemp(pallint)
            pout = p1;
        else
            pout = gpcmex('xor', p1, pallint);
        end
        
        if isemp(pout)
            indout = [];
        else
            indout = ipoly;
        end
        
        ptemp = [pint pxor {pout}];
        indtemp = [indint indxor {indout}];
        
        isbad = cellfun(@(a,b) isempty(a) & isempty(b), ptemp, indtemp);
        pnew = ptemp(~isbad);
        indices = indtemp(~isbad);
        
    end
    
    % Convert back to cell arrays
    
    [xnew, ynew] = cellfun(vectorsFromGPC, pnew, 'uni', 0);
    xnew = xnew(:);
    ynew = ynew(:);
  
end
    

%------------------------------
% Subfunction: stripped down
% polybool
%------------------------------

% Note: This speeds things up just a bit, but much more speed could be
% gained by skipping vectors[To/From]GPC entirely, because it's the
% clockwise checks in there that eat up a lot of unecessary time.

function [x3, y3] = polyboolfast(op, x1, y1, x2, y2, vectorsToGPC, vectorsFromGPC, gpcmex)

% Handle empties

if  all(isnan(x1)) && ~isempty(x1)
    x1(1:end) = [];
    y1(1:end) = [];
end

if all(isnan(x2)) && ~isempty(x2)
    x2(1:end) = [];
    y2(1:end) = [];
end

if isempty(x2)
    if strcmp(op,'int')
        % Intersection is empty, but preserve shape
        % by using x2 and y2 rather than [].
        x3 = x2;
        y3 = y2;
    else
        % Union, exclusive or, or difference with
        % empty leaves x1 and y1 unaltered.
        x3 = x1;
        y3 = y1;
    end
    emptyInput = true;
elseif isempty(x1)
    if any(strcmp(op,{'int','diff'}))
        % Intersection or difference is empty, but preserve
        % shape by using x1 and y1 rather than [].
        x3 = x1;
        y3 = y1;        
    else
        % Union or exclusive or with empty leaves x2 and y2 unaltered.
        x3 = x2;
        y3 = y2;
    end
    emptyInput = true;
else
    x3 = [];
    y3 = [];
    emptyInput = false;
end

% Calculate

if ~emptyInput
    p1 = vectorsToGPC(x1, y1);
    p2 = vectorsToGPC(x2, y2);
    if (length(p1)==1 && p1.ishole) || (length(p1)>1 && (p1(1).ishole || ~all([p1(2:end).ishole]))) || ...
       (length(p2)==1 && p2.ishole) || (length(p1)>1 && (p2(1).ishole || ~all([p1(2:end).ishole])))     
        blah
    end
    p3 = gpcmex(op, p1, p2);
    [x3, y3] = vectorsFromGPC(p3);
end



      
   