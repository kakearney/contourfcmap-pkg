function [xnew, ynew, indices] = multiplepolyint(x, y, varargin)
%MULTIPLEPOLYINT Multiple polygon intersection
%
% [xnew, ynew, indices] = multiplepolyint(x, y)
% [xnew, ynew, indices] = multiplepolyint(x, y, p1, v1, ...)
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
% Optional input parameters (passed as parameter/value pairs)
%
%   method:     method to use for polygon boolean operations:
%               'polybool': use polybool function, default for pre-2017b
%                   (correponds to older fast = false option)
%               'gpcmex': call gpcmex directly.  Faster than polybool
%                   option in some older versions of Matlab (corresponds to
%                   older fast = true option)
%               'polyshape': use polyshape object methods, default for
%                   R2017b+
%
%   v2gpcpath:  path to your local copy of vectorsToGPC.m.  The default is 
%               [matlabroot]/toolbox/map/map/private/vectorsToGPC.m, where
%               [matlabroot] is the directory where the MATLAB software is
%               installed.
%
%   gpcmexpath: path to your local copy of the gpcmex mex file.  The
%               default is
%               [matlabroot]/toolbox/map/map/private/gpcmex.[mex], where
%               [matlabroot] is the directory where the MATLAB software is
%               installed, and [mex] is the appropriate mex file extension
%               for your operating system. 
%
%   vfgpcpath:  path to your local copy of vectorsFromGPC.m.  The default
%               is [matlabroot]/toolbox/map/map/private/vectorsFromGPC.m,
%               where [matlabroot] is the directory where the MATLAB
%               software is installed.   
%
%   fast:       This is a mostly-deprecated option.  If no method is
%               specified and Matlab detects pre-2017b is running, this
%               logical scalar will switch between the polybool (false) and
%               gpcmex (true) methods.
%
%               Note: Can also be passed via the older syntax:
%                 [xnew,ynew,ind] = multiplepolyint(x,y,fastflag)
%               This option maintains back-compatibility with older
%               versions of this function. 
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

if nargin == 3 % For old syntax
    Opt.fast = varargin{1}; 
    Opt.v2gpcpath  = '';
    Opt.gpcmexpath = '';
    Opt.vfgpcpath  = '';
    validateattributes(Opt.fast, {'logical'}, {'scalar'});
else
    p = inputParser;

    p.addParameter('fast',       false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    p.addParameter('v2gpcpath',  '',    @(x) validateattributes(x, {'char'}, {}));
    p.addParameter('gpcmexpath', '',    @(x) validateattributes(x, {'char'}, {}));
    p.addParameter('vfgpcpath',  '',    @(x) validateattributes(x, {'char'}, {}));
    p.addParameter('method', '', @(x) validateattributes(x, {'char'}, {}));
    p.parse(varargin{:});

    Opt = p.Results;
end

if ~iscell(x) || ~iscell(y) || ~isequal(size(x), size(y))
    error('x and y must be cell arrays with the same dimensions');
end
if ~all(cellfun(@isvector, x)) || ~all(cellfun(@isvector, y))
    error('Contents of x and y must be vectors');
end

if isempty(Opt.method)
    if verLessThan('matlab', '9.3.0') % 2017b
        if Opt.fast
            Opt.method = 'gpcmex';
        else
            Opt.method = 'polybool';
        end
    else
        Opt.method = 'polyshape';
    end
end 
validatestring(Opt.method, {'gpcmex', 'polybool', 'polyshape'});
    
x = x(:);
y = y(:);

% If the fast flag is on, we need access to two functions hidden in a
% private directory of the mapping toolbox.

if strcmp(Opt.method, 'gpcmex')
    
    msg = ['Could not find copies of vectorsToGPC.m (%s) and/or ', ...
           'mex function gpcmex (%s).  Please verify these paths ', ...
           '(if you passed them as inputs) or provide the ', ...
           'appropriate paths (if the default is not ', ...
           'working). The proper files should be found in the ', ...
           'toolbox/map/map/private folder under the ', ...
           'directory where Matlab is installed'];

    mappath = fullfile(matlabroot, 'toolbox', 'map', 'map', 'private');
    if isempty(Opt.v2gpcpath)
        Opt.v2gpcpath = fullfile(mappath, 'vectorsToGPC.m');
    end
    if isempty(Opt.vfgpcpath)
        Opt.vfgpcpath = fullfile(mappath, 'vectorsFromGPC.m');
    end
    
    if isempty(Opt.gpcmexpath)
        Gpc = dir(fullfile(mappath, 'gpcmex*'));
        if length(Gpc) < 1
            error('multiplepolyint:gpcmex', 'Could not find gpcmex in default location (%s); please pass as input', mappath);
        end
        if length(Gpc) > 1
            warning('Mutiple gpcmex files found; using %s; to change, include path as input parameter', Gpc(1).name);
            Gpc = Gpc(1);
        end
        Opt.gpcmexpath = fullfile(mappath, Gpc.name);
            
    end
    
    if ~exist(Opt.vfgpcpath, 'file') || ~exist(Opt.gpcmexpath, 'file')
        error('multiplepolyint:privatepath', msg, Opt.v2gpcpath, Opt.gpcmexpath);
    end
    
    vectorsToGPC = function_handle(Opt.v2gpcpath);
    vectorsFromGPC = function_handle(Opt.vfgpcpath);
    gpcmex = function_handle(Opt.gpcmexpath);
end

%------------------------------
% Find intersections
%------------------------------

switch Opt.method
    case 'polybool' % The original way, using polybool

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

    case 'gpcmex' % The riskier way, going straight to gpcmex
    
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
    case 'polyshape'
        
        % Convert all polygons to GPC-style structures
        
        W = warning('off'); % turn off duplicate vertices warning
        try
            p = polyshape;
            for ii = length(x):-1:1
                p(ii) = polyshape(x{ii}, y{ii});
            end
            warning(W);
        catch ME
            warning(W);
            rethrow(ME);
        end

        % Start with first polygon

        pnew = p(1);
        indices = {1};

        isemp = @(p) isempty([p.Vertices]);

        pint = polyshape;
        pxor = polyshape;
        for ipoly = 2:length(x)

            p1 = p(ipoly);
            for icomp = 1:length(pnew)

                p2 = pnew(icomp);

                % Intersecting

                pint(icomp) = intersect(p1, p2);
                noint = isemp(pint(icomp));

                % Only in 2

                if noint
                    pxor(icomp) = p2;
                else
                    pxor(icomp) = xor(p2, pint(icomp));
                end

                noxor = isemp(pxor(icomp));

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

            
            pallint = pint(1); 
            for ii = 2:length(pint)
                pallint = union(pallint, pint(ii));
            end
            if isemp(pallint)
                pout = p1;
            else
                pout = xor(p1, pallint);
            end

            if isemp(pout)
                indout = [];
            else
                indout = ipoly;
            end

            ptemp = [pint pxor pout];
            indtemp = [indint indxor {indout}];

            isbad = cellfun(@(a,b) isempty(a.Vertices) & isempty(b), num2cell(ptemp), indtemp);
            pnew = ptemp(~isbad);
            indices = indtemp(~isbad);

        end

        % Convert back to cell arrays

        [xnew, ynew] = deal(cell(length(pnew),1));
        for ii = 1:length(pnew)
            [xnew{ii}, ynew{ii}] = boundary(pnew(ii));
        end
%         [xnew, ynew] = cellfun(vectorsFromGPC, pnew, 'uni', 0);
%         xnew = xnew(:);
%         ynew = ynew(:);
  
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



      
   