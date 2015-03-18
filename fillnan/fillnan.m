function [v, A] = fillnan(v, coords, varargin)
%FILLNAN Fill in missing values in an array with nearest neighbor
%
% [v, A] = fillnan(v)
% [v, A] = fillnan(v, coords, p1, v1, ...)
%
% This function fills NaN (or otherwise masked) points in a gridded array
% with the value from the nearest non-masked grid point.  It is much faster
% than using a scattered interpolant for this, or inpainting functions like
% inpaint_nans (which perform more complex fitting calculations).  The
% output also allows for resuse of the nearest-neighbor calculation on
% different datasets that share a coordinate system.
%
% Input variables:
%
%   v:      array of values, 2D or 3D array, nd-grid style
%
%   coords: 1 x 2 or 1 x 3 cell array of grid coordinates, corresponding to
%           x, y, and z coordinates, respectively.  Cells should either be
%           the same size as v, or vectors following ndgrid conventions.
%           If empty or not included, coordinates will correspond to array
%           indices. 
%
% Optional input variables (passed as parameter/value pairs):
%
%   mask:   mask indicating where data is missing (1) and present (0).  If
%           not included, defaults to isnan(v)
%
%   geo:    true if input coordinates represent geographic coordinates
%           {longitude, latitude, height (m)}, rather than cartesian
%           {x,y,z} coordinates. [false]
%
%   perim:  true to fill only values along the edge of data/no-data
%           boundaries.  This option is useful when preparing to upsample a
%           dataset when you only need to prevent propagation of NaNs along
%           the edges of a mask (for example, upsampling ocean data where
%           land is masked), without wasting computation time on filling
%           inland grid cells. [false]     
%
% Output variables:
%
%   v:      array matching input v, with all NaNs (or NaNs bordering data)
%           now filled in. 
%
%   A:      1 x 1 structure with the following fields.  Can be used to
%           repeat interpolation without recalculating nearest neighbors,
%           e.g. v2(A.idxnn) = v2(A.idxfill).
%
%           idxfill:    indices to the missing-data cells that are
%                       filled in 
%
%           idxnn:      indices to the nearest-neighbor data cells used to
%                       fill in each of the cells in the above array,
%                       respectively.

% Copyright 2013-2015 Kelly Kearney

%------------------------
% Parse input
%------------------------

Opt.mask = isnan(v);
Opt.geo  = false;
Opt.perim = false;

Opt = parsepv(Opt, varargin);

if ndims(v) > 3
    error('Only 2D and 3D arrays supported at this time');
end
[nx, ny, nz] = size(v);

% Convert coordinates if necessary

if nargin < 2 || isempty(coords)
    coords = {1:nx, 1:ny, 1:nz};
end

if cellfun(@isvector, coords) % all implicit
    [coords{:}] = ndgrid(coords{:});
end

if Opt.geo

    E = referenceEllipsoid('earth'); % WGS84, m

    phi    = degtorad(coords{1});
    lambda = degtorad(coords{2});
    if length(coords) == 3
        h = coords{3};
    else
        h = zeros(size(phi));
    end
    [x,y,z] = geodetic2ecef(phi, lambda, h, E);
else
    x = coords{1};
    y = coords{2};
    if length(coords) == 3
        z = coords{3};
    else
        z = zeros(size(x));
    end
end

%------------------------
% Match no-data points to
% data points
%------------------------

if Opt.perim

    % Find perimeter values

    mask = padarray(Opt.mask, [1 1 1], 'both', 'replicate'); % Pad so box isn't considered boundary
    bwp = bwperim(mask);
    bwp = bwp(2:end-1, 2:end-1, 2:end-1);

    % Dilate perimeter in all directions

    bwd = imdilate(bwp, ones(3,3,3));

    % Match data-perimeter values to missing-perimeter values

    is1 = bwd &  Opt.mask; % perimeter w/o data
    is2 = bwd & ~Opt.mask; % perimeter w/ data
    
else
    
    is1 =  Opt.mask;
    is2 = ~Opt.mask;
    
end

idx = knnsearch([x(is2) y(is2) z(is2)], [x(is1) y(is1) z(is1)]);

A.idxfill = find(is1);
tmp = find(is2);
A.idxnn = tmp(idx);

% Fill with nearest neighbor

v(A.idxfill) = v(A.idxnn);








