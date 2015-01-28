function [v, A] = fillnan(v, coords, varargin)
%FILLNAN Fill in missing values in an array with nearest neighbor
%
% [v, A] = fillnan(v, coords, p1, v1, ...)
%
% When upsampling geographic data, I often need to first perform a
% nearest-neighbor interpolation to fill in land-masked data cells and
% prevent replication of these NaNs into now-resolved regions.  This
% function is much faster than inpainting via a scattered interpolant.
%
% Input variables:
%
%   v:      array of values, 2D or 3D array, nd-grid style
%
%   coords: 1 x 2 or 1 x 3 array of coordinates
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
%           the edges of a dataset, without wasting computation time on
%           inland grid cells. [false] 
%
% Output variables:
%
%   v:      array matching input v, with all NaNs (or NaNs bordering data)
%           now filled in. 
%
%   A:      1 x 1 structure with the following fields.  Can be used to
%           repeat interpolation without recalculating nearest neighbors.
%
%           idxfill:    indices to the missing-data cells that are
%                       filled in 
%
%           idxnn:      indices to the nearest-neighbor data cells used to
%                       fill in each of the cells in the above array,
%                       respectively.

% Copyright 2013 Kelly Kearney

%------------------------
% Parse input
%------------------------

Opt.mask = isnan(v);
Opt.geo  = false;
Opt.perim = false;

Opt = parsepv(Opt, varargin);

[nx, ny, nz] = size(v);

% Convert coordinates if necessary

if cellfun(@isvector, coords)
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








