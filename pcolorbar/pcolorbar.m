function h = pcolorbar(clev, cmap, varargin)
%PCOLORBAR Create a pseudo-colorbar with static color blocks
%
% h = pcolorbar(clev, cmap, ...)
%
% This function creates an axis that mimics a colorbar.  However, rather
% than using data from the axis' actual colormap, it plots the specific
% color intervals and colors specified by the user.
%
% Input variables:
%
%   clev:       vector of length n defining the edges of the color blocks
%
%   cmap:       n-1 x 3 colormap array, colors corresponding to the
%               intervals between each level in the clev array 
%
% Optional input variables (passed as parameter/value pairs):
%
%   peer:       axis to attach colorbar to [gca]
%
%   location:   location of colorbar (same options as regular colorbar
%               command) ['east']
%
%   even:       logical scalar, if true, create evenly-sized color blocks
%               regardless of the values assigned.  Ticks labels will
%               reflect the values in clev. [false]
%
%   tkint:      interval used to label ticks, with tick 1:tkint:end
%               receiving a label.
%
%   cmapt:      n-1 x 3 colormap array corresponding to the top/right of
%               each color block.  If a row in cmapt is different than the
%               corresponding row in cmap, then the block will show a
%               linear gradient between the two colors.  By default,
%               solid-colored blocks are drawn.
%
% Output variables:
%
%   h:          structure with the following fields:
%
%               cb: handle to the non-visible, real colorbar.  The
%                   pseudocolorbar will shadow this one for position.
%                   (Note: if you want to change the position of the
%                   pseudocolorbar, change this handle's position
%                   property).
%
%               ax: handle to pseudocolorbar axis
%
%               p:  handle to pseudocolorbar patch object
%
%               el: event listener that links h.ax position to h.cb
%                   position.

% Copyright 2018 Kelly Kearney

% Parse input

p = inputParser;
p.addParameter('peer', gca, @(x) validateattributes(x, {}, {'scalar'}));
p.addParameter('location', 'east', @(x) validateattributes(x, {'char','string'}, {'scalartext'}));
p.addParameter('even', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('tkint', 1, @(x) validateattributes(x, {'numeric'}, {'integer'}));
p.addParameter('cmapt', zeros(0,3), @(x) validateattributes(x, {'numeric'}, {'ncols',3,'>=',0,'<=',1}));

p.parse(varargin{:});

Opt = p.Results;

validateattributes(cmap, {'numeric'}, {'ncols',3,'>=',0,'<=',1}, 'pcolorbar', 'cmap');
validateattributes(clev, {'numeric'}, {'vector', 'increasing'}, 'pcolorbar', 'clev');

nblock = length(clev) - 1;

if isempty(Opt.cmapt)
    Opt.cmapt = cmap;
end

if ~isequal(nblock, size(cmap,1), size(Opt.cmapt,1))
    error('Colormap array(s) should include one less color than clev array');
end

clev = clev(:);

% Create a real colorbar, used for position info

h.cb = colorbar('peer', Opt.peer, 'location', Opt.location);

loc = get(h.cb, 'Location'); % Easier than lower/upper
pos = get(h.cb, 'position');

% Calculate y-values for each color interval, and labels to go along

if Opt.even
    tk = linspace(0,1,nblock+1);
    y1 = tk(1:end-1);
    y2 = tk(2:end);
else
    y1 = clev(1:end-1);
    y2 = clev(2:end);
    tk = clev;
end
tklbl = strtrim(cellstr(num2str(clev)));

islbl = ismember(1:length(tk), 1:Opt.tkint:length(tk));
[tklbl{~islbl}] = deal(' ');

% Coordinates and values for patch vertices

ypatch = [y1 y2 y2 y1 y1]';
xpatch = repmat([0 0 1 1 0], size(ypatch,2), 1)';

cpatch = cat(3, cmap, Opt.cmapt, Opt.cmapt, cmap, cmap);
cpatch = permute(cpatch, [3 1 2]);

% Create pseudo-colorbar

currax = get(gcf, 'CurrentAxes');
h.ax = axes('position', pos, 'box', 'on');

switch lower(loc)
    case {'east', 'west', 'eastoutside', 'westoutside'}
        h.p = patch(xpatch, ypatch, cpatch);
        set(h.ax, 'ytick', tk, 'yticklabel', tklbl, 'ylim', minmax(tk), ...
            'xlim', [0 1], 'xtick', []);
    otherwise
        h.p = patch(ypatch, xpatch, cpatch);   
        set(h.ax, 'xtick', tk, 'xticklabel', tklbl, 'xlim', minmax(tk), ...
            'ylim', [0 1], 'ytick', []);
end

switch lower(loc)
    case {'eastoutside', 'west'}
        set(h.ax, 'yaxislocation', 'right');
    case {'south', 'northoutside'}
        set(h.ax, 'xaxislocation', 'top');
end

set(h.ax, 'layer', 'top');

set(h.p, 'edgecolor', 'none');
set(h.cb, 'visible', 'off');

% Create listener so axis size updates with colorbar

try
    h.el = addlistener(h.cb, 'MarkedClean', @(he,ed) set(h.ax, 'position', get(he, 'position')));
catch
    warning('Unable to add colorbar resizing listener');
end