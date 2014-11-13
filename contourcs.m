function Cout = contourcs(varargin)
%CONTOURCS   Wrapper to CONTOURS to Obtain Structure Output
%   S = CONTOURCS(...) takes the exact same input arguments as the default
%   CONTOURC but its output is a struct array with fields:
%
%     Level  - contour line value
%     Length - number of contour line points
%     X      - X coordinate array of the contour line
%     Y      - Y coordinate array of the contour line
%
%   See also contourc.

% Version 1.0 (Aug 11, 2010)
% Written by: Takeshi Ikuma
% Revision History:
%  - (Aug. 11, 2010) : initial release

% Run CONTOURC and get output matrix
try
   C = contourc(varargin{:});
catch ME
   throwAsCaller(ME);
end

% Count number of contour segments found (K)
K = 0;
n0 = 1;
while n0<=size(C,2)
   K = K + 1;
   n0 = n0 + C(2,n0) + 1;
end

% initialize output struct
el = cell(K,1);
Cout = struct('Level',el,'Length',el,'X',el,'Y',el);

% fill the output struct
n0 = 1;
for k = 1:K
   Cout(k).Level = C(1,n0);
   idx = (n0+1):(n0+C(2,n0));
   Cout(k).Length = C(2,n0);
   Cout(k).X = C(1,idx);
   Cout(k).Y = C(2,idx);
   n0 = idx(end) + 1; % next starting index
end

% Copyright (c)2010, Takeshi Ikuma
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   * Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer. 
%   * Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%   documentation and/or other materials provided with the distribution.
%   * Neither the names of its contributors may be used to endorse or
%   promote products derived from this software without specific prior
%   written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
