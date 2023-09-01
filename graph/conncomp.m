function varargout =  conncomp(varargin)
%CONNCOMP finds the connected components in graph.
%
%  CONNCOMP method for BIOGRAPH objects extracts the adjacency matrix and
%  calls GRAPHCONNCOMP function. All other input and output arguments are
%  analogous to GRAPHCONNCOMP.
% 
%  See also graphconncomp.

% Copyright 2006 The MathWorks, Inc.


if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[varargout{1:nargout}] = graphconncomp(getmatrix(varargin{1}),varargin{2:end});
