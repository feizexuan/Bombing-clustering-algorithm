function varargout =  maxflow(varargin)
%MAXFLOW calculates the maximum flow in a directed graph.
%
%  MAXFLOW method for BIOGRAPH objects extracts the adjacency matrix and
%  calls GRAPHMAXFLOW function. All other input and output arguments are
%  analogous to GRAPHMAXFLOW.
% 
%  See also graphmaxflow.

% Copyright 2006 The MathWorks, Inc.


if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[varargout{1:nargout}] = graphmaxflow(getweightmatrix(varargin{1}),varargin{2:end});



