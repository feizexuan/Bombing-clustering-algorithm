function varargout =  shortestpath(varargin)
%SHORTESTPATH solves the shortest path problem in graph.
%
%  SHORTESTPATH method for BIOGRAPH objects extracts the adjacency matrix
%  and calls GRAPHSHORTESTPATH function. All other input and output
%  arguments are analogous to GRAPHSHORTESTPATH.
% 
%  See also graphshortestpath.

% Copyright 2006 The MathWorks, Inc.


if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[varargout{1:nargout}] = graphshortestpath(getweightmatrix(varargin{1}),varargin{2:end});
