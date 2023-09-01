function varargout =  allshortestpaths(varargin)
%ALLSHORTESTPATHS finds all the shortest paths in graph.
%
%  ALLSHORTESTPATHS method for BIOGRAPH objects extracts the adjacency
%  matrix and calls GRAPHALLSHORTESTPATHS function. All other input and
%  output arguments are analogous to GRAPHALLSHORTESTPATHS.
% 
%  See also graphallshortestpaths.

% Copyright 2006 The MathWorks, Inc.


[varargout{1:nargout}] = graphallshortestpaths(getweightmatrix(varargin{1}),varargin{2:end});



