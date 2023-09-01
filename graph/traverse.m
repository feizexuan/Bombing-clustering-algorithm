function varargout =  traverse(varargin)
%TRAVERSE performs traversal of the graph by following adjacent nodes.
%
%  TRAVERSE method for BIOGRAPH objects extracts the adjacency matrix and
%  calls GRAPHTRAVERSE function. All other input and output arguments are
%  analogous to GRAPHTRAVERSE.
% 
%  See also graphtraverse.

% Copyright 2006 The MathWorks, Inc.


[varargout{1:nargout}] = graphtraverse(getmatrix(varargin{1}),varargin{2:end});
