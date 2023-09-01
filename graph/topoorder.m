function varargout =  topoorder(varargin)
%TOPOORDER performs topological sort of a directed acyclic graph.
%
%  TOPOORDER method for BIOGRAPH objects extracts the adjacency matrix and
%  calls GRAPHTOPOORDER function.
% 
%  See also graphtopoorder.

% Copyright 2006 The MathWorks, Inc.


[varargout{1:nargout}] = graphtopoorder(getmatrix(varargin{1}),varargin{2:end});
