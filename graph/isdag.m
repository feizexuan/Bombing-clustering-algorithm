function varargout =  isdag(varargin)
%ISDAG tests for cycles in a directed graph.
%
%  ISDAG method for BIOGRAPH objects extracts the adjacency matrix and
%  calls GRAPHISDAG function.
% 
%  See also graphisdag.

% Copyright 2006 The MathWorks, Inc.


[varargout{1:nargout}] = graphisdag(getmatrix(varargin{1}),varargin{2:end});
