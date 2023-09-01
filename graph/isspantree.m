function varargout =  isspantree(varargin)
%ISSPANTREE  test a spanning tree.
%
%  ISSPANTREE method for BIOGRAPH objects extracts the adjacency matrix and
%  calls GRAPHISSPANTREE function.
%
%  Remark: ISSPANTREE ignores the direction of the edges in the BIOGRAPH
%  object.
% 
%  See also graphisspantree.

% Copyright 2006 The MathWorks, Inc.


m = getmatrix(varargin{1});
m = max(m,m');
[varargout{1:nargout}] = graphisspantree(m,varargin{2:end});
