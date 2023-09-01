function varargout =  minspantree(varargin)
%MINSPANTREE finds the minimal spanning tree in graph.
%
%  MINSPANTREE method for BIOGRAPH objects extracts the adjacency matrix
%  and calls GRAPHMINSPANTREE function. All other input and output
%  arguments are analogous to GRAPHMINSPANTREE.
%
%  Remark: MINSPANTREE ignores the direction of the edges in the BIOGRAPH
%  object.
% 
%  See also graphminspantree.

% Copyright 2006 The MathWorks, Inc.


if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

m = getweightmatrix(varargin{1});   
m = max(m,m');
[varargout{1:nargout}] = graphminspantree(m,varargin{2:end});








