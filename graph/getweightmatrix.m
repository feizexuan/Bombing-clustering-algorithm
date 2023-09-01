function [cm, lab] = getweightmatrix(bg)
%GETWEIGHTMATRIX gets a connection matrix with weights.
%
%   [MATRIX, ID] = GETWEIGHTMATRIX(BIOGRAPH) converts the BIOGRAPH object
%   into a double sparse matrix, where non-zeros entries indicate the
%   weight from the source node (row index) to the destination node (column
%   index). ID is a list of the nodes 'ID' property and corresponds to the
%   rows and columns of MATRIX.
%
%   Example:
%
%      cm = [0 1 1 0 0;2 0 0 4 4;4 0 0 0 0;0 0 0 0 2;4 0 5 0 0];
%      bg = biograph(cm);
%      [cm,IDs] = getweightmatrix(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID, BIOGRAPH.BIOGRAPH/GETMATRIX
%   BIOGRAPH.BIOGRAPH/GETNODESBYID, BIOGRAPH.BIOGRAPH/VIEW,
%   BIOGRAPH.NODE/GETANCESTORS, BIOGRAPH.NODE/GETRELATIVES.

%   Copyright 2006-2016 The MathWorks, Inc.


w = cell2mat(get(bg.edges,'Weight'));
[i,j] = find(bg.from);
n = numel(bg.Nodes);
cm = sparse(j,i,w,n,n);  

if nargout>1
    lab = get(bg.nodes,'ID');
end
