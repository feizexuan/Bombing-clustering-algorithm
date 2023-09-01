function [cm, lab, dist] = getmatrix(bg, varargin)
%GETMATRIX gets a connection matrix.
%
%   [MATRIX, ID, DISTANCES] = GETMATRIX(BIOGRAPH) converts the BIOGRAPH
%   object into a logical sparse matrix, where 1's indicate that a node
%   (row index) is connected to another node (column index). ID is a list
%   of the nodes 'ID' property and corresponds to the rows and columns of
%   MATRIX. DISTANCES is a column vector with one entry for every nonzero
%   entry in MATRIX traversed columnwise and representing the respective
%   'Weight' property for each edge. 
%
%   Example:
%
%      cm = [0 1 1 0 0;2 0 0 4 4;4 0 0 0 0;0 0 0 0 2;4 0 5 0 0];
%      bg = biograph(cm);
%      [cm,IDs,dist] = getmatrix(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID, BIOGRAPH.BIOGRAPH/GETNODESBYID,
%   BIOGRAPH.BIOGRAPH/VIEW, BIOGRAPH.NODE/GETANCESTORS,
%   BIOGRAPH.NODE/GETRELATIVES.

%   Copyright 2006-2016 The MathWorks, Inc.



 %%% check arguments
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:biograph:getmatrix:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'no_input_arguments'};
    for j = 1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        if(isstruct(pval))
             error(message('bioinfo:biograph:getmatrix:StructParamError'));
        end
        k = find(strncmpi(pname, okargs, length(pname)));
        if isempty(k)
            error(message('bioinfo:biograph:getmatrix:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:biograph:getmatrix:AmbiguousParameterName', pname));
        else
            switch(k)
            end
        end
    end
end

cm  = bg.to>0;
if nargout>1
    lab = get(bg.nodes, 'ID');
end
if nargout>2
    dist = get(bg.edges(nonzeros(bg.to)),'Weight');
    if numel(dist)>1
        dist = cell2mat(dist);
    end
end
