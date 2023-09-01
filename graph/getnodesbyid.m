function nodes = getnodesbyid(h,ids)
%GETNODESBYID gets handles to nodes.
%
%   NODES = GETNODESBYID(BIOGRAPH,NODE_IDS) gets the node handles for the
%   queried ID's. NODE_IDS can be a string or cell string.
%
%   Example:
%
%      % Creates a BIOGRAPH object
%      species = {'Homo','Pan','Gorilla','Pongo','Baboon','Macaca','Gibbon'};
%      cm = magic(7)>25 & ~eye(7);
%      bg = biograph(cm,species)
%
%      % Find the handles to the members of the Cercopithecidae family and
%      % members of the Hominidae family and colors them differently
%      Cercopithecidae = {'Macaca','Baboon'};
%      Hominidae = {'Homo','Pan','Gorilla','Pongo'};
%      CercopithecidaeNodes = getnodesbyid(bg,Cercopithecidae);
%      HominidaeNodes = getnodesbyid(bg,Hominidae);
%      set(CercopithecidaeNodes,'Color',[.7 1 .7])
%      set(HominidaeNodes,'Color',[1 .7 .7])
%      view(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID, BIOGRAPH.BIOGRAPH/VIEW,
%   BIOGRAPH.NODE/GETANCESTORS, BIOGRAPH.NODE/GETDESCENDANTS,
%   BIOGRAPH.NODE/GETRELATIVES, GET, SET.

% Copyright 2003-2004 The MathWorks, Inc.


switch nargin
    case 1 % returns all nodes (undocumented)
        nodes = h.Nodes;
    otherwise 
        if isempty(ids)
           nodes = handle(zeros(0,1));
        else
            switch numel(h.Nodes)
                case 0
                    nodes = handle(zeros(0,1));
                case 1
                    nodes = h.Nodes(isequal(h.Nodes.ID,ids));
                otherwise
                    [tf,loc] = ismember(ids,get(h.Nodes,'ID'));
                    nodes = h.Nodes(nonzeros(loc));
            end
        end
end

