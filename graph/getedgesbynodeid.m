function edges = getedgesbynodeid(h,from,to)
%GETEDGESBYNODEID gets handles to edges
%
%   EDGES = GETEDGESBYNODEID(BIOGRAPH,SOURCE_IDS,SINK_IDS) gets the edge
%   handles that connect any of the nodes' SOURCE_IDS to any of the nodes'
%   SINK_IDS. SOURCE_IDS and SINK_IDS can be either strings, cell strings,
%   or empty strings (gets all).
%
%   EDGES = GETEDGESBYNODEID(BIOGRAPH,NODE_IDS) gets the edge handles that
%   connect the nodes of the first column of NODE_IDS to the nodes in the
%   second column. NODE_IDS is a cell array of strings with two columns and
%   any number of row.
%
%   EDGES = GETEDGESBYNODEID(BIOGRAPH,PATH_IDS) get the edge handles that
%   connect contiguous nodes in a path. PATH_IDS is a cell string.
%
%   Examples:
%
%      % Create a BIOGRAPH object.
%      species = {'Homo','Pan','Gorilla','Pongo','Baboon','Macaca','Gibbon'};
%      cm = magic(7)>25 & ~eye(7);
%      bg = biograph(cm,species)
%
%      % Find all the edges that connect to the homo sapiens node.
%      homoEdgesIn = getedgesbynodeid(bg,[],'Homo');
%      homoEdgesOut = getedgesbynodeid(bg,'Homo',[]);
%      set(homoEdgesIn,'LineColor',[0 1 0])
%      set(homoEdgesOut,'LineColor',[1 0 0])
%      view(bg)
%
%      % Find all edges that connect members of the Cercopithecidae family
%      % to members of the Hominidae family.
%      Cercopithecidae = {'Macaca','Baboon'};
%      Hominidae = {'Homo','Pan','Gorilla','Pongo'};
%      edgesSel = getedgesbynodeid(bg,Cercopithecidae,Hominidae);
%      set(bg.edges,'LineColor',[.5 .5 .5])
%      set(edgesSel,'LineColor',[0 0 1])
%      view(bg)
%
%      % Find edges between specific pairs of nodes.
%      GroupA = {'Macaca','Baboon','Homo','Homo'};
%      GroupB = {'Gorilla','Pongo','Pan','Gorilla'};
%      edgesSel = getedgesbynodeid(bg,[GroupA;GroupB]');
%      set(bg.edges,'LineColor',[.5 .5 .5])
%      set(edgesSel,'LineColor',[0 0 1])
%      view(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETNODESBYID, BIOGRAPH.BIOGRAPH/VIEW,
%   BIOGRAPH.NODE/GETANCESTORS, BIOGRAPH.NODE/GETDESCENDANTS,
%   BIOGRAPH.NODE/GETRELATIVES, GET, SET.


% Copyright 2003-2016 The MathWorks, Inc.


switch nargin
    case 1 % returns all edges
        edges = h.Edges;
    case 2 % variable 'from' is either PATH_IDS or NODE_IDS
        if ~iscellstr(from)
            error(message('bioinfo:biograph:getedgesbynodeid:NoCellStr'));
        end
        if size(from,2) == 2 % NODE_IDs with two columns
            to = from(:,2);
            from = from(:,1);
        elseif numel(from) == 1 % PATH_IDs with a single node
            from = from(:);
            to = from;
        else % PATH_IDS
            from = from(:);
            to = from(2:end);
            from = from(1:end-1);
        end
        if isempty(from)
            fromNodes = mycell2mat(get(h.getnodesbyid,'idx'));
        else
            fromNodes = mycell2mat(get(h.getnodesbyid(from),'idx'));
        end
        if isempty(to)
            toNodes = mycell2mat(get(h.getnodesbyid,'idx'));
        else
            toNodes = mycell2mat(get(h.getnodesbyid(to),'idx'));
        end
        if numel(fromNodes) ~= numel(toNodes) || numel(fromNodes) ~= size(from,1) || numel(toNodes) ~= size(to,1)
            error(message('bioinfo:biograph:getedgesbynodeid:NoValidNodeID'));
        end
        idxs = sub2ind(repmat(numel(h.Nodes),1,2),fromNodes,toNodes);
        edges = h.edges(nonzeros(h.to(idxs)));
    case 3
        if isempty(from)
            fromNodes = mycell2mat(get(h.getnodesbyid,'idx'));
        else
            fromNodes = mycell2mat(get(h.getnodesbyid(from),'idx'));
            if (iscell(from) && numel(fromNodes) ~= numel(from)) || ...
                    (ischar(from) && numel(fromNodes) ~= size(from,1))
                error(message('bioinfo:biograph:getedgesbynodeid:NoValidNodeID'));
            end
        end
        if isempty(to)
            toNodes = mycell2mat(get(h.getnodesbyid,'idx'));
        else
            toNodes = mycell2mat(get(h.getnodesbyid(to),'idx'));
             if (iscell(to) && numel(toNodes) ~= numel(to)) || ...
                    (ischar(to) && numel(toNodes) ~= size(to,1))
                error(message('bioinfo:biograph:getedgesbynodeid:NoValidNodeID'));
             end
        end
        edges = h.edges(nonzeros(h.to(fromNodes,toNodes)));        
    otherwise
        error(message('bioinfo:biograph:getedgesbynodeid:IncorrectNumberOfArguments', mfilename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = mycell2mat(c)
if numel(c)==1
    m = c;
else
    m = cell2mat(c);
end
    
