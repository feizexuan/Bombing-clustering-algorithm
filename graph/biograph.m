function h = biograph(cm,ids,pvp)
%BIOGRAPH class constructor.
%
%  h = biograph.biograph(cm,ids)

warning(message("bioinfo:biograph:biograph:ClassToBeRemoved"));

% Copyright 2003-2016 The MathWorks, Inc.

% object constructor
if nargin > 1
    if isstring(ids)
        ids = cellstr(ids);
    end
end

if nargin > 2
    if isstring(pvp)
        pvp = cellstr(pvp);
    end
end

h = biograph.biograph;

% add nodes
numNodes = size(cm,1);
for i = numNodes:-1:1
    nodes(i) = biograph.node(h,ids{i},i); 
end

% save handles to nodes in the biograph object
h.nodes = nodes;

% remove self_loops (not valid in a biograph object)
cm = cm - diag(diag(cm));

% find all nonzeros in the conexion matrix
[kk,jj] = find(cm'); % transpose to make it backwards compatible

% add edges
numEdges = numel(jj);
for i = numEdges:-1:1
     j = jj(i); k = kk(i);
     str = [h.Nodes(j).ID ' -> ' h.Nodes(k).ID];
     edges(i) = biograph.edge(h,h.Nodes(j),h.Nodes(k),str,cm(j,k),i);
end

% save handles to edges in the biograph object
if numEdges>0
    h.edges = edges;
else
    h.edges = [];
end

% save sparse pointers for efficient access
h.to   = sparse(jj,kk,1:numEdges,numNodes,numNodes);
h.from = sparse(kk,jj,1:numEdges,numNodes,numNodes);

% modify some object properties (undocumented)
npvp = numel(pvp);
if npvp
    for j = 1:2:npvp-1
        h.(pvp{j}) = pvp{j+1};
    end
end

