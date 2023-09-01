function bh = deepCopy(bho)
% DEEPCOPY 

%   Copyright 2003-2006 The MathWorks, Inc.


bh = bho.copy;
for k = 1:numel(bho.Nodes)
     bh.Nodes(k) = bho.Nodes(k).copy;
     bh.Nodes(k).connect(bh,'up')
end
[kk,jj] = find(bh.from);
for i = 1:numel(bho.Edges)
    bh.Edges(i) = bho.Edges(i).copy;
    bh.Edges(i).FromNode = bh.Nodes(jj(i));
    bh.Edges(i).ToNode   = bh.Nodes(kk(i));
    bh.Edges(i).connect(bh,'up')
end
