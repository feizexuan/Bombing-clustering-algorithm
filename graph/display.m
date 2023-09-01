function display(h)
%DISPLAY command window display.

% Copyright 2003-2004 The MathWorks, Inc.


disp(sprintf('Biograph object with %d nodes and %d edges.',...
     numel(h.Nodes),numel(h.Edges)))
