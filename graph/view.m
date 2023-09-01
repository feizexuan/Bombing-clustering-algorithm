function bh = view(h)
% VIEW draw the figure for a BIOGRAPH object
%
%    VIEW(BIOGRAPH) opens a new figure window and draws the graph
%    represented by the BIOGRAPH object. When the BIOGRAPH object is
%    already represented by an existing figure, this function only updates
%    the graph properties. 
%
%    BIOGRAPH_COPY = VIEW(BIOGRAPH) returns a handle to a deep copy of the
%    BIOGRAPH object and which is contained by the figure. When updating an
%    existing figure, you can use the returned handle to change object
%    properties programmatically or from the command line. The original
%    BIOGRAPH object is left unchanged. When the figure is closed the
%    handle is no longer valid.
%
%    Changing some properties (either programmatically or using the GUI)
%    updates automatically the graph. However, in some cases it is
%    necessary to use call again the layout engine (DOLAYOUT) to correct
%    the position and shape of the nodes and edges.
%
%   Examples:
%
%      % Create a BIOGRAPH object.
%      cm = [0 1 1 0 0;1 0 0 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 0 0];
%      bg = biograph(cm)
%
%      % Render it into a HG figure and get the figure handle.
%      h = view(bg)
%
%      % Change the color of all nodes and all edges.
%      set(h.Nodes,'Color',[.5 .7 1])
%      set(h.Edges,'LineColor',[0 0 0])
%
%      % Changing some properties requires to rerun the layout engine to
%      % correct the position and shape of the graph elements.
%      set(h.Nodes,'Shape','Diamond')
%      dolayout(h)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/DOLAYOUT,
%   BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID, BIOGRAPH.BIOGRAPH/GETNODESBYID,
%   BIOGRAPH.NODE/GETANCESTORS, BIOGRAPH.NODE/GETDESCENDANTS,
%   BIOGRAPH.NODE/GETRELATIVES, GET, SET.


% Copyright 2003-2004 The MathWorks, Inc.

if isempty(h.up) 
    % call constructor for gui
    g = biograph.bggui(h);
    if nargout
       bh = g.down;
    end
elseif isequal(class(h.up),'biograph.bggui')
    h.hgReDraw
    bh = h;
else
    error(message('bioinfo:biograph:view:nobggui'));
end

