function hgReDraw(h,e)
%HGREDRAW


% Copyright 2003-2006 The MathWorks, Inc.

if isempty(h.hgAxes)
    error(message('bioinfo:biograph:hgUpdate:nohgaxes'));
end

if nargin>1 % special case when called with two arguments then it means 
            % that a property has been updated, then we only redraw what is
            % needed...
    switch e.Source.Name
        case {'ShowArrows','ShowWeights','EdgeType','ArrowSize','EdgeFontSize','EdgeTextColor'}
            h.edges.hgUpdate;
        case {'Scale'}
            h.nodes.hgUpdate;
            h.edges.hgUpdate;
            h.hgUpdate
        case {'ShowTextInNodes','CustomNodeDrawFcn'}
            h.nodes.hgUpdate;
        otherwise
            % do nothing
    end
else
    h.nodes.hgUpdate;
    h.edges.hgUpdate;
end

