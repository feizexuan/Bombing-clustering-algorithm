function hgUpdate(h)
%HGUPDATE


% Copyright 2003-2010 The MathWorks, Inc.

if isempty(h.hgAxes)
    error(message('bioinfo:biograph:hgUpdate:nohgaxes'));
end

% get figure size
parent = get(h.hgAxes,'Parent');
fs = hgconvertunits(parent,get(parent,'position'),get(parent,'units'),'points',0);

% required points for graph
rs = h.BoundingBox*h.Scale + [ 0 0 8 8]; % gives extra 8 points for margins

% scale to be used for data
sc = min(fs(3:4)./rs(3:4));

% axes position (in points) (to cover all the figure)
as = [0 0 fs(3:4)];

% axes limits
al = [0 0 fs(3:4)/sc]-4; %compensates on each side for the 8 points for margins
al = al - (fs([3 4 3 4])/sc-rs([3 4 3 4]))/2;

setappdata(h.hgAxes,'Scale',sc)
set(h.hgAxes,'position',as,'XLim',al([1,3]),'YLim',al([2,4]),'Visible','on')

if ~isempty(get(h.hgAxes,'Children'))
    h.hgCorrectFontSize
end

% saves current view so 'Reset View' in pan and zoom mode will restore
% appropriately
resetplotview(h.hgAxes,'SaveCurrentView');



                  
            
