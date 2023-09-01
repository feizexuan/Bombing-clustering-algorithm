function hgCorrectFontSize(h)
%HGCORRECTFONTSIZE

% Copyright 2003-2015 The MathWorks, Inc.

minAllowableFontSize = 1.5;

if isempty(h.hgAxes)
    error(message('bioinfo:biograph:hgUpCorrectFontSize:nohgaxes'));
end

% get the current scale from appdata stored in the axes
scl = getappdata(h.hgAxes,'Scale');

scl = min(1,scl);

for i =1:numel(h.Nodes)
    
    if ~isempty(h.Nodes(i).hgText)
        nfs = h.Nodes(i).FontSize*scl;
        if minAllowableFontSize>nfs
            set(h.Nodes(i).hgText,'Visible','off');
            set(h.Nodes(i).hgText,'FontSize',minAllowableFontSize);
        else
            set(h.Nodes(i).hgText,'Visible','on');
            set(h.Nodes(i).hgText,'FontSize',nfs);
        end
    end
end

if strcmp(h.showWeights,'on')
   edgeFontSize =  max(h.edgeFontSize * scl,minAllowableFontSize);
   lineObjs = get(h.Edges, 'hgline');
   textObjs = get([lineObjs{:}], 'UserData');
   set([textObjs{:}], 'FontSize', edgeFontSize);   
end
