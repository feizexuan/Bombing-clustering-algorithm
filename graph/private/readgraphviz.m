function  readgraphviz(h,outFilename)
% READGRAPHVIZ parses GRAPHVIZ output into a BIOGRAPH object.

% Copyright 2003-2012 The MathWorks, Inc.


% open file for reading
fid = fopen(outFilename,'r');
if fid < 0
    error(message('bioinfo:biograph:readgraphviz:FileOpenError', outFilename)); 
end

% get file line by line
i = 1;
fileLine = fgetl(fid);
while (fileLine ~= -1)    
    u = regexp(fileLine,{'\[','\]'});
    
    % process line
    if numel(u{1})
        key{i}   = fileLine(2:u{1}-2); %#ok
        attrib{i} = fileLine(u{1}+1:u{2}-1); %#ok
        i = i+1;
    end
    
    % get next line
    fileLine = fgetl(fid);
    while fileLine(end)=='\'
        % remove line continuations
        fileLine = [fileLine(1:end-1) fgetl(fid)];
    end
end
fclose(fid);

% remove '"' placed to help graphviz with difficult chars
key = regexprep(key,'"','');

% separate keys and attribs by type
graphAttrib = attrib(ismember(key,'graph'));
edgePtrs = find(~cellfun('isempty',strfind(key,'->')));
nodePtrs = setdiff(find(~(ismember(key,{'graph','node','edge'}))),edgePtrs);

% get Bounding Box
bb = sscanf(graphAttrib{strmatch('bb',graphAttrib)},'bb="%f,%f,%f,%f')';
h.BoundingBox = bb.*h.LayoutScale;
h.IsLaidout = true;

hash_table = java.util.HashMap;
for i = 1:numel(h.Nodes)
   hash_table.put(num2str(i),i);
end

% get Position and Extent for every node
for i = nodePtrs
    n = h.Nodes(hash_table.get(key{i}));
    n.Position = parseAttribs(attrib{i},'pos')*h.LayoutScale;
    n.Size = 72 * h.LayoutScale * [parseAttribs(attrib{i},'width') ...
                                   parseAttribs(attrib{i},'height')];
    % 72 ppi seems to be the default unless post-scaled
end

ns = regexp(key(edgePtrs),'(.+) -> (.+)','tokens','once');
for i = 1:numel(ns) 
    e = h.edges(h.to(hash_table.get(ns{i}{1}),hash_table.get(ns{i}{2})));
    e.ControlPoints = h.LayoutScale * parseEdgeControlPoints(attrib{edgePtrs(i)});
end

clear hash_table

% ======================================================
function attrib = parseAttribs(str,prop)
%PARSEATTRIBS parse the property from the given string
[s,f,t] = regexp(str,[prop '=(\w+),??']);  
attrib = '';
if isempty(s) % not a word, might be a number
    [s,f,t] = regexp(str,[prop '="([^"]*)",??']); 
    if ~isempty(t)
        attrib = str2num(str(t{1}(1):t{1}(2))); %#ok
    end
else % word
    attrib = str(t{1}(1):t{1}(2));
end

% ======================================================
function  cp = parseEdgeControlPoints(str)
%PARSEEDGECONTROLPOINTS parse control points from string
str = regexp(str,'pos="([^"]+)','tokens','once');
startStr = regexp(str,'s,[\d]+,[\d]+','match','once');
startPt = sscanf(startStr{:},'s,%d,%d');
str = regexprep(str,'s,[\d]+,[\d]+','');
endStr = regexp(str,'e,[\d]+,[\d]+','match','once');
endPt = sscanf(endStr{:},'e,%d,%d');
str = regexprep(str,'e,[\d]+,[\d]+','');
cp = sscanf(str{:},'%d,%d');
cp = reshape(cp,2,numel(cp)/2);
if isempty(startPt)&&~isempty(endPt)
    cp = [cp endPt];
elseif ~isempty(startPt)&&isempty(endPt)
    cp = [fliplr(cp) startPt];
elseif ~isempty(startPt)&&~isempty(endPt)
    cp = [startPt cp endPt] ;
end
