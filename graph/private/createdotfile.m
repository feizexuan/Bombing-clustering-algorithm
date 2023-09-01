function inFilename = createdotfile(h,doOnlyPaths)
% CREATEDOTFILE creates the input for GRAPHVIZ given a BIOGRPH object.

% Copyright 2003-2012 The MathWorks, Inc.


% generate temporary filenames
inFilename = [tempname '.dot'];

%open input file to write
fid = fopen(inFilename,'w');
if fid < 0
    error(message('bioinfo:biograph:createdotfile:FileOpenError', inFilename)); 
end

% start graph 
fprintf(fid,'digraph mygraph{\n');
fprintf(fid,'	maxiter=99;\n');
fprintf(fid,'	nslimit=99;\n');
fprintf(fid,'	nslimit1=99;\n');
% splines and overlap mess twopi up in some unix platforms
if doOnlyPaths || ispc || ~isequal(h.LayoutType,'radial')
   fprintf(fid,'	splines=true;\n');
   fprintf(fid,'	overlap=false;\n');
end

nodeOrder = 1:numel(h.Nodes); 
% for twopi we need to put the unconnected nodes at the end, otherwise
% twopi errors out.
if isequal(h.LayoutType,'radial') && ~doOnlyPaths
    NonConnectedNodes = ~(full(any(h.to))|full(any(h.from)));
    nodeOrder = [nodeOrder(~NonConnectedNodes) nodeOrder(NonConnectedNodes)];
end

% keep track of maximum width and height used
maxWidth = 0; maxHeight = 0;

if strcmp(h.ShowTextInNodes,'label') 
    if numel(nodeOrder)>1
        if all(cellfun('isempty',get(h.Nodes,'Label')))
            h.ShowTextInNodes = 'id';
        end
    elseif isempty(h.Nodes.Label)
        h.ShowTextInNodes = 'id';
    end
end
    
if isequal(h.NodeAutoSize,'on')
    % create a phantom figure to calculate the extents for the labels
    hf = figure('visible','off');
    ha = newplot;
    set(ha,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1]);
    ht = text(0,0,'','Interpreter','none','Units', 'inches','VerticalAlignment','bottom');
    % nodes
    for k = 1:numel(nodeOrder)
        i = nodeOrder(k);
        shape =  h.Nodes(i).Shape;
        switch shape
            case {'rect','rectangle'}
                shape = 'box';
            case 'diamond'
                shape = 'ellipse';
        end
        
		switch h.ShowTextInNodes
			case 'label'
				label = h.Nodes(i).Label;
			case 'id'
				label = h.Nodes(i).id;
			case 'none'
				label = '';
		end
		label = regexprep(label,'^(default|remove|factory)$','\\$0');
		
        % Get the extents of the label and add a little bit of space around the text.
        % Thanks to Brad Gulko for his suggestion for this code.
        set(ht,'string',label,'FontSize',h.Nodes(i).Fontsize);
        extents = get(ht,'Extent');    
        height  = extents( 4 ); 
        width   = extents( 3 ) + height*2/5; 

        switch h.Nodes(i).Shape;
            case 'circle'
                height = max(height,width);
                width = height;
            case {'house','invhouse'}
                height = height * 6/5;
            case {'trapezium','invtrapezium','parallelogram'}
                width =  width + height * 2/5;
            case 'diamond'
                width =  width + height * 2/5;
                height = height * 7/5;
        end
        width = width/h.LayoutScale;
        height = height/h.LayoutScale;
        maxWidth = max(maxWidth,width); maxHeight = max(maxHeight,height);
        switch doOnlyPaths
            case false
                if isequal(h.LayoutType,'radial') && NonConnectedNodes(i)
                    fprintf(fid,'    "%s" [shape=%s,height=%f,width=%f,fixedsize=true,fontsize=%d];\n',...
                        num2str(i),shape,maxHeight,maxWidth,h.Nodes(i).Fontsize);
                else
                    fprintf(fid,'    "%s" [shape=%s,height=%f,width=%f,fixedsize=true,fontsize=%d];\n',...
                        num2str(i),shape,height,width,h.Nodes(i).Fontsize);
                end
            case true
                fprintf(fid,'    "%s" [pos="%f,%f",shape=%s,height=%f,width=%f,fixedsize=true,fontsize=%d];\n',...
                    num2str(i),h.Nodes(i).Position(1)/h.LayoutScale,...
                               h.Nodes(i).Position(2)/h.LayoutScale,...
                    shape,height,width,h.Nodes(i).Fontsize);
        end
    end
    close(hf)

else % ~isequal(h.NodeAutoSize,'on')
    for k = 1:numel(nodeOrder)
        i = nodeOrder(k);
        shape =  h.Nodes(i).Shape;
        switch shape
            case {'rect','rectangle'}
                shape = 'box';
            case 'diamond'
                shape = 'ellipse';
        end     
        width = h.Nodes(i).Size(1)/h.LayoutScale/72;
        height = h.Nodes(i).Size(2)/h.LayoutScale/72;
        maxWidth = max(maxWidth,width); maxHeight = max(maxHeight,height);
        switch doOnlyPaths
            case false
                if isequal(h.LayoutType,'radial') && NonConnectedNodes(i)
                    fprintf(fid,'    "%s" [shape=%s,height=%f,width=%f,fixedsize=true,fontsize=%d];\n',...
                       num2str(i),shape,maxHeight,maxWidth,h.Nodes(i).Fontsize);
                else
                   fprintf(fid,'    "%s" [shape=%s,height=%f,width=%f,fixedsize=true,fontsize=%d];\n',...
                       num2str(i),shape,height,width,h.Nodes(i).Fontsize);
                end
            case true
                fprintf(fid,'    "%s" [pos="%f,%f",shape=%s,height=%f,width=%f,fixedsize=true,fontsize=%d];\n',...
                    num2str(i),h.Nodes(i).Position(1)/h.LayoutScale,...
                               h.Nodes(i).Position(2)/h.LayoutScale,...
                    shape,height,width,h.Nodes(i).Fontsize);
        end
    end
end

% edges
% other directions: 'back','both','none'
for i = 1:numel(h.Edges)
    name1 = num2str(find(h.Nodes == h.Edges(i).FromNode));
    name2 = num2str(find(h.Nodes == h.Edges(i).ToNode));
    %len = h.Edges(i).Weight;
    %fprintf(fid,'    "%s" -> "%s" [dir=forward,len=%f];\n',name1,name2,20*len);
    fprintf(fid,'    "%s" -> "%s" [dir=forward];\n',name1,name2);
end

% finish graph
fprintf(fid,'}\n');

fclose(fid);
