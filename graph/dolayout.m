function dolayout(h,varargin)
%DOLAYOUT calculates node positions and edge trajectories
%
%   DOLAYOUT(BIOGRAPH) calls the layout engine to calculate the optimal
%   position for each node so that its 2-D rendering is clean and
%   uncluttered, and then calculates the best curves to represent the
%   edges. The following BIOGRAPH properties interact with the layout
%   engine:
%
%       'LayoutType'   Selects the layout engine from three available
%                      options; hierarchical, radial and equilibrium. The
%                      default option is hierarchical.
%       'LayoutScale'  Rescales the size of the node sizes before calling
%                      the layout engine. This gives more space to the
%                      layout and reduces the overlapping of nodes.
%       'NodeAutoSize' When 'on' it uses the node properties ('FontSize'
%                      and 'Shape') and 'LayoutScale' to precalculate the
%                      actual size of every node. When 'off' the layout
%                      engine uses the node 'Size' property.
%
%   DOLAYOUT(BIOGRAPH,'PathsOnly', TRUE) calls the layout engine; leaving
%   the nodes at their current positions and only calculating new curves
%   for the edges. Default is FALSE.
%
%   Examples:
%
%      % Create a BIOGRAPH object.
%      cm = [0 1 1 0 0;1 0 0 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 0 0];
%      bg = biograph(cm)
%      bg.nodes(1).Position   % <== nodes do not have a position yet
%
%      % Call the layout engine and render the graph.
%      dolayout(bg)
%      bg.nodes(1).Position
%      view(bg)
%
%      % Manually modify the node position and recalculate the paths.
%      bg.nodes(1).Position = [150 150];
%      dolayout(bg,'PathsOnly',true)
%      view(bg)
%
%   See also BIOGRAPH/BIOGRAPH, BIOGRAPH.BIOGRAPH/VIEW,
%   BIOGRAPH.BIOGRAPH/GETNODESBYID, BIOGRAPH.BIOGRAPH/GETEDGESBYNODEID.

% Copyright 2003-2014 The MathWorks, Inc.


% default
doOnlyPaths = false;

% get input arguments
nvarargin = numel(varargin);
if  nvarargin
    if rem(nvarargin,2)
        error(message('bioinfo:biograph:dolayout:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'pathsonly',''};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, length(pname)));
        if isempty(k)
            error(message('bioinfo:biograph:dolayout:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:biograph:dolayout:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % 'pathsonly'
                    doOnlyPaths = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            end
        end
    end
end

if doOnlyPaths && ~h.isLaidout
    error(message('bioinfo:biograph:dolayout:NoLayout'))
end

debugMakeGraph = false;
debugViewImage = false; % <= if true, it requires a local Graphviz install

% generate temporary for output
outFilename =  [tempname '.dot'];

% prepare Graphviz input file
inFilename = createdotfile(h,doOnlyPaths);

% open input file for debugging
if debugMakeGraph
    edit(inFilename) %#ok<UNRCH>
end

change_dir = false;
try
    % In case we have a UNC path and running from PC we need to change to
    % the temp directory.
    current_dir = pwd;
    if ispc && strncmp(current_dir,'\\',2)
        change_dir = true;
        cd(tempdir); %#ok<*MCCD>
    end
    
    % call Graphviz
    if doOnlyPaths
        stat = callgraphviz('neato','-nop -s -Tdot',inFilename,'-o',outFilename);
    else
        switch h.LayoutType
            case 'hierarchical'
                layoutEngine = 'dot';
            case 'radial'
                layoutEngine = 'twopi';
            case 'equilibrium'
                layoutEngine = 'neato';
        end
        stat = callgraphviz(layoutEngine,'-Tdot',inFilename,'-o',outFilename);
    end
    
catch allExceptions
    if change_dir
        cd(current_dir);
    end
    error(message('bioinfo:biograph:dolayout:graphvizHardError'))
end

if change_dir
    cd(current_dir);
end


% check Graphviz execution output
if stat ~= 0
    if exist(which('callgraphviz'),'file')
        disp('Starting GraphViz Output:')
        [~,theMessage] = callgraphviz(layoutEngine,'-Tdot',inFilename,'-o',outFilename,'-v'); %#ok<NASGU>
        disp('End GraphViz Output.')
    end
    error(message('bioinfo:biograph:dolayout:graphvizError'))
end

% open output file for debugging
if debugMakeGraph
    edit(outFilename) %#ok<UNRCH>
end

% parse Graphviz output
readgraphviz(h,outFilename)

% view the image generated by graphviz
if debugViewImage
    % need to install the latest graphviz distribution
    outImage = [tempname '.jpeg']; %#ok<UNRCH>
    localGraphvizPath = 'D:\Applications\ATT\Graphviz\bin\';
    str = sprintf('"%s%s" -Tjpeg %s -o %s',localGraphvizPath,...
        layoutEngine, inFilename, outImage);
    [stat, msg] = system(str); 
    imtool(outImage)
    delete(outImage)
end

if ~debugMakeGraph
    delete(inFilename)
    delete(outFilename)
end

if ~isempty(h.hgAxes)
    h.hgUpdate
end

