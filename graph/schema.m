function schema
% schema for biograph

% Copyright 2003-2016 The MathWorks, Inc.

% set layout types
if isempty(findtype('BiographLayouts'))
	schema.EnumType('BiographLayouts',{'hierarchical','equilibrium','radial'});
end
% set layout types
if isempty(findtype('BiographEdgeType'))
	schema.EnumType('BiographEdgeType',{'straight','curved','segmented'});
end
% set text in Nodes
if isempty(findtype('BiographTextInNodes'))
	schema.EnumType('BiographTextInNodes',{'none','id','label'});
end

% add biograph class to package
bgpk = findpackage('biograph');
cls = schema.class(bgpk,'biograph');

% === public properties
p = schema.prop(cls,'ID','ustring');
p.FactoryValue =  '';

p = schema.prop(cls,'Label','ustring');
p.FactoryValue =  '';

p = schema.prop(cls,'Description','ustring');
p.FactoryValue =  '';

p = schema.prop(cls,'LayoutType','BiographLayouts');
p.FactoryValue = 'hierarchical';

p = schema.prop(cls,'LayoutScale','double');
p.FactoryValue = 1;
p.SetFunction = @verifyNonnegativeValue;

p = schema.prop(cls,'Scale','double');
p.FactoryValue = 1;
p.SetFunction = @verifyNonnegativeValue;

p = schema.prop(cls,'NodeAutoSize','on/off');
p.FactoryValue = 'on';

p = schema.prop(cls,'ShowTextInNodes','BiographTextInNodes');
p.FactoryValue = 'label';

p = schema.prop(cls,'EdgeType','BiographEdgeType');
p.FactoryValue = 'curved';

p = schema.prop(cls,'EdgeTextColor','color');
p.FactoryValue = [0 0 0];

p = schema.prop(cls,'ShowArrows','on/off');
p.FactoryValue = 'on';

p = schema.prop(cls,'ArrowSize','double');
p.FactoryValue = 8;
p.SetFunction = @verifyArrowSizeValue;

p = schema.prop(cls,'ShowWeights','on/off');
p.FactoryValue = 'off';

p = schema.prop(cls,'EdgeFontSize','double');
p.FactoryValue = 8;
p.SetFunction = @verifyNonnegativeValue;

p = schema.prop(cls,'NodeCallbacks','MATLAB array');
p.FactoryValue = @(node) inspect(node);
p.SetFunction = @verifyCallbackValue;

p = schema.prop(cls,'EdgeCallbacks','MATLAB array');
p.FactoryValue = @(edge) inspect(edge);
p.SetFunction = @verifyCallbackValue;

p = schema.prop(cls,'CustomNodeDrawFcn','MATLAB array');
p.FactoryValue = [];
p.SetFunction = @verifyCallbackValue;

p = schema.prop(cls,'Nodes','handle vector');
p.FactoryValue = handle([]);

p = schema.prop(cls,'Edges','handle vector');
p.FactoryValue = handle([]);

% === private properties
p = schema.prop(cls,'IsLaidout','bool'); 
p.FactoryValue = false;
p.Visible = 'off';

p = schema.prop(cls,'BoundingBox','MATLAB array');   
p.FactoryValue = [];
p.Visible = 'off';

p = schema.prop(cls,'DestroyListener','handle');
p.FactoryValue = [];
p.Visible = 'off';

p = schema.prop(cls,'to','MATLAB array');            
p.FactoryValue = [];
p.Visible = 'off';

p = schema.prop(cls,'from','MATLAB array');          
p.FactoryValue = [];
p.Visible = 'off';

% === handles to HG
p = schema.prop(cls,'hgAxes','MATLAB array');        
p.FactoryValue = [];
p.Visible = 'off';

%--------------------------------------------------------------------------
function valueStored = verifyArrowSizeValue(~, valueProposed)
% Validate the arrow size is integer and non-negative

% check it is a (scalar) integer and nonnegative
if isscalar(valueProposed) && valueProposed >= 0 && rem(valueProposed,1)==0
	valueStored = valueProposed;
else
	error(message('bioinfo:biographschema:InvalidArrowSize'))
end
%--------------------------------------------------------------------------
function valueStored = verifyNonnegativeValue(~, valueProposed)
% Validate the input value to set is non-negative

% check it is a (scalar) nonnegative
if isscalar(valueProposed) && valueProposed >= 0
	valueStored = valueProposed;
else
	error(message('bioinfo:biographschema:InvalidProperty'))
end
%--------------------------------------------------------------------------
function valueStored = verifyCallbackValue(~, valueProposed)
% Validate the custom callbacks for nodes and edges.

% transform to cell array
if ~iscell(valueProposed)
    valueProposed = {valueProposed};
end

% check each value in the cell array
for i = 1:numel(valueProposed)
    if ~isa(valueProposed{i}, 'function_handle')
        if ischar(valueProposed{i})
           if ~exist(valueProposed{i},'file') && ~exist(valueProposed{i},'builtin') && ~exist(valueProposed{i}, 'class')
			   error(message('bioinfo:biographschema:InvalidCallback'))
		   end
		else
            error(message('bioinfo:biographschema:InvalidCallback'))
        end
    end
end

if numel(valueProposed) == 1
    valueStored = valueProposed{1};
else
    valueStored = valueProposed;
end
			
