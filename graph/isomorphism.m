function varargout =  isomorphism(varargin)
%ISOMORPHISM finds an isomorphism between two graphs.
%
%  ISOMORPHISM method for BIOGRAPH objects extracts the adjacency matrix
%  and calls GRAPHISOMORPHISM function. All other input and output
%  arguments are analogous to GRAPHISOMORPHISM.
% 
%  See also graphisomorphism.

% Copyright 2006-2016 The MathWorks, Inc.

 if nargin > 0
     [varargin{:}] = convertStringsToChars(varargin{:});
 end
 
 if numel(varargin) < 2
     error(message('bioinfo:biograph:isomorphism:IncorrectNumberOfArguments'));
 elseif ~isa(varargin{1},'biograph.biograph') ||  ~isa(varargin{2},'biograph.biograph')
      error(message('bioinfo:biograph:isomorphism:InvalidInput'));
 end
      
[varargout{1:nargout}] = graphisomorphism(getmatrix(varargin{1}),getmatrix(varargin{2}),varargin{3:end});
