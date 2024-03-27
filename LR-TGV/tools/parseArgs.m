function out = parseArgs(args,defaults)
% Parses name/value argument pairs. Generally used when a function accepts
% a variable number of arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUTS:
%       args: A list of argument names, default values, and (optional)
%             argument types for each input. It can either be a (nx1),
%             (nx2), or (nx3) cell array corresponding to one of the
%             following forms:
%                   {Name; ...}
%                   {Name, DefaultValue; ...}
%                   {Name, DefaultValue, Type, ...}
%             A structure may also be passed, which is first converted to
%             the form above, or a cell row vector containing name value
%             pairs as:
%                   {Name1, DefaultValue, Name2, DefaultValue2, ...}
%        
%   defaults: The default name/value pairs, which will vary by application.
%             Can be in any of the formats listed above
%--------------------------------------------------------------------------
%   OUTPUTS:
%        out: A structure containing the parsed arguments
%--------------------------------------------------------------------------
% EXAMPLE USAGE: Generally used at the beginning of a function call
%       myfunction(arg1,arg2,varargin)
%   or  myfunction(arg1,arg2,optParams)
%               
%       defaults = {Name1, defVal1, type1; Name2, defVal2, type2;...}
%       out      = parseArgs(varargin(:),defaults);
%   or  out      = parseArgs(optParams, defaults);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args     = structurize(args);
defaults = parseDefaults(defaults);
out      = mergeArgs(args,defaults);

end
%--------------------------------------------------------------------------
function out = parseDefaults(x)
% Returns an (nx3) cell array in the form {Name,DefaultValue,Type; ...}
if isstruct(x)
    if ~isscalar(x)
        error('Error: Structure defaults must be scalar');
    end
    x = [fieldnames(x) struct2cell(x)];
end

if ~iscell(x)
    error('Invalid Defaults. Must be a structure or cell array');
end

if (size(x,1) == 1) && (size(x,2) > 1)
    if floor(numel(x)/2) == numel(x)/2
        x = transpose(reshape(x,[2 numel(x)/2]));
    else
        error('Error: Single row cell arrays must contain key/value pairs');
    end
end

if size(x,2) < 2
    x(:,2) = {[]};
end
if size(x,2) < 3
    x(:,3) = {[]};
end

out = x;    
end

%--------------------------------------------------------------------------
function out = structurize(x)
% Convert a structure or name/value list or record list to structure format

if isempty(x)
    out = struct;
elseif iscell(x)
    if (size(x,1) == 1) && (size(x,2) > 1)
        if floor(numel(x)/2) == numel(x)/2
            x = transpose(reshape(x,[2 numel(x)/2]));
        else
            error('Error: Single row cell arrays must contain key/value pairs');
        end
    end
    if size(x,2) ~= 2
        error('Invalid arguments: cells must be (nx2) {Name, Val; ...} or row vector {Name, Val, ...} list');
    end
    x   = x.';
    x   = x(:);
    out = struct(x{:});
elseif isstruct(x)
    if ~isscalar(x)
        error('Error: Structural arguments must be scalar');
    end
    out = x;
end

end

%--------------------------------------------------------------------------
function out = mergeArgs(args,defaults)
% Applies user arguments, checking the validity of input types
out   = structurize(defaults(:,[1 2]));
names = fieldnames(args);

for l = 1:numel(names)
    out.(names{l}) = args.(names{l});
end

% Check for valid argument types
for l = 1:size(defaults,1)
    [name,~,type] = defaults{l,:};
    if ~isempty(type)
        out.(name) = needa(type, out.(name), name);
    end
end

end

%--------------------------------------------------------------------------
function out = needa(type,value,name)    
% Check that a value is of a given type
switch type
    case 'cellstr'
        isType = iscellstr(value);
    case 'datenum'
        isType = isnumeric(value);
    otherwise
        isType = isa(value,type);
end

if isType
    out = value;
else
    error('Error: Argument %s must be a %s; found a %s', name, type, class(value));
end

end

%--------------------------------------------------------------------------