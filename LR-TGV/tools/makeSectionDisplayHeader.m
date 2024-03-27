function dispHeader = makeSectionDisplayHeader(headerStr, varargin)
% makeSectionDisplayHeader - 
% Formats a section display header for console output

if (nargin == 1)
    
    char = '*';   

elseif (nargin == 2)
    
    char = varargin{1};
    
elseif (nargin < 1)
    
    error('makeSectionDisplayHeader:TooFewArguments', ...
        'Too few arguments');
    
else
    
    error('makeSectionDisplayHeader:TooManyArguments', ...
        'Too many arguments');
    
end

if ~ischar(headerStr)
    
    error('dispHeader:badInputString', ...
        'First input must be a valid string');
end

if ~(isscalar(char) && ischar(char))
    
    error('dispHeader:badCharacter', ...
        'Style must be a single character');
end

% Determine the current width of the MATLAB console
consoleWidth = matlab.desktop.commandwindow.size;

% Create the display header string
dividerStr  = sprintf([repmat(char, [1, consoleWidth(1)]), '\n']);
dispHeader  = sprintf([dividerStr, ' ', headerStr, '\n', dividerStr]);

end

