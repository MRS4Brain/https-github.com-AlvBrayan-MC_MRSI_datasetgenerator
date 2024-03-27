function x_gpu = transferToGPU(x)
% transferToGPU -
% Transfers the input to the currently selected GPU
% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   INPUTS:
%           x: The numerical array to transfer to the GPU. If x is a 
%              container object such as a MATLAB struct or cell array, 
%              the function is called recursively on each item or element
%              until a numerical array is found
%--------------------------------------------------------------------------
%   OUTPUTS:
%       x_gpu: The output gpuArray object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% Check the class of the input and act accordinly

if isnumeric(x)
    
    x = gpuArray(x);
    
elseif isstruct(x)
    % Loop through each item of the input structure
    
    fldNames = fieldnames(x);
    
    for elem = 1 : numel(x)
    
        for name = 1 : numel(fldNames)
            
            x(elem).(fldNames{name}) = transferToGPU(x(elem).(fldNames{name}));
            
        end
    
    end
            
elseif iscell(x)
    % Loop through each element of the cell array
    
    for elem = 1 : numel(x)
        
        x{elem} = transferToGPU(x{elem});
        
    end
    
elseif islogical(x)
    % Do nothing
    
else
    fprintf(['Warning: Objects of type ''%s'' are not currently supported ', ...
        'on the GPU.\nSuch objects were not transferred.\n'], class(x));

end

x_gpu = x;

end

