function x_new = gatherFromGPU(x)
% sendParamsToGPU -
% Gathers a variable from the currently selected GPU
% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   INPUTS:
%           x: The object to gather that is currently sitting on the
%              selected GPU. X_GPU can also be a container class, such as a
%              MATLAB struct or cell array, that contains gpuArrays as
%              elements. In this case, the function is called recursively
%              on each item or element until a gpuArray object is found
%--------------------------------------------------------------------------
%   OUTPUTS:
%       x_new: The output object, which should only contain items that
%              exist on the CPU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmpi(class(x), 'gpuArray'))
    
    x = gather(x);
    
elseif isstruct(x)
    % Loop through each item of the input structure
    
    fldNames = fieldnames(x);
    
    for elem = 1 : numel(x)
        
        for name = 1 : numel(fldNames)
            
            x(elem).(fldNames{name}) = gatherFromGPU(x(elem).(fldNames{name}));
            
        end
        
    end
    
elseif iscell(x)
    % Loop through each element of the cell array
    
    for elem = 1 : numel(x)
        
        x{elem} = gatherFromGPU(x{elem});
        
    end
    
end

x_new = x;

end

