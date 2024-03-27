function filtMRSI = mrsiExpFilter(mrsiData,tempSamplerate,FilterFreq)
% mrsiExpFilter -
%Apodizes the MRSI data with a exponential filter.
%
% INPUTS:
%           mrsiData: the input MRSI dataset ( T),
%                  T is the number of temporal
%                     samples                                 (double)
%     tempSamplerate: the temporal sampling rate in s              (double)
%                 Tw: the kernel length time in s           (double)
%        
%--------------------------------------------------------------------------
% OUTPUTS:
%           filtMRSI: the filtered MRSI dataset               (double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(size(mrsiData)>1)>1
    error('Error in mrsiHSVDFilter, mrsiData cointains more than one dimension.') 
end
T    = max(size(mrsiData));

DwellTime=1.0/tempSamplerate;
Tw=1.0/FilterFreq;
% Kernel
TimeKernel=exp(-([1:T]'*DwellTime)/Tw);
%TimeKernel=exp(-0.5*(([1:T]'*DwellTime)/Tw).^2);


filtMRSI     = squeeze(mrsiData).*TimeKernel ;


end

