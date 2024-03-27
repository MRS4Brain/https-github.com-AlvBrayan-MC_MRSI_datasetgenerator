function [filtMRSI ,  WaterMRSI, WaterFreq , WaterAmp] = mrsiHSVDFilter(mrsiData,scannerFrequency,tempSamplerate,ppmOffset,modelOrder,Water_minfreq,Water_maxfreq,varargin)
% mrsiHSVDFilter -
% Estimates and removes frequency components from each spectrum in an MRSI 
% dataset.
%
% INPUTS:
%           mrsiData: the input MRSI dataset ( T), T is the number of temporal
%                     samples                                 (double)
%   scannerFrequency: frequency of the spectrometer           (double)
%     tempSamplerate: the temporal sampling rate              (double)
%          ppmOffset: the offset for the ppm axis             (double)
%         modelOrder: the number of expected complex exponentials
%                     comprising the input signal             (double)
%--------------------------------------------------------------------------   
% OPTIONAL INPUTS:
%   A series of name/value pairs as described below:
%         'ppmRange': a 2-element vector containing the minimum and maximum 
%                     ppm values circumscribing the range of frequencies
%                     to remove                               (double)
%          'verbose': flag indicating whether processing information should
%                     be printed to the command window        (logical)
%--------------------------------------------------------------------------
% OUTPUTS:
%           filtMRSI: the new MRSI dataset with the indicated frequencies
%                     removed                                 (double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(size(mrsiData)>1)>1
    error('Error in mrsiHSVDFilter, mrsiData cointains more than one dimension.') 
end
    
%[M,T]    = size(mrsiData);
T    = max(size(mrsiData));
ppmScale = buildppmscale(T, tempSamplerate, scannerFrequency, ppmOffset);

defaults = {'ppmRange', [min(ppmScale) max(ppmScale)], 'double';
            'verbose' , true                         , 'logical'};

inParams = parseArgs(varargin,defaults);

% Determine the search range in Hz
%freqSearchRange = -(inParams.ppmRange - ppmOffset) * scannerFrequency;
%minFreq         = min(freqSearchRange);
%maxFreq         = max(freqSearchRange);
minFreq         = Water_minfreq;%in hz
maxFreq         = Water_maxfreq;%in hz
filtMRSI = zeros(size(mrsiData));

% Loop over all spectra in the MRSI dataset
%for spec = 1:M
    %origSignal           = squeeze(mrsiData(spec, :));
   % [freqs,~,basis,amps] = hsvd(transpose(origSignal), tempSamplerate, modelOrder);
     [freqs,~,basis,amps] = HSVD(mrsiData, tempSamplerate, modelOrder);
      %  freqs
    indx                 = find((freqs >= minFreq) & (freqs <= maxFreq));
   
  %  newSignal            = sum(basis(:, indx) * diag(amps(indx), 0), 2);
    %filtMRSI(spec,:)     = origSignal - newSignal.';
    
    % filtMRSI(spec,:) = squeeze(mrsiData(spec, :))-sum(basis(:, indx) * diag(amps(indx), 0), 2).';
 
   filtMRSI = mrsiData-sum(basis(:, indx) * diag(amps(indx), 0), 2);
   WaterMRSI = sum(basis(:, indx) * diag(amps(indx), 0), 2);
   WaterAmp = sum(abs(WaterMRSI));
   if(numel(indx)>0)
        WaterFreq = sum(abs(amps(indx)).*freqs(indx))/sum(abs(amps(indx)));
   else
        WaterFreq = 0;
   end
   % if inParams.verbose
    %    fprintf('Processing spectrum %d of %d\n', spec, M);
    %end
%end

end

