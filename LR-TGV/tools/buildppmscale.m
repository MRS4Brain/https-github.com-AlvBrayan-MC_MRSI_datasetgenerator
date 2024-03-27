function ppmscale = buildppmscale(T,Fs,specfreq,ppmoffset)
%Builds the ppm scale to be used for plotting NMR data
% PARAMETERS --------------------------------------------------------------
%   INPUTS: 
%             T: Number of sampled points
%            Fs: temporal sampling frequency (samp/sec)
%      specfreq: Frequency of the spectrometer (MHz)
%     ppmoffest: The desired ppm value corresponding to the
%                spectrometer frequency
%
%   OUTPUTS:
%      ppmscale: The ppm scale
%--------------------------------------------------------------------------
%NOTES:
%Typical usage with NMR spectra would be using the following plot command:
% plot(ppmscale,real(fftshift(fft(signal)))), set(gca,'XDir','reverse')
% which plots the Fourier transform of a given FID (signal), using the
% convention that ppm values should increase from right to left
%--------------------------------------------------------------------------

ppmscale = (-floor(T/2):ceil(T/2-1))*(Fs/T)*(1/specfreq) + ppmoffset;

end

