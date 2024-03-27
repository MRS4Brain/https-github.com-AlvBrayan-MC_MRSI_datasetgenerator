function [ FreqMap,CCoefMap,RefSpectrumShort ] = MeasureFreqMap(RefTimeSerie, Data_trr,mrsiReconParams)


[~,MaxPPM_pt]=min(abs(mrsiReconParams.MaxPPM - mrsiReconParams.ppm));
[~,MaxPPM_pt_B0Corr]=min(abs(mrsiReconParams.LRTGVModelParams.CorrB0Map_MaxPPM - mrsiReconParams.ppm));
[~,MinPPM_pt_B0Corr]=min(abs(mrsiReconParams.LRTGVModelParams.CorrB0Map_MinPPM  - mrsiReconParams.ppm));

FreqRange=20;
FreqPrec=0.2;

Data_frr=fft(Data_trr,[],1);
Data_trr=ifft(Data_frr(MinPPM_pt_B0Corr:(end-MaxPPM_pt+MaxPPM_pt_B0Corr),:,:),[],1);


RefSpectrumShort=fft(RefTimeSerie,[],1);
RefSpectrumShort=RefSpectrumShort(MinPPM_pt_B0Corr:(end-MaxPPM_pt+MaxPPM_pt_B0Corr));
RefTimeSerie=ifft(RefSpectrumShort,[],1);

NbT=size(Data_trr,1);
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;
Time=([0 :(NbT-1)]'/Fs);
 FreqMap=zeros(size( mrsiReconParams.BrainMask2D));
 CCoefMap=zeros(size( mrsiReconParams.BrainMask2D));
 
for a=1:size(Data_trr,2)
    for b=1:size(Data_trr,3)
       % c=13;
            if mrsiReconParams.BrainMask2D(a,b)
                [FreqMap(a,b), CCoefMap(a,b)] = MeasureFreqShift( RefTimeSerie,Data_trr(:,a,b),Time,FreqRange,FreqPrec );
            end
    end   
end


FreqMap=FreqMap.*abs(CCoefMap);%/max(abs(CCoefMap(mrsiReconParams.BrainMask2D>0)));% weighting by correlation coeficient 

%smoothing of Maps
[NX,NY]=size(FreqMap);
sigma=mrsiReconParams.GaussianSigma;
[X,Y] = ndgrid(1:NX, 1:NY);
xc=floor(NX/2)+1;yc=floor(NY/2)+1;
exponent = -((X-xc).^2 + (Y-yc).^2)./(2*sigma^2);
Kernel = exp(exponent)/sum(exp(exponent(:))); % no need to normalize


FreqMap=FreqMap-mrsiReconParams.BrainMask2D.*mean(FreqMap(mrsiReconParams.BrainMask2D>0));

FreqMap=conv2(FreqMap,Kernel,'same');

CCoefMap=conv2(CCoefMap,Kernel,'same');

end

