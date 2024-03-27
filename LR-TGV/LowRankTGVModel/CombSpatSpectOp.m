function DataOut_rrt = CombSpatSpectOp(Data_rrt, reconParams)

DataOut_rrt   = zeros(size(Data_rrt), class(Data_rrt));
nDims = ndims(Data_rrt);
NbT=size(Data_rrt,3);
index = repmat({':'}, 1, nDims-1);
FreqMap=reconParams.WaterFreqMap ;
Fs=Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;;
DelayT=reconParams.DelayT;
NbT = size(Data_rrt,3);
%{
for l = 1 :  size(Data_rrt,3)
    DataOut_rrt(index{:}, l) = exp(-2*pi*1i*((l-1)/ Fs+DelayT) * FreqMap ) .* ...
                        CombEncodOp(exp(2*pi*1i*((l-1)/ Fs+DelayT) * FreqMap ) .* ...
                        squeeze(Data_rrt(:,:,l)),reconParams.kmask,reconParams.SENSE);
end
%}
Time_rrt=permute(repmat(([0 :(NbT-1)]'/Fs+DelayT),[1 size(Data_rrt,1) size(Data_rrt,2)]),[2,3,1]);

Freqshift_rrt=exp(2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT]));

RepData_crrt=permute(repmat(Freqshift_rrt.*Data_rrt,[1 1 1 size(reconParams.SENSE,1)]),[4,1,2,3]);
kmask_ckkt=permute(repmat(reconParams.kmask,[1 1 size(reconParams.SENSE,1) NbT]),[3,1,2,4]);
SENSE_crrt=repmat(reconParams.SENSE,[1 1 1 NbT]);
Data_ckkt = fft(fft(SENSE_crrt.*RepData_crrt,[],2),[],3).*kmask_ckkt;
DataOut_rrt=Freqshift_rrt.*squeeze(sum(conj(SENSE_crrt).*ifft(ifft(Data_ckkt,[],2),[],3),1)) ;
end
