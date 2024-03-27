function Data_rrt = AdjSpatSpectOp(Data_ckkt, reconParams)

T=size(Data_ckkt,4);
Data_rr   = zeros(size(squeeze(Data_ckkt(1,:,:,:))), class(Data_ckkt));
nDims = ndims(Data_ckkt)-1;


FreqMap=reconParams.WaterFreqMap ;
Fs=Fs=reconParams.mrProt.samplerate*NbT/reconParams.mrProt.VSize;
DelayT=reconParams.DelayT;

for l = 1 : T
   temp=AdjEncodOp(squeeze(Data_ckkt(:,:,:,l)), reconParams.kmask,reconParams.SENSE);
   Data_rrt(:,:,l) = exp(-2*pi*1i*((l-1)/ Fs+DelayT) *  FreqMap) .*temp;
end

end
