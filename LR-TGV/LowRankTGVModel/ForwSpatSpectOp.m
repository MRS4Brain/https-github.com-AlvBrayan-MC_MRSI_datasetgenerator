function Data_ckkt = ForwSpatSpectOp(Data_rrt, reconParams)



NbT=size(Data_rrt,3);
nDims = ndims(Data_rrt);
%spa_index = repmat({':'}, 1, nDims-1);
FreqMap=reconParams.WaterFreqMap ;
Fs=reconParams.mrProt.samplerate*NbT/reconParams.mrProt.VSize;
DelayT=reconParams.DelayT;


Time_rrt=permute(repmat(([0 :(NbT-1)]'/Fs+DelayT),[1 size(Data_rrt,1) size(Data_rrt,2)]),[2,3,1]);

%Freqshift_rrt=exp(2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT]));

RepData_crrt=permute(repmat(exp(2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT])).*Data_rrt,[1 1 1 size(reconParams.SENSE,1)]),[4,1,2,3]);
kmask_ckkt=permute(repmat(reconParams.kmask,[1 1 size(reconParams.SENSE,1) NbT]),[3,1,2,4]);
%SENSE_crrt=repmat(reconParams.SENSE,[1 1 1 NbT]);
Data_ckkt = fft(fft(repmat(reconParams.SENSE,[1 1 1 NbT]).*RepData_crrt,[],2),[],3).*kmask_ckkt;

%{
Data_ckkt   = zeros([  size(reconParams.SENSE,1) size(Data_rrt)] , class(Data_rrt));
for l = 1 : NbT
    Data_ckkt(:,:,:,l) = ForwEncodOp(exp(2*pi*1i*((l-1)/ Fs+DelayT) *  FreqMap ) .* squeeze(Data_rrt(:,:,l)), reconParams.kmask,reconParams.SENSE);
end
%}

end
