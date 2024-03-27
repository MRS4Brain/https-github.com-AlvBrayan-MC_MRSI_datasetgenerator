function Data_stack_rf=make_CNN_StackData_frrr( mrsiData_frrr,Mask,Bound,mrsiReconParams,NameStackData)
N1=Bound(1);
N2=Bound(2);
N = size(mrsiData_frrr);
HzpP=mrsiReconParams.mrProt.samplerate/N(1);
Npt=N(1);
NMRfreq=mrsiReconParams.mrProt.NMRFreq*1E6;
%lipid_mask=mrsiReconParams.SkMask;

Mask_rrrf=repmat(Mask,[1 1 1 (N2-N1+1)]);
data_rrrf=permute(mrsiData_frrr,[2,3,4,1]);
data_rrrf=data_rrrf(:,:,:,N1:N2);
Data_stack_rf=data_rrrf(Mask_rrrf>0);
Data_stack_rf=reshape(Data_stack_rf,[],(N2-N1+1));
hdf5write(NameStackData,'realData',real(Data_stack_rf),'imagData',imag(Data_stack_rf),...
    'samplerate',mrsiReconParams.mrProt.samplerate,...
    'Npt',Npt,'N1',N1-1,'N2',N2,'NMRfreq',NMRfreq);

end
