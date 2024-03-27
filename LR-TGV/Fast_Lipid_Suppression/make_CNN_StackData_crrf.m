function Data_stack_rf=make_CNN_StackData_crrf( mrsiData_crrf,Bound,mrsiReconParams,NameStackData)
N1=Bound(1);
N2=Bound(2);
N = size(mrsiData_crrf);

Npt=N(1);
NMRfreq=mrsiReconParams.mrProt.NMRFreq*1E6;


mrsiData_crrf=mrsiData_crrf(:,:,:,N1:N2);
Data_stack_rf=mrsiData_crrf;
clear mrsiData_crrrf
Data_stack_rf=reshape(Data_stack_rf,[],(N2-N1+1));

hdf5write(NameStackData,'realData',real(Data_stack_rf),'imagData',imag(Data_stack_rf),...
    'samplerate',mrsiReconParams.mrProt.samplerate,...
    'Npt',Npt,'N1',N1-1,'N2',N2,'NMRfreq',NMRfreq);

end
