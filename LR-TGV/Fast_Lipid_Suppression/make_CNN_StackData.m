function Data_stack_rf=make_CNN_StackData( mrsiData_ctkk,Mask,Bound,mrsiReconParams,NameStackData)
N1=Bound(1);
N2=Bound(2);
N = size(mrsiData_ctkk);
HzpP=mrsiReconParams.mrProt.samplerate/N(2);
Npt=N(2);
NMRfreq=mrsiReconParams.mrProt.NMRFreq*1E6;
%lipid_mask=mrsiReconParams.SkMask;

Mask_rrf=repmat(Mask,[1 1 (N2-N1+1)]);
NbptInMask=sum(Mask_rrf(:)>0);
Data_stack_rf=zeros([N(1),NbptInMask]);
for c=1:N(1)
   % c = c
    data_rrf=permute(ifft(ifft(fft(mrsiData_ctkk(c,:,:,:),[],2),[],3),[],4),[3,4,2,1]);
    data_rrf=data_rrf(:,:,N1:N2);
   
    %Data_stack_rf=[Data_stack_rf ; reshape(data_rrrf(Mask_rrrf>0),[],(N2-N1+1))];
    Data_stack_rf(c,:)=data_rrf(Mask_rrf>0);
end
Data_stack_rf=reshape(Data_stack_rf,[],(N2-N1+1));


hdf5write(NameStackData,'realData',real(Data_stack_rf),'imagData',imag(Data_stack_rf),...
    'samplerate',mrsiReconParams.mrProt.samplerate,...
    'Npt',Npt,'N1',N1-1,'N2',N2,'NMRfreq',NMRfreq);

end
