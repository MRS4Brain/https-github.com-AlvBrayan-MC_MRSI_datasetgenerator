function Data_stack_rf=make_CNN_StackData( mrsiData_ctkkk,Mask,Bound,mrsiReconParams,NameStackData)
N1=Bound(1);
N2=Bound(2);
N = size(mrsiData_ctkkk);
HzpP=mrsiReconParams.mrProt.samplerate/N(2);
Npt=N(2);
NMRfreq=mrsiReconParams.mrProt.NMRFreq*1E6;
%lipid_mask=mrsiReconParams.SkMask;

Mask_kkkf=repmat(Mask,[1 1 1 (N2-N1+1)]);
NbptInMask=sum(Mask_kkkf(:)>0);
Data_stack_rf=zeros([N(1),NbptInMask]);
for c=1:N(1)
   % c = c
    data_kkkf=permute(fft(mrsiData_ctkkk(c,:,:,:,:),[],2),[3,4,5,2,1]);
    data_kkkf=data_kkkf(:,:,:,N1:N2);
   
    %Data_stack_rf=[Data_stack_rf ; reshape(data_kkkf(Mask_rrrf>0),[],(N2-N1+1))];
    Data_stack_rf(c,:)=data_kkkf(Mask_kkkf>0);
end
Data_stack_rf=reshape(Data_stack_rf,[],(N2-N1+1));

%Data_stack_rf=permute(ifft(ifft(ifft(fft(mrsiData_ctkkk,[],2),[],3),[],4),[],5),[1,3,4,5,2]);%crrrf
%Mask_crrrf=logical(permute(repmat(logical(Mask),[1 1 1 N(1) N(2)]),[4,1,2,3,5]));
%Data_stack_rf=reshape(Data_stack_rf(Mask_crrrf>0),[],N(2));

hdf5write(NameStackData,'realData',real(Data_stack_rf),'imagData',imag(Data_stack_rf),...
    'samplerate',mrsiReconParams.mrProt.samplerate,...
    'Npt',Npt,'N1',N1-1,'N2',N2,'NMRfreq',NMRfreq);

end
