%% ------------------------------------------------------------------------
function Data_ckk =  ForwEncodOp(Data_rr,kmask,SENSE)

%imDims = size(Data_rr);
%for c=1:size(SENSE,1)
 %  Data_ckk(c,:,:) = fft(fft(squeeze(SENSE(c,:,:)).*Data_rr,[],1),[],2).*kmask;
%end

%same but compact format
Data_crr=permute(repmat(Data_rr,[1 1 size(SENSE,1)]),[3,1,2]);
kmask_ckk=permute(repmat(kmask,[1 1 size(SENSE,1)]),[3,1,2]);
Data_ckk = fft(fft(SENSE.*Data_crr,[],2),[],3).*kmask_ckk;
end