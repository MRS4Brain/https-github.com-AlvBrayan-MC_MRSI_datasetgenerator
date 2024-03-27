function Data_rr =  AdjEncodOp(Data_ckk, kmask,SENSE)

%Data_rr = zeros(size(kmask));
%for c=1:size(SENSE,1)
    %u=real(ifft(ifft(datakk.*kmask,[],1),[],2));
 %   Data_rr = Data_rr + conj(squeeze(SENSE(c,:,:))).*ifft(ifft(squeeze(Data_ckk(c,:,:)).*kmask,[],1),[],2); % Adjoint Solution
%end
%same but compact format
kmask_ckk=permute(repmat(kmask,[1 1 size(SENSE,1)]),[3,1,2]);
Data_rr = conj(SENSE).*ifft(ifft(Data_ckk.*kmask_ckk,[],2),[],3);
Data_rr=squeeze(sum(Data_rr,1));
end