function Data_ckkc = ForwSkullProjOp(Data_rrc,SkullMask,SENSE)
imDims=size(Data_rrc);
nDims = ndims(Data_rrc);
Data_kkc=zeros(imDims);
Data_crrc=permute(repmat(Data_rrc,[1 1 1 size(SENSE,1)]),[4,1,2,3]);
SENSE_crrc= repmat(SENSE,[1 1 1 imDims(end)]);
SkullMask_crrc=permute(repmat(SkullMask,[1 1 imDims(end) size(SENSE,1)]),[4,1,2,3]);
%SkullMask_rrc=repmat(SkullMask,[1 1 imDims(end)]);
Data_ckkc = fft(fft(SkullMask_crrc.*SENSE_crrc.*Data_crrc,[],2),[],3);

end