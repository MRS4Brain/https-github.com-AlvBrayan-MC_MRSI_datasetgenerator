function Resized_raw_tckk=ResizeDataToRightSize(raw_tckk,NPhase,NRead);

[NbT NbC M N]=size(raw_tckk);
data_kk=squeeze(sum(sum(raw_tckk,1),2));

norm_k1=sum(abs(data_kk),2);
min_k1=min(find([norm_k1>0]));
max_k1=max(find([norm_k1>0]));

norm_k2=sum(abs(data_kk),1);
min_k2=min(find([norm_k2>0]));
max_k2=max(find([norm_k2>0]));

raw_tckk=raw_tckk(:,:,min_k1:max_k1,min_k2:max_k2);
data_kk=squeeze(sum(sum(abs(raw_tckk),1),2));
[X,Y]=meshgrid(1:M,1:N);
CMX=sum(data_kk.*X)/sum(data_kk(:);
CMY=sum(data_kk.*Y)/sum(data_kk(:);

Resized_raw_tckk=zeros([NbT, NbC,NRead,NPhase]);

end