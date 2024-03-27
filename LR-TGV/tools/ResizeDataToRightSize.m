function Resized_raw_tckk=ResizeDataToRightSize(raw_tckk,NPhase,NRead);

[NbT NbC M N]=size(raw_tckk);

Resized_raw_tckk=zeros([NbT, NbC,2*M,2*N]);
Resized_raw_tckk(1:NbT, 1:NbC,1:M,1:N)=raw_tckk;

data_kk=squeeze(sum(sum(abs(Resized_raw_tckk).^2,1),2));
maxValue=max(data_kk(:));
[CM_r CM_p] = find(data_kk == maxValue);

P_CM_r=log2(CM_r-1);
P_CM_p=log2(CM_p-1);
norm_conv_cm=abs(P_CM_r-round(P_CM_r))/(P_CM_r) +abs(P_CM_p-round(P_CM_p))/(P_CM_p);
if(norm_conv_cm > 1e-4)
	warning('Original k-space center of mass (2^n+1 x 2^m +1 expected) is unusual in ResizeDataToRightSize.m!');
        fprintf(['Original center of mass is:[',num2str(CM_r),' ',num2str(CM_p),'].\n']);
        CM_r=2^round(P_CM_r)+1; CM_p=2^round(P_CM_p)+1;
        fprintf(['Center of mass was corrected to:[',num2str(CM_r),' ',num2str(CM_p),'].\n']);  
end

GoalCM_r=round((NRead +1)/2);
GoalCM_p=round((NPhase+1)/2);
shift_r=(CM_r-GoalCM_r);
shift_p=(CM_p-GoalCM_p);

Resized_raw_tckk=Resized_raw_tckk(:,:,(1+shift_r):(NRead+shift_r),(1+shift_p):(NPhase+shift_p));




%norm_k1=sum(abs(data_kk),2);
%min_k1=min(find([norm_k1>0]))
%max_k1=max(find([norm_k1>0]))

%norm_k2=sum(abs(data_kk),1);
%min_k2=min(find([norm_k2>0]))
%max_k2=max(find([norm_k2>0]))

%raw_tckk=raw_tckk(:,:,min_k1:max_k1,min_k2:max_k2);
%[NbT NbC M N]=size(raw_tckk)
%data_kk=squeeze(sum(sum(abs(raw_tckk),1),2));
%[Nr,Np]=meshgrid(1:N,1:M);
%CM_r=round(sum(sum(data_kk.*Nr,1),2)/sum(data_kk(:)))
%CM_p=round(sum(sum(data_kk.*Np,1),2)/sum(data_kk(:)))

%GoalCM_r=round((NRead +1)/2)
%GoalCM_p=round((NPhase+1)/2)
%Resized_raw_tckk=zeros([NbT, NbC,NRead,NPhase]);

%shift_r=(GoalCM_r-CM_r)
%shift_p=(GoalCM_p-CM_p)
%Resized_raw_tckk(:,:,(1+shift_r):(M+shift_r),(1+shift_p):(N+shift_p))=raw_tckk;
end
