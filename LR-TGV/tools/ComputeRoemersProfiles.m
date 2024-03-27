function SENSE = ComputeRoemersProfiles(mrsiReconParams,Water_ctkk)
reduction = 2^(-8);     % usually there is no need to change this
alpha = 10;%mrsiReconParams.UndersamplingF*mrsiReconParams.mu_tv/reduction; % usually there is no need to change this
maxit = 500;            % use 1000 Iterations for optimal image quality
innerIter = 20;         % usually there is no need to change this
%kmask = mrsiReconParams.kmask ;
%Water_crr=squeeze(sum(ifft(ifft(mrsiReconParams.Water_ctkk(:,1:round(end/10),:,:),[],3),[],4),2));

MSize_data=size(mrsiReconParams.mrsiData);
WSize_data=size(Water_ctkk);

 RowSet=[1:WSize_data(3)/2, (MSize_data(3)-round(WSize_data(3)/2)+1):MSize_data(3)];
 ColSet= [1:WSize_data(4)/2, (MSize_data(4)-round(WSize_data(4)/2)+1):MSize_data(4)];
%kmask = 0*mrsiReconParams.kmask;
%kmask(RowSet,ColSet) = mrsiReconParams.kmask ;


%Gaussian Kernel
%sigma=size(mrsiReconParams.Water_ctkk,3)/4;
%sigma=size(Water_ctkk,3)/3; %Optimized with comparison to the Full Res Water Data
sigma=mrsiReconParams.GaussianSigma*2;

[X,Y] = ndgrid(1:MSize_data(3), 1:MSize_data(4));
xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;
exponent = -((X-xc).^2 + (Y-yc).^2)./(2*sigma^2);
%Kernel = fftshift(1 / (2 * sqrt(2*pi)).*exp(exponent));  
%Kernel = fftshift(fft2(exp(exponent))); % no need to normalize
Kernel = fft2(fftshift(exp(exponent))); % no need to normalize

% Zero Padding Water and filtering
WaterZeroPad_ctkk=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4));

WaterZeroPad_ctkk(:,:,RowSet,ColSet)=Water_ctkk;
Kernel_ctkk=permute(repmat(Kernel,[1,1,WSize_data(1),WSize_data(2)]),[3 4 1 2]);
WaterZeroPad_ctrr  = ifft(ifft(WaterZeroPad_ctkk.*Kernel_ctkk,[],3),[],4);
%WaterZeroPad_ctrr  = ifft(ifft(WaterZeroPad_ctkk,[],3),[],4);
%WaterAmp_crr=squeeze(sum(abs(fft(ifft(ifft(mrsiReconParams.Water_ctkk,[],3),[],4),[],2)),2));%Doesn't work with data with too sharp peak in frequency


%[Uorig,Sorig,Vorig] = svd(reshape(permute(WaterZeroPad_ctrr,[1 3 4 2]),[],WSize_data(2)),0);
%US_crrc=reshape(Uorig*Sorig, WSize_data(1),MSize_data(3),MSize_data(4),[]);
WaterAmp_crr=squeeze(mean(abs(WaterZeroPad_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
WaterPh_crr=squeeze(angle(sum(WaterZeroPad_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));

clear WaterZeroPad_ctrr Kernel_ctkk WaterZeroPad_ctkk
%%

for k=1:size(mrsiReconParams.Water_ctkk,1)
    %[USreal(:,:,k), funCost] = TGV2_Denoise_MultiMetabVol(im, mrsiReconParams.mu_tv*datanorm);
   SENSE(k,:,:)=exp(1i*squeeze(WaterPh_crr(k,:,:)-WaterPh_crr(mrsiReconParams.mrProt.Main_coil_element,:,:))).*squeeze(WaterAmp_crr(k,:,:)).^(0.5);%./squeeze(sum(WaterAmp_crr,1).^(0.5)); 
% SENSE(k,:,:)=exp(1i*squeeze(WaterPh_crr(k,:,:))).*squeeze(WaterAmp_crr(k,:,:)).^(0.5);%./squeeze(sum(WaterAmp_crr,1).^(0.5)); 

end

%SENSE=US_crrc(:,:,:,1);% 1st SVD component represents water at best


Data_name=mrsiReconParams.NameData;
% s1=['./',mrsiReconParams.Log_Dir,'/','RoemersProfiles_Abs_',Data_name,'.ps'];
% s2=['./',mrsiReconParams.Log_Dir,'/','RoemersProfiles_Phase_',Data_name,'.ps'];
% s3=['./',mrsiReconParams.Log_Dir,'/','RoemersOrigSignal_Abs_',Data_name,'.ps'];
% s4=['./',mrsiReconParams.Log_Dir,'/','RoemersOrigSignal_Phase_',Data_name,'.ps'];
% if exist(s1);delete(s1);end
% if exist(s2);delete(s2);end
% if exist(s3);delete(s3);end
% if exist(s4);delete(s4);end
% figs=figure('visible', 'off');

SENSESQ = sqrt(sum(abs(SENSE).^2,1));
for k=1:size(Water_ctkk,1)
    SENSE(k,:,:)=SENSE(k,:,:)./SENSESQ;
  

%         imagesc(abs(squeeze(SENSE(k,:,:))));%,[ 0, 10*mean(image2plot(:))] )
%         colormap default;colorbar;
%         print(figs, '-append', '-dpsc2', s1);
%         imagesc(angle(squeeze(SENSE(k,:,:))));%,[ 0, 10*mean(image2plot(:))] )
%         colormap default;colorbar;
%         print(figs, '-append', '-dpsc2', s2);
%         
%         
%         imagesc(abs(squeeze(US_crrc(k,:,:,1))));%,[ 0, 10*mean(image2plot(:))] )
%         colormap default;colorbar;
%         print(figs, '-append', '-dpsc2', s3);
%         imagesc(squeeze( angle(US_crrc(k,:,:,1))));%,[ 0, 10*mean(image2plot(:))] )
%         colormap default;colorbar;
%         print(figs, '-append', '-dpsc2', s4);
  
end 
SENSE(find(isnan(SENSE)))=0;
clear Water_crr Kernel exponent
close all;
end
