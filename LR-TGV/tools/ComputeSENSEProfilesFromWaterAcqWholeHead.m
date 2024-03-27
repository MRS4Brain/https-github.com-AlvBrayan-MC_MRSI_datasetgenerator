function [SENSE, WaterFreqMap_rr,Water_trr,Water_rr]=ComputeSENSEProfilesFromWaterAcqWholeHead(mrsiReconParams)

%Water_crr=squeeze(sum(ifft(ifft(mrsiReconParams.Water_ctkk(:,1:round(end/10),:,:),[],3),[],4),2));

MSize_data=size(mrsiReconParams.mrsiData);
WSize_data=size(mrsiReconParams.Water_ctkk);


nDimsOri= ndims(mrsiReconParams.Water_ctkk);

%}
% *************************************************************************
% Compute Profile in Original Resolution
% *************************************************************************

Water_ctrr = ifft(ifft(mrsiReconParams.Water_ctkk,[],3),[],4);

[Uorig,Sorig,Vorig] = svd(reshape(permute(Water_ctrr,[1 3 4 2]),[],WSize_data(2)),0);
US_crrc=reshape(Uorig*Sorig, WSize_data(1),WSize_data(3),WSize_data(4),[]);
WaterAmpOrig_crr=abs(US_crrc(:,:,:,1));
WaterPhOrig_crr=angle(US_crrc(:,:,:,1));


% *************************************************************************
% Compute Profile in High Resolution
% *************************************************************************

RowSet=[1:WSize_data(3)/2, (MSize_data(3)-round(WSize_data(3)/2)+1):MSize_data(3)];
 ColSet= [1:WSize_data(4)/2, (MSize_data(4)-round(WSize_data(4)/2)+1):MSize_data(4)];


%Hanning Kernel
sigma=(size(mrsiReconParams.Water_ctkk,3)+size(mrsiReconParams.Water_ctkk,4))/6;%3; %Optimized with comparison to the Full Res Water Data
[X,Y] = ndgrid(1:MSize_data(3), 1:MSize_data(4));
xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;
temp = (2*(X-xc)/WSize_data(3)).^2 + (2*(Y-yc)/WSize_data(4)).^2;
HKernel = fftshift(0.5*(1+cos(pi*temp)));  

% Zero Padding Water and filtering
WaterZeroPad_ctkk=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4));
WaterZeroPad_ctkk(:,:,RowSet,ColSet)=mrsiReconParams.Water_ctkk;
HKernel_ctkk=permute(repmat(HKernel,[1,1,WSize_data(1),WSize_data(2)]),[3 4 1 2]);
WaterZeroPad_ctrr  = ifft(ifft(WaterZeroPad_ctkk.*HKernel_ctkk,[],3),[],4);


WaterAmp_crr=squeeze(mean(abs(WaterZeroPad_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
WaterPh_crr=squeeze(angle(sum(WaterZeroPad_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));
US_crr=WaterAmp_crr.*exp(1j*WaterPh_crr);

SENSE=US_crr;

if ndims(US_crr)==3
    SENSESQ = squeeze(sqrt(sum(abs(US_crr).^2,1)));
else
    SENSESQ = abs(US_crr); 
end

%%


%SignalThres=mean(SENSESQ(mrsiReconParams.ImMask2D==0))*2;
%SignalThres=quantile(SENSESQ(mrsiReconParams.ImMask2D==0),0.95)*2;
%if isnan(SignalThres)
    SignalThres=0.5*mean(SENSESQ(:));
%end
for k=1:WSize_data(1)   
  %SENSE(k,:,:)=squeeze(SENSE(k,:,:))./SENSESQ.*(SENSESQ>SignalThres);
  SENSE(k,:,:)=squeeze(SENSE(k,:,:))./SENSESQ;
  SENSE(k,:,:)=squeeze(SENSE(k,:,:)).*mrsiReconParams.ImMask2D;
  %SENSE(k,:,:)=squeeze(SENSE(k,:,:)).*mrsiReconParams.BrainMask2D;

end
SENSE(find(isnan(SENSE)))=0;
clear Water_crr exponent

%%%%% compute Freq Maps

% Simple FFT and Coil Combination
kmask=ones(MSize_data(3),MSize_data(4));
%reduction = 2^(-8);     % usually there is no need to change this
alpha = 1E-32; %no Regularization
maxit = 1000;            % use 1000 Iterations for optimal image quality

%for k=1:size(US_ckkc,4)
fprintf('Start the Water reconstruction...\n');


s5=['./',mrsiReconParams.Log_Dir,'/','Water_Homogenization_',mrsiReconParams.NameData,'.ps'];
if exist(s5);delete(s5);end
figs=figure('visible', 'off');

%smoothing of Maps
NX=MSize_data(3);NY=MSize_data(4);
sigma=mrsiReconParams.GaussianSigma;; 
[X,Y] = ndgrid(1:NX, 1:NY);
xc=floor(NX/2)+1;yc=floor(NY/2)+1;
exponent = -((X-xc).^2 + (Y-yc).^2)./(2*sigma^2);
%Kernel = exp(exponent)/sum(exp(exponent(:))); % no need to normalize

Kernel = fftshift(fftshift(exp(exponent)/sum(exp(exponent(:))),1),2); % no need to normalize
Kernel=fft(fft(Kernel,[],1),[],2);

%WaterAmp_ckk=squeeze(mean(abs(mrsiReconParams.Water_ctkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
%WaterPh_ckk=squeeze(angle(sum(mrsiReconParams.Water_ctkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));
WaterAmp_ckk=squeeze(mean(abs(WaterZeroPad_ctkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
WaterPh_ckk=squeeze(angle(sum(WaterZeroPad_ctkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));

Water_ckk=WaterAmp_ckk.*exp(1j*WaterPh_ckk);
clear WaterAmp_ckk WaterPh_ckk

skmask=ones(MSize_data(3),MSize_data(4));
WaterHomo=1;
WH_it=1;
Threshold=1E-8;
RadiusErode=1;%2;
MaskWaterAdjust=1+imerode(round(mrsiReconParams.ImMask2D),offsetstrel(ones(2*RadiusErode+1)));


while WaterHomo>1E-2 & WH_it<3 %25
    [Water_rr e] = tgv2_l2_2D_multiCoil(Water_ckk,SENSE, skmask, 2*alpha, alpha, maxit,Threshold);
  
    MeanWater=mean(abs(Water_rr(MaskWaterAdjust>0)));
    Water_rr=MaskWaterAdjust.*abs(Water_rr)/MeanWater+(~MaskWaterAdjust);
  
    %SmoWater_rr=conv2(Water_rr,Kernel,'same');
    SmoWater_rr=abs(ifft(ifft(Kernel.*fft(fft(Water_rr,[],1),[],2),[],1),[],2));
    %SmoWater_rr=SmoWater_rr.*mrsiReconParams.ImMask2D+(~mrsiReconParams.ImMask2D);
    
    subplot(1,2,1);
    imagesc(  Water_rr);%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;title(['Water Signal It:',num2str(WH_it)]);
    subplot(1,2,2);
    imagesc( SmoWater_rr);%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;title(['Water Signal Smoothed It:',num2str(WH_it)]);
    print(figs,'-bestfit', '-append', '-dpsc2', s5);
    
    for k=1:WSize_data(1)
        SENSE(k,:,:)=squeeze(SENSE(k,:,:)).*SmoWater_rr; %No correction for Synthetic Data 
    end
    WaterHomo=std(SmoWater_rr(:))/mean(SmoWater_rr(:)) %0 for Synthetic Data 
    WH_it=WH_it+1
end


for k=1:WSize_data(1)   

  SENSE(k,:,:)=squeeze(SENSE(k,:,:)).*mrsiReconParams.BrainMask2D;

end

%[Uorig,Sorig,Vorig] = svd(reshape(permute(mrsiReconParams.Water_ctkk,[1 3 4 2]),[],WSize_data(2)),0);
[Uorig,Sorig,Vorig] = svd(reshape(permute(WaterZeroPad_ctkk.*HKernel_ctkk,[1 3 4 2]),[],WSize_data(2)),0);
V=Vorig(:,1:mrsiReconParams.modelOrder);
%U_ckkc=reshape(Uorig, WSize_data(1),WSize_data(3),WSize_data(4),[]);
U_ckkc=reshape(Uorig, MSize_data(1),MSize_data(3),MSize_data(4),[]);

Recon_US_rrc=zeros(MSize_data(3),MSize_data(4),mrsiReconParams.modelOrder);
 
parfor k=1:mrsiReconParams.modelOrder
        [Recon_US_rrc(:,:,k) e] = tgv2_l2_2D_multiCoil(U_ckkc(:,:,:,k),SENSE, skmask, 2*alpha, alpha, maxit,Threshold);
        Recon_US_rrc(:,:,k)=Sorig(k,k)*Recon_US_rrc(:,:,k);
end
    
 Water_trr=formTensorProduct(Recon_US_rrc,V,2);
 Water_trr=permute(Water_trr,[3,1,2]);
 
 WaterAmp_rr=squeeze(mean(abs(Water_trr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),1));
 WaterPh_rr=squeeze(angle(sum(Water_trr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:),1)));
 Water_rr=WaterAmp_rr.*exp(1j*WaterPh_rr);

WaterFreqMap_rr=DetermineSingleLowFreq( Water_trr,60,mrsiReconParams);
WaterFreqMap_rr=WaterFreqMap_rr-mean(WaterFreqMap_rr(mrsiReconParams.BrainMask2D>0));%remove intercept


close all;figs=figure('visible', 'off');
Data_name=mrsiReconParams.NameData;
s1=['./',mrsiReconParams.Log_Dir,'/','WHB_SENSECoilProfiles_FromWaterMeas_Abs_',Data_name,'.ps'];
s2=['./',mrsiReconParams.Log_Dir,'/','WHB_SENSECoilProfiles_FromWaterMeas_Phase_',Data_name,'.ps'];
s3=['./',mrsiReconParams.Log_Dir,'/','WHB_HeadWater_perCoil_Abs_',Data_name,'.ps'];
s4=['./',mrsiReconParams.Log_Dir,'/','WHB_HeadWater_perCoil_Phase_',Data_name,'.ps'];
if exist(s1);delete(s1);end
if exist(s2);delete(s2);end
if exist(s3);delete(s3);end
if exist(s4);delete(s4);end

for k=1:WSize_data(1)  
    imagesc(abs(squeeze(SENSE(k,:,:))));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s1);
    imagesc(angle(squeeze(SENSE(k,:,:))));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s2);
    
    
    imagesc(abs(squeeze(WaterAmpOrig_crr(k,:,:))));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s3);
    imagesc(squeeze( WaterPhOrig_crr(k,:,:)));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s4);
end


s=['./',mrsiReconParams.Log_Dir,'/',mrsiReconParams.NameData, '_WB_Water_Amp_Phase_FreqMap.ps'];
if exist(s);delete(s);end
figs=figure('visible', 'off');

imagesc(abs(Water_rr),[0 max(abs(Water_rr(:)))]);%,[ 0, 10*mean(image2plot(:))] )
colormap default;colorbar;title('Water Signal Amplitude')
print(figs, '-append', '-dpsc2', s);

imagesc(angle(Water_rr));%,[ 0, 10*mean(image2plot(:))] )
colormap default;colorbar;title('Water Signal Phase')
print(figs, '-append', '-dpsc2', s);

imagesc(WaterFreqMap_rr);%,[ 0, 10*mean(image2plot(:))] )
colormap default;colorbar;title('Freq Map Based on Water Signal')
print(figs, '-append', '-dpsc2', s);

%imagesc(mrsiReconParams.b0map2D,[min(WaterFreqMap_rr(:)) max(WaterFreqMap_rr(:))] );%,[ 0, 10*mean(image2plot(:))] )
%colormap default;colorbar;title('Measured B0 Field map')
%print(figs, '-append', '-dpsc2', s);

close all;
fprintf('Coil profile computation and water reconstruction finished.\n');
end
