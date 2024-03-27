function [SENSE, WaterFreqMap_rr,Water_trr,Water_rr]=ComputeSENSEProfilesFromHeadBodyWaterAcq(mrsiReconParams)

%Water_crr=squeeze(sum(ifft(ifft(mrsiReconParams.Water_ctkk(:,1:round(end/10),:,:),[],3),[],4),2));

MSize_data=size(mrsiReconParams.mrsiData);
BWSize_data=size(mrsiReconParams.BodyWater_ctkk);
WSize_data=size(mrsiReconParams.HeadWater_ctkk);
if(sum(BWSize_data(2:end)==WSize_data(2:end))<(size(BWSize_data)-1))
   error('Size of Body Water Data differs from Head Water Data!') ;
end
nDimsOri= ndims(mrsiReconParams.BodyWater_ctkk);

%}
% *************************************************************************
% Compute Profile in Original Resolution
% *************************************************************************

HWater_ctrr = ifft(ifft(mrsiReconParams.HeadWater_ctkk,[],3),[],4);

[Uorig,Sorig,Vorig] = svd(reshape(permute(HWater_ctrr,[1 3 4 2]),[],WSize_data(2)),0);
USHead_crrc=reshape(Uorig*Sorig, WSize_data(1),WSize_data(3),WSize_data(4),[]);
WaterAmpOrig_crr=abs(USHead_crrc(:,:,:,1));
WaterPhOrig_crr=angle(USHead_crrc(:,:,:,1));




% *************************************************************************
% Compute Profile in High Resolution
% *************************************************************************

BWater_ctrr = ifft(ifft(mrsiReconParams.BodyWater_ctkk,[],3),[],4);
%BWater_ctrr = ifft(ifft(mrsiReconParams.HeadWater_ctkk,[],3),[],4);

WaterAmp_crr=squeeze(mean(abs(BWater_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
WaterPh_crr=squeeze(angle(sum(BWater_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));
USBody_crr=WaterAmp_crr.*exp(1j*WaterPh_crr);
BodyNorm=norm(USBody_crr(:));
if ndims(USBody_crr)==3
    SENSESQ = squeeze(sqrt(sum(abs(USBody_crr).^2,1)));
else
    SENSESQ = abs(USBody_crr); 
end
%%

Data_name=mrsiReconParams.NameData;
s1=['./',mrsiReconParams.Log_Dir,'/','WHB_SENSECoilProfiles_FromWaterMeas_Abs_',Data_name,'.ps'];
s2=['./',mrsiReconParams.Log_Dir,'/','WHB_SENSECoilProfiles_FromWaterMeas_Phase_',Data_name,'.ps'];
s3=['./',mrsiReconParams.Log_Dir,'/','WHB_HeadWater_perCoil_Abs_',Data_name,'.ps'];
s4=['./',mrsiReconParams.Log_Dir,'/','WHB_HeadWater_perCoil_Phase_',Data_name,'.ps'];
if exist(s1);delete(s1);end
if exist(s2);delete(s2);end
if exist(s3);delete(s3);end
if exist(s4);delete(s4);end


WaterAmp_crr=squeeze(mean(abs(HWater_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
WaterPh_crr=squeeze(angle(sum(HWater_ctrr (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));
USHead_crrc=WaterAmp_crr.*exp(1j*WaterPh_crr);
HeadNorm=norm(USHead_crrc(:));
SENSE=BodyNorm/HeadNorm*USHead_crrc(:,:,:,1);

%SignalThres=mean(SENSESQ(mrsiReconParams.ImMask2D==0))*2;
SignalThres=quantile(SENSESQ(mrsiReconParams.ImMask2D==0),0.95)*2;
if isnan(SignalThres)
    SignalThres=quantile(SENSESQ(:),0.50);
end
for k=1:size(mrsiReconParams.HeadWater_ctkk,1)   
  SENSE(k,:,:)=squeeze(SENSE(k,:,:))./SENSESQ.*(SENSESQ>SignalThres);
  SENSE(k,:,:)=squeeze(SENSE(k,:,:)).*mrsiReconParams.BrainMask2D;
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

%WaterZeroPad_ckkt=fft(fft(permute(WaterZeroPad_ctrr,[1 3 4 2]),[],2),[],3);
%[Uorig,Sorig,Vorig] = svd(reshape(WaterZeroPad_ckkt,[],WSize_data(2)),0);
%V=Vorig(:,1:mrsiReconParams.modelOrder);
%U_ckkc=reshape(Uorig, MSize_data(1),MSize_data(3),MSize_data(4),[]);
%Recon_US_rrc=zeros(MSize_data(3),MSize_data(4),20);
skmask=ones(WSize_data(3),WSize_data(4));
[Uorig,Sorig,Vorig] = svd(reshape(permute(mrsiReconParams.HeadWater_ctkk,[1 3 4 2]),[],WSize_data(2)),0);
V=Vorig(:,1:mrsiReconParams.modelOrder);
U_ckkc=reshape(Uorig, WSize_data(1),WSize_data(3),WSize_data(4),[]);

Recon_US_rrc=zeros(WSize_data(3),WSize_data(4),mrsiReconParams.modelOrder);
Threshold=1E-8;
WaterHomo=1;
%smoothing kernel
NX=WSize_data(3);NY=WSize_data(4);
sigma=mrsiReconParams.GaussianSigma;;%1 was not enough
[X,Y] = ndgrid(1:NX, 1:NY);
xc=floor(NX/2)+1;yc=floor(NY/2)+1;
exponent = -((X-xc).^2 + (Y-yc).^2)./(2*sigma^2);
Kernel = exp(exponent)/sum(exp(exponent(:))); % no need to normalize
    

while WaterHomo>1E-3;
    parfor k=1:mrsiReconParams.modelOrder
        [Recon_US_rrc(:,:,k) e] = tgv2_l2_2D_multiCoil(U_ckkc(:,:,:,k),SENSE, skmask, 2*alpha, alpha, maxit,Threshold);
        Recon_US_rrc(:,:,k)=Sorig(k,k)*Recon_US_rrc(:,:,k);
    end
    %clear U_ckkc Uorig
    
    Water_trr=formTensorProduct(Recon_US_rrc,V,2);
    Water_trr=permute(Water_trr,[3,1,2]);   
    %Water_rr=(Recon_US_rrc(:,:,1));
    
    WaterAmp_rr=squeeze(mean(abs(Water_trr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),1));
    WaterPh_rr=squeeze(angle(sum(Water_trr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:),1)));
    Water_rr=WaterAmp_rr.*exp(1j*WaterPh_rr);
    MeanWater=mean(abs(Water_rr(mrsiReconParams.BrainMask2D>0)));
    Water_rr=mrsiReconParams.BrainMask2D.*abs(Water_rr)/MeanWater+(~mrsiReconParams.BrainMask2D);
 
    %smoothing of Maps
    
    SmoWater_rr=conv2(Water_rr,Kernel,'same');
    SmoWater_rr=mrsiReconParams.BrainMask2D.*SmoWater_rr+(~mrsiReconParams.BrainMask2D);
    
    for k=1:size(mrsiReconParams.HeadWater_ctkk,1)
        SENSE(k,:,:)=squeeze(SENSE(k,:,:)).*SmoWater_rr; %No correction for Synthetic Data 
    end
    WaterHomo=std(SmoWater_rr(:))/mean(SmoWater_rr(:)) %0 for Synthetic Data 
end

%RoemersCoefs = mrsiReconParams.RoemersCoefs;% ComputeRoemersProfiles(mrsiReconParams);
%RoemersCoefs_ctrr=permute(repmat(RoemersCoefs,[1 1 1 size(mrsiReconParams.Water_ctkk,2)]),[1,4,2,3]);

%WaterZeroPad_trr= conj(RoemersCoefs_ctrr).*WaterZeroPad_ctrr;
%WaterZeroPad_trr =squeeze(sum(WaterZeroPad_trr,1));
%clear RoemersCoefs_ctrr RoemersCoefs


WaterFreqMap_rr=DetermineSingleLowFreq( Water_trr,60,mrsiReconParams);
WaterFreqMap_rr=WaterFreqMap_rr-mean(WaterFreqMap_rr(mrsiReconParams.BrainMask2D>0));%remove intercept

figs=figure('visible', 'off');

for k=1:size(mrsiReconParams.HeadWater_ctkk,1)  
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
