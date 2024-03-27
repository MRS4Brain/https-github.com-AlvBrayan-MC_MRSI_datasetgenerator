function [SENSE, WaterFreqMap_rr,Water_trr,Water_rr]=ComputeSENSEProfilesWithESPIRiT(mrsiReconParams)
HomogenizeWithWaterSignal=0;

% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.02;
% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;
ksize = mrsiReconParams.ESPIRIT_kernel;%[5,5]; [6,6]; % kernel size

MSize_data=size(mrsiReconParams.mrsiData);
WSize_data=size(mrsiReconParams.Water_ctkk);

RowSet=[1:WSize_data(3)/2, (MSize_data(3)-round(WSize_data(3)/2)+1):MSize_data(3)];
 ColSet= [1:WSize_data(4)/2, (MSize_data(4)-round(WSize_data(4)/2)+1):MSize_data(4)];

nDimsOri= ndims(mrsiReconParams.Water_ctkk);

DATA=permute(sum(mrsiReconParams.Water_ctkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2),[3,4,1,2]);
DATA=fft(fft(fftshift(fftshift(ifft(ifft(DATA,[],1),[],2),1),2),[],1),[],2);
DATA=fftshift(fftshift(DATA,1),2);
Check_kmask=fftshift(fftshift(mrsiReconParams.Water_kmask,1),2);

[Nc,~,sx,sy] = size(mrsiReconParams.mrsiData);

%ncalib = round([sx/3,sy/3]);%round(min(sx,sy)/sqrt(2.0)/2); % FOR ELLIPTICAL ENCODING
%{
 ncalib=0;
 filling_callib=1;
while filling_callib==1
    ncalib= ncalib+1;
    filling_callib = mean(vectorizeArray( crop(Check_kmask,[ncalib,ncalib])));
end
ncalib=ncalib-1
%}






%% Compute ESPIRiT EigenVectors
% Here we perform calibration in k-space followed by an eigen-decomposition
% in image space to produce the EigenMaps. 


% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
% crop a calibration area
%calib = crop(DATA,[ncalib,ncalib,Nc]);
%[k,S] = dat2Kernel(calib,ksize);
% explicit kernel calculation with more extended calibration area (not only square center of kspace)

count=0;A=[];
%Used_Kernel=0*Check_kmask;
for y=1:(size(DATA,2)-ksize(2)+1)
    for x=1:(size(DATA,1)-ksize(1)+1)
        Check_Kernel=Check_kmask(x:(x+ksize(1)-1),y:(y+ksize(2)-1) ) ;
        if ( mean(Check_Kernel(:))==1) %then kernel is acquired at this location
            count = count+1;
            %Used_Kernel(x:(x+ksize(1)-1),y:(y+ksize(2)-1))=1;
            A(count,:,:) = reshape(DATA(x:(x+ksize(1)-1),y:(y+ksize(2)-1),:),1,ksize(1)*ksize(2),Nc);
        end
    end
end
%SA=size(A)
fprintf([num2str(count), ' kernels of size ',num2str(ksize(1)),'x',num2str(ksize(2)) , ' found in the Water measurement data.\n']);

A = reshape(A,size(A,1),size(A,2)*size(A,3));
[~,S,V] = svd(A);%,'econ');

k = reshape(V,ksize(1),ksize(2),Nc,size(V,2));
S = diag(S);S = S(:);
clear V
%


idx = max(find(S >= S(1)*eigThresh_1));


%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

%%
% crop sensitivity maps 
%maps = M(:,:,:,:,end).*repmat(W(:,:,:,end)>eigThresh_2,[1,1,1,Nc]);
maps = M(:,:,:,end).*repmat(mrsiReconParams.ImMask2D,[1,1,Nc]);

% figure, imshow3(abs(maps),[],[4,ceil(Nc/4)]); 
% title('Absolute sensitivity maps');
% colormap((gray(256))); colorbar;
% 
% figure, imshow3(angle (maps),[],[4,ceil(Nc/4)]); 
% title('Phase of sensitivity maps');
% colormap((jet(256))); colorbar;
SENSE=permute(maps,[3,1,2]);
%%


figs=figure('visible', 'off');

%for k=1:size(SENSE,1) 
 %   SENSE(k,:,:)=squeeze(SENSE(k,:,:)).*mrsiReconParams.ImMask2D;%+~mrsiReconParams.BrainMask*1E-6);  
%end
SENSE(find(isnan(SENSE)))=0;
clear Water_crrr exponent


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
%{
NX=MSize_data(3);NY=MSize_data(4);
sigma=mrsiReconParams.GaussianSigma;; 
[X,Y] = ndgrid(1:NX, 1:NY);
xc=floor(NX/2)+1;yc=floor(NY/2)+1;
exponent = -((X-xc).^2 + (Y-yc).^2)./(2*sigma^2);
%Kernel = exp(exponent)/sum(exp(exponent(:))); % no need to normalize

Kernel = fftshift(fftshift(exp(exponent)/sum(exp(exponent(:))),1),2); % no need to normalize
Kernel=fft(fft(Kernel,[],1),[],2);

%}
%Hanning Kernel
%{
[X,Y] = ndgrid(1:MSize_data(3), 1:MSize_data(4));
xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;
temp = (2*(X-xc)/WSize_data(3)).^2 + (2*(Y-yc)/WSize_data(4)).^2;
HKernel = fftshift(0.5*(1+cos(pi*temp))); 
%}

LowResWater_ctrr=ifft(ifft(mrsiReconParams.Water_ctkk,[],3),[],4);
% Zero Padding Water and filtering
%{
WaterZeroPad_ctkk=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4));
WaterZeroPad_ctkk(:,:,RowSet,ColSet)=mrsiReconParams.Water_ctkk;


WaterZPadKmask=zeros(MSize_data(3),MSize_data(4));
WaterZPadKmask(RowSet,ColSet)=mrsiReconParams.Water_kmask;

HKernel_ctkk=permute(repmat(HKernel,[1,1,WSize_data(1),WSize_data(2)]),[3 4 1 2]);
WaterZeroPad_ctkk = WaterZeroPad_ctkk.*HKernel_ctkk;
%}
%{
WaterAmp_ckk=squeeze(mean(abs(WaterZeroPad_ctkk (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:)),2));
WaterPh_ckk=squeeze(angle(sum(WaterZeroPad_ctkk (:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2)));
Water_ckk=WaterAmp_ckk.*exp(1j*WaterPh_ckk);
%}

clear WaterAmp_ckk WaterPh_ckk

if HomogenizeWithWaterSignal
Water_ckk=squeeze(mean(WaterZeroPad_ctkk(:,1:mrsiReconParams.NbPtForWaterPhAmp,:,:),2));
FullkSpace=ones(MSize_data(3),MSize_data(4));
%skmask=ones(MSize_data(3),MSize_data(4));
WaterHomo=1;
WH_it=1;
Threshold=1E-8;
RadiusErode=1;%2;
MaskWaterAdjust=1+imerode(round(mrsiReconParams.ImMask2D),offsetstrel(ones(2*RadiusErode+1)));


while WaterHomo>1E-2 & WH_it<5 
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
end




%{
[Uorig,Sorig,Vorig] = svd(reshape(permute(WaterZeroPad_ctkk,[1 3 4 2]),[],WSize_data(2)),0);

if size(WaterZeroPad_ctkk,2)<mrsiReconParams.modelOrder
	WaterOrder=size(WaterZeroPad_ctkk,2);
else
	WaterOrder=mrsiReconParams.modelOrder;
end

V=Vorig(:,1:WaterOrder);
U_ckkc=reshape(Uorig, WSize_data(1),MSize_data(3),MSize_data(4),[]);
Recon_US_rrc=zeros(MSize_data(3),MSize_data(4),WaterOrder);
%alpha = 5E-3;
alpha = 5E-6;
parfor k=1:WaterOrder
        [Recon_US_rrc(:,:,k) e] = tgv2_l2_2D_multiCoil(Sorig(k,k)*U_ckkc(:,:,:,k),SENSE, WaterZPadKmask, 2*alpha, alpha, maxit,Threshold);
        %Recon_US_rrc(:,:,k)=Sorig(k,k)*Recon_US_rrc(:,:,k);
end

 Water_trr=formTensorProduct(Recon_US_rrc,V,2);
 Water_trr=permute(Water_trr,[3,1,2]);
%} 

%MSize_data
%WSize_data
[Xw,Yw] = meshgrid(linspace(MSize_data(4)/WSize_data(4),MSize_data(4),WSize_data(4))-0.5*MSize_data(4)/WSize_data(4),linspace(MSize_data(3)/WSize_data(3),MSize_data(3),WSize_data(3))-0.5*MSize_data(3)/WSize_data(3));
[Xm,Ym] = meshgrid((1:MSize_data(4))-0.5,(1:MSize_data(3))-0.5);

Water_ctrr=zeros(MSize_data);
for c=1:WSize_data(1) 
    for t=1:WSize_data(2)  
        Water_ctrr(c,t,:,:) = interp2(Xw,Yw,squeeze(LowResWater_ctrr(c,t,:,:)),Xm,Ym,'spline');
    end
end
SENSE_c1rr=reshape(SENSE,[MSize_data(1) 1  MSize_data(3)  MSize_data(4)]  ); %sum(abs(SENSE).^2)=1!
Water_trr=squeeze(sum(Water_ctrr.*conj(SENSE_c1rr),1)).*reshape(mrsiReconParams.ImMask2D,[1 MSize_data(3)  MSize_data(4)]);
Water_rr=squeeze(mean(Water_trr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:),1));
WaterFreqMap_rr=DetermineSingleLowFreq( Water_trr,60,mrsiReconParams);
WaterFreqMap_rr=mrsiReconParams.BrainMask2D.*WaterFreqMap_rr-mean(WaterFreqMap_rr(mrsiReconParams.BrainMask2D>0));%remove intercept and set Freq to 0 outside Brain


close all;figs=figure('visible', 'off');
Data_name=mrsiReconParams.NameData;
s1=['./',mrsiReconParams.Log_Dir,'/','WHB_ESPIRIT_CoilProfiles_FromWaterMeas_Abs_',Data_name,'.ps'];
s2=['./',mrsiReconParams.Log_Dir,'/','WHB_ESPIRIT_CoilProfiles_FromWaterMeas_Phase_',Data_name,'.ps'];

if exist(s1);delete(s1);end
if exist(s2);delete(s2);end

for k=1:WSize_data(1)  
    imagesc(abs(squeeze(SENSE(k,:,:))));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s1);
    imagesc(angle(squeeze(SENSE(k,:,:))));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s2);
    
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
