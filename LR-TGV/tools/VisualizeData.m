function VisualizeData( mrsiData, HeadWater ,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Visualization_RawData.ps',name_fig);
if exist(s);delete(s);end
Size_data=size(mrsiData);

HWater_ctrr = ifft(ifft(HeadWater,[],3),[],4);
%[Uorig,Sorig,Vorig] = svd(reshape(permute(HWater_ctrr,[1 3 4 2]),[],size(HWater_ctrr,2)),0);
%USHead_crrc=reshape(Uorig*Sorig, Size_data(1),Size_data(3),Size_data(4),[]);

WaterAmp_crr=squeeze(mean(abs(fft(HWater_ctrr,[],2)),2));
%WaterPh_crr=squeeze(angle(USHead_crrc(:,:,:,1)));
WaterPh_crr=squeeze(angle(sum(HWater_ctrr(:,1:5,:,:),2)));
Head_crr=WaterAmp_crr.*exp(1j*WaterPh_crr);
if size(HWater_ctrr,1)==1
  Head_crr= reshape(Head_crr,[1 size(Head_crr)]) ;
end

MRSIData_ctrr = ifft(ifft(mrsiData,[],3),[],4);
%[Uorig,Sorig,Vorig] = svd(reshape(permute(MRSIData_ctrr,[1 3 4 2]),[],Size_data(2)),0);
%USData_crrc=reshape(Uorig*Sorig, Size_data(1),Size_data(3),Size_data(4),[]);

WaterAmp_crr=squeeze(mean(abs(fft(MRSIData_ctrr,[],2)),2));
%WaterPh_crr=squeeze(angle(USData_crrc(:,:,:,1)));
WaterPh_crr=squeeze(angle(sum(MRSIData_ctrr(:,1:5,:,:),2)));
Data_crr=WaterAmp_crr.*exp(1j*WaterPh_crr);
if size(MRSIData_ctrr,1)==1
  Data_crr= reshape(Data_crr,[1 size(Data_crr)]) ;
end

for C=1:Size_data(1);
   figs(C)=figure('visible', 'off'); 
        subplot(2,2,1),imagesc(abs(squeeze(Data_crr(C,:,:)))),colorbar;  
        axis('off');
        colormap default ;title('MRSI Data')
        subplot(2,2,3),imagesc(angle(squeeze(Data_crr(C,:,:)))),colorbar;  
        axis('off');
        colormap default 
        subplot(2,2,2),imagesc(abs(squeeze(Head_crr(C,:,:)))),colorbar;  
        axis('off');
        colormap default; title('Head Water Data')
        subplot(2,2,4),imagesc(angle(squeeze(Head_crr(C,:,:)))),colorbar;  
        axis('off');
        colormap default 
        print(figs(C), '-append', '-dpsc2', s); 
    end;
close all;




end

