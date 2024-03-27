function VisualizeCombinedData( mrsiData, Water ,mrsiReconParams,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Visualization_CombinedRawData.ps',name_fig);
if exist(s);delete(s);end

Size_data=size(mrsiData);

Sdata=size(Water);
SMask=size(squeeze(mrsiReconParams.SENSE(1,:,:)));
    
[Xm,Ym] = meshgrid(linspace(1,Sdata(4),SMask(2)),linspace(1,Sdata(3),SMask(1)));
[Xd,Yd] = meshgrid(1:Sdata(4),1:Sdata(3));
%for c=1:size(mrsiReconParams.SENSE,1)
%	WSENSE(c,:,:,:) = interp3(Xm,Ym,Zm,squeeze(mrsiReconParams.SENSE(c,:,:,:)),Xd,Yd,Zd,'linear');
%end
%WImMask=interp3(Xm,Ym,Zm,mrsiReconParams.ImMask,Xd,Yd,Zd,'nearest');

SENSE=mrsiReconParams.SENSE./sum(abs(mrsiReconParams.SENSE).^2,1);
%WSENSE=WSENSE./sum(abs(WSENSE).^2,1);


Water_ctrr = ifft(ifft(Water,[],3),[],4);

%WaterAmp_crrr=squeeze(mean(abs(fft(HWater_ctrrr,[],2)),2));
%WaterPh_crrr=squeeze(angle(sum(HWater_ctrrr(:,1:5,:,:,:),2)));
%Head_crrr=WaterAmp_crrr.*exp(1j*WaterPh_crrr);
Water_crr=squeeze(mean(Water_ctrr(:,1:5,:,:),2));
if size(Water_ctrr,1)==1
  Water_crr= reshape(Water_crr,[1 size(Water_crr)]) ;
end
HRWater_crr=[];
for c=1:size(mrsiReconParams.SENSE,1)
	HRWater_crr(c,:,:) = interp2(Xd,Yd,squeeze(Water_crr(c,:,:)),Xm,Ym,'spline');
end

%{

MSize_data=size(mrsiData);
WSize_data=size(Water);

RowSet=[1:WSize_data(3)/2, (MSize_data(3)-round(WSize_data(3)/2)+1):MSize_data(3)];
ColSet= [1:WSize_data(4)/2, (MSize_data(4)-round(WSize_data(4)/2)+1):MSize_data(4)];

[X,Y] = ndgrid(1:MSize_data(3), 1:MSize_data(4));
xc=floor(MSize_data(3)/2)+1;yc=floor(MSize_data(4)/2)+1;
temp = (2*(X-xc)/WSize_data(3)).^2 + (2*(Y-yc)/WSize_data(4)).^2;
HKernel = fftshift(0.5*(1+cos(pi*temp))); 

% Zero Padding Water and filtering
WaterZeroPad_ctkk=zeros(WSize_data(1),WSize_data(2), MSize_data(3),MSize_data(4));
WaterZeroPad_ctkk(:,:,RowSet,ColSet)=Water;

HRWater_crr=squeeze(mean(ifft(ifft(WaterZeroPad_ctkk(:,1:5,:,:),[],3),[],4),2));
%}

CombWater_rr=squeeze(sum(HRWater_crr.*conj(SENSE),1)).*mrsiReconParams.ImMask2D;
clear HWater_ctrrr  WaterAmp_crrr WaterPh_crrr
MRSIData_ctrr =ifft(ifft(mrsiData,[],3),[],4);

%WaterAmp_crrr=squeeze(mean(abs(fft(MRSIData_ctrrr,[],2)),2));
%WaterPh_crrr=squeeze(angle(sum(MRSIData_ctrrr(:,1:5,:,:,:),2)));
%Data_crrr=WaterAmp_crrr.*exp(1j*WaterPh_crrr);
Data_crr=squeeze(mean(MRSIData_ctrr(:,1:5,:,:,:),2));
if size(MRSIData_ctrr,1)==1
  Data_crr= reshape(Data_crr,[1 size(Data_crr)]) ;
end
CombData_rr=squeeze(sum(Data_crr.*conj(SENSE),1)).*mrsiReconParams.ImMask2D;

clear MRSIData_ctrrr WaterAmp_crrr  WaterPh_crrr


   figs=figure('visible', 'off'); 
         
       image2plot=abs(CombData_rr);
       imagesc(abs(image2plot));
       title('MRSI Data Amp.'); axis('off');colormap default ;
       print(figs, '-append','-bestfit', '-dpsc2', s); 
       
        image2plot=abs(CombWater_rr);
       imagesc(abs(image2plot)),
       title('Water Data Amp.'); axis('off');colormap default ;      
        print(figs, '-append','-bestfit', '-dpsc2', s); 

close all;


   figs=figure('visible', 'off'); 
   
        image2plot=(angle(CombData_rr));
       imagesc(angle(image2plot));
       title('MRSI Data Phase'); axis('off');colormap default ;
        print(figs, '-append','-bestfit', '-dpsc2', s); 
        
        image2plot=(angle(CombWater_rr));
       imagesc(angle(image2plot)),
       title('Water Data Phase'); axis('off');colormap default ;
        
        print(figs, '-append','-bestfit', '-dpsc2', s); 

close all;


end

