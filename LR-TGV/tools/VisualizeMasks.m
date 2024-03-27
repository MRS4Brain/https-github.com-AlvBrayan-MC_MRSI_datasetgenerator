function VisualizeMasks( mrsiData ,HMask, BMask,SkMask,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Visualization_Masks.ps',name_fig);
if exist(s);delete(s);end
Size_data=size(mrsiData);

MRSIData_ctrr = ifft(ifft(mrsiData,[],3),[],4);

WaterAmp_crr=squeeze(mean(abs(fft(MRSIData_ctrr,[],2)),2));
%WaterPh_crr=squeeze(angle(USData_crrc(:,:,:,1)));
%WaterPh_crr=squeeze(angle(sum(MRSIData_ctrr(:,1:5,:,:),2)));
if size(MRSIData_ctrr,1)==1
  WaterAmp_crr= reshape(WaterAmp_crr,[1 size(WaterAmp_crr)]) ;
end
Data_rr=squeeze(sum(WaterAmp_crr.^2,1));

if(sum(size(Data_rr)==size(HMask))<3)
	Data_rr=imresize(Data_rr,size(HMask));
end

figs=figure('visible', 'off');

imagesc(abs(Data_rr)),colorbar,daspect([ 1 1 1] );
axis('off');
colormap default ;title(' Data')
print(figs, '-append', '-dpsc2', s);

imagesc(abs(Data_rr).*HMask),colorbar,daspect([ 1 1 1] );
axis('off');
colormap default; title('Data Head Masked')
print(figs, '-append', '-dpsc2', s);

imagesc(abs(Data_rr).*BMask),colorbar,daspect([ 1 1 1] );
axis('off');
colormap default; title('Data Brain Masked')
print(figs, '-append', '-dpsc2', s);

imagesc(abs(Data_rr).*SkMask),colorbar,daspect([ 1 1 1] );
axis('off');
colormap default; title('Data Skull Masked')
print(figs, '-append', '-dpsc2', s);

close all;




end

