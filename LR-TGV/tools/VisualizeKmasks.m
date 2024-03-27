function VisualizeKmasks( mrsiData_ctkk,kmask_data, HeadWater_ctkk ,kmask_Head,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

s=sprintf('%s_Visualization_Kmasks.ps',name_fig);
if exist(s);delete(s);end
Size_data=size(mrsiData_ctkk);

mrsiData_kk = squeeze(sum(sum(abs(fft(mrsiData_ctkk,[],2)),1),2));
HWater_kk = squeeze(sum(sum(abs(fft(HeadWater_ctkk,[],2)),1),2));


figs=figure('visible', 'off');
subplot(2,2,1),imagesc(mrsiData_kk,[0 0.5*max(mrsiData_kk(:))]),colorbar;
axis('off');
colormap default ;title('MRSI Data k-space Combined Amplitude')
subplot(2,2,3),imagesc(kmask_data),colorbar;
axis('off');
colormap default; title('MRSI Data k-mask')
subplot(2,2,2),imagesc(HWater_kk,[0 0.5*max(HWater_kk(:))]),colorbar;
axis('off');
colormap default; title('Head Water k-space Combined Amplitude')
subplot(2,2,4),imagesc(kmask_Head),colorbar;
axis('off');
colormap default ; title('Head Water k-mask')

print(figs, '-append', '-dpsc2', s);

close all;




end

