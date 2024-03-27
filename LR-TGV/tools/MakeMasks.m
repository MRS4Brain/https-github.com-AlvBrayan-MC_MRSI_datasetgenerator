function [ ImMask, BrainMask,SkMask ] = MakeMasks(mrsiReconParams )

FiltParam.FilterFreq=0;%5;%hz
FiltParam.Water_minFreq=mrsiReconParams.FiltParam.Water_minFreq;%-60;%hz
FiltParam.Water_maxFreq=mrsiReconParams.FiltParam.Water_maxFreq;%60hz
FiltParam.Comp=16;
mrProt=mrsiReconParams.Water_mrProt;
kmask=mrsiReconParams.Water_kmask;

%SENSE_ctrr=permute(repmat(mrsiReconParams.RoemersCoefs,[1 1 1 size(mrsiReconParams.Water_ctkk,2)]),[1,4,2,3]);
%raw_tkk = fft(fft(conj(SENSE_ctrr).*ifft(ifft(mrsiReconParams.Water_ctkk,[],3),[],4),[],3),[],4);
%SENSE_ctrr=permute(repmat(mrsiReconParams.RoemersCoefs,[1 1 1 size(mrsiReconParams.mrsiData,2)]),[1,4,2,3]);
%raw_tkk = fft(fft(conj(SENSE_ctrr).*ifft(ifft(mrsiReconParams.mrsiData,[],3),[],4),[],3),[],4);
%raw_tkk = squeeze(sum(raw_tkk,1));

raw_tkk=squeeze(sum(mrsiReconParams.mrsiData,1));


clear SENSE_ctrr


if mrsiReconParams.L2SVDparams.PercentThres==0
Energy_rr=squeeze(sum(abs(ifft(ifft(raw_tkk(1:10,:,:),[],2),[],3)).^2,1));

ImMask_temp=Energy_rr>5*quantile(Energy_rr(:),0.05);;
SkMask_temp=0*ImMask_temp;

else 

Lipid_tkk=zeros(size(raw_tkk));
Water_tkk=zeros(size(raw_tkk));


[Lipid_tkk,Water_tkk, ~, ~] = FilterMRSIData(raw_tkk,mrProt,FiltParam,kmask);

Lipid_frr=fft(ifft(ifft(Lipid_tkk,[],3),[],2),[],1);
Water_frr=fft(ifft(ifft(Water_tkk,[],3),[],2),[],1);



[~,high_bnd_L]=min(abs(mrsiReconParams.LipidMaxPPM - mrsiReconParams.ppm));
[~,low_bnd_L]=min(abs(mrsiReconParams.LipidMinPPM  - mrsiReconParams.ppm));

size(Lipid_frr)
%HzpP=mrProt.samplerate/size(Lipid_frr,1);
%low_bnd_L=round(600/HzpP);%round(300/HzpP);
%high_bnd_L=round(2000/HzpP);%round(600/HzpP);

Lipid_Vol=squeeze(sum(abs(Lipid_frr(low_bnd_L:high_bnd_L,:,:)),1));
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
LThr=mrsiReconParams.Threshold_LipMask;%0;
SkMask_temp=(Lipid_Vol>LThr);

Water_Vol=squeeze(sum(abs(Water_frr),1));

Energy_rr=Lipid_Vol+Water_Vol;

% LThr=mrsiReconParams.Threshold_LipMask;
% Total_Vol=(Lipid_Vol+Water_Vol);
% Total_Vol=Total_Vol/max(Total_Vol(:));
% ImMask_temp=(Total_Vol>LThr);
ImMask_temp = imfill(SkMask_temp,'holes');
% Dilate Masks
%{
ImMask=zeros(size(ImMask_temp));
SkMask=zeros(size(ImMask_temp));
NB_Vox=zeros(3,3);
for a=2:(size(ImMask_temp,1)-1)
    for b=2:(size(ImMask_temp,2)-1)
        NB_Vox=ImMask_temp((a-1):(a+1),(b-1):(b+1));
        ImMask(a,b)=(sum(NB_Vox(:))>0);
       %NB_Vox=SkMask_temp((a-1):(a+1),(b-1):(b+1));
        %SkMask(a,b)=(sum(NB_Vox(:))>0);
    end
end
SkMask=SkMask_temp;
%}
end
SkMask=SkMask_temp;
ImMask=ImMask_temp;
BrainMask=ImMask-SkMask;
BrainMask(BrainMask<0)=0;

s=['./',mrsiReconParams.Log_Dir,'/','AutoMasks_Lipid_Head',mrsiReconParams.NameData,'.ps'];
if exist(s);delete(s);end

figs=figure('visible','off');
imagesc(Energy_rr);colorbar;%,[ 0, 10*mean(image2plot(:))] );
title('Original Data');
print(figs, '-append', '-dpsc2', s);
imagesc(SkMask);colorbar;
title('Lipid Mask');
print(figs, '-append', '-dpsc2', s); 
imagesc(ImMask);colorbar;
title('Head Mask');
print(figs, '-append', '-dpsc2', s); 
imagesc(BrainMask);colorbar;
title('BrainMask');
print(figs, '-append', '-dpsc2', s); 
end

