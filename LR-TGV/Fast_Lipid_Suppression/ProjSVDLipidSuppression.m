function [mrsiDataLR_tkk,LipidFree_frr,Lipid_rrf, NBasis,LipidProj]  = ProjSVDLipidSuppression( mrsiData_tkk,LipData_tkk, mrsiReconParams ,NameData,NBasis)
 % mrsiReconParams.mrsiData dims: time-k-k
Data_kkf=fft(permute(mrsiData_tkk,[2,3,1]),[],3);
Data_rrf=ifft(ifft(Data_kkf,[],1),[],2);


LipData_rrt=ifft(ifft(permute(LipData_tkk,[2,3,1]),[],1),[],2);

N = size(LipData_rrt);
HzpP=mrsiReconParams.mrProt.samplerate/N(3);
 
LipData_rrf=fft(LipData_rrt,[],3);

lipid_mask=mrsiReconParams.SkMask2D;


[~,high_bnd_L]=min(abs(mrsiReconParams.LipidMaxPPM - mrsiReconParams.ppm));
[~,low_bnd_L]=min(abs(mrsiReconParams.LipidMinPPM  - mrsiReconParams.ppm));


Lipid_Vol=squeeze(sum(abs(LipData_rrf(:,:,low_bnd_L:high_bnd_L)),3));
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
SortedLip=sort(Lipid_Vol(:),'descend');
LThr=0;%SVD & Thrs=0 seems to be the best option for 2D-MRSI

Noise=Lipid_Vol(:);%=Lipid_Vol(mrsiReconParams.SkMask2D==1);
if mrsiReconParams.Threshold_LipMask>=0
	LThr=mean(Noise)+mrsiReconParams.Threshold_LipMask*std(Noise);%+%1.0*std(Noise); %quantile(Lipid_Vol(:),0.5);%0; 
else
	LThr=0;
end

[~,NumVox]=min(abs(LThr-SortedLip(:)));
fprintf(["NumVox="+num2str(NumVox)]);
lipid_mask=((Lipid_Vol.*lipid_mask)>LThr);

meta_mask=mrsiReconParams.BrainMask2D;
Brainmask_rr1=reshape(meta_mask,[size(meta_mask,1) size(meta_mask,2) 1]);
%% apply L2 lipid-basis recon to dual-density image

if(nargin==4)% no Nbasis given
     [LipidProj, Slipid, NBasis,Lipid ]  = make_SVD_LipidBasis(LipData_rrf,Data_rrf, lipid_mask,mrsiReconParams.L2SVDparams.PercentThres,mrsiReconParams); 
elseif (nargin==5)
    [LipidProj, Slipid, NBasis,Lipid ]  = make_SVD_LipidBasis(LipData_rrf,Data_rrf, lipid_mask,mrsiReconParams.L2SVDparams.PercentThres,mrsiReconParams,NBasis);    
end
%[LipidProj,NBasis] = make_L2_LipidBasis( LipData_rrf,Data_rrf, lipid_mask ,mrsiReconParams.L2SVDparams.PercentThres,mrsiReconParams);

Lipid_kkf= reshape(Data_kkf,[],N(end)) * LipidProj;
Lipid_kkf= reshape( Lipid_kkf,[N(1),N(2),N(3)]);

Lipid_kkt=ifft(Lipid_kkf,[],3);
Lipid_rrf=ifft(ifft(Lipid_kkf,[],1),[],2);
mrsiDataLR_tkk =(mrsiData_tkk-permute(Lipid_kkt,[3,1,2]));
LipidFree_frr = fft(ifft(ifft(mrsiDataLR_tkk,[],2),[],3),[],1);

[~,high_bnd_L]=min(abs(0 - mrsiReconParams.ppm));
[~,low_bnd_L]=min(abs(-4.7  - mrsiReconParams.ppm));

if  ~isempty(NameData)
    if exist('Lipid')
	s=['./',mrsiReconParams.Log_Dir,'/',NameData, '_Lipid'];
	VisualizeSpectral( ifft(Lipid(low_bnd_L:high_bnd_L,:),[],1),Slipid ,s)
    end
    
    s=['./',mrsiReconParams.Log_Dir,'/',NameData, '_Lipid_Images.ps'];
    if exist(s);delete(s);end
    figs=figure('visible', 'off');
    imagesc( ~lipid_mask.*sum(abs(Lipid_rrf),3));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    title('Filtered out Lipids in Brain & Outside head')
    print(figs, '-append', '-dpsc2', s);
   
    imagesc( ~lipid_mask.*sum(abs(Data_rrf),3));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    title('Original Data in Brain & Outside head')
    print(figs, '-append', '-dpsc2', s);   
   
    imagesc( ~lipid_mask.*squeeze(sum(abs( LipidFree_frr),1)));%,[ 0, 10*mean(image2plot(:))] )
    title('Lipid-Free Data in Brain & Outside head');
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
    
    
    imagesc( sum(abs(Lipid_rrf),3));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    title('Filtered out Lipids in Image')
    print(figs, '-append', '-dpsc2', s);
   
    imagesc( sum(abs(Data_rrf),3));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    title('Original Data in Image')
    print(figs, '-append', '-dpsc2', s);   
   
    imagesc( squeeze(sum(abs( LipidFree_frr),1)));%,[ 0, 10*mean(image2plot(:))] )
    title('Lipid-Free Data in Image');
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
    
    imagesc( lipid_mask);%,[ 0, 10*mean(image2plot(:))] )
    title('Lipid masks');
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
end
end

