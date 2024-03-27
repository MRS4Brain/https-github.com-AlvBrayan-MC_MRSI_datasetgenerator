function [mrsiDataLR_ctkk, NBasis,LipidProj]  = ProjSVDLipidSuppressionWithRecon( mrsiData_ctkk, mrsiReconParams ,NameData)
 % mrsiReconParams.mrsiData dims: time-k-k


OrderLip=128;
N = size(mrsiData_ctkk);
IOP=diag(ones(N(2),1));

Brainmask_1rr1=reshape(mrsiReconParams.BrainMask2D,[1 N(3) N(4) 1]);

%Determine approximatif Data_rrf
SENSE_c1rr=reshape(mrsiReconParams.SENSE,[N(1) 1 N(3) N(4)]);
Data_rrf=squeeze(sum(conj(SENSE_c1rr).*ifft(ifft(mrsiData_ctkk,[],3),[],4),1)); % t-r-r
Data_rrf=permute(fft(Data_rrf,[],1),[2,3,1]);
%Determine Lipids_rrf
[Uorig,Sorig,Vorig] = svd(reshape(permute(mrsiData_ctkk,[1 3 4 2]),[],N(2)),0);
V_tc=Vorig(:,1:OrderLip);
S=Sorig(1:OrderLip,1:OrderLip);
U_ckkc=reshape(Uorig(:,1:OrderLip), N(1),N(3),N(4),[]);

kmask=mrsiReconParams.kmask;
%alpha = mean(kmask(:))*mrsiReconParams.mu_tv/reduction; %  % correted for undersampling
%reduction = 2^(-8);     % usually there is no need to change this
alpha = 1E-8;%/reduction; %  % correted for undersampling
maxit = 1000 ;          % use 1000 Iterations for optimal image quality        % usually there is no need to change this
Threshold=1E-6;
U_rrc=zeros(N(3),N(4),OrderLip);
fprintf('Start Lipid Reconstruction...\n');
parfor k=1:OrderLip
    % fprintf([ 'Processing Component ', num2str(k), ' ...\n']);
    [U_rrc(:,:,k) e] = tgv2_l2_2D_multiCoil(U_ckkc(:,:,:,k),mrsiReconParams.SENSE, kmask, 2*alpha, alpha, maxit,Threshold);
end

Lipids_rrt=formTensorProduct(U_rrc, V_tc*S,2);

N = size(Lipids_rrt);
HzpP=mrsiReconParams.mrProt.samplerate/N(end);
 
Lipids_rrf=fft(Lipids_rrt,[],3);

lipid_mask=mrsiReconParams.SkMask2D;

[~,high_bnd_L]=min(abs(mrsiReconParams.LipidMaxPPM - mrsiReconParams.ppm));
[~,low_bnd_L]=min(abs(mrsiReconParams.LipidMinPPM  - mrsiReconParams.ppm));

Lipid_Vol=squeeze(sum(abs(Lipids_rrt(:,:,low_bnd_L:high_bnd_L)),3));
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
SortedLip=sort(Lipid_Vol(:),'descend');
LThr=0;%SVD & Thrs=0 seems to be the best option for 2D-MRSI

%Noise=Lipid_Vol(:);%=Lipid_Vol(mrsiReconParams.SkMask2D==1);
%LThr=mean(Noise)+0.5*std(Noise);%+%1.0*std(Noise); %quantile(Lipid_Vol(:),0.5);%0; 

[~,NumVox]=min(abs(LThr-SortedLip(:)));
%fprintf(["NumVox="+num2str(NumVox)]);
lipid_mask=((Lipid_Vol.*lipid_mask)>LThr);

meta_mask=mrsiReconParams.BrainMask2D;

[LipidProj, Slipid, NBasis,Lipid ]  = make_SVD_LipidBasis(Lipids_rrf,Data_rrf,lipid_mask,mrsiReconParams.L2SVDparams.PercentThres,mrsiReconParams); 
%[LipidProj,NBasis] = make_L2_LipidBasis( LipData_rrf, lipid_mask ,mrsiReconParams.L2SVDparams.PercentThres,mrsiReconParams);
%[LipidProj,NBasis] = make_L2MetabCorr_LipidBasis( LipData_rrf, lipid_mask ,mrsiReconParams.L2SVDparams.PercentThres,mrsiReconParams);

mrsiDataLR_ctkk=0*mrsiData_ctkk;
for c=1:size(mrsiData_ctkk,1)
	OrigData_rrf=squeeze(permute(fft(mrsiData_ctkk(c,:,:,:),[],2),[3,4,2,1]));
	OrigData_rrf=reshape(OrigData_rrf,[],N(end)) * (IOP-LipidProj);
	mrsiDataLR_ctkk(c,:,:,:)=permute(ifft(reshape(OrigData_rrf,N),[],3),[3,1,2]);
end
clear OrigData_rrf
 
Lipid_rrf= reshape(Data_rrf,[],N(end)) * LipidProj;
Lipid_rrf= reshape( Lipid_rrf,[N(1),N(2),N(3)]);

LipidFree_rrf= reshape(Data_rrf,[],N(end)) * (IOP-LipidProj);
LipidFree_rrf= reshape( LipidFree_rrf,[N(1),N(2),N(3)]);


if  ~isempty(NameData)
    s=['./',mrsiReconParams.Log_Dir,'/',NameData, '_Lipid'];
    VisualizeSpectral( ifft(Lipid(low_bnd_L:high_bnd_L,:),[],1),Slipid ,s)
    
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
   
    imagesc( ~lipid_mask.*squeeze(sum(abs( LipidFree_rrf),3)));%,[ 0, 10*mean(image2plot(:))] )
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
   
    imagesc( squeeze(sum(abs( LipidFree_rrf),3)));%,[ 0, 10*mean(image2plot(:))] )
    title('Lipid-Free Data in Image');
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
    
    imagesc( lipid_mask);%,[ 0, 10*mean(image2plot(:))] )
    title('Lipid masks');
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
end
end

