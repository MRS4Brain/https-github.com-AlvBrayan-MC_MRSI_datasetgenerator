function [mrsiDataLR_tkk,LipidFree_frr,Lipid_Vol]  = ProjSVDLipidSuppression_fullKSpace( mrsiData_tkk, mrsiReconParams, Nbasis, Nfit ,NameData)
 % mrsiReconParams.mrsiData dims: time-k-k
%{ 
Data_kkt=permute(mrsiData_tkk,[2,3,1]);
% *************************************************************************
% LOW RANK TGV recon
%*************************************************************************
[mrsiReconParams.US, mrsiReconParams.V]=LowRankTGV(Data_tkk,mrsiReconParams);
TempRes_rrt=formTensorProduct(mrsiReconParams.US,mrsiReconParams.V,nDims-1);
%mrsiReconParams.FilteredData_trr=permute( TempRes_rrt,[ndims( TempRes_rrt),1:(ndims( TempRes_rrt)-1)]);
Data_rrt=permute( TempRes_rrt,[ndims( TempRes_rrt),1:(ndims( TempRes_rrt)-1)]);
%}
Data_rrt=ifft(ifft(permute(mrsiData_tkk,[2,3,1]),[],1),[],2);

N = size(Data_rrt);
HzpP=mrsiReconParams.mrProt.samplerate/N(3);
%b0map2D=imresize(mean(mrsiReconParams.b0map,3),[N(1),N(2)]);
 
%Phase correction ... does not improve results
for l = 1 : size(Data_rrt,3)
 %   Data_rrt(:,:,l)=exp(-2*pi*1i*(l-1) * mrsiReconParams.WaterFreqMap / mrsiReconParams.mrProt.samplerate).*squeeze(Data_rrt(:,:,l));
end

Data_rrf=fft(Data_rrt,[],3);

%combined_2avg = Data_rrt;%fft(Data_rrt,[],3);%rrf% if suppresion is done in time
 
%lipid_mask=imresize(sum(mrsiReconParams.SkMask,3),[N(1),N(2)]);
%lipid_mask=round(lipid_mask/max(lipid_mask(:)));
lipid_mask=mrsiReconParams.SkMask2D;

low_bnd_L=round(200/HzpP);
high_bnd_L=round(500/HzpP);


Lipid_Vol=squeeze(sum(abs(Data_rrf(:,:,low_bnd_L:high_bnd_L)),3));
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
LThr=0;%quantile(Lipid_Vol(:),0.5);
lipid_mask=((Lipid_Vol.*lipid_mask)>LThr);


meta_mask=mrsiReconParams.BrainMask2D;

Brain_rrf = Data_rrf.*repmat(meta_mask,[1 1 N(end)]);
OutsideBrain_rrf=Data_rrf.*repmat(~meta_mask,[1 1 N(end)]);

Lipid = make_SVD_LipidBasis(Data_rrf, lipid_mask,Nbasis); 
   
    
[Uorig,Sorig,Vorig] = svd(reshape(Brain_rrf,[],N(end)),0);
LipidProj=Lipid *Lipid';
 V(:,1:Nfit) =  LipidProj * Vorig(:,1:Nfit);
 
 Uorig_rsh=reshape(Uorig,N(1),N(2),[]);   
  
 Lipid_rrf= ReformSVDProduct(Uorig_rsh(:,:,1:Nfit),Sorig(1:Nfit,1:Nfit), (V(:,1:Nfit)) ,numel(N)-1) + OutsideBrain_rrf ;
   %EpsSVD=Brain_rrf-SVD_img;
   %LipidFree_rrf = (EpsSVD + ReformSVDProduct(Uorig_rsh(:,:,1:Nfit),Sorig(1:Nfit,1:Nfit), (Vorig(:,1:Nfit)-V(:,1:Nfit)) ,numel(N)-1));
   
    %EpsSVD=img_combined_L2-SVD_img;

Lipid_rrt=ifft(Lipid_rrf,[],3);


%Phase correction ... does not improve results
for l = 1 : size(Lipid_rrt,3)
   %Lipid_rrt(:,:,l)=exp(-2*pi*1i*(l-1) * -mrsiReconParams.WaterFreqMap / mrsiReconParams.mrProt.samplerate).*squeeze(Lipid_rrt(:,:,l));
end

Data_rrt=ifft(ifft(permute(mrsiData_tkk,[2,3,1]),[],1),[],2);
LipidFree_frr=permute(fft(Data_rrt-Lipid_rrt,[],3),[3,1,2]);

%Lipid_kkt=fft(fft(ifft(Lipid_rrf,[],3),[],1),[],2);
%mrsiDataLR_tkk =mrsiData_tkk-permute(Lipid_kkt,[3,1,2]);
%LipidFree_frr = fft(ifft(ifft(mrsiDataLR_tkk,[],2),[],3),[],1);
mrsiDataLR_tkk =ifft(fft(fft(LipidFree_frr,[],2),[],3),[],1);

 s=['./',mrsiReconParams.Log_Dir,'/',NameData, '_Lipid_Image.ps'];      
   delete(s);  
   figs=figure('visible', 'off'); 
   
  imagesc( sum(abs(Lipid_rrf),3));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
    
   s=['./',mrsiReconParams.Log_Dir,'/',NameData, '_Lipid_Free_Image.ps'];      
   delete(s);  
   figs=figure('visible', 'off'); 
   
  imagesc( squeeze(sum(abs( LipidFree_frr),1)));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
    
    s=['./',mrsiReconParams.Log_Dir,'/',NameData, '_Lipid_Mask.ps'];      
   delete(s);  
   figs=figure('visible', 'off'); 
   
  imagesc( lipid_mask);%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
end

