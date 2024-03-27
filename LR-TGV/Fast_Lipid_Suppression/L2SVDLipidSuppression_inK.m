function [mrsiDataLR_tkk,LipidFree_img,Lipid_Vol,VariableBeta]  = L2SVDLipidSuppression_inK( mrsiData_tkk, mrsiReconParams, Nbasis, Nfit, Beta )
 % mrsiReconParams.mrsiData dims: time-k-k
 Data_kkt=permute(mrsiData_tkk,[2,3,1]);
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
b0map2D=imresize(mean(mrsiReconParams.b0map,3),[N(1),N(2)]);
 
%{
for l = 1 : size(Data_rrt,3)
  % Data_rrt(:,:,l)=exp(-2*pi*1i*(l-1) * b0map2D / mrsiReconParams.mrProt.samplerate).*squeeze(Data_rrt(:,:,l));
    Data_rrt(:,:,l)=exp(-2*pi*1i*(l-1) * mrsiReconParams.WaterFreqMap  / mrsiReconParams.mrProt.samplerate).*squeeze(Data_rrt(:,:,l));
end
%}

Data_rrf = fft(Data_rrt,[],3);%rrf
Data_kkf = fft(Data_kkt,[],3);%rrf

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
  
%meta_mask=imresize(mean(mrsiReconParams.BrainMask,3),[N(1),N(2)]);
%meta_mask=round(meta_mask/max(meta_mask(:)));
meta_mask=mrsiReconParams.BrainMask2D;
 %% apply L2 lipid-basis recon to dual-density image

% Beta=100;

%Nbasis=30;

%img_combined_L2 = combined_2avg.*repmat(meta_mask,[1 1 N(end)]);

%Lipid_Vol=squeeze(sum(abs(img_combined_L2(:,:,low_bnd_L:high_bnd_L)),3));
Lipid_Vol=squeeze(sum(abs(Data_kkf(:,:,low_bnd_L:high_bnd_L)).^2,3));
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
FatThres =[1 1E-1 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 1E-8 1E-9];
%FatThres =[1 3E-1 1E-1 3E-2 1E-2 3E-3 1E-3 3E-4 1E-4 3E-5];
MinLip=min(Lipid_Vol(Lipid_Vol>0));
%Beta=0.1/FatThres(sum(MinLip<FatThres));
Beta=Beta/FatThres(sum(MinLip<FatThres));

Lipid = make_QR_LipidBasis_inK(Data_rrf, lipid_mask,Nbasis); 

LipidFree_img =zeros(N);
VariableBeta =zeros(N(1),N(2));

 %data_dir=fullfile(pwd, ['reconResults_',mrsiReconParams.NameData]);
 %s=[data_dir, '/Lipid_Contamination_Masks.ps'];      
 s=[mrsiReconParams.NameData, '_Lipid_Contamination_Masks.ps'];      
   delete(s);  
   figs=figure('visible', 'off'); 
for Fatind= 2:numel(FatThres)
    
    ContamMask= (Lipid_Vol>FatThres(Fatind)) & (Lipid_Vol<=FatThres(Fatind-1));
    %img_LipThres_L2=img_combined_L2.*repmat(ContamMask,[1 1 N(end)]);
    img_LipThres_L2=Data_kkf.*repmat(ContamMask,[1 1 N(end)]);
    
    imagesc(ContamMask);%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
    
    AdaptBeta=Beta*FatThres(Fatind-1);
    VariableBeta = VariableBeta + AdaptBeta*ContamMask;
    Lipid_inv = inv( eye(N(3)) + AdaptBeta * (Lipid * Lipid') );
    
    [Uorig,Sorig,Vorig] = svd(reshape(img_LipThres_L2,[],N(end)),0);
   %{ 
    for a = 1:N(1)
        for b = 1:N(2)
            if(ContamMask(a,b))
             LipidFree_img(a,b,:)=squeeze(Lipid_inv)*squeeze(img_LipThres_L2(a,b,:));

            % plot(1:N(3),squeeze(real(LipidFree_img(a,b,:))),1:N(3),squeeze(real(img_LipThres_L2(a,b,:))));
           
            end
        end
    end
%}
    
for comp = 1:Nfit
     %SVD
      V(:,comp) = (squeeze(Lipid_inv) * Vorig(:,comp));
  end
    
   Uorig_rsh=reshape(Uorig,N(1),N(2),[]);   
   SVD_img=ReformSVDProduct(Uorig_rsh(:,:,1:Nfit),Sorig(1:Nfit,1:Nfit), Vorig(:,1:Nfit),numel(N)-1);
   EpsSVD=img_LipThres_L2-SVD_img;
   LipidFree_img = LipidFree_img+ repmat(ContamMask,[1 1 N(end)]).*(EpsSVD + ReformSVDProduct(Uorig_rsh(:,:,1:Nfit),Sorig(1:Nfit,1:Nfit), (V(:,1:Nfit)) ,numel(N)-1));
   
    %EpsSVD=img_combined_L2-SVD_img;

end


%Data_rrt=ifft(LipidFree_img,[],3);
Data_kkt=ifft(LipidFree_img,[],3);
%Data_rrt=LipidFree_img; % if suppresion is done in time

%{
for l = 1 : size(Data_rrt,3)
   %Data_rrt(:,:,l)=exp(-2*pi*1i*(l-1) * -b0map2D / mrsiReconParams.mrProt.samplerate).*squeeze(Data_rrt(:,:,l));
    Data_rrt(:,:,l)=exp(+2*pi*1i*(l-1) * mrsiReconParams.WaterFreqMap  / mrsiReconParams.mrProt.samplerate).*squeeze(Data_rrt(:,:,l));
end
%}
%LipidFree_img=fft(Data_rrt,[],3);

%mrsiDataLR_tkk =fft(fft(permute(Data_rrt,[3,1,2]),[],2),[],3);
mrsiDataLR_tkk = permute(Data_kkt,[3,1,2]);

end

