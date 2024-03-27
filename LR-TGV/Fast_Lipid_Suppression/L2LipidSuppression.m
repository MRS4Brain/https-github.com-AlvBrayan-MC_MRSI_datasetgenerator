function [mrsiDataLR_tkk,mrsiData_LipidSup_rrf, mrsiData_LipidPenality_rrf]  = L2LipidSuppression( mrsiData_tkk, mrsiReconParams, Beta ,NameData,meta_mask)
 % mrsiReconParams.mrsiData dims: time-k-k
Data_rrt=ifft(ifft(permute(mrsiData_tkk,[2,3,1]),[],1),[],2);
Data_kkf=fft(permute(mrsiData_tkk,[2,3,1]),[],3);
N = size(Data_rrt);
HzpP=mrsiReconParams.mrProt.samplerate/N(3);
% b0map2D=imresize(mean(mrsiReconParams.b0map,3),[N(1),N(2)]);


combined_2avg = fft(Data_rrt,[],3);%rrf
 
%lipid_mask=imresize(sum(mrsiReconParams.SkMask,3),[N(1),N(2)]);
%lipid_mask=round(lipid_mask/max(lipid_mask(:)));
lipid_mask=mrsiReconParams.SkMask2D;

%Skull_leaking_approx=Conv_PSF_2D(mrsiReconParams.SkMask,mrsiReconParams.kmask);

low_bnd_L=round(200/HzpP);
high_bnd_L=round(500/HzpP);

Lipid_Vol=squeeze(sum(abs(combined_2avg(:,:,low_bnd_L:high_bnd_L)),3));
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
%LThr=quantile(Lipid_Vol(:),0.5);
LThr=0;
lipid_mask=((Lipid_Vol.*lipid_mask)>LThr);
  

 %% apply L2 lipid-basis recon to dual-density image
mrsiData_LipidSup_kkf = Data_kkf;%.*repmat(meta_mask,[1 1 N(end)]);


betamin = 1;%100;%6.5e-1;      % beta chosen to match the data consistency of iterative L1-lipid-basis recon
betamax = 1E7;
%betamax = 0;
%betamin = 1000;
Lipid = make_LipidBasis(combined_2avg, lipid_mask); 

Lipid_Vol=squeeze(sum(abs(combined_2avg(:,:,low_bnd_L:high_bnd_L)).^2,3));
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));

betaMap=Lipid_Vol*betamax+betamin;
%betaMap=Skull_leaking_approx.^2*betamax+betamin;
%BetaVals=(2.^linspace(0,log2(max(betaMap(:))/betamin),30) )*betamin;
%BetaVals=(2.^linspace(0,log2(max(betaMap(:))/betamin),10) )*betamin;
BetaVals=Beta;

%for val=1:numel(BetaVals)
 %  Lipid_inv(val,:,:) = inv( eye(N(3)) + BetaVals(val) * (Lipid * Lipid') );
%end

 Lipid_inv = inv( eye(N(3)) + BetaVals * (Lipid * Lipid') );

mrsiData_LipidPenality_kkf=zeros(size(mrsiData_LipidSup_kkf));
ActualbetaMap=0*betaMap;

if ~exist('meta_mask','var') || isempty(meta_mask)
    meta_mask=ones(size(Data_kkf,1),size(Data_kkf,2));
end

for a = 1:size(Data_kkf,1)
    for b= 1:size(Data_kkf,2)
        %if meta_mask(a,b)  
           % [~, index_Inv] = min(abs(betaMap(a,b) -BetaVals) );
           index_Inv=1;
            %ActualbetaMap(a,b)=BetaVals(index_Inv);
           % xi = Data_kkf(a,b,:);
           %mrsiData_LipidSup_kkf(a,b,:) = squeeze(Lipid_inv(index_Inv,:,:)) * xi(:);
           mrsiData_LipidSup_kkf(a,b,:) = Lipid_inv * squeeze(Data_kkf(a,b,:));
           mrsiData_LipidPenality_kkf(a,b,:)=(Lipid * Lipid')*squeeze(Data_kkf(a,b,:));
            
       % end
        
    end
end

mrsiData_LipidSup_rrf=ifft(ifft(mrsiData_LipidSup_kkf,[],2),[],1);
mrsiData_LipidPenality_rrf=ifft(ifft(mrsiData_LipidPenality_kkf,[],2),[],1);

mrsiDataLR_tkk =ifft(permute(mrsiData_LipidSup_kkf,[3,1,2]),[],1);%tkk



if ~isempty(NameData)
 
    s=['./',mrsiReconParams.Log_Dir,'/',NameData, '_LipidL2Sup_Images.ps'];
    delete(s);
    figs=figure('visible', 'off');
    imagesc( sum(abs(mrsiData_LipidPenality_rrf),3));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    title('Filtered out Lipids')
    print(figs, '-append', '-dpsc2', s);
   
    imagesc( sum(abs(combined_2avg),3));%,[ 0, 10*mean(image2plot(:))] )
    colormap default;colorbar;
    title('Original Data')
    print(figs, '-append', '-dpsc2', s);   
   
    imagesc( squeeze(sum(abs( mrsiData_LipidSup_rrf),3)));%,[ 0, 10*mean(image2plot(:))] )
    title('Lipid-Free Data');
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
    
    imagesc( lipid_mask);%,[ 0, 10*mean(image2plot(:))] )
    title('Lipid masks');
    colormap default;colorbar;
    print(figs, '-append', '-dpsc2', s);
end

end

