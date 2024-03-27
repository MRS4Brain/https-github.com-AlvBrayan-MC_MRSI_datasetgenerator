function [LipidProj, Slipid, Nbasis,Lipids ] = make_WdSVD_LipidBasis( lipid_rrf,data_rrf, lipid_mask, SVratio ,ReconParams)
%GET_LIPIDBASIS Summary of this function goes here
%   Detailed explanation goes here


count = 0;
N = size(data_rrf);
IOP=diag(ones(N(end),1));
Lipid_Basis = zeros(N(3), sum(lipid_mask(:)));

%HzpP=ReconParams.mrProt.samplerate/N(end);
%low_bnd_L=round(100/HzpP);
%high_bnd_L=round(600/HzpP);

[~,high_bnd_L]=min(abs(ReconParams.LipidMaxPPM - ReconParams.ppm));
[~,low_bnd_L]=min(abs(ReconParams.LipidMinPPM  - ReconParams.ppm));

[~,high_LipId]=min(abs(-1.9 - ReconParams.ppm));
[~,low_LipId]=min(abs(0  - ReconParams.ppm));
%high_LipId=1;
low_LipId=1056;

% SVD basis of the skull lipids 
%InitEnergyMap=sum(abs(data_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);

LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
[~,Sorig,Vorig] = svd(reshape(lipid_rrf(LipidMask_rrf>0),[],N(end)),'econ');%0);
lipid_stack_rf=reshape(lipid_rrf(LipidMask_rrf>0),[], size(lipid_rrf,3));

lipid_Wdw_rrf=lipid_rrf(:,:,high_LipId:low_LipId);
LipidMask_rrf=repmat(lipid_mask,[1 1 size(lipid_Wdw_rrf,3)]);
lipid_Wdw_stack_rf=reshape(lipid_Wdw_rrf(LipidMask_rrf>0),[], size(lipid_Wdw_rrf,3));
[~,SWdw,VWdw] = svd(reshape(lipid_Wdw_rrf(LipidMask_rrf>0),[], size(lipid_Wdw_rrf,3)),'econ');%0);

Window=zeros(N(end),1);
Window(high_LipId:low_LipId)=1;
WindowOPSqr=diag(Window)*N(end)/(low_LipId-high_LipId);
WindowOP=IOP(1:N(end),high_LipId:low_LipId);


N1=4;N2=64;
LipidProj=Vorig(high_LipId:low_LipId,1:N1)*Vorig(:,1:N1)';
WdwLipidProj=VWdw(:,1:N2)*VWdw(:,1:N2)';
%LipFree_rrf=reshape(reshape(data_rrf,[],N(3))*(IOP-WindowOP*WdwLipidProj*LipidProj),[N(1) N(2) N(3)]);
Nbasis=N1;
LipidProj=WindowOP*WdwLipidProj*LipidProj;
Lipids=Vorig(:,1:Nbasis);
Slipid=Sorig(1:Nbasis,1:Nbasis);
%LipFree_rrf = data_rrf-Lip_rrf;
%{
LipFree1_rrf=LipFree_rrf;
lipid_rrf=LipFree_rrf;
LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
[~,Sorig,Vorig] = svd(reshape(lipid_rrf(LipidMask_rrf>0),[],N(end)),'econ');%0);

Nbasis=64;%32
Nstep=32;%16
while Nstep>=1
    LipidProj=Vorig(:,1:Nbasis)*Vorig(:,1:Nbasis)';
    LipFree_rrf=reshape(reshape(LipFree1_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3)]);
% WdwLipidProj=VWdw(:,1:N2)*VWdw(:,1:N2)';
   % LipidProj=Vorig(high_LipId:low_LipId,1:Nbasis)*Vorig(high_LipId:low_LipId,1:Nbasis)';
   %  LipWdw_rrf=reshape(reshape(data_rrf(:,:,high_LipId:low_LipId),[],low_LipId-high_LipId+1)*(WdwLipidProj),[N(1) N(2) low_LipId-high_LipId+1]);
% LipFree_rrf=reshape(reshape(data_rrf(:,:,:),[],N(3))*(IOP-WindowOP*LipidProj),[N(1) N(2) N(3)]);
% LipFree_rrf=reshape(reshape(data_rrf(:,:,:),[],N(3))*(IOP-LipidProj),[N(1) N(2) N(3)]);
 %Lip_rrf=reshape(reshape(data_rrf(:,:,high_LipId:low_LipId),[],low_LipId-high_LipId+1)*(WdwLipidProj*LipidProj),[N(1) N(2) N(3)]);
   %  LipFree_rrf = data_rrf-Lip_rrf;
     %LipFree_rrf = data_rrf-Lip_rrf;
  %plot(1:1056,real(squeeze(data_rrf(a,b,:))),high_LipId:low_LipId,real(squeeze(LipWdw_rrf(a,b,:))),1:1056,real(squeeze(Lip_rrf(a,b,:))),1:1056,real(squeeze(LipFree_rrf(a,b,:))))
   EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
    
    SkullEnergy=mean(EnergyMap(lipid_mask>0))
    BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0))
    %BrainEnergy/SkullEnergy;
    if SkullEnergy>(BrainEnergy/SVratio);
        Nbasis=Nbasis+Nstep
    else
        Nbasis=Nbasis-Nstep
    end
    Nstep=Nstep/2 ;
end
%[~,Sorig,Vorig] = svds(reshape(Lipid_img,[],N(end)),Nbasis);%0);%crashing some times
%SingVal=diag(Sorig)./cumsum(diag(Sorig));%(abs(Sorig(1,1));
%Nbasis=min(find(abs(diff(SingVal))<SVratio)); % too much variation with diff

%SingVal=log(diag(Sorig))./cumsum(log(diag(Sorig)));
%Nbasis=min(find(abs(SingVal)<SVratio)); 
%figure;
%imagesc(sum(abs(reshape(LipFree_rf,[N(1) N(2) N(3)])),3));

if Nbasis>ReconParams.L2SVDparams.NBasisMax
    Nbasis=ReconParams.L2SVDparams.NBasisMax;
end
if Nbasis<ReconParams.L2SVDparams.NBasisMin
    Nbasis=ReconParams.L2SVDparams.NBasisMin;
end

Lipids=Vorig(:,1:Nbasis);
Slipid=Sorig(1:Nbasis,1:Nbasis);
%{ 
%Do some HSVD filtering on the basis (Didn't work well)
if ReconParams.L2SVDparams.LipidT2star>0
    for comp=1:Nbasis
        [freqs,dampings,basis,amps] = HSVD(fft(Lipids(:,comp)), ReconParams.mrProt.samplerate, 128);
        indx                 = find(-dampings >= 1/(1E-3*ReconParams.L2SVDparams.LipidT2star));
        Lipids(:,comp) = ifft(sum(basis(:, indx) * diag(amps(indx), 0), 2));
    end
end
%}

%Lipid_img=reshape(data_rrf.*repmat(lipid_mask,[1 1 N(end)]),[],N(end));
%LipidVol=sum(abs(Lipid_img).^2,2);
%[~,Ind]=sort(LipidVol(:));
%Lipid_img_sorted=Lipid_img(flip(Ind),:);
%[Q,R]=qr(Lipid_img_sorted(1:Nbasis,:));


%[Q,R]=qr(Lipids'); %for SVD basis of the skull lipids 
%useless already orthonormal
%Lipids=R';
%Lipids=Lipids*sqrt(inv(diag(diag(Lipids'*Lipids))));

LipidProj=Lipids*Lipids';
%}
end

