function [LipidProj,Lambda]  = make_L2MetabCorr_LipidBasis( LipData_rrf,data_rrf, lipid_mask ,SVratio,ReconParams)


count = 0;
N = size(data_rrf);
%Lipid_Basis = zeros(N(4), sum(lipid_mask(:)));

HzpP=ReconParams.mrProt.samplerate/N(end);
[~,high_bnd_L]=min(abs(ReconParams.LipidMaxPPM - ReconParams.ppm));
[~,low_bnd_L]=min(abs(ReconParams.LipidMinPPM  - ReconParams.ppm));
%SimData=h5read('Dataset_1E4Train_1E3Test_25PkW_5-20SNR_0NbBL_AcqDel_1056Npt_17Metab.h5','/train/spectra' );
%DataClean_rf=transpose(SimData.r+1j*SimData.i);
%[~,Sorig,VorigMetab] = svd(DataClean_rf,'econ');%0);
%MetabProj=VorigMetab(:,1:16)*VorigMetab(:,1:16)';
%load('DataCleanSpectra');

LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
Lipid_stack_rf=reshape(LipData_rrf(LipidMask_rrf>0),[],N(end));

LipidProj=zeros(N(end),N(end));
%[LipidProj, Slip, NBasis]  = make_SVD_LipidBasis(data_rrf, lipid_mask,SVratio,ReconParams); 
IOP=diag(ones(N(end),1));

Lambda=1;
Ratio=999;
step=0;
LoopRatio=SVratio;%*SVratio;%*2;
damp=1;

[X,Y] = ndgrid(1:N(1), 1:N(2));
xc=floor(N(1)/2)+1;yc=floor(N(2)/2)+1;
temp = (2*(X-xc)/N(1)).^2 + (2*(Y-yc)/N(2)).^2;
HKernel = fftshift(0.5*(1+cos(pi*temp)));
HF_data_rrf=ifft(ifft(fft(fft(data_rrf,[],1),[],2).*HKernel,[],1),[],2);


while (abs(LoopRatio-Ratio)/LoopRatio>0.05) & (step<1000)
  step=step+1;
   
     LipidRM1=inv( eye(N(end)) + Lambda * Lipid_stack_rf'*Lipid_stack_rf );
     %LipFree_rrf=reshape(reshape(HF_data_rrf,[],N(end))*LipidRM1,[N(1) N(2) N(3)]);
     LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*LipidRM1,[N(1) N(2) N(3)]);
     
    EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
    
   SkullEnergy=mean(EnergyMap(lipid_mask>0));
    BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0));
    
    Ratio=(BrainEnergy/SkullEnergy);
    % imagesc(Vol2Image(EnergyMap.*lipid_mask.*(EnergyMap>SkullEnergy)));drawnow;
    if (LoopRatio)>Ratio;
        Lambda=Lambda*(4)^(1/damp);
        %T2Star = T2Star +Nstep
    else
        damp=damp+1;
        Lambda=Lambda/(4)^(1/damp);
        % T2Star = T2Star -Nstep
    end
end
LipidProj1=IOP-LipidRM1;
%LipFree1_rrf=LipFree_rrf;

Nbasis=16;

RadiusErode=round((size(ReconParams.BrainMask2D,1)+size(ReconParams.BrainMask2D,2))/16);
CenterBrain=1+imerode(round(ReconParams.BrainMask2D),offsetstrel(ones(2*RadiusErode+1)));
%CenterBrain=ReconParams.BrainMask2D;

BrainMask_rrf=repmat(CenterBrain,[1 1 N(end)]);
%Brain_stack_rf=reshape(HF_data_rrf(BrainMask_rrf>0),[],N(end))*LipidRM1;
Brain_stack_rf=reshape(data_rrf(BrainMask_rrf>0),[],N(end))*LipidRM1;
  
[~,Sb,Vb] = svd(Brain_stack_rf,'econ');%0);
MetabProj=Vb(:,1:Nbasis)*Vb(:,1:Nbasis)';
MetabRM=IOP-MetabProj;
    figure
    plot(1:1056,abs(Sb(1,1)*Vb(:,1)),1:1056,abs(Sb(2,2)*Vb(:,2)));drawnow

Lambda=1;
damp=1;
step=0;
Ratio=999;
LoopRatio=SVratio;
while (abs(LoopRatio-Ratio)/LoopRatio>0.05) & (step<100)
    step=step+1;
     
    Lipid_stack_rf=reshape(LipData_rrf(LipidMask_rrf>0),[],N(end));%*MetabRM;
   % LipidRM=inv( eye(N(end)) + (MetabProj+MetabRM*Lambda * Lipid_stack_rf'*Lipid_stack_rf) );
     LipidRM=inv( eye(N(end)) + (Lambda * Lipid_stack_rf'*Lipid_stack_rf) );
    
    LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*(MetabProj+MetabRM*LipidRM),[N(1) N(2) N(3)]);
    EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
    SkullEnergy=mean(EnergyMap(lipid_mask>0));
    BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0));
    Ratio=(BrainEnergy/SkullEnergy);
    
        if (LoopRatio)>Ratio;
        Lambda=Lambda*(4)^(1/damp);
        %T2Star = T2Star +Nstep
    else
        damp=damp+1;
        Lambda=Lambda/(4)^(1/damp);
        % T2Star = T2Star -Nstep
    end
    
end

Ratio=Ratio;
LipidProj=IOP-(MetabProj+MetabRM*LipidRM);


%%%%%
%Checking with original code
%{
count = 0;
N = size(data_rrf);
[~,high_bnd_L]=min(abs(ReconParams.LipidMaxPPM - ReconParams.ppm));
[~,low_bnd_L]=min(abs(ReconParams.LipidMinPPM  - ReconParams.ppm));

LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
Lipid_stack_rf=reshape(LipData_rrf(LipidMask_rrf>0),[],N(end));

LipidProj=zeros(N(end),N(end));
LipidRM=zeros(N(end),N(end));
IOP=diag(ones(N(end),1));

norm_Data= norm(data_rrf(:));

Lambda=1E5*norm_Data^(-2);
Ratio=999;
LoopRatio=SVratio;
step=0;
while (abs(LoopRatio-Ratio)/LoopRatio>0.05) & (step<100)
  step=step+1;
%      LipidProj=([zeros(size(DataClean_rf))',Lambda*Lipid_stack_rf'/norm(Lipid_stack_rf)]...
%          /[DataClean_rf'/norm(DataClean_rf),Lambda*Lipid_stack_rf'/norm(Lipid_stack_rf)])';
%      LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3) N(4)]);
%    
     LipidRM2=inv( eye(N(end)) + Lambda * Lipid_stack_rf'*(Lipid_stack_rf) );
     %LipidRM=inv( eye(N(end)) + Lambda * (Lipid_stack_rf.*abs(Lipid_stack_rf).^(1/2))'*(Lipid_stack_rf) );
     LipFree2_rrf=reshape(reshape(data_rrf,[],N(end))*LipidRM2,[N(1) N(2) N(3)]);
     
    EnergyMap=sum(abs(LipFree2_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
    SkullEnergy=mean(EnergyMap(lipid_mask>0));
    BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0));
    %SkullEnergy=mean(EnergyMap((EnergyMap.*lipid_mask)>quantile(EnergyMap(lipid_mask>0),0.95)));
    %BrainEnergy=mean(EnergyMap((EnergyMap.*ReconParams.BrainMask2D)>quantile(EnergyMap(ReconParams.BrainMask2D>0),0.95)));

    Ratio=(BrainEnergy/SkullEnergy);
    % imagesc(Vol2Image(EnergyMap.*lipid_mask.*(EnergyMap>SkullEnergy)));drawnow;
    if (LoopRatio)>Ratio;
        Lambda=Lambda*1.2;
        %T2Star = T2Star +Nstep
    else
        Lambda=Lambda/1.5;
        % T2Star = T2Star -Nstep
    end
end

BrainMask_rrf=repmat(ReconParams.BrainMask2D,[1 1 N(end)]);
Brain_stack_rf=reshape(data_rrf(BrainMask_rrf>0),[],N(end))*LipidRM1;
[~,Sl,Vl] = svd(Brain_stack_rf,'econ');%0);
Brain_stack_rf=reshape(data_rrf(BrainMask_rrf>0),[],N(end))*LipidRM2;
[~,Sr,Vr] = svd(Brain_stack_rf,'econ');%0);
Brain_stack_rf=reshape(data_rrf(BrainMask_rrf>0),[],N(end))*(MetabProj+MetabRM*LipidRM);
[~,Sm,Vm] = svd(Brain_stack_rf,'econ');%0);
%}

end

