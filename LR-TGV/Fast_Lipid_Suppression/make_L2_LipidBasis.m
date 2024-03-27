function [LipidProj,Lambda]  = make_L2_LipidBasis( LipData_rrf, data_rrf,lipid_mask ,SVratio,ReconParams)
%GET_LIPIDBASIS Summary of this function goes here
%   Detailed explanation goes here

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
     LipidRM=inv( eye(N(end)) + Lambda * Lipid_stack_rf'*(Lipid_stack_rf) );
     %LipidRM=inv( eye(N(end)) + Lambda * (Lipid_stack_rf.*abs(Lipid_stack_rf).^(1/2))'*(Lipid_stack_rf) );
     
    LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*LipidRM,[N(1) N(2) N(3)]);
     
    EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
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
LipidProj=IOP-LipidRM;
step=step;
Lambda=Lambda;

%{
BrainMask_rrf=repmat(ReconParams.BrainMask2D,[1 1 N(end)]);
%Brain_stack_rf=reshape(data_rrf(BrainMask_rrf>0),[],N(end))*LipidRM;
Brain_stack_rf=reshape(HF_data_rrf(BrainMask_rrf>0),[],N(end))*LipidRM;
[~,Sl,Vl] = svd(Brain_stack_rf,'econ');%0);
figure; plot(1:1056,abs(Vl(:,1:2)*Sl(1:2,1:2)));drawnow
%}
end

