function [LipidRMOP,Lambda]  = make_WdL2_LipidRMOP( LipData_rrf, data_rrf,lipid_mask ,SVratio,ReconParams)
%GET_LIPIDBASIS Summary of this function goes here
%   Detailed explanation goes here

count = 0;
N = size(data_rrf);
IOP=eye(N(end));
[~,high_bnd_L]=min(abs(ReconParams.LipidMaxPPM - ReconParams.ppm));
[~,low_bnd_L]=min(abs(ReconParams.LipidMinPPM  - ReconParams.ppm));


[~,high_LipId]=min(abs(-1.9 - ReconParams.ppm));
[~,low_LipId]=min(abs(0  - ReconParams.ppm));
low_LipId=1056;


LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
[~,Sorig,Vorig] = svd(reshape(LipData_rrf(LipidMask_rrf>0),[],N(end)),'econ');%0);

lipid_Wdw_rrf=LipData_rrf(:,:,high_LipId:low_LipId);
LipidMask_rrf=repmat(lipid_mask,[1 1 size(lipid_Wdw_rrf,3)]);
[~,SWdw,VWdw] = svd(reshape(lipid_Wdw_rrf(LipidMask_rrf>0),[], size(lipid_Wdw_rrf,3)),'econ');%0);

Window=zeros(N(end),1);
Window(high_LipId:low_LipId)=1;
WindowOPSqr=diag(Window)*N(end)/(low_LipId-high_LipId);
WindowOP=IOP(1:N(end),high_LipId:low_LipId);


N1=4;N2=64;
LipidProj=Vorig(high_LipId:low_LipId,1:N1)*Vorig(:,1:N1)';
WdwLipidProj=VWdw(:,1:N2)*VWdw(:,1:N2)';
LipFree1_rrf=reshape(reshape(data_rrf,[],N(3))*(IOP-WindowOP*WdwLipidProj*LipidProj),[N(1) N(2) N(3)]);
Nbasis=N1;
LipidProj1=WindowOP*WdwLipidProj*LipidProj;


LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
Lipid_stack_rf=reshape(LipFree1_rrf(LipidMask_rrf>0),[],N(end));

LipidProj=zeros(N(end),N(end));
LipidRM=zeros(N(end),N(end));
IOP=diag(ones(N(end),1));

Window=zeros(N(end),1);
Window(high_LipId:low_LipId)=1;
WindowOP=diag(Window)*N(end)/(low_LipId-high_LipId);

norm_Data= norm(data_rrf(:));

Lambda=1000*norm_Data^(-2);
Ratio=999;
LoopRatio=SVratio;
step=0;
while (abs(LoopRatio-Ratio)/LoopRatio>0.05) & (step<100)
  step=step+1;
%      LipidProj=([zeros(size(DataClean_rf))',Lambda*Lipid_stack_rf'/norm(Lipid_stack_rf)]...
%          /[DataClean_rf'/norm(DataClean_rf),Lambda*Lipid_stack_rf'/norm(Lipid_stack_rf)])';
%      LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3) N(4)]);
%     
     LipidRM=inv( eye(N(end)) + Lambda * Lipid_stack_rf'*Lipid_stack_rf );
     %LipidProj=inv( eye(N(end))- Lambda*0*WindowOP - Lambda* Lipid_stack_rf'*Lipid_stack_rf );
     LipFree_rrf=reshape(reshape(LipFree1_rrf,[],N(end))*LipidRM,[N(1) N(2) N(3)]);
     %Lip_rrf=reshape(reshape(data_rrf,[],N(end))*LipidProj,[N(1) N(2) N(3)]);
    EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
    SkullEnergy=mean(EnergyMap(lipid_mask>0));
    BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0));
    %SkullEnergy=mean(EnergyMap((EnergyMap.*lipid_mask)>quantile(EnergyMap(lipid_mask>0),0.95)));
    %BrainEnergy=mean(EnergyMap((EnergyMap.*ReconParams.BrainMask2D)>quantile(EnergyMap(ReconParams.BrainMask2D>0),0.95)));

    Ratio=(BrainEnergy/SkullEnergy)
    % imagesc(Vol2Image(EnergyMap.*lipid_mask.*(EnergyMap>SkullEnergy)));drawnow;
    if (LoopRatio)>Ratio;
        Lambda=Lambda*1.2
        %T2Star = T2Star +Nstep
    else
        Lambda=Lambda/1.5
        % T2Star = T2Star -Nstep
    end
end




LipidProj=IOP-LipidRM;
step=step;
Lambda=Lambda;

LipidRMOP=LipidRM*LipidProj1;

LipFree2_rrf=LipFree_rrf;

LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
Lipid_stack_rf=reshape(LipData_rrf(LipidMask_rrf>0),[],N(end));

LipidProj=zeros(N(end),N(end));
LipidRM=zeros(N(end),N(end));
IOP=diag(ones(N(end),1));

Window=zeros(N(end),1);
Window(high_LipId:low_LipId)=1;
WindowOP=diag(Window)*N(end)/(low_LipId-high_LipId);

norm_Data= norm(data_rrf(:));

Lambda=1000*norm_Data^(-2);
Ratio=999;
LoopRatio=SVratio;
step=0;
while (abs(LoopRatio-Ratio)/LoopRatio>0.05) & (step<100)
  step=step+1;
%      LipidProj=([zeros(size(DataClean_rf))',Lambda*Lipid_stack_rf'/norm(Lipid_stack_rf)]...
%          /[DataClean_rf'/norm(DataClean_rf),Lambda*Lipid_stack_rf'/norm(Lipid_stack_rf)])';
%      LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3) N(4)]);
%     
     LipidRM=inv( eye(N(end)) + Lambda * Lipid_stack_rf'*Lipid_stack_rf );
     %LipidProj=inv( eye(N(end))- Lambda*0*WindowOP - Lambda* Lipid_stack_rf'*Lipid_stack_rf );
     LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*LipidRM,[N(1) N(2) N(3)]);
     %Lip_rrf=reshape(reshape(data_rrf,[],N(end))*LipidProj,[N(1) N(2) N(3)]);
    EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
    SkullEnergy=mean(EnergyMap(lipid_mask>0));
    BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0));
    %SkullEnergy=mean(EnergyMap((EnergyMap.*lipid_mask)>quantile(EnergyMap(lipid_mask>0),0.95)));
    %BrainEnergy=mean(EnergyMap((EnergyMap.*ReconParams.BrainMask2D)>quantile(EnergyMap(ReconParams.BrainMask2D>0),0.95)));

    Ratio=(BrainEnergy/SkullEnergy)
    % imagesc(Vol2Image(EnergyMap.*lipid_mask.*(EnergyMap>SkullEnergy)));drawnow;
    if (LoopRatio)>Ratio;
        Lambda=Lambda*1.2
        %T2Star = T2Star +Nstep
    else
        Lambda=Lambda/1.5
        % T2Star = T2Star -Nstep
    end
end


end

