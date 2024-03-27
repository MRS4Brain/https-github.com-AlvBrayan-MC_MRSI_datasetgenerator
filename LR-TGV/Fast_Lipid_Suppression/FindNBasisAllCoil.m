function  Nbasis  = FindNBasisAllCoil(Lipids_tkk,mrsiReconParams )
% mrsiReconParams.mrsiData dims: time-k-k


lipid_mask=mrsiReconParams.SkMask2D;
assignin("base","SkMask",mrsiReconParams.SkMask) %ADDED TO CHECK
MaxNbasis=64; %Must be a power of 2 

[~,high_bnd_L]=min(abs(mrsiReconParams.LipidMaxPPM - mrsiReconParams.ppm));
[~,low_bnd_L]=min(abs(mrsiReconParams.LipidMinPPM  - mrsiReconParams.ppm));

[NCoil, Nt, Nx,Ny] = size(mrsiReconParams.mrsiData);

assignin("base","N",[NCoil, Nt, Nx,Ny]) %ADDED TO CHECK

Lipid_SingSpect_fsc=zeros(Nt,MaxNbasis,NCoil);
assignin("base","Lipid_SingSpect_fsc",Lipid_SingSpect_fsc) %ADDED TO CHECK
LipidMask_rrf=zeros(Nx,Ny,Nt);


for coil = 1:NCoil
    CoilLipids_rrf=permute(fft(mrsiReconParams.SENSE(coil,:,:).*ifft(ifft(Lipids_tkk,[],2),[],3),[],1),[2,3,1]);
    
    Lipid_Vol=squeeze(sum(abs(CoilLipids_rrf),3));
    
    Signal=Lipid_Vol(:);
    if mrsiReconParams.Threshold_LipMask>=0
        LThr=mean(Signal)+mrsiReconParams.Threshold_LipMask*std(Signal);
    else
        LThr=0;
    end
    
    CoilLipid_mask=((Lipid_Vol.*lipid_mask)>LThr);
    assignin("base","CoilLipid_mask",CoilLipid_mask) %ADDED TO CHECK
    fprintf(['Coil=' num2str(coil) ,', Mean Signal =' num2str(mean(Signal)) ', STD Signal=' num2str(std(Signal)) ', LThr=' num2str(LThr) ', NumVox=' num2str(sum(CoilLipid_mask(:))) '\n']);
    
    LipidMask_rrf=repmat(squeeze(CoilLipid_mask),[1 1 Nt]);
    assignin("base","LipidMask_rrf",LipidMask_rrf) %ADDED TO CHECK
    Lipid_stack_rf=reshape(CoilLipids_rrf(LipidMask_rrf>0),[],Nt);
    assignin("base","Lipid_stack_rf",Lipid_stack_rf) %ADDED TO CHECK
    [~,Sorig,Vorig] = svd(Lipid_stack_rf,'econ');%0);

    assignin("base","Vorig",Vorig) %ADDED TO CHECK

    Lipid_SingSpect_fsc(:,:,coil) = Vorig(:,1:MaxNbasis);
    
end

Nbasis=MaxNbasis/2;
Nstep=MaxNbasis/4;
IOP=diag(ones(Nt,1));
LipidProj=0*IOP;
LipFree_rrf=zeros(Nx,Ny,Nt);
%CoilData_rrf=zeros(Nx,Ny,Nt);
CombLipFreeData_rrf=zeros(Nx,Ny,Nt);
lipid_mask=mrsiReconParams.SkMask2D;
BrainMask=mrsiReconParams.BrainMask2D.*~lipid_mask;

CoilData_crrf=permute(mrsiReconParams.mrsiData,[1,3,4,2]);%ckkt
CoilData_crrf=fft(ifft(ifft(CoilData_crrf,[],2),[],3),[],4);%crrf

LipSupStep=1;
figs=figure('visible', 'off');
while  Nstep>=1
    CombLipFreeData_rrf=CombLipFreeData_rrf*0;
    fprintf("coil = ");
    for coil = 1:NCoil
        fprintf([num2str(coil)+", "]);
                
        LipidProj=Lipid_SingSpect_fsc(:,1:Nbasis,coil)*Lipid_SingSpect_fsc(:,1:Nbasis,coil)';
    
        LipFree_rrf=reshape(reshape(CoilData_crrf(coil,:,:,:),[],Nt)*(IOP-LipidProj),[Nx Ny Nt]); 
        
        CombLipFreeData_rrf = CombLipFreeData_rrf + squeeze(conj(mrsiReconParams.SENSE(coil,:,:))).*LipFree_rrf;       
    end
   
    EnergyMap=sum(abs(CombLipFreeData_rrf(:,:,low_bnd_L:high_bnd_L).^2),3);
    %SkullEnergy=mean(EnergyMap(lipid_mask>0));% Not sensitive enough to inhomogenious Skull Lipids.... We want it homogenously removed
    %BrainEnergy=mean(EnergyMap(BrainMask>0));
    %SkullEnergy=max(EnergyMap(lipid_mask>0));
    %BrainEnergy=max(EnergyMap(BrainMask>0));
    SkullEnergy=quantile(EnergyMap(lipid_mask>0),0.97);
    BrainEnergy=quantile(EnergyMap(BrainMask>0),0.97);
    Ratio=(BrainEnergy/SkullEnergy);
   
    s=['./',mrsiReconParams.Log_Dir,'/LipidSuppression/Energy_After_LipSup_N',num2str(Nbasis),'.ps'];
    if exist(s);delete(s);end
 
    plotImage= EnergyMap;
    imagesc(plotImage,[0 , quantile(EnergyMap(:),0.97) ]);
    colormap default;colorbar;
    title('Energy Map');    
    print(figs, '-append', '-dpsc2', s);    

    plotImage= sqz(EnergyMap.*BrainMask) ;
    imagesc(plotImage,[0 , quantile(EnergyMap(:),0.97) ]);
    colormap default;colorbar;
    title(['Energy Map in Brain, Brain Energy = ',num2str(BrainEnergy)]);    
    print(figs, '-append', '-dpsc2', s);    

    plotImage= EnergyMap.*lipid_mask ;
    imagesc(plotImage,[0 , quantile(EnergyMap(:),0.97) ]);
    colormap default;colorbar;
    title(['Energy Map in Lipid Mask, Skull Energy = ',num2str(SkullEnergy)]);    
    print(figs, '-append', '-dpsc2', s);    

    fprintf(["Ratio = "+num2str(Ratio)+", Nbasis = "+num2str(Nbasis)+"\n"]);	    
	if (mrsiReconParams.L2SVDparams.PercentThres)>Ratio;
		Nbasis=Nbasis+Nstep;
	else
		Nbasis=Nbasis-Nstep;
    end
	Nstep=Nstep/2;

    LipSupStep=LipSupStep+1;
end

if Nbasis>mrsiReconParams.L2SVDparams.NBasisMax
   Nbasis=mrsiReconParams.L2SVDparams.NBasisMax;
end
if Nbasis<mrsiReconParams.L2SVDparams.NBasisMin
    Nbasis=mrsiReconParams.L2SVDparams.NBasisMin;
end

fprintf([" NBasis = "+num2str(Nbasis)+"\n"]);

end
