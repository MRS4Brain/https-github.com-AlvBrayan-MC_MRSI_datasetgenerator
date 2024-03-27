function [LipidProj, Slipid, Nbasis,Lipids ] = make_SVD_LipidBasis( lipid_rrf,data_rrf, lipid_mask, SVratio ,ReconParams,Nbasis)
%GET_LIPIDBASIS Summary of this function goes here
%   Detailed explanation goes here


count = 0;
N = size(data_rrf);
IOP=diag(ones(N(end),1));
Lipid_Basis = zeros(N(3), sum(lipid_mask(:)));

[~,high_bnd_L]=min(abs(ReconParams.LipidMaxPPM - ReconParams.ppm));
[~,low_bnd_L]=min(abs(ReconParams.LipidMinPPM  - ReconParams.ppm));

% SVD basis of the skull lipids 

LipidMask_rrf=repmat(lipid_mask,[1 1 N(end)]);
[~,Sorig,Vorig] = svd(reshape(lipid_rrf(LipidMask_rrf>0),[],N(end)),'econ');%0);

if(nargin==5)% no Nbasis given
	Nbasis=64;%32
	Nstep=32;%16
	while Nstep>=1
	    
	    LipidProj=Vorig(:,1:Nbasis)*Vorig(:,1:Nbasis)';
	    
	    LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3)]);
	    EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
	    
	    SkullEnergy=mean(EnergyMap(lipid_mask>0));
	    BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0));
	    %BrainEnergy/SkullEnergy;
	    if SkullEnergy>(BrainEnergy/SVratio);
		Nbasis=Nbasis+Nstep;
	    else
		Nbasis=Nbasis-Nstep;
	    end
	    Nstep=Nstep/2 ;
	end

	if Nbasis>ReconParams.L2SVDparams.NBasisMax
	    Nbasis=ReconParams.L2SVDparams.NBasisMax;
	end
	if Nbasis<ReconParams.L2SVDparams.NBasisMin
	    Nbasis=ReconParams.L2SVDparams.NBasisMin;
	end

elseif (nargin==6)
	LipidProj=Vorig(:,1:Nbasis)*Vorig(:,1:Nbasis)';
	    
	LipFree_rrf=reshape(reshape(data_rrf,[],N(end))*(IOP-LipidProj),[N(1) N(2) N(3)]);
	EnergyMap=sum(abs(LipFree_rrf(:,:,low_bnd_L:high_bnd_L)).^2,3);
	    
	SkullEnergy=mean(EnergyMap(lipid_mask>0));
	BrainEnergy=mean(EnergyMap(ReconParams.BrainMask2D>0));


end	
Nbasis = Nbasis
Lipids=Vorig(:,1:Nbasis);
Slipid=Sorig(1:Nbasis,1:Nbasis);

LipidProj=Lipids*Lipids';
end

