function [FilteredData_trr,U_rrc,V_tc,S] = LowRankTGV_RecurModel(mrsiData_ckkt,mrsiReconParams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the Low-Rank TGV Recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist([mrsiReconParams.Log_Dir '/LowRankTGV_Recon']);mkdir([mrsiReconParams.Log_Dir '/LowRankTGV_Recon']);end

fprintf('Compute initial SVD Adjoint solution...\n');
SiOri=size(mrsiData_ckkt);
nDimsOri= ndims(mrsiData_ckkt);

NbT=SiOri(4);
FreqMap=mrsiReconParams.WaterFreqMap; % WaterFreqMap ;
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;

Time_rrt=permute(repmat(([0 :(NbT-1)]'/Fs),[1 SiOri(2) SiOri(3)]),[2,3,1]);
Freqshift_rrt=exp(-2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT])); %-2pi to go from k -> r
Freqshift_crrt= permute(repmat(Freqshift_rrt,[1 1 1 SiOri(1)]),[4 1 2 3]);
clear  Time_rrt Freqshift_rrt
Brainmask_crrt=permute(repmat(mrsiReconParams.BrainMask2D,[1 1 size(mrsiData_ckkt,1) size(mrsiData_ckkt,4)]),[3 1 2 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute initial condition for the iteration process (from a SVD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Uorig,Sorig,Vorig] = svd(reshape(fft(fft(ifft(ifft(mrsiData_ckkt,[],2),[],3).*Freqshift_crrt.*Brainmask_crrt,[],2),[],3) ,[],SiOri(4)),0);
if (size(Vorig,2)<mrsiReconParams.modelOrder)
	mrsiReconParams.modelOrder = size(Vorig,2);
end

V_tc=Vorig(:,1:mrsiReconParams.modelOrder);
S=Sorig(1:mrsiReconParams.modelOrder,1:mrsiReconParams.modelOrder);
U_ckkc=reshape(Uorig(:,1:mrsiReconParams.modelOrder), SiOri(1),SiOri(2),SiOri(3),[]);

clear Brainmask_crrt Freqshift_crrt

kmask=mrsiReconParams.kmask;

alpha = 1E-12*mrsiReconParams.mu_tv;
maxit = 500 ;          % use 1000 Iterations for optimal image quality        % usually there is no need to change this
Threshold=1E-12;
Brainmask_1rr=reshape(mrsiReconParams.BrainMask2D,[1 size(mrsiReconParams.BrainMask2D,1) size(mrsiReconParams.BrainMask2D,2)]);
Init_U_rrc=zeros(SiOri(2),SiOri(3),mrsiReconParams.modelOrder);
HomCorr_1rr = reshape(mrsiReconParams.HomCorr,[1 size(mrsiReconParams.HomCorr) ]);
fprintf('Start initial reconstruction...\n');
for k=1:mrsiReconParams.modelOrder
    [Init_U_rrc(:,:,k) e] = tgv2_l2_2D_multiCoil(U_ckkc(:,:,:,k),mrsiReconParams.SENSE.*Brainmask_1rr.*HomCorr_1rr, kmask, 2*alpha, alpha, maxit,Threshold);
end


VisualizeSpectral( V_tc,S, [ './',mrsiReconParams.Log_Dir,'/LowRankTGV_Recon/',mrsiReconParams.NameData,'_Initial']);
IFFTData=zeros(size(U_ckkc,2),size(U_ckkc,3),mrsiReconParams.modelOrder);
for comp=1:mrsiReconParams.modelOrder
    IFFTData(:,:,comp)= squeeze(sum(conj(mrsiReconParams.SENSE).*ifft(ifft(U_ckkc(:,:,:,comp),[],2),[],3),1));
end
VisualizeTGV( IFFTData ,Init_U_rrc,['./',mrsiReconParams.Log_Dir,'/LowRankTGV_Recon/',mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_Initial']);

clear  US_ckkc U_ckkc Uorig Vorig OriginalData_trr


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load initial condition in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_rrc=Init_U_rrc;

SizeVol = size(U_rrc);


fprintf(makeSectionDisplayHeader('Start the recurrence reconstruction...\n'));

alpha = mrsiReconParams.mu_tv;
maxit= mrsiReconParams.LRTGVModelParams.Maxit;
minit= mrsiReconParams.LRTGVModelParams.Minit;
Threshold=1;% threshold for stopping the process. Not used Anymore


[U_rrc,V_tc,S, costFunVal,DynamicFreqMap] = tgv2_l2_2D_multiCoil_LowRank_CombinedConv(mrsiData_ckkt,U_rrc,V_tc,S, 2*alpha, alpha, maxit,minit,mrsiReconParams,Threshold);

VisualizeTGV( Init_U_rrc ,U_rrc,['./',mrsiReconParams.Log_Dir,'/',mrsiReconParams.NameData,'_Final']);
VisualizeSpectral( V_tc,S, [ './',mrsiReconParams.Log_Dir,'/',mrsiReconParams.NameData,'_Final'])
s=[ './',mrsiReconParams.Log_Dir,'/',mrsiReconParams.NameData,'_Final_FrequencyMap.ps'];
figs=figure('visible','off');imagesc(DynamicFreqMap);
colorbar; title('Final Adaptative Frequency Map');print(figs, '-dpsc2',s);

fprintf('Recombination of the Data...\n');
TempRes_rrt=formTensorProduct(U_rrc,V_tc*S,nDimsOri-2);
FilteredData_trr=permute(TempRes_rrt,[ndims( TempRes_rrt),1:(ndims( TempRes_rrt)-1)]);

end
