function Recon_OriginalSize_Data_trr = MRSI_LRTGV(Metab_matrix,B0drift_map,model_mask,Param_Filename,Data_name,mu_tv,PercentThres,WaterSuppComp,Comp,UndersamplingF)

% Some system variables initialization ...
[SCRIPT_DIR, ~, ~]=fileparts(mfilename('fullpath'));
addpath(genpath (SCRIPT_DIR));
CURRENT_DIR=pwd;
warning('off','MATLAB:DELETE:FileNotFound')

% Storing reconstruction input in a single structure:mrsiReconParams
mrsiReconParams=struct('NameData',Data_name,'mu_tv',mu_tv,'modelOrder',Comp);
mrsiReconParams.Log_Dir = ['Log_Files_', Data_name];
mrsiReconParams.Results_Dir = ['Results_',Data_name];

% Read the parameters written in the parameter file: Param_Filename
run(Param_Filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mrsiReconParams.InitialB0MapFromWater= 1; %yes / no, Compute a initial B0 map from water signal (generally 1 but might be incorrect in phantoms)

%mrsiReconParams.AcqDelay = 0.0013; %0 for a Spin-echo. For FID-MRSI sequence, put the acquisition delay in second. First missing point of the FID will be predicted to correct the 1st order phase
mrsiReconParams.AcqDelay = 0; %0 for a Spin-echo. For FID-MRSI sequence, put the acquisition delay in second. First missing point of the FID will be predicted to correct the 1st order phase


mrsiReconParams.MinPPM   = -10; % begining of the frequency window for low-rank Recon. Put to -1000 for full window
mrsiReconParams.MaxPPM   = 10; % End of the frequency window for low-rank Recon.Put to 1000 for full window
%mrsiReconParams.MaxPPM   = 1000.0; % End of the frequency window for low-rank Recon.Put to 1000 for full window


mrsiReconParams.DoHomogeneityCorrection= 0; %yes / no % use the Water signal to uniform the metabolite signal during reconstruction. This factor is also applied to the resulting Water signal so the final quantification is not affected by this option  .DoHomogeneityCorrection leads to instabilities with high B1 inhomogeneity and enhance lowest slice signal by strong factor. Not recommended
mrsiReconParams.HammingFidFilter=0;% Filter applied on the recon data fidelity terms. 0 - no filter , 1 max filtering. It helps the convergence in case of very low SNR Data

mrsiReconParams.Threshold_LipMask=0.5; % Threshold to reduce the size of the skull mask for each coil indivudually (Used for lipid suppression only)

mrsiReconParams.UndersamplingF = UndersamplingF; % Retrospective undesampling (for acceleration simulation purpose)

mrsiReconParams.ZeroPaddingF = 0;% Pad the end of the data: 1 add 100%, 0.5 add 50%, 0 no padding .(LCModel works somehow better with padded data)

%Parameters for the lipid suppression 
mrsiReconParams.L2SVDparams.PercentThres = PercentThres; % strength of the Lipid suppresion
mrsiReconParams.L2SVDparams.Nfit = Comp;%
mrsiReconParams.L2SVDparams.NBasisMax =128;%Max allowed components for Lipid suppression
mrsiReconParams.L2SVDparams.NBasisMin =5;%Min allowed components for Lipid suppression

mrsiReconParams.LipidMinPPM=-3;%Range of lipid in ppm (used to compute lipid contamination in lipid suppression)
mrsiReconParams.LipidMaxPPM=0.0;%Range of lipid in ppm (used to compute lipid contamination in lipid suppression)
mrsiReconParams.ESPIRIT_kernel=[6,6]; % Size of the ESPIRIT kernel for coil sensitivity profile computation

%Parameters for the HSVD of Water
mrsiReconParams.FiltParam.Water_minFreq=-NMRFreq*2;%hz %Min freq for water removal by HSVD (0Hz = 4.7ppm)
mrsiReconParams.FiltParam.Water_maxFreq=NMRFreq/2;%hz %Max freq for water removal by HSVD (0Hz = 4.7ppm)
mrsiReconParams.FiltParam.Comp=WaterSuppComp; %Nb of component for the HSVD water removal (advised: 16 at 3T and 32 at 7T)

mrsiReconParams.GaussianSigma=3;% The kernel width for B0 fieldmap correction

mrsiReconParams.LRTGVModelParams.Maxit= 1500; % Maximum number of iteration
mrsiReconParams.LRTGVModelParams.Minit= 600; % Minimum number of iteration
mrsiReconParams.LRTGVModelParams.check_it=150; % Step where the iterative recon check the different convergence parameters (output in log file)
mrsiReconParams.LRTGVModelParams.Plot_it=2000; % Step where the iterative recon make plot of the temporary results (output in Log_Files.../LowRankTGV_Recon/)
mrsiReconParams.LRTGVModelParams.CorrB0Map_it=25; % Step where the iterative recon estimatea a correction to the B0 Field Map
mrsiReconParams.LRTGVModelParams.CorrB0Map_Maxcount=100;% Max amount of  B0 Field Map correction
mrsiReconParams.LRTGVModelParams.CorrB0Map_MaxPPM=-0.5;
mrsiReconParams.LRTGVModelParams.CorrB0Map_MinPPM=-3.5;
mrsiReconParams.LRTGVModelParams.Orthogonalize_it=25;% Step where the iterative recon forces orthogonality in the temporaal & spatial component
mrsiReconParams.LRTGVModelParams.SpecItFact=2;
mrsiReconParams.LRTGVModelParams.reduction=100; %100 invivo %1000 phantom % Nb of iterations for soft introduction of the regularization in the recon
mrsiReconParams.LRTGVModelParams.min_SpectStep=1/32;% minimum spectral update step
mrsiReconParams.LRTGVModelParams.max_SpectStep=1/2;%1/2 invivo, 1/4 if diverge 1/8 phantom % maximum spectral update step
mrsiReconParams.LRTGVModelParams.min_taup=1/32;% minimal spatial update step
mrsiReconParams.LRTGVModelParams.max_taup=1/8; % %1/8 invivo , 1/16 if diverge, 1/16 for Synthetic Data, 1/32 phantom  % maximal spatial update step

%Make directory for log file
if ~exist(mrsiReconParams.Log_Dir);mkdir(mrsiReconParams.Log_Dir);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and filter the Field map Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(makeSectionDisplayHeader('Starting Data Preprocessing'));

%These masks will be determined later 
BrainMask=ones(MatSize);
ImMask=ones(MatSize);
SkMask=0; % SkullMask

mrsiReconParams.BrainMask        = BrainMask;
mrsiReconParams.ImMask           = ImMask;
mrsiReconParams.SkMask           = SkMask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load MRSI Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_tkk=Metab_matrix;
mrProt.samplerate=Bandwidth;
mrProt.DwellTime = 1.0/mrProt.samplerate;
mrsiReconParams.mrProt = mrProt;
mrsiReconParams.mrProt.VSize = FidPoints;
mrsiReconParams.mrProt.NMRFreq = NMRFreq;
NRG=squeeze(sum(abs(raw_tkk).^2,1));
mrsiReconParams.kmask = log(NRG)>-20; %must be centered in (1,1,1)

raw_ctkk = reshape(raw_tkk,[1 size(raw_tkk)]); % only one coil in the data (SO FAR)
clear raw_tkk;
Size_data=size(raw_ctkk);

Nb_Coil=Size_data(1);
mrsiReconParams.mrsiData = single(raw_ctkk);
clear raw_ctkk


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Anatomical masks at the mrsiData size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Size_data=size(mrsiReconParams.mrsiData);

mrsiReconParams.BrainMask2D=imresize(mean(mrsiReconParams.BrainMask,3),[Size_data(3),Size_data(4)]);
mrsiReconParams.BrainMask2D=round(mrsiReconParams.BrainMask2D/max(mrsiReconParams.BrainMask2D(:)));

mrsiReconParams.ImMask2D=imresize(mean(mrsiReconParams.ImMask,3),[Size_data(3),Size_data(4)]);
mrsiReconParams.ImMask2D=round(mrsiReconParams.ImMask2D/max(mrsiReconParams.ImMask2D(:)));

mrsiReconParams.SkMask2D=mrsiReconParams.ImMask2D-mrsiReconParams.BrainMask2D;
mrsiReconParams.SkMask2D(mrsiReconParams.SkMask2D<0)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute usefull variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDims = ndims(mrsiReconParams.mrsiData)-2;
Size_data=size(mrsiReconParams.mrsiData);
mrsiReconParams.spa_index = repmat({':'}, 1, nDims);
mrsiReconParams.ppm=(-4.7+((1:mrsiReconParams.mrProt.VSize)*mrsiReconParams.mrProt.samplerate/(mrsiReconParams.mrProt.VSize*mrsiReconParams.mrProt.NMRFreq)));
[~,mrsiReconParams.MaxPPM_pt]=min(abs(mrsiReconParams.MaxPPM  - mrsiReconParams.ppm));
[~,mrsiReconParams.MinPPM_pt]=min(abs(mrsiReconParams.MinPPM  - mrsiReconParams.ppm));

mrsiReconParams.NbMissingPoints=round(mrsiReconParams.mrProt.samplerate*mrsiReconParams.AcqDelay);% Number of missing points at the beginning of the FID


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply the UnderSampling pattern k-space mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mrsiReconParams.Originalkmask=mrsiReconParams.kmask;
if mrsiReconParams.UndersamplingF < 1
    [US_MASK]=Make_undersampled_mask(1.0,mrsiReconParams.UndersamplingF,mrsiReconParams.kmask,6);%Weighted Random distribution
    mrsiReconParams.kmask=mrsiReconParams.kmask.*US_MASK;
end

s=['./',mrsiReconParams.Log_Dir,'/','Kmasks_Original_US_CS_',Data_name,'.ps'];
if exist(s);delete(s);end
figs=figure('visible','off');
imagesc(fftshift(mrsiReconParams.Originalkmask));%,[ 0, 10*mean(image2plot(:))] );
title('Original K-mask');
print(figs, '-append', '-dpsc2', s);
figs=figure('visible', 'off');
imagesc(fftshift(mrsiReconParams.kmask));
title('Random undersampled K-mask');
print(figs, '-append', '-dpsc2', s);
close all;

if mrsiReconParams.UndersamplingF < 1
    fprintf(makeSectionDisplayHeader('Apply the Under-sampling pattern k-space mask ...'));
    kmask_ctkk=permute(repmat(mrsiReconParams.kmask,[1,1,Size_data(1),Size_data(2)]),[3,4,1,2]);
    mrsiReconParams.mrsiData=mrsiReconParams.mrsiData.*kmask_ctkk;
    clear kmask_ctkk   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Roemers Weights and Make Anatomical masks at the mrsiData size
% & compute SENSE operator (WITH ONE COIL ELEMENT, THEIR ARE JUST =1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf( 'Compute Roemers Weights & Anatomical Masks...\n');
mrsiReconParams.RoemersCoefs = ones(MatSize);
mrsiReconParams.SENSE = ones([1 MatSize]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate anatomical mask based on the water signal energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

 mrsiReconParams.BrainMask2D=model_mask;

 mrsiReconParams.SkMask2D =  mrsiReconParams.ImMask2D -  mrsiReconParams.BrainMask2D;
 mrsiReconParams.SkMask2D(mrsiReconParams.SkMask2D (:)<0)=0;
 
 mrsiReconParams.WaterFreqMap = B0drift_map;


if ~mrsiReconParams.InitialB0MapFromWater
	mrsiReconParams.WaterFreqMap=0*mrsiReconParams.WaterFreqMap; % set intial FreqMap to zero
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Uniformization Factor and Correct Water Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mrsiReconParams.DoHomogeneityCorrection 
        fprintf( 'Compute Homogeneity Correction map...\n');
	mrsiReconParams.HomCorr = abs(mrsiReconParams.Water_rr);
	
	mrsiReconParams.HomCorr = mrsiReconParams.HomCorr/mean(mrsiReconParams.HomCorr(mrsiReconParams.BrainMask2D>0));
	mrsiReconParams.HomCorr(mrsiReconParams.ImMask2D==0) = 1;
	mrsiReconParams.Water_trr = mrsiReconParams.Water_trr./reshape(mrsiReconParams.HomCorr,[ 1 size(mrsiReconParams.HomCorr)]);
else
	mrsiReconParams.HomCorr = ones(size(mrsiReconParams.WaterFreqMap));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restrict MRSI data to the desired ppm span & bB0map correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Time_ctkk=permute(repmat(([0 :(Size_data(2)-1)]'/mrsiReconParams.mrProt.samplerate),[1 Size_data(1) Size_data(3) Size_data(4)]),[2 1 3 4]);
% OriginalFullData_cfkk=fft(mrsiReconParams.mrsiData.*exp(2*pi*1i*Time_ctkk*mrsiReconParams.mrProt.NMRFreq),[],2); %Shift of 1 ppm
% clear Time_ctkk
OriginalFullData_cfkk=fft(mrsiReconParams.mrsiData,[],2);

mrsiData_ckkf=permute(OriginalFullData_cfkk(:,1:mrsiReconParams.MaxPPM_pt,:,:), [1,3,4,2]);
mrsiData_ckkf(:,:,:,1:mrsiReconParams.MinPPM_pt)=0;
mrsiData_ckkt=ifft(mrsiData_ckkf,[],4);

SENSE_ctrr=permute(repmat(mrsiReconParams.SENSE,[1 1 1 size(OriginalFullData_cfkk,2)]),[1,4,2,3]);
OriginalFullData_trr= conj(SENSE_ctrr).*ifft(ifft(ifft(OriginalFullData_cfkk,[],2),[],3),[],4);
OriginalFullData_trr =squeeze(sum(OriginalFullData_trr,1));

SENSE_ctrr=permute(repmat(mrsiReconParams.SENSE,[1 1 1 size(mrsiData_ckkt,4)]),[1,4,2,3]);
OriginalData_trr= conj(SENSE_ctrr).*ifft(ifft(ifft(OriginalFullData_cfkk(:,1:mrsiReconParams.MaxPPM_pt,:,:),[],2),[],3),[],4);
OriginalData_trr =squeeze(sum(OriginalData_trr,1));

fprintf( 'No Hamming Filter is computed neither used...\n');
mrsiReconParams.HKernel=ones(Size_data(3),Size_data(4));

clear FullmrsiData_cfkk OriginalFullData_cfkk SENSE_ctrr mrsiData_ckkf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOW RANK TGV Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(makeSectionDisplayHeader('TGV regularization...'));

OrigName=mrsiReconParams.NameData;

[mrsiReconParams.ReconDataShort_trr,mrsiReconParams.U_rrc,mrsiReconParams.V_tc,mrsiReconParams.S] = LowRankTGV_RecurModel(mrsiData_ckkt,mrsiReconParams);

clear TempRes_rrt Recon_US Recon_V mrsiData_ckkt;

OriginalSize_ReconData_frr=fft(OriginalFullData_trr,[],1); %frrr
OriginalSize_ReconData_frr(1:mrsiReconParams.MaxPPM_pt,:,:,:)=fft(mrsiReconParams.ReconDataShort_trr,[],1); %frr

Time_trr=repmat(([0 :(Size_data(2)-1)]'/mrsiReconParams.mrProt.samplerate),[1 Size_data(3) Size_data(4)]);
Recon_OriginalSize_Data_trr=ifft(OriginalSize_ReconData_frr,[],1);%trrr

clear OriginalSize_ReconData_frr;