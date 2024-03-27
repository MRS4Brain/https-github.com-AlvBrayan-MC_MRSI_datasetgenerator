function MRSI_CSSENSELR_Recon(Metab_Filename,Water_Filename,Param_Filename,Data_name,mu_tv,PercentThres,WaterSuppComp,Comp,UndersamplingF)

% Some system variables initialization ...
[SCRIPT_DIR, ~, ~]=fileparts(mfilename('fullpath'));
addpath(genpath (SCRIPT_DIR));
CURRENT_DIR=pwd;
warning('off','MATLAB:DELETE:FileNotFound')

% Storing reconstruction input in a single structure:mrsiReconParams
mrsiReconParams=struct('NameData',Data_name,'mu_tv',mu_tv,'modelOrder',Comp,'Metab_Filename',Metab_Filename,'Water_Filename',Water_Filename);
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


mrsiReconParams.MinPPM   = -5.3; % begining of the frequency window for low-rank Recon. Put to -1000 for full window
mrsiReconParams.MaxPPM   = 2.0; % End of the frequency window for low-rank Recon.Put to 1000 for full window
%mrsiReconParams.MaxPPM   = 1000.0; % End of the frequency window for low-rank Recon.Put to 1000 for full window


mrsiReconParams.DoHomogeneityCorrection= 0; %yes / no % use the Water signal to uniform the metabolite signal during reconstruction. This factor is also applied to the resulting Water signal so the final quantification is not affected by this option  .DoHomogeneityCorrection leads to instabilities with high B1 inhomogeneity and enhance lowest slice signal by strong factor. Not recommended
mrsiReconParams.HammingFidFilter=0;% Filter applied on the recon data fidelity terms. 0 - no filter , 1 max filtering. It helps the convergence in case of very low SNR Data

mrsiReconParams.Threshold_LipMask=0.5; % Threshold to reduce the size of the skull mask for each coil indivudually (Used for lipid suppression only)

mrsiReconParams.UndersamplingF = UndersamplingF; % Retrospective undesampling (for acceleration simulation purpose)

mrsiReconParams.ZeroPaddingF = 0;% Pad the end of the data: 1 add 100%, 0.5 add 50%, 0 no padding .(LCModel works somehow better with padded data)

mrsiReconParams.L2SVDparams.PercentThres = PercentThres; % strength of the Lipid suppresion
mrsiReconParams.L2SVDparams.Nfit = Comp;%
mrsiReconParams.L2SVDparams.NBasisMax =128;%Max allowed components for Lipid suppression
%mrsiReconParams.L2SVDparams.NBasisMax =21;%Max allowed components for Lipid suppression
mrsiReconParams.L2SVDparams.NBasisMin =5;%Min allowed components for Lipid suppression
%mrsiReconParams.L2SVDparams.NBasisMin =21;%Min allowed components for Lipid suppression

mrsiReconParams.LipidMinPPM=-3;%Range of lipid in ppm (used to compute lipid contamination in lipid suppression)
mrsiReconParams.LipidMaxPPM=0.0;%Range of lipid in ppm (used to compute lipid contamination in lipid suppression)
mrsiReconParams.NbPtForWaterPhAmp= 5; %First pt of the time series taken for Amplitude and Phase calculation
mrsiReconParams.ESPIRIT_kernel=[6,6]; % Size of the ESPIRIT kernel for coil sensitivity profile computation


mrsiReconParams.FiltParam.Water_minFreq=-NMRFreq*2;%hz %Min freq for water removal by HSVD (0Hz = 4.7ppm)
mrsiReconParams.FiltParam.Water_maxFreq=NMRFreq/2;%hz %Max freq for water removal by HSVD (0Hz = 4.7ppm)
mrsiReconParams.FiltParam.Comp=WaterSuppComp; %Nb of component for the HSVD water removal (advised: 16 at 3T and 32 at 7T)

mrsiReconParams.GaussianSigma=3;% The kernel width for B0 fieldmap correction

mrsiReconParams.LRTGVModelParams.Maxit= 1500; % Maximum number of iteration
mrsiReconParams.LRTGVModelParams.Minit= 600; % Minimum number of iteration
mrsiReconParams.LRTGVModelParams.check_it=25; % Step where the iterative recon check the different convergence parameters (output in log file)
mrsiReconParams.LRTGVModelParams.Plot_it=100; % Step where the iterative recon make plot of the temporary results (output in Log_Files.../LowRankTGV_Recon/)
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


name_Preproc=['PreProcReconResults_HBRef_WS',num2str(WaterSuppComp) , 'C.mat'];

%Make directory for log file
if ~exist(mrsiReconParams.Log_Dir);mkdir(mrsiReconParams.Log_Dir);end


% if( exist(fullfile(CURRENT_DIR,name_Preproc),'file')~=2) % Enable the if statement to save intermediate data after water suppression.
    
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
    raw_tkk=Read_brucker_CSI(Metab_Filename,FidPoints,MatSize,NA);
    mrProt.samplerate=Bandwidth;
    mrProt.DwellTime = 1.0/mrProt.samplerate;
    mrsiReconParams.mrProt = mrProt;
    mrsiReconParams.mrProt.VSize = FidPoints;
    mrsiReconParams.mrProt.NMRFreq = NMRFreq;
    NRG=squeeze(sum(abs(raw_tkk).^2,1)); % Energy Mask
    mrsiReconParams.kmask = log(NRG)>-20; %must be centered in (1,1,1)
    
    raw_ctkk = reshape(raw_tkk,[1 size(raw_tkk)]); % only one coil in the data (SO FAR)
    clear raw_tkk;
    Size_data=size(raw_ctkk);
    
    Nb_Coil=Size_data(1);
    mrsiReconParams.mrsiData = single(raw_ctkk);
    clear raw_ctkk

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load MRSI Data of the Water
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    Water_tkk=Read_brucker_CSI(Water_Filename,FidPoints,MatSize,NA);
    Water_mrProt.samplerate=Bandwidth;
    Water_mrProt.DwellTime = 1.0/Water_mrProt.samplerate;
    mrsiReconParams.Water_mrProt = Water_mrProt;
    
    NRG=squeeze(sum(abs(Water_tkk).^2,1));
    mrsiReconParams.Water_kmask = log(NRG)>-20; %must be centered in (1,1,1)
    
    mrsiReconParams.Water_ctkk= reshape(Water_tkk,[1 size(Water_tkk)]);% only one coil in the data (SO FAR)
    clear Water_tkk 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shift data in Frequence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %metabolites
    Fs = mrsiReconParams.mrProt.samplerate;
    Time_t=[0 :(FidPoints-1)]'/Fs;
    Freqshift_1t=exp(2*pi*1i*Time_t*AcqFreqShift); %-2pi to go from k -> r
    Freqshift_1t= reshape(Freqshift_1t,[1 size(Freqshift_1t)]);
    
    mrsiReconParams.mrsiData = mrsiReconParams.mrsiData.*Freqshift_1t;
   
    %water
    Fs = mrsiReconParams.Water_mrProt.samplerate;
    Time_t=[0 :(FidPoints-1)]'/Fs;
    Freqshift_1t=exp(2*pi*1i*Time_t*WaterAcqFreqShift); %-2pi to go from k -> r
    Freqshift_1t= reshape(Freqshift_1t,[1 size(Freqshift_1t)]);
    
    mrsiReconParams.Water_ctkk = mrsiReconParams.Water_ctkk.*Freqshift_1t;
   

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
    % Print a Rough estimate of Water signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	s=['./',mrsiReconParams.Log_Dir,'/','RoughWaterSignalEnergy_',Data_name,'.ps'];
	if exist(s);delete(s);end
	Temp_cfrr=fftshift(ifft(ifft(ifft(mrsiReconParams.mrsiData,[],2),[],3),[],4),2); % c t k k
	HzpPt=mrsiReconParams.mrProt.samplerate/(mrsiReconParams.mrProt.VSize)	;
	Energy_rr=squeeze(sum(sum(abs(Temp_cfrr(:,round(end/2+mrsiReconParams.FiltParam.Water_minFreq/HzpPt):round(end/2+mrsiReconParams.FiltParam.Water_maxFreq/HzpPt),:,:)).^2,2),1));
	clear Temp_cfrr ;
	figs=figure('visible','off');
	imagesc(Energy_rr);colorbar;%,[ 0, 10*mean(image2plot(:))] );
	title('Water Energy Signal');
	print(figs, '-bestfit', '-dpsc2', s);
	close all;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do Water Suppression by HSVD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    if mrsiReconParams.FiltParam.Comp>0
        fprintf(makeSectionDisplayHeader('Starting Water Removal by HSVD...'));
        Filtered_ctkk=zeros(size(mrsiReconParams.mrsiData));
        
        fprintf([ 'Processing Coil ']);

	for coil=1:Size_data(1);
            fprintf([ ' ', num2str(coil), ',']);drawnow('update');
            [Filtered_ctkk(coil,:,:,:)] = WaterSuppression(squeeze(mrsiReconParams.mrsiData(coil,:,:,:)),mrsiReconParams.mrProt,        mrsiReconParams.FiltParam, mrsiReconParams.kmask);
        end
        fprintf('\n');
        mrsiReconParams.mrsiData= single( Filtered_ctkk);
        clear Filtered_ctkk
    end
    
%{
 % Enable the if statement to save intermediate data after water suppression.
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save initial- intermediate data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(makeSectionDisplayHeader(['Saving exisiting preprocessed data in file: ',fullfile(CURRENT_DIR,name_Preproc)]));
    save( fullfile(CURRENT_DIR,name_Preproc), 'mrsiReconParams','-v7.3');
    

else

    NotSaved_mrsiReconParams=mrsiReconParams;
    fprintf(makeSectionDisplayHeader(['Loading exisiting preprocessed data in file: ',fullfile(CURRENT_DIR,name_Preproc)]));
    load( fullfile(CURRENT_DIR,name_Preproc),'-mat');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Replace only the needed preprocessed data fields
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    NotSaved_mrsiReconParams.kmask = mrsiReconParams.kmask  ;
    
    NotSaved_mrsiReconParams.SkMask = mrsiReconParams.SkMask;
    NotSaved_mrsiReconParams.BrainMask = mrsiReconParams.BrainMask;
    NotSaved_mrsiReconParams.ImMask = mrsiReconParams.ImMask;
    
    NotSaved_mrsiReconParams.SkMask2D = mrsiReconParams.SkMask2D;
    NotSaved_mrsiReconParams.BrainMask2D = mrsiReconParams.BrainMask2D;
    NotSaved_mrsiReconParams.ImMask2D = mrsiReconParams.ImMask2D;
    
    if(isfield(mrsiReconParams,'Water_ctkk') )
        NotSaved_mrsiReconParams.Water_ctkk = mrsiReconParams.Water_ctkk;
        NotSaved_mrsiReconParams.Water_kmask = mrsiReconParams.Water_kmask;
        NotSaved_mrsiReconParams.Water_mrProt = mrsiReconParams.Water_mrProt;
        NotSaved_mrsiReconParams.Water_Filename = mrsiReconParams.Water_Filename;
    end
    if(isfield(mrsiReconParams,'WaterOfMRSIData_ctkk') )
        NotSaved_mrsiReconParams.WaterOfMRSIData_ctkk = mrsiReconParams.WaterOfMRSIData_ctkk;
    end
    if(isfield(mrsiReconParams,'HeadWater_ctkk') )
        NotSaved_mrsiReconParams.HeadWater_ctkk = mrsiReconParams.HeadWater_ctkk;
        NotSaved_mrsiReconParams.HeadWater_kmask = mrsiReconParams.HeadWater_kmask;
        NotSaved_mrsiReconParams.HeadWater_mrProt = mrsiReconParams.HeadWater_mrProt;
        NotSaved_mrsiReconParams.HeadWater_Filename = mrsiReconParams.HeadWater_Filename;
        NotSaved_mrsiReconParams.Water_ctkk = mrsiReconParams.HeadWater_ctkk;
        NotSaved_mrsiReconParams.Water_kmask = mrsiReconParams.HeadWater_kmask;
        NotSaved_mrsiReconParams.Water_mrProt = mrsiReconParams.HeadWater_mrProt;
        NotSaved_mrsiReconParams.Water_Filename = mrsiReconParams.HeadWater_Filename;
    end
    if(isfield(mrsiReconParams,'BodyWater_ctkk') )
        NotSaved_mrsiReconParams.BodyWater_ctkk = mrsiReconParams.BodyWater_ctkk;
        NotSaved_mrsiReconParams.BodyWater_kmask = mrsiReconParams.BodyWater_kmask;
        NotSaved_mrsiReconParams.BodyWater_mrProt = mrsiReconParams.BodyWater_mrProt;
        NotSaved_mrsiReconParams.BodyWater_Filename = mrsiReconParams.BodyWater_Filename;
    end
    
    NotSaved_mrsiReconParams.Metab_Filename = mrsiReconParams.Metab_Filename;
    NotSaved_mrsiReconParams.mrProt = mrsiReconParams.mrProt;
    NotSaved_mrsiReconParams.mrsiData = mrsiReconParams.mrsiData;
    
    mrsiReconParams = NotSaved_mrsiReconParams;
end

clear NotSaved_mrsiReconParams;
%}

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
 
mrsiReconParams.Water_trr=squeeze(  ifft(ifft(mrsiReconParams.Water_ctkk,[],3),[],4));

mrsiReconParams.Water_rr = squeeze(mean(mrsiReconParams.Water_trr(1:mrsiReconParams.NbPtForWaterPhAmp,:,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate anatomical mask based on the water signal energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Energy_rr=squeeze(abs(mrsiReconParams.Water_rr).^2);


 mrsiReconParams.BrainMask2D=Energy_rr>0.5*mean(Energy_rr(:));;
 %mrsiReconParams.BrainMask2D=Energy_rr>mean(Energy_rr(:));;
 mrsiReconParams.SkMask2D =  mrsiReconParams.ImMask2D -  mrsiReconParams.BrainMask2D;
 mrsiReconParams.SkMask2D(mrsiReconParams.SkMask2D (:)<0)=0;
 
 WaterFreqMap_rr=DetermineSingleLowFreq( mrsiReconParams.Water_trr,200,mrsiReconParams);
 WaterFreqMap_rr=mrsiReconParams.BrainMask2D.*WaterFreqMap_rr-mean(WaterFreqMap_rr(mrsiReconParams.BrainMask2D>0));%remove intercept and set Freq to 0 outside Brain
 mrsiReconParams.WaterFreqMap = WaterFreqMap_rr;
 
 VisualizeMasks( mrsiReconParams.Water_ctkk ,mrsiReconParams.ImMask2D, mrsiReconParams.BrainMask2D,mrsiReconParams.SkMask2D,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData])
 VisualizeKmasks( mrsiReconParams.mrsiData,mrsiReconParams.kmask, mrsiReconParams.Water_ctkk ,mrsiReconParams.Water_kmask,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData]);
 fprintf( 'Compute Sensitivity maps based on Head Water only...\n');


 VisualizeData( mrsiReconParams.mrsiData, mrsiReconParams.Water_ctkk ,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_Initial']);
 VisualizeCombinedData( mrsiReconParams.mrsiData, mrsiReconParams.Water_ctkk ,mrsiReconParams,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_Initial'])


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
	mrsiReconParams.HomCorr = ones(size(mrsiReconParams.Water_rr));
end

s=[ './',mrsiReconParams.Log_Dir,'/',mrsiReconParams.NameData,'_HomogeneityCorrection.ps'];
figs=figure('visible','off');imagesc(mrsiReconParams.HomCorr);
colorbar; title('Signal Homogeneity Estimated with Water Signal');print(figs, '-dpsc2',s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and apply Hann Filter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mrsiReconParams.HammingFidFilter>0

	Ham=mrsiReconParams.HammingFidFilter;
	fprintf( ['Compute Hamming Filter with param =',num2str(Ham) ,' for data fidelity term...\n']);
	[X,Y] = ndgrid(1:Size_data(3), 1:Size_data(4));
	xc=floor(Size_data(3)/2)+1;yc=floor(Size_data(4)/2)+1;
	temp = (2*(X-xc)/Size_data(3)).^2 + (2*(Y-yc)/Size_data(4)).^2;

	mrsiReconParams.HKernel = fftshift( (1-Ham/2)+Ham/2*cos(pi*temp) );
	clear X Y xc temp
else
	fprintf( 'No Hamming Filter is computed neither used...\n');
	mrsiReconParams.HKernel=ones(Size_data(3),Size_data(4));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do Lipid Supression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(mrsiReconParams.L2SVDparams.PercentThres>0)
    fprintf(makeSectionDisplayHeader('Starting Lipid Removal by Skull Lipid Projection...'));
    Lipids_ctkk = mrsiReconParams.mrsiData;
    SiSE=size(mrsiReconParams.SENSE);
    Lipids_tkk = sqz(ifft(ifft( sum( conj(reshape(mrsiReconParams.SENSE,[SiSE(1),1,SiSE(2),SiSE(3)])).*ifft(ifft(Lipids_ctkk,[],3),[],4) ,1) ,[],3),[],4) );

 if ~exist([mrsiReconParams.Log_Dir '/LipidSuppression']);mkdir([mrsiReconParams.Log_Dir '/LipidSuppression']);end
  
if(mrsiReconParams.L2SVDparams.PercentThres<2)
     Nbasis  = FindNBasisAllCoil(Lipids_tkk,mrsiReconParams )
else
    Nbasis  = round(mrsiReconParams.L2SVDparams.PercentThres)
end
    mrsiDataLr=zeros(size(mrsiReconParams.mrsiData));
    LipidProj_cff=zeros(Size_data(1),Size_data(2),Size_data(2));

    fprintf([ 'Processing Coil']);
    for coil = (1:Size_data(1))
        fprintf([ ' ', num2str(coil), ',']);drawnow('update');
        [mrsiDataLr(coil,:,:,:),~,~,Nbasis_c(coil),LipidProj_cff(coil,:,:)]   = ProjSVDLipidSuppression( squeeze(mrsiReconParams.mrsiData(coil,:,:,:)),squeeze(mrsiReconParams.mrsiData(coil,:,:,:)), mrsiReconParams ,['LipidSuppression/',mrsiReconParams.NameData,'_Coil',num2str(coil) ],Nbasis);
    end

    fprintf('\n');
    mrsiReconParams.L2SVDparams.Nbasis=Nbasis;
    mrsiReconParams.mrsiData=mrsiDataLr;
    clear mrsiDataLr LipidProj mrsiData LipidsParams mrsiData_Lipids Lipids_ctkk
    VisualizeData( mrsiReconParams.mrsiData, mrsiReconParams.Water_ctkk ,[mrsiReconParams.Log_Dir '/' mrsiReconParams.NameData '_AfterLipidSupp']);
    
else
    fprintf('No Lipid Suppression.\n');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restrict MRSI data to the desired ppm span & bB0map correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Time_ctkk=permute(repmat(([0 :(Size_data(2)-1)]'/mrsiReconParams.mrProt.samplerate),[1 Size_data(1) Size_data(3) Size_data(4)]),[2 1 3 4]);
OriginalFullData_cfkk=fft(mrsiReconParams.mrsiData.*exp(2*pi*1i*Time_ctkk*mrsiReconParams.mrProt.NMRFreq),[],2); % Empirical 1 ppm first order phase correction (CHECK k-space)
clear Time_ctkk

mrsiData_ckkf=permute(OriginalFullData_cfkk(:,1:mrsiReconParams.MaxPPM_pt,:,:), [1,3,4,2]); % Selection of the data to apply the LR denoising in the frequency domain
mrsiData_ckkf(:,:,:,1:mrsiReconParams.MinPPM_pt)=0;
mrsiData_ckkt=ifft(mrsiData_ckkf,[],4);

SENSE_ctrr=permute(repmat(mrsiReconParams.SENSE,[1 1 1 size(OriginalFullData_cfkk,2)]),[1,4,2,3]);
OriginalFullData_trr= conj(SENSE_ctrr).*ifft(ifft(ifft(OriginalFullData_cfkk,[],2),[],3),[],4);
OriginalFullData_trr =squeeze(sum(OriginalFullData_trr,1));

SENSE_ctrr=permute(repmat(mrsiReconParams.SENSE,[1 1 1 size(mrsiData_ckkt,4)]),[1,4,2,3]);
OriginalData_trr= conj(SENSE_ctrr).*ifft(ifft(ifft(OriginalFullData_cfkk(:,1:mrsiReconParams.MaxPPM_pt,:,:),[],2),[],3),[],4);
OriginalData_trr =squeeze(sum(OriginalData_trr,1));

clear FullmrsiData_cfkk OriginalFullData_cfkk SENSE_ctrr mrsiData_ckkf;

if ~exist(mrsiReconParams.Results_Dir);mkdir(mrsiReconParams.Results_Dir);end

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
mrsiReconParams.Recon_OriginalSize_Data_trr=ifft(OriginalSize_ReconData_frr,[],1);%trrr

clear OriginalSize_ReconData_frr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero padding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(makeSectionDisplayHeader('Zero padding & Phase correction...'));

VSize=round((mrsiReconParams.ZeroPaddingF+1)*size(OriginalFullData_trr,1))-mrsiReconParams.NbMissingPoints;
%VSize=round((mrsiReconParams.ZeroPaddingF+1)*size(OriginalFullData_trr,1))+mrsiReconParams.NbMissingPoints; %we concluded that we need to sum the orig points with the missing points
mrsiReconParams.ppm=(-4.7+((1:VSize)*mrsiReconParams.mrProt.samplerate/(VSize*mrsiReconParams.mrProt.NMRFreq)));

[~,mrsiReconParams.MaxPPM_pt]=min(abs(mrsiReconParams.MaxPPM - mrsiReconParams.ppm));
[~,mrsiReconParams.MinPPM_pt]=min(abs(mrsiReconParams.MinPPM - mrsiReconParams.ppm));

PadSizeSmall=(mrsiReconParams.MaxPPM_pt)-size(mrsiReconParams.ReconDataShort_trr,1);
PadSizeFull=round((mrsiReconParams.ZeroPaddingF)*size(OriginalFullData_trr,1))-mrsiReconParams.NbMissingPoints;
%PadSizeFull=0; %We don't intend to pad the data, however due to missing points we need to fix it to 0

BrainMask2D=mrsiReconParams.BrainMask2D;

mrsiReconParams.ReconDataShort_trr=padarray(mrsiReconParams.ReconDataShort_trr,PadSizeSmall,0,'post');
OriginalData_trr=padarray(OriginalData_trr,PadSizeSmall,0,'post');
ZeroPaddedV_tc=padarray(conj(mrsiReconParams.V_tc),PadSizeSmall,0,'post');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OriginalData_frr=fft(padarray(OriginalFullData_trr,PadSizeFull,0,'post'),[],1);

clear OriginalFullData_trr OriginalData_trr

FilteredData_frr=OriginalData_frr;
FilteredData_frr(1:size(mrsiReconParams.ReconDataShort_trr,1) ,:,:)=fft(mrsiReconParams.ReconDataShort_trr,[],1);

ZeroPaddedV_fc=zeros(size(OriginalData_frr,1),size(ZeroPaddedV_tc,2));
ZeroPaddedV_fc(1:size(ZeroPaddedV_tc,1) ,:)=fft(ZeroPaddedV_tc,[],1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add first missing points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Size_data=size(FilteredData_frr);
NaNData=[repmat(NaN,[mrsiReconParams.NbMissingPoints,Size_data(2)*Size_data(3)]); reshape(ifft(FilteredData_frr,[],1),Size_data(1),[]) ];       
FilteredData_frr = fft(reshape(fillgaps(double(NaNData)),[ size(NaNData,1) Size_data(2:3)] ),[],1);

Size_data=size(OriginalData_frr);
NaNData=[repmat(NaN,[mrsiReconParams.NbMissingPoints,Size_data(2)*Size_data(3)]); reshape(ifft(OriginalData_frr,[],1),Size_data(1),[]) ];       
OriginalData_frr = fft(reshape(fillgaps(double(NaNData)),[ size(NaNData,1) Size_data(2:3)] ),[],1);
  
Size_data=size(ZeroPaddedV_fc);
NaNData=[repmat(NaN,[mrsiReconParams.NbMissingPoints,Size_data(2)]); reshape(conj(ifft(ZeroPaddedV_fc,[],1)),Size_data(1),[]) ];       
ZeroPaddedV_fc = fft(reshape(fillgaps(double(NaNData)),[ size(NaNData,1) Size_data(2)] ),[],1);

mrsiReconParams.mrProt.VSize=size(FilteredData_frr,1);

mrsiReconParams.ppm=(-4.7+((1:mrsiReconParams.mrProt.VSize)*mrsiReconParams.mrProt.samplerate/(mrsiReconParams.mrProt.VSize*mrsiReconParams.mrProt.NMRFreq)));
mrsiReconParams.MaxPPM_pt=size(mrsiReconParams.ReconDataShort_trr,1);%min(abs(mrsiReconParams.MaxPPM- mrsiReconParams.ppm));
[~,mrsiReconParams.MinPPM_pt]=min(abs(mrsiReconParams.MinPPM - mrsiReconParams.ppm));
PaddedRange_max=mrsiReconParams.MaxPPM_pt;
PaddedRange_min=mrsiReconParams.MinPPM_pt;

Time_tc=repmat(([0 :(size(OriginalData_frr,1)-1)]'/mrsiReconParams.mrProt.samplerate),[1 size(ZeroPaddedV_tc,2)]);
mrsiReconParams.TimeComp=ifft(ZeroPaddedV_fc,[],1);%tc;

Size_data=size(OriginalData_frr);
Time_trr=repmat(([0 :(size(FilteredData_frr,1)-1)]'/mrsiReconParams.mrProt.samplerate),[1 Size_data(2) Size_data(3)]);

Recon_ZP_Data_trr=ifft( FilteredData_frr,[],1);%tre;
Original_RePhased_Data_trr=ifft( OriginalData_frr,[],1);%trr;
clear Time_trr Time_tc


PadSize=size(Recon_ZP_Data_trr,1)-size(mrsiReconParams.Water_trr,1);
mrsiReconParams.Water_trr=padarray(squeeze(mrsiReconParams.Water_trr),PadSize,0,'post');

mrsiReconParams.Recon_ZP_Data_trr = Recon_ZP_Data_trr;

mrsiReconParams.Original_RePhased_Data_trr=Original_RePhased_Data_trr;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save  data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reconResults= mrsiReconParams;
clear mrsiReconParam

reconResults=rmfield(reconResults, {'mrsiData','Water_ctkk'});

reconResults.MainRDA_Filename=[reconResults.Results_Dir,'/MRSI_CSSENSELR_Recon_',reconResults.NameData,'.rda'];
reconResults.OrigRDA_Filename=[reconResults.Results_Dir,'/MRSI_LipidandWaterRemoved_NoRecon_',reconResults.NameData,'.rda'];
reconResults.WaterRDA_Filename=[reconResults.Results_Dir,'/MRSI_Water_',reconResults.NameData,'.rda'];

if ~exist(reconResults.Results_Dir);mkdir(reconResults.Results_Dir);end

fprintf('Saving Data. \n');
save( [reconResults.Results_Dir,'/MRSI_MatlabReconResult_',Data_name,'.mat'], 'reconResults', '-v7.3');


