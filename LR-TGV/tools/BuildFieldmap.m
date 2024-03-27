function [fieldmap MagVol]=BuildFieldmap(fNameMagImg,fNamePhaseImg)
%% ************************************************************************
% Load the phase and magnitude images into the workspace
%**************************************************************************

%fNameMagImg   = '/Users/Jeff/HUG Siemens 3T Data/EXPERIMENT_1H_05_11_2015/gre_field_mapping_240_492_12_184216';
%fNamePhaseImg = '/Users/Jeff/HUG Siemens 3T Data/EXPERIMENT_1H_05_11_2015/gre_field_mapping_240_492_14_184217';

%fNameMagImg   = '/data/MR_data/Prisma/Data_Klauser/Data-2015-12-15-AK/P-LOGAN_PHANTOM_15_12_15-19_35_45-STD-1_3_12_2_1107_5_2_43_67018/RECHERCHE_A_JOUR_PRISMA_ANTOINE_20151215_193633_297000/GRE_FIELD_MAPPING_240_492_0024/20151215_193633grefieldmapping240492s024a1001.nii.gz';
%fNamePhaseImg = '/data/MR_data/Prisma/Data_Klauser/Data-2015-12-15-AK/P-LOGAN_PHANTOM_15_12_15-19_35_45-STD-1_3_12_2_1107_5_2_43_67018/RECHERCHE_A_JOUR_PRISMA_ANTOINE_20151215_193633_297000/GRE_FIELD_MAPPING_240_492_0026/20151215_193633grefieldmapping240492s026a2001.nii.gz';


magImg   = load_nifti(fNameMagImg);
phaseImg = load_nifti(fNamePhaseImg);

%figure
%imagesc(squeeze(magImg.vol(:,:,1)));
%figure
%imagesc(squeeze(phaseImg.vol(:,:,1)));
% *************************************************************************
% Compute the phase maps
%**************************************************************************

% Get the bit depth from the corresponding DICOM headers
bitDepth = 12;

% Rescale the phase image (assuming Siemens data scaled between 1 and 4096)
%phaseImg.vol = (phaseImg.vol / 2^(bitDepth - 1)) - 1;
max_ph=max(max(max(phaseImg.vol)));
min_ph=min(min(min(phaseImg.vol)));
phaseImg.vol = (phaseImg.vol) / (max_ph-min_ph+1);

% Form the complex image, and unwrap the phase
%Z        = exp(2*pi*1i * (phaseImg.vol));
%fieldmap = atan2(imag(Z), real(Z));

for z=1:size(phaseImg.vol,3)
    fieldmap(:,:,z) = QualityGuidedUnwrap2D(squeeze(exp(2*pi*1i * phaseImg.vol(:,:,z))));
end
% Scale by the difference in echo times
TE1      = 4.92;
TE2      = 7.38;
dTE      = (TE2 - TE1) / 1000;      % seconds
fieldmap = fieldmap / (2*pi*dTE);   % Hz

%figure
%imagesc(squeeze(fieldmap(:,:,1)), [-pi , pi]  );
% ************************************************************************
% Set an appropriate threshold using the magnitude images
%**************************************************************************
%[nRows, nCols, nSlices] = size(fieldmap);


%figure, hIm = montage(reshape(magImg, [nRows, nCols, 1 nSlices]), ...
            %  'DisplayRange', [min(magImg(:)), max(magImg(:))]);        
%interactiveImageThresholder(hIm);          

%% ************************************************************************
% Threshold the result using the corresponding magnitude images
%**************************************************************************

%tHold = 296.746;
MagVol=squeeze(magImg.vol(:,:,:,1));
tHold=0.05*max(MagVol(:));

fieldmap(mean(MagVol,3)< tHold) = 0;


%% ************************************************************************
% Empirical reorientation following .rda orientation (and MRSIData out CombineDatFiles)
%**************************************************************************

fieldmap=flip(fieldmap,1);
MagVol=flip(MagVol,1);
% *************************************************************************
 %B0 needs to take opposite sign. Empirical observation
 %*************************************************************************
fieldmap=-fieldmap;

%figure
%imagesc(squeeze(fieldmap(:,:,1))  );
