% Code Authors: Alves B
% Script that will go over the simulated sets (generated from
% Create_MRSI_datasets.m) and denoise them using MP-PCA (MP.m) and LR-TGV (MRSI-LRTGV.m).
% Output will be RAW files for all 3 methodologies, ready for LCModel
% quantification


% If you this code is of any use for your research project, please cite one of the following papers in your manuscript :

%% CODE THAT WILL SCROLL THROUGH A MC STUDY AND APPLY THE LR RECONSTRUCTION FROM ANTOINE

close all;
clear all;
clc;

start_path = ''; %Path where the simulated sets are located (TO FILL)
cd(start_path)
simu_name = ['']; %Names of the simulated sets (can be an array)
water_name = [""]; %Name of the water reference (used as as surrogate, only one set necessary)
load('border_rapid_coil_vmask.mat'); %Name of the file containing the map geometry (such as rodent brain, in a MATLAB file)

for loop=1:length(simu_name)
    
    cd(start_path)
    name = simu_name(loop);
    load([convertStringsToChars(name) '.mat']);

    water = load(convertStringsToChars(water_name(loop)));
    water_nrrt = water.study.data.real+sqrt(-1)*water.study.data.imag;
    
    power_mask = border.power_mask;
    
    data_nrrt = study.data.real+sqrt(-1)*study.data.imag;
    data_nkkt = fft(fft(data_nrrt,[],2),[],3);
    B0map_nrr = study.shift;
    studyname = study.filename;
    
    dims = size(data_nrrt);
    
    clear study;
    %% Computing useful data and checking the number of components of MPPCA
    
    
    time_start = datetime;
    
    Nx = dims(2);
    Ny = dims(3);
    Nel = dims(1);
    
    tot_npars = zeros(1,Nel);
    
    
    %% DENOISING LOOP
    for ii = 1:Nel
        cd(start_path)
        num = num2str(ii);
        metab_data_tkk = permute(squeeze(data_nkkt(ii,:,:,:)),[3,1,2]);
        water_signal = permute(squeeze(water_nrrt(ii,:,:,:)),[3,1,2]);
        B0map = squeeze(B0map_nrr(ii,:,:));
    
        %%%% MP-PCA Denoising %%%%
        MPPCA_signal = zeros(size(metab_data_tkk));
        metab_data_trr = ifft(ifft(metab_data_tkk,[],2),[],3);
        number_spectra = sum(sum(power_mask));
        temp2 = [];
        vox_pos = [];
        for x=1:Nx
            for y=1:Ny
                if(power_mask(x,y))
                    temp2 = [temp2;real(metab_data_trr(:,x,y))';imag(metab_data_trr(:,x,y))'];
                    vox_pos = [vox_pos;x,y];
                end
            end
        end
    
        [dn_full_fid,sigma,npars]=MP_savingFig(temp2,50,true,[convertStringsToChars(name) '_newSlice_N' num]);
        tot_npars(ii) = npars;
        dn_full_fid = dn_full_fid';
        count = 1;
        for p = 1:2:2*number_spectra-1
            real_n = dn_full_fid(:,p);
            imag_n = dn_full_fid(:,p + 1);
            pos = vox_pos(count,:);
            MPPCA_signal(:,pos(1),pos(2)) = real_n + 1i*imag_n;
            count = count+1;
        end
        
        %%%% LR Denoising %%%%
        LR_signal = MRSI_LRTGV(metab_data_tkk,B0map,power_mask,'AcqParam.m',[convertStringsToChars(name) '_' num],1E-3,0.5,0,npars,1);
    
    
        % CREATING RAW FILES TO BE PROCESSED BY LCMODEL 
        
    
        DenoisingToRaw_Maker(metab_data_trr,MPPCA_signal,LR_signal,water_signal,'C:\Data Processing\',studyname,['newSlice_N' num],border);
    
    end
    
    cd(start_path)
    save(['MPPCA_Components_' convertStringsToChars(name) '.mat'],"tot_npars");
    clear test test_mppca test_lr water_nrrt metab_data_trr MPPCA_signal LR_signal
    clear data_nkkt B0map B0map_nrr study water data_struct

    disp(['It is well and truly over ! You can take a break now :) '])
    time_end = datetime;
end