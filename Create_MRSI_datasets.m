% Code Authors: Mosso J & Alves B
% Synthetic Dataset generator for Monte-Carlo Simulations - from basis set files
% This code generates a Nx X Ny MRSI from basis-set elements simulated in
% NMR ScopeB.

% The code was inspired by the work from Mosso et al. (GitHub :
% https://github.com/jessie-mosso/DWMRS-MPPCA)

% If you this code is of any use for your research project, please cite one of the following papers in your manuscript :
% Mosso J, Simicic D, Şimşek K, Kreis R, Cudalbu C, Jelescu IO. MP-PCA denoising for diffusion MRS data: promises and pitfalls. NeuroImage. 2022;263:119634. doi:10.1016/j.neuroimage.2022.119634

clear; clc; close all;

%% to choose: 

b0drift=false; %add b0 drift distortions between repetitions
singleavSNR=1000; %single average SNR on real and imaginary
MM=false; %include MM in the spectra

saveMnoiseless=false; %
saveMtot=true; 

%be careful here: need to change the CreateStudyStruc.m depending on your application
%this "save" saves the output matrices for Nico's GUI

%Matrix dimension
Nx = 31;
Ny = 31;

ncomp = 3; % Number of compartiments 
na=30; % Number of MC realization
np=2048; % Number of time points in the FID
nMC=1; % generate nMC times the input matrix for Monte Carlo studies
folder_name = ['nMc=' convertStringsToChars(num2str(na)) '_SNR=' convertStringsToChars(num2str(singleavSNR))];


dw=1/7142.85714; % dwell time (1/Bandwidth) [s]
sfrq1H=599.0787; % 1H resonance frequency at 14.1T [Hz] (to change if necessary)
receiveroffset_ppm = 4.7; % Add an offset similar to in vivo acquisition

%Load the power map used as a reference for the simulation of the MRSI
%slice

t=[0:dw:(np-1)*dw]; 
fmax=1/(2*dw);
f=[fmax:-2*fmax/(np-1):-fmax];
scale_ppm=f/sfrq1H+receiveroffset_ppm;

load('border_rapid_coil_vmask.mat'); % Brain map (to change if necessary)
brain_mask = border.power_mask;

%% B_0 Shift computation from the Water signal

if b0drift
    water_file = load('MRSI_MatlabReconResult_rapidVmask_1Em3TGC_0.5LipSup_16WaterSup_20LR.mat'); %Change the name of the text file containing the frequency if necessary
    WaterFreqMap_rr=water_file.reconResults.WaterFreqMap;

    % WaterFreqMap_rr=load(''); % Load a .mat matrix containing the B0 shift 
end


%% Name & paths to find the simulated metabolites
%You need to make sure that all simulated metabolites that are in the
%variable 'metnames' can be found in .RAW format. Change 'metnames' and
%'pathfidraw' if necessary 

path_code = ''; %CHANGE TO MAKE IT FIT FOR YOUR PATH
pathfidraw=[path_code '\fid_simul\'];

if MM
    metnames=["Ala","Asc","Asp","Cr","PCr","GABA","Gln","Glu","GSH","Ins","Lac","NAA","Tau","Glc-A","Glc-B","NAAG","PE","GPC","PCho","MM_steam_clean_11trunc_try1"]; % bhb replaced by MM here for simplicity
else
    metnames=["Ala","Asc","Asp","Cr","PCr","GABA","Gln","Glu","GSH","Ins","Lac","NAA","Tau","Glc-A","Glc-B","NAAG","PE","GPC","PCho"]; % bhb replaced by MM here for simplicity
end
    
nbptsrawfiles=2048; %number of points in the .RAW files

metfids=zeros(nbptsrawfiles*2,length(metnames)); %x2 for real and imaginary (which are interleaved in the .RAW file)

%% Load metabolites FID (.RAW)
for nmet=1:length(metnames)
    fid = fopen([pathfidraw + metnames(nmet) + ".RAW"]);
    C = textscan(fid, '%f', 8*512, 'headerlines', 2);
    metfids (:,nmet)= C{1,1};
    fclose(fid);
end

%% simulation parameters

% 3 compartiment separation, following the geometry  :

lim_x = 19;
lim_y = 15;
lim_radius = 7;
center_x = 19;
center_y = 16;

pfactor = 0; % To use only if you want to add a tumor region, otherwise put this at 0


% Hippocampus with tumour in here 
conc(4,:)= [ 0.64, 2.49, 1.74,4.26, 4.99, 1.28 , 3.22 ,10.05 ,0.89 ,7.80*pfactor ,1.33 , 9.52 , 6.93 , 1.00  , 1.00  , 0.80 , 2.49, 0.50, 0.41 , 0.7];%5];%0.2 for bhb
% Hippocampus
conc(1,:)= [ 0.64, 2.49, 1.74,4.26, 4.99, 1.28 , 3.22 ,10.05 ,0.89 ,7.80 ,1.33 , 9.52 , 6.93 , 1.00  , 1.00  , 0.80 , 2.49, 0.50, 0.41 , 0.7];%5];%0.2 for bhb
% Striatum
conc(2,:)= [ 0.67, 2.04, 1.60,4.03, 4.46, 1.44 , 4.79 ,9.68  ,1.06 ,5.19 ,1.44 , 8.10 , 9.74 , 1.10  , 1.10  , 0.74 ,3.42, 0.53, 0.83 , 1.9];
% Cortex
conc(3,:)= [ 0.84, 2.69, 1.20,5.95, 6.67, 1.14 , 3.94 ,10.44  ,0.62 ,6.95 ,1.30 , 9.28 , 5.50 , 1.26  , 1.26  , 1.03 ,1.41, 0.45, 0.53 , 1.82];



Mnoiseless_t=zeros(ncomp,np);
linewidths=zeros(ncomp,1);

for MCit=1:nMC
    disp(["MCit=" + num2str(MCit)])
    
    Mtot_t=zeros(na,Nx,Ny,np);  %will contain complex FID
        
    % random fluctuations definition
    randomphases=rand(na,Nx,Ny)./2; %change amplitude if needed

    %% Construct the ncomps 1D spectra
    for comps=1:ncomp

        temp0=0;
        for c=1:size(metfids,2)
            data = metfids(:,c);

            if c==20 %flip the MM spectrum
                metbs_r=data(1:2:end);
                metbs_i=data(2:2:end);
                flipmetbs_r=metbs_r;
                flipmetbs_i=-metbs_i;%- on imaginary part of FID ==> flip symetrically on the spectral dimension
                metbs_new=zeros(np,1);
                metbs_new=flipmetbs_r + sqrt(-1)*flipmetbs_i;
                data_c=[metbs_new(1)/2;metbs_new(2)/2;metbs_new(3:end)];
                data_c = data_c'; % transpose of the data
            else

                data=[data(1)/2; data(2)/2; data(3:np*2)]; 
                a=1:2:(length(data)-1);
                b=2:2:length(data);
                
                data_c=data(a)+1i*data(b);
                
                data_c = data_c'; % transpose of the data
                data_c=conj(data_c); % conjugate of the data
            end

            %Take out the TMS using HSVD
    
            hue =  data_c';
            hue =  Fast_HSVD_Filter(hue,7143,16,-1700,-1600)';

            list_hue(c,:) = hue;
        
            
            %ponderated by its corresponding concentration
            temp0=temp0+hue*conc(comps,c);
        end
        
        %Correction term for the receiver offset (12.55 was the optimal
        %value for correcting)

        time_pts = [0:1:np-1]*dw;
        shift_spec = 12.55*sfrq1H;
        corr = exp(sqrt(-1).*shift_spec.*time_pts);
        spectrum_c_t(comps,:)=temp0.*corr; %complex - fid tot at 1 bval

        % add LB on noiseless spectrum 
        tt=0.000200:0.000200:np*0.000200;
        spectrum_c_t(comps,:)=spectrum_c_t(comps,:).*exp(-tt*15);

        %Calculate the linewidth
        spect = fftshift(fft(squeeze(spectrum_c_t(comps,:)),[],2),2);
        [val,~] = max(abs(spect));
        int_spect = abs(spect)-ones(size(spect,1),1)*val*0.5;
        list = find(int_spect(1475:1500)>=0);

        left_ind = list(1)+1474;
        right_ind = list(end)+1474;

        %Then, we translate those indexes into frequencies (Hz)
        left_Hz = f(left_ind);
        left_Hz_under = f(left_ind-1);
        right_Hz = f(right_ind);
        right_Hz_under = f(right_ind+1);

        %Finally, we use a linear interpolation to find the expected
        %frequency where we are at half maximum

        left_alpha = int_spect(left_ind-1)/(int_spect(left_ind)-int_spect(left_ind-1));
        left_ans = (1+left_alpha)*left_Hz_under - left_alpha*left_Hz;

        right_alpha = int_spect(right_ind)/(int_spect(right_ind+1)-int_spect(right_ind));
        right_ans = (1+right_alpha)*right_Hz - right_alpha*right_Hz_under;

        linewidths(comps) = left_ans-right_ans;


    end
   

    %% Store noiseless matrix
    if MCit==1 
        Mnoiseless_t=spectrum_c_t;
        if saveMnoiseless
            CreateStudyStruc_BA("noiseless",['M_noiseless'],[pwd],Mnoiseless_t,na,np,Nx,Ny);
        end 
    end


    %% Define noise level from the first shell, first iteration 
    if MCit==1 %only once
        sigmanoise=max(real(spectrum_c_t(1,:)))/singleavSNR; 
        disp(sigmanoise);
    end

    %% add water spectrum with random phase
    ph_water=-pi+randn(1)*2*pi; %randon phase for water residual
    fid_water=exp(-tt.*50).*exp(i*ph_water); %lorenzian shape with random phase


    %% Construction of the 2D MRSI slice

    for nx=1:Nx
        for ny=1:Ny
            %% construct the na iterations of noise
            for n=1:na  
                %% add noise
                spectrum_c_tnoisy=spectrum_c_t+sigmanoise*(randn(ncomp,np)+1i*randn(ncomp,np)); 
                
                %% add B0 drift - add small random shift
                if b0drift
                    pts=[0:1:np-1]*dw;    
                    shift_spec=WaterFreqMap_rr(nx,ny); %change amplitude if needed
                    shift_spec_tot(n,nx,ny)=shift_spec; %save the b0 drift distortions applied
                    spectrum_c_tnoisy=spectrum_c_tnoisy.*exp(i.*shift_spec.*pts);
                end

                % Circle
                if (sqrt((nx-center_x)^2+(ny-center_y)^2) > lim_radius)
                    Mav(n,nx,ny,:)=spectrum_c_tnoisy(3,:)*brain_mask(nx,ny); %Cortex Concentration

                elseif nx < lim_x
                    Mav(n,nx,ny,:)=spectrum_c_tnoisy(2,:)*brain_mask(nx,ny); %Striatum Concentration
                else
                    Mav(n,nx,ny,:)=spectrum_c_tnoisy(1,:)*brain_mask(nx,ny); % Hippocampus Concentration
                end

            end
        end
    end


    %% Noisy matrix na x nx x ny x np for 1 MC
    Mtot_t=Mav; 
    
    %% Store noisy matrix
    if saveMtot
        CreateStudyStruc_BA("noisy",['testPaper_SNR' num2str(singleavSNR) '_' num2str(MCit)],['C:\' folder_name],Mtot_t,na,np,Nx,Ny,shift_spec_tot,linewidths);
    end 

end