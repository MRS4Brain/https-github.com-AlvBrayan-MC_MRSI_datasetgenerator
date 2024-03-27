%% ESPIRiT Maps Demo
% This is a demo on how to generate ESPIRiT maps. It is based on the paper
% Uecker et. al, MRM 2013 DOI 10.1002/mrm.24751. ESPIRiT is a method that
% finds the subspace of multi-coil data from a calibration region in
% k-space using a series of eigen-value decompositions in k-space and image
% space. 

%%
% Set parameters

%load PreProcReconResults_HBRef_WS16C.mat

DATA=permute(mrsiReconParams.Water_ctkk(:,5,:,:),[3,4,1,2]);
DATA=fft(fft(fftshift(fftshift(ifft(ifft(DATA,[],1),[],2),1),2),[],1),[],2);
DATA=fftshift(fftshift(DATA,1),2);

[sx,sy,Nc] = size(DATA);

%ncalib = 24; % use 24 calibration lines to compute compression
ncalib = round(min(sx,sy)/sqrt(2.0)); % FOR ELLIPTICAL ENCODING


%size = [6,6]; % kernel size
ksize = [12,12]; % kernel size

% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.

eigThresh_1 = 0.02;

% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);

%%
% Display coil images: 
% im = ifft2c(DATA);
% figure, imshow3(abs(im),[],[4,ceil(Nc/4)]); 
% title('magnitude of physical coil images');
% colormap((gray(256))); colorbar;
% 
% figure, imshow3(angle(im),[],[4,ceil(Nc/4)]); 
% title('phase of physical coil images');
% colormap('default'); colorbar;

%% Compute ESPIRiT EigenVectors
% Here we perform calibration in k-space followed by an eigen-decomposition
% in image space to produce the EigenMaps. 


% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels

[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));

%% 
% This shows that the calibration matrix has a null space as shown in the
% paper. 
% 
% kdisp = reshape(k,[ksize(1)*ksize(2)*Nc,ksize(1)*ksize(2)*Nc]);
% figure, subplot(211), plot([1:numel(S)],S,'LineWidth',2);
% hold on, 
% plot([1:ksize(1)*ksize(2)*Nc],S(1)*eigThresh_1,'r-','LineWidth',2);
% plot([idx,idx],[0,S(1)],'g--','LineWidth',2)
% legend('signular vector value','threshold')
% title('Singular Vectors')
% subplot(212), imagesc(abs(kdisp)), colormap(gray(256));
% xlabel('Singular value #');
% title('Singular vectors')


%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

%%
% show eigen-values and eigen-vectors. The last set of eigen-vectors
% corresponding to eigen-values 1 look like sensitivity maps

% figure, imshow3(abs(W),[],[4,ceil(Nc/4)]); 
% title('Eigen Values in Image space');
% colormap((gray(256))); colorbar;
% 
% figure, imshow3(abs(M),[],[Nc,Nc]); 
% title('Magnitude of Eigen Vectors');
% colormap(gray(256)); colorbar;
% 
% figure, imshow3(angle(M),[],[Nc,Nc]); 
% title('Magnitude of Eigen Vectors');
% colormap(jet(256)); colorbar;
% 

%%
% project onto the eigenvectors. This shows that all the signal energy
% lives in the subspace spanned by the eigenvectors with eigenvalue 1.
% These look like sensitivity maps. 


% alternate way to compute projection is:
% ESP = ESPIRiT(M);
% P = ESP'*im;

% P = sum(repmat(im,[1,1,1,Nc]).*conj(M),3);
% figure, imshow3(abs(P),[],[4,ceil(Nc/4)]); 
% title('Magnitude of Eigen Vectors');
% colormap(sqrt(gray(256))); colorbar;
% 
% figure, imshow3(angle(P),[],[4,ceil(Nc/4)]); 
% title('Magnitude of Eigen Vectors');
% colormap((jet(256))); colorbar;
% 


%%
% crop sensitivity maps 
maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);

figure, imshow3(abs(maps),[],[4,ceil(Nc/4)]); 
title('Absolute sensitivity maps');
colormap((gray(256))); colorbar;

figure, imshow3(angle (maps),[],[4,ceil(Nc/4)]); 
title('Phase of sensitivity maps');
colormap((jet(256))); colorbar;


