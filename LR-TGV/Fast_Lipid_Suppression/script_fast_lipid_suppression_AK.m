%% load invivo high res spectroscopy data (20 avers) and also 2 aver data 
%% 2 aver data : used for estimating high res lipid image in dual-density recon
%% 20 aver data : used for generating low resolution (0.56 cc) 20 aver image
addpath(genpath('/data/MR_data/CombineDatFiles'))
addpath(genpath('/data/MR_data/MRSI_Recon_V1'))
load PrePreProcReconResults;          % 2 average .16cc data %r-r-f
     % 20 average .16cc data %r-r-f
 mrsiReconParams.mrsiData=PG_extrap(mrsiReconParams.mrsiData, mrsiReconParams,0.10,'Combined' );
 
 Data_rrt=ifft(ifft(permute(mrsiReconParams.mrsiData,[2,3,1]),[],1),[],2);
 N = size(Data_rrt);
 HzpP=mrsiReconParams.mrProt.samplerate/N(3);
 b0map2D=imresize(mean(mrsiReconParams.b0map,3),[N(1),N(2)]);

for l = 1 : size(Data_rrt,3)
  % Data_rrt(:,:,l)=exp(-2*pi*1i*(l-1) * -b0map2D / mrsiReconParams.mrProt.samplerate).*squeeze(Data_rrt(:,:,l));
end

img_2avg=fft(Data_rrt,[],3);
img=img_2avg;



%load lipid_mask;      % lipid layer mask
%load meta_mask;       % brain mask

lipid_mask=imresize(sum(mrsiReconParams.SkMask,3),[N(1),N(2)]);
lipid_mask=round(lipid_mask/max(lipid_mask(:)));

low_bnd_L=round(200/HzpP);
high_bnd_L=round(500/HzpP);
 Lipid_Vol=squeeze(sum(abs(img(:,:,low_bnd_L:high_bnd_L)),3));
 
Lipid_Vol=Lipid_Vol/max(Lipid_Vol(:));
LThr=quantile(Lipid_Vol(:),0.5);
lipid_mask=((Lipid_Vol.*lipid_mask)>LThr);
  



meta_mask=imresize(sum(mrsiReconParams.BrainMask,3),[N(1),N(2)]);
%
meta_mask=round(meta_mask/max(meta_mask(:)));

plot_projection([img, img_2avg]), title('Sum over frequency, Left: 20 aver -- Right: 2 aver MRSI data')


%% generate dual density image by combining 2 aver high-resolution and 20 aver low-resolution images

%csi_rad = 16;                       % size of low res image, has lipid ringing
%mask_csi = repmat(circular_mask([64,64], csi_rad), [1,1,N(3)]); %kspace mask
                      % size of low res image, has lipid ringing
                      
mask_csi = repmat(mrsiReconParams.kmask, [1,1,N(3)]); %kspace mask

mask_high = ones(N) - mask_csi; 
FT_high = FT_v2(mask_high);         % outer k-space
FT_csi  = FT_v2(mask_csi);          % inner k-space
FT_full = FT_v2(ones(N));           % full k-space

tic
    data1 = FT_high * (img_2avg .* repmat(lipid_mask,[1,1,N(3)]));      % high-res k-space from 2avg data
    %data1 = FT_csi * (img_2avg .* repmat(lipid_mask,[1,1,N(3)]));
    data2 = FT_csi * img;                                               % low-res k-space from 20avg data
toc

%combined_2avg = FT_full' * (data1+data2);           % dual-density image   
combined_2avg = img;           % dual-density image                        
plot_projection( log(abs([combined_2avg, img])),2), title('Sum over frequency in dB, Left: Dual-density recon -- Right: 20 aver MRSI data')


%% apply iterative L1 lipid-basis recon to dual-density image

param = init;
param.xfmWeight = 1e-3;       % Lipid basis penalty chosen by L-curve heuristic

param.FT = 1;
param.data = param.FT * combined_2avg; 
param.Itnlim = 10;

param.Lipid = get_LipidBasis(combined_2avg, lipid_mask);    % Lipid basis functions
param.Bmask = meta_mask;      % Brain mask 
img_combined_L1 = combined_2avg;

tic
%{
for t = 1:10
    img_combined_L1 = fnlCg_LipidBasis(img_combined_L1, param);
    figure(1), imagesc(  [ sum(abs(img_combined_L1),3).*meta_mask , sum(abs(combined_2avg),3).*meta_mask ] ), axis image, colorbar, drawnow
end
%}
toc

img_combined_L1_meta = img_combined_L1 .* repmat(meta_mask,[1,1,N(3)]);

% compare spectra from dual-density and L1-lipid-basis recon
overplot_spectra(img_combined_L1_meta, combined_2avg, 30, 40, 5, 5)


%% apply L2 lipid-basis recon to dual-density image

beta = 100;%6.5e-1;      % beta chosen to match the data consistency of iterative L1-lipid-basis recon

tic

Lipid = make_LipidBasis(combined_2avg, lipid_mask); 
Lipid_inv = inv( eye(N(3)) + beta * (Lipid * Lipid') );


img_combined_L2 = combined_2avg;

for ay = 1:size(img_combined_L2,1)
    for cey = 1:size(img_combined_L2,2)
        
        if meta_mask(ay,cey)
            
            xi = img_combined_L2(ay,cey,:);
            xi_L2 = Lipid_inv * xi(:);
            img_combined_L2(ay,cey,:) = xi_L2;
            
        end
        
    end
end
toc

img_combined_L2_meta = img_combined_L2 .* repmat(meta_mask,[1,1,N(3)]);
combined_2avg_meta = combined_2avg .* repmat(meta_mask,[1,1,N(3)]);

plot_projection([img_combined_L2_meta, img_combined_L1_meta, combined_2avg_meta], 1, [0,.3])
title('Sum over frequency inside brain mask, Left: L2-lipid-basis  --  Middle: L1-lipid-basis  --  Right: Dual-density recon')

disp(['L2-lipid-basis data consistency: ', num2str(norm( img_combined_L2(:) - combined_2avg(:) )^2)])        % L2 consistency
disp(['L1-lipid-basis data consistency: ', num2str(norm( img_combined_L1(:) - combined_2avg(:) )^2)])               % L1 consistency



%% compare dual-density, L1 and L2 lipid-basis recon spectra
   
overplot_spectra( combined_2avg_meta, img_combined_L2_meta,  46, 38, 4, 4, 1, N(3)-45)
overplot_spectra( combined_2avg_meta, img_combined_L2_meta,  38, 18, 6, 6, 1, N(3)-45, 3, 4)

overplot_spectra( img_combined_L1_meta, img_combined_L2_meta,  46, 38, 4, 4, 1, N(3)-45, 5, 6)
overplot_spectra( img_combined_L1_meta, img_combined_L2_meta,  38, 18, 6, 6, 1, N(3)-45, 7, 8)



%% generate lipid maps

lipid_range = 55:195;

lipid_img = 20*log10(sum(abs([img(:,:,lipid_range), combined_2avg(:,:,lipid_range), img_combined_L1(:,:,lipid_range), img_combined_L2(:,:,lipid_range)]),3));

h = figure(1); close(h);
figure(1), imagesc(lipid_img, [-35, 15]), axis image off, colormap jet
title('Lipid maps in dB, Left: no suppression, Mid-Left: Dual-density, Mid-Right: L1-lipid-basis, Right: L2-lipid-basis')


% plot section through the lipid image 
x_section = 40;

lipid1 = 20*log10(sum(abs(img(x_section,:,lipid_range)),3));
lipid2 = 20*log10(sum(abs(combined_2avg(x_section,:,lipid_range)),3));
lipid3 = 20*log10(sum(abs(img_combined_L1(x_section,:,lipid_range)),3));
lipid4 = 20*log10(sum(abs(img_combined_L2(x_section,:,lipid_range)),3));

h = figure(2); close(h);
figure(2), hold on, plot(lipid1(:)), axis tight
figure(2), hold on, plot(lipid2(:), 'k'), axis tight
figure(2), hold on, plot(lipid3(:), 'r'), axis tight
figure(2), hold on, plot(lipid4(:), 'm'), axis tight
title('Section through the lipid image')


