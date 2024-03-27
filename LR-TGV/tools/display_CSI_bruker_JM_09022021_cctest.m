%% Open and display Bruker CSI 
%JM - 09/02/2021 

clear;close all;

%% infos 
Nav=1;
MatSize=[31,31]; 
FidPoints=1024; %complex points
grpdly=77;  % these first points of FID are not ok -need to be taken out
%% open FID
expnb=22; % Bruker data are stored in folders with numbers based on order they were acquired. for teh data I sent you I modified the names of these folders 
fullpath="Z:\DATA\CSI\14T-FIDCSI-Data\20201207_100404_rat_07122020_rat_07122020_1_1\" +num2str(expnb)+"\ser"; % if fid then write fid if ser then write ser here

fileid=fopen(fullpath,'r','ieee-le'); %read binary format
if fileid == -1
    disp('Cannot open file');
    return
end

%% Read data		
buffer=fread(fileid,'double'); %note: CSI Bruker format is double 
buffer_c=buffer(1:2:end)+1i*buffer(2:2:end);

fid_mat_c=reshape(buffer_c, FidPoints,MatSize(1)*MatSize(2),Nav);

Nav_todisplay=1; %could also do mean here if Nav>1
fid_mat_c=fid_mat_c(:,:,Nav_todisplay)';

% group delay - circular shift
fid_mat_c_shift=[fid_mat_c(:,grpdly:end),fid_mat_c(:,1:grpdly-1)];

% saves if needed 
save('fid_mat_c.mat')
save('fid_mat_c_shift.mat')

%% Plot 
pix_to_display=602; %out of MatSize(1)*MatSize(2)

figure() %fid
plot(real(fid_mat_c_shift(pix_to_display,:)))
figure() %spectrum
ft_mat_c_shift=fftshift(fft(fid_mat_c_shift(pix_to_display,:),[],2),2);
plot(real(ft_mat_c_shift))
