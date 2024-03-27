function [SENSE SENSE_Adj]=TWIX_GRE_2_SENS(twix_filename,mrsiReconParams);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read TWIX data file and put them in the same format as reading_RDA after coil-element combination and truncation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(exist(twix_filename)~=2 ) error(strcat('ERROR: Twix .dat file ',twix_filename,' was not found! in TWIX_GRE_2_SENS.m')); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run mapVBVD but muted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_dump=evalc('raw_data=mapVBVD(twix_filename)');
twix.hdr=raw_data{2}.hdr;

raw_kcke=squeeze(mean(mean(raw_data{2}.image(),5),6));%%% average the slices and averages and squeeze
raw_ckk = permute(squeeze(raw_kcke(:,:,:,1)),[2,1,3]);%pick the first echo only
twix.RawSize=size(raw_ckk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interpolate the to its actual  Size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(twix.hdr.Config.RawCol*2==raw_data{2}.image.NCol)
    raw_crr = ifftshift(ifftshift(ifft(ifft(raw_ckk,[],2),[],3),2),3);
	temp=raw_crr(:,(end*1/4+1):(end*3/4),:); 
	raw_ckk=fft(fft(temp,[],2),[],3);
    raw_crr=temp;
    twix.RawSize(2)=twix.RawSize(2)/2;
else
	error(strcat('ERROR: GRE Twix .dat Column number is unsual in TWIX_GRE_2_SENS.m'));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
senseIter = 1000;   
SENSE=raw_crr;
for c = 1:twix.RawSize(1);
   % [SENSE(c,:,:) e] = h1_l2_2D_pd(squeeze(temp(c,:,:)),squeeze(raw_ckk(c,:,:)), 1000, senseIter,0.001);
  % smooth sensitivity estimate
    mask = [0 1 0; 1 4 1; 0 1 0]/8;
    for j=1:100
        SENSE(c,:,:) = conv2(squeeze(abs(SENSE(c,:,:))),mask,'same');
    end
end

 SENSE_Adj=pinv(reshape(SENSE,twix.RawSize(1),twix.RawSize(2)*twix.RawSize(3)));
SENSE_Adj=reshape(SENSE_Adj,twix.RawSize(2),twix.RawSize(3),twix.RawSize(1));
end
