function [combined_raw_trr,kmask,header,mrProt]=CombineDatFiles(DirPath,UndersamplingF)


[SCRIPT_DIR, trash1, trash2]=fileparts(mfilename('fullpath'));
CURRENT_DIR=pwd;
search_str = '/*.dat';

search_str = fullfile(DirPath, search_str);
filelist = dir(search_str);
if min(size(filelist))==0
    error(sprintf('ERROR in CombinedDatFile: could not find any *.dat files in %s \n',DirPath));
end

fprintf('reading file: %s \n', fullfile(DirPath,filelist(1).name));
comb_twix=reading_twix_small( fullfile(DirPath,filelist(1).name));

NDim=sum([comb_twix.hdr.Config.PhaseMatrix,comb_twix.hdr.Config.SliceMatrix,comb_twix.hdr.Config.ReadMatrix]>1);

if NDim==2
    comb_twix.raw_tckk(end,end,comb_twix.hdr.Config.PhaseMatrix,comb_twix.hdr.Config.ReadMatrix)=0;
elseif NDim==3
    comb_twix.raw_tckk(end,end,comb_twix.hdr.Config.PhaseMatrix,comb_twix.hdr.Config.SliceMatrix,comb_twix.hdr.Config.ReadMatrix)=0;
end


header=comb_twix.hdr;

for k=2:numel(filelist)
    fprintf('combine file: %s \n', fullfile(DirPath,filelist(k).name));
    twixk=reading_twix_small( fullfile(DirPath,filelist(k).name));
   % size(twixk.raw_tckk)
  if NDim==2
       twixk.raw_tckk(:,:,header.Config.PhaseMatrix,header.Config.ReadMatrix)=0;
    elseif NDim==3
        twixk.raw_tckk(:,:,header.Config.PhaseMatrix,header.Config.SliceMatrix,header.Config.ReadMatrix)=0;
   end

   if~all(size(comb_twix.raw_tckk)==size(twixk.raw_tckk))
      error('Raw data files to combine have different size (CombineDatFile.m, line 28)'); 
   end
    comb_twix.raw_tckk=comb_twix.raw_tckk+twixk.raw_tckk;
end
 clear  twixk
 
 comb_twix.raw_fckk=fft(conj( comb_twix.raw_tckk),[ ],1);
kmask=squeeze((sum(sum(abs(comb_twix.raw_tckk),1),2)>0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make undersampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(UndersamplingF<1)
    fprintf('Undersamplig data with factor: %s \n',UndersamplingF);
    if NDim==2
        [US_MASK]=Make_undersampled_mask(-1.0,UndersamplingF,kmask,8);
        kmask=kmask.*US_MASK;
       
        for i=1:size( comb_twix.raw_fckk,1);
            for j=1:size( comb_twix.raw_fckk,2); 
        	comb_twix.raw_fckk(i,j,:,:)=squeeze(comb_twix.raw_fckk(i,j,:,:)).*US_MASK;
            comb_twix.raw_tckk(i,j,:,:)=squeeze(comb_twix.raw_tckk(i,j,:,:)).*US_MASK;
        end;end
    elseif NDim>2
        error('Undersampling not supported for Ndim>2');
    end
end
%imagesc(kmask)
%figure()
%imagesc(squeeze(abs(sum(sum(comb_twix.raw_fckk,1),2))))
    
if NDim==2
    %comb_twix.raw_tcrr=fftshift(fft(fftshift(fft( comb_twix.raw_tckk,[ ],3), 3),[ ],4), 4);
    comb_twix.raw_tckk=fftshift(fftshift(comb_twix.raw_tckk,3),4);
    comb_twix.raw_tcrr=ifft(ifft(comb_twix.raw_tckk,[ ],3),[ ],4);
   comb_twix.raw_tcrr=fftshift(fftshift( comb_twix.raw_tcrr, 3), 4);
   
    comb_twix.raw_tcrr=flip(permute( comb_twix.raw_tcrr, [1,2,4,3]),3);
     kmask=permute( kmask,[2,1]);
      kmask=fftshift(fftshift( kmask,1),2);
    comb_twix.raw_tcrr=size( comb_twix.raw_tcrr,1)^2*conj( comb_twix.raw_tcrr);
    raw_fcrr=fft( comb_twix.raw_tcrr,[ ],1);

    for i=1:size( raw_fcrr,3);for j=1:size( raw_fcrr,4); 
        for coil_el=1:size( raw_fcrr,2) 
            rel_phase(coil_el,i,j)=angle(sum( raw_fcrr(:,coil_el,i,j),1));
            rel_amp(coil_el,i,j)=sum( abs(raw_fcrr(:,coil_el,i,j)).^2,1);
        end
        rel_amp(:,i,j)=rel_amp(:,i,j).^(0.5)/sum(rel_amp(:,i,j).^(0.5),1); 
    end;end;
   % clear raw_fcrr;
    
     comb_twix.rel_phase=rel_phase;
     comb_twix.rel_amp=rel_amp;

     %%%WSVD
     %Range_noise=[size(comb_twix.raw_fcrr,1)/2-250 size(comb_twix.raw_fcrr,1)/2+250];
     %for i=1:size( comb_twix.raw_fcrr,3);for j=1:size( comb_twix.raw_fcrr,4); 
      %  [MRSIdata_f(:,i,j) CCoef Q]=WSVD(squeeze(comb_twix.raw_fcrr(:,:,i,j)), Range_noise);
     %end;end

    [Trash Main_coil_element]=max(sum(sum(rel_amp(:,:,:),2),3));
    comb_twix.Main_coil_element=Main_coil_element;

    
    %showing the coil element volumes
    %{ 
for coil_el=1:size(comb_twix.raw_tcrr,2) 
       %phasemap(:,:,coil_el)= squeeze(sum(abs(comb_twix.raw_fcrr(215:240,coil_el,:,:)),1))/max(max(max(squeeze(sum(abs(comb_twix.raw_fcrr(215:240,:,:,:)),1)))));
    phasemap(:,:,coil_el)= (squeeze(angle(sum(raw_fcrr(:,coil_el,:,:),1)))+pi)/(2*pi);
       subplot(4,4,coil_el),subimage(phasemap(:,:,coil_el));
        colormap default
     end 
    %}

  %  combined_raw_frr=zeros(size( comb_twix.raw_fcrr,1),size( comb_twix.raw_fcrr,3),size( comb_twix.raw_fcrr,4));  
    combined_raw_trr=zeros(size( comb_twix.raw_tcrr,1),size( comb_twix.raw_tcrr,3),size( comb_twix.raw_tcrr,4)); 
    for i=1:size( comb_twix.raw_tcrr,3);for j=1:size( comb_twix.raw_tcrr,4)
        for coil_el=1:size( comb_twix.raw_tcrr,2) 
          %  combined_raw_frr(:,i,j)=combined_raw_frr(:,i,j)+ rel_amp(coil_el,i,j)*comb_twix.raw_fcrr(:,coil_el,i,j).*exp(-1i*(rel_phase(coil_el,i,j)-rel_phase(Main_coil_element,i,j)));
            combined_raw_trr(:,i,j)=combined_raw_trr(:,i,j)+ rel_amp(coil_el,i,j)*comb_twix.raw_tcrr(:,coil_el,i,j).*exp(-1i*(rel_phase(coil_el,i,j)-rel_phase(Main_coil_element,i,j)));
        end
    end;end;

    %at this point, the data are ordered and oriented following the .rda
    %conventions.

    %combined_raw_tkk=fftshift(fft(fftshift(fft( combined_raw_trr,[ ],2), 2),[ ],3), 3);
    %clear combined_raw_trr;

    
elseif NDim==3    
    
    %comb_twix.raw_tcrr=fftshift(fft(fftshift(fft(fftshift(fft( comb_twix.raw_tckk,[ ],3), 3),[ ],4), 4),[ ],5), 5);
    comb_twix.raw_tckk=fftshift(fftshift(fftshift(comb_twix.raw_tckk,3),4),5);
    comb_twix.raw_tcrr=ifft(ifft(ifft(comb_twix.raw_tckk,[ ],3),[ ],4),[ ],5);
   comb_twix.raw_tcrr=fftshift(fftshift(fftshift( comb_twix.raw_tcrr, 3), 4), 5);
   
     comb_twix.raw_tcrr= fftshift(fftshift(fftshift(comb_twix.raw_tcrr,3),4),5);
   comb_twix.raw_tcrr=flip(permute( comb_twix.raw_tcrr, [1,2,5,3,4]),3);
    kmask=permute( kmask,[3,1,2]);
    kmask=fftshift(fftshift(fftshift( kmask,1),2),3);
    comb_twix.raw_tcrr=size( comb_twix.raw_tcrr,1)^2*conj( comb_twix.raw_tcrr);
    raw_fcrr=fft( comb_twix.raw_tcrr,[ ],1);

    for i=1:size( raw_fcrr,3);for j=1:size( raw_fcrr,4);for k=1:size( raw_fcrr,5)    
        for coil_el=1:size( raw_fcrr,2) 
            rel_phase(coil_el,i,j,k)=angle(sum( raw_fcrr(:,coil_el,i,j,k),1));
            rel_amp(coil_el,i,j,k)=sum( abs(raw_fcrr(:,coil_el,i,j,k)).^2,1);
        end
        rel_amp(:,i,j,k)=rel_amp(:,i,j,k).^(0.5)/sum(rel_amp(:,i,j,k).^(0.5),1);
    end;end;end
    clear raw_fcrr;

    comb_twix.rel_phase=rel_phase;
    comb_twix.rel_amp=rel_amp;
    [Trash Main_coil_element]=max(sum(sum(sum(rel_amp(:,:,:,:),2),3),4));
    comb_twix.Main_coil_element=Main_coil_element;
    %{
    figure
    z_dim=4;
    for coil_el=1:size(comb_twix.raw_tcrr,2) 
    phasemap(:,:,coil_el)= squeeze(sum(abs(comb_twix.raw_fcrr(:,coil_el,:,z_dim,:)),1))/max(max(max(squeeze(sum(abs(comb_twix.raw_fcrr(:,:,:,z_dim,:)),1)))));
    subplot(4,4,coil_el),subimage(phasemap(:,:,coil_el));
    colormap default
    end  
    %}
    
   % combined_raw_frr=zeros(size( comb_twix.raw_fcrr,1),size( comb_twix.raw_fcrr,3),size( comb_twix.raw_fcrr,4),size( comb_twix.raw_fcrr,5));  
    combined_raw_trr=zeros(size( comb_twix.raw_tcrr,1),size( comb_twix.raw_tcrr,3),size( comb_twix.raw_tcrr,4),size( comb_twix.raw_tcrr,5)); 
    for i=1:size( comb_twix.raw_tcrr,3);for j=1:size( comb_twix.raw_tcrr,4);for k=1:size( comb_twix.raw_tcrr,5)
    for coil_el=1:size( comb_twix.raw_tcrr,2) 
        %combined_raw_frr(:,i,j,k)=combined_raw_frr(:,i,j,k)+ rel_amp(coil_el,i,j,k)*comb_twix.raw_fcrr(:,coil_el,i,j,k).*exp(-1i*(rel_phase(coil_el,i,j,k)-rel_phase(Main_coil_element,i,j,k)));
        combined_raw_trr(:,i,j,k)=combined_raw_trr(:,i,j,k)+ rel_amp(coil_el,i,j,k)*comb_twix.raw_tcrr(:,coil_el,i,j,k).*exp(-1i*(rel_phase(coil_el,i,j,k)-rel_phase(Main_coil_element,i,j,k)));
    end
    end;end;end

    %at this point, the data are ordered and oriented following the .rda
    %conventions.

    %combined_raw_tkk=fftshift(fft(fftshift(fft(fftshift(fft( combined_raw_trr,[ ],2), 2),[ ],3), 3),[ ],4), 4);
   % clear combined_raw_trr;
else
    error('the dimensionality of .dat files in CombineDatFile.m is different than 2 or 3. (line 127)');
end

mrProt.FoVHeight=header.Config.PeFOV;
mrProt.FoVWidth=header.Config.RoFOV;
mrProt.FoV3D=header.Config.VoI_SliceThickness;
mrProt.DwellTime=header.Config.DwellTime*1e-9;% in s
mrProt.NMRFreq=header.Config.Frequency*1E-6; %42.576*3; %in Mhz
mrProt.VSize=size(combined_raw_trr,1);

% save('CombinedDat-CSI.m','combined_raw_tkk','kmask','header','-mat');
% clear all;
end
