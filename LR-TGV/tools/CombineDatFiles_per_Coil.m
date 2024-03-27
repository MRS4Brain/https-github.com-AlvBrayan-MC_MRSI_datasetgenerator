function [raw_ctkk,kmask,header,MrProt]=CombineDatFiles_per_Coil(FilePath,UndersamplingF,DoPermute)


[SCRIPT_DIR, trash1, trash2]=fileparts(mfilename('fullpath'));
CURRENT_DIR=pwd;
VoxelShift=0.5;

%{
search_str = '/*.dat';
search_str = fullfile(DirPath, search_str);
filelist = dir(search_str);
if min(size(filelist))==0
    error(sprintf('ERROR in CombinedDatFile: could not find any *.dat files in %s \n',DirPath));
end
fprintf('reading file: %s \n', fullfile(DirPath,filelist(1).name));
comb_twix=reading_twix_small(fullfile(DirPath,filelist(1).name));
%}
fprintf('reading file: %s \n', FilePath);
comb_twix=reading_twix_small( FilePath);
header=comb_twix.hdr;

if isfield(header.MeasYaps.sProtConsistencyInfo,'tMeasuredBaselineString')
    SoftVersionString=header.MeasYaps.sProtConsistencyInfo.tMeasuredBaselineString;
elseif isfield(header.MeasYaps.sProtConsistencyInfo,'tBaselineString')
   SoftVersionString=header.MeasYaps.sProtConsistencyInfo.tBaselineString;
elseif isfield(header.Meas,'SoftwareVersions');
    SoftVersionString=header.Meas.SoftwareVersions;
else
    error('Could not determine the software version from the .dat header');
end


if(~isempty(strfind(SoftVersionString,'N4_VD13'))| ~isempty(strfind(SoftVersionString,'N4_VE11')))
    MrProt.Nlines=header.Config.PhaseEncodingLines;
    if isfield(header.Config,'ReadResolution')
        MrProt.Ncol=header.Config.ReadResolution;
    else
    	MrProt.Ncol=header.Config.BaseResolution;  
    end
    MrProt.Nslc=header.Config.NPar;
    
    MrProt.FoVHeight=header.Config.PeFOV;
    MrProt.FoVWidth=header.Config.RoFOV;

    if isfield(header.Config,'VoI_SliceThickness')
        MrProt.FoV3D=header.Config.VoI_SliceThickness;
    elseif isfield(header.Meas,'VoiThickness')
        MrProt.FoV3D=header.Meas.VoiThickness;
    else  
         MrProt.FoV3D=header.Spice.VoiThickness;
    end

    if isfield(header.Config,'DwellTime')
        MrProt.DwellTime=header.Config.DwellTime*1e-9;% in s
    elseif isfield(header.Meas,'alDwellTime')
        MrProt.DwellTime= header.Meas.alDwellTime(1)*1e-9;;
    else  
        MrProt.DwellTime= header.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;
    end

    if isfield(header.Config,'ReadResolution')
        MrProt.Sequence=header.Config.SequenceString;
    else
        MrProt.Sequence=header.Meas.SequenceString;
    end
    if isfield(header.Config,'Frequency')
        MrProt.NMRFreq = header.Config.Frequency*1E-6; %42.576*3; %in Mhz
    elseif isfield(header.Meas,'lFrequency')
        MrProt.NMRFreq = header.Meas.lFrequency*1E-6;
    else  
        MrProt.NMRFreq = header.Dicom.lFrequency*1E-6;
    end
    
    
    NDim=sum([ MrProt.Nlines,MrProt.Ncol,MrProt.Nslc]>1);
    if NDim==2
        comb_twix.raw_tckk(end,end,MrProt.Nlines,MrProt.Ncol)=0;
    elseif NDim==3
        comb_twix.raw_tckk(end,end,MrProt.Nlines,MrProt.Nslc,MrProt.Ncol)=0;
    end
    
elseif(~isempty(strfind(SoftVersionString,'B17')))
    MrProt.Nlines=header.Config.PhaseEncodingLines;
    MrProt.Ncol= header.Config.BaseResolution;
    MrProt.Nslc=header.Meas.Partitions;
    
    %check if data have the right size and otherwise force the right size
    comb_twix.raw_tckk=permute(comb_twix.raw_tckk,[1,2,4,3]);%2D VB17 Raw DATA are apparently switched
    comb_twix.raw_tckk=ResizeDataToRightSize(comb_twix.raw_tckk,MrProt.Ncol,MrProt.Nlines );

    MrProt.FoVHeight=header.Config.PeFOV;
    MrProt.FoVWidth=header.Config.RoFOV;
    MrProt.FoV3D=header.Meas.VoiThickness;
    if isfield(header.Config,'ReadResolution')
        MrProt.Sequence=header.Config.SequenceString;
    else
        MrProt.Sequence=header.Meas.SequenceString;
    end
    MrProt.DwellTime=header.Meas.alDwellTime(1)*1e-9;% in s
    MrProt.NMRFreq=header.Meas.lFrequency*1E-6; %42.576*3; %in Mhz
    
else
    error(strcat('Twix file header format: ',SoftVersionString,' is unknown!'));
end
MrProt.header=header;
    
NDim=sum([ MrProt.Nlines,MrProt.Ncol,MrProt.Nslc]>1);

%Now resized above for VB17 already 
% if NDim==2 
%     comb_twix.raw_tckk(end,end,MrProt.Nlines,MrProt.Ncol)=0;
% elseif NDim==3
%     comb_twix.raw_tckk(end,end,MrProt.Nlines,MrProt.Nslc,MrProt.Ncol)=0;
% end

%{
if(mod(size(comb_twix.raw_tckk,3),2)==1 & size(comb_twix.raw_tckk,3)>2 ) 
        fprintf('Making %s dimensions even.\n Original: %g ,Now: %g \n', FilePath,size(comb_twix.raw_tckk,3),size(comb_twix.raw_tckk,3)+1);
	comb_twix.raw_tckk(:,:,end+1,:,:)=0;			
end
if(mod(size(comb_twix.raw_tckk,4),2)==1 & size(comb_twix.raw_tckk,4)>2 ) 
	fprintf('Making %s dimensions even.\n Original: %g ,Now: %g \n', FilePath,size(comb_twix.raw_tckk,4),size(comb_twix.raw_tckk,4)+1);
	comb_twix.raw_tckk(:,:,:,end+1,:)=0;			
end
if(mod(size(comb_twix.raw_tckk,5),2)==1 & size(comb_twix.raw_tckk,5)>2 ) 
	fprintf('Making %s dimensions even.\n Original: %g ,Now: %g \n', FilePath,size(comb_twix.raw_tckk,5),size(comb_twix.raw_tckk,5)+1);
	comb_twix.raw_tckk(:,:,:,:,end+1)=0;		
end
%}



if NDim==2

	[Kr,Kc] = meshgrid(0:(size(comb_twix.raw_tckk,4)-1),0:(size(comb_twix.raw_tckk,3)-1));
	Kr=2*pi*Kr/size(comb_twix.raw_tckk,4)*double(~mod(size(comb_twix.raw_tckk,4),2));
	Kc=2*pi*Kc/size(comb_twix.raw_tckk,3)*double(~mod(size(comb_twix.raw_tckk,3),2));
	comb_twix.raw_tckk=comb_twix.raw_tckk.* reshape(exp(1i*VoxelShift*Kr).*exp(1i*VoxelShift*Kc),[1 1 size(Kr)]);


	comb_twix.raw_tckk=fftshift(fftshift(comb_twix.raw_tckk,3),4);
	comb_twix.raw_tcrr=ifft(ifft(comb_twix.raw_tckk,[ ],3),[ ],4);
	comb_twix.raw_tcrr=fftshift(fftshift( comb_twix.raw_tcrr, 3), 4);

   
   % comb_twix.raw_tcrr=circshift(comb_twix.raw_tcrr,-1,3);
   % comb_twix.raw_tcrr=circshift(comb_twix.raw_tcrr,-1,4);  


   if(strcmp(DoPermute,'permute'))
    comb_twix.raw_tcrr=flip(permute( comb_twix.raw_tcrr, [1,2,4,3]),3);
   end

    comb_twix.raw_tcrr=conj( comb_twix.raw_tcrr);%*size( comb_twix.raw_tcrr,1)^2;
    raw_fcrr=fft( comb_twix.raw_tcrr,[ ],1);

    for i=1:size( raw_fcrr,3);for j=1:size( raw_fcrr,4); 
        for coil_el=1:size( raw_fcrr,2) 
            rel_phase(coil_el,i,j)=angle(sum( raw_fcrr(:,coil_el,i,j),1));
            rel_amp(coil_el,i,j)=sum( abs(raw_fcrr(:,coil_el,i,j)).^2,1);
        end
        rel_amp(:,i,j)=rel_amp(:,i,j).^(0.5)/sum(rel_amp(:,i,j).^(0.5),1); 
    end;end;
    clear raw_fcrr;
    
     comb_twix.rel_phase=rel_phase;
     comb_twix.rel_amp=rel_amp;

     %%%WSVD
     %Range_noise=[size(comb_twix.raw_fcrr,1)/2-250 size(comb_twix.raw_fcrr,1)/2+250];
     %for i=1:size( comb_twix.raw_fcrr,3);for j=1:size( comb_twix.raw_fcrr,4); 
      %  [MRSIdata_f(:,i,j) CCoef Q]=WSVD(squeeze(comb_twix.raw_fcrr(:,:,i,j)), Range_noise);
     %end;end

    [Trash Main_coil_element]=max(sum(sum(rel_amp(:,:,:),2),3));
    comb_twix.Main_coil_element=Main_coil_element;
    
    comb_twix=rmfield(comb_twix, {'raw_tckk'});
     raw_ctkk=squeeze(fft(fft(permute(comb_twix.raw_tcrr, [2,1,3,4]),[],3),[],4));
elseif NDim==3    
    
    %comb_twix.raw_tcrr=fftshift(fft(fftshift(fft(fftshift(fft( comb_twix.raw_tckk,[ ],3), 3),[ ],4), 4),[ ],5), 5);
    comb_twix.raw_tckk=fftshift(fftshift(fftshift(comb_twix.raw_tckk,3),4),5);
    comb_twix.raw_tcrr=ifft(ifft(ifft(comb_twix.raw_tckk,[ ],3),[ ],4),[ ],5);
   comb_twix.raw_tcrr=fftshift(fftshift(fftshift( comb_twix.raw_tcrr, 3), 4), 5);
   
  % kmask=fftshift(fftshift(fftshift( kmask,1),2),3);
   
   
   if(strcmp(DoPermute,'permute'))
    comb_twix.raw_tcrr=flip(permute( comb_twix.raw_tcrr, [1,2,5,3,4]),3);
  %  kmask=permute( kmask,[3,1,2]);%flippinf in real space doesn't flip k-space (of course but cumun mistake)

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1Voxel shift measured empirically   
    comb_twix.raw_tcrr=circshift(circshift(circshift(comb_twix.raw_tcrr,1,3),-1,4),-1,5);
   else
    comb_twix.raw_tcrr=circshift(circshift(circshift(comb_twix.raw_tcrr,-1,3),-1,4),-1,5); 
   end
    
    
    comb_twix.raw_tcrr=conj( comb_twix.raw_tcrr);%*size( comb_twix.raw_tcrr,1)^2;
   
    
     [Trash Main_coil_element]=max(squeeze(sum(sum(sum(sum(abs(comb_twix.raw_tcrr(:,:,:,:)),1),3),4),5)));
    comb_twix.Main_coil_element=Main_coil_element;
    
    comb_twix=rmfield(comb_twix, {'raw_tckk'});
    raw_ctkk=fft(fft(fft(permute(comb_twix.raw_tcrr, [2,1,3,4,5]),[],3),[],4),[],5);
else
    error('the dimensionality of .dat files in CombineDatFile.m is different than 2 or 3. (line 127)');
end




MrProt.VSize=size(raw_ctkk,2);
MrProt.OriginalVSize=MrProt.VSize;
%MrProt.rel_phase=comb_twix.rel_phase;
%MrProt.rel_amp=comb_twix.rel_amp;
MrProt.Main_coil_element=comb_twix.Main_coil_element;
kmask = log(squeeze(sum(sum(abs(raw_ctkk).^2,1),2)))>-20;
% save('CombinedDat-CSI.m','combined_raw_tkk','kmask','header','-mat');
% clear all;
end
