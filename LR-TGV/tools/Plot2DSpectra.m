function Plot2DSpectra(reconResults)


% ************************************************************************
% Plot Diagonal spectra 
%*************************************************************************

%Reference DataSet
if isfield(reconResults,'mrsiData_wLip_trr')
    DatawL_frr=fft(reconResults.mrsiData_wLip_trr,[],1);
elseif isfield(reconResults,'Original_RePhased_Data_trr')
    DatawL_frr=fft(reconResults.Original_RePhased_Data_trr,[],1);
else
    error('No reference dataset found in reconResults!')
end

%Resulting DataSet
if isfield(reconResults,'Recon_ZP_Data_trr')
    'Use Recon Data'
    Data_frr=fft(reconResults.Recon_ZP_Data_trr,[],1);
elseif isfield(reconResults,'mrsiDataLipRem_trr')
    'Use LipRem Data'
  Data_frr=fft(reconResults.mrsiDataLipRem_trr,[],1);
 else
    error('No resulting dataset found in reconResults!')
end   

TimeSize=size(DatawL_frr,1);
Final_ppm=reconResults.ppm; %(-4.7+((1:TimeSize)*reconResults.mrProt.samplerate/(TimeSize*reconResults.mrProt.NMRFreq)));
[~,maxppm]=min(abs(-1.5 - reconResults.ppm));
pts=reconResults.MinPPM_pt:maxppm;

marginx=10;
marginy=5;
xstep=1.1;
ystep=2*max(abs(squeeze(Data_frr(pts,size(Data_frr,2)/2,size(Data_frr,2)/2))));% 0.01 Braino %0.1 Invivo
jump=2;
s=['./',reconResults.Log_Dir,'/','2DSpectra_plot_', reconResults.NameData, '.ps'];
if exist(s);delete(s);end
figs=figure('units','inches');
hold on;
for x=(1+marginx):jump:(size(Data_frr,2)-marginx)
    for y=(1+marginy):jump:(size(Data_frr,3)-marginy)
   % plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),abs(squeeze(Data_frr(pts,x,y)))+x/jump*ystep,'b-' )
     plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),imag(squeeze(Data_frr(pts,x,y)))+x/jump*ystep,'r-' )
      plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),real(squeeze(Data_frr(pts,x,y)))+x/jump*ystep,'b-' )
    end
end 
%axis([Final_ppm(pts(1)) Final_ppm(pts(end))+2*(Final_ppm(pts(end))-Final_ppm(pts(1)))+(Final_ppm(pts(end))-Final_ppm(pts(1)))*size(Data_frr,2)/jump 0 size(Data_frr,2)/jump*ystep])
set(gca,'xtick',[])
set(gca,'ytick',[])

pos = get(gcf,'pos');
%set(gcf,'pos',[pos(1) pos(2) 20 20])

print(figs,'-append', '-dpsc2',s );
close all 

ystep=2*max(abs(squeeze(DatawL_frr(pts,size(DatawL_frr,2)/2,size(DatawL_frr,2)/2))));% 0.01 Braino %0.1 Invivo
figs=figure('units','inches');
hold on;
for x=(1+marginx):jump:(size(Data_frr,2)-marginx)
    for y=(1+marginy):jump:(size(Data_frr,3)-marginy)
   % plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),abs(squeeze(DatawL_frr(pts,x,y)))+x/jump*ystep,'g-' )
     plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),imag(squeeze(DatawL_frr(pts,x,y)))+x/jump*ystep,'r-' )
      plot(xstep*(Final_ppm(pts(end))-Final_ppm(pts(1)))*(y-1)/jump+Final_ppm(pts),real(squeeze(DatawL_frr(pts,x,y)))+x/jump*ystep,'b-' )
    end
end 

%axis([Final_ppm(pts(1)) Final_ppm(pts(end))+2*(Final_ppm(pts(end))-Final_ppm(pts(1)))+(Final_ppm(pts(end))-Final_ppm(pts(1)))*size(Data_frr,2)/jump 0 size(Data_frr,2)/jump*ystep])
set(gca,'xtick',[])
set(gca,'ytick',[])
pos = get(gcf,'pos');
%set(gcf,'pos',[pos(1) pos(2) 20 20])

print(figs,'-append', '-dpsc2',s );
 
 close all;
 
end
