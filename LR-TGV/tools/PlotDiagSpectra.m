function PlotDiagSpectra(reconResults)


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
    Data_frr=fft(reconResults.Recon_ZP_Data_trr,[],1);
elseif isfield(reconResults,'mrsiDataLipRem_trr')
  Data_frr=fft(reconResults.mrsiDataLipRem_trr,[],1);
 else
    error('No resulting dataset found in reconResults!')
end   

%2nd Resulting DataSet    
if isfield(reconResults,'mrsiDataLipRem2_trr')
    Data2_frr=fft(reconResults.mrsiDataLipRem2_trr,[],1);
    Data2_frr=permute(Data2_frr,[1,3,2]);
end
Data_frr=permute(Data_frr,[1,3,2]);
DatawL_frr=permute(DatawL_frr,[1,3,2]);


TimeSize=size(DatawL_frr,1);
Final_ppm=reconResults.ppm; %(-4.7+((1:TimeSize)*reconResults.mrProt.samplerate/(TimeSize*reconResults.mrProt.NMRFreq)));
[~, min_ppm_pt]=min(abs(Final_ppm+4.2))
[~, max_ppm_pt]=min(abs(Final_ppm-0))
pts=min_ppm_pt:max_ppm_pt;
xstep=0.3;
ystep=0.25*max(abs(squeeze(Data_frr(pts,size(Data_frr,2)/2,size(Data_frr,2)/2))));% 0.01 Braino %0.1 Invivo
figs=figure(); 
hold on;
size(Data_frr)
for x=14:(size(Data_frr,2)-14)
    y=(size(Data_frr,3))-x+1;
	if (x>0 & y>0)
	    plot(Final_ppm(pts(end))-Final_ppm(pts(1))+Final_ppm(pts)+x*xstep,abs(squeeze(Data_frr(pts,x,y)))+(size(Data_frr,2)-x)*ystep,'g-' )
	    plot(Final_ppm(pts)+x*xstep,abs(squeeze(DatawL_frr(pts,x,y)))+(size(Data_frr,2)-x)*ystep, 'r-')
	    if isfield(reconResults,'mrsiDataLipRem2_trr')
		plot(2*Final_ppm(pts(end))-2*Final_ppm(pts(1))+Final_ppm(pts)+x*xstep,abs(squeeze(Data2_frr(pts,x,y)))+(size(Data_frr,2)-x)*ystep, 'b-')
	    end
	end
end 

axis([Final_ppm(pts(1)) Final_ppm(pts(end))+2*(Final_ppm(pts(end))-Final_ppm(pts(1)))+size(Data_frr,2)*xstep 0 size(Data_frr,2)*ystep])

 print(figs, '-dpsc2', ['./',reconResults.Log_Dir,'/','DiagonalSpectra_plot_', reconResults.NameData, '.ps']);
 close all;
 
end
