function PlotDiagSpectraData(reconResults,Data1_frr,Data2_frr,name)


% ************************************************************************
% Plot Diagonal spectra 
%*************************************************************************



TimeSize=size(Data1_frr,1);
Final_ppm=reconResults.ppm; %(-4.7+((1:TimeSize)*reconResults.mrProt.samplerate/(TimeSize*reconResults.mrProt.NMRFreq)));

pts=1:size(Data1_frr,1);%reconResults.MinPPM_pt:reconResults.MaxPPM_pt;

xstep=0.1*pts(end);
ystep=max(abs(squeeze(Data1_frr(pts,size(Data1_frr,2)/2,size(Data1_frr,2)/2))));% 0.01 Braino %0.1 Invivo
figs=figure(); 
hold on;
for x=1:size(Data1_frr,2)
    y=size(Data1_frr,2)-x+1;
    plot(pts(end)+pts+x*xstep,abs(squeeze(Data1_frr(pts,x,y)))+(size(Data1_frr,2)-x)*ystep,'g-' )
    plot(pts+x*xstep,abs(squeeze(Data2_frr(pts,x,y)))+(size(Data1_frr,2)-x)*ystep, 'r-')

end 

axis([pts(1) pts(end)+2*(pts(end)-pts(1))+size(Data1_frr,2)*xstep 0 size(Data1_frr,2)*ystep])

 print(figs, '-dpsc2', name);
 close all;
 
end