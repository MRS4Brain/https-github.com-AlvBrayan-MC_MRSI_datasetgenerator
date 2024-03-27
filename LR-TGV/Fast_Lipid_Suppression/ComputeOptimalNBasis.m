function [Nbasis,SVratCoil] = ComputeOptimalNBasis( mrsiData_Lipids,mrsiData,mrsiReconParams );

%Local copy of the variable to allow parallelization
LipidsParams.mrProt.samplerate = mrsiReconParams.mrProt.samplerate;
LipidsParams.SkMask2D = mrsiReconParams.SkMask2D;
LipidsParams.BrainMask2D = mrsiReconParams.BrainMask2D;
LipidsParams.kmask = mrsiReconParams.kmask;
LipidsParams.L2SVDparams =  mrsiReconParams.L2SVDparams;
LipidsParams.Log_Dir = mrsiReconParams.Log_Dir;
LipidsParams.WaterFreqMap = mrsiReconParams.WaterFreqMap;
Size_data=size(mrsiData_Lipids);
Nfit = mrsiReconParams.L2SVDparams.Nfit;
%mrsiDataLr=zeros(size(mrsiData));


HzpP=mrsiReconParams.mrProt.samplerate/size(mrsiData,2);
low_bnd_L=round(100/HzpP);
high_bnd_L=round(500/HzpP);

Max_pts=round((high_bnd_L-low_bnd_L)*0.10);

BrainMask_frr=permute(repmat(LipidsParams.BrainMask2D ,[ 1, 1,(high_bnd_L-low_bnd_L+1)]),[3 1 2]);
NBasisMax=LipidsParams.L2SVDparams.NBasisMax;
NBasisMin=LipidsParams.L2SVDparams.NBasisMin;
NBasis=NBasisMin:LipidsParams.L2SVDparams.NBasisStep:NBasisMax;

for coil = 1 :  Size_data(1)
    fprintf(['Coil:',num2str(coil),','] );
    
    mrsiData_Lipids_frr=squeeze(fft(ifft(ifft(mrsiData_Lipids(coil,:,:,:),[],4),[],3),[],2));
   % [SortedLipid_frr , PeakPositions]=sort(abs(mrsiData_Lipids_frr(low_bnd_L:high_bnd_L,:,:)).^2,'descend');
     [SortedLipid_frr , PeakPositions]=sort(abs(mrsiData_Lipids_frr(low_bnd_L:high_bnd_L,:,:)),'descend');
    
     SortedLipid_frr=SortedLipid_frr.*BrainMask_frr;
    SumLipidOrig(coil)=sum(sum(sum(SortedLipid_frr(1:Max_pts,:,:),1),2),3);
    SumPeakSpreadOrig(coil)=sum(sum(std(PeakPositions(1:Max_pts,:,:),[],1),2),3);
    SumPeaksAreaOrig(coil)=sum(sum( sum(SortedLipid_frr(1:Max_pts,:,:),1)./(eps+std(PeakPositions(1:Max_pts,:,:),[],1)) ,2),3);
    
    clear mrsiData_Lipids_frr SortedLipid_frr;
    
    mrsiData_Lipids_tkk=squeeze(mrsiData_Lipids(coil,:,:,:));
    mrsiData_tkk=squeeze(mrsiData(coil,:,:,:));

    mrsiDataBrain_frr=fft(ifft(ifft(mrsiData_tkk,[],2),[],3),[],1);
    mrsiDataBrain_frr=mrsiDataBrain_frr(low_bnd_L:high_bnd_L,:,:).*BrainMask_frr;
    
    parfor Nb = 1:numel(NBasis)
    %    for Nb = 1:NBasisMax
        [~,LipidFree_frr,Lipid_rrf,Slip] = ProjSVDLipidSuppression( mrsiData_tkk,mrsiData_Lipids_tkk, LipidsParams,NBasis(Nb), Nfit,''); 
        SVratio(Nb)=Slip(end,end)/Slip(1,1);
       % Lipid_frr=sort(abs(permute(Lipid_rrf(:,:,low_bnd_L:high_bnd_L),[3 1 2])),'descend');
       % [SortedLipidFree_frr , PeakPositions]=sort(abs(LipidFree_frr(low_bnd_L:high_bnd_L,:,:)).^2,'descend'); 
    	[SortedLipidFree_frr , PeakPositions]=sort(abs(LipidFree_frr(low_bnd_L:high_bnd_L,:,:)),'descend'); 
        SortedLipidFree_frr=SortedLipidFree_frr.*BrainMask_frr;
    
        PeakPositions=PeakPositions.*BrainMask_frr;
        
        if(coil==mrsiReconParams.mrProt.Main_coil_element | coil==1 )
            Ptx_plot_data(Nb,:,:,:) = PeakPositions(1:Max_pts,:,:);
            Pty_plot_data(Nb,:,:,:) = SortedLipidFree_frr(1:Max_pts,:,:);
            Spect_plot_data(Nb,:,:,:) = LipidFree_frr(low_bnd_L:high_bnd_L,:,:).*BrainMask_frr;
        end
        
      %  SumLipid(coil,Nbasis)=sum(sum(sum(Lipid_frr(1:Max_pts,:,:),1),2),3);%squeeze(sum(sum(quantile(abs(Lipid_rrf),0.95,3),1),2));
        SumLipidFree(coil,Nb)=sum(sum(sum(SortedLipidFree_frr(1:Max_pts,:,:),1),2),3);%squeeze(sum(sum(quantile(abs(LipidFree_frr),0.95,1),2),3));
        SumPeakSpread(coil,Nb)=sum(sum(std(PeakPositions(1:Max_pts,:,:),[],1),2),3);
        SumPeaksArea(coil,Nb)=sum(sum( sum(SortedLipidFree_frr(1:Max_pts,:,:),1)./(eps+std(PeakPositions(1:Max_pts,:,:),[],1)) ,2),3);
        
        Lipid_frr=permute(Lipid_rrf(:,:,low_bnd_L:high_bnd_L),[3,1,2]).*BrainMask_frr;
        LipidFree_frr= LipidFree_frr(low_bnd_L:high_bnd_L,:,:).*BrainMask_frr;
        Consist_Sq(coil,Nb)=sumsqr(abs( LipidFree_frr(:)-mrsiDataBrain_frr(:)));
        Lipid_Sq(coil,Nb)=sumsqr(abs(Lipid_frr(:)-mrsiDataBrain_frr(:)));
        
        %SumPeakPositions(coil,Nb)=sum(sum(sum(PeakPositions(1:Max_pts,:,:),1),2),3);
        %figure
       % plot(low_bnd_L:high_bnd_L,abs(LipidFree_frr(low_bnd_L:high_bnd_L,45,45)),low_bnd_L-1+PeakPositions(1:Max_pts,45,45),abs(SortedLipidFree_frr(1:Max_pts,45,45)),'*')
    end
    SumLipidFree(coil,:)=smooth(squeeze(SumLipidFree(coil,:)),5);
    %RelDiffSumLipid(coil,:)=gradient(squeeze(SumLipid(coil,:)))./SumLipidGrad(coil,:);
    RelDiffSumLipidFree(coil,:)=diff(squeeze(SumLipidFree(coil,:)))./SumLipidOrig(coil);
     
    %RelDiffPeakPositions(coil,:)=gradient(squeeze(SumPeakPositions(coil,:)))./SumPeakPositions(coil,:);
    RelDiffPeakSpread(coil,:)=diff(squeeze(SumPeakSpread(coil,:)))./SumPeakSpreadOrig(coil);
    RelDiffPeaksArea(coil,:)=diff(squeeze(SumPeaksArea(coil,:)))./SumPeaksAreaOrig(coil);
  %  RelDiffProduct(coil,:)=gradient(squeeze(SumLipidFree(coil,:)./SumPeakSpread(coil,:)))./(SumLipidFree(coil,:)./SumPeakSpread(coil,:));

    RelSumLipidFree(coil,:)=(squeeze(SumLipidFree(coil,:)))/SumLipidOrig(coil);
    RelSumPeakSpread(coil,:)=(squeeze(SumPeakSpread(coil,:)))./SumPeakSpreadOrig(coil);
    RelSumPeaksArea(coil,:)=(squeeze(SumPeaksArea(coil,:)))./SumPeaksAreaOrig(coil);
    %Nbasis(coil)=min(find(abs(RelDiffSumLipidFree(coil,2:end))<mrsiReconParams.L2SVDparams.PercentThres));
    
    pp_SumLipidFree=spline(NBasis,squeeze(SumLipidFree(coil,:)));
    Interp_N=linspace(NBasisMin,NBasisMax,1E3);
    HD_DDSumLipidFree(coil,:)=ppval(fnder(pp_SumLipidFree,2),Interp_N);
    pp_SumLipidFree=spline(NBasis,squeeze(log(SumLipidFree(coil,:))));
    Interp_N=linspace(NBasisMin,NBasisMax,1E3);
    HD_DDLogSumLipidFree(coil,:)=ppval(fnder(pp_SumLipidFree,2),Interp_N);
    
    [M(coil),I(coil)]=max(gradient(gradient(squeeze(SumLipidFree(coil,:)))));
    
    PtUnderThres=find(abs(RelDiffSumLipidFree(coil,:))<mrsiReconParams.L2SVDparams.PercentThres);
   % PtUnderThres=find(gradient(gradient(squeeze(SumLipidFree(coil,:))))<M(coil)*mrsiReconParams.L2SVDparams.PercentThres);
    PtUnderThres(PtUnderThres<(NBasisMin+2))=999; % Avoid false minimum Values
    
    if ~isempty(PtUnderThres)
       Nbasis(coil)=NBasis(min(PtUnderThres));
       SVratCoil(coil)=SVratio(min(PtUnderThres));
    else
        fprintf(['\nWarning , unable to find a Nbasis value to go under the Lipid Suppression Threshold for coil:',num2str(coil),'. Value set to NBasisMax.\n'] );  
        Nbasis(coil)=NBasisMax;
        SVratCoil(coil)=SVratio(NBasisMax);
      %error('Unable to find a Nbasis value to go under the Lipid Suppression Threshold. ComputeOptimalNBasis.m') ;
    end
    
    %{
    [K(coil,:),B(coil,:),pp_Lip{coil},pp_Consist{coil}] = ComputeLCurvature( 1:numel(NBasis),squeeze(Consist_Sq(coil,:)) ,squeeze(Lipid_Sq(coil,:)), 1E5);
    [~,I]=max(abs(squeeze(K(coil,:))));
    Optimized_NBasis_L2_Supp(coil)=B(coil,I);
   %}
end

ind_Reg=1;
for coil = 1 : Size_data(1)
    for Nb = 1:numel(NBasis)
           if(Nb<7)
            X_Reg(ind_Reg)=Nb;
            Y_Reg(ind_Reg)=log(RelSumLipidFree(coil,Nb));
            ind_Reg=ind_Reg+1;
        end
    end
end


%plotting the optimization in a .ps
s=['./',mrsiReconParams.Log_Dir,'/','LipidNbasis_Optimization_', mrsiReconParams.NameData,'.ps'];
delete(s);
figs=figure('visible', 'off');

semilogy(NBasis,abs(RelSumPeaksArea)','r',NBasis,abs(RelSumLipidFree)','b',NBasis,abs(RelSumPeakSpread)','g',NBasis,mrsiReconParams.L2SVDparams.PercentThres,'-')
legend('RelSumPeaksArea','RelSumLipidFree','RelSumPeakSpread')
print(figs, '-append', '-dpsc2', s);

semilogy(NBasis,mrsiReconParams.L2SVDparams.PercentThres,'-',NBasis,abs(RelSumPeaksArea)','r');
title('RelSumPeaksArea');
print(figs, '-append', '-dpsc2', s);

semilogy(NBasis,mrsiReconParams.L2SVDparams.PercentThres,'-',NBasis,abs(RelSumLipidFree)','b');
title('RelSumLipidFree');
print(figs, '-append', '-dpsc2', s);

semilogy(NBasis,mrsiReconParams.L2SVDparams.PercentThres,'-',NBasis,abs(RelSumPeakSpread)','g');
title('RelSumProduct');
print(figs, '-append', '-dpsc2', s);

semilogy(NBasis(2:end),mrsiReconParams.L2SVDparams.PercentThres,'-',NBasis(2:end),abs(RelDiffPeaksArea)','r');
title('RelDiffPeaksArea');
print(figs, '-append', '-dpsc2', s);

semilogy(NBasis(2:end),mrsiReconParams.L2SVDparams.PercentThres,'-',NBasis(2:end),abs(RelDiffSumLipidFree)','b');
title('RelDiffSumLipidFree');
print(figs, '-append', '-dpsc2', s);

semilogy(NBasis(2:end),mrsiReconParams.L2SVDparams.PercentThres,'-',NBasis(2:end),abs(RelDiffPeakSpread)','g');
title('RelDiffPeakSpread');
print(figs, '-append', '-dpsc2', s);

%Regression:
X=[ones(length(X_Reg),1) X_Reg'];
b1 = X\Y_Reg';
yCalc = X*b1;
scatter(X_Reg',Y_Reg')
hold on
plot(X_Reg',yCalc')
hold off
xlabel('N');xlabel('log RelSumLip');legend([ num2str(b1(2)) '*x + ' num2str(b1(1)) ] )
print(figs, '-append', '-dpsc2', s);


semilogy(NBasis,squeeze(SumLipidFree),'-*b', NBasis,5*gradient(gradient(squeeze(SumLipidFree))),'-*r')
hold on
for coil = 1 :  Size_data(1)
    semilogy(I(coil),5*M(coil),'ok');
    semilogy(NBasis,(5*M(coil)*mrsiReconParams.L2SVDparams.PercentThres),'k');
end
hold off
title('Double derivative data');
print(figs, '-append', '-dpsc2', s);

clf
plot(NBasis,squeeze(SumLipidFree),'-*b', NBasis,5*gradient(gradient(squeeze(SumLipidFree))),'-*r')
hold on
for coil = 1 :  Size_data(1)
    plot(I(coil),5*M(coil),'ok');
    plot(NBasis,(5*M(coil)*mrsiReconParams.L2SVDparams.PercentThres),'k');
end
hold off
title('Double derivative data');
print(figs, '-append', '-dpsc2', s);


%{
plot(Consist_Sq',Lipid_Sq','-*');
hold on
for c=1:size(B,1)
    plot( ppval(pp_Consist{c},squeeze(B(c,:))),ppval(pp_Lip{c},squeeze(B(c,:))),'-r');
end
hold off
title('L-Curve');xlabel('Log Consist. Sq.');ylabel('Log Lipid Sq.');
print(figs, '-append', '-dpsc2', s);

semilogx(squeeze(B'),abs(squeeze(K')));
title('Curvature');xlabel('NBasis');ylabel('Curvature');
print(figs, '-append', '-dpsc2', s);
%}
for x=1:size(Spect_plot_data,3);
    if mrsiReconParams.ImMask2D(x,x)==1
        plot(low_bnd_L:high_bnd_L,abs(squeeze(Spect_plot_data(:,:,x,x))),low_bnd_L-1+squeeze(Ptx_plot_data(:,:,x,x)),abs(squeeze(Pty_plot_data(:,:,x,x))),'*r')
        title(['Peak Signal Suppression, Coil ',num2str(mrsiReconParams.mrProt.Main_coil_element),', Diagonal element:',num2str(x)]);
        print(figs, '-append', '-dpsc2', s);
    end
end
close all;


end

