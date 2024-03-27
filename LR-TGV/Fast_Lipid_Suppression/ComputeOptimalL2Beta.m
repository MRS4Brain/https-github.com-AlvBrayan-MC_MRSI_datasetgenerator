function Betas = ComputeOptimalL2Beta( mrsiData_Lipids_ctkk,mrsiReconParams );


%Local copy of the variable to allow for parallelization


NbPtOpt=mrsiReconParams.L2Sup.NbPtOpt;
LipidsParams.mrProt.samplerate = mrsiReconParams.mrProt.samplerate;
LipidsParams.SkMask2D = mrsiReconParams.SkMask2D;
LipidsParams.BrainMask2D = mrsiReconParams.BrainMask2D;
LipidsParams.ImMask2D = mrsiReconParams.ImMask2D;
LipidsParams.kmask = mrsiReconParams.kmask;
LipidsParams.L2SVDparams =  mrsiReconParams.L2SVDparams;
LipidsParams.Log_Dir = mrsiReconParams.Log_Dir;
LipidsParams.WaterFreqMap = mrsiReconParams.WaterFreqMap;
Size_data=size(mrsiData_Lipids_ctkk);

%mrsiData=mrsiReconParams.mrsiData;
%mrsiDataLr=zeros(size(mrsiData));

HzpP=mrsiReconParams.mrProt.samplerate/Size_data(2);
low_bnd_L=round(100/HzpP);
high_bnd_L=round(500/HzpP);

Max_pts=round((high_bnd_L-low_bnd_L)*0.10);

BrainMask_frr=permute(repmat(LipidsParams.BrainMask2D ,[ 1, 1,(high_bnd_L-low_bnd_L+1)]),[3 1 2]);

Betas(1:Size_data(1))=1000;
BetaValues=1E2*10.^(0.5*(0:1:(NbPtOpt-1)));

Ptx_plot_data=zeros(NbPtOpt,Max_pts,size(BrainMask_frr,2),size(BrainMask_frr,3));
Pty_plot_data=zeros(NbPtOpt,Max_pts,size(BrainMask_frr,2),size(BrainMask_frr,3));
Spect_plot_data=zeros(NbPtOpt,size(BrainMask_frr,1),size(BrainMask_frr,2),size(BrainMask_frr,3));

for coil = 1 :   Size_data(1)
    fprintf(['Coil:',num2str(coil),','] );
 
    mrsiData_Lipids_frr=squeeze(fft(ifft(ifft(mrsiData_Lipids_ctkk(coil,:,:,:),[],4),[],3),[],2));
    %[SortedLipid_frr , PeakPositions]=sort(abs(mrsiData_Lipids_frr(low_bnd_L:high_bnd_L,:,:)).^2,'descend');
    [SortedLipid_frr , PeakPositions]=sort(abs(mrsiData_Lipids_frr(low_bnd_L:high_bnd_L,:,:)),'descend');
    clear mrsiData_Lipids_frr;
     SortedLipid_frr=SortedLipid_frr.*BrainMask_frr;
    SumLipidOrig(coil)=sum(sum(sum(SortedLipid_frr(1:Max_pts,:,:),1),2),3);
    SumPeakSpreadOrig(coil)=sum(sum(std(PeakPositions(1:Max_pts,:,:),[],1),2),3);
    SumPeaksAreaOrig(coil)=sum(sum( sum(SortedLipid_frr(1:Max_pts,:,:),1)./(eps+std(PeakPositions(1:Max_pts,:,:),[],1)) ,2),3);
   
    mrsiData_Lipids_tkk=squeeze(mrsiData_Lipids_ctkk(coil,:,:,:));   
    mrsiDataBrain_frr=fft(ifft(ifft(mrsiData_Lipids_tkk,[],2),[],3),[],1);
    mrsiDataBrain_frr=mrsiDataBrain_frr(low_bnd_L:high_bnd_L,:,:).*BrainMask_frr;
    
    
    
    parfor Nb = 1:NbPtOpt
        %ifor Nb = 1:NbPtOpt
        [~,mrsiData_LipidSup_rrf, Lipid_rrf]  = L2LipidSuppression( mrsiData_Lipids_tkk, LipidsParams, BetaValues(Nb),'');
        
        % Lipid_frr=sort(abs(permute(Lipid_rrf(:,:,low_bnd_L:high_bnd_L),[3 1 2])),'descend');
        LipidFree_frr=permute(mrsiData_LipidSup_rrf,[3,1,2]);
        %[SortedLipidFree_frr , PeakPositions]=sort(abs(LipidFree_frr(low_bnd_L:high_bnd_L,:,:)).^2,'descend');
        [SortedLipidFree_frr , PeakPositions]=sort(abs(LipidFree_frr(low_bnd_L:high_bnd_L,:,:)),'descend');
        SortedLipidFree_frr=SortedLipidFree_frr.*BrainMask_frr;
        
        PeakPositions=PeakPositions.*BrainMask_frr;
        if(coil==1);%mrsiReconParams.mrProt.Main_coil_element)
            Ptx_plot_data(Nb,:,:,:) = PeakPositions(1:Max_pts,:,:);
            Pty_plot_data(Nb,:,:,:) = SortedLipidFree_frr(1:Max_pts,:,:);
            Spect_plot_data(Nb,:,:,:)= LipidFree_frr(low_bnd_L:high_bnd_L,:,:).*BrainMask_frr;
        end
        SumLipidFree(coil,Nb)=sum(sum(sum(SortedLipidFree_frr(1:Max_pts,:,:),1),2),3);%squeeze(sum(sum(quantile(abs(LipidFree_frr),0.95,1),2),3));
        SumPeakSpread(coil,Nb)=sum(sum(std(PeakPositions(1:Max_pts,:,:),[],1),2),3);
        SumPeaksArea(coil,Nb)=sum(sum( sum(SortedLipidFree_frr(1:Max_pts,:,:),1)./(eps+std(PeakPositions(1:Max_pts,:,:),[],1)) ,2),3);
        
        Lipid_frr=permute(Lipid_rrf(:,:,low_bnd_L:high_bnd_L),[3,1,2]).*BrainMask_frr;
        LipidFree_frr= LipidFree_frr(low_bnd_L:high_bnd_L,:,:).*BrainMask_frr;
        Consist_Sq(coil,Nb)=sumsqr(abs( LipidFree_frr(:)-mrsiDataBrain_frr(:)));
        Lipid_Sq(coil,Nb)=sumsqr(abs(Lipid_frr(:)));
       
    end
 %   SumLipidOrig(coil)=SumLipidFree(coil,1);
  %  SumPeakSpreadOrig(coil)=SumPeakSpread(coil,1);
  %  SumPeaksAreaOrig(coil)=SumPeaksArea(coil,1);
    
    RelDiffSumLipidFree(coil,:)=diff(squeeze(SumLipidFree(coil,:)))./SumLipidOrig(coil);
    RelDiffPeakSpread(coil,:)=diff(squeeze(SumPeakSpread(coil,:)))./diff(BetaValues(:))';
    RelDiffPeaksArea(coil,:)=diff(squeeze(SumPeaksArea(coil,:)))./diff(BetaValues(:))';
    
    RelSumLipidFree(coil,:)=(squeeze(SumLipidFree(coil,:)))/SumLipidOrig(coil);
    RelSumPeakSpread(coil,:)=(squeeze(SumPeakSpread(coil,:)))./SumPeakSpreadOrig(coil);
    RelSumPeaksArea(coil,:)=(squeeze(SumPeaksArea(coil,:)))./SumPeaksAreaOrig(coil);
    
    pp_SumLipidFree=spline(BetaValues(1:end),squeeze(SumLipidFree(coil,:)));
    Interp_Beta=exp(linspace(log(BetaValues(1)),log(BetaValues(end)),1E3));
    HD_DSumLipidFree(coil,:)=ppval(fnder(pp_SumLipidFree,1),Interp_Beta);
    HD_DDSumLipidFree(coil,:)=ppval(fnder(pp_SumLipidFree,2),Interp_Beta);
    
    [M(coil),I(coil)]=max(gradient(gradient(squeeze(SumLipidFree(coil,:)))));
    
    [Max_RelDiff,IMax_RelDiff]=max(squeeze(abs(RelDiffSumLipidFree(coil,:))));
    PtUnderThres=find(abs(RelDiffSumLipidFree(coil,:))<mrsiReconParams.L2SVDparams.PercentThres);
    PtUnderThres(PtUnderThres<IMax_RelDiff)=Max_RelDiff;
    
   % PtUnderThres=find(gradient(gradient(squeeze(SumLipidFree(coil,:))))<M(coil)*mrsiReconParams.L2SVDparams.PercentThres);  
   
    
    % [~, PtatThres]=min(abs(abs(HD_RelDiffSumLipidFree(coil,:))-mrsiReconParams.L2SVDparams.PercentThres)));
   % PtatThres=find( abs(RelDiffSumLipidFree(coil,:))<mrsiReconParams.L2SVDparams.PercentThres );
    %[~,Imax]=max(abs(RelDiffSumLipidFree(coil,:)));
    if isempty(PtUnderThres)
        error('Unable to find a Beta Value value to go under the Lipid Suppression Threshold. ComputeOptimalL2Beta.m') ;
    end
    
    Betas(coil)=BetaValues(min(PtUnderThres))
    
   % [K(coil,:),B(coil,:),pp_Lip{coil},pp_Consist{coil}] = ComputeLCurvature(BetaValues,squeeze(Lipid_Sq(coil,:)),squeeze(Consist_Sq(coil,:)) , 1E5);
   % [~,I]=max(abs(squeeze(K(coil,:))));
   
end

ind_Reg=1;
for coil = 1 :  size(RelSumLipidFree,1)
    for Nb = 1:NbPtOpt
        if(BetaValues(Nb)>=1E3 & BetaValues(Nb)<1E6)
            X_Reg(ind_Reg)=Nb;%log(BetaValues(Nb));
            Y_Reg(ind_Reg)=log(RelSumLipidFree(coil,Nb));
            ind_Reg=ind_Reg+1;
        end
    end
end


%plotting the optimization in a .ps
s=['./',mrsiReconParams.Log_Dir,'/','L2LipSupp_Optimization_', mrsiReconParams.NameData,'.ps'];
delete(s);
figs=figure('visible', 'off');

loglog(BetaValues,abs(RelSumPeaksArea),'r',BetaValues,abs(RelSumLipidFree),'b',BetaValues,abs(RelSumPeakSpread),'g',Interp_Beta,mrsiReconParams.L2SVDparams.PercentThres,'-')
legend('RelSumPeaksArea','RelSumLipidFree','RelSumPeakSpread')
print(figs, '-append', '-dpsc2', s);

loglog(BetaValues,mrsiReconParams.L2SVDparams.PercentThres,'-',BetaValues,abs(RelSumPeaksArea),'-*r');
title('RelSumPeaksArea');
print(figs, '-append', '-dpsc2', s);

loglog(BetaValues,mrsiReconParams.L2SVDparams.PercentThres,'-',BetaValues,abs(RelSumLipidFree),'-*b');
title('RelSumLipidFree');
print(figs, '-append', '-dpsc2', s);

loglog(BetaValues,mrsiReconParams.L2SVDparams.PercentThres,'-',BetaValues,abs(RelSumPeakSpread),'-*g');
title('RelSumProduct');
print(figs, '-append', '-dpsc2', s);

loglog(BetaValues(2:end),mrsiReconParams.L2SVDparams.PercentThres,'-',BetaValues(2:end),abs(RelDiffPeaksArea),'r');
title('RelDiffPeaksArea');
print(figs, '-append', '-dpsc2', s);

loglog(BetaValues(2:end),mrsiReconParams.L2SVDparams.PercentThres,'-',BetaValues(2:end),abs(RelDiffSumLipidFree),'b');
title('RelDiffSumLipidFree');
print(figs, '-append', '-dpsc2', s);

loglog(BetaValues(2:end),mrsiReconParams.L2SVDparams.PercentThres,'-',BetaValues(2:end),abs(RelDiffPeakSpread),'g');
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
xlabel('Beta nb');ylabel('log RelSumLip');legend([ num2str(b1(2)) '*x + ' num2str(b1(1)) ] )
print(figs, '-append', '-dpsc2', s);

semilogy(log(BetaValues),squeeze(SumLipidFree),'-*b', log(BetaValues),5*gradient(gradient(squeeze(SumLipidFree))),'-*r')
hold on
for coil = 1 :  numel(I)
    semilogy(log(BetaValues(I(coil))),5*M(coil),'ok');
    semilogy(log(BetaValues),(5*M(coil)*mrsiReconParams.L2SVDparams.PercentThres),'k');
end
hold off
title('Double derivative data');
print(figs, '-append', '-dpsc2', s);

clf
plot(log(BetaValues),squeeze(SumLipidFree),'-*b', log(BetaValues),5*gradient(gradient(squeeze(SumLipidFree))),'-*r')
hold on
for coil = 1 : numel(I)
    plot(log(BetaValues(I(coil))),5*M(coil),'ok');
    plot(log(BetaValues),(5*M(coil)*mrsiReconParams.L2SVDparams.PercentThres),'k');
end
hold off
title('Double derivative data');
print(figs, '-append', '-dpsc2', s);

%{
plot(Consist_Sq',Lipid_Sq','-*')
hold on
for c=1:size(B,1)
    plot( ppval(pp_Consist{c},squeeze(B(c,:))),ppval(pp_Lip{c},squeeze(B(c,:))),'-r');
end
hold off
title('L-Curve');xlabel('Log Consist. Sq.');ylabel('Log Lipid Sq.');
print(figs, '-append', '-dpsc2', s);

semilogx(squeeze(B'),abs(squeeze(K')));
title('Curvature');xlabel('Beta');ylabel('Curvature');
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

