function [U_rrc,V_tc,S,costFunVal,FreqMap] = tgv2_l2_2D_multiCoil_LowRank_CombinedConv(Data_ckkt,U_rrc,V_tc,S, alpha0, alpha1, maxits,minits,mrsiReconParams,Threshold)
% Primal dual TGV2 algorithm adapated for spatio-spectral reconstruction.

% As described in the TGV paper :
% Knoll, F.; Bredies, K.; Pock, T.; Stollberger, R.: Second Order Total
% Generalized Variation (TGV) for MRI: Magnetic Resonance in Medicine, 

% Original code given in: 
% https://cai2r.net/resources/second-order-tgv-reconstruction-for-undersampled-radial-mri/

Data_ckkt=(single(Data_ckkt));
U_rrc=(single(U_rrc));
V_tc=(single(V_tc));
S=(single(S));

[dxm,dym,dzm,dxp,dyp,dzp] = defineDiffOperators();
check_it =mrsiReconParams.LRTGVModelParams.check_it;%25;%25;
Plot_it=mrsiReconParams.LRTGVModelParams.Plot_it;%=50;%50;
CorrB0Map_it=mrsiReconParams.LRTGVModelParams.CorrB0Map_it;%=25 % 1E9 for Synthetic Data  %25 %invivo
CorrB0Map_Maxcount=mrsiReconParams.LRTGVModelParams.CorrB0Map_Maxcount;%=100;%6

Orthogonalize_it=mrsiReconParams.LRTGVModelParams.Orthogonalize_it;%=999
SpecItFact=mrsiReconParams.LRTGVModelParams.SpecItFact;%=5
reduction =mrsiReconParams.LRTGVModelParams.reduction;%= 10^(-2)%10^(-3)

min_SpectStep=mrsiReconParams.LRTGVModelParams.min_SpectStep;%=1E-3
max_SpectStep=mrsiReconParams.LRTGVModelParams.max_SpectStep;%=1/2

min_taup=mrsiReconParams.LRTGVModelParams.min_taup;%=1E-3
max_taup=mrsiReconParams.LRTGVModelParams.max_taup;%=1/8 % %1/8 invivo , 1/16 if diverge, 1/16 for Synthetic Data

CorrB0Map_count=1;

BMask_rr1=single(reshape(mrsiReconParams.BrainMask2D,[size(mrsiReconParams.BrainMask2D) 1]));
SENSE_crr1=single(reshape(mrsiReconParams.SENSE ,[size(mrsiReconParams.SENSE) 1]));
kmask_1kk1=single(reshape(mrsiReconParams.kmask,[1 size(mrsiReconParams.kmask) 1]));
HannF_1kk1=single(reshape(mrsiReconParams.HKernel,[1 size(mrsiReconParams.HKernel) 1]));


SENSE_crr1=SENSE_crr1.* reshape(BMask_rr1,[1 size(BMask_rr1)]); % Mask the Sensitivity profiles to the Brain
U_rrc=BMask_rr1.*U_rrc;% Mask the Spatial components to the Brain

Init_U_rrc= U_rrc;
[ M N NbComp] = size(U_rrc); % numSamplesOnSpoke, numSamplesOnSpoke, nCh
NbCoil = size(Data_ckkt,1);
NbT= size(Data_ckkt,4);
SizeVol = size(U_rrc);
DimVol = ndims(U_rrc)-1;
UIndcs  = repmat({':'}, [1, numel(SizeVol)]);
numSpatialPts=size(U_rrc, 1)*size(U_rrc,2);

%for Forward Transform

FreqMap=mrsiReconParams.WaterFreqMap ;
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;

Time_rrt=single(permute(repmat(([0 :(NbT-1)]'/Fs),[1 M N]),[2,3,1]));
Freqshift_rrt=single(exp(2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT])));
RepData_crrt=single(0*Data_ckkt);
CompData_ckkt=single(0*Data_ckkt);

HomCorr_1rr1=( reshape(mrsiReconParams.HomCorr,[1 size(mrsiReconParams.HomCorr) 1]));
G=single(0*formTensorProduct(U_rrc, V_tc*S,DimVol));

v = single(0*Data_ckkt);%zeros([NbCoil SizeVol]); % coil-k-k-c
p = single(zeros([SizeVol,2]));
q = single(zeros([SizeVol,3]));
xi = single(zeros([SizeVol,2]));
xi_ = single(xi); % v in article

StepNormGrad = [];
StepDiff=Threshold;
StepDiffu = [];
StepDiffr = [];
StepDiffxi = [];
StepDiffq = [ ];
StepDiffp = [];
StepDiffww = [ ];
StepDiffdivp = [];

tau_p = min_taup;%Synthetic Data *1/4 %1/128;%64;%%1/64;%1/16; %the lowest the most stable against singular point
tau_d = tau_p*2;%;1/8;%1/64;%1/8;
stepSize = min_SpectStep; %Synthetic Data *1/4

PrevNormGradU=1E12;
NormGradWentDown=0;
NbDivergSteps=0;

RelDiffSq_V=1;
RelDiffSq_U=1;

t_old    = 2;
PrevNormGrad=0;
step_noConv=0;

u = U_rrc;
uold = u; xiold = xi;
V_old    = V_tc;
V_proj    = V_tc;

G = formTensorProduct(U_rrc, V_tc*S,DimVol); %r r t   
RepData_crrt = reshape(Freqshift_rrt.*G,[1,size(G)]);% c-r-r-t
RepData_crrt = fft(fft(SENSE_crr1.*HomCorr_1rr1.*RepData_crrt,[],2),[],3).*kmask_1kk1;% c-k-k-t
RepData_crrt = (HannF_1kk1.*(Data_ckkt-RepData_crrt));
dataConsistencyCost = norm(RepData_crrt(:));
costFunVal = dataConsistencyCost;

k=-1;

alpha00 = alpha0/reduction;
alpha10 = alpha1/reduction;
alpha002 = alpha0*reduction;
alpha102 = alpha1*reduction;
alpha01 = alpha0;
alpha11 = alpha1;

fftw('wisdom',[]);
fftw('dwisdom',[]);
fftw('planner','measure');

while (k<maxits) %( abs(StepDiff(end))>Threshold  & (k<maxits) ) | k<(minits); 
    k=k+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPATIAL CONVERGENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   if k<=(minits/2)
  	 alpha0 = exp(k/(minits/2)*log(alpha01) + ((minits/2)-k)/(minits/2)*log(alpha00));
  	 alpha1 = exp(k/(minits/2)*log(alpha11) + ((minits/2)-k)/(minits/2)*log(alpha10));
  elseif k<=(minits)
   	 alpha0 = exp((k-minits/2)/(minits/2)*log(alpha01) + ((minits/2)-(k-minits/2))/(minits/2)*log(alpha002));
  	 alpha1 = exp((k-minits/2)/(minits/2)*log(alpha11) + ((minits/2)-(k-minits/2))/(minits/2)*log(alpha102));
  else   
	 alpha0 = alpha01;
  	 alpha1 = alpha11;
  end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE VARIABLES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uold = u;
    xiold = xi;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DUAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    G  = formTensorProduct(U_rrc, V_tc*S,DimVol); %rrt
    
    %explicit Forward Transformation
    RepData_crrt=reshape(Freqshift_rrt.*G,[1,size(G)]);% c-r-r-t
    CompData_ckkt = fft(fft(SENSE_crr1.*HomCorr_1rr1.*RepData_crrt,[],2),[],3).*kmask_1kk1 ;
    
    r =  HannF_1kk1.*(CompData_ckkt -  Data_ckkt); % Data Fidelity gradient % Dimensions % coil-k-k-c
    v = r;
    
    % gradient
    ux = dxp(U_rrc);
    uy = dyp(U_rrc);
    
    
    p(UIndcs{:},1) = p(UIndcs{:},1) - tau_d*(ux + xi_(UIndcs{:},1));
    p(UIndcs{:},2) = p(UIndcs{:},2) - tau_d*(uy + xi_(UIndcs{:},2));

    % projection
    
    for comp=1:NbComp
        absp = sqrt(abs(p(UIndcs{1:(end-1)},comp,1)).^2 + abs(p(UIndcs{1:(end-1)},comp,2)).^2 );
        denom = max(1,absp/(alpha1*max(diag(S))/S(comp,comp)));
        p(UIndcs{1:(end-1)},comp,1) = p(UIndcs{1:(end-1)},comp,1)./denom;
        p(UIndcs{1:(end-1)},comp,2) = p(UIndcs{1:(end-1)},comp,2)./denom;
    end
    
    % symmetrized gradient
    gradxi1 = dxm(xi_(UIndcs{:},1));
    gradxi2 = dym(xi_(UIndcs{:},2));
    gradxi3 = (dym(xi_(UIndcs{:},1)) + dxm(xi_(UIndcs{:},2)))/2;
    
    q(UIndcs{:},1) = q(UIndcs{:},1) - tau_d*gradxi1; % line
    q(UIndcs{:},2) = q(UIndcs{:},2) - tau_d*gradxi2;
    q(UIndcs{:},3) = q(UIndcs{:},3) - tau_d*gradxi3;
    
    % projection
    for comp=1:NbComp
        absq = sqrt(abs(q(UIndcs{1:(end-1)},comp,1)).^2 + abs(q(UIndcs{1:(end-1)},comp,2)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,3)).^2);
        denom = max(1,absq/(alpha0*max(diag(S))/S(comp,comp)));
        q(UIndcs{1:(end-1)},comp,1) = q(UIndcs{1:(end-1)},comp,1)./denom;
        q(UIndcs{1:(end-1)},comp,2) = q(UIndcs{1:(end-1)},comp,2)./denom;
        q(UIndcs{1:(end-1)},comp,3) = q(UIndcs{1:(end-1)},comp,3)./denom;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRIMAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % dual operator
    
    Diff_crrt=ifft(ifft( kmask_1kk1.*v,[],2),[],3);
    ww =  conj(SENSE_crr1).*Diff_crrt; % coil -r - r -c
    ww = squeeze(sum(ww,1)); %r - r -t
    
    ww= formTensorProduct(conj(Freqshift_rrt).*ww.*BMask_rr1 ,(V_tc*inv(S))', DimVol ) ;
    
    % divergence
    divp = dxm(p(UIndcs{:},1)) + dym(p(UIndcs{:},2));
    
    u = u - tau_p*(ww + divp); %lines 8
  
    % divergence
    divq1 = dxp(q(UIndcs{:},1)) + dyp(q(UIndcs{:},3));
    divq2 = dxp(q(UIndcs{:},3)) + dyp(q(UIndcs{:},2));
    
    xi(UIndcs{:},1) = xi(UIndcs{:},1) - tau_p*(divq1 - p(UIndcs{:},1));%line 11
    xi(UIndcs{:},2) = xi(UIndcs{:},2) - tau_p*(divq2 - p(UIndcs{:},2));%line 11
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUXILIARY UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


    U_rrc = (2*u - uold).*BMask_rr1; %Force the values to zero outside the Brain mask (avoid wave behavior outside brain)
    xi_ = 2*xi - xiold;
    

   if  mod(k+1,Orthogonalize_it)==0
            U_rrc = reshape(OrthogonalizeComponents(reshape(U_rrc,[size(U_rrc,1)*size(U_rrc,2) , size(U_rrc,3)])),size(U_rrc));
	    u = reshape(OrthogonalizeComponents(reshape(u,[size(u,1)*size(u,2) , size(u,3)])),size(u));
   end

    for c=1:NbComp
        NormUrrc=sqrt(sum(sum(abs(U_rrc(:,:,c)).^2,1),2));
        S(c,c)=S(c,c)*NormUrrc;
        U_rrc(:,:,c) = U_rrc(:,:,c)/NormUrrc;
        u(:,:,c) = u(:,:,c)/NormUrrc;
    end

    RelDiffSq_U = norm(U_rrc(:)-u(:))^2/norm(u(:))^2 ;%+ norm(xi_(:)-xi(:))^2/norm(xi(:))^2;
    norm_p=norm(p(:),1);
    norm_q=norm(q(:),1);
    RelDiffSq_old=RelDiffSq_U;
    
    NormGradU=norm(ww(:) + divp(:));
      
    if k>10 & NormGradU<PrevNormGradU;
        NormGradWentDown=1;
        tau_p = tau_p*1.10;
        if tau_p>max_taup; tau_p=max_taup;end;
    elseif NormGradWentDown==1;
        tau_p=tau_p*0.85;%0.2;
        if tau_p<min_taup;tau_p=min_taup;end;
    end
    tau_d=tau_p*2;
    PrevNormGradU=NormGradU;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPECTRAL CONVERGENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    for SpecInIt=1:SpecItFact

        G = formTensorProduct(U_rrc, V_tc*S,DimVol); %r r t
        
	RepData_crrt=reshape(Freqshift_rrt.*G,[1,size(G)]);% c-r-r-t
	RepData_crrt = HannF_1kk1 .*(fft(fft(SENSE_crr1.*HomCorr_1rr1.*RepData_crrt,[],2),[],3).*kmask_1kk1- Data_ckkt);% c-k-k-t
	G=BMask_rr1.*conj(Freqshift_rrt).*squeeze(sum(conj(SENSE_crr1).*ifft(ifft(RepData_crrt,[],2),[],3),1)); % r-r-t;
	gradStep = ((reshape(U_rrc, [], NbComp)*inv(S))' * reshape(G , [numSpatialPts, NbT]))' ;

	[ gradStep] = OrthogonalizeComponents(gradStep);

	V_proj   = V_tc -stepSize * gradStep;

	for c=1:size(V_tc,2)
		V_proj(:,c) = V_proj(:,c) ./ norm(V_proj(:,c));
	end

	if k>0 
		RelDiffSq_V=norm(V_proj(:)-V_old(:))/norm(V_proj(:));
	end
	NormGrad=norm(gradStep(:));
	if (NormGrad<PrevNormGrad || k==0)
		step_noConv=0;
		stepSize = 1.2 * stepSize;
		if stepSize>max_SpectStep;stepSize=max_SpectStep;end
	else
		step_noConv=step_noConv+1;
		stepSize = 0.75 * stepSize;
		if stepSize<min_SpectStep;stepSize=min_SpectStep;end
	end

	PrevNormGrad=NormGrad;
	StepNormGrad = [StepNormGrad NormGrad];

	t_new    = (1 + sqrt(1 + 4*t_old.^2)) / 2;
	fact=((t_old - 1) / t_new);
	% V_update = V + ((t_old - 1) / t_new) * (V - V_old);
	V_tc = (1-fact)*V_old +(fact)* V_proj;

	for c=1:size(V_tc,2)
	V_tc(:,c) = V_tc(:,c) ./ norm(V_tc(:,c));
	end
	V_old    = V_proj;
	t_old    = t_new;
          
end
    
if  mod(k+1,Orthogonalize_it)==0
	V_tc = OrthogonalizeComponents(V_tc);
	V_old= OrthogonalizeComponents(V_old);
end


  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FREQUENCY MAP DYNAMIC CORRECTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
     if mod(k+1,CorrB0Map_it) == 5 & CorrB0Map_count<=CorrB0Map_Maxcount & k>100
	CorrB0Map_count=CorrB0Map_count+1;
        A=6;
        if(size(V_tc,2)<A) %% ADDED TO AVOID MISTAKES WHEN NCOMPS IS LOWER THAN 6 (BA-02/11/2022)
            A=size(V_tc,2);
        end
        RefTimeSerie=conj(sum((V_tc(:,1:A)*S(1:A,1:A)),2));
        [FreqMapCorr, CCoefMap,RefSpectrumShort ] = MeasureFreqMap(RefTimeSerie, permute(formTensorProduct(U_rrc, V_tc*S,DimVol),[3,1,2]),mrsiReconParams);
        MFreqMapCorr=mean(abs(FreqMapCorr(mrsiReconParams.BrainMask2D>0)));
        MFreqMap=mean(abs(FreqMap(mrsiReconParams.BrainMask2D>0)));
        RelFreqCorr=MFreqMapCorr/MFreqMap;
	fprintf('Frequency Map correction Nb.%g : RelFreqCorr = %g, MeanFreqCorr= %g, MeanFreqMap= %g, it = %g \n', CorrB0Map_count-1,RelFreqCorr,MFreqMapCorr,MFreqMap,k+1);
      
        if(RelFreqCorr>0.001)

            FreqMap=FreqMap-FreqMapCorr.*mrsiReconParams.BrainMask2D;
            FreqMap=FreqMap-mrsiReconParams.BrainMask2D.*mean(FreqMap(mrsiReconParams.BrainMask2D>0));% avoid freq. shifting
            Freqshift_rrt=exp(2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT]));
            
            s=[ './',mrsiReconParams.Log_Dir,'/LowRankTGV_Recon/',mrsiReconParams.NameData,'_FreShiftCorr_step', num2str(k+1), '.ps'];
            if exist(s);delete(s);end
            figs=figure('visible','off');
            imagesc(FreqMapCorr);colorbar;
            title('Frequency Map correction');
            print(figs, '-append', '-dpsc2',s);
            imagesc(FreqMap,[-30 30]);colorbar;
            title('Resulting Frequency Map');
            print(figs, '-append', '-dpsc2',s);
            
            imagesc(CCoefMap);colorbar;
            title('Correlation Map');
            print(figs, '-append', '-dpsc2',s);
             plot(1:numel(RefSpectrumShort),real(RefSpectrumShort),1:numel(RefSpectrumShort),imag(RefSpectrumShort));
            title('Reference Time serie');
            print(figs, '-append', '-dpsc2',s);
            close
        end
     end
    

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keep other norms for monitoring
    
    norm_xi=norm(xi_(:));
    norm_u=norm(U_rrc(:));
    norm_r=norm(NormGradU(:));
    norm_divp=norm(divp(:));
    norm_ww=norm(ww(:));
    
    
    
    StepDiffu = [ StepDiffu norm_u];
    StepDiffr = [ StepDiffr norm_r];
    StepDiffxi = [ StepDiffxi norm_xi];
    StepDiffq = [ StepDiffq norm_q];
    StepDiffp = [ StepDiffp norm_p];
    StepDiffww = [ StepDiffww norm_ww ];
    StepDiffdivp = [ StepDiffdivp norm_divp ];
    
    if mod(k+1,check_it) == 0
        
	G = formTensorProduct(U_rrc, V_tc*S,DimVol); %r r t   
	RepData_crrt = reshape(Freqshift_rrt.*G,[1,size(G)]);% c-r-r-t
	RepData_crrt = fft(fft(SENSE_crr1.*HomCorr_1rr1.*RepData_crrt,[],2),[],3).*kmask_1kk1;% c-k-k-t
	RepData_crrt = HannF_1kk1.*(Data_ckkt-RepData_crrt);
	dataConsistencyCost = norm(RepData_crrt(:));
        StepDiff = [ StepDiff, (dataConsistencyCost-costFunVal(end))/(dataConsistencyCost*check_it)];
        costFunVal = [costFunVal, dataConsistencyCost];

        %fprintf('\nTGV2-L2-2D-PD: it = %g, costFunVal = %g,tau_p = %g, RelDiffSq = %g ', k+1,costFunVal(end) ,taup_c(1), RelDiffSq_U);
        fprintf('\nTGV2-L2-2D-PD: it = %g, costFunDiff = %g,tau_p = %g, RelDiffSq_U = %g ', k+1,StepDiff(end) ,tau_p, RelDiffSq_U);
        fprintf('\nSpectral Iteration: %d , step_size = %d, RelDiffSq_V = %d, step_noConv = %d', k+1, stepSize,RelDiffSq_V,step_noConv);
        CovV=V_tc'*V_tc;
        U_rc =reshape(U_rrc,[],size(U_rrc,3));
        CovU=(U_rc'*U_rc);
        fprintf('\nTime component independance: %d , Spatial component independance: %d\n', trace(abs(CovV))/sum(abs(CovV(:))),  trace(abs(CovU))/sum(abs(CovU(:))) );

    
        if mod(k+1,Plot_it) == 0
		s=[ './',mrsiReconParams.Log_Dir,'/LowRankTGV_Recon/',mrsiReconParams.NameData,'_TGV_Diagnostic.ps'];
		if exist(s);delete(s);end

		imagesc(abs(U_rrc(:,:,2)),[0 5*nanmean(nanmean(abs(U_rrc(:,:,2))))]);drawnow;
		print(s,'-bestfit', '-append', '-dpsc2'); 
		if (NbComp>3);imagesc(abs(U_rrc(:,:,4)),[0 5*nanmean(nanmean(abs(U_rrc(:,:,4))))]);drawnow;end;
		print(s,'-bestfit', '-append', '-dpsc2'); 
		if (NbComp>7); imagesc(abs(U_rrc(:,:,8)),[0 5*nanmean(nanmean(abs(U_rrc(:,:,8))))]);drawnow;end;
		print(s,'-bestfit', '-append', '-dpsc2'); 
		if (NbComp>15); imagesc(abs(U_rrc(:,:,16)),[0 5*nanmean(nanmean(abs(U_rrc(:,:,16))))]);drawnow;end;
		print(s,'-bestfit', '-append', '-dpsc2'); 
		if (NbComp>19); imagesc(abs(U_rrc(:,:,20)),[0 5*nanmean(nanmean(abs(U_rrc(:,:,20))))]);drawnow;end;
		print(s,'-bestfit', '-append', '-dpsc2');         
		plot(StepDiffu);title('norm U');drawnow;
		print(s,'-bestfit', '-append', '-dpsc2'); 
		plot(StepDiffxi);title('norm Xi');drawnow;
		print(s,'-bestfit', '-append', '-dpsc2'); 
		plot(StepDiffp);title('norm p');drawnow;
		print(s,'-bestfit', '-append', '-dpsc2'); 
		plot(StepDiffq);title('norm q');drawnow;
		print(s,'-bestfit', '-append', '-dpsc2');
		plot(costFunVal);title('Data Fidelity');drawnow;
		print(s,'-bestfit', '-append', '-dpsc2');
		plot(StepDiff);title('Diff Data Fidelity');drawnow;
		print(s,'-bestfit', '-append', '-dpsc2'); 
        end
        
    end
    
    if  mod(k+1,Plot_it)==0
        VisualizeTGV( Init_U_rrc ,U_rrc,['./',mrsiReconParams.Log_Dir,'/LowRankTGV_Recon/',mrsiReconParams.NameData,'_muTV', num2str(mrsiReconParams.mu_tv),'_step', num2str(k+1)]);
        VisualizeSpectral( V_tc,S, [ './',mrsiReconParams.Log_Dir,'/LowRankTGV_Recon/',mrsiReconParams.NameData,'_step', num2str(k+1)])
            end
   close all;
end

%Reorder Component following signular value:

[SDiag, DescOrder]=sort(diag(S),'descend');
U_rrc=U_rrc(:,:,DescOrder).*repmat(mrsiReconParams.BrainMask2D,[1 1 NbComp]);
V_tc=V_tc(:,DescOrder);
S=diag(SDiag);
fprintf([ '\nTGV Recon done in ', num2str(k),' steps. Relative Step Diff U =',  num2str(RelDiffSq_U), ' steps. Relative Step Diff V =',  num2str(RelDiffSq_V),'\n']);

