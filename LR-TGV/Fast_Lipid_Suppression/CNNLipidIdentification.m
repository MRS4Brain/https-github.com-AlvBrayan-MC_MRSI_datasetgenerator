function [ mrsiDataLR_ctkk, mrsiDataLipIDOp_ckkt ]  = CNNLipidIdentification(  LipidProj_cff,Lipids_tkk, mrsiReconParams ,NameData)

% mrsiReconParams.mrsiData dims: time-k-k

MinPPM_Lip= -4.7; %mrsiReconParams.MinPPM ;
MaxPPM_Lip=+2; %mrsiReconParams.MaxPPM ;
NbTrainEx=1000000;
MaxLipScaling=500;
NbTrainEpoch=60;
CNNPadding = 8;
dropout = 0.08; %0.05 for NFilter=32
NFilters = 32;
Nbatch = 256;
MaxPeakW=50;

curDir = pwd;
NameStackLipDataPath=[curDir,'/',mrsiReconParams.Log_Dir,'/LipidSuppression/Lipid_Stack_File.h5'];
NameLipOpPath=[curDir,'/',mrsiReconParams.Log_Dir,'/LipidSuppression/LipidProjOpFile.h5'];

CNNTrainingDataPath=[mrsiReconParams.Log_Dir,'/LipidSuppression/CNNLipTrainingData.h5'];
CNNLipModelPath=[mrsiReconParams.Log_Dir,'/LipidSuppression/CNNLipIDModel'];


NameStackInputAllDataPath=[mrsiReconParams.Log_Dir,'/LipidSuppression/MRSIDataWithLip_Stack_File.h5'];
NameStackInputLipOpIDDataPath=[mrsiReconParams.Log_Dir,'/LipidSuppression/MRSIDataLipIDByProjOp_Stack_File.h5'];
NameStackOutputDataPath=[mrsiReconParams.Log_Dir,'/LipidSuppression/MRSIDataLipRemovedByCNN_Stack_File.h5'];

[~,MaxPPM_pt_In]=min(abs(MaxPPM_Lip- mrsiReconParams.ppm));
[~,MinPPM_pt_In]=min(abs(MinPPM_Lip  - mrsiReconParams.ppm));
MaxPPM_pt_In=MaxPPM_pt_In - mod((MaxPPM_pt_In-MinPPM_pt_In+1),16)+16;
WindowSizeIn=(MaxPPM_pt_In-MinPPM_pt_In+1);

MaxPPM_pt_Out=MaxPPM_pt_In -CNNPadding;
MinPPM_pt_Out=MinPPM_pt_In+CNNPadding;
WindowSizeOut=(MaxPPM_pt_Out-MinPPM_pt_Out+1);
fprintf([ 'Lipid Suppression by CNN over the Spectral Range: ' num2str(MinPPM_pt_Out), ' - ', num2str(MaxPPM_pt_Out), ' pts\n']);
[SCRIPT_DIR, ~, ~]=fileparts(mfilename('fullpath'));



mrsiDataLR_ctkk = single(0*mrsiReconParams.mrsiData);

if ~isempty(NameData)
    N=size(mrsiReconParams.mrsiData);
    HzpPt=mrsiReconParams.mrProt.samplerate/(mrsiReconParams.mrProt.VSize);
    HistLR_r=[];Hist_r=[];
    HistNoise_r=[];HistFreq_r=[];

    for coil=1:N(1)   
       Data_rrf=fft(squeeze(permute(single(mrsiReconParams.mrsiData(coil,:,:,:)),[3,4,2,1])),[],3);% k -k -k-f
       Data_rrf = ifft(ifft(Data_rrf,[],1),[],2);
       Hist_r=[Hist_r ; reshape(std(Data_rrf(:,:,50:300),0,4),[],1)];
       [~, I ] = max(Data_rrf(:,:,50:300),[],4);
       HistFreq_r =[HistFreq_r ; HzpPt*reshape(I,[],1)];
       Data_rrf= reshape(reshape(Data_rrf,[],N(2)) * single(eye(N(2))-squeeze(LipidProj_cff(coil,:,:))),[N(3),N(4),N(2)]);
       HistLR_r=[HistLR_r ; reshape(std(Data_rrf(:,:,50:250),0,4),[],1)];
       HistNoise_r=[HistNoise_r ; reshape(std(Data_rrf(:,:,600:700),0,4),[],1)];
       mrsiDataLR_ctkk(coil,:,:,:)= permute(ifft(fft(fft(Data_rrf,[],1),[],2),[],3),[3,1,2]);    
    end

      
    %save( 'Data_LinLipRM_ctkkk', 'mrsiDataLR_ctkkk','-v7.3');
    HistFreq_r = HistFreq_r(abs(HistNoise_r)>0);
    HistFreq_r = HistFreq_r-median(HistFreq_r(:));
    HistLR_r = HistLR_r(abs(HistNoise_r)>0);
    Hist_r = Hist_r(abs(HistNoise_r)>0);
    HistNoise_r = HistNoise_r(abs(HistNoise_r)>0);
    
    s=['./',mrsiReconParams.Log_Dir,'/',NameData,'_Histograms.ps'];
    delete(s);
    
    figs=figure('visible', 'off');
    hist(abs(HistLR_r)./HistNoise_r,256);
    title('Metabolite SNR');
    print(figs, '-append', '-dpsc2', s);

    hist(abs(Hist_r)./HistNoise_r,256);
    title('Lipids SNR');
    print(figs, '-append', '-dpsc2', s);

    hist(abs(Hist_r)./abs(HistLR_r),256);
    title('Lipid/Metab');
    print(figs, '-append', '-dpsc2', s);

    hist(HistFreq_r,256);
    title('Lipid Freq. shift in Hz');
    print(figs, '-append', '-dpsc2', s);
     
    close(figs);
    
    MRSIDataOrig_rrf=squeeze(sum( conj(mrsiReconParams.SENSE).*fft(ifft(ifft(permute(mrsiDataLR_ctkk,[1,3,4,2]),[],2),[],3),[],4),1));


    clear mrsiDataLR_ctkkk
    SizeD=size(MRSIDataOrig_rrf);
    s1=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpecComp_BeforeCNNSup.ps'];
    s2=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpatComp_BeforeCNNSup.ps'];
    delete(s1);delete(s2);
    [Uorig,Sorig,Vorig] = svd(reshape(MRSIDataOrig_rrf,[],SizeD(3)),0);
    U_rrc=reshape(Uorig(:,1:10), SizeD(1),SizeD(2),[]);
    figs=figure('visible', 'off');
    for comp=1:10
        
        plot(MinPPM_pt_Out:MaxPPM_pt_Out,real(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)),...
            MinPPM_pt_Out:MaxPPM_pt_Out,imag(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)),...
            MinPPM_pt_Out:MaxPPM_pt_Out,abs(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)));
        title(['Spectral comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s1);
        
        plotImage= abs(U_rrc(:,:,comp)) ;
        imagesc(plotImage);
        colormap default;colorbar;
        title(['Spatial comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s2);
        
    end
end   % if ~isempty(NameData)         
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if exist([CNNLipModelPath,'/model.json'])==0 % if the Model doesn't exist, then make it
    
    N = size(mrsiReconParams.mrsiData);
    Npt=N(2);
    N1=MinPPM_pt_In;
    N2=MaxPPM_pt_In;
    NMRfreq=mrsiReconParams.mrProt.NMRFreq*1E6;
    hdf5write(NameLipOpPath,'realLipidProj',real(LipidProj_cff),'imagLipidProj',imag(LipidProj_cff),...
        'samplerate',mrsiReconParams.mrProt.samplerate,...
        'Npt',Npt,'N1',N1-1,'N2',N2,'NMRfreq',NMRfreq);

   [~] = make_CNN_StackData( reshape(Lipids_tkk,[1 size(Lipids_tkk)]),mrsiReconParams.SkMask2D,[1,Npt],mrsiReconParams,NameStackLipDataPath);  

    clear Lipids_ctkkk
    fprintf([ 'Generation of the training dataset ...\n']);
    
   command= [mrsiReconParams.PrefixKeras, ' ',SCRIPT_DIR,'/LipIDCNN/GenerationFID/Gene_FID_withLipids.py -o ',CNNTrainingDataPath, ' -l ',NameStackLipDataPath, ' -lop ',NameLipOpPath,' --ntrain ',num2str(NbTrainEx),' --maxLipSc ',num2str(MaxLipScaling), ' --maxPkW ',num2str(MaxPeakW) , ' -vf'];
   [status,cmdout] = system(command,'-echo');
    if(status>0);fprintf(cmdout);end
    
    fprintf([ 'Training dataset finished. Starting Model training ...\n']);
    
    %command= [mrsiReconParams.PrefixKeras, ' ',SCRIPT_DIR,'/LipIDCNN/Training/train_LipID_UNet.py -o ',CNNLipModelPath, ' -i ',CNNTrainingDataPath, ' --Nepochs ', num2str(NbTrainEpoch)];
    command= [mrsiReconParams.PrefixKeras, ' ',SCRIPT_DIR,'/LipIDCNN/Training/train_LipID_UNet.py -o ',CNNLipModelPath, ' -i ',CNNTrainingDataPath, ' --Nepochs ', num2str(NbTrainEpoch), ' --dropOut ', num2str(dropout), ' --nFilters ', num2str(NFilters), ' --Nbatch ', num2str(Nbatch)];
    [status,cmdout] = system(command);
    if(status>0);fprintf(cmdout);end
    
    if exist([CNNLipModelPath,'/model.json'])==0
        error('CNN Training Failed. Stopping Here!')
    end
    fprintf([ 'Model training finished.\n']);
    delete(NameStackLipDataPath);
    delete(CNNTrainingDataPath);
    delete(NameLipOpPath);
else
    fprintf([ 'LipReModel exists! Skipping the training and proceeding to the Lipid Removal.\n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf([ 'Starting Lipid Removal by CNN ...\n']);

% Store raw data  with lipids in a .h5 file...
%[~] = make_CNN_StackData_InKSpace( mrsiReconParams.mrsiData_ctkkk ,DataMask,[MinPPM_pt_In,MaxPPM_pt_In],mrsiReconParams,NameStackInputAllDataPath);
[~] = make_CNN_StackData_crrf( fft(ifft(ifft(permute(mrsiReconParams.mrsiData,[1,3,4,2]),[],2),[],3),[],4) ,[MinPPM_pt_In,MaxPPM_pt_In],mrsiReconParams,NameStackInputAllDataPath);

% Apply the LipRM operator

[Nc,Nt,Nk1,Nk2] =size(mrsiReconParams.mrsiData);
mrsiDataLipIDOp_ckkt = zeros(Nc,Nk1,Nk2,Nt);
for coil=1:N(1)
    Data_kkf=fft(squeeze(permute(single(mrsiReconParams.mrsiData(coil,:,:,:)),[3,4,2,1])),[],3);
    %Data_kkkf= reshape(reshape(Data_kkf,[],Nt) * single(eye(Nt)-squeeze(LipidProj_cff(coil,:,:))),[N(3),N(4),N(5),N(2)]);
    Data_kkf= reshape(reshape(Data_kkf,[],Nt) * single(squeeze(LipidProj_cff(coil,:,:))),[Nk1,Nk2,Nt]);
    mrsiDataLipIDOp_ckkt(coil,:,:,:)= ifft(Data_kkf,[],3);
end

DataMask = mrsiReconParams.kmask;

% Store ProjOpLipRMData  with lipids in a .h5 file...
%[~] = make_CNN_StackData_InKSpace( permute(mrsiDataLipIDOp_ckkkt,[1,5,2,3,4]) ,DataMask,[MinPPM_pt_In,MaxPPM_pt_In],mrsiReconParams,NameStackInputLipOpIDDataPath);
[~] = make_CNN_StackData_crrf( fft(ifft(ifft(mrsiDataLipIDOp_ckkt,[],2),[],3),[],4) ,[MinPPM_pt_In,MaxPPM_pt_In],mrsiReconParams,NameStackInputLipOpIDDataPath);

% [U_rc,S,V_tc] = svd(reshape(mrsiDataLipIDOp_ckkkt,[],Nt),0);
% V_tc=V_tc(:,1:mrsiReconParams.modelOrder*5);
% S=S(1:mrsiReconParams.modelOrder*5,1:mrsiReconParams.modelOrder*5);
% U_rc = U_rc(:,1:mrsiReconParams.modelOrder*5);
% [~] = make_CNN_StackData( permute(reshape( (U_rc * S * (V_tc') ),[Nc,Nk1,Nk2,Nk3,Nt]),[1,5,2,3,4]) ,DataMask,[MinPPM_pt_In,MaxPPM_pt_In],mrsiReconParams,NameStackInputLipOpIDDataPath);
% ResidueData_ckkkt = mrsiDataLipIDOp_ckkkt - reshape( (U_rc * S * (V_tc') ),[Nc,Nk1,Nk2,Nk3,Nt]);

%clear mrsiDataLipIDOp_ckkkt

fprintf([ 'Using model ',CNNLipModelPath, ' for lipid ID\n']);

command= [mrsiReconParams.PrefixKeras, ' ',SCRIPT_DIR,'/LipIDCNN/Training/Run_model_on_LipData.py -m ',CNNLipModelPath, ' -i1 ',NameStackInputAllDataPath, ' -i2 ', NameStackInputLipOpIDDataPath,' -o ',NameStackOutputDataPath];

% command= [mrsiReconParams.PrefixKeras, ' ',SCRIPT_DIR,'/LipIDCNN/Training/Run_LinModel_on_LipData.py -m ',CNNLipModelPath, ' -i1 ',NameStackInputAllDataPath, ' -i2 ', NameStackInputLipOpIDDataPath,' -o ',NameStackOutputDataPath];

[status,cmdout] = system(command,'-echo');
if(status>0);fprintf(cmdout);end

ReData_rf = h5read(NameStackOutputDataPath,'/LipData_rf/realData' );
ImData_rf = h5read(NameStackOutputDataPath,'/LipData_rf/imagData' );
ReData_rf = ReData_rf(:,(1+CNNPadding):(end-CNNPadding));
ImData_rf = ImData_rf(:,(1+CNNPadding):(end-CNNPadding));
SizeD=size(mrsiReconParams.mrsiData);

%Mask_crrr=permute(repmat(DataMask,[ 1 1 1 SizeD(1)]),[4 1 2 3]);
%TempShort_crrrf(Mask_crrr(:)>0,:)=ReData_rf+1j*ImData_rf;
TempShort_crrf=zeros([SizeD(1)*SizeD(3)*SizeD(4),WindowSizeOut]);
TempShort_crrf=ReData_rf+1j*ImData_rf;

TempShort_crrf =reshape(TempShort_crrf,[SizeD(1),SizeD(3),SizeD(4),WindowSizeOut]);
LipShort_cfkk = fft(fft( permute(TempShort_crrf,[1,4,2,3]) ,[],3),[],4); %Now really : cfkkk

%Data are in k-space already

delete(NameStackInputAllDataPath);
delete(NameStackInputLipOpIDDataPath);
delete(NameStackOutputDataPath)

data_fkk=zeros(SizeD(2:4));

mrsiDataLR_ctkk = single(0*mrsiReconParams.mrsiData);
N=size(mrsiReconParams.mrsiData);
for c=1:SizeD(1)

    AllInShort_fkk=squeeze(fft(single(mrsiReconParams.mrsiData(c,:,:,:)),[],2));
    AllInShort_fkk = AllInShort_fkk(MinPPM_pt_Out:MaxPPM_pt_Out,:,:);
    Data_kkf=fft(squeeze(permute(single(mrsiReconParams.mrsiData(c,:,:,:)),[3,4,2,1])),[],3);
    Data_kkf= reshape(reshape(Data_kkf,[],N(2)) * single(eye(N(2))-squeeze(LipidProj_cff(c,:,:))),[N(3),N(4),N(2)]);
    mrsiDataLR_ctkk(c,:,:,:)= permute(ifft(Data_kkf,[],3),[3,1,2]);
    
    fprintf([ 'Storing Lipid-Free data in coil element ',num2str(c) , '.\n']);
    
    data_fkk=squeeze(fft(mrsiDataLR_ctkk(c,:,:,:),[],2));
  
    data_fkk(MinPPM_pt_Out:MaxPPM_pt_Out,:,:) = AllInShort_fkk - squeeze(LipShort_cfkk(c,:,:,:)) ;
  
    mrsiDataLR_ctkk(c,:,:,:) = ifft(data_fkk,[],1);
    
    %for plotting only:
    DataLRShort_frr=squeeze(ifft(ifft(AllInShort_fkk - squeeze(LipShort_cfkk(c,:,:,:)),[],2),[],3));
 
 

    if ~isempty(NameData)
        LipidShort_frr=squeeze( ifft(ifft( fft( single(mrsiReconParams.mrsiData(c,:,:,:)),[],2),[],3),[],4) );
        LipidShort_frr = LipidShort_frr(MinPPM_pt_Out:MaxPPM_pt_Out,:,:);
        s=['./',mrsiReconParams.Log_Dir,'/',NameData,'_coil',num2str(c) , '_Lipid_Images.ps'];
        if exist(s);delete(s);end
        close all
        figs=figure('visible', 'off');
        
        
        plotImage= squeeze(sum(abs(LipidShort_frr),1));
        imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
        colormap default;colorbar;
        title('Data with Lipids in Whole Volume');
        print(figs, '-append', '-dpsc2', s);
        
        plotImage= squeeze(sum(abs(DataLRShort_frr),1));
        imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
        colormap default;colorbar;
        title('Lipid Removed Data in Whole Volume');
        print(figs, '-append', '-dpsc2', s);
        
        plotImage= squeeze(sum(abs(LipidShort_frr-DataLRShort_frr),1));
        imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
        colormap default;colorbar;
        title('Lipid-Identified Data in Whole Volume');
        print(figs, '-append', '-dpsc2', s);
        
        
        plotImage= (~mrsiReconParams.SkMask2D).*squeeze(sum(abs(LipidShort_frr),1));
        imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
        colormap default;colorbar;
        title('Data with Lipids in Brain & Outside Head');
        print(figs, '-append', '-dpsc2', s);
        
        plotImage=((~mrsiReconParams.SkMask2D).*squeeze(sum(abs(DataLRShort_frr),1)) );
        imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
        colormap default;colorbar;
        title('Lipid Removed Data in Brain & Outside Head');
        print(figs, '-append', '-dpsc2', s);
        
        plotImage= ((~mrsiReconParams.SkMask2D).*squeeze(sum(abs(LipidShort_frr-DataLRShort_frr),1)));
        imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
        colormap default;colorbar;
        title('Lipid-Identified Data in Brain & Outside Head');
        print(figs, '-append', '-dpsc2', s);
        
        plotImage= (mrsiReconParams.SkMask2D );
        imagesc(plotImage);%,[ 0, 10*mean(image2plot(:))] )
        colormap default;colorbar;
        title('Lipid mask');
        print(figs, '-append', '-dpsc2', s);
    end
    
end

   %save( 'Data_NNLipRM_ctkkk', 'mrsiDataLR_ctkkk','-v7.3');


MRSIDataLR_rrf=squeeze(sum( conj(mrsiReconParams.SENSE).*fft(ifft(ifft(permute(mrsiDataLR_ctkk,[1,3,4,2]),[],2),[],3),[],4),1));
MRSIDataLR_rrf = MRSIDataLR_rrf.*reshape(mrsiReconParams.BrainMask2D,[size(mrsiReconParams.BrainMask2D) 1]);
SizeD=size(MRSIDataLR_rrf)
if ~isempty(NameData)
    s1=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpecComp_AfterCNNSup.ps'];
    s2=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpatComp_AfterCNNSup.ps'];
    delete(s1);delete(s2);
    [Uorig,Sorig,Vorig] = svd(reshape(MRSIDataLR_rrf,[],SizeD(3)),0);
    U_rrc=reshape(Uorig(:,1:10), SizeD(1),SizeD(2),[]);
    figs=figure('visible', 'off');
    for comp=1:10
        
        plot(MinPPM_pt_Out:MaxPPM_pt_Out,real(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)),...
            MinPPM_pt_Out:MaxPPM_pt_Out,imag(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)),...
            MinPPM_pt_Out:MaxPPM_pt_Out,abs(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)))
            ;
        title(['Spectral comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s1);
        
        plotImage= (abs(U_rrc(:,:,comp)) );
        imagesc(plotImage);
        colormap default;colorbar;
        title(['Spatial comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s2);
        
    end
end     
fprintf([ 'Lipid Removing by CNN finished.\n']);





end

