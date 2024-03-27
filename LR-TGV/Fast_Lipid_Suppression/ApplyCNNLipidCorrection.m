function mrsiDataLR_frrr = ApplyCNNLipidCorrection( Data_frrr, mrsiReconParams ,NameData)

% mrsiReconParams.mrsiData dims: time-k-k

MinPPM_Lip = mrsiReconParams.MinPPM% -4.4; %mrsiReconParams.MinPPM ;
MaxPPM_Lip = mrsiReconParams.MaxPPM%+1; %mrsiReconParams.MaxPPM ;
CNNPadding = 8;

[~,MaxPPM_pt_In]=min(abs(MaxPPM_Lip- mrsiReconParams.ppm));
[~,MinPPM_pt_In]=min(abs(MinPPM_Lip  - mrsiReconParams.ppm));
MaxPPM_pt_In=MaxPPM_pt_In - mod((MaxPPM_pt_In-MinPPM_pt_In+1),16)+16;
WindowSizeIn=(MaxPPM_pt_In-MinPPM_pt_In+1);

MaxPPM_pt_Out=MaxPPM_pt_In -CNNPadding;
MinPPM_pt_Out=MinPPM_pt_In + CNNPadding;
WindowSizeOut=(MaxPPM_pt_Out-MinPPM_pt_Out+1);

fprintf([ 'Lipid Suppression Correction by CNN over the Spectral Range: ' num2str(MinPPM_pt_Out), ' - ', num2str(MaxPPM_pt_Out), ' pts\n']);
[SCRIPT_DIR, ~, ~]=fileparts(mfilename('fullpath'));



MRSIDataOrig_rrrf=permute(Data_frrr,[2,3,4,1]);
MRSIDataOrig_rrrf = MRSIDataOrig_rrrf(:,:,:,MinPPM_pt_In:MaxPPM_pt_In);
SizeD=size(MRSIDataOrig_rrrf);
if ~isempty(NameData)
    
    s1=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpecComp_BeforeCNNCorr.ps'];
    s2=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpatComp_BeforeCNNCorr.ps'];
    delete(s1);delete(s2);
    [Uorig,Sorig,Vorig] = svd(reshape(MRSIDataOrig_rrrf.*mrsiReconParams.BrainMask,[],SizeD(4)),0);
    U_rrrc=reshape(Uorig(:,1:10), SizeD(1),SizeD(2),SizeD(3),[]);
    figs=figure('visible', 'off');
    for comp=1:10
        
        plot(MinPPM_pt_In:MaxPPM_pt_In,real(Vorig(:,comp)),...
            MinPPM_pt_In:MaxPPM_pt_In,imag(Vorig(:,comp)),...
            MinPPM_pt_In:MaxPPM_pt_In,abs(Vorig(:,comp)));
        title(['Spectral comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s1);
        
        plotImage= squeeze(abs(U_rrrc(:,:,:,comp))); %plotImage= Vol2Image(abs(U_rrrc(:,:,:,comp)) );
        imagesc(plotImage);
        colormap default;colorbar;
        title(['Spatial comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s2);
        
    end
end



curDir = pwd;
CNNLipModelPath=[curDir,'/',mrsiReconParams.Log_Dir,'/LipidSuppression/CNNLipCorrModel'];




if exist([CNNLipModelPath,'/model.json'])==0 % if the Model doesn't exist
	error(['CNN model: ',CNNLipModelPath , ' needed for lipid suppression correction doesnt exist. Processing is aborded!'] )
end

fprintf([ 'Starting Correction by CNN ...\n']);

fprintf([ 'Writing the data files ...\n']);

NameStackInputDataPath=[curDir,'/',mrsiReconParams.Log_Dir,'/LipidSuppression/MRSIDataAfterLipRem_Stack_File.h5'];
NameStackOutputDataPath=[curDir,'/',mrsiReconParams.Log_Dir,'/LipidSuppression/MRSIDataLipCorrected_Stack_File.h5'];
DataMask = ones(size(mrsiReconParams.ImMask));

[~] = make_CNN_StackData_frrr(Data_frrr ,DataMask,[MinPPM_pt_In,MaxPPM_pt_In],mrsiReconParams,NameStackInputDataPath);

command= [mrsiReconParams.PrefixKeras, ' ',SCRIPT_DIR,'/LipCorrCNN/Training/Run_model_on_LipData.py -m ',CNNLipModelPath, ' -i ',NameStackInputDataPath,' -o ',NameStackOutputDataPath];
%[status,cmdout] = system(command,'-echo');
[status,cmdout] = system(command);

if(status>0);fprintf(cmdout);end



ReData_rf = h5read(NameStackOutputDataPath,'/LipCorrData_rf/realData' );
ImData_rf = h5read(NameStackOutputDataPath,'/LipCorrData_rf/imagData' );
ReData_rf = ReData_rf(:,(1+CNNPadding):(end-CNNPadding));
ImData_rf = ImData_rf(:,(1+CNNPadding):(end-CNNPadding));
SizeD=size(Data_frrr);
if numel(SizeD)==3;SizeD(4)=1;end
TempShort_rrrf=zeros([SizeD(2)*SizeD(3)*SizeD(4),WindowSizeOut]);
TempShort_rrrf(DataMask(:)>0,:)=ReData_rf+1j*ImData_rf;
TempShort_rrrf =reshape(TempShort_rrrf,[SizeD(2),SizeD(3),SizeD(4),WindowSizeOut]);

delete(NameStackInputDataPath);
delete(NameStackOutputDataPath);

mrsiDataLR_frrr=Data_frrr;

mrsiDataLR_frrr(MinPPM_pt_Out:MaxPPM_pt_Out,:,:,:)=permute(TempShort_rrrf,[4,1,2,3]);

MRSIDataLR_rrrf=permute(mrsiDataLR_frrr,[2,3,4,1]);


SizeD=size(MRSIDataLR_rrrf);
if ~isempty(NameData)
    fprintf([ 'Making figures...\n']);
    s1=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpecComp_AfterCNNCorr.ps'];
    s2=['./',mrsiReconParams.Log_Dir,'/',NameData,'_SVD_SpatComp_AfterCNNCorr.ps'];
    delete(s1);delete(s2);
    [Uorig,Sorig,Vorig] = svd(reshape(MRSIDataLR_rrrf.*mrsiReconParams.BrainMask,[],SizeD(4)),0);
    U_rrrc=reshape(Uorig(:,1:10), SizeD(1),SizeD(2),SizeD(3),[]);
    figs=figure('visible', 'off');
    for comp=1:10
        
        plot(MinPPM_pt_Out:MaxPPM_pt_Out,real(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)),...
            MinPPM_pt_Out:MaxPPM_pt_Out,imag(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)),...
            MinPPM_pt_Out:MaxPPM_pt_Out,abs(Vorig(MinPPM_pt_Out:MaxPPM_pt_Out,comp)));
        
        title(['Spectral comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s1);
        
        plotImage= squeeze(abs(U_rrrc(:,:,:,comp))); %Vol2Image(abs(U_rrrc(:,:,:,comp)) );
        imagesc(plotImage);
        colormap default;colorbar;
        title(['Spatial comp ',num2str(comp)]);
        print(figs, '-append', '-dpsc2', s2);
        
    end


    MRSIDataOrig_rrrf=permute(Data_frrr,[2,3,4,1]);
    MRSIDataOrig_rrrf = MRSIDataOrig_rrrf(:,:,:,MinPPM_pt_Out:MaxPPM_pt_Out);
    s=['./',mrsiReconParams.Log_Dir,'/',NameData,'_Diagonal_Spectra_BeforeAndAfterCorrection.ps'];
    delete(s);     
	for c=1:size(MRSIDataLR_rrrf,3);
		for a=1:min(size(MRSIDataLR_rrrf,1),size(MRSIDataLR_rrrf,2));if(mrsiReconParams.BrainMask(a,a,c)==1)
            
                subplot(2,2,1);title('Real part')
		sp1=squeeze(MRSIDataLR_rrrf(a,a,c,MinPPM_pt_Out:MaxPPM_pt_Out));
		sp2=squeeze(MRSIDataOrig_rrrf(a,a,c,:));

    		 plot(MinPPM_pt_Out:MaxPPM_pt_Out,real(sp1),...
           	 MinPPM_pt_Out:MaxPPM_pt_Out,real(sp2),...
            	MinPPM_pt_Out:MaxPPM_pt_Out,real(sp1-sp2));

		subplot(2,2,2);title('Imag part')
    		 plot(MinPPM_pt_Out:MaxPPM_pt_Out,imag(sp1),...
           	 MinPPM_pt_Out:MaxPPM_pt_Out,imag(sp2),...
            	MinPPM_pt_Out:MaxPPM_pt_Out,imag(sp1-sp2));

		subplot(2,2,3);title('Magnitude')
    		 plot(MinPPM_pt_Out:MaxPPM_pt_Out,abs(sp1),...
           	 MinPPM_pt_Out:MaxPPM_pt_Out,abs(sp2),...
            	MinPPM_pt_Out:MaxPPM_pt_Out,abs(sp1)-abs(sp2));
                suptitle(['x=',num2str(a) , 'y=',num2str(a) ,'z=',num2str(c)]);
    		print(figs, '-append', '-dpsc2', s);
		end;end
	end
end
fprintf([ ' Correction by CNN finished.\n']);



end

