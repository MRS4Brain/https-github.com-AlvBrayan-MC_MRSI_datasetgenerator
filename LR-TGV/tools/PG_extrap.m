function  Results_tkk= PG_extrap( mrsiData_tkk,mrsiReconParams,init_LThr,Vid_name)
%kmask=fftshift(fftshift(mrsiReconParams.kmask,1),2); %it must be centered at N/2
kmask=mrsiReconParams.kmask;
numSlices=size(mrsiReconParams.SkMask,3);
%SkullMask=sum(mrsiReconParams.SkMask,3);
SkullMask=imresize(sum(mrsiReconParams.SkMask,3),0.5);
SkullMask=squeeze(round(SkullMask/max(SkullMask(:))));

ImDim=size(SkullMask);

Data_fkk_orig=fft(mrsiData_tkk,[],1); 
%Data_fkk_orig=fftshift(fftshift(Data_fkk_orig,2),3); %it must be centered at N/2
NbTp=size(Data_fkk_orig,1);
diff_size=(ImDim-size(kmask));
Data_fkk_orig_zp=fftshift(fftshift(padarray(fftshift(fftshift(Data_fkk_orig,2),3),[0 ,diff_size(1)/2,diff_size(2)/2 ]),2),3);
kmask_zp=fftshift(fftshift(padarray(fftshift(fftshift(kmask,1),2),[diff_size(1)/2,diff_size(2)/2 ]),1),2);



Freq_kmask_zp=zeros(NbTp,size(kmask_zp,1),size(kmask_zp,2));
Freq_kmask=zeros(NbTp,size(kmask,1),size(kmask,2));
for(a=1:NbTp)
    Freq_kmask_zp(a,:,:)=kmask_zp(:,:);
    Freq_kmask(a,:,:)=kmask(:,:);
 end
Freq_kmask_zp=logical(Freq_kmask_zp);
Freq_kmask=logical(Freq_kmask);

SR=mrsiReconParams.mrProt.samplerate;
HzpP=SR/NbTp;


Data_trr=ifft(ifft(ifft(Data_fkk_orig_zp,[],1),[],2),[],3);
%for l = 1 : NbTp
%   Data_trr(l,:,:)=exp(-2*pi*1i*(l-1) * b0map / SR).*squeeze(Data_trr(l,:,:));
%end
Data_frr=fft(Data_trr,[],1);




low_bnd_L=round(200/HzpP);
high_bnd_L=round(500/HzpP);
low_bnd_M=round(100/HzpP);
high_bnd_M=round(200/HzpP);
LThr=init_LThr;
mu=0.95;%0.95;

plot_bound=max( abs(Data_frr(round(100/HzpP):round(400/HzpP),round(size(Data_frr,2)/3),round(size(Data_frr,3)/2))) );

Energy(1)=0;
Energy(2)=1;
n=2;

while(abs(Energy(n)-Energy(n-1))/Energy(n)>1E-9)
 
        
     Metabo_Vol=squeeze(sum(abs(Data_frr(low_bnd_M:high_bnd_M,:,:)),1));
  Lipid_Vol=squeeze(sum(abs(Data_frr(low_bnd_L:high_bnd_L,:,:)),1));
 
  LipMask=Lipid_Vol/max(Lipid_Vol(:));
  LipMask=((LipMask.*SkullMask)>LThr);
 
%{
  Freq_LipMask=zeros(NbTp,size(LipMask,1),size(LipMask,2));
  for(a=1:size(Data_frr,1))
      Freq_LipMask(a,:,:)=LipMask;
  end
  Freq_LipMask=logical( Freq_LipMask);
  %}
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Make video frames
 if(~isempty(Vid_name))
  % [pathstr,name,ext] = fileparts([ filename, '.test'] ) ;
     ImSiC=2048;ImSiR=2048;
    frame=figure('Visible','off','Position', [100, 100, ImSiC,ImSiR]);
  
    image2plot=Lipid_Vol;
    subplot(2,2,1),imagesc(image2plot/(max(image2plot(:))));  
    image2plot=LipMask;
    subplot(2,2,2),imagesc(image2plot);  
    image2plot=Metabo_Vol.*abs(LipMask-1);
    subplot(2,2,3),imagesc(image2plot); 
    subplot(2,2,4),plot(HzpP*(round(100/HzpP):round(400/HzpP)),abs(Data_frr(round(100/HzpP):round(400/HzpP),round(size(Data_frr,2)*0.4),round(size(Data_frr,3)*0.25)))); 
    axis([-inf inf 0 plot_bound  ])
    %axis('off');
    colormap default 
    saveas(frame,sprintf('%s/Lip_fr%g.png','.',n-1),'png');

    close all;
  end
  
  
  
  Data_frr_masked=zeros(size(Data_frr)); 
  
  for(a=1:size(Data_frr,2))
      for(b=1:size(Data_frr,3))
          if(LipMask(a,b)==1)
            Data_frr_masked(:,a,b)=Data_frr(:,a,b);
          end 
      end
  end
  
%   Data_frr_masked(Freq_LipMask)=Data_frr(Freq_LipMask);
   
 Data_tkk=ifft(fft(fft(Data_frr_masked,[],2),[],3),[],1);
%for l = 1 : NbTp
%   Data_tkk(l,:,:)=exp(2*pi*1i*(l-1) * b0map / SR).*squeeze(Data_tkk(l,:,:));
%end
Data_fkk=fft(Data_tkk,[],1);

 for(a=1:size(Data_frr,2))
      for(b=1:size(Data_frr,3))
          if(kmask_zp(a,b)==1)
            Data_fkk(:,a,b)=Data_fkk_orig_zp(:,a,b);
          end 
      end
  end
%Data_fkk(Freq_kmask)=Data_fkk_orig_zp(Freq_kmask_zp);

 
 Data_trr=ifft(ifft(ifft(Data_fkk,[],1),[],2),[],3); 
 %for l = 1 : NbTp
 %  Data_trr(l,:,:)=exp(-2*pi*1i*(l-1) * b0map / SR).*squeeze(Data_trr(l,:,:));
%end
Data_frr=fft(Data_trr,[],1);
 
 
 
 LThrmask=mu*min(Lipid_Vol(LipMask))/max(Lipid_Vol(LipMask));
    if(LThr<LThrmask & mod(n,1)==0)
        LThr=LThrmask;
    end
    
 n=n+1;
 temp=squeeze(sum(abs(Data_frr).^2,1)).*abs(LipMask-1);
 Energy(n)=sum(temp(:))/numel(temp(:));
% abs(Energy(n)-Energy(n-1))/Energy(n);
 Lipid_Energy= sum(sum((Lipid_Vol.*abs(LipMask-1)).^2,1),2) /sum(sum(abs(LipMask-1),1),2) ;
 
 LThr;
  
end
%Results_rrrt=permute(repmat( ifft(Data_frr,[],1),[1,1, 1, numSlices]),[2 3 4 1]);
Results_tkk=fft(fft(ifft(Data_frr,[],1),[],2),[],3);

fprintf('Lipid Energy: %g \n',Lipid_Energy)
fprintf('Lipid Threshold: %g \n',LThr)

if(~isempty(Vid_name))
    Fig2Video( './Lip_fr',n-2,['PGLipidSupp_',Vid_name] )
    delete(  './Lip_fr*.png');
end
end
