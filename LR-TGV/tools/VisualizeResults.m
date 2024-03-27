function VisualizeResults(reconResults,filename)
mrsiOrig=fft(fft(reconResults.mrsiData,[],3),[],2);


%mrsiOrig=reconResults.mrsiData;
[Uorig,Sorig,Vorig] = svd(reshape(mrsiOrig,size(mrsiOrig, 1),[])',0);

mrsiResults=formTensorProd(reconResults.U,reconResults.V);
%mrsiResults=fft(mrsiResults,[],4);

[Ures,Sres,Vres] = svd(reshape(mrsiResults,[],size(mrsiResults, 4)),0);

SiOri=size(mrsiOrig);
SiRes=size(mrsiResults);

reshUres=reshape(Ures, SiRes(1),SiRes(2),SiRes(3),SiRes(4));
reshUorig=reshape(Uorig, SiOri(2),SiOri(3),SiOri(1));


	
figure('visible', 'off'); 
for k=1:18
    image2plot=squeeze(sum(abs(reshUres(:,:,2,k)),3));
    subplot(6,6,k),subimage(image2plot/max(image2plot(:)));
    colormap default

    image2plot=abs(squeeze(reshUorig(:,:,k)));
    subplot(6,6,k+18),subimage(image2plot/max(image2plot(:)));
    colormap default
end;
s=sprintf('%s_spatial.jpeg',filename);
	print(s  ,'-djpeg')
	s=sprintf('%s_spatial.eps',filename);
	print(s  ,'-depsc2')
close;
	
figure('visible', 'off'); 
for k=1:18
    subplot(6,6,k),plot(real(fft(squeeze(Vres(k,:)))));
    subplot(6,6,k+18),plot(real(fft(squeeze(Vorig(k,:)))));
     %subplot(6,6,k),plot(real((squeeze(Vres(k,:)))));
   % subplot(6,6,k+18),plot(real((squeeze(Vorig(k,:)))));
end;
s=sprintf('%s_spectral.jpeg',filename);
	print(s  ,'-djpeg')
	s=sprintf('%s_spectral.eps',filename);
	print(s  ,'-depsc2')
    close;
end