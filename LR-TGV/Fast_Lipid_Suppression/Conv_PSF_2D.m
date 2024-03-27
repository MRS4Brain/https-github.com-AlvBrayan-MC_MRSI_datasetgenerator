function SkMask2D_ktrunc=Conv_PSF_2D(SkMask,kmask)
SkMask2D=sum(SkMask,3)/size(SkMask,3);

OverResFact=10;
Filt_sig=1.5;
SkMask2D=imresize(SkMask2D,OverResFact);

ImSi=size(SkMask2D);
SpecSi=size(kmask);


%PSF=@(x)sin(pi*x/OverResFact)./(pi*x/OverResFact);
%FiltKern2=@(x)abs(sin(pi*x/OverResFact)./(pi*x/OverResFact));

FiltKern=@(x)normpdf(x,0,Filt_sig) ;
x_axis=(1:SpecSi(1))-round(SpecSi(1)/2);
y_axis=(1:SpecSi(2))-round(SpecSi(2)/2);

[XFilt,YFilt] = meshgrid(FiltKern(x_axis),FiltKern(y_axis));
Filt2d=XFilt.*YFilt;
Filt2d=fftshift(Filt2d/max(Filt2d(:)));

%{
[XFilt2,YFilt2] = meshgrid(FiltKern2(x_axis),FiltKern2(y_axis));
Filt2d2=XFilt2.*YFilt2;
Filt2d2=fftshift(Filt2d2/max(Filt2d2(:)));
%}

DiffSi=ImSi-SpecSi;
%kmask_zp=fftshift(fftshift(padarray(fftshift(fftshift(kmask,1),2),[ round(DiffSi(1)/2),round(DiffSi(2)/2) ]),1),2);
Center_k_SkMask=fftshift(fftshift(fft2(SkMask2D),1),2);
LowResk_SkMask=fftshift(fftshift(Center_k_SkMask((1+round(DiffSi(1)/2)):(end-round(DiffSi(1)/2)),(1+round(DiffSi(2)/2)):(end-round(DiffSi(2)/2))),1),2);
%SkMask2D_ktrunc=ifft2(kmask_zp.*fft2(SkMask2D));

SkMask2D_ktrunc=abs(ifft2(kmask.*LowResk_SkMask));
SkMask2D_ktrunc=ifft2(fft2(SkMask2D_ktrunc).*fft2(Filt2d));
SkMask2D_ktrunc=SkMask2D_ktrunc/max(SkMask2D_ktrunc(:));

end