function [ CData_KKCorr ] = KKRecurs( CData_t)

%CData_t=fftshift(CData_t);
orig_CData_f=fft(CData_t);
%Re_CData_f=real(fft(CData_t));
%Im_CData_f=imag(fft(CData_t));
%CData_f=abs(fft(CData_t));

CData_f=fft(CData_t);
Norm_CData_f=abs(fft(CData_t));
N=numel(CData_f);
Im_it_diff=1;

iter=1;
%figure();
%diffiter(1)=0;
while Im_it_diff>5e-4 & iter< 1E3%Im_it_diff>1e-4  5e-4 & iter< 1E3
    if(iter<500)
		p=0.5;
	elseif(iter<1000)
		p=0.5;
	else
		p=0.2;	
    end
    if(mod(iter,2)==1)
     Reph_CData_t=ifft(Rephase_poly(CData_f));
     else
     Reph_CData_t= ifft(Rephase_poly_minImag(CData_f)); 
    end
    
    Reph_CData_t=ifft(CData_f);
    
    %Re_CData_t_KK=[real(Reph_CData_t(1:floor(N*0.5)) ) -real(Reph_CData_t(ceil(N*0.5):end ))];%the data go crescendo
	%Re_CData_t_KK=[imag(Reph_CData_t(1:floor(N*0.5)) ) -imag(Reph_CData_t(ceil(N*0.5):end ))];
     New_CData_f=fft([Reph_CData_t(1:round(N*0.5));  -Reph_CData_t(round(N*0.5+1):end)]);
 
  % New_CData_f=fft(Reph_CData_t);

    %New_CData_f=-1i*imag(hilbert(real(Reph_CData_f)));
    %New_CData_f=New_CData_f+imag(hilbert(imag(Reph_CData_f)));
    

    New_CData_f=Norm_CData_f./abs(New_CData_f).*New_CData_f;
 
		
	%plot([1:250],real(New_CData_f(1:250)),[1:250],imag(New_CData_f(1:250)),[1:250],real(orig_CData_f(1:250)));
   % plot([1:N],real(New_CData_f(1:N)),[1:N],imag(New_CData_f(1:N)),[1:N],abs(Norm_CData_f(1:N)));
	%drawnow;
     
    New_CData_f=(New_CData_f*p+CData_f*(1-p));

    Im_it_diff=sum(abs(New_CData_f(:)-CData_f(:)))/sum(abs(CData_f(:)));
 %   diffiter(end+1)=Im_it_diff;
    CData_f=New_CData_f;
   % Re_CData_f=real(CData_f);
  %  Im_CData_f=imag(CData_f);
	iter=iter+1;
    if(mod(iter,100)==0)
  %      iter
    %    Im_it_diff
    end

%drawnow;
end
%{
figure
Aorig_CData_f=abs(orig_CData_f) .* exp (-1j * imag (hilbert (log(abs(orig_CData_f))))); 
plot([1:N],real(New_CData_f(1:N)),[1:N],real(Aorig_CData_f(1:N)),[1:N],real(orig_CData_f(1:N)));
drawnow;
%}
%figure
%semilogy(diffiter);
%drawnow;
iter=iter
%CData_KKCorr=fftshift(ifft(CData_f));
CData_KKCorr=ifft(CData_f);



end

