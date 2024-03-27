function [ Reph_CData_f ] = Rephase_poly(  CData_f)

    N=numel(CData_f);
   % Nmax=sum(mrsiReconParams.ppm<0);
  %  Ph=[0;0;0];
   %  Ph=[0;0;0];
    step=1;  
    Reph_CData_f=CData_f;
    
  % CData_f=fft( Reph_CData_t);

    CorrSpectr=@(Ph)1.0/EvalExcursRatio(CData_f, Ph,0.10); 
 
   %f=fit((1:Nmax)',unwrap(angle(CData_f(1:Nmax))),'poly1')
   
  optnew = optimset('MaxFunEvals',1e4,'MaxIter',1e4);
 
  %optnew = optimset('Display','iter','PlotFcns',@optimplotfval);
  
    Ph=fminsearch(CorrSpectr,[0,0],optnew );
   % Ph=fminsearch(CorrSpectr,[0,0,0],optnew );
    
  % options.Algorithm = 'levenberg-marquardt';
  % Ph=lsqnonlin(CorrSpectr,[0,0,0] ,[],[], options);

   % Reph_CData_t=ifft(fft(Reph_CData_t).*exp(-1j*(Ph(3)+(1:N)'*Ph(2)+(1:N).^2'*Ph(1))));
   %  Reph_CData_f=Reph_CData_f.*exp(-1j*(Ph(3)+(1:N)'*Ph(2)+(1:N).^2'*Ph(1)));
     Reph_CData_f=Reph_CData_f.*exp(-1j*(Ph(2)+(1:N)'*Ph(1)));
    step=step+1;

 %  figure();
  % ppm=(-4.7+((1:N)*4000/(N*3*42.58)));
 %   plot(ppm,real(fft(CData_t)),ppm,real(fft(Reph_CData_t)))
end

