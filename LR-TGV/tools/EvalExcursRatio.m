function [ Ratio ] = EvalExcursRatio(CData_f, Ph,Ratio_pt)
%EVALEXCURSRATIO Summary of this function goes here
%   Detailed explanation goes here
N=numel(CData_f);
MeasN=round(N*Ratio_pt);
R_Corr_CData_f=real(CData_f.*exp(-1j*(Ph(2)+(1:N)'*Ph(1))));
%R_Corr_CData_f=real(CData_f.*exp(-1j*(Ph(3)+(1:N)'*Ph(2)+(1:N).^2'*Ph(1))));
[SortR_Corr_CData_f,I]=sort(R_Corr_CData_f);
Ratio=sum(abs(SortR_Corr_CData_f((end-MeasN+1):end)).^2)/sum(abs(SortR_Corr_CData_f(1:MeasN)).^2);

%Ratio=sum(imag(R_Corr_CData_f));
%ppm=(-4.7+((1:N)*4000/(N*3*42.58)));

%plot(ppm,real(CData_f), ppm(I((end-MeasN+1):end)),real(SortR_Corr_CData_f((end-MeasN+1):end)),'*',ppm(I(1:MeasN)),real(SortR_Corr_CData_f(1:MeasN)),'.')

end

