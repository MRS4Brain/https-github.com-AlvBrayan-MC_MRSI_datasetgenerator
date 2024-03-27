function [FreqShift,MaxCCoef] = MeasureFreqShift( TS1,TS2, Time,FreqRange,FreqPrec )
ind=0;
ftTS1=fft(TS1);
ftTS2=fft(TS2);
range=1:numel(TS2);%20:110;
PNSq=(norm(ftTS1(range))*norm(ftTS2(range)))^2;
for FreqS=-FreqRange:FreqPrec:FreqRange
    ind=ind+1;
    %CM=corrcoef(TS1,TS2.*exp(2*pi*1i*Time*FreqS));
    %CorCoef(ind)=abs(CM(1,2))^2;
   ftTS2=fft(TS2.*exp(2*pi*1i*Time*FreqS));
    CorCoef(ind)=abs(ftTS1(range)'*ftTS2(range))^2/PNSq;
    FreqInd(ind)=FreqS;
end
%[MaxCCoef, I]=max(CorCoef);
%FreqShift=FreqInd(I);

ThresCC=max(CorCoef(:))/2;
pt=(CorCoef>ThresCC);
 FreqShift=(sum(CorCoef(pt).*(FreqInd(pt)))/sum(CorCoef(pt)));
MaxCCoef=max(CorCoef(:));

% plot(FreqInd,CorCoef)
%plot(1:numel(TS1),real(fft(TS1)),1:numel(TS1),norm(TS1)/norm(TS2)*real(fft(TS2.*exp(2*pi*1i*Time*FreqShift))))
 
end

