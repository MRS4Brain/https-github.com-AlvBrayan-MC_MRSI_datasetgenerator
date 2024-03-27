function [filtMRSI , WaterMRSI, WaterFreq, WaterAmp ]=FilterMRSIData(MRSIData_tkk,mrProt,FiltParam,kmask)

HzpP = mrProt.samplerate/size(MRSIData_tkk,1);
%load('./CombinedDat-CSI.m','-mat');

NDim=sum(size(MRSIData_tkk)>1)-1;


if(NDim==2)
    step=0;
    tot_steps= sum(kmask(:));%(size(MRSIData_tkk,2)*size(MRSIData_tkk,3));
    for a=1:size(MRSIData_tkk,2);
        for b=1:size(MRSIData_tkk,3);
            if(kmask(a,b))
                step=step+1;
                if(mod(step,round(tot_steps/5))==0)
                    fprintf(sprintf('filter step %g out of %g. \n',step,tot_steps));
                end
                [filtMRSI(:,a,b), WaterMRSI(:,a,b), WaterF1 , WaterA1] = mrsiHSVDFilter(MRSIData_tkk(:,a,b),mrProt.NMRFreq,1.0/mrProt.DwellTime,-4.7, FiltParam.Comp,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
               %Time reversed filtering
               [TempfiltMRSI, TempWaterMRSI, WaterF2  , WaterA2] = mrsiHSVDFilter(flip(filtMRSI(:,a,b)),mrProt.NMRFreq,1.0/mrProt.DwellTime,-4.7, FiltParam.Comp,-FiltParam.Water_maxFreq,-FiltParam.Water_minFreq);
                 %filtMRSI(:,a,b,c) = flip( mrsiHSVDFilter(flip(filtMRSI(:,a,b,c)),mrProt.NMRFreq,1.0/mrProt.DwellTime,-4.7, FiltParam.Comp,-FiltParam.Water_maxFreq,-FiltParam.Water_minFreq));
                  filtMRSI(:,a,b) = flip(TempfiltMRSI);
                 WaterMRSI(:,a,b) = WaterMRSI(:,a,b) + flip(TempWaterMRSI);
                 WaterAmp(a,b) = sum(abs(fftshift(fft(WaterMRSI(:,a,b)))),1);
                 WaterFreq(a,b) = HzpP* (sum(abs(fftshift(fft(WaterMRSI(:,a,b)))).*(1:size(WaterMRSI,1))',1)/sum(abs(fftshift(fft(WaterMRSI(:,a,b)))),1) - size(WaterMRSI,1)/2);

               if FiltParam.FilterFreq>0
                    filtMRSI(:,a,b)= mrsiExpFilter(filtMRSI(:,a,b),1.0/mrProt.DwellTime,FiltParam.FilterFreq);
               end

               %  plot(1:numel(squeeze(MRSIData_tkk(:,a,b))),real(fft(squeeze(filtMRSI(:,a,b)),[],1)),1:numel(squeeze(MRSIData_tkk(:,a,b))),real(fft(MRSIData_tkk(:,a,b),[],1)));
              % drawnow;
                 % figure();
                % plot(1:numel(reflect_coeffs),reflect_coeffs,1:numel(reflect_coeffs),cumsum(reflect_coeffs)./(1:numel(reflect_coeffs))' );
            end
        end;
    end
elseif(NDim==3)
    step=0;
    tot_steps=(size(MRSIData_tkk,2)*size(MRSIData_tkk,3)*size(MRSIData_tkk,4));
    for a=1:size(MRSIData_tkk,2);for b=1:size(MRSIData_tkk,3);for c=1:size(MRSIData_tkk,4);
        step=step+1;
        if(mod(step,round(tot_steps/100))==0)
            fprintf(sprintf('filter step %g out of %g. \n',step,tot_steps));
        end          
         [filtMRSI(:,a,b), WaterMRSI(:,a,b) ] = mrsiHSVDFilter(MRSIData_tkk(:,a,b,c),mrProt.NMRFreq,1.0/mrProt.DwellTime,-4.7,FiltParam.Comp,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
        %Time reversed filtering
        [TempfiltMRSI, TempWaterMRSI ] = mrsiHSVDFilter(flip(filtMRSI(:,a,b,c)),mrProt.NMRFreq,1.0/mrProt.DwellTime,-4.7, FiltParam.Comp,-FiltParam.Water_maxFreq,-FiltParam.Water_minFreq);
       filtMRSI(:,a,b,c) = flip(TempfiltMRSI);
       WaterMRSI(:,a,b,c) = WaterMRSI(:,a,b,c) + flip(TempWaterMRSI);
        if FiltParam.FilterFreq>0
            filtMRSI(:,a,b,c) = mrsiExpFilter(filtMRSI(:,a,b,c),1.0/mrProt.DwellTime,FiltParam.FilterFreq);
        end
       
    end;end;end
else
    error('the dimensionality of MRSIData_tkk in FilterMRSIData is different than 2 or 3.');
end

end