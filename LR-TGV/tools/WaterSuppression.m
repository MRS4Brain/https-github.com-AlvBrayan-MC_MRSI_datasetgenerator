function [filtMRSI ]=WaterSuppression(MRSIData_tkk,mrProt,FiltParam,kmask)


filtMRSI=zeros(size(MRSIData_tkk));
step=0;
MRSIData_trr=ifft(ifft(MRSIData_tkk,[],2),[],3);
tot_steps= sum(kmask(:));
for b=1:size(MRSIData_tkk,3);
    for a=1:size(MRSIData_tkk,2);
      %  if(kmask(a,b))
            step=step+1;
%             if(mod(step,round(tot_steps/5))==0)
%                 fprintf(sprintf('filter step %g out of %g. \n',step,tot_steps));
%             end
            MRSIData_trr(:,a,b) =  Fast_HSVD_Filter(MRSIData_trr(:,a,b),1.0/mrProt.DwellTime,FiltParam.Comp,FiltParam.Water_minFreq,FiltParam.Water_maxFreq);
      %  end
    end
end
filtMRSI=fft(fft(MRSIData_trr,[],2),[],3);
end