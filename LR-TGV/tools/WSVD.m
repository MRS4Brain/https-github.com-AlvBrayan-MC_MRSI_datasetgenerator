function [MRSIdata_trr CCoef Q]=WSVD(MRSIdata_ctrr, Range_noise);
SData=size(MRSIdata_ctrr);
MRSIdata_trr=zeros(SData(2:4));
for a=1:SData(3);for b=1:SData(4)
        Noise_tc=transpose(MRSIdata_ctrr(:,Range_noise(1):Range_noise(2),a,b));
        Noise_tc=detrend(Noise_tc,'constant');% remove the mean (normally zero but in case of bad baseline)
        
        for c=1:size(Noise_tc,2)for d=1:size(Noise_tc,2)
                Noise_corr(c,d)=mean(squeeze(Noise_tc(:,c).*conj(Noise_tc(:,d))));
            end;end
        
        [X,D] = eig(Noise_corr);% following the notation of Rodgers et al. 2016
        W=D^(-0.5)*X';
        if(isnan(sum(W(:))))
            W=eye(size(W));
        end
        S=transpose(W*MRSIdata_ctrr(:,:,a,b));
        
        [V,Sig,U]=svd(S);
        
        Q=Sig(1,1)*V(:,1);
        
        CCoef=inv(W)*U(:,1);
        
        for c=1:size(CCoef,1)
            MRSIdata_trr(:,a,b)=MRSIdata_trr(:,a,b)+transpose(CCoef(c)*MRSIdata_ctrr(c,:,a,b));
        end
        
    end;end;
