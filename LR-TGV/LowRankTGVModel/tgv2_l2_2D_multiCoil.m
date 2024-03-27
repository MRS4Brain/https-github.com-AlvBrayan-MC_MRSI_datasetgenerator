function [u e] = tgv2_l2_2D_multiCoil(datackk,sense, kmask, alpha0, alpha1, maxits,Threshold)
% Primal dual TGV2 algorithm, as described in the TGV paper

% alpha0, alpha1, maxits, TGV convergeance parameters
% 
% (c) 31.8.2010
% Florian Knoll (florian.knoll@tugraz.at)
% Kristian Bredies (kristian.bredies@uni-graz.at
% Thomas Pock (pock@icg.tugraz.at)
% Rudolf Stollberger (rudolf.stollberger@tugraz.at)
%
% If you consider this code to be useful for your research, please cite:
% Knoll, F.; Bredies, K.; Pock, T.; Stollberger, R.: Second Order Total
% Generalized Variation (TGV) for MRI: Magnetic Resonance in Medicine, 
% to appear (2010)

[dxm,dym,dzm,dxp,dyp,dzp] = defineDiffOperators();
check_it =1;

alpha00 = alpha0;
alpha10 = alpha1;
%alpha01 = alpha0*reduction;
%alpha11 = alpha1*reduction;

[CN M N] = size(datackk); % numSamplesOnSpoke, numSamplesOnSpoke, nCh


u = zeros(M,N);
for c=1:CN
    %u=real(ifft(ifft(datakk.*kmask,[],1),[],2));
    u = u + conj(squeeze(sense(c,:,:))).*ifft(ifft(squeeze(datackk(c,:,:)).*kmask,[],1),[],2); % Adjoint Solution
end

%u = max(0, real(u)); %combined, positive data

v = single(zeros(size(datackk)));
p = single(zeros(M,N,2));
q = single(zeros(M,N,3));
xi = single(zeros(M,N,2));
u=single(u);
datackk=single(datackk);
sense=single(sense);
kmask=single(kmask); 
single(alpha0) ;
single(alpha1);
kmask_1kk=reshape(kmask,[1 size(kmask)]);

u_ = u;
xi_ = xi; % v in article

e = [];
 StepDiff=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiindices of the spatial derivatives

% derivatives
% | uxx uxy |
% | uyx uyy |

% multiindices
% | 1 3 |
% | 3 2 |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = sqrt(64);
tau_p = 1/16;
tau_d = 1/8;

PrevNormDiv=1E12;
PrevNormGrad=1E12;

k=-1;
RelDiffSq=1;
%for k=0:maxits
while k<100 | (RelDiffSq>Threshold & k < maxits)
k=k+1;

% update alpha's
  alpha0 = alpha0;
  alpha1 = alpha1;
 
 % alpha0 = exp(k/maxits*log(alpha01) + (maxits-k)/maxits*log(alpha00));
 % alpha1 = exp(k/maxits*log(alpha11) + (maxits-k)/maxits*log(alpha10));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SAVE VARIABLES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  uold = u; xiold = xi;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DUAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % operator

 
  %Ku0 = fft(fft(u_,[],1),[],2).*kmask;
  %for c=1:CN
   %Ku0(c,:,:) = fft(fft(squeeze(sense(c,:,:)).*u_,[],1),[],2).*kmask;
  %end
  Ku0 = fft(fft(sense.*reshape(u_,[1,size(u_)]),[],2),[],3).*kmask_1kk;
  
  Ku_ = Ku0;  
  r =  Ku_ - datackk.*kmask_1kk;% Data Fidelity gradient
  
  v = (v + tau_d*r)/(1+tau_d);%line 6 ('v' = r, 'r' = Ku - g)
  
  % gradient
  ux = dxp(u_);
  uy = dyp(u_);
  
  p(:,:,1) = p(:,:,1) - tau_d*(ux + xi_(:,:,1));
  p(:,:,2) = p(:,:,2) - tau_d*(uy + xi_(:,:,2));
  
  % projection
  absp = sqrt(abs(p(:,:,1)).^2 + abs(p(:,:,2)).^2);
  denom = max(1,absp/alpha1);
  p(:,:,1) = p(:,:,1)./denom;
  p(:,:,2) = p(:,:,2)./denom;  
  
  % symmetrized gradient
  gradxi1 = dxm(xi_(:,:,1));
  gradxi2 = dym(xi_(:,:,2));
  gradxi3 = (dym(xi_(:,:,1)) + dxm(xi_(:,:,2)))/2;
  
  q(:,:,1) = q(:,:,1) - tau_d*gradxi1; % line
  q(:,:,2) = q(:,:,2) - tau_d*gradxi2;
  q(:,:,3) = q(:,:,3) - tau_d*gradxi3;
  
 
  % projection
  absq = sqrt(abs(q(:,:,1)).^2 + abs(q(:,:,2)).^2 + 2*abs(q(:,:,3)).^2);
  denom = max(1,absq/alpha0);
  q(:,:,1) = q(:,:,1)./denom;
  q(:,:,2) = q(:,:,2)./denom;
  q(:,:,3) = q(:,:,3)./denom;  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PRIMAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  % dual operator
  ww = zeros(size(u));
 
  
  %for c=1:CN
   %ww = ww +  conj(squeeze(sense(c,:,:))).*ifft(ifft(squeeze(v(c,:,:)).*kmask,[],1),[],2);
  %end

   ww =  conj(sense).*ifft(ifft(v.*kmask_1kk,[],2),[],3);
  ww=squeeze(sum(ww,1));
  
  % divergence
  divp = dxm(p(:,:,1)) + dym(p(:,:,2));
  
  u = u - tau_p*(ww + divp); %lines 8
  
  % divergence
  divq1 = dxp(q(:,:,1)) + dyp(q(:,:,3));
  divq2 = dxp(q(:,:,3)) + dyp(q(:,:,2));
  
  xi(:,:,1) = xi(:,:,1) - tau_p*(divq1 - p(:,:,1));%line 11
  xi(:,:,2) = xi(:,:,2) - tau_p*(divq2 - p(:,:,2));%line 11
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % AUXILIARY UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  u_ = 2*u - uold;
  xi_ = 2*xi - xiold;
% tau_p = tau_p
% tau_d = tau_d
%RelDiffSq = norm(uold-u)^2/norm(uold)^2;
  if mod(k,check_it) == 0
    primal_dual_energy = 0;
    dual_energy = sum(sum(abs(u).^2))/2.0;
    e = [e dual_energy];
     RelDiffSq = norm(uold-u)^2/norm(uold)^2;
    %fprintf('TGV2-L2-2D-PD: it = %g, alpha0 = %g, energy = %g\n', k, alpha0, e);
    %imagesc(abs(u_),[0 5*mean(abs(u_(:)))]);
    %drawnow;
    %it = k
    %RelDiffSq =   RelDiffSq
   % StepDiff = [ StepDiff RelDiffSq];
   % alpha0 = alpha0 
  end
  
end
fprintf([ 'TGV Recon done in ', num2str(k),' steps. Relative Step Diff =',  num2str(RelDiffSq), '\n']);

