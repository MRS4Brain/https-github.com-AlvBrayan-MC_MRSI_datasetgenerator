function [frequencies, dampings, basis, ahat] = HSVD(y, fs, K)
%
% Decompose the signal y using the method of Barkhuijsen, et al. (1987)
%
% Obligatory arguments:
% 'y' is the FID (a linear combination of complex exponentials).
% 'fs' is the sampling frequency (bandwidth).
% 'K' is the desired model order (expected number of complex
% exponentials comprising the signal)
%
% Outputs:
% 'frequencies' - frequencies of components in signal
% 'dampings' - damping of components in signal
% 'basis' - basis vectors, one for each component
% 'ahat' - amplitudes of each basis in signal
%
% Author: Greg Reynolds (remove.this.gmr001@bham.ac.uk)
% Date: August 2006 

N = length(y);
L = floor(0.5*N);
M = N+1-L;

% H is the LxM Hankel LP matrix
H       = hankel(y(1:L), y(L:N));

% compute the singular value decomposition
[U,~,~] = svd(H);

% construct H of rank K
Uk = U(:,1:K);

% find Ukt and Ukb
Ukt = Uk(2:end,:);
Ukb = Uk(1:end-1,:);
Zp  = pinv(Ukb)*Ukt;

% find the poles
[Q,~]       = eig(Zp);
Z           = inv(Q)*Zp*Q;
q           = log(diag(Z));
dt          = 1/fs;
dampings    = real(q)/dt;
frequencies = imag(q)/(2*pi)/dt;

%dampings(dampings>10)
dampings(dampings>10)    = 10;


% construct the basis
t     = (0:dt:(length(y)-1)*dt);
basis = exp(t.'*(dampings.' + 2*pi*1i*frequencies.'));


% compute the amplitude estimates
ahat = pinv(basis(1:length(y),:))*y;

end
