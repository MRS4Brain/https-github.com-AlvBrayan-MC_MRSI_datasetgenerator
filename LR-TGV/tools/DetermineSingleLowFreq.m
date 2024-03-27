function FreqMusic=DetermineSingleLowFreq(data_trr,fMax,mrsiWaterParams)
% Les points dans la time serie sont acquis a 4000 Hz et il y a 160 pts.
% La résolution du spectre calculé avec FFT est donc de 25Hz. 
% Cependant, je voudrais calculer la position du pique en fréquence avec une meilleure précision (si possible ~ a 2-3 Hz près).
% Use EIG and or MUSIC


%% Parameters
Fe = mrsiWaterParams.Water_mrProt.samplerate; % [Hz]
p = 1; % Number of pole pairs
N = size(data_trr,1);


t = (0 : N-1) / Fe;

f = linspace (-fMax, fMax, 1024);

%% Fourier analysis
% The expected frequency (12.5 Hz) is the half of df
% df = 1 / t(end);

%% Autocorrelation matrix
% rx = xcorr (x, p, 'biased');
% rx = rx (p+1:end);
% Rx = transpose (toeplitz (rx));

%% Eigen decomposition, sorted in eigenvalue decreasing order
% [V, D] = eig (Rx);
% vp = diag(D);
% [vp, ordre] = sort (vp);
% ordre = ordre (end : - 1 : 1);
% D = D (ordre, ordre);
% V = V (:, ordre);
% Lambda = diag (D);

%% Noise variance computation
% VarNoise = Lambda (end);

%% PSD Estimate, only use the noise subspace, the one starting at 2, first is for the sine wave
% for iter = 2 : p+1
%     % [Den(:, iter-1), ~] = freqz (V(:, iter), 1, 1024);
%     Den(:, iter-1) = fft (V(:, iter), NFFT); % Do not use freqz because the signal is complex and the spectrum not even
% end;
% w = linspace (0, 2*pi, NFFT);
% Px = 1./(abs(Den).^2);
% PxSum = 1./(sum(abs(Den).^2, 2));
FreqMusic=zeros(size(data_trr,2),size(data_trr,3));
%FreqEIG=zeros(size(data_trr,2),size(data_trr,3));

for a=1:size(data_trr,2)
    parfor b=1:size(data_trr,3)
        %for b=1:size(data_trr,3)
        xn=data_trr(:,a,b);
        %% Music
        % [S, f] = pmusic (xn, 2, NFFT, Fe);
        [S, freq] = pmusic (xn, p, f, Fe);
        [Max, IxMax] = max (S);
        FreqMusic(a,b) = freq(IxMax);

        %% EIG
       % % [S1, f1] = peig (xn, 2, NFFT, Fe);
       %[S1, freq] = peig (xn, 2, f, Fe);
       % [Max1, IxMax1] = max (S1);
      %  FreqEIG(a,b) = freq(IxMax1);

    end
end


end
