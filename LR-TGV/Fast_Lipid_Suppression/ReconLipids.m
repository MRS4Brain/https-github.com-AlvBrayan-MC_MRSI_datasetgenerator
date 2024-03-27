function [Lipids_ctkk,Lipids_tkk]  = ReconLipids(  mrsiReconParams )
 % mrsiReconParams.mrsiData dims: time-k-k


OrderLip=mrsiReconParams.OrderLipRecon;%256;%128;
N = size(mrsiReconParams.mrsiData);
IOP=diag(ones(N(end),1));


%Determine Lipids_rrf
[Uorig,Sorig,Vorig] = svd(reshape(permute(mrsiReconParams.mrsiData,[1 3 4 2]),[],N(2)),0);
V_tc=Vorig(:,1:OrderLip);
S=Sorig(1:OrderLip,1:OrderLip);
U_ckkc=reshape(Uorig(:,1:OrderLip), N(1),N(3),N(4),[]);

kmask=mrsiReconParams.kmask;

%alpha = mean(kmask(:))*mrsiReconParams.mu_tv/reduction; %  % correted for undersampling
%reduction = 2^(-8);     % usually there is no need to change this
alpha = 1E-8;%/reduction; %  % correted for undersampling
maxit = 1000 ;          % use 1000 Iterations for optimal image quality        % usually there is no need to change this
Threshold=1E-6;
U_rrc=zeros(N(3),N(4),OrderLip);
fprintf('Start Lipid Reconstruction...\n');
parfor k=1:OrderLip
    % fprintf([ 'Processing Component ', num2str(k), ' ...\n']);
    [U_rrc(:,:,k) e] = tgv2_l2_2D_multiCoil(U_ckkc(:,:,:,k),mrsiReconParams.SENSE, kmask, 2*alpha, alpha, maxit,Threshold);
end

Lipids_rrt=formTensorProduct(U_rrc, V_tc*S,2);

N = size(Lipids_rrt);
HzpP=mrsiReconParams.mrProt.samplerate/N(end);
 
Temp_trr=zeros(N(3),N(1),N(2));
Lipids_tkk = fft(fft(permute(Lipids_rrt,[3,1,2]),[],2),[],3);
for c=1:size(mrsiReconParams.mrsiData,1)
	Temp_trr=mrsiReconParams.SENSE(c,:,:).*permute(Lipids_rrt,[3,1,2]);
    Lipids_ctkk(c,:,:,:)=fft(fft(Temp_trr,[],2),[],3);
end


end

