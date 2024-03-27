function [Lipids ] = make_QR_LipidBasis_inK( data_rrf, lipid_mask, Nbasis )
%GET_LIPIDBASIS Summary of this function goes here
%   Detailed explanation goes here


count = 0;
N = size(data_rrf);
Lipid_Basis = zeros(N(3), sum(lipid_mask(:)));


%{
Lipid_img=data_rrf.*repmat(lipid_mask,[1 1 N(end)]);
Lipid_img=fft(fft(Lipid_img,[],1),[],2);%in K space
Lipid_img=reshape(Lipid_img,[],N(end));

LipidVol=sum(abs(Lipid_img).^2,2);
[~,Ind]=sort(LipidVol(:));
Lipid_img_sorted=Lipid_img(flip(Ind),:);
[Q,R]=qr(Lipid_img_sorted(1:Nbasis,:));
%}

% SVD basis of the skull lipids 
Lipid_img=data_rrf.*repmat(lipid_mask,[1 1 N(end)]);
Lipid_img=fft(fft(Lipid_img,[],1),[],2);%in K space
[~,~,Vorig] = svd(reshape(Lipid_img,[],N(end)),0);
Lipids=conj(Vorig(:,1:Nbasis));

%{
[Q,R]=qr(Lipids'); %for SVD basis of the skull lipids 

Lipids=R';

Lipids=Lipids*inv(diag(diag(Lipids'*Lipids)));
%}
%Lipids=transpose(Lipid_img_sorted(1:Nbasis,:));
end

