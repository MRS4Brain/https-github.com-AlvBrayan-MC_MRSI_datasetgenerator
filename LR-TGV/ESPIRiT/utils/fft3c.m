function res = fft3c(x)

S = size(x);
fctr = S(1)*S(2)*S(3);

x = reshape(x,S(1),S(2),S(3),prod(S(4:end)));

res = zeros(size(x));
for n=1:size(x,4)
    res(:,:,:,n) = fft(fft(fft(ifftshift(ifftshift(ifftshift(x(:,:,:,n),1),2),3),[],1),[],2),[],3);
	res(:,:,:,n) = 1/sqrt(fctr)*fftshift(fftshift(fftshift(res(:,:,:,n),1),2),3);
end

res = reshape(res,S);



