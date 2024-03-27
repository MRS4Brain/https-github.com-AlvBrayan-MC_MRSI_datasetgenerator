function Xnew = formTensorProd(X,Y)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%X spatial x nb components
% Y spectral x nb components

nDimsX = ndims(X);
sizeX  = size(X);

if (nDimsX < 4)
    sizeX = [sizeX, ones(1, 4 - nDimsX)];
end

Xnew = reshape(reshape(X, [], size(Y, 2)) * Y', [sizeX(1:(end-1)), size(Y,1)]);

end




