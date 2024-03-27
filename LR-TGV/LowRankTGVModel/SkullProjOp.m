function out = SkullProjOp(Data_rrc,SkullMask)
imDims=size(x);
nDims = ndims(x);
out=zeros(imDims);
spa_index_im = repmat({':'},[ 1, nDims-1]);

for comp=1:imDims(end)
    out(spa_index_im{:},comp) = SkullMask.*x(spa_index_im{:},comp);
end
end