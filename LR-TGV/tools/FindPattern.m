function [ position ] = FindPattern(pat,mask)
%pat is the pattern to be found in the binary matrix 2D mask
%It is a vector or a submatrix made of 1 or 0

Spat=size(pat);
Smask=size(mask);
position=[];
for a=0:(Smask(1)-Spat(1))
    for b=0:(Smask(2)-Spat(2))
        match=(pat==mask(a+[1:Spat(1)],b+[1:Spat(2)]));        
        if prod(match(pat>0))
            position(end+1,1)=a+1;
            position(end,2)=b+1;
        end
    end
end
end

