function [ c ] = circular_mask( mesh_size, rad )
%CIRCULAR_MASK Summary of this function goes here
%   Detailed explanation goes here

[x,y] = meshgrid(1:mesh_size(1),1:mesh_size(2));


c = (x-1-mesh_size(1)/2).^2 + (y-1-mesh_size(2)/2).^2 <= rad^2;

end

