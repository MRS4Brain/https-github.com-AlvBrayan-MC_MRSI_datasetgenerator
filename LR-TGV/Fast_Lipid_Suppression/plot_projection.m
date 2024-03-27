function [  ] = plot_projection( img, fig_no, axis_scale )
%PLOT_PROJECTION Summary of this function goes here
%   Detailed explanation goes here


proj = sum(abs(img),3);

if nargin < 2
    figure(1), imagesc( proj ) , axis image off, drawnow
else
    figure(fig_no), imagesc( proj ) , axis image off, drawnow
end


if nargin == 3
    caxis(axis_scale)
end


end

