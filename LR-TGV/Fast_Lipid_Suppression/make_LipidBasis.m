function [ Lipid_Basis ] = make_LipidBasis( csi, lipid_mask )
%GET_LIPIDBASIS Summary of this function goes here
%   Detailed explanation goes here


count = 0;
N = size(csi);
Lipid_Basis = zeros(N(3), sum(lipid_mask(:)));

for ay = 1:N(1)
    for cey = 1:N(2)
        
        if lipid_mask(ay,cey) ~= 0
            % invivo lipid voxel
            count = count + 1;
            
            Lipid_Basis(:,count) = csi(ay,cey,:);
        end
        
    end
end




end

