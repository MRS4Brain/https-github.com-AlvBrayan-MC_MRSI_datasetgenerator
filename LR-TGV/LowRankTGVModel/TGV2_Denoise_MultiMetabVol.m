function [imOut, funCost] = TGV2_Denoise_MultiMetabVol(im,mu_tv)
% tgv2_denoise_multivol - 
% Solves the total generalized variation denoising problem for a collection
% of 3D volumes. For details see:
% Bredies, K., Kunisch, K., and Pock, T. "Total Generlized Variation" SIAM
% Journal on Imaging Sciences 3(3), 2010, pp. 492-526
% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUTS:
%       im: The original, noisy images [(Nd) x numVols], where (Nd) denotes 
%           the dimensions of each image volume
%   alpha0: Regularization parameter controlling the weighting of first
%           order information
%   alpha1: Regularization parameter controlling the weighting of second
%           order information
%   params: A paramter structure containing any of the following fields:
%           
%           init: An initial guess for the solution [(Nd) x numVols]
%        verbose: Determines whether verbose output is desired (logical)
%      reduction: Constant reduction factor for alpha0 and alpha1 over the
%                 course of the iterations (scalar double)
%        maxIter: The maximum number of allowed iterations (scalar)
%       stepSize: A length Nd vector containing step sizes to use in each
%                 dimension for finite difference calculations
%--------------------------------------------------------------------------
%   OUTPUTS:
%       imOut: The denoisied images
%     funCost: A vector containing the values of the cost functional at
%              each iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the input dimensions
[nRows, nCols, nSlices, nVols] = size(im);
imDims = [nRows, nCols, nSlices, nVols];

        
maxIter   = 1000;
minIter   = 50;
u         = zeros(imDims, class(im));
verbose   = 1;
stepSize  = [1, 1, 1];
reduction = 1;
alpha0    = 2 ;
alpha1    = 1 ;
ConvThres = 1.0E-5;
  
% Initialize all other variables 

alpha00   = alpha0*mu_tv;
alpha10   = alpha1*mu_tv;
alpha01   = alpha0*mu_tv * reduction;
alpha11   = alpha1*mu_tv * reduction;

p         = zeros([imDims, 3], class(im));
q         = zeros([imDims, 6], class(im));
xi        = zeros([imDims, 3], class(im));
u_update  = u;
xi_update = xi;
colIndcs  = repmat({':'}, [1, numel(imDims)]); 
step_p    = 1 / sqrt(12);   % Step size for primal update
step_d    = 1 / sqrt(12);   % Step size for dual update

funCost   = [];

% Define the differential operators 
%--------------------------------------------------------------------------
% Indexing scheme:
%   Derivatives      Multi-indeces
% | uxx uxy uxz |     | 1 4 5 |
% | uyx uyy uyz |     | 4 2 6 |
% | uzx uzy uzz |     | 5 6 3 |
%--------------------------------------------------------------------------
[fhdxm, fhdym, fhdzm, fhdxp, fhdyp, fhdzp] = defineDiffOperators(stepSize);


% Begin the iterations
iter=0;
CostConv=1;
while (iter <= minIter || (iter <= maxIter && CostConv>ConvThres ))
%for iter = 0 : maxIter
    
    alpha0 = exp((iter / maxIter) * log(alpha01) + ...
              ((maxIter - iter) / maxIter) * log(alpha00));
          
    alpha1 = exp((iter / maxIter) * log(alpha11) + ...
              ((maxIter - iter) / maxIter) * log(alpha10));
          
    %----------------------------------------------------------------------
    % Store the results of the previous iteration
    % ---------------------------------------------------------------------
    
    u_old  = u;
    xi_old = xi;
    
    %----------------------------------------------------------------------
    % Dual Update
    % ---------------------------------------------------------------------
    
    % Compute the gradient
    du_x = fhdxp(u_update);
    du_y = fhdyp(u_update);
    du_z = fhdzp(u_update);
    
    % Take a step 
    p(colIndcs{:}, 1) = p(colIndcs{:}, 1) - step_d * (du_x + xi_update(colIndcs{:}, 1));
    p(colIndcs{:}, 2) = p(colIndcs{:}, 2) - step_d * (du_y + xi_update(colIndcs{:}, 2));
    p(colIndcs{:}, 3) = p(colIndcs{:}, 3) - step_d * (du_z + xi_update(colIndcs{:}, 3));
    
    % Compute the projection
    p_norm            = sqrt(abs(p(colIndcs{:}, 1)).^2 + ...
                             abs(p(colIndcs{:}, 2)).^2 + ...
                             abs(p(colIndcs{:}, 3)).^2);
              
    denom             = max(1, (p_norm / alpha1));
    p(colIndcs{:}, 1) = p(colIndcs{:}, 1) ./ denom;
    p(colIndcs{:}, 2) = p(colIndcs{:}, 2) ./ denom;
    p(colIndcs{:}, 3) = p(colIndcs{:}, 3) ./ denom;
    
    % Compute the symmetrized gradient
    dxi_1 = fhdxm(xi_update(colIndcs{:}, 1));
    dxi_2 = fhdym(xi_update(colIndcs{:}, 2));
    dxi_3 = fhdzm(xi_update(colIndcs{:}, 3));
    dxi_4 = (fhdxm(xi_update(colIndcs{:}, 2)) + fhdym(xi_update(colIndcs{:}, 1))) / 2;
    dxi_5 = (fhdxm(xi_update(colIndcs{:}, 3)) + fhdzm(xi_update(colIndcs{:}, 1))) / 2;
    dxi_6 = (fhdym(xi_update(colIndcs{:}, 3)) + fhdzm(xi_update(colIndcs{:}, 2))) / 2;
    
    % Take a step
    q(colIndcs{:}, 1) = q(colIndcs{:}, 1) - step_d * dxi_1;
    q(colIndcs{:}, 2) = q(colIndcs{:}, 2) - step_d * dxi_2;
    q(colIndcs{:}, 3) = q(colIndcs{:}, 3) - step_d * dxi_3;
    q(colIndcs{:}, 4) = q(colIndcs{:}, 4) - step_d * dxi_4;
    q(colIndcs{:}, 5) = q(colIndcs{:}, 5) - step_d * dxi_5;
    q(colIndcs{:}, 6) = q(colIndcs{:}, 6) - step_d * dxi_6;
    
    % Compute the projection
    q_norm            = sqrt(abs(q(colIndcs{:}, 1)).^2 + ...
                             abs(q(colIndcs{:}, 2)).^2 + ...
                             abs(q(colIndcs{:}, 3)).^2 + ...
                             2 * abs(q(colIndcs{:}, 4)).^2 + ...
                             2 * abs(q(colIndcs{:}, 5)).^2 + ...
                             2 * abs(q(colIndcs{:}, 6)).^2);
              
    denom             = max(1, (q_norm / alpha0));
    q(colIndcs{:}, 1) = q(colIndcs{:}, 1) ./ denom;
    q(colIndcs{:}, 2) = q(colIndcs{:}, 2) ./ denom;
    q(colIndcs{:}, 3) = q(colIndcs{:}, 3) ./ denom;
    q(colIndcs{:}, 4) = q(colIndcs{:}, 4) ./ denom;
    q(colIndcs{:}, 5) = q(colIndcs{:}, 5) ./ denom;
    q(colIndcs{:}, 6) = q(colIndcs{:}, 6) ./ denom;
    
    %----------------------------------------------------------------------
    % Primal Update
    % ---------------------------------------------------------------------
    
    % Compute the first order divergence
    divp = fhdxm(p(colIndcs{:}, 1)) + fhdym(p(colIndcs{:}, 2)) + ...
            fhdzm(p(colIndcs{:}, 3));
    
    % Take a step    
    w     = u - step_p * divp;
    u     = (w + (step_p * im)) ./ (1 + step_p);
    
    % Compute the second order divergence
    divq_1 = fhdxp(q(colIndcs{:}, 1)) + fhdyp(q(colIndcs{:}, 4)) + fhdzp(q(colIndcs{:}, 5));
    divq_2 = fhdxp(q(colIndcs{:}, 4)) + fhdyp(q(colIndcs{:}, 2)) + fhdzp(q(colIndcs{:}, 6));
    divq_3 = fhdxp(q(colIndcs{:}, 5)) + fhdyp(q(colIndcs{:}, 6)) + fhdzp(q(colIndcs{:}, 3));
    
    % Take a step
    xi(colIndcs{:}, 1) = xi(colIndcs{:}, 1) - step_p * (divq_1 - p(colIndcs{:}, 1));
    xi(colIndcs{:}, 2) = xi(colIndcs{:}, 2) - step_p * (divq_2 - p(colIndcs{:}, 2));
    xi(colIndcs{:}, 3) = xi(colIndcs{:}, 3) - step_p * (divq_3 - p(colIndcs{:}, 3));
    
    %----------------------------------------------------------------------
    % Update the auxiliary variables
    % ---------------------------------------------------------------------
    
    u_update  = 2 * u - u_old;
    xi_update = 2 * xi - xi_old;
    
    %----------------------------------------------------------------------
    % Cost calculations
    % ---------------------------------------------------------------------
    
    dataConsistencyCost = sum(abs(im(:) - u(:)).^2);
    tgv1Cost            = sum(abs(vectorizeArray(cat(5, du_x, du_y, du_z) - xi)).^2);
    tgv0Cost            = sum(vectorizeArray(abs(dxi_1) + abs(dxi_2) + abs(dxi_3) + ...
                           abs(dxi_4) + abs(dxi_5) + abs(dxi_6)));
    totalCost           = dataConsistencyCost + alpha1 * tgv1Cost + ...
                           alpha0 * tgv0Cost;
    funCost             = [funCost, totalCost];
    
    if iter>1
        CostConv=abs((funCost(end-1)-totalCost)/totalCost);
    end
    
    if verbose && mod(iter,100)==0
            fprintf('TGV2: iter = %4d, cost = %d, CostConv = %d\n', iter, totalCost,CostConv);
        
    end
   iter=iter+1;   
end %while (iter <= maxIter && CostConv>ConvThres )

imOut = u;

end
