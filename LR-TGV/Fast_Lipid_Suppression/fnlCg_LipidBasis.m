function x = fnlCg_LipidBasis(x0,params)
%-----------------------------------------------------------------------
%
% given k-space measurements d, and a fourier operator F the function 
% finds the image m that minimizes:
%
% Phi(x) = || F * m - k_space ||^2 +       
%          lambda * sum_over((x,y) in brain_mask){ ||L * m(x,y)||_1 }
%
% the optimization method used is non linear conjugate gradient with 
% backtracking line-search.
% 
% params.FT : DFT operator 
%
% params.data :  k-space data of lipid contaminated image
%
% params.Bmask : brain mask
% params.Lipid : lipid basis functions
%
% Original code: Michael Lustig, sparse MRI
%-------------------------------------------------------------------------

x = x0;

% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;     beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;


% compute g0  = grad(Phi(x))
g0 = wGradient(x,params);
dx = -g0;


% iterations
while(1)

	[FT_x, FT_dx] = preobjective(x, dx);
    
	f0 = objective(FT_x, FT_dx, x, dx, 0, params);
	t = t0;
        [f1,consis,lip]  =  objective(FT_x, FT_dx, x, dx, t, params);
	
	lsiter = 0;

	while ( f1 > f0 - alpha * t * abs( dot(g0(:),dx(:)) ) ) ^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1,consis,lip]  =  objective(FT_x, FT_dx, x, dx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	if lsiter > 2
		t0 = t0 * beta;
	end 	
	if lsiter < 1
		t0 = t0 / beta;
	end

	x = (x + t*dx);

	%---------------------------------------------------------------
	disp(sprintf('%d   , obj: %f , consistency: %f ,  lipid: %f ', k, f1, consis, lip));
	%---------------------------------------------------------------
	
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
	bk = dot(g1(:),g1(:))/(dot(g0(:),g0(:))+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
		break;
	end

end

disp(sprintf('L2 of consistency: %f ,  L1 of lipid: %f ', sqrt(consis), lip/params.xfmWeight));

return;


% -------------------------------------------------------------------------
function [FT_x, FT_dx] = preobjective(x, dx)
% precalculates transforms to make line search cheap

FT_x = params.FT * x;
FT_dx = params.FT * dx;

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [res,obj,XFM] = objective(FT_x, FT_dx, x, dx, t, params)
%calculates the objective function

p = params.pNorm;

temp = (FT_x + t* FT_dx) - params.data;
obj = norm(temp(:))^2;


G = 0;

if params.xfmWeight > 0
    
    for ay = 1:size(x,1)
        for cey = 1:size(x,2)
            if params.Bmask(ay,cey) == 0
               continue; 
            end

            % gradient of L1 norm at each (x,y)
            Dx = LipidXFM( reshape( x(ay,cey,:) + t*dx(ay,cey,:) ,[size(x,3),1]), params );
            G = G + sum((Dx.*conj(Dx) + params.l1Smooth).^(p/2));
        end
    end
    
end

XFM = G*params.xfmWeight;

res = obj + XFM;

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function grad = wGradient(x,params)

gradXFM = 0;

gradObj = gOBJ(x,params);

if params.xfmWeight
    gradXFM = gLipid(x,params);
end

grad = gradObj +  params.xfmWeight.*gradXFM; 

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

F_x = params.FT * x;

gradObj =  params.FT'*( F_x - params.data) ;
gradObj = 2 * gradObj;

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function grad = gLipid(x,params)
% computes gradient of lipid basis operator

grad = zeros(size(x));

for ay = 1:size(x,1)
    for cey = 1:size(x,2)
        if params.Bmask(ay,cey) == 0
           continue; 
        end
        
        % gradient of L1 norm at each (x,y)
        Dx = LipidXFM( reshape(x(ay,cey,:), [size(x,3),1]), params );       
        grad(ay,cey,:) = sign(Dx)'*params.Lipid;        
    end
end

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function XFM = LipidXFM(spect,params)
% compute the inner product of given spectrum with the lipid basis functions
% 'Lipid Transform'
% each row of params.Lipid : lipid spectrum

XFM = params.Lipid * conj(spect);

end
% -------------------------------------------------------------------------

end