function [K,Interp_Beta,pp_Lip,pp_Consist] = ComputeLCurvature( Betas,Lipid_Sq,Consist_Sq , NbPt)
%COMPUTELCURVATURE Summary of this function goes here
%   Detailed explanation goes here
pp_Consist=pchip(Betas,Consist_Sq);%spline(Betas,Consist_Sq);
pp_Lip=pchip(Betas,Lipid_Sq);%spline(Betas,Lipid_Sq);
Interp_Beta=exp(linspace(min(log(Betas)),max(log(Betas)),NbPt));
% in the notetation of Hansen 2000
rho=ppval(pp_Consist,Interp_Beta);
rhop=ppval(fnder(pp_Consist,1),Interp_Beta);
rhopp=ppval(fnder(pp_Consist,2),Interp_Beta);
eta=ppval(pp_Lip,Interp_Beta);
etap=ppval(fnder(pp_Lip,1),Interp_Beta);
etapp=ppval(fnder(pp_Lip,2),Interp_Beta);

K= 2*eta.*rho./etap.*(Interp_Beta.^2.*etap.*rho + 2*Interp_Beta.*eta.*rho+Interp_Beta.^4.*eta.*etap)./(Interp_Beta.^2.*eta.^2 + rho.^2).^(3.0/2);
%K=2*(rhop.*etapp-rhopp.*etap)./((rhop.^2+etap.^2).^(3.0/2));

end

