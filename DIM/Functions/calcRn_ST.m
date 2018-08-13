function Rn = calcRn_ST(H0,L0,Beta,n)
% Rn = calcRn(H0,L0,Beta,n)
% 
% calculates the runup elevation using equation 9 in PH12
% 
% inputs 
% H0 : deep water significant wave height
% L0 : deep water wavelength
% Beta : foreshore slope
% n : number of standard deviations of swash that are used to calculate
% runup. Use n = 1 for R16 and n = 2 for R2.
%
% outputs
% Rn : time series of runup elevation
% 
%
% Kilian Vos 2017 WRL

c1 = 0.35*Beta.*sqrt(H0.*L0); % contribution of mean water level to runup
c2 = sqrt(H0.*L0.*(0.563*Beta.^2 + 0.004))*n/4; % contribution of n stds of swash to runup

Rn = 1.1*(c1 + c2);

end