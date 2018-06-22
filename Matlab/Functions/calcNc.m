function [Nc p] = calcNc(zb,zSWL,dt,Tp,H0,L0,Beta)
% Nc = calcNc(zb,zR,zSWL,dt,Tp,H0,L0,Beta)
% 
% calculates the number of collisions Nc using equation 12 in PH12
% 
%
% inputs 
% zb = elevation of the dune base
% zSWL = elevation of the Still Water Line (tide + surge)
% dt = time spacing
% Tp = peak wave period
% H0 : deep water significant wave height
% L0 : deep water wavelength
% Beta : foreshore slope
%
% outputs
% Nc : number of collisions during time dt
% 
%
% Kilian Vos 2017 WRL

eta = 0.35*Beta*sqrt(H0*L0); % mean water level
sigma_s = sqrt(H0*L0*(0.563*Beta^2 + 0.0004))/2; % 1 standard deviation of swash
% 
% figure(100)
% hold on
% box on
% grid on
% x = eta + zSWL - 3*sigma_s:0.01:eta + zSWL + 3*sigma_s;
% y = pdf('Normal',x,eta + zSWL,sigma_s);
% plot(x,y,'r-')
% plot([zb zb],[0 max(y)],'b-')
%
p =1-cdf('norm',zb,eta + zSWL,sigma_s);
Nc = p*(dt/Tp);

end