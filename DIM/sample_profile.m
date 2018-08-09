% Calculate functional relationship between x, z, dv

% Clean
clear
close all

% Load 
load june2016.mat
load GP_query_B0.mat

%% Profile Data

% Settings
dx = 0.5;
padding = 10; % Padding in m for model over-erosion
zMHW = 1.45 - 0.925; % Fixed MHW (1.45 = MHW, 0.925

% Extract Data
sample = 209;
pfx = june2016.profiles(sample).pfx;
pf1 = flipud(june2016.profiles(sample).pf1);
pf2 = flipud(june2016.profiles(sample).pf2);
xb1 = june2016.profiles(sample).xb1;
zb1 = june2016.profiles(sample).zb1;
xb2 = june2016.profiles(sample).xb2;
zb2 = june2016.profiles(sample).zb2;

% Calculate functional form
x = [0:dx:dx*(length(pf1)-1)]'; % X grid with dx spacing
x_dv = x(x >= xb1 & x <= xb2+padding); % X coordinates between pre/post dune base + padding for over erosion
z1_dv = pf1(x >= xb1 & x <= xb2+padding); % Z coordinates between pre/post dune base for pf1
m_dv = (zb2-zb1)/(xb2-xb1); % Slope between pre/post dune base
zb = m_dv*(0:0.5:range(x_dv))' + zb1*ones(length(z1_dv),1); % Coordinates of recession slope (i.e, dune base)
dv = cumtrapz(x_dv,z1_dv-zb); % Observed eroded dune volume

% Plot
%%%% Profile
% subplot(121)
% plot(pfx,pf1,'b-')
% hold on
% plot(pfx,pf2,'r-')
% plot(xb1,zb1,'b.','MarkerSize',20)
% plot(xb2,zb2,'r.','MarkerSize',20)
% grid on
% xlabel('Cross-shore distance (m)')
% ylabel('Elevation (m)')
% title(['Profile ' num2str(sample)])
% xlim([100 250])
% set(gca,'FontSize',16)
% set(gcf,'units','centimeters','position',[3 6 40 15])
%%%% Functional Form
% subplot(122)
% plot(m_z,dV)
% hold on
% plot([zb1 zb1],ylim,'k--')
% plot([zb2 zb2],ylim,'k--')
% grid on
% xlabel('Dune Base Elevation (m)')
% ylabel('Cumulative Dune Erosion Volume (m^3/m)')
% title('Dune Base Elevation vs \Delta Dune Volume')
% set(gca,'FontSize',16)

%% Wave Data
% Waves
Hs   = june2016.forcing(sample).Hm0_offshore;
Tp   = june2016.forcing(sample).Tp_offshore;
L    = june2016.forcing(sample).L0_offshore;
zSWL = june2016.forcing(sample).waterLevel;
% Slope
idx_xb1 = find(x == xb1); % Index of xb1
idx_beachFace = pf1 > 0 & pf1 < 1 & x < xb1; % Indices of beach face
xMHW = interp1(pf1(idx_beachFace),x(idx_beachFace),zMHW); % x location of MHW
b0 = (zMHW-zb1)/(xMHW-xb1); % Slope of dune toe to MHW
% Stockdon Runup
R_st = calcRn(Hs,L,b0,2) + zSWL;
% GP Runup
index = sample*108-107:1:sample*108;
R_gp = calcRn_GP(GP_mean,GP_sigma,index,'mean') + zSWL;

% Save
save sample_profile.mat zb dv Tp R_st R_gp