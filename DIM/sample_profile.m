% Calculate functional relationship between x, z, dv

% Clean
clear
close all

% Load 
load june2016.mat

% Settings
dx = 0.5;
padding = 10; % Padding in m for model over-erosion
zMHW = 0.52; % Fixed MHW for the entire coastline

% Extract Data
sample = 209;
pfx = june2016.profiles(sample).pfx;
pf1 = flipud(june2016.profiles(sample).pf1);
pf2 = flipud(june2016.profiles(sample).pf2);
xb1 = june2016.profiles(sample).xb1;
zb1 = june2016.profiles(sample).zb1;
xb2 = june2016.profiles(sample).xb2;
zb2 = june2016.profiles(sample).zb2;

% Plot sample
subplot(121)
plot(pfx,pf1,'b-')
hold on
plot(pfx,pf2,'r-')
plot(xb1,zb1,'b.','MarkerSize',20)
plot(xb2,zb2,'r.','MarkerSize',20)
grid on
xlabel('Cross-shore distance (m)')
ylabel('Elevation (m)')
title(['Profile ' num2str(sample)])
xlim([100 250])
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[3 6 40 15])

% Calculate functional form
x = 0:dx:dx*(length(pf1)-1); % X grid with dx spacing
x_dV = x(x >= xb1 & x <= xb2+padding); % X coordinates between pre/post dune base + padding for over erosion
z1_dV = pf1(x >= xb1 & x <= xb2+padding); % Z coordinates between pre/post dune base for pf1
m_dV = (zb2-zb1)/(xb2-xb1); % Slope between pre/post dune base
m_z = m_dV*(0:0.5:range(x_dV))' + zb1*ones(length(z1_dV),1); % Coordinates of recession slope
dV = cumtrapz(x_dV,z1_dV-m_z); % Observed eroded dune volume

% Plot
subplot(122)
plot(m_z,dV)
hold on
plot([zb1 zb1],ylim,'k--')
plot([zb2 zb2],ylim,'k--')
grid on
xlabel('Dune Base Elevation (m)')
ylabel('Cumulative Dune Erosion Volume (m^3/m)')
title('Dune Base Elevation vs \Delta Dune Volume')
set(gca,'FontSize',16)

% Save
% zb = m_z;
% dv = dV;
% save sample_profile.mat zb dv