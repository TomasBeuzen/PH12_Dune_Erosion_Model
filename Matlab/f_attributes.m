function [attr] = f_attributes(profile,forcing,options_plot)
% attr = f_attributes(june2016,j,options_plot)
%
% This function calculates attributes (slope, volumes, total water level..etc) 
% 
% Inputs: 
%   -profile : structure containing all the info about a given profile
%   -forcing : structure containing all the forcing conditions at a given
%   profile
%  
% KV WRL 10.2017


dx = 0.5;                         % grid spacing
zMHW = 1.45 - 0.925;              % from MHL chart
idx_nan = isnan(profile.pf1);
z1 = flipud(profile.pf1(~idx_nan)); % pre-storm
z2 = flipud(profile.pf2(~idx_nan)); % post-storm
x = (0:dx:dx*(length(z1)-1))';

% Dune toes
xb1 = profile.xb1;
zb1 = profile.zb1;
xb2 = profile.xb2;
zb2 = profile.zb2;

% Calculate slopes

% Dune retreat slope
B1 = (zb2-zb1)/(xb2-xb1); % end-point slope calculation

% Mean high water (MHW) profile intersection
idx_beachFace = z1 > 0 & z1 < 1 & x < xb1;
xMHW = interp1(z1(idx_beachFace),x(idx_beachFace),zMHW);

% Compute the slope between the MHW intersection and the dune toe
B0 = (zMHW-zb1)/(xMHW-xb1); % end-point slope calculation

% Compute the slope between the MHW intersection and the 1.2 m contour
zMitch = 1.2;
idx_beachFace2 = z1 > 1 & z1 < 1.5 & x < xb1;
xMitch = interp1(z1(idx_beachFace2),x(idx_beachFace2),zMitch);
B0m = (zMitch-zMHW)/(xMitch-xMHW); % end-point slope calculation

if options_plot
    hfig1 = figure;
    hfig1.Units = 'normalized';
    hfig1.Position = [0 0 1 1];
    box on
    hold on
    grid on
    plot(x,z1,'b-','displayname','pre-storm')
    plot(x,z2,'r-','displayname','post-storm')
    plot(x,zMHW*ones(length(x),1),'g--','displayname','mean high water line')
    plot(xb1,zb1,'bo','markerfacecolor','b','markersize',5,'displayname','dune toe pre-storm')
    plot(xb2,zb2,'ro','markerfacecolor','r','markersize',5,'displayname','dune toe post-storm')
    hLegend = legend('show');
    hLegend.Location = 'best';
    plot(xMHW,zMHW,'go','markerfacecolor','g','markersize',5,'displayname','mean high water')
    plot([xMHW x(end)],[zMHW B0*(x(end)-xMHW) + zMHW],'k--')
    plot([xMHW x(end)],[zMHW B0m*(x(end)-xMHW) + zMHW],'k-')
    plot([xb1 x(end)],[zb1 B1*(x(end)-xb1) + zb1],'k:')
    text(x(end),B0*(x(end)-xMHW) + zMHW,'\beta_0')
    text(x(end),B0m*(x(end)-xMHW) + zMHW,'\beta_{0m}')
    text(x(end),B1*(x(end)-xb1) + zb1,'\beta_{1}')
    title([profile.site{1} ' profile # ' num2str(j)])
    xlabel('cross-shore distance [m]')
    ylabel('elevation [m]')
%     xlim([xMHW-5 xb2 + 5])
%     ylim([zMHW-0.5 zb2 + 0.5])
end

% Save attributes in a structure ATTR

% Dune toes
attr.x = x;
attr.z1 = z1;
attr.z2 = z2;
attr.dt1  = [xb1 zb1];
attr.dt2  = [xb2 zb2];
attr.MHW = [xMHW zMHW];

% Slopes
attr.B0 = B0;
attr.B0m = B0m;
attr.B1 = B1;

% Water levels and runup
H0 = forcing.H0_10m_max;
Tp = forcing.Tp_10m;
L0 = forcing.L0_10m;
zSWL = forcing.waterLevel;
attr.R2 = calcRn(H0,L0,B0m,2);
attr.SWL = zSWL;
attr.TWL = attr.R2 + zSWL;

% max(total water level) - initial dune toe elevation
attr.maxTWL = max(attr.TWL);
attr.peak_tide = max(zSWL);

%  Volumes
% Volume above retreat slope dV_o1
% Volume above flat slope dV_o2
if abs(xb2 - xb1) <= 1
    attr.dV_o1 = 0;
    attr.dV_o2 = 0;
else
    xV = x(x >= xb1 & x <= xb2);
    z1V = z1(x >= xb1 & x <= xb2);
    z2V = z2(x >= xb1 & x <= xb2);
    mV = (zb2-zb1)/(xb2-xb1);
    if mV < 0
        mV = 0;
    end
    zVm = mV*(0:0.5:range(xV))' + zb1*ones(length(z1V),1);
    attr.dV_o1 = trapz(xV,z1V-zVm);
    attr.dV_o2 = trapz(xV,z1V-zb1*ones(length(z1V),1));
end

% Profile and wave orientations
attr.orientation = profile.orientation;
attr.waveDm = median(forcing.Dm_10m);

end