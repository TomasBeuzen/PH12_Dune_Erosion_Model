 function out = f_PH12model_single(profile,forcing,time,Cs,rB,plot_option,folder)
% out = f_PH12model_single(june2016,Cs,ratio_B,plot_option)
%
% This function runs the PH12 model on a single profile. 
%
% Inputs: 
%   -profile : structure containing all the info about a given profile
%   -forcing : structure containing all the forcing conditions at a given
%   profile
%   -time : time vector in datenums
%   -Cs : Cs parameter
%   -rB : ratio between dune retreat slope and foreshore slope
%   - plot_option : to plot the output
%   - folder : where to save the plot
%
%
% KV WRL 10.2017

dx = 0.5;                         % grid spacing
zMHW = 1.45 - 0.925;              % from MHL chart
options.plot = plot_option;       % 1 if you want to save the plots

out.t_c = 0;
out.nans = 0;
out.too_steep = 0;
out.xskip = 0;
thresh_skip = 10;
out.break = 0;

    
% Load profile coordinates (x-z)
idx_nan = isnan(profile.pf1);
z1 = flipud(profile.pf1(~idx_nan)); % pre-storm
z2 = flipud(profile.pf2(~idx_nan)); % post-storm
x = [0:dx:dx*(length(z1)-1)]';

% Dune toes (xb1,zb1 and xb2,zb2)
idx_xb1 = find(x == profile.xb1);
xb1 = profile.xb1;
zb1 = profile.zb1;
xb2 = profile.xb2;
zb2 = profile.zb2;
    
% Calculate observed volume difference
if abs(xb2 - xb1) <= 1
    out.dV_o = 0;
    out.dV_o2 = 0;
else
    xV = x(x >= xb1 & x <= xb2);
    z1V = z1(x >= xb1 & x <= xb2);
    mV = (zb2-zb1)/(xb2-xb1);
    % if negative slope set slope to 0 (only volume above dune toe)
    if mV < 0
        mV = 0;
    end
    zVm = mV*(0:0.5:range(xV))' + zb1*ones(length(z1V),1);
    out.dV_o = trapz(xV,z1V-zVm);
    out.dV_o2 = trapz(xV,z1V-zb1*ones(length(z1V),1));
end
    
% Mean high water (MHW) profile intersection
idx_beachFace = z1 > 0 & z1 < 1 & x < xb1;
xMHW = interp1(z1(idx_beachFace),x(idx_beachFace),zMHW);

% Compute the slope between the dune toe and the MHW
B0 = (zMHW-zb1)/(xMHW-xb1); % end-point slope calculation

% Compute the slope between the 1.2m contour and the MHW
z12 = 1.2;
idx_beachFace = z1 > 1 & z1 < 1.5 & x < xb1;
x12 = interp1(z1(idx_beachFace),x(idx_beachFace),z12);
B0m = (z12-zMHW)/(x12-xMHW); % end-point slope calculation

% Forcing conditions
t = time;
t = (t - t(1))*24;          % convert to hours
dt = diff(t(1:2))*3600;     % dt in seconds

H0 = forcing.H0_10m_max;    % deep-water Hs [m]
Tp = forcing.Tp_10m;        % deep-water Tp [s]
L0 = forcing.L0_10m;        % deep-water Lo [m]
zSWL = forcing.waterLevel;  % water-level   [m]
    
zb(1) = zb1;            % initial dune toe elevation
idx_start = idx_xb1;     % initial grid point, where dune toe is
xb(1) = x(idx_start);   % initial dune toe position
    
% pre-computed trajectory (using Bt)
Bt = rB * B0;        % slope for pre-computed dune trajectory
zb_traj = [nan(1,idx_start -1)' ; (Bt*(x(idx_start:end) - x(idx_start)) + zb(1)*ones(length(x(idx_start:end)),1))];
    
% Initialise vectors to run model
zR = zeros(length(t),1);         % runup [m]
zTotal = zeros(length(t),1);     % total water level [m]
Nc = zeros(length(t),1);         % number of collisions in 1h [-]
dV = zeros(length(t),1);         % eroded volume in 1h [m^3/m]
xb = [xb1 ; zeros(length(t),1)]; % x dune base [m]
zb = [zb1 ; zeros(length(t),1)]; % z dune base [m]
deltax = zeros(length(t),1);     % dune toe recession in 1h [m]
Vt = 0;                          % total eroded volume (sum of dV)
dVt = 0;                         % cumulated volume when not enough to erode a grid point
    
ti = 0;

% Time loop during storm
while ti < length(t)
    
    ti = ti + 1;
    
    % Compute cumulative volume between the initial profile and dune trajectory
    clear Vc
    Vc = [nan(idx_start-1,1) ; cumtrapz(x(idx_start:end),z1(idx_start:end)-zb_traj(idx_start:end))];
    Vc = Vc - Vc(idx_start);
    dVc = [diff(Vc) ; NaN]; % compute the difference of the cumulative volume
    
    % if the dune retreat trajectory is above the initial profile
    if dVc(idx_start) < 0
        
        if any(dVc(idx_start:end) > 0) % if there is intersection
            % skip until the dune retreat trajectory is above the initial profile
            idx_skip = find(dVc(idx_start+1:end) > 0,1,'first');
            % count how many meters are skipped
            out.xskip = out.xskip + x(idx_start + idx_skip) - x(idx_start);
        
            idx_start = idx_start + idx_skip;
        
            % re-compute the cumulative volume difference at the new position
            clear Vc
            Vc = [nan(idx_start-1,1) ; cumtrapz(x(idx_start:end),z1(idx_start:end)-zb_traj(idx_start:end))];
            Vc = Vc - Vc(idx_start);
            dVc = [diff(Vc) ; NaN];
            
        else % if no intersection stop there
            
            % update the variables zb, xb and deltax
            zb(ti+1:end) = zb_traj(idx_start)*ones(length(zb(ti+1:end)),1);
            xb(ti+1:end) = x(idx_start)*ones(length(zb(ti+1:end)),1);
            deltax(ti+1) = xb(ti+1)-xb(ti);
            
            % flag that the loop was broken
            out.break = 1;
            break;
            
        end
            
        
    end
        
    % Calculate zR, zTot, Nc and dV at time t
    zR(ti) = calcRn(H0(ti),L0(ti),B0m,2);
    zTotal(ti) = zR(ti) + zSWL(ti);
    [Nc(ti) p(ti)] = calcNc(zb(ti),zSWL(ti),dt,Tp(ti),H0(ti),L0(ti),B0m);
    dV(ti) = 4*Cs*(max(zTotal(ti)-zb(ti),0))^2*Nc(ti);
    Vt = dV(ti) + Vt; % sum up the eroded volume
    
    % Compute the time in collision regime
    if zTotal(ti) >= zb(ti)
        out.t_c = out.t_c + 1;
    end
    
    % auxiliary variable for the eroded volume
    dVt = dV(ti) + dVt;
    
    % if dune got eroded
    if dVt > 0 && dVt >= Vc(idx_start+1)
               
        if ~any(Vc - dVt > 0) % if dVt is always larger than Vc
            
            idx_start = find(Vc == max(Vc),1,'first');
            
            % update the variables zb, xb and deltax
            zb(ti+1:end) = zb_traj(idx_start)*ones(length(zb(ti+1:end)),1);
            xb(ti+1:end) = x(idx_start)*ones(length(zb(ti+1:end)),1);
            deltax(ti+1) = xb(ti+1)-xb(ti);
            
            % flag that the loop was broken
            out.break = 1;
            break;
        
        else
            % erode until interception between dVt and Vc
            idx_erosion = find(Vc - dVt > 0,1,'first'); % erode until Vc > dVt
            
            if any(dVc(idx_start:idx_erosion) < 0) % in case the profile goes below the dune retreat trajectory
                idx_erosion = idx_start + find(dVc(idx_start+1:idx_erosion) < 0,1,'first');
            end
            
            % Go to grid point before and keep remaining eroded volume for
            % the next hit
            if idx_erosion - idx_start > 1
                idx_erosion = idx_erosion - 1;
            end
            
            dVt = dVt - Vc(idx_erosion); % keep the extra erosion
            idx_start = idx_erosion;    % update starting point for next iteration
            
        end
                
    end
    
    zb(ti+1) = zb_traj(idx_start);
    xb(ti+1) = x(idx_start);
    deltax(ti) = xb(ti+1)-xb(ti);
end
    
% Save results (out structure)
out.xb_m = xb(end);              % dune toe x-position MODEL
out.zb_m = zb(end);              % dune toe z-position MODEL
out.dV_m = Vt - dVt;             % eroded volume MODEL
out.dx_m = sum(deltax);          % dune toe recession MODEL
out.dz_m = zb(end) - zb1;        % dune toe elevation difference MODEL
out.xb_o = [xb1 xb2];            % dune toe x-position DATA
out.zb_o = [zb1 zb2];            % dune toe z-position DATA
out.dx_o = xb2 - xb1;            % dune toe recession DATA
out.dz_o = zb2 - zb1;            % dune toe elevation difference DATA

out.error_dx = xb2 - xb(end);    % error in dune toe x-position
out.error_dz = zb2 - zb(end);    % error in dune toe z-position
out.error_dV = out.dV_o - out.dV_m; % error in eroded volume above dune toe

% Compute BSS for deltaX
out.BSS = 1 - (abs(xb(end) - xb2))/((xb2 - xb1));
out.ABSS = 1 - (abs(xb(end) - xb2))/((xb2 - (xb1 + out.xskip + 0.01)));
 
% Compute normilised volume difference
out.dVnorm = 1 - (abs(out.dV_m - out.dV_o))/((out.dV_o));
    
% flag if too many grid points were skipped
if out.xskip >= thresh_skip
    out.too_steep = 1;
end
    
% print a figure
if options.plot %|| out.break %|| out.too_steep
    
    str = sprintf('%s profile #%d',profile.site{1},j);
    
    hfig1 = figure;
    hfig1.Units = 'normalized';
    hfig1.Position = [0 0 0.6 0.9];
    
    % Subplot on the left
    hfig1sub1 = subplot(8,1,[1 2 3]);
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
    plot([xMHW x(end)],[zMHW B0*(x(end)-xMHW) + zMHW],'--','color',[0.5 0.5 0.5])
    %     plot([xMHW x(end)],[zMHW B0r*(x(end)-xMHW) + zMHW],'--','color',[0.3 0.3 0.3])
    plot([xMHW x(end)],[zMHW B0m*(x(end)-xMHW) + zMHW],'--','color','k')
    text(x(end),B0*(x(end)-xMHW) + zMHW,'\beta_0')
    %     text(x(end),B0r*(x(end)-xMHW) + zMHW,'\beta_{0r}')
    text(x(end),B0m*(x(end)-xMHW) + zMHW,'\beta_{0} mitch')
    %     plot(xb,zb,'ko','markersize',3,'markerfacecolor','k')
    plot([xb1 x(end)],[zb1 Bt*(x(end)-xb1)+zb1],'k:')
    plot(xb(end),zb(end),'ko','markerfacecolor','k','markersize',3)
    xlabel('cross-shore distance [m]')
    ylabel('elevation [m]')
    %         title(sprintf('%s (Cs = %.4f, rB = %.2f)    tc = %d h \n Edx = %.1f m (xjump = %.1f m), BSSdx = %.2f (ABSS = %.2f) \n EdV = %.2f m^3/m (%.d %% error)',...
    %             str,Cs,rB,out.t_c,out.error_dx,out.xskip,out.BSS,out.ABSS,out.error_dV, round(100*out.error_dV/out.dV_o)),'FontWeight','normal');
    title(sprintf('%s (Cs = %.4f, rB = %.2f)    tc = %d h \n Edx = %.1f m (dxobs = %.1f, dxmod = %.1f), BSSdx = %.2f \n EdV = %.2f m^3/m (dVobs = %.2f, dVmod = %.2f), BSSdV = %.2f ',...
        str,Cs,rB,out.t_c,out.error_dx,out.dx_o,out.dx_m,out.BSS,out.error_dV,out.dV_o,out.dV_m,out.dVnorm),'FontWeight','normal');
    
    % xlim([xb1 - 10 xb2 + 10])
    % Subplot on the right
    hfig1sub2 = subplot(8,1,5);
    hold on;grid on;box on;
    title('zTotal')
    ylabel('zTotal [m]')
    plot(t,zTotal,'b-')
    plot(t,zb(2:end),'r-')
    hfig1sub2.XTickLabel = [];
    hfig1sub4 = subplot(8,1,6);
    hold on;grid on;box on;
    title('N_C')
    ylabel('N_C [-]')
    plot(t,Nc,'b-')
    hfig1sub4.XTickLabel = [];
    hfig1sub6 = subplot(8,1,7);
    hold on;grid on;box on;
    title('\Delta V')
    ylabel('dV [m^3/m]')
    plot(t,dV,'b-')
    hfig1sub6.XTickLabel = [];
    hfig1sub8 = subplot(8,1,8);
    hold on;grid on;box on;
    title('\Delta x')
    ylabel('dx [m]')
    plot(t,deltax,'b-')
    xlabel('time [hrs]')
    %         print(hfig1,fullfile(folder,str),'-djpeg','-r300')
    %         close(hfig1)
end
        
end
