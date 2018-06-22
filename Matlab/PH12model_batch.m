% Runs PH12 model on all the profiles and plots performance statistics
%
% KV WRL 2017


clear
close all
clc
format bank
format compact
addpath(genpath('.\functions'))
addpath(genpath('.\data'))

% Load variables
load('june2016.mat')
load('calib_rB05_B0m.mat'); % calibration data
idx_keep = calib.idx_keep; % 600 profiles to use for data analysis

% Remove profiles that had no erosion
june2016.profiles(calib.idx_noerosion) = [];
june2016.forcing(calib.idx_noerosion) = [];

% Remove profiles that could not be calibrated
% june2016.profiles(~calib.idx_keep) = [];
% june2016.forcing(~calib.idx_keep) = [];
june2016.profiles(~idx_keep) = [];
june2016.forcing(~idx_keep) = [];

Np = length(june2016.profiles);

% Cs vs maxTWL - zb regression line
a = calib.regression.dx.p1;
b = calib.regression.dx.p2;
Cs_reg = @(x) 10.^((x-b)/a);

rB = calib.rB;
time = june2016.info.time;

hW = waitbar(0,'Please wait...');
for i = 1:Np
    
    % Update waitbar
    hW = waitbar(i/Np);
    hW.Children.Title.String = [num2str(i) '/' num2str(Np)];
    
    profile = june2016.profiles(i);
    forcing = june2016.forcing(i);
    
    % Compute Cs from regression
    attr = f_attributes(profile,forcing,0);
    var = attr.maxTWL - attr.dt1(2);
    Cs = Cs_reg(var);
    
    % run model
    out = f_PH12model_single(profile,forcing,time,Cs,rB,0); 
    
    % store outputs
    idx_toosteep(i) = out.too_steep;
    idx_break(i) = out.break;
    tc(i) = out.t_c;
    xo(i) = out.xb_o(2);
    zo(i) = out.zb_o(2);
    xi(i) = out.xb_o(1);
    zi(i) = out.zb_o(1);
    xm(i) = out.xb_m;
    zm(i) = out.zb_m;
    dVm(i) = out.dV_m;
    dVo(i) = out.dV_o;
    xskip(i) = out.xskip;
    BSS(i) = out.BSS;
    dV_norm(i) = out.dVnorm;
    
end
delete(hW);

% Calculate statistics
retreat_obs = xo - xi;
retreat_model = xm - xi;
dV = dVm - dVo;

% Remove outliers
idx1 = find(abs(retreat_obs - retreat_model) > mean(retreat_obs - retreat_model)+3*std(retreat_obs - retreat_model));
retreat_model(idx1) = []; 
retreat_obs(idx1) = [];
dV(idx1) = [];
dV_norm(idx1) = [];
BSS(idx1) = [];

% Percentage of positive BSS and dVnorm
nBSS = sum(BSS >= 0)/length(BSS); 
ndVnorm = sum(dV_norm >= 0)/length(dV_norm);

%% Figure 1 
figure('units','normalized','position',[0 0 1 1])
% subplot 1:1 plot dx
h(1) = subplot(221);
hold on; grid on; box on ;axis equal;
rho = corrcoef(retreat_model,retreat_obs);
title(['1:1 plot \Deltax,   \rho = ' num2str(rho(2,1),2)])
xlabel('\Deltax_{obs}');ylabel('\Deltax_{model}');
plot(retreat_obs,retreat_model,'b.');
plot([0 max(retreat_model)], [0 max(retreat_model)],'r--');
% subplot error distribution dx
h(2) = subplot(222);
hold on; grid on; box on ;
h1 = histogram(retreat_obs - retreat_model);
xlabel('error in x [m]');ylabel('counts')
h1.BinWidth = 0.5;
h1.FaceColor = [.5 .5 .5];
h1.LineWidth = 0.5; 
str = sprintf('RMSE = %.1f m\n sigma = %.1f m\n median = %.1f m\n q_{90} = %.1f m\n n = %d',rms(retreat_obs - retreat_model),std(retreat_obs - retreat_model),median(retreat_obs - retreat_model),quantile(abs(retreat_obs - retreat_model),0.9),Np);
text(0.5*max(retreat_obs - retreat_model),0.75*max(h1.Values),str,'fontweight','normal')
title('Error distribution (\Deltax_{obs} - \Deltax_{model})')
% subplot histogram relative error dx
h(3) = subplot(223);
hold on; grid on; box on ;
xlabel('dx relative error');ylabel('counts');
title(sprintf('Relative error dx \n %d %% >= 0 , mean = %.2f , std = %.2f',round(nBSS*100),mean(BSS(BSS >=0)),std(BSS(BSS >=0))));
histogram(BSS(BSS >=0),'BinWidth',0.05)

% subplot histogram relative error dV
h(4) = subplot(224);
hold on; grid on; box on ;
xlabel('dV relative error ');ylabel('counts');
title(sprintf('Relative error dV \n %d %% >=0 , mean = %.2f , std = %.2f',round(ndVnorm*100),mean(dV_norm(dV_norm >=0)),std(dV_norm(dV_norm >=0))));
histogram(dV_norm(dV_norm >=0),'BinWidth',0.05)

%% Figure 2
sites = [june2016.profiles.site]';
sites(idx1) = [];
sites_unique = unique(sites,'stable');

r2 = mean((retreat_model-mean(retreat_model)).*(retreat_obs-mean(retreat_obs)))/(std(retreat_obs)*std(retreat_model));
rmse = sqrt(mean((retreat_model-retreat_obs).^2));
rmsm = sqrt(mean(retreat_obs.^2));
sci     = rmse/max(rmsm,abs(mean(retreat_obs)));
relbias = mean(retreat_model-retreat_obs)/max(rmsm,abs(mean(retreat_obs)));
str = sprintf('Dune toe retreat: R^2 = %.2f , RMSE = %.2f m    (sci = %.2f , relbias = %.2f)',r2,rmse,sci,relbias);

width = 16;
height = 6;
margHor = [1.5 0.5];
margVer = [0.5 1.5];
gaps = [0 1.5];


figure
colors = colormap(lines(length(sites_unique)));
counter = 1;
subplot(3,1,1)
hold on; box on;axis tight;
title(str)
plot(0,100,'o','color','k','markerfacecolor','k','markersize',3,'displayname','modelled retreat')
hLeg = legend('show');
ylim([0 30])
ylabel('observed retreat [m]')
set(gca,'XTickLabel',[])
bool1 = 1; bool2 = 1;
counterTicks = 1;
for i = 1 :length(sites_unique)
    if counter > 200 && bool1
        set(gca,'XTick',xticks,'XTickLabel',xticklabels);
        xtickangle(45)
        subplot(3,1,2)
        hold on; box on;axis tight;
        ylim([0 30])
        ylabel('retreat distance [m]')
        set(gca,'XTickLabel',[])
        bool1 = 0;
        counterTicks = 1;
    end
    if counter > 400 && bool2
        set(gca,'XTick',xticks,'XTickLabel',xticklabels);
        xtickangle(45)
        subplot(3,1,3)
        hold on; box on;axis tight;
        ylim([0 30])
        ylabel('retreat distance [m]')
        set(gca,'XTickLabel',[])
        bool2 = 0;
        counterTicks = 1;
    end
    idx = find(strcmpi(sites,sites_unique{i}));
    xticks(counterTicks) = counter + round(length(idx)/2);
    xticklabels{counterTicks} = sites_unique{i};
    counterTicks = counterTicks + 1;
    
    for j = 1:length(idx)
       hbar = bar(counter,retreat_obs(idx(j)),1);
       hbar.FaceColor = colors(i,:);
       plot(counter,retreat_model(idx(j)),'o','color','k','markerfacecolor','k','markersize',3)
       counter = counter + 1;
    end
    counter = counter + 1;
end
set(gca,'XTick',xticks,'XTickLabel',xticklabels);
xtickangle(45)

%% Figure 3
% figure
% hold on; grid on; box on ;
% bar(retreat_obs,1);
% plot(retreat_model,'ro-','MarkerFaceColor','r','linewidth',0.5,'markersize',5)
% r2 = mean((retreat_model-mean(retreat_model)).*(retreat_obs-mean(retreat_obs)))/(std(retreat_obs)*std(retreat_model));
% rmse = sqrt(mean((retreat_model-retreat_obs).^2));
% rmsm = sqrt(mean(retreat_obs.^2));
% sci     = rmse/max(rmsm,abs(mean(retreat_obs)));
% relbias = mean(retreat_model-retreat_obs)/max(rmsm,abs(mean(retreat_obs)));
% str = sprintf('Dx: R^2 = %.2f , rms = %.2f , sci = %.2f , relbias = %.2f ',r2,rmse,sci,relbias);
% title(str);ylabel('\Deltax [m]')


