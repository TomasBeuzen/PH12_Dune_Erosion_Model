%% Main 4: Calibrate Cs parameter with the best dx and dV

clear
close all force
clc
format bank
format compact
addpath(genpath('./Functions'))
addpath(genpath('./Data'))

load('june2016.mat')

Np = length(june2016.profiles);

% Index for profiles that did not erode during the storm (less than a
% meter)
idx_noerosion = abs([june2016.profiles.xb1] - [june2016.profiles.xb2]) <= 1;

% Generate N different Cs values in the range specified
N = 20;
Cs_range = [1e-5 5e-1];
Cs_values = logspace(log10(Cs_range(1)), log10(Cs_range(2)), N);

% Choose slope ratio
rB = 0.5;

idx_skip = nan(Np,N);
dx = nan(Np,N);
ABSS = nan(Np,N);
BSS = nan(Np,N);
dV = nan(Np,N);
dV_norm = nan(Np,N);
M = nan(Np,12);

hW = waitbar(0,'Please wait...');

time = june2016.info.time;

for i = 1:Np

    hW = waitbar(i/Np);
    hW.Children.Title.String = [num2str(i) '/' num2str(Np)];
    
    if idx_noerosion(i)
        continue
    end
    
    profile = june2016.profiles(i);
    forcing = june2016.forcing(i);

    for ii = 1:N
        % run model
         out = f_PH12model_single(profile,forcing,time,Cs_values(ii),rB,0);
        
        % save performance
        dx(i,ii) = out.error_dx;
        BSS(i,ii) = out.BSS;
        ABSS(i,ii) = out.ABSS; %adjusted brier skill score (if grid points were skipped)
        dV(i,ii) = out.error_dV;
        dV_norm(i,ii) = out.dVnorm;
        
        if out.too_steep || out.break
            idx_skip(i,ii) = 1;
        else
            idx_skip(i,ii) = 0;
        end
        
    end
    
    % Compute attributes and put into matrix M
    attr = f_attributes(profile,forcing,0);
    M(i,:) = [i attr.dt1(2) attr.maxTWL attr.peak_tide attr.orientation, ...
        attr.waveDm attr.B0 attr.B0m attr.dV_o1 attr.dV_o2 attr.maxTWL-attr.dt1(2),...
        attr.dt1(2)-attr.peak_tide];

end
delete(hW);
%% Save outputs into structure
calib.Cs_values = Cs_values;
calib.Cs_range = Cs_range;
calib.rB = rB;
calib.slopeR2 = 'B0m';

calib.dx = dx;
calib.BSS = BSS;
calib.ABSS = ABSS;
calib.dV = dV;
calib.dV_norm = dV_norm;

calib.idx_skip = idx_skip;
calib.idx_noerosion = idx_noerosion;

calib.M = M;

% save('calib1_final.mat','calib');