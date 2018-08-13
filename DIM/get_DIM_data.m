% Get profile and wave data to run Dune Impact Model
% There are 1000 profiles available in the june2016.mat dataset. This code
% filters through these to find which are applicable for the DIM model. It
% extracts required morphological data (dune base, dune erosion) and wave
% data (wave period, runup).
%
% TB 2018
%
% ISSUES
% Evaluate all profiles but NaN bad ones and flag

%% Pre-processing
% Clean
clear
close all

% Paths
addpath('Data')
addpath('Functions')

% Load 
load june2016.mat
load GP_query_B0.mat

% Hard-wired Settings/Parameters
Np = length(june2016.profiles); % Number of profiles
dx = 0.5;                       % Grid spacing
dt = 108;                       % Timesteps in data (i.e. storm = 108hrs)
padding = 10;                   % Padding in m to account for model over-erosion
zMHW = 1.45 - 0.925;            % Fixed MHW (1.45 = MHW, 0.925
draws = 10;                     % # of draws to take from GP to form an ensemble

% Pre-allocation
flag = false(Np, 1);        % Flag for good profiles
R_gp_draws = nan(dt,draws); % Initialise GP draws matrix
data = struct(...           % Output structure
    'zb',{},...
    'zb_final',{},...
    'dv',{},...
    'dv_obs',{},...
    'Tp',{},...
    'R_st',{},...
    'R_gp',{},...
    'R_gp_draws',{});

%% Data Loop
for i=161:Np

    %%%% Morphological Data
    % Extract Data
    pfx = june2016.profiles(i).pfx;        
    pf1 = flipud(june2016.profiles(i).pf1); % Pre-storm profile
    pf2 = flipud(june2016.profiles(i).pf2); % Post-storm profile
    xb1 = june2016.profiles(i).xb1;         % Pre-storm dune base x
    zb1 = june2016.profiles(i).zb1;         % Pre-storm dune base z
    xb2 = june2016.profiles(i).xb2;         % Post-storm dune base x
    zb2 = june2016.profiles(i).zb2;         % Post-storm dune base z
    m   = (zb2-zb1)/(xb2-xb1);              % Slope between pre/post dune base
    
    if abs(xb2 - xb1)<=1 || m < 0
        data(i,1).zb = NaN;
        data(i,1).zb_final = NaN;
        data(i,1).dv = NaN;
        data(i,1).dv_obs = NaN;
        data(i,1).Tp = NaN;
        data(i,1).R_st = NaN;
        data(i,1).R_gp = NaN;
    else
        flag(i) = true;
        % Calculate Dune Base and Volume Erosion
        x = [0:dx:dx*(length(pf1)-1)]';                        % X grid with dx spacing
        xdv = x(x >= xb1 & x <= xb2+padding);                  % X coordinates between pre/post dune base + padding for over erosion
        zdv = pf1(x >= xb1 & x <= xb2+padding);                % Z coordinates between pre/post dune base for pf1
        zb = m*(0:0.5:range(xdv))' + zb1*ones(length(zdv),1);  % Coordinates of recession slope (i.e, dune base)
        dv = cumtrapz(xdv,zdv-zb);                             % Cumulative eroded dune volume along projected slope
        zb_final = zb2;
        dv_obs = dv(xdv==xb2);                                 % Observed dune volume (no padding)
        % Calcuate Beach Slope
        idx_xb1 = find(x == xb1);                                 % Index of xb1
        idx_beachFace = pf1 > 0 & pf1 < 1 & x < xb1;              % Indices of beach face
        xMHW = interp1(pf1(idx_beachFace),x(idx_beachFace),zMHW); % x location of MHW
        b0 = (zMHW-zb1)/(xMHW-xb1);                               % Slope of dune toe to MHW
        
        %%% Wave Data
        % Waves
        Hs = june2016.forcing(i).Hm0_offshore; % Offshore significant wave height
        Tp = june2016.forcing(i).Tp_offshore;  % Peak wave period
        L = june2016.forcing(i).L0_offshore;   % Offshore wave length
        zSWL = june2016.forcing(i).waterLevel; % Still water level
        % Runup elevation using Stockdon
        R_st = calcRn_ST(Hs,L,b0,2) + zSWL;
        % Runup elevation using Gaussian Process
        index = i*dt-(dt-1):1:i*dt;
        R_gp = calcRn_GP(GP_mean,GP_sigma,index,'mean') + zSWL;
        % Runup elevation using Gaussian Process
        for k=1:draws
            R_gp_draws(:,k) = calcRn_GP(GP_mean,GP_sigma,index,'draw') + zSWL;
        end
        
        % Store
        data(i,1).zb = zb;
        data(i,1).zb_final = zb_final;
        data(i,1).dv = dv;
        data(i,1).dv_obs = dv_obs;
        data(i,1).Tp = Tp;
        data(i,1).R_st = R_st;
        data(i,1).R_gp = R_gp;
        data(i,1).R_gp_draws = R_gp_draws;
    end
    
end

%% Save Output
% save DIM_data.mat data flag