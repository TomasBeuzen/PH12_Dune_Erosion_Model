%Code to plot a model run for a given profile

clear all
close all
clc

%load data
load DIM_data.mat
data=data_subset;

%profile #
p=72;
%number of ensembles
ens=10;

%LEH parameter
Cs=8e-4;

%retrieve the data
%run it through GP
[SigDuneErosionGP,zbmGP] = LEH04ensembles(data_subset(p).dv,data_subset(p).zb,data_subset(p).R_gp,data_subset(p).Tp,Cs);
%run it through 'draws' from GP
[SigDuneErosionGPD,zbmGPD] = LEH04ensembles(data_subset(p).dv,data_subset(p).zb,data_subset(p).R_gp_draws(:,1:ens),data_subset(p).Tp,Cs);

%erosion and runup timeseries plots
figure
subplot(2,1,1)
plot(data(p).R_gp_draws(:,1:ens))
hold on
plot(data(p).R_gp,'k','LineWidth',5)
plot(108,data(p).zb_final,'*')
plot(zbmGPD)
plot(zbmGP,'k','LineWidth',2)
xlabel('time (hours)')
ylabel('elevation of runup and dune toe height')

subplot(2,1,2)
plot(SigDuneErosionGPD)
hold on
plot(SigDuneErosionGP,'k','LineWidth',5)
plot(108,data(p).dv_obs,'*')
xlabel('time (hours)')
ylabel('total delV')


