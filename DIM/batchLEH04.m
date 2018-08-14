%This script loops through the structure and calculates zb and dune erosion
%volumes for all profiles in the structure. It outputs the time series of
%erosion and profile eovlution, as well as the final volumes
%
%"It is easier to write a new code than to understand an old one"
%-John von Neumann to Marston Morse, 1952
%
%EBG Aug. 9 2018

clear all
close all
clc

%load data
load DIM_data.mat
data=data_subset;
profile_indicies = 1:508;

Cs = logspace(-5,-2,50); %LEH04 param
%best is Cs= 8e-4

%number of ensembles
ens=100;

%preallocate
errorGP=zeros(length(profile_indicies),length(Cs));
errorST=zeros(length(profile_indicies),length(Cs));
errorGPmean=zeros(length(profile_indicies),length(Cs));
errorGPmode=zeros(length(profile_indicies),length(Cs));
errorGPmedian=zeros(length(profile_indicies),length(Cs));
MaxGP = zeros(length(profile_indicies),length(Cs),ens);
MinGP = zeros(length(profile_indicies),length(Cs),ens);
MeanGP = zeros(length(profile_indicies),length(Cs),ens);
ModeGP = zeros(length(profile_indicies),length(Cs),ens);
MedianGP = zeros(length(profile_indicies),length(Cs),ens);

%loop through those indicies
for i = 1:1:length(profile_indicies)
    i
    for j = 1:length(Cs) %LEH04 param
    %pull all relevant data out of the structure
    dv = data(profile_indicies(i)).dv;
    zb = data(profile_indicies(i)).zb;
    T = data(profile_indicies(i)).Tp;
    R_st = data(profile_indicies(i)).R_st;
    R_gp = data(profile_indicies(i)).R_gp;
    R_gp_draws = data(profile_indicies(i)).R_gp_draws;
    
    %run it through the St model
    [SigDuneErosionST,zbmST] = LEH04ensembles(dv,zb,R_st,T,Cs(j));
    
    %run it through GP
    [SigDuneErosionGP,zbmGP] = LEH04ensembles(dv,zb,R_gp,T,Cs(j));

    %run it through 'ens' number of 'draws' from GP
    [SigDuneErosionGPD,zbmGPD] = LEH04ensembles(dv,zb,R_gp_draws,T,Cs(j));
    
    %record the max and min from GP draws
    for k=1:ens
        MaxGP(i,j,k) = max(SigDuneErosionGPD(end,1:k));
        MinGP(i,j,k) = min(SigDuneErosionGPD(end,1:k));
        MeanGP(i,j,k) = mean(SigDuneErosionGPD(end,1:k));
        ModeGP(i,j,k) = mode(SigDuneErosionGPD(end,1:k));
        MedianGP(i,j,k) = median(SigDuneErosionGPD(end,1:k));
    end
          
    %record error
    errorST(i,j)=abs([data(profile_indicies(i)).dv_obs]-SigDuneErosionST(end));     
    errorGP(i,j)=abs([data(profile_indicies(i)).dv_obs]-SigDuneErosionGP(end));
    
    %GP ensemble mean is done for max number of ensembles.
    errorGPmean(i,j)=abs([data(profile_indicies(i)).dv_obs]-MeanGP(i,j,ens));
    errorGPmode(i,j)=abs([data(profile_indicies(i)).dv_obs]-ModeGP(i,j,ens));
    errorGPmedian(i,j)=abs([data(profile_indicies(i)).dv_obs]-MedianGP(i,j,ens));
    
    
    %Record output in the structure
     data(profile_indicies(i)).SigDuneErosionST = SigDuneErosionST;
     data(profile_indicies(i)).zbmST = zbmST;
     data(profile_indicies(i)).SigDuneErosionGP = SigDuneErosionGP;
     data(profile_indicies(i)).zbmGP = zbmGP;
     data(profile_indicies(i)).SigDuneErosionGPD = SigDuneErosionGPD;
     data(profile_indicies(i)).zbmGPD = zbmGPD;
         
    end
end

save data.mat
%save('DIM_data_model_ensemble_Cs','data')

%capture stats
PercentWithin=zeros(length(Cs),ens);

obs=[data.dv_obs];
for m = 1:length(Cs)
    for n = 1:ens
    within=obs(obs<=MaxGP(:,m,n)' & obs>=MinGP(:,m,n)');
    PercentWithin(m,n)=100*(length(within)/508);
    end
end


figure
subplot(2,1,1)
semilogx(Cs,mean(errorST),'r.');
hold on
semilogx(Cs,mean(errorGP),'b.');
semilogx(Cs,mean(errorGPmean),'k.');
semilogx(Cs,mean(errorGPmedian),'k--');
semilogx(Cs,mean(errorGPmode),'k.-');
xlabel('Cs')
ylabel('MAE (all profiles)')
legend('ST','GP', 'E. Mean', 'E. Median', 'E. Mode')

subplot(2,1,2)
semilogx(Cs,PercentWithin)
xlabel('Cs')
ylabel('percent within ensemble')


%%%%%%%%%%%
%find index of value Cs=7e-4; in the Cs vector... 
%Hard coded for now

j=31;

Mean_GP=squeeze(MeanGP(:,j,:));
Mode_GP=squeeze(ModeGP(:,j,:));
Median_GP=squeeze(MedianGP(:,j,:));

a=[data.dv_obs]';
b=repmat(a,1,100);

MAE_GPmean = mean(abs(b-Mean_GP));
MAE_GPmode = mean(abs(b-Mode_GP));
MAE_GPmedian = mean(abs(b-Median_GP));

%CAPTURE STATS
Capture=squeeze(PercentWithin(j,:));

figure
subplot(2,1,1)
plot(MAE_GPmean)
hold on
plot(MAE_GPmode)
plot(MAE_GPmedian)
xlabel('# of ensembles')
ylabel('MAE')
legend('Mean','Mode','Median')

subplot(2,1,2)
plot(Capture)
xlabel('# of ensembles')
ylabel('Captured %')




