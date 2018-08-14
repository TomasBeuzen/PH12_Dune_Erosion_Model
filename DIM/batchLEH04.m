%This script loops through the structure and calculates zb and dune erosion
%volumes for all profiles in the structure. It outputs the time series of
%erosion and profile eovlution, as well as the final volumes

%to do: 
%1. compare predicted dv and zb with the observations 
%    (which struc. field are the final zb and dv for each profile?)
%3. do the ensembles
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

Cs = linspace(2e-3,1e-5,100); %LEH04 param
%best is Cs= 8e-4

%number of ensembles
ens=100;

%preallocate
errorGP=zeros(length(profile_indicies),length(Cs));
errorST=zeros(length(profile_indicies),length(Cs));
MaxGP = zeros(length(profile_indicies),length(Cs),ens);
MinGP = zeros(length(profile_indicies),length(Cs),ens);


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
    
    %error
     errorGP(i,j)=abs([data(profile_indicies(i)).dv_obs]-SigDuneErosionGP(end));
     errorST(i,j)=abs([data(profile_indicies(i)).dv_obs]-SigDuneErosionST(end));     
    
    %run it through 10 'draws' from GP
    [SigDuneErosionGPD,zbmGPD] = LEH04ensembles(dv,zb,R_gp_draws,T,Cs(j));
    %record the max and min.
    for k=1:ens
        MaxGP(i,j,k) = max(SigDuneErosionGPD(end,1:k));
        MinGP(i,j,k) = min(SigDuneErosionGPD(end,1:k));
    end
    
  
%     %put time series back in the structure
%     data(profile_indicies(i)).SigDuneErosionST = SigDuneErosionST;
%     data(profile_indicies(i)).zbmST = zbmST;
%     data(profile_indicies(i)).SigDuneErosionGP = SigDuneErosionGP;
%     data(profile_indicies(i)).zbmGP = zbmGP;
%     data(profile_indicies(i)).SigDuneErosionGPD = SigDuneErosionGPD;
%     data(profile_indicies(i)).zbmGPD = zbmGPD;
%     
%     %make a new column for total erosion and final zb
%     data(profile_indicies(i)).dVst = SigDuneErosionST(end);
%     data(profile_indicies(i)).zbst = zbmST(end); 
%     data(profile_indicies(i)).dVGP = SigDuneErosionGP(end);
%     data(profile_indicies(i)).zbGP = zbmGP(end);
%     
%     data(profile_indicies(i)).dVGPD = SigDuneErosionGPD(end,:)';
%     data(profile_indicies(i)).zbGPD = zbmGPD(end,:)';
    
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
plot(Cs,mean(errorST),'r.');
hold on
plot(Cs,mean(errorGP),'b.');
xlabel('Cs')
ylabel('MAE all profiles')
legend('ST','GP')

subplot(2,1,2)
plot(Cs,PercentWithin)
xlabel('Cs')
ylabel('percent within ensemble')



