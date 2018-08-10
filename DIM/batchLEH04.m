%This script loops through the structure and calculates zb and dune erosion
%volumes for all profiles in the structure. It outputs the time series of
%erosion and profile eovlution, as well as the final volumes

%to do: 
%1. compare predicted dv and zb with the observations 
%    (which struc. field are the final zb and dv for each profile?)
%2. calibration
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

%find the indicies of structure that are non-empty
profile_indicies = find(~cellfun(@isempty,{data.zb}))';

%loop through those indicies
for i = 1:1:length(profile_indicies)
    
    %pull all relevant data out of the structure
    dv = data(profile_indicies(i)).dv;
    zb = data(profile_indicies(i)).zb;
    T = data(profile_indicies(i)).Tp;
    R_st = data(profile_indicies(i)).R_st;
    R_gp = data(profile_indicies(i)).R_gp;
    
    %run it through the 2 models
    [SigDuneErosionGP,zbmGP] = LEH04(dv,zb,R_gp,T);
    [SigDuneErosionST,zbmST] = LEH04(dv,zb,R_st,T);
    
    %put time series back in the structure
    data(profile_indicies(i)).SigDuneErosionGP = SigDuneErosionGP;
    data(profile_indicies(i)).zbmGP = zbmGP;
    data(profile_indicies(i)).SigDuneErosionST = SigDuneErosionST;
    data(profile_indicies(i)).zbmST = zbmST;
    
    %make a new column for total erosion and final zb
    data(profile_indicies(i)).dVGP = SigDuneErosionGP(end);
    data(profile_indicies(i)).dVst = SigDuneErosionST(end);
    data(profile_indicies(i)).zbGP = zbmGP(end);
    data(profile_indicies(i)).zbst = zbmST(end); 
end

%plot the comparison on dune erosion metrics
figure 
subplot(1,2,1);
plot([data.dVGP],[data.dVst],'.')
xlabel('dV GP')
ylabel('dV st')
hold on
plot(0:40,0:40,'k-')
axis([0 40 0 40])

subplot(1,2,2)
plot([data.zbGP],[data.zbst],'.')
xlabel('zb GP')
ylabel('zb st')
hold on
plot(2:5,2:5,'k-')
axis([2 5 2 5])


%there are nans, so these are not correct... and some errors are large.
MSE_zb_GP= mean(nansum(([data.zb_final]-[data.zbGP]).^2))
MSE_zb_st= mean(nansum(([data.zb_final]-[data.zbst]).^2))

MAE_dv_GP= mean(abs(([data.dv_obs]-[data.dVGP])))
MAE_dv_st= mean(abs(([data.dv_obs]-[data.dVst])))
