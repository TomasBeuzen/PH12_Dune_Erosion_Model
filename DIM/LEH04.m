function [SigDuneErosion,zbm] = LEH04(dv,zb,R,T)
%for a given set of points x, z, delV, and a given Runup and period timeseries
%erode the dune

%counter for cumulative dune erosion, which is length of forcing
SigDuneErosion=0*R;

%counter for dune foot elevation, which is length of forcing
zbm=0*R;

%iteratively calculate the dune erosion fore the entire runup timeseries
for i=1:1:length(R)
    
    Cs = 1.4e-3; %LEH04 param
    t=60*60; %hourly timesteps
        
    %LEH model
    if i==1 %use the initial value from the curve, which is the pre storm
        DuneErosion=4*Cs*(R(i)-zb(1))*(t/T(i));
    else %use the model value from previosu time step
        DuneErosion=4*Cs*(R(i)-zbm(i-1))*(t/T(i));
    end
    
    %if the runup is higher than the dune toe, then erode the dune
    if DuneErosion>0
        %for the dune erosion calculated, find the new zb.
        %1. find cumulative dune erosion
        if i==1
            SigDuneErosion(i) = DuneErosion;
        else
            SigDuneErosion(i) = DuneErosion + SigDuneErosion(i-1);
        end
        
        %2. Use to find the new zb, using inteprp1
        zbm(i) = interp1(dv,zb,SigDuneErosion(i));
        
        %if the new zb is higher than R(i), but it started out lower than
        %(Ri), then set new zb to R(i); i.e., no over-erosion
        if zbp < R(i) && zbm(i) > R(i)
            zbm(i) = R(i);
            %and then only erode the amount of sand to get to that
            %readjusted height and overprint the DuneErosion number
            SigDuneErosion(i) = interp1(zb,dv,newzbm);
        end
    end
end
end
