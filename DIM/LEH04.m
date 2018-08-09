% for a given set of points x, z, delV, and a given Runup and period timeseries 
%erode the dune

%counter for cumulative dune erosion, which is length of forcing
SigDuneErosion=0*R;

%counter for cumulative dune erosion, which is length of forcing
SigDuneErosion=0*R;

%start at the beginning of the dune toe
zb=;

%iteratively calculate the dune erosion fore the entire runup timeseries
for i=1:1:length(R)
    
    Cs = 1.4e-3; %LEH04 param
    t=60*60; %hourly timesteps
    
    %LEH model
    DuneErosion=4*Cs*(R(i)-zb)*(t/T(i));
    
    %if the runup is higher than the dune toe, then erode the dune
    if DuneErosion>0
        %for the dune erosion calculated, find the new zb.
        %1. find DuneToe at beginning fo timestep, 
        %and match that with a dunevolume previously eroded 
        
        %2. add together to find the index of cumulative volume lost, to
        %then interp1 to find the new zb
        
        if newzb>R(i)
            %if zb is greater than R(i), then set zb to R(i)
        
            %and only erode the amount of sand to get to that height
        
            %and overprint the DuneErosion number
        end
        
        %set new Zb for future runs
        
        %Count the cumulative erosion
        if i==1
            SigDuneErosion(i)=DuneErosion(i);
        else
            SigDuneErosion(i)=DuneErosion(i)+SigDuneErosion(i-1);
        end  
    end
    
end
