function [xd, zd] = findDuneToe(x,z)
% 
%function [xd zd] = findDuneToe(x,z);
%
% function to locate dune toe location and elevation to be used to find
% beach face slope and in R2 calculations.  Basic stuff as a first attack,
% based on derivatives.
% ASSUMES you're starting at the shore
% Inputs:
% x - cross-shore locations
% z - bed elevation, z+ is above AHD (dry beach)
% Outputs:
% xd,zd - location and elevation of dune toe
%
% Kristen Splinter, 2009.
zFull=z;
xFull=x;
clear z
dx = diff(x(1:2));
for ii=1:size(zFull,1)
    z = zFull(ii,:);
    % keep stuff above AHD =0
    clear dry high
    [val high] = max(z);
    z=z(high:end);
    x=x(high:end); clear high
    %high = min(high,4);
    high = find(z>=min(max(z),3.0),1,'last');
    dry = find(z<0.1,1,'first');
    figure(1)
    clf
       if(isempty(dry)==0)
        z=z(high:dry);
        x=x(high:dry);
        plot(x,z,'k'), hold on
        %p=polyfit(x,z,6);
        %z=polyval(p,x);
        %hold on, plot(x,z,'g')
        dzdx = [0 diff(z)./dx];  %slope
        d2zdx2 = [0 diff(z,2)./(dx.^2) 0]; %gradient of slope
        %smooth it out a bit
        filt = hann(3)./sum(hann(3));
        dz2filt = conv2(d2zdx2(:),filt(:),'same');
        %ch = find(diff(sign(d2zdx2))==-2)
        %d3zdx3 = gradient(dz2filt,dx);
        neg = find(dzdx<0);
    %large = find(d2zdx2(neg)>0.005);
     %   if(isempty(neg)==0)
    %    if (large(1)==1)
    %        large=large(2:end);
    %    end
            [foo,ixd]=max(d2zdx2(neg));
            if(isempty(ixd)==0)
                xd(ii) = x(neg(ixd));
                zd(ii) = z(neg(ixd));
            else
                xd(ii) = NaN;
                zd(ii) = NaN;
            end
        else
            xd(ii) = NaN;
            zd(ii) = NaN;
    end
        plot(xd,zd,'ks'), hold on, ylim([-1 4])
        pause(1)
    end
end