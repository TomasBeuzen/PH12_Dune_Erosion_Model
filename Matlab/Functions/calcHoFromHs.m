 function Ho = calcHoFromHs(Hs,T,h)
 %
 %  function Ho = calcHoFromHs(Hs,T,h)
 % function to reverse shoal  wave height to deep water equivalent.
 % Assumes no refraction. Used in R2 (Stockdon, 06) where Ho is needed.
 % Under variable bathy, using an inshore wave height and reverse shoaling
 % may be more accurate than using the offshore one registered at wave
 % gauge
 %
 % Hs = sig wave height in water depth (h)
 % T = wave period
 %
 % kristen, 09
 %
gamma = 0.78;
g = 9.81;
Cgo = 1./4.*g.*T./pi;
k = dispsol(h, 1./T);
L=2.*pi./k;
C = L(:)./T(:);
Cg = 0.5.*(1+ 4*pi.*h(:)./L(:)./sinh(4*pi.*h./L(:))).*C(:);

 Ho = Hs(:).*sqrt(Cg(:)./Cgo(:));
 
 