function Lo = calcLo(T)
%
%function Lo=calcLo(T)
% function to calculate offshore (deep water) wave length based on known
% wave period (T).
%
% Kristen Splinter, 2009.
g = 9.81;
Lo = g.*T.^2./(2.*pi);