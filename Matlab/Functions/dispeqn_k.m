function y = dispeqn_k(k,h,f)
% FUNCTION y = dispeqn(k,h,f)
% radial k, h, and f are all assumed scalars.

sig = 2*pi*f;
sigsq=sig.^2;
g = 9.80665;
y = rms((g*k*tanh(k*h))./sigsq-1);


