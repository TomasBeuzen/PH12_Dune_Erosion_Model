function R = calcRn_GP(GP_mean,GP_sigma,index,type)
% INPUTS
% GP_mean  = pre-computed runup values
% GP_sigma = associated pre-computed sigma values
% index    = index of inquiry
% type     = required sample, can be "mean" or "draw"
%
% OUTPUTS
% R        = runup elevation (m)
%
% TB 2018

if strcmp(type,'mean')
    R = GP_mean(index);
elseif strcmp(type,'draw')
    R = normrnd(GP_mean(index),GP_sigma(index));
end