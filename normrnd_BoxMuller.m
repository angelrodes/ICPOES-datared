function [ r ] = normrnd_BoxMuller(mu,sigma,varargin)
% NORMRND_BOXMULLER returns 
%   an array of random numbers chosen from a normal distribution with mean
%   MU and standard deviation SIGMA.
%   Normally distributed random numbers are generated using the Box-Muller transform.
%   https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
%   The purpose of this function is replacing the NORMRND function if the
%   NORMRND function is not avialable. E.g. if you are using Octave without
%   the statistics package.


% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;

% decide the size of the output array
if numel(sigma)>1 % if more than one SIGMA is given
    sizeout=size(sigma); % ignore the size given after SIGMA
elseif nargin > 2 % else if size is given after SIGMA
    sizeout=cell2mat(varargin); % use it
else % if not
    sizeout=1; % return just one random number
end


u1=rand(sizeout);
u2=rand(sizeout);
z=(-2*log(u1)).^0.5.*cos(2*pi*u2);

r=mu+sigma.*z;
end

