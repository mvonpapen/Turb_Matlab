function [mu, sigma] = meanangle(in,dim,sens);

% MEANANGLE will calculate the mean of a set of angles (in degrees) based
% on polar considerations.
%
% Usage: [mu, sigma] = meanangle(in,dim)
%
% in is a vector or matrix of angles (in degrees)
% out is the mean of these angles mu along the dimension dim and its
% standard deviation sigma
%
% If dim is not specified, the first non-singleton dimension is used.
%
% A sensitivity factor is used to determine oppositeness, and is how close
% the mean of the complex representations of the angles can be to zero
% before being called zero.  For nearly all cases, this parameter is fine
% at its default (1e-12), but it can be readjusted as a third parameter if
% necessary:
%
% [mu, sigma] = meanangle(in,dim,sensitivity)
%
% Written by J.A. Dunne, 10-20-05
% added sigma by M. von Papen, 2015

if nargin<3
    sens = 1e-12;
end

if nargin<2
    ind = min(find(size(in)>1));
    if isempty(ind)
        %This is a scalar
        mu = in;
        return
    end
    dim = ind;
end

in2 = in * pi/180;

in2 = exp(1i*in2);
mid = mean(in2,dim);
mu  = atan2(imag(mid),real(mid))*180/pi;
mu(abs(mid)<sens) = nan;

% % wrap angles so that mu is in the middle and compute standard deviation
% in    = mod(in - mu + 180, 360);
% sigma = std(in);