%% Calculate the longitude in Saturn's magnetic phase system based on
%% Provan et al., 2013 (JGR) for a given time/date and given local time of
%% the spacecraft (e.g. Cassini).
%%
%% Author: M. von Papen
%% Date: 24.08.2015
%%
%% Input: t     -   numeric date in **unit of days** with t=0 at year 0 (e.g. generated with datenum)
%%        lt    -   local time of spacecraft in dec.h
%%        r     -   radial distance of Cassini in Rs
    
function [Phin, Phis, k] = magphase(t,lt,r)


if nargin<2; lt=12; end

% T0=1January2004
t0=datenum(2004,1,1);
t = t-t0;


% Find correct interval parameters (see Table B1 in Provan et al., 2013)
if t >= 2025 && t < 2125 %E1-1
    Phi0n = 316.8;
    Tn    = 10.655;
    Phi0s = 266.4;
    Ts    = 10.745;
    k     = 0.9;
elseif t >= 2125 && t < 2175 %E1-1
    Phi0n = 144;
    Tn    = 10.655;
    Phi0s = 180;
    Ts    = 10.74;
    k     = 1.05;
elseif t<2025 || t > 2175
    error('There are no magphases for given times t.')
end
    
phicas  = 180-lt*15;
Phin    = Phi0n+(360*24/Tn)*(t-t0)+phicas;
Phis    = Phi0s+(360*24/Ts)*(t-t0)+phicas;
% plot(t-t0,cosd(Phi))


%% Radial phase delay (see Provan et al., 2014, eq. 2b)
if nargin>2
    Phin(r>12) = Phin(r>12) - 3 * ( r(r>12)-12 );
    Phis(r>12) = Phis(r>12) - 3 * ( r(r>12)-12 );
end


%% Modulus 360
Phin = mod(Phin,360);
Phis = mod(Phis,360);