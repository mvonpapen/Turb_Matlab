%% Calculates Saturn's (axisymmetric) internal Magnetic Field
%% based on Burton et al., 2010, GRL
%% 
%% 
%%      INPUT:
%%              r:      radial distance to Saturn in R_S
%%              z:      vertical distance to geographic equator in R_S
%%
%%      OUTPUT:
%%              Bint:   internal magnetic field in nT
%%
%%
%% Date: 07.03.2014
%% Author: M. von Papen
%%


function [B,Br,Bt,Bp]=intbsat(r,z)

%% Parameters of model
g10=21136;
g20=1526;
g30=2219;
G10=0; %-11.6 only valid inside 4R_S
Rs=60268000;

%% Theta and Legendre polynomes
costheta=z./r;
P1=legendre(1,costheta,'sch');
P2=legendre(2,costheta,'sch');
P3=legendre(3,costheta,'sch');
P10=P1(1,:);
P20=P2(1,:);
P30=P3(1,:);
% Derivatives dP/dtheta
dP10 = (costheta.*P10 - 1) ./ sqrt(1-costheta.^2.);
dP20 = (2*costheta.*P20 - 2*P10) ./ sqrt(1-costheta.^2.);
dP30 = (3*costheta.*P30 - 3*P20) ./ sqrt(1-costheta.^2.);    


%% Model
% Eq. (1)
Br =  2*(1./r).^3*g10.*P10 ...
    + 3*(1./r).^4*g20.*P20 ...
    + 4*(1./r).^5*g30.*P30 ...
    + G10.*P10;
% Eq. (2)
Bt = -(1./r).^3*g10.*dP10 ...
    - (1./r).^4*g20.*dP20 ...
    - (1./r).^5*g30.*dP30 ...
    - G10.*dP10;
% Eq. (3)
Bp = 0;

B=sqrt(Br.^2+Bt.^2);