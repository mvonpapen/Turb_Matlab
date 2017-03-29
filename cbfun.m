%%
%% Author: Michael von Papen
%% Date: 26.06.2013

function Y = cbfun(ky,kz,f)

si=-10/3;
cb=2/3;
L=1e9;
kx=2*pi*f/600e3;
psi=90;

Y = ((kx.*sind(psi)-kz.*cosd(psi)).^2+ky.^2).^(si/2).*exp(-L^(1-cb).*sqrt((kx.*cosd(psi)+kz.*sind(psi)).^2)./((kx.*sind(psi)-kz.*cosd(psi)).^2+ky.^2).^(cb/2));