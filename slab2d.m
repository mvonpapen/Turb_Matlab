%% Calculates PSD for Slab + 2D Turbulence

function [P,Slab_percentage]=slab2d(f,theta,v,C,si,TPratio)

if nargin<6; TPratio=0; end
if nargin<5; si=[5/3 2]; end
if nargin<4; C=[1 1]; end
if nargin<3; v=600e3; end
if nargin<2; theta=[0 90]; end

C2D=C(1); CS=C(2);

if TPratio~=0
    for i=1:length(theta)
        % Pxx, after Bieber1996 but with factor 1/2 b/c we handle P seperately
        P(1,:,i)=C2D*2/(1+si(1))*(v*abs(sind(theta(i))))^(si(1)-1)*(2*pi*f).^(-si(1))...
                 +CS*(v*abs(cosd(theta(i))))^(si(2)-1)*(2*pi*f).^(-si(2));
        % Pyy=si*Pxx, after Bieber1996
        P(2,:,i)=C2D*2*si(1)/(1+si(1))*(v*abs(sind(theta(i))))^(si(1)-1)*(2*pi*f).^(-si(1))...
                 +CS*(v*abs(cosd(theta(i))))^(si(2)-1)*(2*pi*f).^(-si(2));
        % Pzz, here only P in integrand <- insert TPratio
        % if T=P, then Pxx+Pyy=Pzz: (2-(kx^2+ky^2)/kp^2)=1=kp^2/k^2
        P(3,:,i)=C2D*(v*abs(sind(theta(i))))^(si(1)-1)*(2*pi*f).^(-si(1))/TPratio;
        % Psum
        P(4,:,i)=sum(P(1:3,:,i));
    end
    Slab_percentage=100-100./(P(4,:,1)./P(4,:,length(theta))+1);
else
    for i=1:length(theta)
        P(:,i)=2*C2D*(2*pi/v/abs(sind(theta(i))))^(1-si(1))*f.^(-si(1))...
               +2*CS*(2*pi/v/abs(cosd(theta(i))))^(1-si(2))*f.^(-si(2));
    end
    Slab_percentage=100-100./(P(:,1)./P(:,length(theta))+1);
end