%% Calculates PSD for Slab + 2D Turbulence

function P=slab2d_test(f,theta,v,C,si1,rho,si2,TPratio)

if nargin<8; TPratio=0; end
if nargin<7; si2=[7/3 5]; end
if nargin<6; rho=1e5; end
if nargin<5; si=[5/3 2]; end
if nargin<4; C=[1 1]; end
if nargin<3; v=600e3; end
if nargin<2; theta=[0 90]; end

C2D=C(1); CS=C(2);

if TPratio~=0
    for i=1:length(theta)
        
        kpe=2*pi*f(:)/v/sind(theta(i));
        kpa=2*pi*f(:)/v/cosd(theta(i));
        c1=rho^(si1(1)-si2(1));
        c2=(rho*cosd(theta(i))/sind(theta(i)))^(si1(2)-si2(2));
                
        %% Sort MHD range and kinetic range
        % Slab+2D, after Bieber1996 but with factor 1/2 b/c we handle P seperately
        oz=zeros(length(kpe),1); oz(kpe*rho<=1)=1;
        P2D=C2D*(v*sind(theta(i)))^(-1)/(1+si1(1))*(...
                kpe.^(-si1(1)).*oz+c1*kpe.^(-si2(1)).*~oz );
        PS=CS*(v*cosd(theta(i)))^(-1)*(...
                kpa.^(-si1(2)).*oz+c2*kpa.^(-si2(2)).*~oz );
            
        if theta(i)==0; P2D=0; end
        if theta(i)==90; PS=0; end
        
        %Pxx
        P(1,:,i)=P2D+PS;
        % Pyy=si*Pxx, after Bieber1996
        P(2,:,i)=si1(1)*P2D+PS;
        % Pzz, here only P in integrand <- insert TPratio
        % if T=P, then Pxx+Pyy=Pzz: (2-(kx^2+ky^2)/kp^2)=1=kp^2/k^2
        P(3,:,i)=(1+si1(1))*P2D/TPratio;
        % Psum
        P(4,:,i)=sum(P(1:3,:,i));
    end
else
    for i=1:length(theta)
        kpe=2*pi*f(:)/v/sind(theta(i));
        kpa=2*pi*f(:)/v/cosd(theta(i));
        c1=rho^(si1(1)-si2(1));
        c2=(rho*cosd(theta(i))/sind(theta(i)))^(si1(2)-si2(2));
        
        oz=zeros(length(kpe),1); oz(kpe*rho<=1)=1;
        P(:,i)=2*C2D*(v*sind(theta(i)))^(-1)*(...
                    kpe.^(-si1(1)).*oz+c1*kpe.^(-si2(1)).*~oz )...
                +2*CS*(v*cosd(theta(i)))^(-1)*(...
                    kpa.^(-si1(2)).*oz+c2*kpa.^(-si2(2)).*~oz );
    end
end