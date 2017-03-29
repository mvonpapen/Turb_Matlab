%% Numerically integrates P_i(f,\theta)=\int dk^3 P_i(k) \delta(2\pi f - k*v)
%% for given kx(=2pi f/v, kx is in rotated frame with kx||v) and \theta
%% Model of P_i(k) is critically balanced with k_para~L^(1/3)k_perp^(2/3) in ion
%% inertial range and k_\para~L^(1/3)\rho^(1/3)k_\perp^(1/3) in electron inertial range
%% P_i(k) as described in Cho et al., 2002 (ApJ)
%%
%% Author: Michael von Papen
%% Date: 26.06.2013

function [P,kern] = cbspec2 (f,theta,n,fun,L0,rho,v,va,B,bounds,ratioTP,rhoe,si_in,cb_in,mirror)


%% Check Input
if nargin<15; mirror=0; end
if nargin<14; cb_in=[2/3 1/3]; end
if nargin<13; si_in=[-10/3 -11/3]; end
if nargin<11; ratioTP=0; end
if nargin<10; bounds=[-10 1]; end
if nargin<9; B=1e-9; end
if nargin<8; va=60e3; end
if nargin<7; v=600e3; end
if nargin<6; rho=1e5; end
if nargin<5; L0=1e9; end
if nargin<4; fun=1; end
if nargin<3; n=500; end
if nargin<2; theta=[0 90]; end
if nargin<12; rhoe=rho/42.85; end


%% Wave vector in rotated coordinate system
kix=2*pi*f/v;

%% Basic Parameter
si=si_in(1); %-10/3->k^{-5/3}
si2=si_in(2); %-11/3->k^{-7/3}
cb=cb_in(1); % 2/3->alfven
cb2=cb_in(2); %1/3->KAW


%% K-space gridpoints
% Set boundaries a little bit wider than kmin,kmax for numerical reasons.
% Later everything outside [kmin,kmax] will be disregarded
kmin=bounds(1);
kmax=bounds(2);

    
ky1=10.^[kmin+(0:n-2)*(kmax-kmin)/(floor(n)-1), kmax];
dky=[ky1(2:end) 2*10^kmax]-ky1;
ky=ky1(ones(1,2*n),:)'; % <=> ky=repmat(ky',1,2*n);
dky=dky(ones(1,2*n),:)';

%% Set output variable
if ratioTP ~=0
    P=zeros(4,length(kix),length(theta));
else
    P=zeros(length(kix),length(theta));
end


% %% Ion-cyclotron frequency cut-off
% % All k_para ~ w_ic/V_A are subject to ion-cyclotron damping
% % Thus, parallel scales cannot reach k_para > w_ic/V_A
% if nargin<6;
    kz_max=kmax;
% else
%     kz_max=min([kmax log10(1.6e-19*B/1.67e-27/va)]);
% end

%% Begin with loop over theta
for k=1:length(theta)
   
    
    %% Begin loop over frequency
    for i=1:length(kix)
        
        % Find kiz, so that kz~0
        kiz=[-ky1(end:-1:1) ky1];
        kz=kix(i)*cosd(theta(k))+kiz.*sind(theta(k));
        [i1,i2]=min(abs(kz));
        %nz log verteilt auf pos UND neg Achse
        if theta(k)==0;
            kiz=[-ky1(end:-1:1) ky1];
        else
            kiz=[-ky1(end:-1:1) ky1]+kiz(i2);    
        end
        
        dkiz=[kiz(2:end) 2*10^kmax]-kiz;
        kiz=kiz(ones(n,1),:); % <=> kiz=repmat(kiz',1,2*n);
        dkiz=dkiz(ones(n,1),:);

        %% Calculate PSD at z
        kx=kix(i)*sind(theta(k))-kiz.*cosd(theta(k)); %=kx in unrotated system
        %ky=kiy in unrotated system
        kz=kix(i)*cosd(theta(k))+kiz.*sind(theta(k)); %=kz in unrotated system

        kern=zeros(n,2*n);
            
        kp2=ky.^2+kx.^2; % k_perp^2
        kabs2=kp2+kz.^2; % |k|
%         cbetaT=cbetaT

        %% Equations written in unprimed coordinates for the sake of
        %% brevity, but integration is done over primed variables,
        %% which is why dkiy and dkiz is used.
        
        for L=L0
            %% Single components Alfven cascade
            i1=find(kp2 <= 1/rho^2 & kp2 > 1/L^2 & abs(kz) <= 10^kz_max ...
                & kabs2 <= 10^(2*kmax));
            switch fun
                case 1 %'exp'
                    kern( i1 ) = kp2(i1).^(si/2)...
                        .*exp(-L^(1-cb).*abs(kz(i1))./kp2(i1).^(cb/2))...
                        .*dky(i1).*dkiz(i1);
                case 2 %'expdamp' for Ti=Te => rho_e=sqrt(me/mi)*rho_i
                    kern( i1 ) = kern(i1) + kp2(i1).^(si/2)...
                        .*exp(-L^(1-cb).*abs(kz(i1)./kp2(i1).^(cb/2))...
                              -sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1);                    
                case 3 %'gauss'
                    kern( i1 ) = kern(i1) + kp2(i1).^(si/2)...
                        .*exp(-(L^(1-cb)*abs(kz(i1))./kp2(i1).^(cb/2)-1).^2)...
                        .*dky(i1).*dkiz(i1)/sqrt(pi);
                case 4 %'gaussdamp'
                    kern( i1 ) = kern(i1) + kp2(i1).^(si/2)...
                        .*exp(-(L^(1-cb)*abs(kz(i1))./kp2(i1).^(cb/2)-1).^2 ...
                              -sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1)/sqrt(pi); 
                case 5 %'heavi'
                    i2=find(L^(1-cb)*abs(kz(i1)./kp2(i1).^(cb/2))<=1);
                    kern( i1(i2) ) = kern(i1(i2)) + kp2(i1(i2)).^(si/2).*dky(i1(i2)).*dkiz(i1(i2));
                case 6 %'delta'
                    [tmp,i2]=min((L^(1-cb)*abs(kz(i1))-kp2(i1).^(cb/2)).^2);
                    kern( i1(i2) ) = kern(i1(i2)) + L^(1-cb).*kp2(i1(i2)).^(si/2)...
                        .*dky(i1(i2)).*dkiz(i1(i2));
                case 7 %'expisodamp'
                    kern( i1 ) = kern(i1) + kp2(i1).^(si/2)...
                        .*exp(-sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1);
            end


            %% Single components KAW cascade
            i1=find(kp2 > 1/rho^2 & kp2 > 1/L^2 & abs(kz) <= 10^kz_max ...
                & kabs2 <= 10^(2*kmax) & kp2 <= 1/rhoe^2);
            switch fun
                case 1 %'exp'
                    kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                        .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2) )...
                        .*dky(i1).*dkiz(i1);
                case 2 %'expdamp' for Ti=Te => rho_e=sqrt(me/mi)*rho_i
                    kern( i1 ) = kern(i1) + rho^(si2-si).*kp2(i1).^(si2/2)...
                        .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)...
                               -sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1);                  
                case 3 %'gauss'
                    kern( i1 ) = kern(i1) + rho^(si2-si).*kp2(i1).^(si2/2)...
                        .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)-1).^2)...
                        .*dky(i1).*dkiz(i1)/sqrt(pi);
                case 4 %'gaussdamp'
                    kern( i1 ) = kern(i1) + rho^(si2-si).*kp2(i1).^(si2/2)...
                        .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)-1).^2 ...
                               -sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1)/sqrt(pi);                
                case 5 %'heavi'
                    i2=find(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)<=1);
                    kern( i1(i2) ) = kern(i1(i2)) + rho^(si2-si).*kp2(i1(i2)).^(si2/2)...
                        .*dky(i1(i2)).*dkiz(i1(i2));
                case 6 %'delta'
                    [tmp,i2]=min((L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))-kp2(i1).^(cb2/2)).^2);
                    kern( i1(i2) ) = kern(i1(i2)) + L^(1-cb)*rho^(si2-si).*kp2(i1(i2)).^(si2/2)...
                        .*dky(i1(i2)).*dkiz(i1(i2));
                case 7 %'expisodamp'
                    kern( i1 ) = kern(i1) + rho^(si2-si).*kp2(i1).^(si2/2)...
                        .*exp(-sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1);                    
            end
        
        %% Single components of cascade at electron scales
        i1=find(kp2 > 1/rhoe^2 & abs(kz) <= 10^kz_max ...
            & kabs2 <= 10^(2*kmax));
        switch fun
            case 1 %'exp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2 )...
                    .*dky(i1).*dkiz(i1);
            case 2 %'expdamp' for Ti=Te => rho_e=sqrt(me/mi)*rho_i
                kern( i1 ) = kern(i1) + rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2...
                           -sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1);                  
            case 3 %'gauss'
                kern( i1 ) = kern(i1) + rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2-1).^2)...
                    .*dky(i1).*dkiz(i1)/sqrt(pi);
            case 4 %'gaussdamp'
                kern( i1 ) = kern(i1) + rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2-1).^2 ...
                           -sqrt(kp2(i1))*rhoe).*dky(i1).*dkiz(i1)/sqrt(pi);                
            case 5 %'heavi'
                i2=find(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2<=1);
                kern( i1(i2) ) = kern(i1(i2)) + rho^(si2-si).*kp2(i1(i2)).^(si2/2)...
                    .*dky(i1(i2)).*dkiz(i1(i2));
            case 6 %'delta'
                [tmp,i2]=min((L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))-rhoe^cb2).^2);
                kern( i1(i2) ) = kern(i1(i2)) + L^(1-cb)*rho^(si2-si).*kp2(i1(i2)).^(si/2)...
                    .*dky(i1(i2)).*dkiz(i1(i2));
        end        
        
        end
        
        %% Full version with Toroidal and Poloidal parts
        % Sum up to get power for one ky value
        if ratioTP ~= 0
            Tor=2*ratioTP*kern./kp2;
            Pkern=kern;
            
            % Add mirror modes to Pkern %WIP
            if mirror~=0
                sigma=0.2; sigmaz=1e-7;
                i1=find(kz.^2./kp2<1e-2); %populate Pkern up to certain kz
                % Energy distribution of mirror mode modelled as gaussian
                % in logspace for kp*rho and normal gaussian for kz
                Pkern(i1)=Pkern(i1)+mirror/(2*pi*sqrt(sigma*sigmaz))...
                    *exp(-log10(sqrt(kp2(i1))*rho).^2/2/sigma^2 ...
                         -kz(i1).^2/2/sigmaz^2).*dky(i1).*dkiz(i1);
            end
            
            Pol=2*(1-ratioTP)*Pkern./kabs2;
            P(1,i,k)=sum(sum( ky.^2.*Tor ...
                       +(kx.*kz).^2./kp2.*Pol ));
            P(2,i,k)=sum(sum( kx.^2.*Tor ...
                       +(ky.*kz).^2./kp2.*Pol ));
            P(3,i,k)=sum(sum( kp2.*Pol ));
            P(4,i,k)=sum(P(1:3,i,k));
        else
            P(i,k)=sum(sum(2*kern));
        end
    end
end

P=1/v/L^(1-cb).*P;