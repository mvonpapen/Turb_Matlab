%% Numerically integrates P_i(f,\theta)=\int dk^3 P_i(k) \delta(2\pi f - k*v)
%% for given kx(=2pi f/v, kx is in rotated frame with kx||v) and \theta
%% Model of P_i(k) is critically balanced with k_para~L^(1/3)k_perp^(2/3) in ion
%% inertial range and k_\para~L^(1/3)\rho^(1/3)k_\perp^(1/3) in electron inertial range
%% P_i(k) as described in Cho et al., 2002 (ApJ)
%%
%%
%%  Input:
%%      f       Frequncy in Hz
%%      theta   Field-to-flow angle in degrees
%%      n       Number of nodes in one dimension
%%      fun     function to evaluate ( 1=exp, 2=expdamp, ... )
%%      L       Outer scale in m
%%      rho     Gyro radius (or else) in m
%%      v       Plasma bulk velocity in m/s
%%      va      Alfven velocity in m/s
%%      B       Magnetic field strength in nT
%%      bounds  Boundary in log10 for integration ( bounds = [kmin kmax] )
%%      ratioTP Ratio between toroidal T and poloidal P fluctuations
%%      rhoe    Electron gyro radius (or else)
%%      si_in   Spectral index in [MHD kinetic] range
%%      cb_in   Critical balance exponent in [MHD kinetic] range
%%
%%
%% Author: Michael von Papen
%% Date: 26.06.2013

function [P,kern,Pkern] = cbspec (f,theta,fun,varargin)


%% Check Input
args=struct('B',5,...
            'bounds',[-10 -2],...
            'cb_in',[2/3 1/3],...
            'isopol',0,...
            'L',1e9,...
            'mirror',0,...
            'n',500,...
            'ratioTP',0,...
            'rho',1e5,...
            'rhoe',NaN,...
            'si_in',[-10/3 -11/3],...
            'v',6e5,...
            'va',6e4);
args=parseArgs(varargin,args);
if isnan(args.rhoe); args.rhoe=args.rho/42.85; end
L=args.L;
rho=args.rho;
rhoe=args.rhoe;
va=args.va;


%% Wave vector in rotated coordinate system
kix=2*pi*f/args.v;

%% Basic Parameter
si=args.si_in(1); %-10/3->k^{-5/3}
si2=args.si_in(2); %-11/3->k^{-7/3}
cb=args.cb_in(1); % 2/3->alfven
cb2=args.cb_in(2); %1/3->KAW


%% K-space gridpoints
% Set boundaries a little bit wider than kmin,kmax for numerical reasons.
% Later everything outside [kmin,kmax] will be disregarded
kmin=args.bounds(1);
kmax=args.bounds(2);

    
ky=10.^[kmin+(0:args.n-1)*(kmax-kmin)/args.n];
dky=[ky(2:end) 2*10^kmax]-ky;
%nz log verteilt auf pos UND neg Achse
kiz=repmat([-ky(end:-1:1) ky],args.n,1);
dkiz=repmat([dky(end:-1:1) dky],args.n,1);
ky=ky(ones(1,2*args.n),:)'; % <=> ky=repmat(ky',1,2*n);
dky=dky(ones(1,2*args.n),:)';

% %% Ion-cyclotron frequency cut-off
% % All k_para ~ w_ic/V_A are subject to ion-cyclotron damping
% % Thus, parallel scales cannot reach k_para >> w_ic/V_A
wic=1.6e-19*args.B*1e-9/1.67e-27; %wic for protons

%% Set output variable
if args.ratioTP ~=0
    if args.isopol==1
        P=zeros(6,length(kix),length(theta));
    else
        P=zeros(4,length(kix),length(theta));
    end
else
    P=zeros(length(kix),length(theta));
end

%% Begin with loop over theta
for k=1:length(theta)
%     fprintf('theta=%3.1f°\n',theta(k))
    
    
    %% Begin loop over frequency
    for i=1:length(kix)
        
        %% Calculate PSD at z
        %% CAUTION: Setting kiz=0 on x-axis makes problems when determining
        %% spectral anisotropy on MHD scales for small outer scales and/or
        %% small angles
        if theta(k)==0
            kx=kix(i)*sind(theta(k))-kiz.*cosd(theta(k)); %=kx in unrotated system
            %ky=kiy in unrotated system
            kz=kix(i)*cosd(theta(k))+kiz.*sind(theta(k)); %=kz in unrotated system
        else
            kx=kix(i)*sind(theta(k))-kiz.*cosd(theta(k))+kix(i).*cosd(theta(k))^2/sind(theta(k)); %=kx in unrotated system with kiz=0 on x-axis
            %ky=kiy in unrotated system
            kz=kiz.*sind(theta(k)); %=kz in unrotated system with kiz=0 on x-axis
        end
        

        kern=zeros(args.n,2*args.n);
            
        kp2=ky.^2+kx.^2; % k_perp^2
        kabs2=kp2+kz.^2; % |k|
        
%         %% Get damping factors from Anne, Tmax
%         load('/afs/geo/usr/vonpapen/PhD/publication/cbmodel/anne_data/anne.mat')
%         efi=interp1(kp,damp_factor_chen',sqrt(kp2));


        %% Equations written in unprimed coordinates for the sake of
        %% brevity, but integration is done over primed variables,
        %% which is why dkiy and dkiz is used.

        %% Single components Alfven cascade
        % take out '& kabs2 > 1/L^2' when checking for fmax or showing L-dependence
        i1=find(kabs2 > 1/L^2 & kp2 <= 1/rho^2);
        switch fun
            case 1 %'exp'
                kern( i1 ) = kp2(i1).^(si/2)...
                    .*exp(-L^(1-cb).*abs(kz(i1))./kp2(i1).^(cb/2) )...
                    .*dky(i1).*dkiz(i1);
            case 2 %'expdamp'
                kern( i1 ) = kp2(i1).^(si/2)...
                    .*exp(-L^(1-cb).*abs(kz(i1))./kp2(i1).^(cb/2)...
                          -sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
            case 3 %'gauss'
                kern( i1 ) = kp2(i1).^(si/2)...
                    .*exp(-(L^(1-cb)*abs(kz(i1))./kp2(i1).^(cb/2)-1).^2)...
                    .*dky(i1).*dkiz(i1)/sqrt(pi);
            case 4 %'gaussdamp'
                kern( i1 ) = kp2(i1).^(si/2)...
                    .*exp(-(L^(1-cb)*abs(kz(i1))./kp2(i1).^(cb/2)-1).^2 ...
                          -sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1)/sqrt(pi); 
            case 5 %'heavi'
                i2=find(L^(1-cb)*abs(kz(i1)./kp2(i1).^(cb/2))<=1);
                kern( i1(i2) ) = kp2(i1(i2)).^(si/2).*dky(i1(i2)).*dkiz(i1(i2));
            case 6 %'delta'
                [tmp,i2]=min((L^(1-cb)*abs(kz(i1))-kp2(i1).^(cb/2)).^2);
                kern( i1(i2) ) = L^(1-cb).*kp2(i1(i2)).^(si/2)...
                    .*dky(i1(i2)).*dkiz(i1(i2));
            case 7 %'expisodamp'
                kern( i1 ) = kp2(i1).^(si/2)...
                    .*exp(-sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
            case 8 %'exp anne damp'
                kern( i1 ) = kp2(i1).^(si/2)...
                    .*exp(-L^(1-cb).*abs(kz(i1))./kp2(i1).^(cb/2)-abs(kz(i1))*va/wic)...
                    .*dky(i1).*dkiz(i1).*efi(i1);
            case 9 %'exp' damped only by IC damping
                kern( i1 ) = kp2(i1).^(si/2)...
                    .*exp(-L^(1-cb).*abs(kz(i1))./kp2(i1).^(cb/2)...
                    -abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
        end


        %% Single components KAW cascade
        i1=find(kabs2 > 1/L^2 & kp2 > 1/rho^2 & kp2 <= 1/rhoe^2);
        switch fun
            case 1 %'exp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2) )...
                    .*dky(i1).*dkiz(i1);
            case 2 %'expdamp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)...
                           -sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
            case 3 %'gauss'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)-1).^2)...
                    .*dky(i1).*dkiz(i1)/sqrt(pi);
            case 4 %'gaussdamp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)-1).^2 ...
                           -sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1)/sqrt(pi);                
            case 5 %'heavi'
                i2=find(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)<=1);
                kern( i1(i2) ) = rho^(si2-si).*kp2(i1(i2)).^(si2/2)...
                    .*dky(i1(i2)).*dkiz(i1(i2));
            case 6 %'delta'
                [tmp,i2]=min((L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))-kp2(i1).^(cb2/2)).^2);
                kern( i1(i2) ) = L^(1-cb)*rho^(si2-si).*kp2(i1(i2)).^(si2/2)...
                    .*dky(i1(i2)).*dkiz(i1(i2));
            case 7 %'expisodamp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp(-sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
            case 8 %'exp anne damp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)-abs(kz(i1))*va/wic )...
                    .*dky(i1).*dkiz(i1).*efi(i1);
            case 9 %'exp' damped only by IC damping
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))./kp2(i1).^(cb2/2)...
                    -abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
        end

        %% Single components of cascade at electron scales
        i1=find(kabs2 > 1/L^2 & kp2 > 1/rhoe^2);
        switch fun
            case 1 %'exp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2 )...
                    .*dky(i1).*dkiz(i1);
            case 2 %'expdamp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2...
                           -sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
            case 3 %'gauss'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2-1).^2)...
                    .*dky(i1).*dkiz(i1)/sqrt(pi);
            case 4 %'gaussdamp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2-1).^2 ...
                           -sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1)/sqrt(pi);                
            case 5 %'heavi'
                i2=find(L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2<=1);
                kern( i1(i2) ) = rho^(si2-si).*kp2(i1(i2)).^(si2/2)...
                    .*dky(i1(i2)).*dkiz(i1(i2));
            case 6 %'delta'
                [tmp,i2]=min((L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))-rhoe^cb2).^2);
                kern( i1(i2) ) = L^(1-cb)*rho^(si2-si).*kp2(i1(i2)).^(si/2)...
                    .*dky(i1(i2)).*dkiz(i1(i2));
            case 7 %'expisodamp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp(-sqrt(kp2(i1))*rhoe-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
            case 8 %'exp'
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2-abs(kz(i1))*va/wic )...
                    .*dky(i1).*dkiz(i1).*efi(i1); 
            case 9 %'exp' damped only by IC damping
                kern( i1 ) = rho^(si2-si).*kp2(i1).^(si2/2)...
                    .*exp( -L^(1-cb)*rho^(cb-cb2)*abs(kz(i1))*rhoe^cb2...
                    -abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
        end        
        
        
        %% Full version with Toroidal and Poloidal parts
        % Sum up to get power for one ky value
        if args.ratioTP ~= 0
            Tor=kern./kp2;
            
            % Add mirror modes to Pkern %WIP
            if args.mirror~=0
                sigma=0.2; sigmaz=1e-6;
                i1=find(kz.^2./kp2<1e-2); %populate Pkern up to certain kz
                % Energy distribution of mirror mode modelled as gaussian
                % in logspace for kp*rho and normal gaussian for kz
                kern(i1)=kern(i1)+args.mirror/(2*pi*sqrt(sigma*sigmaz))...
                    *exp(-log10(sqrt(kp2(i1))*rho).^2/2/sigma^2 ...
                         -kz(i1).^2/2/sigmaz^2).*dky(i1).*dkiz(i1);
            end
            
            
            %% Compose Poloidal Spectrum
            if args.isopol==1
                i1=find(kp2 > 1/L^2 & kp2 <= 1/rhoe^2 & kabs2 <= 10^(2*kmax));
                Pkern=zeros(args.n,2*args.n);
                Pkern( i1 ) = kabs2(i1).^(-2)/L^(2/3) ...
                    .*exp(-sqrt(kp2(i1))*rho-abs(kz(i1))*va/wic).*dky(i1).*dkiz(i1);
                Pol=(Pkern+kern)/args.ratioTP./kabs2;
            else
                Pol=kern/args.ratioTP./kabs2;
            end
            
            
            %% Add Spectra to PSD
            P(1,i,k)=sum(sum( ky.^2.*Tor ...
                       +(kx.*kz).^2./kp2.*Pol ));
            P(2,i,k)=sum(sum( kx.^2.*Tor ...
                       +(ky.*kz).^2./kp2.*Pol ));
            P(3,i,k)=sum(sum( kp2.*Pol ));
            P(4,i,k)=sum(P(1:3,i,k));
%             if isopol==1
%                 P(5,i,k)=sum(sum( (1-ratioTP/(1+ratioTP))*kp2.*kern./kabs2 ));
%                 P(6,i,k)=sum(sum( (1-ratioTP/(1+ratioTP))*kp2.*Pkern./kabs2 ));
%             end
        else
            P(i,k)=sum(sum(kern));
        end
    end
end

%% Multiply Spectra with factor to obtain PSD in nT^2/Hz
P=2*args.B^2/args.v/L^(1-cb).*P;