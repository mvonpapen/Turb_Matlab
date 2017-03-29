%% Program calculates the k-space distribution according to Cho et al, 2002
%% Test if calculation really gives k_para~k_perp^2/3 -> k_para~k_perp^1/3

% clear all

%% Isolinien
plotiso=1;
f=logspace(-4,1,5);
psi=1;
v=600e3;

%% Verteilung
type='exp';
mirror=0;
normed=0;

%% Gridpunkte
n=500;
kernmin=-20;
kernmax=18;

%% Turbulence
si=-10/3; %-10/3->k^{-5/3}
si2=-11/3; %-11/3->k^{-7/3}
cb=2/3; % 2/3->alfven
cb2=1/3; %1/3->KAW
B=5; %in nT
va=60e3;


%% Parameter
L=1e9;
kmin=-11;
kmax=-1;
rho=1e5;
rhoe=rho/42.85;
kz_max=kmax;%log10(1.6e-19*B/1.67e-27/va);
wic=1.6e-19*B*1e-9/1.67e-27; %wic for protons

%% K-space gridpoints
kperp=logspace(kmin,kmax,n);
%nz log verteilt auf pos UND neg Achse
kz=logspace(kmin,kmax,n);


%% Berechne k-Verteilung
kern=zeros(n);

switch type
    case 'exp'
        %% Single components Alfven cascade
        i1=find(kperp <= 1/rho & sqrt(kperp.^2+kz.^2)>=1/L);
        for j=1:n
            kern(i1,j) = kperp(i1).^si...
                 .*exp(-L^(1-cb).*kz(j)./kperp(i1).^cb);
        end
        %% Single components KAW cascade
        % Proportionality factor, so that E(AW)=E(KAW) at k=krho: c=(L*krho)^(1/3)*exp( abs(kz)*L^(1/3)/krho^(1/3)*(L^(1/3)-1/krho^(1/3)) );
        i1=find(kperp > 1/rho & kperp <=1/rhoe);
        for j=1:n
            kern(i1,j) = rho^(si2-si).*kperp(i1).^si2...
                 .*exp( -L^(1-cb)*rho^(cb-cb2)*kz(j)./kperp(i1).^cb2 );
        end
        %% Single electron scales cascade
        % Proportionality factor, so that E(AW)=E(KAW) at k=krho: c=(L*krho)^(1/3)*exp( abs(kz)*L^(1/3)/krho^(1/3)*(L^(1/3)-1/krho^(1/3)) );
        i1=find(kperp > 1/rhoe);
        for j=1:n
            kern(i1,j) = rho^(si2-si).*kperp(i1).^si2...
                 .*exp( -L^(1-cb)*rho^(cb-cb2)*kz(j)*rhoe^cb2 );
        end        

    case 'gauss'
        %% Single components Alfven cascade (Gaussian centered at
        %% zi(i1).^2./(ky(j)^2+xi(i1).^2).^cb=1
        i1=find(kperp <= 1/rho);
        for j=1:n
            kern(i1,j) = kperp(i1).^si...
                 .*exp(-(L^(1-cb)*kz(j)./kperp(i1).^cb-1).^2);
        end
        %% Single components KAW cascade
        % Proportionality factors, so that AW1*exp(AW2)=alpha*KAW1*exp(beta*KAW2) at k=1/rho:
        % alpha=(L/rho)^(1/3), beta=alpha^(-1)  (already incorporated in bottom eq.)
        i1=find(kperp > 1/rho);
        for j=1:n
            kern(i1,j) = rho^(si2-si).*kperp(i1).^si2...
                 .*exp( -(L^(1-cb)*rho^(cb-cb2)*kz(j)./kperp(i1).^cb2-1).^2 );
        end
        
    case 'heavi'
        %% Single components Alfven cascade (Heavyside function with
        %% zi(i1).^2. <= (ky(j)^2+xi(i1).^2).^cb
        i1=find(kperp <= 1/rho);
        for j=1:n
            i2=find(L^(1-cb)*kz(j)./kperp(i1).^cb<=1);
            kern(i1(i2),j) = B^2/L^(1-cb).*kperp(i1(i2)).^si;
        end
        %% Single components KAW cascade
        i1=find(kperp > 1/rho);
        for j=1:n
            i2=find(L^(1-cb)*rho^(cb-cb2)*kz(j)./kperp(i1).^cb2<=1);
            kern(i1(i2),j) = B^2/L^(1-cb)*rho^(si2-si).*kperp(i1(i2)).^si2;
        end
    case 'delta'
        %% Single components Alfven cascade (Heavyside function with
        %% zi(i1).^2. <= (ky(j)^2+xi(i1).^2).^cb
        i1=find(kperp <= 1/rho);
        for j=1:n
            [tmp,i2]=min((L^(1-cb)*kz(j)-kperp(i1).^cb).^2);
            kern(i1(i2),j) = B^2.*kperp(i1(i2)).^si;
        end
        %% Single components KAW cascade
        i1=find(kperp > 1/rho);
        for j=1:n
            [tmp,i2]=min((L^(1-cb)*rho^(cb-cb2)*kz(j)-kperp(i1).^cb2).^2);
            kern(i1(i2),j) = B^2*rho^(si2-si).*kperp(i1(i2)).^si2;
        end        
    case 'expdamp'
        %% Single components Alfven cascade
        i1=find(kperp <= 1/rho & sqrt(kperp.^2+kz.^2)>=1/L);
        for j=1:n
            kern(i1,j) = kperp(i1).^si...
                 .*exp(-L^(1-cb).*kz(j)./kperp(i1).^cb...
                 -kperp(i1)*rhoe-abs(kz(j))*va/wic);
        end
        %% Single components KAW cascade
        % Proportionality factor, so that E(AW)=E(KAW) at k=krho: c=(L*krho)^(1/3)*exp( abs(kz)*L^(1/3)/krho^(1/3)*(L^(1/3)-1/krho^(1/3)) );
        i1=find(kperp > 1/rho & kperp <=1/rhoe);
        for j=1:n
            kern(i1,j) = rho^(si2-si).*kperp(i1).^si2...
                 .*exp( -L^(1-cb)*rho^(cb-cb2)*kz(j)./kperp(i1).^cb2...
                 -kperp(i1)*rhoe-abs(kz(j))*va/wic);
        end
        %% Single electron scales cascade
        % Proportionality factor, so that E(AW)=E(KAW) at k=krho: c=(L*krho)^(1/3)*exp( abs(kz)*L^(1/3)/krho^(1/3)*(L^(1/3)-1/krho^(1/3)) );
        i1=find(kperp > 1/rhoe);
        for j=1:n
            kern(i1,j) = rho^(si2-si).*kperp(i1).^si2...
                 .*exp( -L^(1-cb)*rho^(cb-cb2)*kz(j)*rhoe^cb2...
                 -kperp(i1)*rhoe-abs(kz(j))*va/wic);
        end 
end

if mirror==1
    A=1e60*B^2/L^(1-cb); sigma=0.2; sigmaz=1e-7/rho;
    for j=1:n
        i1=find(kz(j)./kperp<0.1);
        kern(i1,j)=kern(i1,j) + A/(2*pi*sqrt(sigma*sigmaz))...
            *exp(-log10(10*kperp(i1)'*rho).^2/2/sigma^2 ...
                 -kz(j).^2/2/sigmaz);
    end
end

kern=B^2/L^(1/3)*kern;
i=find(kern<10.^kernmin);
kern(i)=10.^kernmin;

%% Normalize each column (kperp) to max value
if normed==1
    norm=max(kern,[],2);
    kern=kern./norm(:,ones(500,1));
end

%% Plot results
figure;
subplot(2,1,1)
contourf(log10(kperp*rho),log10(kz*rho),log10(kern)',linspace(kernmin,log10(max(max(kern))),20)), shading flat
hold all
% plot(log10(kperp*rho),2/3*log10(kperp*rho*L^(-1/2)))
% plot([log10(kernmin*rho) log10(kernmax*rho)],[1 1]*kz_max,'--white','linew',3)
% plot([-8 -0.1],[-8 -0.1].*(2/3)+0.7,'--white','linew',3)
% plot([0 1.9],[0 1.9].*(1/3)+0.7,'--white','linew',3)
% plot([2 4],[0 0]+1.3,'--white','linew',3)
subplot(2,1,2)
contourf(log10(kperp*rho),log10(kz*rho),log10(kern)',linspace(kernmin,log10(max(max(kern))),20)), shading flat
hold all, set(gca,'YDir','Reverse')
% plot(log10(kperp),2/3*log10(kperp*L^(-1/2)))
% plot(log10(1/rho*[1e-6 1e-1]),log10((1/rho*[1e-6 1e-1]).^(2/3)/100),'--white','linew',3)
% plot(log10(1/rho*[1 1e8]),log10((1/rho*[1 1e8]).^(1/3)*1e-4),'--white','linew',3)

if plotiso==1
    %% Plot Isocontours in k_perp,k_para plot
    nf=length(f);
    kx=2*pi*f/v;
    k=logspace(-20,10,1000); k=[-fliplr(k) k];
    for i=1:2
        for j=1:nf;
            pl(j,:,:)=[sind(psi)*kx(j)+cosd(psi)*k; cosd(psi)*kx(j)-sind(psi)*k]; 
        end
        i1=find(pl<0); pl(i1)=NaN;
        for j=1:nf;
            j1=~isnan(squeeze(pl(j,1,:)));
            j2=~isnan(squeeze(pl(j,2,:)));
            n=[1:length(j1)];
            j3=intersect(n(j1), n(j2));
            subplot(2,1,i), hold all
            if i==1
                %% Normal log
                plot([-10 log10(squeeze(pl(j,1,j3))*rho)' max(log10(squeeze(pl(j,1,j3))*rho))],...
                    [max(log10(squeeze(pl(j,2,j3))*rho)) log10(squeeze(pl(j,2,j3))*rho)' -10],'--white');
            else
                %% Log für psi>90
                plot([log10(squeeze(pl(j,1,j3))*rho)' min(log10(squeeze(pl(j,1,j3))*rho))],...
                    [log10(squeeze(pl(j,2,j3))*rho)' -10],'--white');
            end
        end
        psi=180-psi;
    end
end