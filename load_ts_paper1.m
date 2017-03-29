%% Everything from scratch: load time series and global wavelet spectra;
%% read out parameters; construct models for density, temperature,
%% scaleheight, plasma beta, gyroradius and inertial length; calculate
%% noise level (for FGM=0) and SNR; find spectral slope for fixed kr frame
%% with SNR>x; calculate heating rate

% clear all

%% Thresholds and fit ranges
dBincmax=1; %maximum allowed value for increments
thr=5; %SNR threshold
fit=[2 50]*sqrt(2); %fit range for kprho, sqrt(2) for correction of gyro radius
dBthr=[0.1 3]; % Threshold for rms


%% Define Variables so that code is faster (max(ind)=n)
n=length(files); %14916; %This value is found by testing lines 53-61 + 89-93
r=NaN(n,1);
ltime=NaN(n,1);
offline=NaN(n,1);
fnam=cell(n,1);
kurt=NaN(n,3);
stdB=NaN(n,3);
vecB=NaN(n,3);
B=NaN(n,1);
stdBgws=NaN(n,3);
vrel=NaN(n,3);
utcnum=NaN(n,2);
zdp=NaN(n,1);
theta_vB=NaN(n,1);
P=NaN(n,156,3);
Pkprho=NaN(n,156,3);
F=NaN(n,156,3);
kappa_kprho=NaN(n,2);
kprho=NaN(n,156);
QLstrong_f=NaN(n,156);
QLstrong=NaN(n,1);
tauc=NaN(n,3);
noise=NaN(n,156);
SNR=NaN(n,156);


%% Basic Parameters
Rs=60268000; % Saturn radius
OmegaS=2*pi/(10.7*3600); %10.7h mean rotation
eV=1.6e-19; % 1 eV = 1.6e-19 J
amu=1.66e-27; % atomic mass unit
mu0=4e-7*pi;
c=2.998e8; %speed of light
eps0=8.8542e-12;

% %% Load all time series
% fold=dir('/home/vonpapen/PhD_archive/Cassini/*');
ind=1;
% for i=3:length(fold);
%     files=dir(['/home/vonpapen/PhD_archive/Cassini/' fold(i).name '/10m*.dat']);
    for j=1:length(file)
%         [dat,str]=importts(['/home/vonpapen/PhD_archive/Cassini/' fold(i).name '/' files(j).name]);
        [dat,str]=importts(file{j});
        tmpoff=str.offline;
        tmpr=str.dist;
        dBinc=max(max(dat(2:end,2:4)-dat(1:end-1,2:4)));
        if tmpoff<=0.14 && tmpr>=6 && tmpr<=20 && dBinc<dBincmax
%             fnam{ind}=['/home/vonpapen/PhD_archive/Cassini/' fold(i).name '/' files(j).name];
            fnam{ind}=file{j};
            offline(ind)=str.offline;
            r(ind)=str.dist;
            ltime(ind)=str.ltime;
            kurt(ind,:)=kurtosis(dat(:,2:4));
%             kurt10s(ind,:)=kurtosis(dat(72:end,:)-dat(1:end-71,:));
            for k=1:3
                % Linearer fit von Zeitreihe abziehen vor Berechnung des rms
                tmp=polyfit(dat(:,1),dat(:,k+1),1);
                stdB(ind,k)=std(dat(:,k+1)-tmp(2)-dat(:,1)*tmp(1));
                % Calculate autocorrelation time
                ac=autocorr(dat(:,k+1),2000);
                tauc(ind,k)=find(ac<=exp(-1),1,'first')*0.14;
            end
            vecB(ind,:)=str.vecB;
            % Load PSD
            [dat,str]=importpsd(strrep(fnam{ind},'.dat','.wti'));
            rev{ind}=[cell2mat(str.rev{1}) cell2mat(str.rev{2})];
            B(ind)=str.meanB;
            stdBgws(ind,:)=str.std;
            vrel(ind,:)=str.vcas*1e3; %in m/s !
            utc(ind,:)=[str.utc{1} str.utc{2}];
            utcnum(ind,1)=datenum(utc(ind,1),'yyyy-mm-ddTHH:MM:SS');
            utcnum(ind,2)=datenum(utc(ind,2),'yyyy-mm-ddTHH:MM:SS');
            zdp(ind)=str.z(1);
            P(ind,:,:)=dat(:,2:4);
            F(ind,:,:)=dat(:,5:7);
            ind=ind+1;
            f=dat(:,1);
        end
    end
% end
r_cyl=sqrt(r.^2-zdp.^2);

clear dat str ind i j k files fold tmpoff tmpr

%% dBrel, nw, H, T, rho_w, lambda_w, V_A & sls
dBrel=sqrt(sum(stdB.^2,2))./B;
scaleheight=sqrt(r_cyl.^2/3/4.56); %scale height from thomsen, 2010, eq.2'
nw=1.38E6*r.^(-5.68)*1e6.*exp(-abs(zdp)./scaleheight); %total ion number density in m^-3
rho=nw*amu*18;
T=3/2*18*amu*(0.6*OmegaS)^2.*(scaleheight.*Rs).^2/eV; %in eV
rhow=sqrt(2*18*amu*T*eV)/eV./B/1e-9;
lambdaw=c./sqrt(eV^2*nw/18/amu/eps0);
beta = nw*eV.*T ./ ( (B*1e-9).^2/2/mu0 );
va=B*1e-9./sqrt(mu0*nw*18*amu);
ficw=eV*B*1e-9/18/amu/2/pi; %eigentlich omega_ic
SLS2=sls(mean(utcnum,2),ltime,2);
SLS3=sls(mean(utcnum,2),ltime,3);
SLS4S_ui=sls(mean(utcnum,2),ltime,4)';
SLS4S_lesia=sls(mean(utcnum,2),ltime,6)';


%% Noise
si=2.4; alias_sn=0;
for n=1:10000
    alias_sn=alias_sn+(f*2*0.14).^(si)./(2*n-f*2*0.14).^(si)+(f*2*0.14).^(si)...
        ./(2*n+f*2*0.14).^(si); %podesta,2006 eq(6)
end
for i=1:length(r)
    if B(i)<40
        noise(i,:)=75e-6./f+3*4.9e-3^2/6*0.14+alias_sn.*squeeze(sum(P(i,:,:),3))'; %only for FGM=0
    else
        noise(i,:)=75e-6./f+3*48.8e-3^2/6*0.14+alias_sn.*squeeze(sum(P(i,:,:),3))'; %only for FGM=1
    end
    SNR(i,:)=sum(P(i,:,:),3)./noise(i,:);
end

%% Theta, kprho, kappa
for i=1:1200 %length(r)
    theta_vB(i)=vecang(vecB(i,:),vrel(i,:));
    if theta_vB(i)>90, theta_vB(i)=180-theta_vB(i); end
    % Gyro radius
    kprho(i,:)=2*pi.*f./norm(vrel(i,:))/sind(theta_vB(i))*rhow(i);
    fi=[find(kprho(i,:)>fit(1) & f'>2e-2,1,'first') : find(SNR(i,:)<thr | kprho(i,:)>fit(2),1,'first')-1];
    kappa_kprho(i,:)=fitkappa(f,sum(P(i,:,:),3),f(fi));
    Pkprho(i,:,:)=P(i,:,:)*norm(vrel(i,:))*sind(theta_vB(i))/2/pi/rhow(i);
    % Inertial length
    kpd(i,:)=2*pi.*f./norm(vrel(i,:))/sind(theta_vB(i))*lambdaw(i);
    fi=[find(kpd(i,:)>fit(1) & f'>2e-2,1,'first') : find(SNR(i,:)<thr | kpd(i,:)>fit(2),1,'first')-1];
    kappa_kpd(i,:)=fitkappa(f,sum(P(i,:,:),3),f(fi));
    % Ion cyclotron freq
    fi=[find(f'/ficw(i)>5 & f'>2e-2,1,'first') : find(SNR(i,:)<thr | f'/ficw(i)>100,1,'first')-1];
    kappa_ficw(i,:)=fitkappa(f,sum(P(i,:,:),3),f(fi));
end
Ekprho=turbint(kprho,sum(Pkprho,3),fit,SNR,thr);

clear i tmp k fi n


%% Heating rates
for i=1:length(rhow)
    QLstrong_f(i,:)=(sum(P(i,:,:),3)*1e-18).^(3/2).*f'.^(7/2)/sqrt(mu0^3*rho(i)) ...
        *(2*pi/norm(vrel(i,:))/sind(theta_vB(i)))^2*rhow(i);
    fi=find(kprho(i,:)'>fit(1) & f>2e-2,1,'first');
    QLstrong(i)=geomean(QLstrong_f(i,...
        fi:fi+find(SNR(i,fi:end)'<thr | kprho(i,fi:end)'>fit(2),1,'first')-2),2);
end
QSstrong=(sqrt(sum(stdB.^2,2))*1e-9).^3.*rhow./sqrt(mu0^3*rho)...
    ./(mean(tauc,2).*sqrt(sum(vrel.^2,2))).^2;

%% Heating between 6-20 Rs
for i=1:14;
%     j=find(min(stdB,[],2)>=dBthr(1) & max(stdB,[],2)<=dBthr(2) & r_cyl>6+(i-1) & r_cyl<6+i ...
%         & ~isnan(QLstrong) & abs(zdp)./scaleheight<=1 & B>5);
    j=find(min(stdB,[],2)>=dBthr(1) & max(stdB,[],2)<=dBthr(2) & ~isnan(QLstrong));
    QLstrong_int(i,1)=log10(geomean(QLstrong(j))*pi*((6+i)^2-(6+(i-1))^2)...
        *2*geomean(scaleheight(j))*Rs^3);
    QLstrong_int(i,2)=sqrt( std(log10(QLstrong(j)))^2+std(log10(scaleheight(j)))^2 );
%     j=find(min(stdB,[],2)>=dBthr(1) & max(stdB,[],2)<=dBthr(2) & r_cyl>6+(i-1) & r_cyl<6+i ...
%         & ~isnan(QSstrong) & abs(zdp)./scaleheight<=1 & B>5);
    j=find(min(stdB,[],2)>=dBthr(1) & max(stdB,[],2)<=dBthr(2) & ~isnan(QSstrong));
    QSstrong_int(i,1)=log10(geomean(QSstrong(j))*pi*((6+i)^2-(6+(i-1))^2)...
        *2*geomean(scaleheight(j))*Rs^3);
    QSstrong_int(i,2)=sqrt( std(log10(QSstrong(j)))^2+std(log10(scaleheight(j)))^2 );
end
QLstrong_heat_6_20=sum(10.^QLstrong_int(~isnan(QLstrong_int(:,1)),1));
QSstrong_heat_6_20=sum(10.^QSstrong_int(~isnan(QSstrong_int(:,1)),1));
clear fi i j

%% Number different orbits
o=find(utcnum(2:end,1)-utcnum(1:end-1,1)>0.2);
o=[o(1) o(2) o(2)+250 o(3:end)']; %In RevA closest approach less than 6Rs, thus truncate manually between in and outbound leg
Orbs{1}=1:o(1);
for i=2:length(o); Orbs{i}=[o(i-1)+1:o(i)]; end
clear ans o