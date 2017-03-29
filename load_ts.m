%% Everything from scratch: load time series and global wavelet spectra;
%% read out parameters; construct models for density, temperature,
%% scaleheight, plasma beta, gyroradius and inertial length; calculate
%% noise level (for FGM=0) and SNR; find spectral slope for fixed kr frame
%% with SNR>x; calculate heating rate

clear all

%% Thresholds and fit ranges
dBincmax=1; %maximum allowed value for increments
thr=5; %SNR threshold
fit=[2 50]; %fit range for kprho, sqrt(2) for correction of gyro radius
dBthr=0.1; % Threshold for rms
minB=2;

%% Define Variables so that code is faster (max(ind)=n)
n=33132; %until end of 2012; 21213 until end of 2010; %This value is found by testing lines 53-62 + 90-94
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
z=NaN(n,1);
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
Te=NaN(n,1);
rev=cell(n,1);
utc=cell(n,2);

%% Basic Parameters
Rs=60268000; % Saturn radius
OmegaS=2*pi/(10.8*3600); %10.8h mean rotation
eV=1.6e-19; % 1 eV = 1.6e-19 J
amu=1.66e-27; % atomic mass unit
mu0=4e-7*pi;
c=2.998e8; %speed of light
eps0=8.8542e-12;

%% Load all time series
fold=dir('/home/vpapenm/PhD_archive/Cassini/*');
ind=1;
for i=3:length(fold);
    files=dir(['/home/vpapenm/PhD_archive/Cassini/' fold(i).name '/10m*.dat']);
    for j=1:length(files)
        [dat,str]=importts(['/home/vpapenm/PhD_archive/Cassini/' fold(i).name '/' files(j).name]);
        tmpoff=str.offline;
        tmpr=str.dist;
        dBinc=max(max(dat(2:end,2:4)-dat(1:end-1,2:4)));
        if tmpoff<=0.14 && tmpr>=6 && tmpr<=20 && dBinc<dBincmax
            fnam{ind}=['/home/vpapenm/PhD_archive/Cassini/' fold(i).name '/' files(j).name];
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
                ac=xcorr(dat(:,k+1)-mean(dat(:,k+1)),2000,'unbiased');
                tauc(ind,k)=find(ac(2e3:end)/ac(2001)<=exp(-1),1,'first')*0.14;
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
            z(ind)=str.z(1);
            P(ind,:,:)=dat(:,2:4);
            F(ind,:,:)=dat(:,5:7);
            ind=ind+1;
            f=dat(:,1);
        end
    end
end

%% Longitude systems
SLS2        = sls(mean(utcnum,2),ltime,2);
SLS3        = sls(mean(utcnum,2),ltime,3);
SLS4S_ui    = sls(mean(utcnum,2),ltime,4)';
SLS4S_lesia = sls(mean(utcnum,2),ltime,6)';
MAGS        = sls(mean(utcnum,2),ltime,8);
MAGN        = sls(mean(utcnum,2),ltime,9);
load MP_out

%% Current sheet flapping model of Arridge et al., 2011, JGR, 116, A11205
r_cyl       = sqrt(r.^2-z.^2);
theta_sun   = satangle(mean(utcnum,2)); %positive at southern summer
theta_tilt  = 12; %deg
r0          = 12;
Rh          = 22; %hinging distance, mean of Table 2
G           = 3; %3deg/Rs from Provan et al., 2014, JGR, 119 (see also Table 2)
% Southern magnetic system
psi_PS      = mod(MAGS-G*(r_cyl-r0),360); %z_cs in-phase with Br: => Br>0 <-> Bt<0 leads to z_cs>0
i           = r_cyl<=12;
z_cs_mag(i) = (r_cyl(i)-Rh*tanh(r_cyl(i)/Rh)).*tand(theta_sun(i));
i           = r_cyl>12;
z_cs_mag(i) = (r_cyl(i)-Rh*tanh(r_cyl(i)/Rh)).*tand(theta_sun(i)) ...
                + (r_cyl(i)-r0).*tand(theta_tilt).*cosd(psi_PS(i));
zdp         = z - z_cs_mag';
% Northern magnetic system
psi_PS_N      = mod(MAGN-G*(r_cyl-r0),360); %z_cs in-phase with Br: => Br>0 <-> Bt<0 leads to z_cs>0
i             = r_cyl<=12;
z_cs_mag_N(i) = (r_cyl(i)-Rh*tanh(r_cyl(i)/Rh)).*tand(theta_sun(i));
i           = r_cyl>12;
z_cs_mag_N(i) = (r_cyl(i)-Rh*tanh(r_cyl(i)/Rh)).*tand(theta_sun(i)) ...
                + (r_cyl(i)-r0).*tand(theta_tilt).*cosd(psi_PS(i));
zdp_N         = z - z_cs_mag';
clear i

%% Correct jumps in local time caused by mean([359° 1°])=180°
i=find(ltime(2:end-1)-ltime(1:end-2)>-13 & ltime(2:end-1)-ltime(1:end-2)<-11 ...
        & ltime(3:end)-ltime(2:end-1)>-13 & ltime(3:end)-ltime(2:end-1)<-11)+1;
ltime(i)=mod(ltime(i)+12,24);
            

clear dat str ind i j k files fold tmpoff tmpr

%% all model parameters: dBrel, nw, H, T, rho_w, lambda_w, V_A, sls, etc.
dBrel=sqrt(sum(stdB.^2,2))./B;
rms=sqrt(sum(stdB.^2,2));
%% Scale heights in Rs & Temperatures in eV
HH=sqrt(r.^2/3/1.6); TH=3/2*amu*(0.6*OmegaS)^2.*(HH.*Rs).^2/eV; %hydrogen
HW=sqrt(r.^2/3/8.7); TW=3/2*18*amu*(0.6*OmegaS)^2.*(HW.*Rs).^2/eV; %water group
HH2=sqrt(r.^2/3/2.75); TH2=3/2*2*amu*(0.6*OmegaS)^2.*(HH2.*Rs).^2/eV; %molecular hydrogen
HT=sqrt(r.^2/3/4.56); TT=3/2*16*amu*(0.6*OmegaS)^2.*(HT.*Rs).^2/eV; %total
%% Densities from Thomsen 2010 in m^-3
nH=1.01e10*r.^(-4.28).*exp(-(zdp./HH).^2); %hydrogen
nH2=79.7e6*r.^(-2.88).*exp(-(zdp./HH2).^2); %molecular hydrogen
nW=8.72e12*r.^(-6.62).*exp(-(zdp./HW).^2); %water group
nT=1.38e12*r.^(-5.68).*exp(-(zdp./HT).^2); %total

% %% Densities from Goldstein et al, 2016, JGR (multiplied by 5)
% load density_goldstein2016.mat
% nW = NaN(size(r));
% j  = find(~isnan(SLS3))';
% for i=j
%     k       = min( [ceil(r(i)-4) 16] ); %index of coefficient-matrix
%     nW(i)   = 1e6 * sum(  Cm(:,k) .* cosd( (0:5)' * mod(SLS3(i)+100,360) ) ...
%                       + Sm(:,k) .* sind( (0:5)' * mod(SLS3(i)+100,360) )  ) ...
%                   *exp( -(abs(zdp(i))./HT(i)).^2);
% end

rho=nW*amu*18;

%% Electron Temperatures and densities after Schippers et al., 2008
Te_cold(r<=8.9)=6.8e-5*r(r<=8.9).^5.9; Te_cold(r>8.9)=4245*r(r>8.9).^(-2.3);
Te_hot(r<=9.3)=0.2*r(r<=9.3).^4.3; Te_hot(r>9.3)=1.2e6*r(r>9.3).^(-2.7);
ne_cold=3.8e10*r.^(-3.9);
ne_hot(r<=8.5)=4.8*r(r<=8.5).^(5.1); ne_hot(r>8.5)=13.9e6*r(r>8.5).^(-1.9);
% Te(r<=9)=6.8e-5*r(r<=9).^5.9;
% Te(r>9 & r<13)=4245*r(r>9 & r<13).^-2.3;
% Te(r>=13)=TT(r>=13)/10;


rhow=sqrt(2*18*amu*TW*eV)/eV./B/1e-9;
lambdaw=c./sqrt(eV^2*nW/18/amu/eps0);
beta = nW*eV.*TW ./ ( (B*1e-9).^2/2/mu0 );
va=B*1e-9./sqrt(mu0*rho);


%% Noise
si=7/3; alias_sn=0;
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
for i=1:length(r)
    theta_vB(i)=vecang(vecB(i,:),vrel(i,:));
    if theta_vB(i)>90, theta_vB(i)=180-theta_vB(i); end
    kprho(i,:)=2*pi*f./norm(vrel(i,:))/sind(theta_vB(i))*rhow(i);
    fi=[find(kprho(i,:)>fit(1) & f'>2e-2,1,'first') : find(SNR(i,:)<thr | kprho(i,:)>fit(2),1,'first')-1];
    kappa_kprho(i,:)=fitkappa(f,sum(P(i,:,:),3),f(fi),5,0);
    Pkprho(i,:,:)=P(i,:,:)*norm(vrel(i,:))*sind(theta_vB(i))/2/pi/rhow(i);
end
Ekprho=turbint(kprho,sum(Pkprho,3),fit,SNR,thr);

clear i tmp k fi n


%% Heating rates
for i=1:length(rhow)
    QLstrong_f(i,:)=(sum(P(i,:,:),3)*1e-18).^(3/2).*f'.^(7/2)/sqrt(mu0^3*rho(i)) ...
        *(2*pi/norm(vrel(i,:))/sind(theta_vB(i)))^2*rhow(i);
%     t_nl_f(i,:) = norm(vrel(i,:))^2*sind(theta_vB(i))^2*sqrt(mu0*rho(i)) ...
%         /4/pi^2/rhow(i)./(sum(P(i,:,:),3)*1e-18).^(1/2)./f'.^(5/2);
    fi=find(kprho(i,:)'>fit(1) & f>2e-2,1,'first');
    QLstrong(i)=geomean(QLstrong_f(i,...
        fi:fi+find(SNR(i,fi:end)'<thr | kprho(i,fi:end)'>fit(2),1,'first')-2),2);
%     t_nl(i) = geomean(t_nl_f(i,...
%         fi:fi+find(SNR(i,fi:end)'<thr | kprho(i,fi:end)'>fit(2),1,'first')-2),2);
end
QSstrong=(sqrt(sum(stdB.^2,2))*1e-9).^3.*rhow./sqrt(mu0^3*rho)...
    ./(mean(tauc,2).*sqrt(sum(vrel.^2,2))).^2;

%% Heating between 6-20 Rs
Hx=HW;
for i=1:14;
    j=find(rms>=dBthr & r_cyl>6+(i-1) & r_cyl<=6+i ...
        & ~isnan(QLstrong) & abs(zdp)./Hx<=1 & B>minB & vecB(:,3)<0 & inside_MP);
    QLstrong_int(i,1)=log10(geomean(QLstrong(j))*pi*((6+i)^2-(6+(i-1))^2)...
        *2*geomean(Hx(j))*Rs^3);
    QLstrong_int(i,2)=sqrt( std(log10(QLstrong(j)))^2+std(log10(Hx(j)))^2 );
    j=find(rms>=dBthr & r_cyl>6+(i-1) & r_cyl<6+i ...
        & ~isnan(QSstrong) & abs(zdp)./Hx<=1 & B>minB & vecB(:,3)<0 & inside_MP);
    QSstrong_int(i,1)=log10(geomean(QSstrong(j))*pi*((6+i)^2-(6+(i-1))^2)...
        *2*geomean(Hx(j))*Rs^3);
    QSstrong_int(i,2)=sqrt( std(log10(QSstrong(j)))^2+std(log10(Hx(j)))^2 );
    j=find(rms>=dBthr & r_cyl>6+(i-1) & r_cyl<6+i ...
        & abs(zdp)./Hx<=1 & B>minB & vecB(:,3)<0 & inside_MP);
    Poynt_int(i,1)=log10(geomean(rms(j).^2*1e-18.*va(j)/mu0)*pi*((6+i)^2 ...
        -(6+(i-1))^2)*Rs^2);
    Poynt_int(i,2)=sqrt( std(log10(rms(j).^2*1e-18.*va(j)/mu0))^2 ...
        +std(log10(Hx(j)))^2 );
end
QLstrong_heat_6_20=sum(10.^QLstrong_int(:,1));
QSstrong_heat_6_20=sum(10.^QSstrong_int(:,1));
Poynt_6_20=sum(10.^Poynt_int(:,1));
clear fi i j Hx

%% Number of different in- and outbound legs
clear Orbs
o=find(utcnum(2:end,1)-utcnum(1:end-1,1)>5); %0.2 ~ 5hours
Orbs{1}=1:o(1);
for i=2:length(o)
    Orbs{i}= [o(i-1)+1:o(i)];
end
Orbs{i+1}  = [o(i)+1:length(utcnum(:,1))];
clear ans i o Leg
k=1;
for i=1:length(Orbs);
    if ismember(i,[23 24])
        Leg{k} = Orbs{i};
        k =k+1;
        continue
    end
    log = sign( diff( r(Orbs{i}) ) ); %-1=in, 1=out
    log = diff([log(1); log; log(end)]);
    if sum(log)==2 && sum(log~=0)==1
        j = find(log==2);
        tmp = abs(diff(ltime(Orbs{i}(j-1:j+1))));
        if tmp(1)>tmp(2)
            j=j-1;
        end
        Leg{k}   = Orbs{i}(1:j);
        Leg{k+1} = Orbs{i}(j+1:end);
        k = k+2;
    elseif sum(log)==2 && sum(log~=0)==2
        j   = find(log~=0);
        Leg{k}   = Orbs{i}(1:j(1));
        Leg{k+1} = Orbs{i}(j(2):end);
        k = k+2;
    elseif sum(log)==0 && sum(log~=0)==0
        Leg{k}   = Orbs{i};
        k = k+1;
    end        
end
clear i in out k tmp

%% Characteristic times
tau_kaw  = vrel(:,2).*mean(tauc,2) ./va./rhow.*HW*Rs;
tau_kaw2 = vrel(:,2).*mean(tauc,2).^(5/3) ./rhow.*(HW*Rs).^(1/3);
tau_nl   = (vrel(:,2).*mean(tauc,2)).^2./dBrel./va./rhow;
tau_nl2  = (vrel(:,2).*mean(tauc,2))./dBrel./va;



save satturb_2004-2012_MAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot sorted by radial distance
i=find(rms>=dBthr & r>6 & r<20 & abs(zdp)./HW<=1 & B>minB & vecB(:,3)<0 & inside_MP);
errorbar(r_cyl(i),kappa_kprho(i,1),kappa_kprho(i,2),'.')



%% Plot sorted by local time
figure
for j=1:4
    subplot(2,2,j)
    if j<4
        i=find(rms>=dBthr & r_cyl>6 & r_cyl<18 & abs(zdp)./HW<1 & inside_MP ...
            & ltime>=(j-1)*6+3 & ltime<mod(j*6+3,24));
    else
        i=find(rms>=dBthr & r_cyl>6 & r_cyl<18 & abs(zdp)./HW<1 & inside_MP ...
            & ltime>=(j-1)*6+3 | ltime<mod(j*6+3,24)); 
    end
    errorbar(r(i),kappa_kprho(i,1),kappa_kprho(i,2),'.');
    ylim([-5 0]);
%     semilogy(r_cyl(i),kurt(i,:),'.');
%     ylim([1 100]);    
    xlim([6 18]);
    title(['LT ', mat2str((j-1)*6+3), 'h-', mat2str(mod(j*6+3,24)), 'h']);
    m=nanmean(kappa_kprho(i(r_cyl(i)>10),1));
    s=nanstd(kappa_kprho(i(r_cyl(i)>10),1));
    legend(['\kappa=', mat2str(m,3), '\pm', mat2str(s,2), 'for r>10Rs']);
    clear m s i
end

%% Plot sorted by MAG-S
figure
for j=1:4
    subplot(2,2,j)
    if j<4
        i=find(rms>=dBthr & r_cyl>6 & r_cyl<18 & abs(zdp)./HW<1 & inside_MP ...
            & MAGS>=(j-1)*90+45 & MAGS<mod(j*90+45,360));
    else
        i=find(rms>=dBthr & r_cyl>6 & r_cyl<18 & abs(zdp)./HW<1 & inside_MP ...
            & MAGS>=(j-1)*90+45 | MAGS<mod(j*90+45,360));
    end
    errorbar(r(i),kappa_kprho(i,1),kappa_kprho(i,2),'.');
    ylim([-5 0]);
%     semilogy(r_cyl(i),kurt(i,:),'.');
%     ylim([1 100]);    
    xlim([6 18]);
    title(['MAG ', mat2str((j-1)*90+45), '^\circ-', mat2str(mod(j*90+45,360)), '^\circ']);
    m=nanmean(kappa_kprho(i,1));
    s=nanstd(kappa_kprho(i,1));
    legend(['\kappa=', mat2str(m,3), '\pm', mat2str(s,2)]);
    clear m s i
end