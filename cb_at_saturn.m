%% Check numerically calculated spectral slope for exponential damping at rhoe
%% inside 11Rs for measured and modelled plasma parameters

% clear Pcb_f fi snr_mod nois kappa comp fn
% clear all

% load Saturn-MAG-Turbulence_data_Diss_vonPapen
% load paper1_index
% k=ind; %Index of paper 1
j=0;
n=200;
fun=2; %expdamp model
fi=logspace(-2.5,0.5,50); %logspace(-4,0.5,20)

%% Basic Parameters
Rs=60268000; % Saturn radius
OmegaS=2*pi/(10.8*3600); %10.8h mean rotation
eV=1.6e-19; % 1 eV = 1.6e-19 J
amu=1.66e-27; % atomic mass unit
mu0=4e-7*pi;

%% Correct height over magnetic dipole
Rh=25; %hinging distance
theta_sun=satangle(mean(utcnum,2));
z_cs=(r-Rh*tanh(r/Rh)).*tand(theta_sun);
z=zdp-z_cs;

%% Scale heights & Temperatures
% HH=sqrt(r.^2/3/1.6); TH=3/2*amu*(0.6*OmegaS)^2.*(HH.*Rs).^2/eV; %hydrogen
HW=sqrt(r.^2/3/8.7); TW=3/2*18*amu*(0.6*OmegaS)^2.*(HW.*Rs).^2/eV; %molecular hydrogen
% HH2=sqrt(r.^2/3/2.75); TH2=3/2*2*amu*(0.6*OmegaS)^2.*(HH2.*Rs).^2/eV; %water group
HT=sqrt(r.^2/3/4.56); TT=3/2*16*amu*(0.6*OmegaS)^2.*(HT.*Rs).^2/eV; %total

%% Electron Temperatures after Schippers et al., 2008
Te_cold(r<=8.9)=6.8e-5*r(r<=8.9).^5.9; Te_cold(r>8.9)=4245*r(r>8.9).^(-2.3);
Te_hot(r<=9.3)=0.2*r(r<=9.3).^4.3; Te_hot(r>9.3)=1.2e6*r(r>9.3).^(-2.7);

%% Densities
% nH=1.01e10*r.^(-4.28).*exp(-(zdp./HH).^2); %hydrogen
% nH2=79.7e10*r.^(-2.88).*exp(-(zdp./HH2).^2); %molecular hydrogen
% nW=8.72e12*r.^(-6.62).*exp(-(zdp./HW).^2); %water group
nT=1.38e12*r.^(-5.68).*exp(-(z./HT).^2); %total
rho=nT*amu*16;
va=B./sqrt(rho*mu0)*1e-9;

%% Gyro radii
rhow=sqrt(18*amu*TW*eV)/eV./B/1e-9;
rhoe_cold=sqrt(amu/1836*Te_cold'*eV)/eV./B/1e-9;
rhoe_hot=sqrt(amu/1836*Te_hot'*eV)/eV./B/1e-9;

%% Calculate Aliasing
si=7/3; alias_sn_tmp=0; dt=0.14;
for ni=1:10000;
    alias_sn_tmp=alias_sn_tmp+(fi*2*dt).^(si)./(2*ni-fi*2*dt).^(si)+(fi*2*dt).^(si)...
        ./(2*ni+fi*2*dt).^(si);
end %podesta,2006 eq(6)
clear ni si

fprintf('%i loops to go!\n', length(r)); %1136
fprintf('Currently processing No.%05d',j)

for i=k; %1:length(r); %1:1136;

    j=j+1;
    fprintf('\b\b\b\b\b%05d',j)

    %% P in fixed frequency frame
    ind=[find(kprho(i,:)>2 & f'>1e-2,1,'first'): ...
        find(SNR(i,:)<5 | kprho(i,:)>50,1,'first')-1];
    TPratio=squeeze(geomean(sum(P(i,ind,1:2),3)./P(i,ind,3)));
    if ~isempty(ind)
        indi=find(fi>=f(ind(1)) & fi<=f(ind(end)));
        Pcb_cold=cbspec(fi(indi),theta_vB(i),n,fun,2*HW(i)*Rs,rhow(i),norm(vrel(i,:)),...
            va(i),B(i),[log10(1/10/HW(i)/Rs) log10(100/rhoe_cold(i))],TPratio,rhoe_cold(i));
        Pcb_cold=sum(Pcb_cold)';
        Pcb_hot=cbspec(fi(indi),theta_vB(i),n,fun,2*HW(i)*Rs,rhow(i),norm(vrel(i,:)),...
            va(i),B(i),[log10(1/10/HW(i)/Rs) log10(100/rhoe_hot(i))],TPratio,rhoe_hot(i));
        Pcb_hot=sum(Pcb_hot)';
        %% Normalize to measured turb intensity
        E=geomean(sum(P(i,ind,:),3));
        cE_hot(j)=geomean(Pcb_hot)/E;
        cE_cold(j)=geomean(Pcb_cold)/E;
        Pcb_cold=Pcb_cold/geomean(Pcb_cold)*E;
        Pcb_hot=Pcb_hot/geomean(Pcb_hot)*E;
        %% Noise in fixed frequency frame
        if B(i)<40; 
            nois_cold=75e-6./fi(indi)'+3*4.9e-3^2/6*dt+alias_sn_tmp(indi)'...
                .*(Pcb_cold+75e-6./fi(indi)');
            nois_hot=75e-6./fi(indi)'+3*4.9e-3^2/6*dt+alias_sn_tmp(indi)'...
                .*(Pcb_hot+75e-6./fi(indi)');
        else
            nois_cold=75e-6./fi(indi)'+3*48.8e-3^2/6*dt+alias_sn_tmp(indi)'...
                .*(Pcb_cold+75e-6./fi(indi)');
            nois_hot=75e-6./fi(indi)'+3*48.8e-3^2/6*dt+alias_sn_tmp(indi)'...
                .*(Pcb_hot+75e-6./fi(indi)');
        end;
        Pcb_cold=Pcb_cold+nois_cold;
        Pcb_hot=Pcb_hot+nois_hot;

        %% Calculate Spectral Slope
        kappa_cold(j,:)=fitkappa(fi(indi),Pcb_cold,minmax(f(ind)));
        kappa_hot(j,:)=fitkappa(fi(indi),Pcb_hot,minmax(f(ind)));
        kappa_cold_nn(j,:)=fitkappa(fi(indi),Pcb_cold-nois_cold,minmax(f(ind)));
        kappa_hot_nn(j,:)=fitkappa(fi(indi),Pcb_hot-nois_hot,minmax(f(ind)));
    else
        kappa_cold(j,:)=[NaN NaN];
        kappa_hot(j,:)=[NaN NaN];
        kappa_cold_nn(j,:)=[NaN NaN];
        kappa_hot_nn(j,:)=[NaN NaN];
    end
end
fprintf('\n')

clear ind E nois_cold nois_hot j i
% clear indi Pcb_cold Pcb_hot