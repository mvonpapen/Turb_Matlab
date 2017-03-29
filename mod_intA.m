%% Model PSD for interval A of Rev A [von Papen et al., 2014]
%%
%% Date: 12.06.2014

clear all

%% Load data
load turb_at_saturn_v5.mat
i=Orbs{3}(39:96); %interval A approx. 10h long

%% Electron Temperatures after Schippers et al., 2008
Te_cold(r<=8.9)=6.8e-5*r(r<=8.9)'.^5.9; Te_cold(r>8.9)=4245*r(r>8.9).^(-2.3);
Te_hot(r<=9.3)=0.2*r(r<=9.3).^4.3; Te_hot(r>9.3)=1.2e6*r(r>9.3).^(-2.7);
Te_cold=Te_cold';
Te_hot=Te_hot';

%% Gyro radii
rhow=sqrt(2*18*amu*TW*eV)/eV./B/1e-9;
rhoe_cold=sqrt(amu/1836*Te_cold*eV)/eV./B/1e-9;
rhoe_hot=sqrt(amu/1836*Te_hot*eV)/eV./B/1e-9;


dat=importpsd('~/IDL/turbana/RevAin.gws'); %load gws of intA

%% Construct 1h averages
for j=1:10;
    k=(j-1)*6+1:min(j*6,58);
    mB(j)=mean(B(i(k)));
    mtheta(j)=meanangle(theta_vB(i(k)));
    mv(j)=norm(mean(vrel(i(k),:),1));
    mrhow(j)=mean(rhow(i(k)));
    mrhoe_cold(j)=mean(rhoe_cold(i(k)));
    mrhoe_hot(j)=mean(rhoe_hot(i(k)));
    mH(j)=mean(HW(i(k)));
    mr(j)=mean(r_cyl(i(k))); 
    mT(j)=mean(TW(i(k)));
    mTe_cold(j)=mean(Te_cold(i(k)));
    mTe_hot(j)=mean(Te_hot(i(k)));
    mPani(j)=mean(geomean(sum(P(i(k),20:120,1:2),3)./P(i(k),20:120,3),2)); 
    mva(j)=mean(va(i(k)));
    mE(j)=geomean(Ekprho(i(k)));
end

%% Calculate synthetic spectra
n=1000;
fi=logspace(-4,0.5,20);
Pmod=NaN(10,4,20);
for j=1:10;
    Pmod_cold(j,:,:)=cbspec(fi,mtheta(j),n,2,mH(j)*Rs,mrhow(j),mv(j),mva(j), ...
        mB(j)*1e-9,log10(1./[(2*mH(j)*Rs) mrhoe_cold(j)/100]),mPani(j),mrhoe_cold(j)); 
    Pmod_hot(j,:,:)=cbspec(fi,mtheta(j),n,2,mH(j)*Rs,mrhow(j),mv(j),mva(j), ...
        mB(j)*1e-9,log10(1./[(2*mH(j)*Rs) mrhoe_hot(j)/100]),mPani(j),mrhoe_hot(j)); 
end
Pperp_cold=squeeze(mean(sum(Pmod_cold(:,1:2,:),2),1));
Ppara_cold=squeeze(mean(Pmod_cold(:,3,:),1));
Pperp_hot=squeeze(mean(sum(Pmod_hot(:,1:2,:),2),1));
Ppara_hot=squeeze(mean(Pmod_hot(:,3,:),1));

%% Calculate noise for synthetic spectra
c=10; %This is a empirically derived factor to correct for high energy of synthetic spectra
si=5/3; alias_sn_tmp=0; dt=0.14;
for ni=1:20000;
    alias_sn_tmp=alias_sn_tmp+(fi*2*dt).^(si)./(2*ni-fi*2*dt).^(si)+(fi*2*dt).^(si)...
        ./(2*ni+fi*2*dt).^(si);
end %podesta,2006 eq(6)
clear ni si
for j=1:10
    for ch=1:3
        nois_cold(j,ch,:)=25e-6./fi+4.9e-3^2/6*dt ...
            +alias_sn_tmp.*(squeeze(Pmod_cold(j,ch,:))'/c+25e-6./fi);
        nois_hot(j,ch,:)=25e-6./fi+4.9e-3^2/6*dt ...
            +alias_sn_tmp.*(squeeze(Pmod_hot(j,ch,:))'/c+25e-6./fi);
    end
end
Pperp_noise_cold=squeeze(mean(sum(Pmod_cold(:,1:2,:)/c+nois_cold(:,1:2,:),2),1));
Ppara_noise_cold=squeeze(mean(Pmod_cold(:,3,:)/c+nois_cold(:,3,:),1));
Pperp_noise_hot=squeeze(mean(sum(Pmod_hot(:,1:2,:)/c+nois_hot(:,1:2,:),2),1));
Ppara_noise_hot=squeeze(mean(Pmod_hot(:,3,:)/c+nois_hot(:,3,:),1));


%% Plot results
figure
subplot(1,2,1) %without noise
loglog(dat(:,1),dat(:,3),dat(:,1),dat(:,2),...
    fi,Pperp_cold,fi,Ppara_cold,fi,Pperp_hot,fi,Ppara_hot)
xlabel('f [Hz]'), ylabel('PSD [nT^2/Hz]')
title('without correction')
subplot(1,2,2) %with noise and energy correction
loglog(dat(:,1),dat(:,3),dat(:,1),dat(:,2),...
    fi,Pperp_noise_cold,fi,Ppara_noise_cold,fi,Pperp_noise_hot,fi,Ppara_noise_hot)
xlabel('f [Hz]'), ylabel('PSD [nT^2/Hz]')
title('with noise and energy correction (c=10)')