%% Check numerically calculated spectral slope for exponential damping at rhoe
%% inside 11Rs for measured and modelled plasma parameters

clear Pcb_f fi snr_mod nois kappa comp fn


load satturb8
OmegaS=2*pi/(10.8*3600); %10.8h mean rotation
j=0;
si1=-10/3; si2=-11/3;
n=100;
%     [tmp,k]=find(r<=20);
k=1:1136; %Interval.A;
fi=logspace(-2.3,0.5,20) %logspace(-4,0.5,20);
indi=find(fi>0.1 & fi<0.6);

%% Temperatures [Thomsen et al., 2010; Schippers et al., 2008]
HW=sqrt(r.^2/3/8.7); TW=3/2*18*amu*(0.6*OmegaS)^2.*(HW.*Rs).^2/eV; %water group
Te_cold(r<=8.9)=6.8e-5*r(r<=8.9).^5.9; Te_cold(r>8.9)=4245*r(r>8.9).^(-2.3);
Te_hot(r<=9.3)=0.2*r(r<=9.3).^4.3; Te_hot(r>9.3)=1.2e6*r(r>9.3).^(-2.7);
rhow=sqrt(2*18*amu*TW*eV)/eV./B/1e-9;
rhoe=rhow.*sqrt(Te./T/1836/18);

%% Calculate Aliasing
si=7/3; alias_sn_tmp=0; dt=0.14;
for ni=1:10000;
    alias_sn_tmp=alias_sn_tmp+(fi*2*dt).^(si)./(2*ni-fi*2*dt).^(si)+(fi*2*dt).^(si)...
        ./(2*ni+fi*2*dt).^(si);
end %podesta,2006 eq(6)
clear ni si

fprintf('%i loops to go!\n', length(k)); %1136
fprintf('Currently processing No.%04d',j)

for i=1:1136; %k

    j=j+1;
    fprintf('\b\b\b\b%04d',j)

    %% PSD in fixed frequency frame
    ind=find(f'>0.1 & snr(i,:)>5);
    comp(j)=mean(sum(PSD(i,ind,1:2),3)./PSD(i,ind,3));
    Pcb_f(j,:,:)=cbspec(fi,theta_vb(i),2*HW(i)*Rs,rhow(i),norm(vrel(i,1:2)),...
        va(i),n,2,norm(stdB(i,:)),[log10(1/HW(i)/Rs) log10(100/rhoe(i))],comp(j),rhoe(i),[si1 si2]);
    %% Normalize to measured turb intensity
    for ch=1:3
        Pcb_f(j,ch,:)=Pcb_f(j,ch,:)/geomean(Pcb_f(j,4,indi),3)*Pavgf(i)/500; %factor 1/500 empirically determinated

        %% Noise in fixed frequency frame
%         SNR=interp1(f,snr(i,:),fi);
        if FGM(i)==0; 
            nois(j,ch,:)=25e-6./fi+4.9e-3^2/6*dt+alias_sn_tmp.*(squeeze(Pcb_f(j,ch,:))'+25e-6./fi);
        else
            nois(j,ch,:)=25e-6./fi+48.8e-3^2/6*dt+alias_sn_tmp.*(squeeze(Pcb_f(j,ch,:))'+25e-6./fi);
        end;
        Pcb_f(j,ch,:)=Pcb_f(j,ch,:)+nois(j,ch,:);
    end
    nois(j,4,:)=sum(nois(j,1:3,:),2);
    Pcb_f(j,4,:)=sum(Pcb_f(j,1:3,:),2);
    snr_mod(j,:)=Pcb_f(j,4,:)./nois(j,4,:); %SNR relativ zu theor. noise level, "echtes" noise=noise*[1,10]

    %% Calculate Spectral Slope
    ind=find(snr_mod(j,:)<5,1,'first');
    fn(j)=min([fi(end) fi(ind)]);
    kappa(j,:)=fitkappa(fi,Pcb_f(j,4,:),[0.1 fn(j)],4); %snr*n damit das gleiche noise (random) verwendet wird
    kprho_mod(j,:)=2*pi*rhow(i)/norm(vrel(i,1:2))/sind(theta_vb(i))*fi;
    ki=[find(kprho_mod(j,:)>2 & fi>2e-2,1,'first') : ...
        find(snr_mod(j,:)<5 | kprho_mod(j,:)>50,1,'first')-1];
    kappa_sn5(j,:)=fitkappa(fi,sum(Pcb_f(j,:,:),2),fi(ki),5,0);
    clear ki x
%         kapinert(j,:)=fitkappa(fi,Pcb_f(j,4,:),[5e-4 5e-3],4); %inertial range kappa

end
fprintf('\n')


%% OLD STUFF:
%
%% PSD in fixed krho frame
%     ind=find(krho(i,:)>=2 & krho(i,:)<=50 & snr(i,:)>=10);
%     tmp=log10(f(minmax(ind)));
%     fcb(j,:)=logspace(tmp(1),tmp(2),5);
%     Pcb(j,:)=cbspec(fcb(j,:),theta_vb(i),HW(i)*Rs,...
%         rhow(i),norm(vrel(i,:)),va(i),n,2,B(i),[-10 0],0,rhow(i)*sqrt(Te(i)/T(i)/1836));
%% Noise in fixed krho frame
%     SNR=interp1(f,snr(i,:),fcb(j,:));
%     nois(j,:)=Pcb(j,:)./(SNR-1);
%     kappa(j,:)=fitkappa(fcb(j,:),Pcb(j,:)+nois(j,:),minmax(fcb(j,:)));
