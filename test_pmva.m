%% Saturn
Rs=60268000;
v=548E3;
L=1e9;
B=10.2e-9;
kmin=-11; %Umfang bei 17Rs: 6.4E9m->k0=-10
rho=76e3; %<rhoW>=3e5m aus meinen daten
kmax=0; %rho_e=rho_W/1840/18 ~ 3
f=logspace(-1,1,15);
fit=[1 3];
va=58e3;
bounds=[-11 0];
n=1000;
% %% Horbury
% L=1e9;
% rho=1e4;
% bounds=[-11 0];
% n=1000;
% load horbury_data

psi=[2];

%% v
P=cbspectra(f,psi,L,rho,v,va,n,bounds,'exp',[-10/3 -11/3],[2/3 1/3],B);
%% v+va
for i=1:length(psi)
    psi_pos(i)=vecang([0 1]*B,[sind(psi(i))*v cosd(psi(i))*v+va]);
end
P_pos=cbspectra(f,psi_pos,L,rho,v,va,n,bounds,'exp',[-10/3 -11/3],[2/3 1/3],B);
%% v-va
for i=1:length(psi)
    psi_neg(i)=vecang([0 1]*B,[sind(psi(i))*v cosd(psi(i))*v-va]);
end
P_neg=cbspectra(f,psi_neg,L,rho,v,va,n,bounds,'exp',[-10/3 -11/3],[2/3 1/3],B);

%% Spectral Index
for i=1:length(psi);
    kappa_va(i,:)=fitkappa(f,squeeze(P_pos(4,:,i))+squeeze(P_neg(4,:,i)),fit,2,0);
    kappa(i,:)=fitkappa(f,squeeze(P(4,:,i)),fit,2,0); 
end
plot(psi,kappa(:,1),psi,kappa_va(:,1))