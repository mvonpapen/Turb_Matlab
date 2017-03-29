clear winkrho winkdi winkfic winkf
mindp=5;
ss=4; % 1/step size
is=5; % interval size [x:is*x]
n=3*ss;

%  & f<0.6 & f>2e-2 & squeeze(sum(PSD(i,:,:),3))'>1e-4
%  & sn>10 & f'>2e-2
for j=0:n-1
    for i=1:1136
        sni=find(SNR(i,:)<5,1,'first')-1;
        x=10.^(j/ss-1);
%         kappa.krho(i,j+1,:)=fitkappa(f,sum(PSD(i,:,:),3),f(krho(i,:)>=x & krho(i,:)<=is*x & f'<f(sni)'),mindp);
        kappa.kprho(i,j+1,:)=fitkappa(f,sum(PSD(i,:,:),3),f(kprho(i,:)>=x & kprho(i,:)<=is*x & f'<f(sni)'),mindp);
        kappa.kdi(i,j+1,:)=fitkappa(f,sum(PSD(i,:,:),3),f(kd(i,:)>=x & kd(i,:)<=is*x & f'<f(sni)'),mindp);
        kappa.fic(i,j+1,:)=fitkappa(f,sum(PSD(i,:,:),3),f(f/ficW(i)>=x & f/ficW(i)<=is*x & f<f(sni)),mindp);
        kappa.f(i,j+1,:)=fitkappa(f,sum(PSD(i,:,:),3),f(f>=x/10 & f<=is*x/10 & f<f(sni)),mindp);
%         kappa.k(i,j+1,:)=fitkappa(f,sum(PSD(i,:,:),3),f(krho(i,:)./rhoi(i)>=x*1e-6...
%             & krho(i,:)./rhoi(i)<=is*x*1e-6 & f'<f(sni)'),mindp);
    end;
end


%% windows
win=10.^([0:n-1]/ss-1);
winf=10.^([0:n-1]/ss-1)/10;
% wink=10.^([0:n-1]/ss-1)*1e-6;


%% Plot Results
for i=1:9
    i1=find(r>5.5+i & r<6.5+i);
    % Frequency
    subplot(2,2,1)
    semilogx(winf,nanmean(kappa.f(i1,:,1)))
    hold all
    % kperp gyro radius
    subplot(2,2,2)
    semilogx(win,nanmean(kappa.kprho(i1,:,1)))
    hold all
    % ion cyclotron freq.
    subplot(2,2,3)
    semilogx(win,nanmean(kappa.fic(i1,:,1)))
    hold all
    % inertial length
    subplot(2,2,4)
    semilogx(win,nanmean(kappa.kdi(i1,:,1)))
    hold all
%     % kprho
%     subplot(3,2,5)
%     semilogx(win,nanmean(kappa.kprho(i1,:,1)))
%     hold all
%     % k
%     subplot(3,2,6)
%     semilogx(wink,nanmean(kappa.k(i1,:,1)))
%     hold all    
end

%% Labels usw.
subplot(2,2,1), xlim([1e-2 0.5]), ylim([-3 -1])
title('[ f'' : 5f'' ]')
xlabel('f'' [Hz]'), ylabel('Spectral Index'), grid on
plot([1e-3 1e2],[-2.6 -2.6],'--black')
subplot(2,2,2), xlim([1e-1 50]), ylim([-3 -1])
title('[ k''_\perp \rho_W : 5k''_\perp \rho_W ]')
xlabel('k''_\perp \rho_W'), ylabel('Spectral Index'), grid on
plot([1e-3 1e2],[-2.6 -2.6],'--black')
subplot(2,2,3), xlim([1e-1 70]), ylim([-3 -1])
title('[ f''/f_{c,W} : 5f''/f_{c,W} ]')
xlabel('f''/f_{c,W}'), ylabel('Spectral Index'), grid on
plot([1e-3 1e2],[-2.6 -2.6],'--black')
subplot(2,2,4), xlim([1e-1 70]), ylim([-3 -1])
title('[ k''_\perp \lambda_W : 5k''_\perp \lambda_W ]')
xlabel('k''_\perp \lambda_W'), ylabel('Spectral Index'), grid on
plot([1e-3 1e2],[-2.6 -2.6],'--black')
% subplot(3,2,5), xlim([1e-1 70]), ylim([-3 -1])
% title('[k_\perp'' \rho_W:5k_\perp'' \rho_W]')
% xlabel('k_\perp'' \rho_W'), ylabel('Spectral Index'), grid on
% plot([1e-3 1e2],[-2.6 -2.6],'--black')
% subplot(3,2,6), xlim([1e-7 1e-4]), ylim([-3 -1])
% title('[k'':5k'']')
% xlabel('k'''), ylabel('Spectral Index'), grid on
% plot([1e-7 1],[-2.6 -2.6],'--black')
legend('7 R_s', '8 R_s', '9 R_s', '10 R_s', '11 R_s', '12 R_s', '13 R_s', '14 R_s', '15 R_s')

clear j i sni is ss mindp n i1