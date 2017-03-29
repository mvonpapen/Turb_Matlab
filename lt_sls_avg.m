%% Compute radially binned averages for LT and SLS

% qqplot(log10(sqrt(sum(stdB(i,:).^2,2))),sqrt(sum(stdB(i,:).^2,2)))


%% LT means
for i=1:24; 
    
    j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./scaleheight<1 ...
        & ltime>=i-1 & ltime<=i & B>=5 & r_cyl>=5 & r_cyl<=10);
    rmsLTmean_Rle10(i,:)=[mean(log10(sqrt(sum(stdB(j,:).^2,2)))) std(log10(sqrt(sum(stdB(j,:).^2,2))))];
    QL_LTmean_Rle10(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
    E_LTmean_Rle10(i,:)=[nanmean(log10(Ekprho(j))) nanstd(log10(Ekprho(j)))];

    j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./scaleheight<1 ...
        & ltime>=i-1 & ltime<=i & B>=5 & r_cyl>=15 & r_cyl<=20);
    rmsLTmean_Rge10(i,:)=[mean(log10(sqrt(sum(stdB(j,:).^2,2)))) std(log10(sqrt(sum(stdB(j,:).^2,2))))];
    QL_LTmean_Rge10(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
    E_LTmean_Rge10(i,:)=[nanmean(log10(Ekprho(j))) nanstd(log10(Ekprho(j)))];

end

%% SLS means of rms
for i=1:36;
    
    j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./scaleheight<1 ...
        & SLS4>=(i-1)*10 & SLS4<=i*10 & B>=5 & r_cyl>=5 & r_cyl<=10);
    rmsSLSmean_Rle10(i,:)=[mean(log10(sqrt(sum(stdB(j,:).^2,2)))) std(log10(sqrt(sum(stdB(j,:).^2,2))))];
    QL_SLSmean_Rle10(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
    E_SLSmean_Rle10(i,:)=[nanmean(log10(Ekprho(j))) nanstd(log10(Ekprho(j)))];
  
    j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./scaleheight<1 ...
        & SLS4>=(i-1)*10 & SLS4<=i*10 & B>=5 & r_cyl>=15 & r_cyl<=20);
    rmsSLSmean_Rge10(i,:)=[mean(log10(sqrt(sum(stdB(j,:).^2,2)))) std(log10(sqrt(sum(stdB(j,:).^2,2))))];
    QL_SLSmean_Rge10(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
    E_SLSmean_Rge10(i,:)=[nanmean(log10(Ekprho(j))) nanstd(log10(Ekprho(j)))];
    
    j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./scaleheight<1 ...
        & SLS4>=(i-1)*10 & SLS4<=i*10 & B>=5 & r_cyl>=10 & r_cyl<=15 & (ltime<=6 | ltime>=18));
    rmsSLSmean_Rge10_night(i,:)=[mean(log10(sqrt(sum(stdB(j,:).^2,2)))) std(log10(sqrt(sum(stdB(j,:).^2,2))))];
    QL_SLSmean_Rge10_night(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
    E_SLSmean_Rge10_night(i,:)=[nanmean(log10(Ekprho(j))) nanstd(log10(Ekprho(j)))];
    
    j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./scaleheight<1 ...
        & SLS4>=(i-1)*10 & SLS4<=i*10 & B>=5 & r_cyl>=10 & r_cyl<=15 & ltime>=6 & ltime<=18);
    rmsSLSmean_Rge10_day(i,:)=[mean(log10(sqrt(sum(stdB(j,:).^2,2)))) std(log10(sqrt(sum(stdB(j,:).^2,2))))];
    QL_SLSmean_Rge10_day(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
    E_SLSmean_Rge10_day(i,:)=[nanmean(log10(Ekprho(j))) nanstd(log10(Ekprho(j)))];

end


% Plot means of rms
figure
subplot(3,1,1)
errorbar([0.5:23.5],rmsLTmean_Rle10(:,1),rmsLTmean_Rle10(:,2),'--+black')
hold all, ylabel('log_1_0(rms [nT])'), xlabel('Local time [h]')
errorbar([0.5:23.5],rmsLTmean_Rge10(:,1),rmsLTmean_Rge10(:,2),'--+red')
legend('r \leq 10 Rs', 'r \geq 10 Rs')
title('Energy in magnetic field fluctuations')
xlim([0 24]), ylim([-0.7 0.2])

subplot(3,1,2)
errorbar([5:10:360],rmsSLSmean_Rle10(:,1),rmsSLSmean_Rle10(:,2),'--+black')
hold all, ylabel('log_1_0(rms [nT])'), xlabel('SLS4')
errorbar([5:10:360],rmsSLSmean_Rge10(:,1),rmsSLSmean_Rge10(:,2),'--+red')
legend('r \leq 10 Rs', 'r \geq 10 Rs')
xlim([0 360]), ylim([-0.7 0.2])
plot([0:360],cosd([0:360]-300)*0.04-0.3)

subplot(3,1,3)
errorbar([5:10:360],rmsSLSmean_Rge10_day(:,1),rmsSLSmean_Rge10_day(:,2),'--+black')
hold all, ylabel('log_1_0(rms [nT])'), xlabel('SLS4')
errorbar([5:10:360],rmsSLSmean_Rge10_night(:,1),rmsSLSmean_Rge10_night(:,2),'--+red')
legend('r \geq 10 Rs, day', 'r \geq 10 Rs, night')
xlim([0 360]), ylim([-0.7 0.2])
plot([0:360],cosd([0:360]-300)*0.04-0.3)



%% Plot means of QL
% figure
% subplot(3,1,1)
% errorbar([0.5:23.5],QL_LTmean_Rle10(:,1),QL_LTmean_Rle10(:,2),'--+black')
% hold all, ylabel('log_1_0(Q_L [W/m^3])'), xlabel('Local time [h]')
% errorbar([0.5:23.5],QL_LTmean_Rge10(:,1),QL_LTmean_Rge10(:,2),'--+red')
% legend('r \leq 10 Rs', 'r \geq 10 Rs')
% title('Heating rate')
% xlim([0 24]), ylim([-17 -14])
% 
% subplot(3,1,2)
% errorbar([5:10:360],QL_SLSmean_Rle10(:,1),QL_SLSmean_Rle10(:,2),'--+black')
% hold all, ylabel('log_1_0(Q_L [W/m^3])'), xlabel('SLS4')
% errorbar([5:10:360],QL_SLSmean_Rge10(:,1),QL_SLSmean_Rge10(:,2),'--+red')
% legend('r \leq 10 Rs', 'r \geq 10 Rs')
% xlim([0 360]), ylim([-17 -14])
% plot([0:360],cosd([0:360]-300)*0.2-15.4)
% 
% subplot(3,1,3)
% errorbar([5:10:360],QL_SLSmean_Rge10_day(:,1),QL_SLSmean_Rge10_day(:,2),'--+black')
% hold all, ylabel('log_1_0(Q_L [W/m^3])'), xlabel('SLS4')
% errorbar([5:10:360],QL_SLSmean_Rge10_night(:,1),QL_SLSmean_Rge10_night(:,2),'--+red')
% legend('r \geq 10 Rs, day', 'r \geq 10 Rs, night')
% xlim([0 360]), ylim([-17 -14])
% plot([0:360],cosd([0:360]-300)*0.2-15.4)



% %% Plot means of Ekprho
% figure
% subplot(3,1,1)
% errorbar([0.5:23.5],E_LTmean_Rle10(:,1),E_LTmean_Rle10(:,2),'--+black')
% hold all, ylabel('log_1_0(E/E_0)'), xlabel('Local time [h]')
% errorbar([0.5:23.5],E_LTmean_Rge10(:,1),E_LTmean_Rge10(:,2),'--+red')
% legend('r \leq 10 Rs', 'r \geq 10 Rs')
% title('Relative spectral power')
% xlim([0 24]), ylim([-0.2 1.8])
% 
% subplot(3,1,2)
% errorbar([5:10:360],E_SLSmean_Rle10(:,1),E_SLSmean_Rle10(:,2),'--+black')
% hold all, ylabel('log_1_0(E/E_0)'), xlabel('SLS4')
% errorbar([5:10:360],E_SLSmean_Rge10(:,1),E_SLSmean_Rge10(:,2),'--+red')
% legend('r \leq 10 Rs', 'r \geq 10 Rs')
% xlim([0 360]), ylim([-0.2  1.8])
% plot([0:360],cosd([0:360]-300)*0.1+0.8)
% 
% subplot(3,1,3)
% errorbar([5:10:360],E_SLSmean_Rge10_day(:,1),E_SLSmean_Rge10_day(:,2),'--+black')
% hold all, ylabel('log_1_0(E/E_0)'), xlabel('SLS4')
% errorbar([5:10:360],E_SLSmean_Rge10_night(:,1),E_SLSmean_Rge10_night(:,2),'--+red')
% legend('r \geq 10 Rs, day', 'r \geq 10 Rs, night')
% xlim([0 360]), ylim([-0.2  1.8])
% plot([0:360],cosd([0:360]-300)*0.1+0.8)