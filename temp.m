L=1e10/1.4^1.1;  %% Wicks 2010
fit=[1e-2 1e-1];
load Wicks_100-150
for i=1:18; k_perp(i,:)=fitkappa(2*pi*f/v*rhoi,squeeze(PThetaXX(:,i)+PThetaYY(:,i)),fit); end
for i=1:9; [siperp(i,1),siperp(i,2)]=wmean([k_perp(i,:);k_perp(19-i,:)]); end
save Wicks_100-150
load Wicks_150-200
for i=1:18; k_perp(i,:)=fitkappa(2*pi*f/v*rhoi,squeeze(PThetaXX(:,i)+PThetaYY(:,i)),fit); end
for i=1:9; [siperp(i,1),siperp(i,2)]=wmean([k_perp(i,:);k_perp(19-i,:)]); end
save Wicks_150-200
load Wicks_200-250
for i=1:18; k_perp(i,:)=fitkappa(2*pi*f/v*rhoi,squeeze(PThetaXX(:,i)+PThetaYY(:,i)),fit); end
for i=1:9; [siperp(i,1),siperp(i,2)]=wmean([k_perp(i,:);k_perp(19-i,:)]); end
save Wicks_200-250
load Wicks_250-300
for i=1:18; k_perp(i,:)=fitkappa(2*pi*f/v*rhoi,squeeze(PThetaXX(:,i)+PThetaYY(:,i)),fit); end
for i=1:9; [siperp(i,1),siperp(i,2)]=wmean([k_perp(i,:);k_perp(19-i,:)]); end
save Wicks_250-300
load Wicks_300-350
for i=1:18; k_perp(i,:)=fitkappa(2*pi*f/v*rhoi,squeeze(PThetaXX(:,i)+PThetaYY(:,i)),fit); end
for i=1:9; [siperp(i,1),siperp(i,2)]=wmean([k_perp(i,:);k_perp(19-i,:)]); end
save Wicks_300-350





%% Plot roots of kp12 for frequencies 0.01, 0.1, 1 & 10 Hz
figure
f=[0.1 0.3];
theta=logspace(0,log10(89),100);
v=6e5;
L=1e9;
rhoi=1e5;
% a=2/3; c=L^(a-1); %MHD
a=1/3; c=L^(a-2/3)*rhoi^(-1/3); %KAW


for i=1:2
    kp0=2*pi*f(i)/v./sind(theta);
    kr_gt_1=find(kp0*rhoi>=0);
    kp12_taylor=[kp0-c.*kp0.^a.*cotd(theta); kp0+c.*kp0.^a.*cotd(theta)];
    switch a
        case 2/3
            for m=1:length(theta)
                x1(:,m)=roots([1,c.*cotd(theta(m)),0,-kp0(m)]);
                x2(:,m)=roots([1,-c.*cotd(theta(m)),0,-kp0(m)]);
            end
%             x1=p3root(1,c.*cotd(theta),0,-kp0);
%             x2=p3root(1,-c.*cotd(theta),0,-kp0);
        case 1/3
            for m=1:length(theta)
                x1(:,m)=roots([1,0,c.*cotd(theta(m)),-kp0(m)]);
                x2(:,m)=roots([1,0,-c.*cotd(theta(m)),-kp0(m)]);
            end
%             x1=p3root(1,0,c.*cotd(theta),-kp0);
%             x2=p3root(1,0,-c.*cotd(theta),-kp0);
    end
    for k=1:length(theta)
        [tmp,j]=min(abs(imag(x1(:,k))));
        kp1(k)=real(x1(j,k)).^3;
        [tmp,j]=min(abs(imag(x2(:,k))));
        kp2(k)=real(x2(j,k)).^3;
    end
    subplot(1,2,i)
    semilogy(theta(kr_gt_1),kp0(kr_gt_1)*rhoi,'-red')
    hold all
    plot(theta(kr_gt_1),kp12_taylor(1,kr_gt_1)*rhoi,'--red','Linew',3)
    plot(theta(kr_gt_1),kp1(kr_gt_1)*rhoi,'-black')
    if i==1; legend('k_\perp', 'k_{\perp,1/2} from Taylor', 'k_{\perp,1/2} (Cardano''s method)'); end
    plot(theta(kr_gt_1),kp12_taylor(2,kr_gt_1)*rhoi,'--red','Linew',3)
    plot(theta(kr_gt_1),kp2(kr_gt_1)*rhoi,'-black')
    title(['f=' mat2str(f(i)) 'Hz'])
    xlabel('\theta [^\circ]')
    ylabel('k_\perp [m^{-1}]')
end
clear tmp i j k


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PLOT ROOTS



c=1;
dt=0.14;
for i=1:1136;
    n=find(acf(i,c,:)<exp(-1),1,'first');
    for j=1:n
        tmp(j)=sum(lag(1:end-j).*squeeze(acf_inc(i,c,1:end-j)));
    end
    tc(i)=1-dt*sum(tmp);
end


L=50; %change this for smooth-ness
t=.01:.01:10;
% y=sin(2*pi*t)+rand(size(t));
win=hanning(L)/sum(hanning(L));
y2=conv(win,y);
y2=y2( L/2 +1 :length(y)+L/2);





%% Versuch Q_turb zu berechnen (bzw. erst mal db^4/B_0
clear all

Orbits={[04182 04183]; %Nr.1
        [04301 04302 04303]; %Nr.2
        [04349 04350 04351]; %Nr.3
        [05015 05016 05017]; %Nr.4
        [05046 05047 05048 05049]; %Nr.5
        [05067 05068 05069]; %Nr.6
        [05087 05088 05089 05090]; %Nr.7
        };
for o=1:length(Orbits)
    files=getfiles('10m*dat',Orbits{o});
    for i=1:length(files)
        [ts,stats]=importts(files{i});
        B0{o}(i)=stats.meanB;
        d{o}(i)=stats.dist;
        b=sqrt(ts(:,2).^2+ts(:,3).^2+ts(:,4).^2);
        b=b-B0{o}(i);
        b4{o}(i)=mean(b.^4);
        clear b
    end;
end
    
    
    
%% Versuch von minimum der ACF
clear all
root='/afs/geo/usr/vonpapen/PhD/Cassini/2004/298-304/';
f=dir([root, '10m*dat']);
[ts,stats]=importts([root, f(77).name]);
maxlag=2000;

for n=1:200
    for i=1:length(ts)-n;
        x(i,n)=ts(i+n,2)-ts(i,2);
        y(i,n)=ts(i+n,3)-ts(i,3);
        z(i,n)=ts(i+n,4)-ts(i,4);
    end
    [cx(:,n), lags]=xcov(x(:,n), maxlag, 'unbiased');
    [cy(:,n), lags]=xcov(y(:,n), maxlag, 'unbiased');
    [cz(:,n), lags]=xcov(z(:,n), maxlag, 'unbiased');
    idx=find(cx(maxlag:end,n)==min(cx(maxlag:end,n)));
    mx(n,:)=[idx(1) cx(maxlag+idx(1)-1,n)];
    idy=find(cy(maxlag:end,n)==min(cy(maxlag:end,n)));
    my(n,:)=[idy(1) cy(maxlag+idy(1)-1,n)];
    idz=find(cz(maxlag:end,n)==min(cz(maxlag:end,n)));
    mz(n,:)=[idz(1) cz(maxlag+idz(1)-1,n)];
end

loglog(mx(:,1),mx(:,2),my(:,1),my(:,2),mz(:,1),mz(:,2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
clear all
for n=1:5
switch n
    case 1
        root='/afs/geo/usr/vonpapen/PhD/Cassini/2004/181-184/';
    case 2
        root='/afs/geo/usr/vonpapen/PhD/Cassini/2004/298-304/';
    case 3
        root='/afs/geo/usr/vonpapen/PhD/Cassini/2004/346-353/';
    case 4
        root='/afs/geo/usr/vonpapen/PhD/Cassini/2005/015-017/';
    case 5
        root='/afs/geo/usr/vonpapen/PhD/Cassini/2005/046-049/';
end
f=dir([root, '10m*.dat']);
for i=1:length(f)
    files{i}=[root, f(i).name];
end
M{n}=moments(files,mat2str(n));
clear files
end



%% Hier wird gecheckt ob das kompensierte PSD einigermaßen waagerecht ist
clear all
root='/afs/geo/usr/vonpapen/PhD/Cassini/2004/181-184/';
f=dir([root, '10m*.wtf']);
for i=1:length(f)
    files{i}=[root, f(i).name];
end
[PSD,n]=stack_wtf('10m',[],[6.5 15],[1e-2 1], files, [0.01 1]);

[PSD,files]=stack_wtf('10m', dnum, [6.5 15],[1e-2 0.5],[],[0.01 1]);
[PSDno,files]=stack_wtf('10m', dnum, [6.5 15],[1e-2 0.5],[],0);

figure
for i=1:3
    subplot(2,3,i), loglog(PSD.f,squeeze(pno(:,i,:)),'o','Markersize', 0.1);
    hold all, loglog(PSD.f,pmno(:,i),'-black','Linewidth', 2);
    xlim([0.004 4]), ylim([1E-8 1E-2]), grid on, title('Compensated PSD (\kappa=-2.5)')
    
    
end
















%%%%%%%%%%%%%%%   SI plot erstellen [SI(f)]
clear ch win j windows x siglog ylog PSD
for ch=1:4
        % Fitting
        win=[5 400 5];
        for j=1:size(all.mean,3)-1
            for i=win(1):win(2);
                windows=[i:win(3)*i];
                x = log10(all.mean(windows,1,1));
                siglog=ones(length(windows),1);
                ylog=log10(all.mean(windows,ch,j+1));
                [PSD.sepfit.ab(i,ch,j,:), PSD.sepfit.chi(i,ch,j), PSD.sepfit.var(i,ch,j,:)]...
                    =lsline(x,ylog,siglog);
            end
        end
        for i=win(1):win(2);
            PSD.sepfit.mean(i,ch,:)=mean(PSD.sepfit.ab(i,ch,:,:),3);
            PSD.sepfit.std(i,ch,:)=std(PSD.sepfit.ab(i,ch,:,:),1,3);
        end;
        clear k i j
        clear windows i tmp siglog ylog x chi vars
end  
label={'B_x', 'B_y', 'B_z'}
for ch=2:4
    subplot(3,1,ch-1)
    hold off
%     plot(win(3)/2*all.mean(win(1):win(2),1,1), sqrt(mean(PSD.sepfit.chi(win(1):win(2),ch,:),3)))
    errorbar(win(3)/2*all.mean(win(1):win(2),1,1), PSD.sepfit.mean(win(1):win(2),ch,1),...
                PSD.sepfit.std(win(1):win(2),ch,1), 'Linewidth', 0.5)
    title(['Spectral Index of ', label{ch-1}])
    if ch==4
        xlabel('Center Frequency in Hz')
    end
    ylabel('\kappa')
%     ylabel('\chi')
    hold all        
    set(gca,'XScale', 'log')
    xlim([0.02 2])
    ylim([-4 0])
    grid on
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure mit labeln versehen
label={'B_x', 'B_y', 'B_z'};
xlab='Zeit in s';
ylab=' in nT';
for i=1:3
    subplot(3,1,i)
    if i==1
        title('Inkrementzeitreihen', 'Fontsize', 20)
    end
    xlabel(xlab), grid on
    ylabel([label{i} ylab])
    set(gca, 'Fontsize', 15)
end

%%%%%%%%%%%%%%%% fits erstellen
fitx1=[find(all.mean(:,1,1)<0.001, 1, 'last') find(all.mean(:,1,1)<0.02, 1, 'last')];
fitx2=[find(all.mean(:,1,1)<0.05, 1, 'last') find(all.mean(:,1,1)<0.8, 1, 'last')];
% fitz=[find(all.mean(:,1,1)<0.02, 1, 'last') find(all.mean(:,1,1)<0.1, 1, 'last')];
x1=log10(all.mean(fitx1(1):fitx1(2),1,1));
x2=log10(all.mean(fitx2(1):fitx2(2),1,1));
% z=log10(all.mean(fitz(1):fitz(2),1,1));

siglog1=ones(length(x1),1);
siglog2=ones(length(x2),1);
% siglogz=ones(length(z),1);
for ch=1:3
    Pm(:,ch)=10.^mean(log10(all.mean(:,ch+1,3:9)),3);
end
for ch=1:3;
    y=log10(Pm(fitx1(1):fitx1(2),ch));
    [ab(1,ch,:), chi(1,ch), var(1,ch,:)]=lsline(x1,y,siglog1);
end
for ch=1:3;
    y=log10(Pm(fitx2(1):fitx2(2),ch));
    [ab(2,ch,:), chi(2,ch), var(2,ch,:)]=lsline(x2,y,siglog2);
end

% y=log10(Pm(fitz(1):fitz(2),ch));
% [ab(1,ch,:), chi(1,ch), var(1,ch,:)]=lsline(z,y,siglogz);


% figure

%%%%%%%%%%%%%%%%  Fits in die allPSD plotten
for ch=1:3
    subplot(1,3,ch)
%     if ch<3
        loglog(all.mean(:,1,1), Pm(:,ch))
        xlim([limx])
        ylim([limy])
        grid on
        hold all
        plot([all.mean(fitx1(1),1,1) all.mean(fitx1(2),1,1)],...
            5.*10.^[(ab(1,ch,2)+log10(all.mean(fitx1(1),1,1))*ab(1,ch,1))...
            (ab(1,ch,2)+log10(all.mean(fitx1(2),1,1))*ab(1,ch,1))],...
            '--black', 'linew', 3)
        text(all.mean(fitx1(2),1,1), 5.*10.^(ab(1,ch,2)+log10(all.mean(fitx1(2),1,1))*ab(1,ch,1)),...
            mat2str(round(100*ab(1,ch,1))/100),...
            'Backgr', 'white')
        plot([all.mean(fitx2(1),1,1) all.mean(fitx2(2),1,1)],...
            5.*10.^[(ab(2,ch,2)+log10(all.mean(fitx2(1),1,1))*ab(2,ch,1))...
            (ab(2,ch,2)+log10(all.mean(fitx2(2),1,1))*ab(2,ch,1))],...
            '--black', 'linew', 3)
        text(all.mean(fitx2(2),1,1), 5.*10.^(ab(2,ch,2)+log10(all.mean(fitx2(2),1,1))*ab(2,ch,1)),...
            mat2str(round(100*ab(2,ch,1))/100),...
            'Backgr', 'white')
%     else
%         plot([all.mean(fitz(1),1,1) all.mean(fitz(2),1,1)],...
%             5.*10.^[(ab(2,ch,2)+log10(all.mean(fitz(1),1,1))*ab(2,ch,1))...
%             (ab(2,ch,2)+log10(all.mean(fitz(2),1,1))*ab(2,ch,1))],...
%             '--black', 'linew', 3)
%         text(all.mean(fitz(2),1,1), 5.*10.^(ab(1,ch,2)+log10(all.mean(fitz(2),1,1))*ab(1,ch,1)),...
%             mat2str(round(100*ab(1,ch,1))/100),...
%             'Backgr', 'white')
%     end
end





        
        
filename='Oa_files_completeOa2.mat';

db_thres=[9e9 0.4 0.2]


load(filename); clear filename
Oa2_complete=files; clear files
n=length(Oa2_complete);

%%CALC \Deltab_max
db=1e-10*ones(n,3);
for i=1:n
    for j=1:length(ts{i})-1;
        for ch=1:3;
            db(i,ch)=max(abs(ts{i}(j,ch+1)-ts{i}(j+1,ch+1)),db(i,ch));
        end;
    end;
end; clear i j ch

for i=1:length(PSD.stats);
    B0(i)=PSD.stats{i}.meanB;
    rms(i,:)=PSD.stats{i}.dB2;
    rtmp(i)=PSD.stats{i}.dist;
    ttmp(i)=datenum(PSD.stats{i}.utc{1}, 'yyyy-mm-ddTHH:MM:SS');
end; clear i

[r, r_ix]=sort(rtmp); clear rtmp
[t, t_ix]=sort(ttmp); clear ttmp

for i=1:length(db_thres)
    clear Oa2_filt
    j=0;
    for k=1:length(db)
        if max(db(k,:))<db_thres(i)
            j=j+1;
            Oa2_filt{j}=Oa2_complete{k};
        end
    end
    [PSD_db{i}, files_db{i}]=stack('10m', [], [6.5 15], 0, Oa2_filt);
    PSD_db{i}.db_thres=db_thres(i);
    files_db{i}.db_thres=db_thres(i);
end



%% PLOT FIGURE
fig=figure('PaperType', 'A4', 'PaperOrientation', 'portrait',...
            'PaperUnits', 'normalized', 'PaperPosition', [0.1 0.1 0.9 0.9], 'PaperPositionMode', 'manual')
subplot(4,1,1)
    plot(t, db(t_ix,:)), hold all
    title('\Deltab_{max} of DOY183')
    for i=1:length(db_thres)
        plot(t, ones(size(t))*db_thres(i), '--black', 'linew', 1.5)
    end
    ylim([0 1]),
    ylabel('\Deltab_{max}=max(|b_{i+1}-b_i|)');
    xlabel('Time')
    datetick_doy('x', 15) %15=HH:MM
    legend('x', 'y', 'z')
    grid on
subplot(4,1,2)
    plot(r, db(r_ix,:)), hold all
    ylim([0 1]),
    xlabel('Distance in R_S')
    ylabel('\Deltab_{max}=max(|b_{i+1}-b_i|)');
    grid on
    for i=1:length(db_thres)
        plot(r, ones(size(r))*db_thres(i), '--black', 'linew', 1.5)
    end

for i=1:length(db_thres)
    subplot(4,length(db_thres),[2*length(db_thres)+i 3*length(db_thres)+i])
    x=[0.001, 4];
    y=[1e-6, 2];
    loglog(PSD.f, PSD_db{i}.mean), hold all
    loglog(x, [10^(log10(1e-4*y(2))-5/3*log10(x(1))), 10^(log10(1e-4*y(2))-5/3*log10(x(2)))], '--black', 'linew', 1)
    xlim(x)
    ylim(y)
    grid on
    if db_thres(i)<1
        title(['PSD with \Deltab_{max}<', mat2str(db_thres(i)), ' nT'])
    else
        title('PSD with \Deltab_{max} free')
    end
    xlabel('f in Hz')
    set(gca, 'XTick', [1E-3, 1E-2, 1E-1, 1E0])
    text(0.05, 0.06, ['#stacks: ', mat2str(length(PSD_db{i}.stats))], 'Units', 'normalized')
    if i==1
        ylabel('PSD in nT^2')
        legend('x', 'y', 'z')
    end
end
saveas(fig, 'pics/183_db.eps', 'psc2')
close(fig)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














% % clear all
% 
% %% set
% dnum(1)=datenum([2004, 6, 30, 15, 19, 0]);
% dnum(2)=datenum([2004, 6, 30, 21, 0, 0]);
% 
% %% load
% [PSD, files]=stack('10m', dnum, [6.5 15])
% for i=1:length(PSD.stats); B0(i)=PSD.stats{i}.meanB; b2(i,:)=PSD.stats{i}.dB2; r(i)=PSD.stats{i}.dist; tnum(i)=datenum(PSD.stats{i}.utc{1}, 'yyyy-mm-ddTHH:MM:SS'); end
% [time, ind]=sort(tnum)
% 
% %% plot
% % close all
% subplot(3,1,1), semilogy(time, B0(ind), '-o'), hold all, grid on
% subplot(3,1,2), semilogy(time, b2(ind,:), '-o'), hold all, grid on
% subplot(3,1,3), semilogy(time, r(ind), '-o'), hold all, grid on
% subplot(3,1,1), datetick_doy('x', 36)
% subplot(3,1,2), datetick_doy('x', 36)
% subplot(3,1,3), datetick_doy('x', 36)




[PSD, files]=stack('10m', [], [6.5 15], 0, Oa2)
for i=1:length(Oa2);
    [ts{i}, stats{i}]=importts(strrep(files{i}, 'psd', 'dat'));
end

%%CALC dB
% Oa2=files;
db=1e-10*ones(length(Oa2),3);
for i=1:length(Oa2)
    for j=1:length(ts{i})-1;
        for ch=1:3;
            db(i,ch)=max(abs(ts{i}(j,ch+1)-ts{i}(j+1,ch+1)),db(i,ch));
        end;
    end;
end

save('Oa_files_Oa2_sort.mat', 'PSD', 'files')


%% db THreshold
dbthres=0.2;
clear Oa2_filt
n=0;
for i=1:length(db)
    if max(db(i,:))<dbthres
        n=n+1;
        Oa2_filt{n}=Oa2{i};
    end
end
[PSD, files]=stack('10m', [], [6.5 15], 0.01, Oa2_filt)

%% PLOT
x=[0.001, 4];
y=[1e-6, 2];
loglog(PSD.f, PSD.mean, x, [10^(log10(1e-4*y(2))-5/3*log10(x(1))), 10^(log10(1e-4*y(2))-5/3*log10(x(2)))], '--black', 'linew', 1)
xlim(x)
ylim(y)
grid on
title('PSD of DOY183 (\Deltab_{max}<0.2 nT)')
xlabel('f in Hz')
% ylabel('PSD in nT^2')
legend('x', 'y', 'z')


%% Variance of b threshold
rmsthres=1;
clear Oa2_filt
n=0;
for i=1:length(ts)
    for ch=1:3
        rms(i,ch)=std(ts{i}(:,ch+1));
    end
    if max(rms(i,:))<rmsbthres
        n=n+1;
        Oa2_filt{n}=Oa2{i};
    end
end
[PSD, files]=stack('10m', [], [6.5 15], 0.01, Oa2_filt)

% %% PLOT
% 
% figure
% x=[0.008, 20]; y=[5e-7, 1];
% loglog(PSD_db01.bin.f, PSD_db01.bin.mean, x, [10^(log10(1e-4*y(2))-5/3*log10(x(1))), 10^(log10(1e-4*y(2))-5/3*log10(x(2)))], '--black', 'linew', 1)
% xlim(x)
% ylim(y)
% grid on
% title('PSD of DOY183 (6.5R_S<r<15R_S)')
% xlabel('f in Hz')
% ylabel('PSD in nT^2')
% legend('x', 'y', 'z')

figure
x=[0.001, 4];
y=[3e-7, 1];
loglog(PSD.f, PSD.mean, x, [10^(log10(1e-4*y(2))-5/3*log10(x(1))), 10^(log10(1e-4*y(2))-5/3*log10(x(2)))], '--black', 'linew', 1)
xlim(x)
ylim(y)
grid on
title(['PSD of DOY183 (RMS<',mat2str(sigbthres), ' nT)'])
xlabel('f in Hz')
ylabel('PSD in nT^2')
legend('x', 'y', 'z')






















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Läd mat files ein, die filenamen eines Zeitraumes enthalten, der genauer
%% untersucht werden soll. Dazu werden die Parameter \deltab, rms^2 und
%% rms^4 variiert und die entwicklung des power spektrums geplottet

clear all

DOY=182;

norm=0;

%% FILENAME % THRESHOLDS
switch DOY
    case 182
        filename='/afs/geo/usr/vonpapen/PhD/Cassini/2004/181-184/Oa1_sort.mat';
        limy1=[0 0.5];
        limy2=[1E-4 1];
        limy3=[1E-8 1];
        db_thres=[0.9];
        rms2_thres=[1E-3];
        rms4_thres=[1E-6];
        range=[6.5 11];
    case 183
        filename='/afs/geo/usr/vonpapen/PhD/Cassini/2004/181-184/Oa2_complete.mat';
        limy1=[0 0.5];
        limy2=[1E-4 1E0];
        limy3=[1E-8 1E0];
        db_thres=[0.9];
        rms2_thres=[1E-3];
        rms4_thres=[1E-6];
        range=[6.5 12];
end

norm=1;

%% LOADING
load(filename); clear filename
complete=files; clear files
n=length(complete);
[PSD, files]=stack('10m', [], [6.5 15], 0, complete, norm);
for i=1:length(complete);
    [ts{i}, stats{i}]=importts(strrep(complete{i}, 'psd', 'dat'));
end

clear files


%%CALC \Deltab_max
db=1e-10*ones(n,3);
for i=1:n
    for j=1:length(ts{i})-1;
        for ch=1:3;
            db(i,ch)=max(abs(ts{i}(j,ch+1)-ts{i}(j+1,ch+1)),db(i,ch));
        end;
    end;
end; clear i j ch

for i=1:length(PSD.stats);
    B0(i)=PSD.stats{i}.meanB;
    rms(i,:)=sqrt(PSD.stats{i}.dB2);
    rtmp(i)=PSD.stats{i}.dist;
    ttmp(i)=datenum(PSD.stats{i}.utc{1}, 'yyyy-mm-ddTHH:MM:SS');
end; clear i

[r, r_ix]=sort(rtmp); clear rtmp
[t, t_ix]=sort(ttmp); clear ttmp

[PSD, files]=stack('10m', [], range, rms2_thres, complete, norm);
 
    if norm==0
        limy=[1E-6 1];
    else
        limy=[5E-4 3E3];
    end

    %% PLOT FIGURE
    fig=figure('PaperType', 'A4', 'PaperOrientation', 'portrait',...
                'PaperUnits', 'normalized', 'PaperPosition', [0.05 0.15 0.9 0.8], 'PaperPositionMode', 'manual')
    subplot(5,1,1), hold off
        plot(r, db(r_ix,[3,2,1])), hold all
        plot(ones(2,1)*range(2), limy1, '--black', 'linew', 3)
        title(['\Deltab_{max} of DOY', mat2str(DOY), ' sorted by distance (R_S):'])
        for i=1:length(db_thres)
            plot(r, ones(size(r))*db_thres(i), '--black', 'linew', 1.5)
        end
        ylim(limy1)
        ylabel('\Deltab_{max} in nT');
        legend('z', 'y', 'x', 'Location', 'NE')
        grid on
    subplot(5,1,2), hold off
        plot(r, rms(r_ix,[3,2,1]).^2), hold all
        plot(ones(2,1)*range(2), limy2, '--black', 'linew', 3)
        ylim(limy2)
        set(gca, 'YScale', 'log', 'YTick', logspace(-4,2,7))
        title('Energy <b^2> (RMS^2)')
        ylabel('<b^2> in nT^2');
        grid on
        for i=1:length(rms2_thres)
            plot(r, ones(size(r))*rms2_thres(i), '--black', 'linew', 1.5)
        end
    subplot(5,1,3), hold off
        plot(r, rms(r_ix,[3,2,1]).^4), hold all
        plot(ones(2,1)*range(2), limy3, '--black', 'linew', 3)
        ylim(limy3)
        set(gca, 'YScale', 'log', 'YTick', logspace(-10,2,7))
        title('Dissipation <b^2>^2 (RMS^4)')
        ylabel('<b^2>^2 in nT^4');
        grid on
        for i=1:length(rms4_thres)
            plot(r, ones(size(r))*rms4_thres(i), '--black', 'linew', 1.5)
        end
    subplot(5,2,[7 9]), hold off
        x=[0.001, 4];
        y=[limy];
        loglog(PSD.f, PSD.mean(:,[1, 2])), hold all
        loglog(x, [10^(log10(1e-4*y(2))-5/3*log10(x(1))), 10^(log10(1e-4*y(2))-5/3*log10(x(2)))], '--black', 'linew', 1)
        xlim(x)
        ylim(y)
        grid on
        title(['parallel PSD with \Deltab<', mat2str(db_thres), 'nT and RMS^2>', mat2str(rms2_thres), ' nT^2'])
        text(0.1, 4*10^(log10(1e-4*y(2))-5/3*log10(0.1)), '-5/3', 'Backgr', 'white', 'rotation', -50)
        text(0.05, 0.06, [' # stacks: ', mat2str(length(PSD.stats))], 'Units', 'normalized', 'Backgr', 'white')
        ylabel('normalized PSD')
        legend('x', 'y')
        xlabel('f in Hz')
        set(gca, 'XTick', [1E-3, 1E-2, 1E-1, 1E0])
    subplot(5,2,[8 10]), hold off
        x=[0.001, 4];
        y=[limy];
        loglog(PSD.f, PSD.mean(:,3)), hold all
        loglog(x, [10^(log10(1e-4*y(2))-5/3*log10(x(1))), 10^(log10(1e-4*y(2))-5/3*log10(x(2)))], '--black', 'linew', 1)
        xlim(x)
        ylim(y)
        grid on
        title(['perp. PSD with \Deltab<', mat2str(db_thres), 'nT and RMS^2>', mat2str(rms2_thres), ' nT^2'])
        text(0.1, 4*10^(log10(1e-4*y(2))-5/3*log10(0.1)), '-5/3', 'Backgr', 'white', 'rotation', -50)
        text(0.05, 0.06, [' # stacks: ', mat2str(length(PSD.stats))], 'Units', 'normalized', 'Backgr', 'white')
        ylabel('normalized PSD')
        xlabel('f in Hz')
        set(gca, 'XTick', [1E-3, 1E-2, 1E-1, 1E0])        

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(fig, ['Oa', mat2str(DOY), '_stack.eps'], 'psc2')








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hier sollen die Spektralen Indizes bestimmt werden
%% PSD's und files müssen schon geladen sein
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fit ranges
f1=[3:42]; %0.05-0.07
f2=[60:240]; %0.1-0.4

%% Fitting
for i=1:length(files)
    for ch=1:3 %% Fit f1
        x = log10(PSD.f(f1));
        siglog=ones(length(f1),1);
        ylog=log10(PSD.single(f1,ch,i));
        [PSD.fit1.ab(:,ch,i), PSD.fit1.chi1(ch,i), PSD.fit1.var(:,ch,i)]...
                    =lsline(x,ylog,siglog);
        clear x ylog
    end
    for ch=1:3 %% Fit f2
        x = log10(PSD.f(f2));
        siglog=ones(length(f2),1);
        ylog=log10(PSD.single(f2,ch,i));
        [PSD.fit2.ab(:,ch,i), PSD.fit2.chi1(ch,i), PSD.fit2.var(:,ch,i)]...
                    =lsline(x,ylog,siglog);
        clear x ylog
    end
end
clear siglog ch i
label={'x', 'y', 'z'}
plot(r,squeeze(PSD.fit1.ab(1,:,r_ix)))
legend(label)
i=1
%% Plotting
label={'x', 'y', 'z'}
for ch=1:3
    subplot(1,3,ch)
    loglog(PSD.f, PSD.single(:,ch,i))
    hold all
    plot([PSD.f(f1(1)), PSD.f(f1(end))],...
        [10^(PSD.fit1.ab(2,ch,i)+log10(PSD.f(f1(1)))*PSD.fit1.ab(1,ch,i)),...
        10^(PSD.fit1.ab(2,ch,i)+log10(PSD.f(f1(end)))*PSD.fit1.ab(1,ch,i))],...
        '-r', 'Linewidth', 3);
    plot([PSD.f(f2(1)), PSD.f(f2(end))],...
        [10^(PSD.fit2.ab(2,ch,i)+log10(PSD.f(f2(1)))*PSD.fit2.ab(1,ch,i)),...
        10^(PSD.fit2.ab(2,ch,i)+log10(PSD.f(f2(end)))*PSD.fit2.ab(1,ch,i))],...
        '-g', 'Linewidth', 3);
    legend(cell2mat(label(ch)), ['slope=', mat2str(round(100*PSD.fit1.ab(1,ch,i))/100)],...
        ['slope=', mat2str(round(100*PSD.fit2.ab(1,ch,i))/100)], 'Location', 'SouthWest');
    hold off
end