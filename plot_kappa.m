%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot der Spektralen indizes der einzelnen Orbits als histogramm und linien
clear all

saveto='/home/vonpapen/PhD/ordnung/';
version='04301_2-8h';

range=[14.25 16.95];
str='10m*.wti';
fitrange=[0.04 0.3];
stdthres=[0.01 20];
bin=[-5.5:0.15:1];
% Orbits={[04182 04183]; %Nr.1
%         [04301 04302 04303]; %Nr.2
%         [04349 04350 04351]; %Nr.3
%         [05015 05016 05017]; %Nr.4
%         [05046 05047 05048 05049]; %Nr.5
%         [05067 05068 05069]; %Nr.6
%         [05087 05088 05089 05090]; %Nr.7
%         };
Orbits={[04301]};
       

for i=1:length(Orbits)
    [H,bin,dfit{i}]=get_kappa(fitrange, str ,Orbits{i}, bin, range, stdthres);
    X1(:,i)=H(:,1,1)+H(:,2,1);
    Z1(:,i)=H(:,3,1);
    if size(fitrange,1)==2
        X2(:,i)=H(:,1,2)+H(:,2,2);
        Z2(:,i)=H(:,3,2);
    end
end

xlab1=['\alpha: ' mat2str(fitrange(1,1)) '-' mat2str(fitrange(1,2)) ' Hz'];
ylab='Amount of blocks';
t2='Seperate histograms';
tx1='Perp. spectral indices \alpha_1';
tz1='Parallel spectral indices  \alpha_1';
leg={'Orbit 1', 'Orbit 2','Orbit 3','Orbit 4','Orbit 5','Orbit 6','Orbit 7'}; %'Orbit 1',
legloc=[40 160 30 30];
if size(fitrange,1)==2
    xlab2=['\alpha: ' mat2str(fitrange(2,1)) '-' mat2str(fitrange(2,2)) ' Hz'];
    tx2='Spectral indices of P_\perp';
    tz2='Spectral indices of P_|_|';
end
fs=12; %Fontsize



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Berechnung & Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(fitrange,1)==2
    pt1='A4';
    pt2='landscape';
else
    pt1='A5';
    pt2='portrait';
end

anz=sum(sum(X1,2));

%% Stacked Plot
fig1=figure('PaperType', pt1, 'PaperOrientation', pt2, 'Visi', 'on',...
            'PaperUnits', 'normalized', 'PaperPosition', [0.05 0.1 0.9 0.85],...
            'PaperPositionMode', 'manual');
if size(fitrange,1)==2
    subplot(2,2,1)
else
    subplot(2,1,1)
end
bar(bin,X1,'stack'); hold all
kappa=bin*sum(X1,2)/anz; kstr=sprintf('%0.2f',kappa);
stdkap=sqrt((bin-kappa).^2*sum(X1,2)/anz); stdstr=sprintf('%0.2f',stdkap);
limy=get(gca,'ylim');
text(0.67, 0.6, ['\alpha=',kstr,'\pm',stdstr],...
    'Units','normalized','Backgr','white', 'Fonts', fs)
plot(ones(2,1)*kappa,limy,'--black', 'Linew', 2)
xlim([-4 1]); title(tx1, 'Fonts', fs); xlabel(xlab1, 'Fonts', fs);
ylabel(ylab, 'Fonts', fs),grid
set(gca, 'Fonts', fs)

if size(fitrange,1)==2
    subplot(2,2,2), bar(bin,X2,'stack'); hold all
    kappa=bin*sum(X2,2)/anz; kstr=sprintf('%0.2f',kappa);
    stdkap=sqrt((bin-kappa).^2*sum(X2,2)/anz); stdstr=sprintf('%0.2f',stdkap);
    limy=get(gca,'ylim');
    text(0.67, 0.6, ['\alpha=',kstr,'\pm',stdstr],...
        'Units','normalized','Backgr','white', 'Fonts', fs)
    plot(ones(2,1)*kappa,limy,'--black', 'Linew', 2)
    xlim([-4 1]); title(tx2, 'Fonts', fs); xlabel(xlab2, 'Fonts', fs);
    ylabel(ylab, 'Fonts', fs),grid
    set(gca, 'Fonts', fs)
end

if size(fitrange,1)==2
    subplot(2,2,3)
else
    subplot(2,1,2)
end
bar(bin,Z1,'stack'); legend(leg,legloc,'Fontsize', fs);
hold all
kappa=bin*sum(Z1,2)*2/anz; kstr=sprintf('%0.2f',kappa);
stdkap=sqrt((bin-kappa).^2*sum(Z1,2)*2/anz);
stdstr=sprintf('%0.2f',stdkap);
limy=get(gca,'ylim');
text(0.67, 0.6, ['\alpha=',kstr,'\pm',stdstr],...
    'Units','normalized','Backgr','white', 'Fonts', fs)
plot(ones(2,1)*kappa,limy,'--black', 'Linew', 2)
xlim([-4 1]); title(tz1, 'Fonts', fs); xlabel(xlab1, 'Fonts', fs);
ylabel(ylab, 'Fonts', fs),grid
set(gca, 'Fonts', fs)

if size(fitrange,1)==2
    subplot(2,2,4), bar(bin,Z2,'stack'); hold all
    kappa=bin*sum(Z2,2)*2/anz; kstr=sprintf('%0.2f',kappa);
    stdkap=sqrt((bin-kappa).^2*sum(Z2,2)*2/anz);
    stdstr=sprintf('%0.2f',stdkap);
    limy=get(gca,'ylim');
    text(0.67, 0.6, ['\alpha=',kstr,'\pm',stdstr],...
        'Units','normalized','Backgr','white', 'Fonts', fs)
    plot(ones(2,1)*kappa,limy,'--black', 'Linew', 2)
    xlim([-4 1]); title(tz2, 'Fonts', fs); xlabel(xlab2, 'Fonts', fs);
    ylabel(ylab, 'Fonts', fs),grid
    set(gca, 'Fontsize', fs)
end


%% Seperate Orbits Plot
clear kstr
fig2=figure('PaperType', pt1, 'PaperOrientation', pt2, 'Visi', 'on',...
            'PaperUnits', 'normalized', 'PaperPosition', [0.05 0.1 0.9 0.85],...
            'PaperPositionMode', 'manual'); clear pt1 pt2
if size(fitrange,1)==2
    subplot(2,2,1)
else
    subplot(2,1,1)
end
% plot(bin,X1); hold all;
bar(bin,X1,'stack'); hold all
kappa=bin*sum(X1,2)/anz; kstr=sprintf('%0.2f',kappa);
stdkap=sqrt((bin-kappa).^2*sum(X1,2)/anz); stdstr=sprintf('%0.2f',stdkap);
text(0.67, 0.6, ['\alpha=',kstr,'\pm',stdstr],...
    'Units','normalized','Backgr','white', 'Fonts', fs, 'Fontw', 'bold')
plot(ones(2,1)*kappa,limy,'--black')
kappa=bin*X1./sum(X1); limy=get(gca,'ylim');
for i=1:length(kappa)
    kstr=sprintf('%0.2f',kappa(i));
    stdkap=sqrt((bin-kappa(i)).^2*X1(:,i)/sum(X1(:,i)));
    stdstr=sprintf('%0.2f',stdkap);
    text(0.7, 1-0.11*i, ['\alpha_',num2str(i),'=',kstr,'\pm',stdstr],...
        'Units','normalized','Backgr','white', 'Fonts', fs-2)
%     plot(ones(2,1)*kappa,limy,'--black')
end
xlim([-4 1]); title(tx1, 'Fonts', fs); xlabel(xlab1, 'Fonts', fs);
ylabel(ylab, 'Fonts', fs),grid
set(gca, 'Fonts', fs)

if size(fitrange,1)==2
    subplot(2,2,2), plot(bin,X2); hold all;
    kappa=bin*X2./sum(X2); limy=get(gca,'ylim');
    for i=1:length(kappa)
        kstr=sprintf('%0.2f',kappa(i));
        stdkap=sqrt((bin-kappa(i)).^2*X2(:,i)/sum(X2(:,i)));
        stdstr=sprintf('%0.2f',stdkap);
        text(0.7, 1-0.11*i, ['\alpha_',num2str(i),'=',kstr,'\pm',stdstr],...
            'Units','normalized','Backgr','white', 'Fonts', fs-2)
        plot(ones(2,1)*kappa,limy,'--black')
    end
    xlim([-4 1]); title(tx2, 'Fonts', fs); xlabel(xlab2, 'Fonts', fs);
    ylabel(ylab, 'Fonts', fs),grid
    set(gca, 'Fonts', fs)
end

if size(fitrange,1)==2
    subplot(2,2,3)
else
    subplot(2,1,2)
end
% plot(bin,Z1);
bar(bin,Z1,'stack'); hold all
kappa=bin*sum(Z1,2)*2/anz; kstr=sprintf('%0.2f',kappa);
stdkap=sqrt((bin-kappa).^2*sum(Z1,2)*2/anz); stdstr=sprintf('%0.2f',stdkap);
text(0.67, 0.6, ['\alpha=',kstr,'\pm',stdstr],...
    'Units','normalized','Backgr','white', 'Fonts', fs, 'Fontw', 'bold')
legend(leg,legloc,'Fontsize', fs);
limy=get(gca,'ylim');
plot(ones(2,1)*kappa,limy,'--black')
kappa=bin*Z1./sum(Z1);
for i=1:length(kappa)
    kstr=sprintf('%0.2f',kappa(i));
    stdkap=sqrt((bin-kappa(i)).^2*Z1(:,i)/sum(Z1(:,i)));
    stdstr=sprintf('%0.2f',stdkap);
    text(0.7, 1-0.11*i, ['\alpha_',num2str(i),'=',kstr,'\pm',stdstr],...
        'Units','normalized','Backgr','white', 'Fonts', fs-2)
%     plot(ones(2,1)*kappa,limy,'--black')
end
xlim([-4 1]); title(tz1, 'Fonts', fs); xlabel(xlab1, 'Fonts', fs);
ylabel(ylab, 'Fonts', fs),grid
set(gca, 'Fonts', fs)

if size(fitrange,1)==2
    subplot(2,2,4), plot(bin,Z2); hold all;
    kappa=bin*Z2./sum(Z2); limy=get(gca,'ylim');
    for i=1:length(kappa)
        kstr=sprintf('%0.2f',kappa(i));
        stdkap=sqrt((bin-kappa(i)).^2*Z2(:,i)/sum(Z2(:,i)));
        stdstr=sprintf('%0.2f',stdkap);
        text(0.7, 1-0.11*i, ['\alpha_',num2str(i),'=',kstr,'\pm',stdstr],...
            'Units','normalized','Backgr','white', 'Fonts', fs-2)
        plot(ones(2,1)*kappa,limy,'--black')
    end
    xlim([-4 1]); title(tz2, 'Fonts', fs); xlabel(xlab2, 'Fonts', fs);
    ylabel(ylab, 'Fonts', fs),grid
    set(gca, 'Fonts', fs)
end

name=[saveto, 'kappa_stack_', version, '.eps'];
saveas(fig1, name, 'psc2')
name=[saveto, 'kappa_sep_', version, '.eps'];
saveas(fig2, name, 'psc2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
