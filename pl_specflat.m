clear zc1a zc1b zc2a zc2b zdat1 zdat2

fig1=figure('Visi', 'off');

%% Set Variables
x=log10(f);
y=d;
z1=log10(Pa'); z2=log10(Pe');
c1=log10(Fa'); c2=log10(Fe');

%% Set Ion cyclotron frequency
amu=1.66e-27;
e=1.6e-19;
fc=mB*e/amu/2/pi*1e-9;


fig1=figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'Visi', 'off',...
            'PaperUnits', 'normalized', 'PaperPosition', [0.05 0.1 0.9 0.85],...
            'PaperPositionMode', 'manual');
        
%% Plot 3D surface object
subplot(2,1,1)
h1=surf(x,y,z1,c1); zdat1=get(h1,'Zdata'); 
shading interp, colormap jet, hold all
view([8 -10 4]), caxis([0 1.3]), daspect([2 8 9])
title('Parallel Fluctuations');
xlabel('Log_1_0(f [Hz])'), ylabel('Distance [R_S]'), zlabel('PSD [nT^2/Hz]')
xlim([min(x) max(x)]), ylim([min(d) max(d)]), zlim([min(min([z1,z2])) max(max([z1,z2]))])

subplot(2,1,2)
h2=surf(x,y,z2,c2); zdat2=get(h2,'Zdata');
shading interp, colormap jet, hold all
view([8 -10 4]), caxis([0 1.3]), daspect([2 8 9])
title('Perpendicular Fluctuations');
xlabel('Log_1_0(f [Hz])'), ylabel('Distance [R_S]'), zlabel('PSD [nT^2/Hz]')
xlim([min(x) max(x)]), ylim([min(d) max(d)]), zlim([min(min([z1,z2])) max(max([z1,z2]))])

% Colorbar und Text
cb=colorbar('Location', 'South'); set(cb,'Position',[0.15,0,0.7,0.02])
text(-0.1,-0.21,'Log_1_0 (Local Intermittency Measure)','Units','norm');


%% Get zdata at ion cyclotron frequencies 
for i=1:length(y)
    zc1a(i)=interp1(x,zdat1(i,:),log10(fc(i)));
    zc1b(i)=interp1(x,zdat1(i,:),log10(fc(i)/18));
    zc2a(i)=interp1(x,zdat2(i,:),log10(fc(i)));
    zc2b(i)=interp1(x,zdat2(i,:),log10(fc(i)/18));
end


%% Plot ion cyclotron f-line on top of surface
subplot(2,1,1)
plot3(log10(fc),d,zc1a,'--white','Linew',2)
plot3(log10(fc/16),d,zc1b,'--white','Linew',2)
subplot(2,1,2)
plot3(log10(fc),d,zc2a,'--white','Linew',2)
plot3(log10(fc/16),d,zc2b,'--white','Linew',2)


% %% Plot f(d) on wall
% ax1n = axes('Position',get(ax1,'Position'),...
%            'XAxisLocation','top',...
%            'Color','none', 'XColor','blue');
% plot3(ones(length(y),1)*-4,y,(mB-mean(mB))./std(mB)*2+4,'--black');
% plot3(ones(length(y),1)*-4,y,varB+3,'-black','Linew',1.5);
% 
% w=1/(10.6*3600)*60000000*0.6;
% rho=1166*exp(-0.607*d)*amu*18*1e6;
% va=mB./sqrt(4e-7*pi*rho)*1e-9;
% plot3(ones(length(y),1)*-4,y,(va-mean(va))./std(va)*2+4,'-black','Linew', 1.5);

saveas(fig1,'test.pdf')
close all