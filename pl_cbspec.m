%% Program calculates the power spectral density (PSD) for critically
%% balanced turbulence according to Goldreich & Sridhar, 1995. The term
%% that handles critical balance is an exponential function as describen in
%% Cho, Lazarian & Vishniac 2002.

% ind=1;
% for os=[1 0.5 0.1 0.05 0.02 0.01 0.009 0.008 0.007 0.006 0.005 0.004 0.003 0.002 0.001 0.0009 0.0008 0.0005 0.0001];


clear P Pani kappa  
%% Einstellungen
where='horbury'; %'solarwind'; % 'saturn'; 'horbury'; 'test';
figout=0;%'P_SI_horb_KAW_0.003AU_v3.pdf' 'P_SI_SatKAW_v0.pdf' 'P_SI_sum_horb.pdf'; %['P_SI_sum_horb' mat2str(os) 'AU.pdf'];

%% Gridpunkte
n=1000;
fun=2; % 'expdamp'
psi=[0 5 15 25 35 45 55 65 75 85 90 95]; %[0:np-1]/(np-1); [0 5 15 25 35 45 55 65 75 85 90 95];

%% Turbulence
si=-10/3; %-10/3->k^{-5/3}
si2=-11/3; %-11/3->k^{-7/3}
cb=2/3; % 2/3->alfven
cb2=1/3; %1/3->KAW


%% Parameter
CS=149597870700; %Characteristic Scale, here: CS=1AU
switch where
    case 'test'
        v=6E5;
        L=0.0005*CS;
        kmin=log10(1/CS);
        rho=1e5;
        kmax=5;
        f=logspace(-4,3,15);
        xtick=10.^[-4:2:2];
        fit=[1e-3 1e-2;1e0 1e1];
        miny=1e-6;
        B=4e-9; %<B>_saturn=17nT
        va=60e3;
        beta_i=1;
    case 'saturn'
        CS=60268000; % CS=1Rs
        unit='R_s';
        v=60E3;
        va=60e3;
        L=1*CS;
        kmin=-10; %Umfang bei 17Rs: 6.4E9m->k0=-10
        kmax=0; %rho_e=rho_W/1840/18 ~ 3
        f=logspace(-6,2,30);
        xtick=10.^[-6:2:2];
        miny=1e-15;
        B=17e-9; %<B>_saturn=17nT
        rho=1e-2/B; %<rhoW>=3e5m aus meinen daten
        fit=[1e-4 1e-2;v/2/pi/rho*[1 10]];
        beta_i=1;
    case 'solarwind'
        v=6E5;
        unit='AU';
        L=0.005*CS;
        kmin=log10(1/CS);
        rho=1e5; %proton inertial length mit daten von sahraoui: 1E5m
        kmax=0;
        f=logspace(-6,2,20);
        xtick=10.^[-6:2:0];
        fit=[1e-5 1e-2;v/2/pi/rho*[1 10]];
        miny=1e-10;
        B=4e-9; %<B>_saturn=17nT
        va=60e3;
        beta_i=1;
    case 'horbury'
        unit='AU';
        L=0.003*CS;
        kmin=log10(1/1.44/CS);
        rho=1e4;
        kmax=0;
        f=fliplr(0.25*(5/8).^[0:9]);
        fit=[0.015 0.098];
        miny=1e-6;
        xtick=10.^[-2:-1];
        beta_i=1;
        load horbury_data
end



%% PSD berechnen
P=cbspectra(f,psi,L,rho,v,va,n,[kmin kmax],fun,[si si2],[cb cb2],B); % PSI in degrees!


%% Beginn der Schleife über Winkel
for k=1:length(psi)
    
    %% Berechne Variance Anisotropy
    Pani(:,k)=(P(1,:,k)+P(2,:,k))./P(3,:,k);
    
    
    %% Berechne spektralen Index
    [a,b]=size(fit);
    for jj=1:a
        for j=1:4
            [a1,b1]=fitkappa(f',P(j,:,k),fit(jj,:),2,0);
            kappa(j,k,jj)=a1(1);
        end
    end
    clear b a1 b1 j jj
        
    %% Berechne wave vector anisotropy after Bieber et al, 1996 & Hamilton et al, 2008
    for jj=1:a
        fn=find(f>=fit(jj,1) & f<=fit(jj,2));
        Pyx=geomean(P(2,fn,k)./P(1,fn,k));
        r(k,jj)=tan(psi(k))^(1+kappa(4,k,jj))*(Pyx-1)/(-kappa(4,k,jj)-Pyx)*(1-kappa(4,k,jj))/2;
    end
    r=1./(1+r);
    clear fn Pyx
    
end

clear k i i1 xi zi kern kx ky kz dky dkz



%% Erstelle visuelle Zusammenfassung
if figout ~=0
    fig=figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'Visible', 'Off',...
            'PaperUnits', 'normalized', 'PaperPosition', [0.05 0.05 0.9 0.9],...
            'PaperPositionMode', 'manual');

    if strcmp(where,'horbury')==1;
        subplot(3,1,1), plot(psi,kappa(1:4,:)','o')
        title(['L=',mat2str(round(L/CS*1e3)/1e3),' AU, v=',mat2str(v/1e3),' km/s'])
        legend('P_x_x', 'P_y_y', 'P_z_z', '|P|','Location', 'SE'), grid on
        xlabel('\psi \angle(v,B)'), ylabel('Spectral index'), xlim([0 95])
        hold all, errorbar(horbury_data(:,1),horbury_data(:,2),horbury_data(:,3));

        subplot(3,1,2), loglog(f,squeeze(P(4,:,:)))
        title(['log_1_0(k_m_i_n)=',mat2str(round(kmin*100)/100),' log_1_0(k_m_a_x)=',mat2str(round(kmax*100)/100)])
        xlim([min(f) max(f)]), ylim([miny max(max(P(4,:,:)))]), grid on
        xlabel('f [Hz]'), ylabel('PSD [au]'), set(gca,'xtick',xtick)

        subplot(3,1,3), plot(psi,squeeze(P(4,7,:)/P(4,7,1))), grid on
        xlabel('\psi \angle(v,B)'), ylabel('rel. Power at 61mHz'), xlim([0 95])

    else

        subplot(3,2,1), plot(psi,kappa(1:3,:,1)','o',psi,kappa(4,:,1)','-o')
        if a==2;
            hold all; A=get(gca,'ColorOrder');
            set(gca,'ColorOrder',A(1:4,:));
            plot(psi,kappa(1:3,:,2)','+',psi,kappa(4,:,2)','-+'); ylim([-6 0]);
        end
        title(['\kappa_1: ',mat2str(round(fit(1,1)*1e5)/1e2), '-',mat2str(round(fit(1,2)*1e5)/1e2),...
            'mHz, \kappa_2: ',mat2str(round(fit(2,1)*1e2)/1e2), '-',mat2str(round(fit(2,2)*1e2)/1e2),'Hz'])
        legend('P_x_x', 'P_y_y', 'P_z_z', '|P|','Location', 'SE'), grid on
        xlabel('\psi \angle(v,B)'), ylabel('Spectral index'), xlim([-5 95])

        subplot(3,2,2), plot(psi,r,'-o')
        title(['Wave vector anisotropy (Bieber et al., 1996)'])
        xlim([-5 95]), ylim([0 1]), grid on
        xlabel('\psi \angle(v,B)'), ylabel('r=C_S/(C_S+C_2)')
        legend('Ion Inertial range', 'Electron Inertial range','Location','NE')
        
        subplot(3,2,[3 4]), loglog(f,squeeze(P(4,:,:))')
        hold all, plot(v/2/pi/rho*[1 1],[miny max(max(P(2,:,:)))],'--black')
        title(['L=',mat2str(round(L/CS*100)/100),' ',unit,', v=',mat2str(v/1e3),' km/s, B=', ...
            mat2str(round(B*1e9)),' nT, \rho_W=',mat2str(round(rho/1e3)), ...
            ' km, [k_m_i_n,k_m_a_x]=10\^[',mat2str(round(kmin*10)/10),',',mat2str(kmax),']'])
        xlim([min(f) max(f)]), ylim([miny max(max(P(4,:,:)))]), grid on
        xlabel('f [Hz]'), ylabel('|P|'), set(gca,'xtick',xtick)

        subplot(3,2,5), loglog(f,squeeze(P(2,:,:)./P(1,:,:)))
        hold all, plot(v/2/pi/rho*[1 1],[0.9 10],'--black')
        xlim([min(f) max(f)]), ylim([0.9 10]), grid on
        xlabel('f [Hz]'), ylabel('P_(_\perp_)/P_(_|_|_)'), set(gca,'xtick',xtick)
        
        subplot(3,2,6), loglog(f,Pani)
        hold all, plot(v/2/pi/rho*[1 1],[0.9 10],'--black')
        xlim([min(f) max(f)]), ylim([0.9 10]), grid on
        xlabel('f [Hz]'), ylabel('P_\perp/P_|_|'), title('Variance Anisotropy'), set(gca,'xtick',xtick)

    end

    if figout ~= 1
        saveas(fig,figout)
        close(fig)
    else 
        set(fig, 'Visible', 'On')
    end
end



%% Folgendes für schleife über L=os*AU
% %% Store data in variable SI
% SI(ind).kappa=kappa(4,:);
% SI(ind).angle=psi*180/pi;
% SI(ind).outerscale=L;
% SI(ind).psd=squeeze(P(4,:,:));
% SI(ind).v=v; SI(ind).f=f; SI(ind).fit=fit;
% SI(ind).chi=sqrt(sum((kappa(4,:)-horbury_data(:,2)').^2./horbury_data(:,3)'.^2)/length(psi));
% SI(ind).rms=sqrt(sum((kappa(4,:)-horbury_data(:,2)').^2./horbury_data(:,2)'.^2)/length(psi));
% SI(ind).ngrid=n;
% 
% %increase index
% ind=ind+1;
%
% end