%% Plot normalized superposition of all PSD

% load satturb

% ki=[0.5,20];
% 
% for i=1:1136
%     fi=[find(f>2e-2,1,'first'):find(sn(i,:)'<10,1,'first')];
%     j=find(krho(i,:)>ki(1) & krho(i,:)<ki(2));
%     Pavg=geomean(sum(PSD(i,j,:),3),2);
%     loglog(krho(i,fi),sum(PSD(i,fi,:),3)/Pavg);
%     hold all
% end

% figure

n=401; %number of reference spectrum (should look normal)
thr=10; %SNR threshold
% snr=SNR;

% %% Normalized to refrence spectrum in f frame
% dx=[0.1 0.6];
% j=find(f>dx(1) & f<dx(2));
% y=sum(PSD,3);
% %% Calculate reference spectrum's energy
% fi0=find(f>2e-2,1,'first'):find(snr(n,:)<thr,1,'first')-1;
% E0=y(n,:);
% for i=1:1136
%     fi=find(f>2e-2,1,'first'):find(snr(i,:)<thr,1,'first')-1;
%     fi=intersect(intersect(fi,j),fi0);
%     E_f(i)=nanmean(y(i,fi)./E0(fi));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Normalized to refrence spectrum in krho frame
dx=[2 100];
xi=logspace(-1,2,50);
j=find(xi>dx(1) & xi<dx(2));
% ii=find(sqrt(sum(stdB.^2,2))'>=0.1 & sqrt(sum(stdB.^2,2))'<=3.1 ...% & offline<=0.14 ...
%     & abs(zdp)./scaleheight<10 & r>6 & r<20);% & kappa_kprho(:,2)'<0.1);
% ii=[1:5933];
y=sum(Pkprho,3);
x=kprho;
%% Calculate reference spectrum's energy
fi=find(f>2e-2,1,'first'):find(SNR(n,:)<thr,1,'first')-1;
E0=interp1(x(n,fi),y(n,fi),xi);
clear p E_krho
for i=1:length(r)
    fi=find(f>2e-2,1,'first'):find(SNR(i,:)<thr,1,'first')-1;
    if length(fi)>=2
        p(i,:)=interp1(x(i,fi),y(i,fi),xi);
    else
        p(i,:)=ones(length(xi),1)*NaN;
    end
    Ekprho(i)=nanmean(p(i,j)./E0(j));
end
i=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & r_cyl>6 & r_cyl<19 ...
    & abs(zdp)./scaleheight<1 & kappa_kprho(:,1)<-1);
for j=1:length(i); loglog(xi,p(i(j),:)/Ekprho(i(j))); hold all, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %% Normalized to refrence spectrum in f/ficW frame
% dx=[5 100];
% xi=logspace(-1,3,100);
% j=find(xi>dx(1) & xi<dx(2));
% clear x y
% for i=1:1136
%     y(i,:)=sum(PSD(i,:,:),3)*ficW(i);
%     x(i,:)=f/ficW(i);
% end
% 
% %% Calculate reference spectrum's energy
% fi=find(f>2e-2,1,'first'):find(snr(n,:)<thr,1,'first')-1;
% E0=interp1(x(n,fi),y(n,fi),xi);
% clear p E_fic
% for i=1:1136
%     fi=find(f>2e-2,1,'first'):find(snr(i,:)<thr,1,'first')-1;
%     if length(fi)>=2
%         p(i,:)=interp1(x(i,fi),y(i,fi),xi);
%     else
%         p(i,:)=NaN;
%     end
%     E_fic(i)=nanmean(p(i,j)./E0(j));
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Normalized to reference spectrum in kp-lambda frame
% dx=[5 100];
% xi=logspace(-1,3,100);
% j=find(xi>dx(1) & xi<dx(2));
% clear x y
% for i=1:1136
%     y(i,:)=sum(PSD(i,:,:),3)*norm(vrel(i,:))*sind(theta_vb(i))/2/pi/di(i);
%     x(i,:)=kpd(i,:);
% end
% 
% %% Calculate reference spectrum's energy
% fi=find(f>2e-2,1,'first'):find(snr(n,:)<thr,1,'first')-1;
% E0=interp1(x(n,fi),y(n,fi),xi);
% clear p E_kd
% for i=1:1136
%     fi=find(f>2e-2,1,'first'):find(snr(i,:)<thr,1,'first')-1;
%     if length(fi)>=2
%         p(i,:)=interp1(x(i,fi),y(i,fi),xi);
%     else
%         p(i,:)=NaN;
%     end
%     E_kd(i)=nanmean(p(i,j)./E0(j));
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Normalized to reference spectrum in k-frame
% dx=[5e-6 2e-5];
% xi=logspace(-7,-3,100);
% j=find(xi>dx(1) & xi<dx(2));
% fi=find(f>2e-2,1,'first'):find(snr(n,:)<thr,1,'first')-1;
% E0=interp1(krho(n,fi)./rhoi(n),sum(Pkrho(n,fi,:),3)*rhoi(n),xi);
% clear p E_k
% for i=1:1136
%     fi=find(f>2e-2,1,'first'):find(snr(i,:)<5,1,'first')-1;
%     if length(fi)>=2
%         p(i,:)=interp1(krho(i,fi)./rhoi(i),sum(Pkrho(i,fi,:),3)*rhoi(i),xi);
%     else
%         p(i,:)=NaN;
%     end
%     E_k(i)=nanmean(p(i,j)./E0(j));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot([3 100],[3 100].^(-2.6)*0.2,'--black','linew',2)
% xlim([0.1 200]), ylim([1e-8 1]), grid on
% ylabel('P/<P/P_0>_{k\rho}')
% xlabel('k_\perp\rho_W')

clear dx xi x y p n thr j i fi