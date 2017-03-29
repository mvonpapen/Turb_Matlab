%% Satistical significance tests for heating rate
% 
% load turb_at_saturn_v4
% x=[5:10:355]';
% % for k=1:7; % & r_cyl>(k-1)*2+6 & r_cyl<k*2+6
%     for i=1:36;
%         j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./HW<1 ...
%             & r_cyl>10 & r_cyl<15 & SLS4S_lesia>10*(i-1) & SLS4S_lesia<10*i ...
%             & ~isnan(QLstrong) & B>minB);
%         QL_sls4lesia_mean(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
%         j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./HW<1 ...
%             & r_cyl>10 & r_cyl<15 & SLS4S_lesia>10*(i-1) & SLS4S_lesia<10*i ...
%             & ~isnan(QSstrong) & B>minB);
%         QS_sls4lesia_mean(i,:)=[nanmean(log10(QSstrong(j))) nanstd(log10(QSstrong(j)))];        
% %         j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./HW<1 ...
% %             & r_cyl>(k-1)*2+6 & r_cyl<k*2+6 & SLS4S_ui>10*(i-1) & SLS4S_ui<10*i ...
% %             & ~isnan(QLstrong) & B>minB);
% %         QL_sls4ui_mean(i,:)=[nanmean(log10(QLstrong(j))) nanstd(log10(QLstrong(j)))];
% %         j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=2 & abs(zdp)./HW<1 ...
% %             & r_cyl>(k-1)*2+6 & r_cyl<k*2+6 & SLS4S_ui>10*(i-1) & SLS4S_ui<10*i ...
% %             & ~isnan(QSstrong) & B>minB);
% %         QS_sls4ui_mean(i,:)=[nanmean(log10(QSstrong(j))) nanstd(log10(QSstrong(j)))];        
%     end
%     plL(k,:)=lsqcurvefit(@myfun,[-16 0.2 300],x,QL_sls4lesia_mean(:,1));
% %     piL(k,:)=lsqcurvefit(@myfun,[-16 0.2 300],x,QL_sls4ui_mean(:,1));
%     plS(k,:)=lsqcurvefit(@myfun,[-16 0.2 300],x,QS_sls4lesia_mean(:,1));
% %     piS(k,:)=lsqcurvefit(@myfun,[-16 0.2 300],x,QS_sls4ui_mean(:,1));
%     figure
%     errorbar([5:10:355],QL_sls4lesia_mean(:,1),QL_sls4lesia_mean(:,2),'+black')
%     hold all
%     errorbar([5:10:355],QS_sls4lesia_mean(:,1),QS_sls4lesia_mean(:,2),'+red')
%     plot([0:360], plL(k,2)*cosd([0:360]-plL(k,3))+plL(k,1),'-black')
%     plot([0:360], plS(k,2)*cosd([0:360]-plS(k,3))+plL(k,1),'-red')
%     xlim([0 360]), ylim([-17 -15])
% end
% figure, plot([7:2:19],plL(:,3),[7:2:19],piL(:,3),[7:2:19],plS(:,3),[7:2:19],piS(:,3))    
% 
% 
% figure
% errorbar([5:10:355],QL_sls4lesia_mean(:,1),QL_sls4lesia_mean(:,2),'+--black')
% hold all
% errorbar([5:10:355],QL_sls4ui_mean(:,1),QL_sls4ui_mean(:,2),'+--red')
% plot([0:360], 0.16*cosd([0:360]-305)-16.1,'-blue')

% clear Ftest_lesia Ftest_ui
% 
minB=2;
mm=1;

%% F-Test
clear k j i leg n p par SSE SST S VM VD Fn F Fc par orbdat m x y X Z ltm
Z=log10(QLstrong);
X=MAGS;
j=find(min(stdB,[],2)>=dBthr(1) & max(stdB,[],2)<=dBthr(2) & B>minB ...
    & abs(zdp)./HT<=1 & ~isnan(X) & ~isnan(Z) & r<20);
leg=ones(182,1)*(-1);
Prob = NaN(106,1);
% i=j; k=1; % In case of ALL data: Fn=24.8 for SLS4_lesia; Fn=49.4/23.6 for SLS3 (QL/QS)
for k=1:182
    i=intersect(Orbs{k},j);
    if ~isempty(i)
        orbdat(k)=mean(utcnum(i,1));
        fnum{k}=i;
        y{k}=Z(i); %QL_sls4lesia_mean(:,1);
        x{k}=X(i); %[5:10:355]';
        n(k)=length(i);
        ltm(k,:)=[mod(meanangle(ltime(i)*15),360)/15, minmax(ltime(i)')];
        p=3;
        
        if n(k)>p

            % Fit parameter and create model
            [par(k,:) r2norm(k)] = lsqcurvefit(@myfun,[-16 0.1 300],x{k},y{k});
            if par(k,2)<0
                par(k,2)=-par(k,2);
                par(k,3)=par(k,3)-180;
            end
            m{k}=par(k,1)+par(k,2)*cosd(mm*(x{k}-par(k,3)));

            % Calculate sum of squares and (critical) F values
%             SST(k)=sum((y{k}-mean(y{k})).^2);
            SSE(k)=sum((y{k}-m{k}).^2);
            SSR(k)=sum((mean(y{k})-m{k}).^2);
            F(k)=SSR(k)./SSE(k)*(n(k)-p)/(p-1);
            Fc(k)=myfinv(0.99,p-1,n(k)-p); %critical F-value for 1%

            % Normalize F value and calculate variances for model and data to check
%             Fn(k)=F(k)/Fc(k);
            Prob(k) = 1-fcdf(F(k),p-1,n(k)-p);
            VM(k)=sum((y{k}-m{k}).^2)/(n(k)-p);
            VD(k)=sum((y{k}-mean(y{k})).^2)/(n(k)-1);

            if length(i)>1
                switch r(i(1))<r(i(end))
                    case 1
                        leg(k)=0; %=outbound
                    case 0
                        leg(k)=1; %=inbound
                end
            end
        end
    end
end
par(:,3)=mod(par(:,3),360);
i0=find(Prob<=1);
i1=find(Prob<=0.01);
clear p k j Z

%% Plot results
figure
[h0,b]=hist(ltm(i0,1),[2:4:22]);
[h1,b]=hist(ltm(i1,1),[2:4:22]);
subplot(1,2,1)
bar(b,h0,'w'), hold all, bar(b,h1,'k')
xlim([0 24])
subplot(1,2,2)
bar(b,h1./h0,'g')
xlim([0 24])
legend(['Avg occ.rate ' mat2str(round(sum(h1)/sum(h0)*100)) '%'])