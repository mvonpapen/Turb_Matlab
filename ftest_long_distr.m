%% Satistical significance tests for heating rate
% 
minB=2;
mm=1;
equi  = datenum(2009,8,11);

%% F-Test
clear k j i leg n p par SSE SST S VM VD Fn F Fc par orbdat m x y X Z ltm Prob
Z=log10(QLstrong);
X=psi_PS;
% Y=psi_PS_N;
norb = length(Leg);
j=find(rms>=dBthr & B>minB & vecB(:,3)<0 ...
    & abs(zdp)./HW<=1 & ~isnan(X) & ~isnan(Z) & r>6 & r<20 & inside_MP); % ...
    %& utcnum(:,1)>equi & (ltime>14 | ltime<10)); % & ~isnan(Y)
Prob   = NaN(norb,1);
orbdat = NaN(norb,1);
ltm    = NaN(norb,3);
n      = NaN(norb,1);
p      = 3;
% i=j; k=1; % In case of ALL data: Fn=24.8 for SLS4_lesia; Fn=49.4/23.6 for SLS3 (QL/QS)
% Example Rev 136 = [104 105];
for k=1:norb
    i=intersect(Leg{k},j);
    n(k) = length(i);
    if n(k)>=10
        orbdat(k)=mean(utcnum(i,1));
        fnum{k}=i;
        y{k}=Z(i);
        x{k}=X(i); %[5:10:355]';
%         x2{k}=Y(i);
        ltm(k,:)=[mod(meanangle(ltime(i)*15),360)/15, ...
            [min(ltime(i)') max(ltime(i)')]];

        % Fit parameter and create model
        % see also: mdl=fitnlm(x{k},y{k},@myfun,[-16.3 0.2 60])
        [par(k,:) r2norm(k)] = lsqcurvefit(@myfun,[-16.3 0.2 60],x{k},y{k});
%         [par(k,:) r2norm(k)] = lsqcurvefit(@model_sine_NS, ...
%             [0.1 300 0.1 300 -16], [x{k} x2{k}], y{k});
        if par(k,2)<0
            par(k,2)=-par(k,2);
            par(k,3)=par(k,3)-180;
        end
        m{k}=par(k,1)+par(k,2)*cosd(mm*(x{k}-par(k,3)));

        % Calculate sum of squares and (critical) F values
%             SST(k)=sum((y{k}-mean(y{k})).^2);
        SSE(k)=sum((y{k}-m{k}).^2);
        SSR(k)=sum((mean(y{k})-m{k}).^2);
        F(k)=SSR(k)*(n(k)-p)./SSE(k)/(p-1);
        Prob(k) = max( [1-fcdf(F(k),p-1,n(k)-p) ...
            fcdf(1/F(k),n(k)-p,p-1)]);
        VM(k)=sum((y{k}-m{k}).^2)/(n(k)-p);
        VD(k)=sum((y{k}-mean(y{k})).^2)/(n(k)-1);
    end
end
par(:,3)=mod(par(:,3),360);
par(par(:,1)==0,:) = NaN;
i0=find(~isnan(Prob));
i1=find(Prob<0.01);
clear p k j Z


%% Estimate standard error for binomial distribution
[h0,b]=hist(ltm(i0,1),[2:4:22]);
[h1,b]=hist(ltm(i1,1),[2:4:22]);
mu    = h1./h0;
mu_0  = sum(h1)/sum(h0);
sigma = sqrt(mu_0*(1-mu_0)./h0); % Eq.(2.7) Bevington2003
% Welch's t-Test
for i=1:6
    for j=1:6
        if i==j
            continue;
        end
        df(i,j) = (sigma(i)^2/h0(i)+sigma(j)^2/h0(j) )^2 / ...
            ( sigma(i)^4/h0(i)^2/(h0(i)-1) + sigma(j)^4/h0(j)^2/(h0(j)-1) );
        t(i,j)   = (mu(i)-mu(j))/sqrt( sigma(i)^2/h0(i)+sigma(j)^2/h0(j) );
    end
end
prob  = 1-tcdf(abs(t),df); % significant if prob<0.025


%% Plot results
figure
subplot(1,2,1)
bar(b,h0,'w'), hold all, bar(b,h1,'k')
legend('Total number (T)', 'Significant modulation (S)')
xlabel('Local time [h]')
ylabel('Number of in- and outbound legs')
xlim([0 24])
ylim([0 60])
subplot(1,2,2)
bar(b, h1./h0)
hold on, errorbar(b, h1./h0, sigma, '.')
xlabel('Local time [h]')
ylabel('Occurence rate S/T')
xlim([0 24])
ylim([0 1])
legend(['Avg occ.rate ' mat2str(round(sum(h1)/sum(h0)*100)) '%'])
[a,b] = meanangle(par(i1,3))