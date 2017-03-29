t=rand(1,10000);
n=length(t);
minsub=100; maxsub=n/5;
a=mean(t);
c=var(t);
maxlag=2000;
[ac,lags]=xcov(t,maxlag,'unbiased');
tc=sum(ac(maxlag/2+1:end-maxlag/5));

for i=minsub:maxsub
    clear m v
    for j=1:n-maxsub
        m(j)=mean(t(j:j+i));
        v(j)=var(t(j:j+i));
    end
    vm(i-minsub+1)=var(m-a);
    vv(i-minsub+1)=var(v-c);
end

loglog([minsub:maxsub],vm)
hold all
loglog([minsub:maxsub],2*c*tc./[minsub:maxsub])