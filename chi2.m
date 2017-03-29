%% Calculates the Chi squared error

function [chi]=chi2(x,dat,err)

x=x(:); dat=dat(:); err=err(:);

chi=sum( sqrt( ( (x-dat)./err ).^2 / length(dat) ) );