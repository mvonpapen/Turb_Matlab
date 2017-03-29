%% x=lsqcurvefit(@myfun,[0.2 240 0.9],SLS4(i),log10(Ekprho(i)))

function F = myfun(x,xdata)
m=1;
F = x(2)*cosd(m*(xdata-x(3)))+x(1);

% % function F = myfun(x,xdata)
% F = x(1)*exp(-(xdata-x(2)).^2/2/x(3)^2);
% F = log(x(1)*xdata.^(-8/3).*exp(-x(2)*xdata));

% function F = myfun(x,xdata)
% F = x(2)*xdata.^x(3)+x(1);