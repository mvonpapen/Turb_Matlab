%% x=lsqcurvefit(@myfun,[0.2 240 0.9],SLS4(i),log10(Ekprho(i)))

function F = model_sine_NS(x,xdata)
F = x(1)*cosd(xdata(:,1)-x(2))+x(3)*cosd(xdata(:,2)-x(4))+x(5);

% % function F = myfun(x,xdata)
% % F = x(1)*exp(-(xdata-x(2)).^2/2/x(3)^2);
% F = log(x(1)*xdata.^(-8/3).*exp(-x(2)*xdata));

% function F = myfun(x,xdata)
% F = x(2)*xdata.^x(3)+x(1);