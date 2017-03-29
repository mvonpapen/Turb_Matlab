%% Calculates weighted mean and standard deviation of weighted mean
%%
%% Input: xvarin=[measured or mean values, standard deviation or measurement error of x]
%% Output: wm=weighted mean
%%         stdwm=weighted error of mean (this is usually strongly underestimated!)
%%
%% Author: M. von Papen, Date: 20.03.2013



function wm=wmean(varin)

[a,b]=size(varin);
if b~=2
    varin=varin';
end

if b==1
    'No std given!'
else
    y=~isnan(varin(:,1)) & ~isnan(varin(:,2));
    
    x=varin(y,1);
    stdx=varin(y,2);
    n=length(x);
        
    wm(1)=sum(x./stdx.^2)/sum(1./stdx.^2);
%     stdwm=sqrt(1/sum(1./stdx.^2));
%     stdwm=sqrt(1/sum(1./stdx.^2)*sum((x-wm).^2./stdx.^2)/(n-1));
    wm(2)=sqrt(var(x,1./stdx.^2));
end