%% Performs a least quares fit to the double logarithmic plot of y=f(x) in
%% the range determined by fit.


function [a,b]=fitkappa(x,y,fit,mindp,noise)


if nargin<5
    noise=0; %Upper boundary for fit, when PSD reaches noise level
end

if nargin<4
    mindp=2;
end
if length(fit)>2, fit=[fit(1) fit(end)]; end

if length(fit)==2 && fit(1)<fit(2)
    %% Bring x and y in the right form. x must be monotonically increasing!
    [i1,i2]=sort(x);
    x=x(i2); y=y(i2);
    x=x(:); y=y(:);
    clear i1 i2

    %% Spectral Indices
    %%Fit ranges bestimmen
    % f=[find(x>=min(fit) & x<=max(fit))];
    f=[find(x>=fit(1), 1, 'first') ...
         : min([find(y<noise,1,'first') ...
                find(x<=fit(2), 1, 'last')])];

    %% Linear least squares fitting yl=a*xl+b            
    if length(f)>=mindp
        xl = log10(x(f));
        yl = log10(y(f));
        yl = yl(~isnan(yl));
        xl = xl(~isnan(yl));
        if length(xl)>=mindp && length(yl)>=mindp
            n=length(xl);
            %% With Matlab Statistics Toolbox (regress)
%             [v, sd]=regress(yl,[xl ones(length(xl),1)]);
%             a=[v(1) sd(1,2)-sd(1,1)];
%             b=[v(2) sd(2,2)-sd(2,1)];
%             c=sqrt( sum( (yl-a(1)*xl-b(1)).^2 ) / (n-2) );
            %% With Gaussian method
            Sx=sum(xl);
            Sy=sum(yl);
            Sxx=sum(xl.^2);
            Sxy=sum(xl.*yl);
            n=length(xl);
            delta=(n*Sxx-Sx^2);
            % Values
            a(1)=(n*Sxy-Sx*Sy)/delta;
            b(1)=Sy/n-a(1)*Sx/n;
            % Standard deviations
            dy=sqrt( sum( (yl-a(1)*xl-b(1)).^2 ) / (n-2) );
            a(2)=sqrt(dy^2*n/delta);
            b(2)=sqrt(dy^2*Sxx/delta);
        else
            a=[NaN NaN];
            b=[NaN NaN];
        end
    else
        a=[NaN NaN];
        b=[NaN NaN];
    end
else
    a=[NaN NaN];
    b=[NaN NaN];
end