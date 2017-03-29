%% Calculate roots of third degree polynomial using two slightly different
%% methods that analytically should give the same result
%%
%%
%% INPUT:
%%          a,b,c,d:    Coefficients of ax^3+bx^2+cx+d=0
%%
%% OUTPUT:
%%          x1:         Three roots from method 1 (Cardano's method)
%%          x2:         Three roots from method 2

function [x1,x2]=p3root(a,b,c,d)

D = 18.*a.*b.*c.*d - 4*b.^3.*d + b.^2.*c.^2 - 4.*a.*c.^3 - 27*a.^2.*d.^2;
D0 = b.^2 - 3*a*c;
D1 = 2*b.^3 - 9*a.*b.*c + 27*a.^2.*d;

C = ( ( D1 + sqrt(D1.^2-4*D0.^3)) / 2 ).^(1/3);
Ci = ( ( D1 - sqrt(D1.^2-4*D0.^3)) / 2 ).^(1/3);

u(1) = 1;
u(2) = -(1-sqrt(-3))/2;
u(3) = -(1+sqrt(-3))/2;

for j=1:3
    x1(j,:)=-1/3./a.*(b+u(j).*C+D0./u(j)./C); % this result is more stable/less numerical error
    x2(j,:)=-1/3./a.*(b+u(j).*C+conj(u(j)).*Ci);
end