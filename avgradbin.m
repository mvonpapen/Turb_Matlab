%% Averages in radial 0.5 Rs bins

function avg=avgradbin(r,x)

for i=1:22
    avg(i)=geomean(x(r>=6.5+(i-1)*0.5 & r<6.5+i*0.5));
end