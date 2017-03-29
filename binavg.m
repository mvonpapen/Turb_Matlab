%% This function averages 'data' at locations 'r' in Rs bins of given length 'binlen'
%%
%% Date: 06.05.2013
%% Author: Michael von Papen


function [binavg,bins]=binavg(data,r,binlen)

imin=floor(min(r)/binlen);
imax=ceil(max(r)/binlen);

for i=imin:imax-1
    binavg(i-imin+1)=geomean(data(r>=i*binlen & r<(i+1)*binlen));
end

bins=binlen*[imin+0.5:1:imax-0.5];