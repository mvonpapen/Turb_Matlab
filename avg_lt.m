%% Local time averaging

function [ltmean,ymean]=avg_lt(ltime,y,rminmax,stdB,zdp,scaleheight,B,r_cyl)

if nargin<3; rminmax=[6 20]; end

for i=1:24
    j=find(min(stdB,[],2)>=0.1 & max(stdB,[],2)<=3 & abs(zdp)./scaleheight<1 ...
        & r_cyl>min(rminmax) & r_cyl<max(rminmax) & ltime>(i-1) & ltime<i ...
        & ~isnan(y) & B>5);
    ymean(i,:)=[mean(y(j)) std(y(j))];
end

ltmean=[0.5:23.5];