%% Loads flatness from wavelets and flatness from time series (or generates
%% it) and then plots the results in 1Rs bins together in one plot
%%
%% Date: 23.01.2014
%% Author: Mitch

load satturb8
n=length(files);

% %% Correct old path to files
% for i=1:1136;
%     files{i}=strrep(files{i},'PhD','PhD_archive');
% end

%% load new flatness from wavelet files
%% Comment out the following lines if you have up-to-date flatness
clear flat
for i=1:n
    dat=importpsd(strrep(files{i},'2wti','wti'));
    Fw(i,:,:)=dat(:,5:7);
end
clear i dat

%% calculate new flatness from time series files
%% Comment out the following lines if you have up-to-date flatness
clear kurt inc Ft
ninc=round(logspace(0.3,3.2,200));
ninc=unique(ninc);
n=length(files);
for i=1:n
    ts=importts(strrep(files{i},'2wti','dat'));
    for j=1:length(ninc)
        inc(:,1:3)=ts(ninc(j):end,2:4)-ts(1:end-ninc(j)+1,2:4);
        Ft(i,j,:)=kurtosis(inc);
        clear inc
    end
    i
end
clear i j ts

%% Average results in 1Rs bins
clear Ftmean Fwmean
rbin=[7:20];
for i=1:length(rbin)
    j=find(r>rbin(i)-0.5 & r<=rbin(i)+0.5);
    Fwmean(i,:,:)=mean(Fw(j,:,:));
    Ftmean(i,:,:)=mean(Ft(j,:,:));
end
%% Plot results in 1Rs bins
figure
i=[1:2:9]; % 7-15Rs
loglog(1./ninc/0.14,geomean(Ftmean(i,:,1:3),3),'-');
xlabel('1/\tau [s^{-1}], f [Hz]')
ylabel('Flatness F')
legend('7 R_s', '9 R_s', '11 R_s', '13 R_s', '15 R_s')
hold all
loglog(f,geomean(Fwmean(i,:,1:3),3),'--')