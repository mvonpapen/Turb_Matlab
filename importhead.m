function stats=importhead(file)

DELIMITER = ' ';
HEADERLINES = 12;

% Import the file
head = importdata(file, DELIMITER, HEADERLINES);

%% UTC
check=0;
stats.utc(1)=textscan(cell2mat(head(3)),...
    '# First DP: %s');
stats.utc(2)=textscan(cell2mat(head(4)),...
    '# Last DP: %s');
for i=1:2
    if isempty(stats.utc{i})
        check=1;
    end
end
if check==1
    fprintf('Error in File %s: stats.utc is empty!\n', fileToRead1)
    stats.utc={{'2000-01-01T00:00:00'}, {'2020-01-01T00:00:00'}};
end
%% Velocity
stats.vcas=cell2mat(textscan(cell2mat(head(5)),...
    '# V_Cassini in km/s [R,PHI,Z]: %f %f %f'));
if isempty(stats.vcas)
    stats.vcas=[0 0 0];
    fprintf('Error in File %s: stats.vcas is empty!\n', fileToRead1)
end
%% Distance
stats.dist=mean(cell2mat(textscan(cell2mat(head(6)),...
    '# Distance: %f-%f R_S')));
if isempty(stats.dist)
    stats.dist=999;
    fprintf('Error in File %s: stats.dist is empty!\n', fileToRead1)
end
%% meanB
stats.meanB=cell2mat(textscan(cell2mat(head(7)),...
    '# Mean field strength: %f'));
if isempty(stats.meanB)
    stats.meanB=1E-99;
    fprintf('Error in File %s: stats.meanB is empty!\n', fileToRead1)
end
%% Time offline
stats.offline=cell2mat(textscan(cell2mat(head(8)),...
    '# Longest time offline [s]: MAX(dtnativ)=%fs'));
if isempty(stats.offline)
    stats.offline=9999;
    fprintf('Error in File %s: stats.offline is empty!\n', fileToRead1)
end
%% Standard deviation STD
tmpi=textscan(cell2mat(head(9)),...
    '# STD of time series (x,y,z) [nT]: %f %f %f');
check=0;
for i=1:3
    if isempty(tmpi{i})
        check=1;
    else
        tmp{i}=tmpi{i};
    end
end
switch check
    case 1
        stats.std=[0 0 0];
        fprintf('Error in File %s: stats.std is empty!\n', fileToRead1)
    case 0
        stats.std=cell2mat(tmp);
end
%% Energy in time domain
tmpi=textscan(cell2mat(head(10)),...
    '# Energy in timedom (x,y,z): %f %f %f');
check=0;
for i=1:3
    if isempty(tmpi{i})
        check=1;
    else
        tmp{i}=tmpi{i};
    end
end
switch check
    case 1
        stats.etd=[0 0 0];
        fprintf('Error in File %s: stats.etd is empty!\n', fileToRead1)
    case 0
        stats.etd=cell2mat(tmp);
end
%% Energy in wtf
tmpi=textscan(cell2mat(head(11)),...
    '# Energy in WTF (x,y,z): %f %f %f');
check=0;
for i=1:3
    if isempty(tmpi{i})
        check=1;
    else
        tmp{i}=tmpi{i};
    end
end
switch check
    case 1
        stats.ewtf=[0 0 0];
        fprintf('Error in File %s: stats.ewtf is empty!\n', fileToRead1)
    case 0
        stats.ewtf=cell2mat(tmp);
end
%% Local time:
stats.lt=cell2mat(textscan(cell2mat(head(12)),...
    '# Local time: %f h'));
if isempty(stats.lt)
    stats.lt=NaN;
    fprintf('Error in File %s: stats.lt is empty!\n', fileToRead1)
end