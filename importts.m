function [data, stats]=importts(fileToRead1)

%% Check ob file FFT oder WTF ist
suffix=fileToRead1(end-3:end);

switch suffix
    %% Wenn es DAT ist: importts
    case {'.dat'}


    DELIMITER = ' ';
    HEADERLINES = 15;

    % Import the file
    newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
    data=newData1.data;   
    % Orbit/Revelation
    stats.rev=textscan(newData1.textdata{1},...
        '# Cassini Orbit: %s %s');
    if isempty(stats.rev)
        stats.rev='';
        fprintf('Error in File %s: stats.rev is empty!\n', fileToRead1)
    end    
    % UTC start/stop
    stats.utc(1)=textscan(newData1.textdata{2},...
        '# First DP: %s');
    stats.utc(2)=textscan(newData1.textdata{3},...
        '# Last DP: %s');
    if isempty(stats.utc{1}) || isempty(stats.utc{1})
        fprintf('Error in File %s: stats.utc is empty!\n', fileToRead1)
    end
    % Distance
    stats.dist=textscan(newData1.textdata{4},...
        '# Distance: %f-%f%s R_S');
    if strcmp(stats.dist{3},'***')==1
        fprintf('File %s: Cassini <1Rs to Moon! (r=%3.1fRs)\n', fileToRead1,mean(cell2mat(stats.dist(1:2))))
        stats.dist=900+mean(cell2mat(stats.dist(1:2)));
    else
        stats.dist=mean(cell2mat(stats.dist(1:2)));
    end
    if isempty(stats.dist)
        stats.dist=999;
        fprintf('Error in File %s: stats.dist is empty!\n', fileToRead1)
    end
    % Z_DP
    stats.zdp=cell2mat(textscan(newData1.textdata{5},...
        '# Z_DP: %f - %f R_S'));
    if isempty(stats.zdp)
        stats.zdp=NaN;
        fprintf('Error in File %s: Z_DP is empty!\n', fileToRead1)
    end
    % Local Time
    stats.ltime=cell2mat(textscan(newData1.textdata{6},...
        '# Local Time: %f h'));
    if isempty(stats.ltime)
        stats.ltime=NaN;
        fprintf('Error in File %s: Local Time is empty!\n', fileToRead1)
    end    
    % meanB
    stats.meanB=cell2mat(textscan(newData1.textdata{7},...
        '# Mean field strength: %f nT'));
    if isempty(stats.meanB)
        stats.meanB=1E-99;
        fprintf('Error in File %s: stats.meanB is empty!\n', fileToRead1)
    end    
    % V_rel
    stats.vrel=cell2mat(textscan(newData1.textdata{8},...
        '# V_rel=V_pl-V_cas [r,phi,z]: %f %f %f km/s'));
    if isempty(stats.vrel)
        stats.vrel=[0,0,0];
        fprintf('Error in File %s: V_rel is empty!\n', fileToRead1)
    end
    %% vecB
    if strcmp(newData1.textdata{9}(32:40), '[r,phi,z]') ~= 1
        error('Error in File %s: vecB is not in RphiZ coordinates!\n', fileToRead1)
    end
    stats.vecB=cell2mat(textscan(newData1.textdata{12},...
        '# e_z: %f %f %f'));
    if isempty(stats.vecB)
        stats.meanB=[0 0 0];
        error('Error in File %s: stats.vecB is empty!\n', fileToRead1)
    end   
    % Time offline
    stats.offline=cell2mat(textscan(newData1.textdata{14},...
        '# Longest time offline [s]: MAX(dtnativ)=%f'));
    if isempty(stats.offline)
        stats.offline=9e9;
        fprintf('Error in File %s: t_offline is empty!\n', fileToRead1)
    end
    

    %% Wenn es 2DAT ist:
    case {'2dat'}

    DELIMITER = ' ';
    HEADERLINES = 14;

    % Import the file
    newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
    data=newData1.data;   
    % Orbit/Revelation
    stats.rev=textscan(newData1.textdata{1},...
        '# Cassini Orbit: %s');
    if isempty(stats.rev)
        stats.rev='';
        fprintf('Error in File %s: stats.rev is empty!\n', fileToRead1)
    end
    % UTC start/stop
    stats.utc(1)=textscan(newData1.textdata{2},...
        '# First DP: %s');
    stats.utc(2)=textscan(newData1.textdata{3},...
        '# Last DP: %s');
    if isempty(stats.utc{1}) || isempty(stats.utc{1})
        fprintf('Error in File %s: stats.utc is empty!\n', fileToRead1)
    end
    % Distance
    stats.dist=mean(cell2mat(textscan(newData1.textdata{4},...
        '# Distance: %f-%f R_S')));
    if isempty(stats.dist)
        stats.dist=999;
        fprintf('Error in File %s: stats.dist is empty!\n', fileToRead1)
    end
    % Z_DP
    stats.zdp=cell2mat(textscan(newData1.textdata{5},...
        '# Z_DP: %f - %f R_S'));
    if isempty(stats.zdp)
        stats.zdp=NaN;
        fprintf('Error in File %s: Z_DP is empty!\n', fileToRead1)
    end
    % meanB
    stats.meanB=cell2mat(textscan(newData1.textdata{6},...
        '# Mean field strength: %f nT'));
    if isempty(stats.meanB)
        stats.meanB=1E-99;
        fprintf('Error in File %s: stats.meanB is empty!\n', fileToRead1)
    end    
    % V_rel
    stats.vrel=cell2mat(textscan(newData1.textdata{7},...
        '# V_rel=V_pl-V_cas [r,phi,z]: %f %f %f km/s'));
    if isempty(stats.vrel)
        stats.vrel=[0,0,0];
        fprintf('Error in File %s: V_rel is empty!\n', fileToRead1)
    end
    %% vecB
    stats.vecB=cell2mat(textscan(newData1.textdata{11},...
        '# e_z: %f %f %f'));
    if isempty(stats.vecB)
        stats.meanB=[0 0 0];
        fprintf('Error in File %s: stats.vecB is empty!\n', fileToRead1)
    end   
    % Time offline
    stats.offline=cell2mat(textscan(newData1.textdata{13},...
        '# Longest time offline [s]: MAX(dtnativ)=%f'));
    if isempty(stats.offline)
        stats.offline=9e9;
        fprintf('Error in File %s: t_offline is empty!\n', fileToRead1)
    end
    
end