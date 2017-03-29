function [data, stats]=importpsd(fileToRead1)

%% Check ob file FFT oder WTF ist
suffix=fileToRead1(end-3:end);

switch suffix
    %% Wenn es WTF ist: importwtf
    case {'.wtf'}

        DELIMITER = ' ';
        HEADERLINES = 17;

        % Import the file
        newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
        data=newData1.data;
        % UTC
        check=0;
        stats.utc(1)=textscan(cell2mat(newData1.textdata(5)),...
            '# First DP: %s');
        stats.utc(2)=textscan(cell2mat(newData1.textdata(6)),...
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
        % Velocity
        stats.vcas=cell2mat(textscan(cell2mat(newData1.textdata(7)),...
            '# V_Cassini in km/s [R,PHI,Z]: %f %f %f'));
        if isempty(stats.vcas)
            stats.vcas=[0 0 0];
            fprintf('Error in File %s: stats.vcas is empty!\n', fileToRead1)
        end
        % Distance
        stats.dist=mean(cell2mat(textscan(cell2mat(newData1.textdata(8)),...
            '# Distance: %f-%f R_S')));
        if isempty(stats.dist)
            stats.dist=999;
            fprintf('Error in File %s: stats.dist is empty!\n', fileToRead1)
        end
        % meanB
        stats.meanB=cell2mat(textscan(cell2mat(newData1.textdata(9)),...
            '# Mean field strength: %f'));
        if isempty(stats.meanB)
            stats.meanB=1E-99;
            fprintf('Error in File %s: stats.meanB is empty!\n', fileToRead1)
        end
        % Time offline
        stats.offline=cell2mat(textscan(cell2mat(newData1.textdata(10)),...
            '# Longest time offline [s]: MAX(dtnativ)=%fs'));
        if isempty(stats.offline)
            stats.offline=9999;
            fprintf('Error in File %s: stats.offline is empty!\n', fileToRead1)
        end
        % Standard deviation STD
        tmpi=textscan(cell2mat(newData1.textdata(11)),...
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
        
%%  Wenn es neues WTI ist:
    case {'.wti', '2wti', '3wti'}

        DELIMITER = ' ';
        HEADERLINES = 17;

        % Import the file
        newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
        data=newData1.data;
        % UTC
        check=0;
        stats.utc(1)=textscan(cell2mat(newData1.textdata(5)),...
            '# First DP: %s');
        stats.utc(2)=textscan(cell2mat(newData1.textdata(6)),...
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
        % Distance
        stats.dist=mean(cell2mat(textscan(cell2mat(newData1.textdata(7)),...
            '# Distance: %f-%f R_S')));
        if isempty(stats.dist)
            stats.dist=999;
            fprintf('Error in File %s: stats.dist is empty!\n', fileToRead1)
        end
        % Z (height over dipole equator)
        stats.z=mean(cell2mat(textscan(cell2mat(newData1.textdata(8)),...
            '# Z_DP: %f - %f R_S')));
        if isempty(stats.z)
            stats.z=999;
            fprintf('Error in File %s: z is empty!\n', fileToRead1)
        end
        % meanB
        stats.meanB=cell2mat(textscan(cell2mat(newData1.textdata(9)),...
            '# Mean field strength: %f nT'));
        if isempty(stats.meanB)
            stats.meanB=1E-99;
            fprintf('Error in File %s: stats.meanB is empty!\n', fileToRead1)
        end
        % Velocity
        stats.vcas=cell2mat(textscan(cell2mat(newData1.textdata(10)),...
            '# V_rel=V_pl-V_cas [r,phi,z]: %f %f %f km/s'));
        if isempty(stats.vcas)
            stats.vcas=[0 0 0];
            fprintf('Error in File %s: stats.vcas is empty!\n', fileToRead1)
        end
        % Time offline
        stats.offline=cell2mat(textscan(cell2mat(newData1.textdata(11)),...
            '# Longest time offline [s]: MAX(dtnativ)=%f'));
        if isempty(stats.offline)
            stats.offline=9999;
            fprintf('Error in File %s: stats.offline is empty!\n', fileToRead1)
        end
        % Standard deviation STD
        tmpi=textscan(cell2mat(newData1.textdata(12)),...
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
        % Orbit/Revelation
        stats.rev=textscan(cell2mat(newData1.textdata(16)),...
            '# Cassini Orbit: %s %s');
        if isempty(stats.rev)
            stats.rev='';
            fprintf('Error in File %s: stats.rev is empty!\n', fileToRead1)
        end
        
        
    %% Wenn es FFT ist: importpsd
    case {'.fft', '.0fft', '.lin', '.inc', '.hinc', '.psd', '.pol'}
        
        DELIMITER = ' ';
        HEADERLINES = 16;
        % Import the file
        newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
        data=newData1.data;
        % UTC
        check=0;
        stats.utc(1)=textscan(cell2mat(newData1.textdata(4)),...
            '# First DP: %s');
        stats.utc(2)=textscan(cell2mat(newData1.textdata(5)),...
            '# Last DP: %s');
        for i=1:2
            if isempty(stats.utc{i})
                check=1;
            end
        end
        if check==1
            fprintf('Error in File %s: stats.utc is empty!\n', fileToRead1)
        end
        % Distance
        stats.dist=mean(cell2mat(textscan(cell2mat(newData1.textdata(6)),...
            '# Distance: %f-%f R_S')));
        if isempty(stats.dist)
            stats.dist=999;
            fprintf('Error in File %s: stats.dist is empty!\n', fileToRead1)
        end
        % meanB
        stats.meanB=cell2mat(textscan(cell2mat(newData1.textdata(7)),...
            '# Mean field strength: %f'));
        if isempty(stats.meanB)
            stats.meanB=1E-99;
            fprintf('Error in File %s: stats.meanB is empty!\n', fileToRead1)
        end
        % var
        tmp=textscan(cell2mat(newData1.textdata(8)),...
            '# STD of time series (x,y,z) [nT]: %f %f %f');
        check=0;
        for i=1:3
            if isempty(tmp{i})
                check=1;
            end
        end
        switch check
            case 1
                stats.std=[0 0 0];
                fprintf('Error in File %s: stats.std is empty!\n', fileToRead1)
            case 0
                stats.std=cell2mat(tmp);
        end        
        % Time offline
        stats.offline=cell2mat(textscan(cell2mat(newData1.textdata(10)),...
            '# Longest time offline [s]: MAX(dtnativ)=%fs'));
        if isempty(stats.offline)
            stats.offline=9999;
            fprintf('Error in File %s: stats.offline is empty!\n', fileToRead1)
        end
        
    %% Wenn es gws ist:
    case {'.gws', '.dat'}

        DELIMITER = ' ';
        HEADERLINES = 9;

        % Import the file
        newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
%         %% Maxdt
%         stats.dt=cell2mat(textscan(cell2mat(newData1.textdata(2)),...
%             '# dt= %f s,  maxdt= %f s'));
%         if isempty(stats.dt)
%             stats.dt=[999 999];
%             fprintf('Error in File %s: dt is empty!\n', fileToRead1)
%         end
        % Distance
        stats.dist=mean(cell2mat(textscan(cell2mat(newData1.textdata(3)),...
            '# Distance= %f - %f')));
        if isempty(stats.dist)
            stats.dist=999;
            fprintf('Error in File %s: dist is empty!\n', fileToRead1)
        end        
        % Z (height over dipole equator)
        stats.z=mean(cell2mat(textscan(cell2mat(newData1.textdata(4)),...
            '# Z= %f - %f')));
        if isempty(stats.z)
            stats.z=999;
            fprintf('Error in File %s: z is empty!\n', fileToRead1)
        end
        % Vcas
        stats.vcas=cell2mat(textscan(cell2mat(newData1.textdata(5)),...
            '# Vcas [r,phi,z]= %f %f %f km/s'));
        if isempty(stats.vcas)
            stats.vcas=999;
            fprintf('Error in File %s: vcas is empty!\n', fileToRead1)
        end
        % MeanB Vector
        stats.Bvec=cell2mat(textscan(cell2mat(newData1.textdata(6)),...
            '# B [r,phi,z]= %f %f %f nT'));
        if isempty(stats.Bvec)
            stats.Bvec=999;
%             fprintf('Error in File %s: Bvec is empty!\n', fileToRead1)
        end
        % meanB
        stats.meanB=cell2mat(textscan(cell2mat(newData1.textdata(7)),...
            '# MeanB= %f'));
        if isempty(stats.meanB)
            stats.meanB=1E-9;
            fprintf('Error in File %s: meanB is empty!\n', fileToRead1)
        end
        % varB
        stats.varB=cell2mat(textscan(cell2mat(newData1.textdata(8)),...
            '# Var(|B|)= %f nT^2'));
        if isempty(stats.varB)
            stats.varB=1E9;
%             fprintf('Error in File %s: varB is empty!\n', fileToRead1)
        end        
        % Read data with textscan, so that NaN won't be ignored
        fid=fopen(fileToRead1,'r');
        a=textscan(fid,'%s %s %s %s %s', 'Headerlines', HEADERLINES);
        for i=1:5
            data(:,i)=str2double(a{i});
        end
        fclose(fid);
end