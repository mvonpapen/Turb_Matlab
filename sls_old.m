%% Calculate the longitude in Saturnian Longitude System 2 (SLS2) based on
%% Kurth et al., 2007 (GRL) for a given time/date and given local time of
%% the spacecraft (e.g. Cassini).
%%
%% Author: M. von Papen
%% Date: 04.10.2013
%%
%% Input: t     -    numeric date in **unit of days** with t=0 at year 0 (e.g. generated with datenum)
%%        lt    -   local time of spacecraft in dec.h
%%        type  -   SLS numer: 2=SLS2, 3=SLS3, 4=SLS4S, 5=SLS4N
%%        w0    -   Saturn's reference rotational period
    
function [phi]=sls_old(t,lt,type,w0,r)


if nargin<2; lt=12; end

% T0=1January2004
t0=datenum(2004,1,1);

% Which SLS? ( 2, 3 or 4 -> SLS-S, 5 -> SLS4-N )
if nargin<3; type=4; end

% Rotational period (synodic) in degrees/day
if nargin<4; w0=360/0.4497; end


if type==2 || type ==3

    % Coefficients in degree per n*day
    switch type
        case 2
            c0=100;
            c1=87.77;
            c2=-2.527;
            c3=3.041e-3;
            c4=-7.913e-7;
            c5=0;
            c6=0;
        case 3
            c0=100;
            c1=86.6681;
            c2=-2.7537;
            c3=4.773e-3;
            c4=-4.8755e-6;
            c5=3.5653e-9;
            c6=-9.1485e-13;
    end
    %% Calculation
    % Phase drift after Eq(4) of Kurth et al, 2008:
    phase=c1+c2.*(t-t0)+c3.*(t-t0).^2+c4.*(t-t0).^3+c5.*(t-t0).^4+c6.*(t-t0).^5;
    % Longitude after Eq(5) using Eq(3) for lambda_sun:
    phi=c0+w0*(t-t0)-phase + (12-lt)*15; %12*15°=180°

else
    %% Load SLS4 data from table
    load 'SLS4NS.mat'
    phi=NaN(1,length(t));
    switch type
        case 4
            % (Gurnett et al., 2011; UIowa) South
            for i=1:length(t)
                j=find(t(i)>=uiowa_south_zero,1,'last');
                if ~isempty(j)==1 && j+1<=length(uiowa_south_zero)
                    t0=uiowa_south_zero(j);
                    t1=uiowa_south_zero(j+1);
                    phi(i)=(t(i)-t0)/(t1-t0)*360 + (12-lt(i))*15;
                end
            end
        case 5
            % (Gurnett et al., 2011; UIowa) North
            for i=1:length(t)
                j=find(t(i)>=uiowa_north_zero,1,'last');
                if ~isempty(j)==1 && j+1<=length(uiowa_north_zero)
                    t0=uiowa_north_zero(j);
                    t1=uiowa_north_zero(j+1);
                    phi(i)=(t(i)-t0)/(t1-t0)*360 + (12-lt(i))*15;
                end
            end
        case 6
            % (Lamy, 2011; LESIA) South
            for i=1:length(t)
                j=find(t(i)>lesia_data(:,1),1,'last');
                if ~isempty(j)==1 && j+1<=length(lesia_data)
                    if lesia_data(j+1,4)-lesia_data(j,4)<0
                        lesia_data(j+1,4)=lesia_data(j+1,4)+360;
                    end
                    phi(i)=mean(lesia_data(j:j+1,4))+100 + (12-lt(i))*15;
                end
            end
        case 7
            % (Lamy, 2011; LESIA) North
            for i=1:length(t)
                j=find(t(i)>lesia_data(:,1),1,'last');
                if ~isempty(j)==1 && j+1<=length(lesia_data)
                    if lesia_data(j+1,4)-lesia_data(j,4)<0
                        lesia_data(j+1,4)=lesia_data(j+1,4)+360;
                    end
                    phi(i)=mean(lesia_data(j:j+1,4))+100 + (12-lt(i))*15;
                end
            end           
    end
end


%% Radial phase delay (see Provan et al., 2014, eq. 2b)
if nargin>3
    phi(r>12) = phi(r>12) - 3 * ( r(r>12)-12 );
end


%% Modulus 360
phi=mod(phi,360);