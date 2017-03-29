%% Calculates Saturn's obliquity for a given day
%%
%% Input:   dnum        =   datenum(date)
%% Output:  theta_sun   =   angle between rotational axis of saturn and sun -90°
%%                          positive means southern summer

function [theta_sun] = satangle ( dnum )

%% Parameter
equi=datenum(2009,11,8);
solst=equi-10759.22/4; %from wikipedia ;)
obliq=26.7;

theta_sun = obliq * cos( (solst-dnum)/(solst-equi) *pi/2 );