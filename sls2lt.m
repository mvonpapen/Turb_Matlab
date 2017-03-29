% Computes local time of a longitudinal location during SKR max,
% when subsolar longitude is 100°

function lt = sls2lt ( sls )

lt = 12 - (sls-100) / 15;
lt = mod(lt,24);