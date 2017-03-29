% Polar plots for coverage (top row), and turbulence heating rates
% (bottom row) for pre- and post-equinox (left and right column, resp.)

Z       = {[],              log10(QLstrong)};
Label   = {{'Pre-equinox'}, {'Post-equinox'}};
Units   = {{''},            {'log$_{10}(q\,$[W/m$^3$])'}};
cbar    = [ 0 30;           -17 -15.5];
X       = psi_PS;
lon     = 'magS';
name    = 'dist_S_prepost';
fs      = 14;
lt_ana  = [14 10]; % Local time sector to analyze, [14 10]: from 14-24LT and 0-10LT

load MP_out
equi  = datenum(2009,8,11);
pre   = rms>0.1 & abs(zdp)./HW<1 & B>minB & r>6 & r<20 & vecB(:,3)<0 & inside_MP ...
    & (utcnum(:,1)<equi);
post  = rms>0.1 & abs(zdp)./HW<1 & B>minB & r>6 & r<20 & vecB(:,3)<0 & inside_MP ...
    & (utcnum(:,1)>equi);


%% Make one figure for each element in Z
PS = [8 8];
fig1 = figure('Visible', 'off', 'PaperPositionmode', 'manual', ...
    'PaperSize', PS, 'PaperPosition', [0 0 PS]);

for i=1:length(Z);

    if isempty(Z{i})
        %Orbital coverage
        subplot(2,2,1)
        polar_plot_2015(r_cyl, ltime, 'orb', [], Z{i}, cbar(i,:), [0 24], pre);
        subplot(2,2,2)
        polar_plot_2015(r_cyl, ltime, 'orb', [], Z{i}, cbar(i,:), [0 24], post);
    else
        %Distribution of specified parameter
        subplot(2,2,3)
        polar_plot_2015(r_cyl, ltime, lon, X, Z{i}, cbar(i,:), lt_ana, pre);
        subplot(2,2,4)
        polar_plot_2015(r_cyl, ltime, lon, X, Z{i}, cbar(i,:), lt_ana, post);
    end
end

ax=gca;
pos=get(gca,'pos');
hc=colorbar('location','northoutside', ...
    'position',[pos(1)-0.25 pos(2)+pos(4)+0.02 pos(3)+0.04 0.02]);
set(hc,'xaxisloc','top');
annotation(fig1,'textbox', [0.16 0.97 0.3 0.03], 'String', Label{1}, ...
    'Fontsize', 18, 'HorizontalAlign', 'center', 'Linestyle', 'none', ...
    'Fontweight', 'bold');
annotation(fig1,'textbox', [0.6 0.97 0.3 0.03], 'String', Label{2}, ...
    'Fontsize', 18, 'HorizontalAlign', 'center', 'Linestyle', 'none', ...
    'Fontweight', 'bold');
annotation(fig1,'textbox', [0.37 0.52 0.3 0.03], 'String', Units{i}, ...
    'Fontsize', fs, 'HorizontalAlign', 'center', 'Linestyle', 'none', ...
    'Interpreter', 'latex');
annotation(fig1,'textbox', [0.03 0.9 0.2 0.03], 'String', 'A', ...
    'Fontsize', fs, 'HorizontalAlign', 'center', 'Linestyle', 'none', ...
    'Fontweight', 'bold');
annotation(fig1,'textbox', [0.82 0.9 0.2 0.03], 'String', 'B', ...
    'Fontsize', fs, 'HorizontalAlign', 'center', 'Linestyle', 'none', ...
    'Fontweight', 'bold');
annotation(fig1,'textbox', [0.03 0.42 0.2 0.03], 'String', 'C', ...
    'Fontsize', fs, 'HorizontalAlign', 'center', 'Linestyle', 'none', ...
    'Fontweight', 'bold');
annotation(fig1,'textbox', [0.82 0.42 0.2 0.03], 'String', 'D', ...
    'Fontsize', fs, 'HorizontalAlign', 'center', 'Linestyle', 'none', ...
    'Fontweight', 'bold');
print(fig1,'-depsc', [name '.eps'],'-painters','-loose')
% close(fig1)