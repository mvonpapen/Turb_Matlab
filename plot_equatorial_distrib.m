%% Plots figures of polar plot for local time (top left), longitude (top
%% right), dayside longitude (bottom left) and nightside longitude
%% distribution (bottom right)

Z       = {[],                                log10(rms),      log10(QLstrong)};
Label   = {{'Orbits in equatorial plane';''}, {'RMS'; 'Southern'; 'Phase'}, {'Heating Rate'; 'Density'; 'Southern'; 'Phase'}};
Units   = {{''},                              {'log$_{10}(\delta B\,$[nT])'}, {'log$_{10}(q\,$[W/m$^3$])'}};
cbar    = [ 0 30;                             -0.5 0;           -17 -15.5];
X       = psi_PS;
lon     = 'magS';
name    = 'dist_S_post';
fs      = 14;

load MP_out
equi  = datenum(2009,8,11);
start = datenum(2009,8,11);
stop  = datenum(2019,8,11);
index = rms>0.1 & abs(zdp)./HW<1 & B>minB & r>6 & r<20 & vecB(:,3)<0 & inside_MP ...
    & (utcnum(:,1)>=start & utcnum(:,2)<=stop);

%% Make one figure for each element in Z
for i=1:length(Z);

    if isempty(Z{i})
        %Orbital coverage
        PS = [8 3];
        fig1 = figure('Visible', 'off', 'PaperPositionmode', 'manual', ...
            'PaperSize', PS, 'PaperPosition', [0 0 PS]);
        subplot(1,3,1)
        polar_plot_2015(r_cyl, ltime, 'orb', [], Z{i}, cbar(i,:), [0 24], index);
        title(Label{i})
        subplot(1,3,2)
        polar_plot_2015(r_cyl, ltime, 'lt',  [], Z{i}, cbar(i,:), [0 24], index);
        subplot(1,3,3)
        polar_plot_2015(r_cyl, ltime, lon, X, Z{i}, cbar(i,:), [0 24], index);
        ax=gca;
        pos=get(gca,'pos');
        set(gca,'pos',[pos(1) pos(2)*1.25 pos(3) pos(4)*0.95]);
        pos=get(gca,'pos');
        hc=colorbar('location','northoutside', ...
            'position',[pos(1)-0.25 pos(2)+pos(4)-0.02 pos(3)+0.25 0.04]);
        set(hc,'xaxisloc','top');
        print(fig1,'-depsc', [name mat2str(i) '.eps'],'-painters','-loose')
        close(fig1)
    else
        %Distribution of specified parameter
        PS = [8 8];
        fig1 = figure('Visible', 'off', 'PaperPositionmode', 'manual', ...
            'PaperSize', PS, 'PaperPosition', [0 0 PS]);
        subplot(2,2,1)
        polar_plot_2015(r_cyl, ltime, 'lt', X, Z{i}, cbar(i,:), [0 24], index);
        subplot(2,2,2)
        polar_plot_2015(r_cyl, ltime, lon, X, Z{i}, cbar(i,:), [0 24], index);
        subplot(2,2,3)
        polar_plot_2015(r_cyl, ltime, lon, X, Z{i}, cbar(i,:), [6 18], index);
        subplot(2,2,4)
        polar_plot_2015(r_cyl, ltime, lon, X, Z{i}, cbar(i,:), [18 6], index);
        ax=gca;
        pos=get(gca,'pos');
        hc=colorbar('location','northoutside', ...
            'position',[pos(1)-0.25 pos(2)+pos(4)+0.02 pos(3)+0.04 0.02]);
        set(hc,'xaxisloc','top');
        annotation(fig1,'textbox', [0.37 0.97 0.3 0.03], 'String', Label{i}, ...
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
        print(fig1,'-depsc', [name mat2str(i) '.eps'],'-painters','-loose')
        close(fig1)
    end
end