%% Polar plot of local time distribution colour coded with spectral index

function h = polar_plot(r, ltime, order, lon, Z, cbar, lt, index)

% figure
%% Choose Parameter to plot and corresponding color bar %%%%%%%%%
% order='sls';
% lon=MAGS';
% Z=kappa_kprho(:,1); %log10(QLstrong);%log10(QLstrong); %log10(QSstrong); %empty Z gives number of points in each bin   %log10(Ekprho);  %log10(QLstrong);
% cbar=[-3 -2];   %[0 2];   %[-17 -15.5]; [0 50]
% lt=[0 24];
% equi=datenum(2009,11,8);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Filter out lobe data and spikes
if lt(1)<lt(2)
    i=find(  ltime>=lt(1) & ltime<=lt(2) & index);
else
    i=find( (ltime>=lt(1) | ltime<=lt(2)) & index);
end
if isempty(Z)~=1
    Z=Z(i);
end

% Plot polar plot without points
% figure
rmax=20;
h = polar(0, rmax,'');
hold all

% Plot data with scatter, so that can be color coded
switch order
    case 'lt' % Local time
        ri=1;ti=0.5;
        xi=[6:ri:rmax];
        yi=[0:ti:24];
        zi=zeros(length(xi),length(yi));
        nzi=zi; X=zi; Y=zi;
        for j=1:length(xi)
            for jj=1:length(yi)-1
                ii=find(r(i)>=xi(j) & r(i)<xi(j)+ri & ...
                    ltime(i)>=yi(jj) & ltime(i)<yi(jj)+ti);
                nzi(j,jj)=length(ii);
                if ~isempty(Z)==1
                    zi(j,jj)=nanmean(Z(ii));
                end
            end
            % Define grid in cartesian coordinates X Y
            X(j,:)=-xi(j)*sin(-yi/12*pi);
            Y(j,:)=-xi(j)*cos(-yi/12*pi);
        end
        if isempty(Z)==1
            zi=nzi;
            zi(nzi==0)=NaN;
        end
        zi(:,end)=zi(:,1); %data at 0h equal to 24h LT
        pcolor(X,Y,zi), shading flat
%         scatter(X,Y,30,Z,'filled');
    case {'sls', 'magS', 'magN'} % SLS
        ri=1;ti=5;
        xi=[6:ri:rmax];
        yi=[0:ti:360];
        zi=zeros(length(xi),length(yi));
        nzi=zi; X=zi; Y=zi;
        for j=1:length(xi)
            for jj=1:length(yi)-1
                ii=find(r(i)>=xi(j) & r(i)<xi(j)+ri & ...
                    lon(i)>=yi(jj) & lon(i)<yi(jj)+ti);
                nzi(j,jj)=length(ii);
                if ~isempty(Z)==1
                    zi(j,jj)=nanmean(Z(ii));
                end
            end
            % Define grid in cartesian coordinates X Y
            X(j,:)=-xi(j)*sind(yi);
            Y(j,:)=-xi(j)*cosd(yi);
        end
        if isempty(Z)==1
            zi=nzi;
            zi(nzi==0)=NaN;
        end
        zi(:,end)=zi(:,1); %data at 0h equal to 24h LT
        pcolor(X,Y,zi), shading flat
%         alpha(nzi)
%         alim([0 10])
    case 'orb' % plot orbits
        X=-r(i).*sin(-ltime(i)/12*pi);
        Y=-r(i).*cos(-ltime(i)/12*pi);
        plot(X,Y,'.k', 'Markersize', 1);
end

%% Configure axes and labels
ph=findall(gca,'type','text');
ps=get(ph,'string');
% disp([num2cell(1:numel(ps)).',ps]);
% ps([1:end])={''}; %theta and rho values
switch order
    case {'lt', 'orb'}
        ps(1:17) = {'6', '18', '4', '16', '2', '14', '0', '12', '22', ...
            '10', '20', '8', '', '', '', 'Local time', ''};
    case {'sls', 'sls4s'}
        ps(1:17) = {'270', '90', '300', '120', '330', '150', '0', '180', '30', ...
            '210', '60', '240', '', '', '', 'SLS', ''};
    case {'magS'}
        ps(1:17) = {'270', '90', '300', '120', '330', '150', '0', '180', '30', ...
            '210', '60', '240', '', '', '', '\Psi_{ps,s}', ''};
    case {'magN'}
        ps(1:17) = {'270', '90', '300', '120', '330', '150', '0', '180', '30', ...
            '210', '60', '240', '', '', '', '\Psi_{ps,n}', ''};
end

set(ph,{'string'},ps);
colormap jet
caxis(cbar)
axis tight


%% Overplot 10Rs and 20Rs circle and cross
X=10*sin([0:0.1:2*pi+0.1]);
Y=10*cos([0:0.1:2*pi+0.1]);
plot(X,Y,'-k', 'Linew', 1);
X=20*sin([0:0.1:2*pi+0.1]);
Y=20*cos([0:0.1:2*pi+0.1]);
plot(X,Y,'-k', 'Linew', 2);
plot([-20 20],[0 0],'-k', 'Linew', 1);
plot([0 0],[-20 20],'-k', 'Linew', 1);

clear ps ph x y i ii j jj lon ans ri ti xi yi X Y rmax nzi zi