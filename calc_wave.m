%% Load data from variable DATA and calculate wavelet transforms and
%% coherency
clear all
load '/afs/geo.uni-koeln.de/usr/neuro/DATA/DATA_257_last13.3.mat'
[t, dat, ndat, rigor, CHlab]=...
    load_dat_from_DATA ( DATA, 'value1', 'Ruhe', 'activity', 2);

FF=[250 350];

%% Set variables
nf=3;
t_smo=3; % number of scales used to smooth
sc_smo=0; % number of neighboring scales to smooth in freq-domain
ns=length(ndat);
datf=cell(ns,1);
WLFP=cell(ns,1);
coiLFP=cell(ns,1);
PLFP=NaN(nf,5,ns);
PEMG=NaN(nf,3,ns);
WEMG=cell(ns,1);
coiEMG=cell(ns,1);
T=cell(ns,1);
f=logspace(log10(FF(1)),log10(FF(2)),nf);
dt=1/2500;

%% Define LFP and EMG channels
LFPch={'C', 'L', 'P', 'A', 'M'};
EMGch={'EDCre', 'EDCli', 'FDIre', 'FDIli', 'FDLre', 'FDLli'};

for n=1:ns
    
    site(n)=DATA(ndat(n)).site;
    
    %% Check time vector
    for i=1:length(t{n})
        ti(i,:)=minmax(t{n}{i});
    end
    ti=[max(ti(:,1)) min(ti(:,2))];
    T{n} = t{n}{i}( t{n}{i}>=ti(1) & t{n}{i}<=ti(2) );
    
    %% LFP data
    iLFP=find(ismember(CHlab{n},LFPch));
    nLFP(n)=length(iLFP);
    for i=iLFP
        ni =  t{n}{i}>=ti(1) & t{n}{i}<=ti(2) ;
        [datf{n}(:,i), WLFP{n}(:,:,i), coiLFP{n}(:,i), PLFP(:,i,n)]...
            = procdata (dat{n}{i}(ni), 'freq', f, 'filter', []); 
    end    
        
%     %% EMG data
%     iEMG=find(ismember(CHlab{n},EMGch));
%     nEMG(n)=length(iEMG);
%     for i=iEMG;
%         ni =  t{n}{i}>=ti(1) & t{n}{i}<=ti(2) ;
%         [datf{n}(:,i),WEMG{n}(:,:,i-nLFP(n)),coiEMG{n}(:,i-nLFP(n)),PEMG(:,i-nLFP(n),n)]...
%             =procdata(dat{n}{i}(ni),'filter', [],...
%             'rect', 0, 'freq', f);
%     end
% 
%     
%     %% Coherency
%     if ~isempty(WLFP{n}) && ~isempty(WEMG{n})
%         [CLE{n}, WLE{n}] = wave_coh(WLFP{n}, WEMG{n},f,dt, t_smo, sc_smo);
%     end

    
end
clear ni nia nib i j n ti ni2

