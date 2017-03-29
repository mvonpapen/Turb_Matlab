%% Calculates Energy in a given frame (e.g. krho=[2,50]) as an estimate of
%% turbulence intensity
%%
%%
%% Input:
%%          x:          vector or matrix (same size as PSD) with x values for PSD
%%          PSD:        power spectral density for variable x
%%          snr:        signal to noise ratio as matrix with same size as PSD
%%          snrthres:   threshold at which signal is believed to be swamped by noise
%%          
%%  date: 19.02.2014



function E=turbint(x,PSD,frame,snr,snrthres)

if nargin<5, snrthres=5; end
if nargin<4, snr=5*ones(size(PSD)); end


%% Check sizes of vectors, matrices
a=size(x);
b=size(PSD);

% for i=1:length(b)
%     ja=find(a==b(i));
% end
% for i=1:length(a)
%     jb=find(b==a(i));
% end
% if ja==1, x=x'; end
% if jb==1, PSD=PSD'; end

%% Calculate Turbulence Intensity
% If x is the same for all PSD it's easy
if min(a)==1
    for i=1:b(1)
        j=[find(x>min(frame),1,'first'):...
            min([find(snr(i,:)<snrthres,1,'first') find(x>max(frame),1,'first')])-1];
        E(i)=nanmean(PSD(i,j)./PSD(1,j));
    end
% If not, we must first interpolate to common grid xi    
else    
    xi=logspace(log10(min(frame)),log10(max(frame)),10);
    for i=1:b(1)
        j=[find(x(i,:)>min(frame),1,'first'):...
            min([find(snr(i,:)<snrthres,1,'first') find(x(i,:)>max(frame),1,'first')])-1];
        if length(x(i,j))<2
            p(i,:)=NaN*ones(1,length(xi));
        else
            p(i,:)=interp1(x(i,j),PSD(i,j),xi);
        end
        E(i)=nanmean(p(i,:)./p(1,:));
    end
end