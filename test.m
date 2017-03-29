for j=1:50
    for i=1:1136
        n2krho(i,j,:)=fitkappa(f,sum(PSD(i,:,:),3),f(krho(i,:)>j*0.1 & krho(i,:)<j*0.5 & f'<0.6 & squeeze(sum(PSD(i,:,:),3))>1e-3));
        n2kdi(i,j,:)=fitkappa(f,sum(PSD(i,:,:),3),f(kdi(i,:)>j*0.1 & kdi(i,:)<j*0.5 & f'<0.6 & squeeze(sum(PSD(i,:,:),3))>1e-3));
        n2kfic(i,j,:)=fitkappa(f,sum(PSD(i,:,:),3),f(f/ficO(i)>j*0.1 & f/ficO(i)<j*0.5 & f<0.6 & squeeze(sum(PSD(i,:,:),3))'>1e-3));
        n2kf(i,j,:)=fitkappa(f,sum(PSD(i,:,:),3),f(f>j*0.01 & f<j*0.05 & f<0.6 & squeeze(sum(PSD(i,:,:),3))'>1e-3));
    end;
end

